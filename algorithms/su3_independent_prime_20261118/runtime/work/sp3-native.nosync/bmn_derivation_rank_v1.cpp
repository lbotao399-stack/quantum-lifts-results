// BMN adapter of the existing native graded product rule; prime set at compile time.
// Evaluate the pinned Q_letter rules with the graded matrix product rule.
// All coefficients and the ordering of the rules are supplied by Python.
#include "trace_frame_bmn_rank_v1.cpp"
#include <array>

struct DifferentialState {
    int count=0;
    Polynomial value={-1,{}};
    std::map<int,Polynomial> images;
};

Polynomial matrix_product(Context& context,const Polynomial& left,const Polynomial& right) {
    int degree=left.degree+right.degree;
    if(degree>context.maximum)throw std::runtime_error("derivation degree overflow");
    int nl=context.masks[left.degree].size(),nr=context.masks[right.degree].size(),no=context.masks[degree].size();
    Polynomial result={degree,std::vector<I>(std::size_t(context.points)*no*matrix_entries)};
    const auto& pairs=context.table(left.degree,right.degree);
    for(int point=0;point<context.points;++point)for(int out=0;out<no;++out) {
        I accum[matrix_entries]={};
        for(const auto& pair:pairs[out]) {
            const I* lhs=left.values.data()+(std::size_t(point)*nl+pair.left)*matrix_entries;
            const I* rhs=right.values.data()+(std::size_t(point)*nr+pair.right)*matrix_entries;
            for(int row=0;row<matrix_size;++row)for(int column=0;column<matrix_size;++column)
                for(int inner=0;inner<matrix_size;++inner) accum[row*matrix_size+column]+=pair.sign*lhs[row*matrix_size+inner]*rhs[inner*matrix_size+column];
        }
        I* destination=result.values.data()+(std::size_t(point)*no+out)*matrix_entries;
        for(int j=0;j<matrix_entries;++j)destination[j]=residue(accum[j]);
    }
    return result;
}

void add_polynomial(std::map<int,Polynomial>& parts,Polynomial value,I coefficient) {
    auto [position,fresh]=parts.try_emplace(value.degree,Polynomial{value.degree,{}});
    if(fresh)position->second.values.resize(value.values.size());
    auto& target=position->second.values;
    if(target.size()!=value.values.size())throw std::runtime_error("incompatible polynomial shapes");
    coefficient=residue(coefficient);
    for(unsigned i=0;i<target.size();++i)target[i]=residue(target[i]+coefficient*value.values[i]);
}

thread_local std::string derivation_error;
extern "C" const char* bmn_derivation_error(){return derivation_error.c_str();}
extern "C" int bmn_derivation_evaluate(int k,int points,int degree,int prime,int letter_count,const int* parity,
    const std::uint64_t* matrix_offsets,const I* matrices,int rule_count,const int* rules,int word_count,
    const std::uint64_t* word_offsets,const unsigned char* encoded,I* output) {
    try {
        if(prime!=modulus||k<degree||k>14||degree>9||degree<0||letter_count>255)
            throw std::runtime_error("unsupported derivation frame");
        Context context(k,points,degree,parity,matrix_offsets,matrices);
        std::vector<std::map<int,Polynomial>> qletters(letter_count+1);
        for(int i=0;i<rule_count;++i) {
            int source=rules[4*i],coefficient=rules[4*i+1],left=rules[4*i+2],right=rules[4*i+3];
            if(source<1||source>letter_count||left<1||left>letter_count||right<1||right>letter_count)
                throw std::runtime_error("invalid rule letter");
            if(parity[left]+parity[right]>degree)continue;
            add_polynomial(qletters[source],matrix_product(context,context.first(left),context.first(right)),coefficient);
        }
        std::vector<std::vector<std::string>> words(word_count);
        std::map<std::string,DifferentialState> traces;
        for(int i=0;i<word_count;++i) {
            std::string trace;int found_degree=0;
            for(auto pos=word_offsets[i];pos<word_offsets[i+1];++pos) {
                int id=encoded[pos];
                if(id){if(id>letter_count)throw std::runtime_error("letter id");trace.push_back(char(id));found_degree+=parity[id];}
                else {if(trace.size()<2)throw std::runtime_error("short trace");words[i].push_back(trace);traces.try_emplace(trace,DifferentialState{});trace.clear();}
            }
            if(!trace.empty()||words[i].empty()||std::abs(found_degree-degree)!=1)
                throw std::runtime_error("derivation word shape or degree");
        }
        DifferentialState matrix_identity;
        matrix_identity.value={0,std::vector<I>(std::size_t(points)*matrix_entries)};
        for(int point=0;point<points;++point)for(int i=0;i<matrix_size;++i)matrix_identity.value.values[std::size_t(point)*matrix_entries+i*(matrix_size+1)]=1;
        std::string previous;std::vector<DifferentialState> stack;stack.reserve(16);
        for(auto& [word,trace]:traces) {
            unsigned common=0;while(common<std::min(word.size(),previous.size())&&word[common]==previous[common])++common;
            stack.resize(common);
            for(unsigned i=common;i<word.size();++i) {
                const auto& state=stack.empty()?matrix_identity:stack.back();
                int letter=(unsigned char)word[i],dl=parity[letter];DifferentialState following;following.count=state.count+dl;
                for(const auto& [current,image]:state.images)if(current+dl<=degree)
                    add_polynomial(following.images,context.matrix_letter(image,letter),1);
                if(state.value.degree>=0) {
                    for(const auto& [current,image]:qletters[letter])if(state.count+current<=degree)
                        add_polynomial(following.images,matrix_product(context,state.value,image),state.count%2?-1:1);
                    if(state.count+dl<=degree)following.value=context.matrix_letter(state.value,letter);
                }
                stack.push_back(std::move(following));
            }
            const auto& state=stack.back();trace.count=state.count;
            if(state.value.degree>=0)trace.value=context.trace(state.value);
            for(const auto& [current,image]:state.images)trace.images[current]=context.trace(image);
            previous=word;
        }
        DifferentialState scalar_identity;scalar_identity.value={0,std::vector<I>(points,1)};
        std::vector<std::string> previous_word;stack.clear();
        auto columns=std::size_t(points)*context.masks[degree].size();
        std::fill(output,output+std::size_t(word_count)*columns,0);
        for(int index=0;index<word_count;++index) {
            const auto& word=words[index];unsigned common=0;
            while(common<std::min(word.size(),previous_word.size())&&word[common]==previous_word[common])++common;
            stack.resize(common);
            for(unsigned i=common;i<word.size();++i) {
                const auto& state=stack.empty()?scalar_identity:stack.back();
                const auto& trace=traces.at(word[i]);DifferentialState following;following.count=state.count+trace.count;
                if(trace.value.degree>=0)for(const auto& [current,image]:state.images)if(current+trace.count<=degree)
                    add_polynomial(following.images,context.scalar_product(image,trace.value),1);
                if(state.value.degree>=0) {
                    for(const auto& [current,image]:trace.images)if(state.count+current<=degree)
                        add_polynomial(following.images,context.scalar_product(state.value,image),state.count%2?-1:1);
                    if(trace.value.degree>=0&&state.count+trace.count<=degree)
                        following.value=context.scalar_product(state.value,trace.value);
                }
                stack.push_back(std::move(following));
            }
            const auto& state=stack.back();auto found=state.images.find(degree);
            if(found!=state.images.end())std::memcpy(output+std::size_t(index)*columns,found->second.values.data(),columns*sizeof(I));
            previous_word=word;
        }
        derivation_error.clear();return 0;
    } catch(const std::exception& error){derivation_error=error.what();return -1;}
}
