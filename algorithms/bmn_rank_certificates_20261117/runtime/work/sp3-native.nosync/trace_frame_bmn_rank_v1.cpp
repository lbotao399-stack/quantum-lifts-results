// Exact modular matrix/exterior evaluation in the pinned subset ordering.
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using I=std::int64_t;
#ifndef TRACE_FRAME_PRIME
#define TRACE_FRAME_PRIME 1000003
#endif
#ifndef TRACE_FRAME_MATRIX_SIZE
#define TRACE_FRAME_MATRIX_SIZE 6
#endif
constexpr int matrix_size=TRACE_FRAME_MATRIX_SIZE;
constexpr int matrix_entries=matrix_size*matrix_size;
constexpr I modulus=TRACE_FRAME_PRIME;
I residue(I x){x%=modulus;return x<0?x+modulus:x;}
struct Pair {int left,right,sign;};
struct Polynomial {int degree;std::vector<I> values;};
struct Context {
    int k,points,maximum;
    const int* parity;const std::uint64_t* offsets;const I* letters;
    std::vector<std::vector<unsigned>> masks;
    std::vector<std::unordered_map<unsigned,int>> indices;
    std::map<std::pair<int,int>,std::vector<std::vector<Pair>>> tables;
    void combinations(int next,int remaining,unsigned mask,std::vector<unsigned>& result) {
        if(!remaining){result.push_back(mask);return;}
        for(int i=next;i<=k-remaining;++i)combinations(i+1,remaining-1,mask|(1u<<i),result);
    }
    Context(int k0,int p0,int d0,const int* par,const std::uint64_t* off,const I* data):k(k0),points(p0),maximum(d0),parity(par),offsets(off),letters(data) {
        masks.resize(maximum+1);indices.resize(maximum+1);
        for(int d=0;d<=maximum;++d){combinations(0,d,0,masks[d]);for(unsigned i=0;i<masks[d].size();++i)indices[d][masks[d][i]]=i;}
    }
    const std::vector<std::vector<Pair>>& table(int left,int right) {
        auto key=std::make_pair(left,right);auto old=tables.find(key);if(old!=tables.end())return old->second;
        int total=left+right;if(total>maximum)throw std::runtime_error("exterior degree overflow");
        std::vector<std::vector<Pair>> result(masks[total].size());
        for(unsigned out=0;out<masks[total].size();++out) {
            unsigned mask=masks[total][out],sub=mask;
            while(true) {
                if(__builtin_popcount(sub)==left) {
                    unsigned complement=mask^sub;int inversions=0;
                    for(int j=0;j<k;++j)if(complement&(1u<<j))inversions+=__builtin_popcount(sub>>(j+1));
                    result[out].push_back({indices[left].at(sub),indices[right].at(complement),inversions%2?-1:1});
                }
                if(!sub)break;sub=(sub-1)&mask;
            }
        }
        return tables.emplace(key,std::move(result)).first->second;
    }
    Polynomial first(int letter) {
        int degree=parity[letter],count=masks[degree].size();const I* begin=letters+offsets[letter];
        return {degree,std::vector<I>(begin,begin+std::size_t(points)*count*matrix_entries)};
    }
    Polynomial matrix_letter(const Polynomial& left,int letter) {
        int dl=left.degree,dr=parity[letter],degree=dl+dr;
        if(degree>maximum)throw std::runtime_error("trace degree overflow");
        int nleft=masks[dl].size(),nright=masks[dr].size(),nout=masks[degree].size();
        const auto& pairs=table(dl,dr);const I* right=letters+offsets[letter];
        Polynomial result={degree,std::vector<I>(std::size_t(points)*nout*matrix_entries)};
        for(int point=0;point<points;++point)for(int out=0;out<nout;++out) {
            I accum[matrix_entries]={};
            for(const auto& pair:pairs[out]) {
                const I* a=left.values.data()+(std::size_t(point)*nleft+pair.left)*matrix_entries;
                const I* b=right+(std::size_t(point)*nright+pair.right)*matrix_entries;
                for(int row=0;row<matrix_size;++row)for(int column=0;column<matrix_size;++column)
                    for(int inner=0;inner<matrix_size;++inner) accum[row*matrix_size+column]+=pair.sign*a[row*matrix_size+inner]*b[inner*matrix_size+column];
            }
            I* destination=result.values.data()+(std::size_t(point)*nout+out)*matrix_entries;
            for(int j=0;j<matrix_entries;++j)destination[j]=residue(accum[j]);
        }
        return result;
    }
    Polynomial trace(const Polynomial& matrix) {
        Polynomial result={matrix.degree,std::vector<I>(std::size_t(points)*masks[matrix.degree].size())};
        for(unsigned i=0;i<result.values.size();++i){const I* a=matrix.values.data()+std::size_t(i)*matrix_entries;I value=0;for(int j=0;j<matrix_size;++j)value+=a[j*(matrix_size+1)];result.values[i]=residue(value);}
        return result;
    }
    Polynomial scalar_product(const Polynomial& left,const Polynomial& right) {
        int degree=left.degree+right.degree;if(degree>maximum)throw std::runtime_error("product degree overflow");
        int nl=masks[left.degree].size(),nr=masks[right.degree].size(),no=masks[degree].size();
        Polynomial result={degree,std::vector<I>(std::size_t(points)*no)};const auto& pairs=table(left.degree,right.degree);
        for(int p=0;p<points;++p)for(int out=0;out<no;++out) {
            I value=0;
            for(const auto& pair:pairs[out])value+=pair.sign*left.values[std::size_t(p)*nl+pair.left]*right.values[std::size_t(p)*nr+pair.right];
            result.values[std::size_t(p)*no+out]=residue(value);
        }
        return result;
    }
};
std::string last_error;
extern "C" const char* trace_frame_error(){return last_error.c_str();}
extern "C" int trace_frame_evaluate(int k,int points,int degree,int prime,int letter_count,const int* parity,
    const std::uint64_t* matrix_offsets,const I* matrices,int word_count,const std::uint64_t* word_offsets,
    const unsigned char* encoded,I* output) {
    try {
        if(prime!=modulus||k<degree||k>14||degree>9||letter_count>255)throw std::runtime_error("unsupported frame");
        Context context(k,points,degree,parity,matrix_offsets,matrices);
        std::vector<std::vector<std::string>> words(word_count);
        std::map<std::string,Polynomial> traces;
        for(int i=0;i<word_count;++i) {
            std::string trace;int found_degree=0;
            for(auto pos=word_offsets[i];pos<word_offsets[i+1];++pos) {
                int id=encoded[pos];
                if(id){if(id>letter_count)throw std::runtime_error("letter id");trace.push_back(char(id));found_degree+=parity[id];}
                else {if(trace.size()<2)throw std::runtime_error("short trace");words[i].push_back(trace);traces.try_emplace(trace,Polynomial{});trace.clear();}
            }
            if(!trace.empty()||words[i].empty()||found_degree!=degree)throw std::runtime_error("word shape or degree");
        }
        std::string previous;std::vector<Polynomial> matrix_stack;matrix_stack.reserve(16);
        for(auto& [word,value]:traces) {
            unsigned common=0;while(common<std::min(word.size(),previous.size())&&word[common]==previous[common])++common;
            matrix_stack.resize(common);
            if(matrix_stack.empty()){matrix_stack.push_back(context.first((unsigned char)word[0]));common=1;}
            for(unsigned i=common;i<word.size();++i)matrix_stack.push_back(context.matrix_letter(matrix_stack.back(),(unsigned char)word[i]));
            value=context.trace(matrix_stack.back());previous=word;
        }
        std::vector<std::string> previous_word;std::vector<Polynomial> scalar_stack;scalar_stack.reserve(16);
        auto columns=std::size_t(points)*context.masks[degree].size();
        for(int index=0;index<word_count;++index) {
            const auto& word=words[index];unsigned common=0;
            while(common<std::min(word.size(),previous_word.size())&&word[common]==previous_word[common])++common;
            scalar_stack.resize(common);
            if(scalar_stack.empty()){scalar_stack.push_back(traces.at(word[0]));common=1;}
            for(unsigned i=common;i<word.size();++i)scalar_stack.push_back(context.scalar_product(scalar_stack.back(),traces.at(word[i])));
            const auto& result=scalar_stack.back();if(result.degree!=degree||result.values.size()!=columns)throw std::runtime_error("output shape");
            std::memcpy(output+std::size_t(index)*columns,result.values.data(),columns*sizeof(I));previous_word=word;
        }
        last_error.clear();return 0;
    } catch(const std::exception& error){last_error=error.what();return -1;}
}
