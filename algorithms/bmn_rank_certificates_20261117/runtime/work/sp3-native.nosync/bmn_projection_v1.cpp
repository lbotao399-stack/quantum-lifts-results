// BMN adapter of the pinned extremal projector; both campaign primes.
#define main q1_trace_driver_main
#include "q1_native.cpp"
#undef main
#include <array>
#include <cstring>
#include <map>
using Expression=std::unordered_map<std::string,std::int64_t>;
std::vector<std::array<int,3>> projection_letters;
std::vector<std::vector<std::pair<int,int>>> transformations;
std::string projection_error,projection_bytes;
std::vector<std::uint64_t> projection_offsets;
std::vector<std::int64_t> projection_coefficients;
Multi decode_projection_word(const std::string& key) {
    Multi result;Word trace;
    for(unsigned char id:key){if(id)trace.push_back(id);else{result.push_back(trace);trace.clear();}}
    if(!trace.empty())throw std::runtime_error("unterminated trace");return result;
}
Expression transformed(const Expression& input,int operation) {
    Expression result;
    for(const auto& [key,c]:input) {
        Multi word=decode_projection_word(key);
        for(unsigned i=0;i<word.size();++i)for(unsigned j=0;j<word[i].size();++j) {
            auto [replacement,factor]=transformations[operation][word[i][j]];if(!replacement)continue;
            Multi changed=word;changed[i][j]=replacement;int sign=canonical(changed);if(!sign)continue;
            auto target=encode(changed);auto& value=result[target];value=norm(value+c*factor*sign);
            if(!value)result.erase(target);
        }
    }
    return result;
}
std::int64_t inverse_projection(std::int64_t value) {
    value=norm(value);if(!value)throw std::runtime_error("zero projector denominator");
    std::int64_t result=1;int power=prime-2;
    while(power){if(power&1)result=result*value%prime;value=value*value%prime;power>>=1;}
    return result;
}
Expression root_project(const Expression& input,int raising,int lowering,int weight,int shift) {
    Expression raised=input,result=input;std::int64_t coefficient=1;
    for(int power=1;power<200;++power) {
        raised=transformed(raised,raising);if(raised.empty())return result;
        coefficient=norm(-coefficient*inverse_projection(power*(weight+shift+power)));
        Expression lowered=raised;
        for(int repeat=0;repeat<power;++repeat)lowered=transformed(lowered,lowering);
        for(const auto& [key,c]:lowered){auto& value=result[key];value=norm(value+coefficient*c);if(!value)result.erase(key);}
    }
    throw std::runtime_error("nonterminating projector");
}
extern "C" const char* full_projection_error(){return projection_error.c_str();}
extern "C" std::uint64_t full_projection_size(){return projection_bytes.size();}
extern "C" void full_projection_copy(unsigned char* bytes,std::uint64_t* offsets,std::int64_t* coefficients) {
    std::memcpy(bytes,projection_bytes.data(),projection_bytes.size());
    std::memcpy(offsets,projection_offsets.data(),projection_offsets.size()*sizeof(std::uint64_t));
    std::memcpy(coefficients,projection_coefficients.data(),projection_coefficients.size()*sizeof(std::int64_t));
}
extern "C" int full_projection(int p,int alphabet,const int* letter_data,const int* sector,int word_count,
    const std::uint64_t* offsets,const unsigned char* bytes,const std::int64_t* coefficients) {
    try {
        if((p!=1000003&&p!=1000033)||alphabet>255)throw std::runtime_error("unsupported projection input");prime=p;
        odd.resize(alphabet+1);projection_letters.resize(alphabet+1);std::map<std::array<int,3>,int> ids;
        for(int i=1;i<=alphabet;++i){for(int j=0;j<3;++j)projection_letters[i][j]=letter_data[(i-1)*3+j];ids[projection_letters[i]]=i;odd[i]=(__builtin_popcount(projection_letters[i][2])+1)%2;}
        transformations.assign(8,std::vector<std::pair<int,int>>(alphabet+1,{0,0}));
        int rows[6]={1,2,0,2,0,1},columns[6]={2,1,2,0,1,0};
        for(int operation=0;operation<8;++operation)for(int id=1;id<=alphabet;++id) {
            auto L=projection_letters[id];int factor=1;
            if(operation<6) {
                int row=rows[operation],column=columns[operation];
                if(!(L[2]&(1<<column))||(L[2]&(1<<row)))continue;
                int exponent=0;for(int j=std::min(row,column)+1;j<std::max(row,column);++j)exponent+=bool(L[2]&(1<<j));
                factor=exponent%2?-1:1;L[2]^=(1<<row)|(1<<column);
            } else if(operation==6){if(!L[1])continue;factor=L[0]+1;++L[0];--L[1];}
            else {if(!L[0])continue;factor=L[1]+1;--L[0];++L[1];}
            auto found=ids.find(L);if(found==ids.end())throw std::runtime_error("alphabet not closed under projector");
            transformations[operation][id]={found->second,factor};
        }
        Expression expr;
        for(int i=0;i<word_count;++i) {
            std::string key(reinterpret_cast<const char*>(bytes+offsets[i]),offsets[i+1]-offsets[i]);
            auto& value=expr[key];value=norm(value+coefficients[i]);if(!value)expr.erase(key);
        }
        for(int i=0;i<3;++i) {
            int row=rows[2*i],column=columns[2*i],difference=sector[row+2]-sector[column+2];
            if(difference<0||difference%2)throw std::runtime_error("not a dominant weight");
            expr=root_project(expr,2*i,2*i+1,difference/2,column-row);
        }
        int difference=sector[0]-sector[1];if(difference<0||difference%2)throw std::runtime_error("not a dominant spin");
        expr=root_project(expr,6,7,difference/2,1);
        for(int operation:{0,4,6})if(!transformed(expr,operation).empty())throw std::runtime_error("highest-weight check failed");
        std::map<std::string,std::int64_t> sorted(expr.begin(),expr.end());
        projection_bytes.clear();projection_offsets.clear();projection_coefficients.clear();projection_offsets.push_back(0);
        for(const auto& [key,c]:sorted){projection_bytes+=key;projection_offsets.push_back(projection_bytes.size());projection_coefficients.push_back(c);}
        projection_error.clear();return projection_coefficients.size();
    } catch(const std::exception& error){projection_error=error.what();return -1;}
}
