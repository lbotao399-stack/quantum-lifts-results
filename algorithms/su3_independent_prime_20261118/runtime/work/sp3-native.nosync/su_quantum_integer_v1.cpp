// Exact signed-integer SU trace reconnection. Rational species rules are input.
// The U(1) completeness subtraction cancels inside the commutators.
#include "su_trace_canonical_v1.hpp"
#include <array>
#include <cstring>
#include <limits>
#include <map>

struct IntegerRule { std::int64_t coefficient; int first,second; };
static std::unordered_map<std::uint32_t,std::vector<IntegerRule>> species[3];
static std::unordered_map<std::string,std::vector<Pattern>> reconnection_cache;
static std::unordered_map<std::string,std::int64_t> output;
static std::vector<unsigned char> exported;
static std::string error;

static std::int64_t checked(__int128 value) {
    if(value<std::numeric_limits<std::int64_t>::min() || value>std::numeric_limits<std::int64_t>::max())
        throw std::runtime_error("exact integer accumulation overflow");
    return std::int64_t(value);
}

static std::uint32_t rule_key(int left,int middle,int right) {
    return (std::uint32_t(left)<<16)|(std::uint32_t(middle)<<8)|std::uint32_t(right);
}

static Multi contract(const Multi& input,int token) {
    std::vector<std::pair<int,int>> positions;
    for(unsigned t=0;t<input.size();++t)for(unsigned u=0;u<input[t].size();++u)
        if(input[t][u]==token)positions.push_back({t,u});
    if(positions.size()!=2)throw std::runtime_error("generator multiplicity");
    auto [p,u]=positions[0];auto [q,v]=positions[1];
    Multi result;
    for(unsigned t=0;t<input.size();++t)if(int(t)!=p && int(t)!=q)result.push_back(input[t]);
    if(p==q) {
        result.push_back(Word(input[p].begin()+u+1,input[p].begin()+v));
        Word right(input[p].begin()+v+1,input[p].end());
        right.insert(right.end(),input[p].begin(),input[p].begin()+u);
        result.push_back(std::move(right));
    } else {
        Word left(input[p].begin()+u+1,input[p].end());
        left.insert(left.end(),input[p].begin(),input[p].begin()+u);
        left.insert(left.end(),input[q].begin()+v+1,input[q].end());
        left.insert(left.end(),input[q].begin(),input[q].begin()+v);
        result.push_back(std::move(left));
    }
    return result;
}

static const std::vector<Pattern>& reconnect(const Multi& input,const Word& par,
        int order,int left,int middle,int right,int firstodd,int secondodd) {
    std::string key;
    for(int x:{order,int(input.size())})key.push_back(char(x));
    for(const auto& word:input)key.push_back(char(word.size()));
    for(int x:par)key.push_back(char(x));
    for(int x:{left,middle,right,firstodd,secondodd})key.push_back(char(x));
    auto found=reconnection_cache.find(key);
    if(found!=reconnection_cache.end())return found->second;
    Multi ids;int size=0;
    for(const auto& word:input) {
        Word row;for(unsigned j=0;j<word.size();++j)row.push_back(size++);
        ids.push_back(std::move(row));
    }
    Word moved{left};if(order==2)moved.push_back(middle);moved.push_back(right);
    for(int i=0;i<size;++i)if(i!=left && i!=middle && i!=right)moved.push_back(i);
    int front=0;
    for(int i=0;i<size;++i)for(int j=i+1;j<size;++j)
        if(moved[i]>moved[j])front+=par[moved[i]]*par[moved[j]];
    Word position(size,-1);position[left]=0;position[right]=1;int index=2;
    for(int i=0;i<size;++i)if(i!=left && i!=middle && i!=right)position[i]=index++;
    Word newpar=par;newpar[left]=firstodd;newpar[right]=secondodd;
    std::map<Multi,int> collected;
    for(int ca=0;ca<2;++ca)for(int cm=0;cm<(order==2?2:1);++cm)for(int cc=0;cc<2;++cc) {
        Multi inserted;
        for(const auto& word:ids) {
            Word row;
            for(int t:word) {
                if(t==left){row.push_back(ca?t:-1);row.push_back(ca?-1:t);}
                else if(order==2 && t==middle){row.push_back(cm?-2:-1);row.push_back(cm?-1:-2);}
                else if(t==right){int token=order==2?-2:-1;row.push_back(cc?t:token);row.push_back(cc?token:t);}
                else row.push_back(t);
            }
            inserted.push_back(std::move(row));
        }
        Multi contracted=contract(inserted,-1);
        if(order==2)contracted=contract(contracted,-2);
        Multi target;Word written;int factor=1;bool zero=false;
        for(const auto& word:contracted) {
            if(word.size()==1){zero=true;break;}
            if(word.empty())factor*=N;
            else {target.push_back(word);written.insert(written.end(),word.begin(),word.end());}
        }
        if(zero)continue;
        if(written.size()!=unsigned(size-(order==2)))throw std::runtime_error("output occurrence count");
        int exponent=front+ca+cm+cc;
        for(unsigned i=0;i<written.size();++i)for(unsigned j=i+1;j<written.size();++j)
            if(position[written[i]]>position[written[j]])exponent+=newpar[written[i]]*newpar[written[j]];
        if(exponent%2)factor=-factor;
        collected[target]+=factor;
    }
    std::vector<Pattern> result;
    for(auto& item:collected)if(item.second)result.push_back({item.second,item.first});
    if(reconnection_cache.size()>8192)reconnection_cache.clear();
    return reconnection_cache.emplace(std::move(key),std::move(result)).first->second;
}

extern "C" const char* su_integer_error(){return error.c_str();}
extern "C" int su_integer_init(int size,int alphabet,const int* parities) {
    try {
        if(size<2 || size>16 || alphabet<1 || alphabet>255)throw std::runtime_error("configuration domain");
        N=size;odd.assign(parities,parities+alphabet+1);
        for(auto& rules:species)rules.clear();reconnection_cache.clear();output.clear();exported.clear();error.clear();
        return 0;
    }catch(const std::exception& e){error=e.what();return -1;}
}
extern "C" int su_integer_rule(int order,int left,int middle,int right,std::int64_t coefficient,int first,int second) {
    try {
        if(order<1 || order>2)throw std::runtime_error("invalid order");
        for(int i:{left,right,first,second})if(i<1 || i>=int(odd.size()))throw std::runtime_error("invalid letter");
        if(order==2 && (middle<1 || middle>=int(odd.size())))throw std::runtime_error("invalid middle letter");
        species[order][rule_key(left,middle,right)].push_back({coefficient,first,second});return 0;
    }catch(const std::exception& e){error=e.what();return -1;}
}
extern "C" int su_integer_apply(int order,const unsigned char* encoded,int length,std::int64_t coefficient) {
    try {
        Multi input;Word row;
        for(int j=0;j<length;++j) {
            int i=encoded[j];
            if(!i){if(row.size()<2)throw std::runtime_error("short trace");input.push_back(row);row.clear();}
            else {if(i>=int(odd.size()))throw std::runtime_error("unknown letter");row.push_back(i);}
        }
        if(!row.empty() || input.empty())throw std::runtime_error("invalid word");
        Word letters,par;for(const auto& trace:input)for(int i:trace){letters.push_back(i);par.push_back(odd[i]);}
        auto apply=[&](int left,int middle,int right) {
            auto selected=species[order].find(rule_key(letters[left],order==2?letters[middle]:0,letters[right]));
            if(selected==species[order].end())return;
            for(const auto& rule:selected->second)
                for(const auto& pattern:reconnect(input,par,order,left,middle,right,odd[rule.first],odd[rule.second])) {
                    Multi target=pattern.traces;
                    for(auto& trace:target)for(int& i:trace)i=i==left?rule.first:i==right?rule.second:letters[i];
                    int sign=canonical(target);if(!sign)continue;
                    std::string key=encode(target);
                    auto& value=output[key];
                    value=checked(__int128(value)+__int128(coefficient)*rule.coefficient*pattern.c*sign);
                    if(!value)output.erase(key);
                }
        };
        if(order==1)for(unsigned i=0;i<letters.size();++i)for(unsigned j=i+1;j<letters.size();++j)apply(i,-1,j);
        else if(order==2)for(unsigned i=0;i<letters.size();++i)for(unsigned j=i+1;j<letters.size();++j)
            for(unsigned k=j+1;k<letters.size();++k)for(int middle:{int(i),int(j),int(k)})
                apply(middle==int(i)?j:i,middle,middle==int(k)?j:k);
        else throw std::runtime_error("invalid order");
        return 0;
    }catch(const std::exception& e){error=e.what();return -1;}
}
extern "C" const unsigned char* su_integer_finish(std::uint64_t* size) {
    exported.clear();
    for(const auto& item:output) {
        std::uint32_t length=item.first.size();std::uint64_t value=std::uint64_t(item.second);
        for(int shift=0;shift<32;shift+=8)exported.push_back((length>>shift)&255);
        for(int shift=0;shift<64;shift+=8)exported.push_back((value>>shift)&255);
        exported.insert(exported.end(),item.first.begin(),item.first.end());
    }
    *size=exported.size();return exported.data();
}
