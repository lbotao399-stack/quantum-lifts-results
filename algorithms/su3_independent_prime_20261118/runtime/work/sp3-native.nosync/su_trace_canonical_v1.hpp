// Graded cyclic SU trace canonicalization; no transpose quotient.
// Species/theta coefficients are supplied by q1_sp.pair_rule; this file
// performs only graded trace reconnection and finite-field accumulation.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using Word=std::vector<int>;
using Multi=std::vector<Word>;
struct Rule { int c,l,r; };
struct Pattern { int c; Multi traces; };
int prime,N;
std::vector<int> odd;
std::vector<std::vector<std::vector<Rule>>> rules;
std::unordered_map<std::string,std::vector<Pattern>> cache;
std::int64_t norm(std::int64_t x) { x%=prime; return x<0?x+prime:x; }
int parity(const Word& w) { int s=0; for(int x:w)s+=odd[x]; return s%2; }
int canonical(Multi& traces) {
    int coefficient=1;
    for(auto& w:traces) {
        if(w.size()<2)return 0;
        Word best; int bestsign=0,total=0;
        for(int x:w)total+=odd[x];
        for(int reverse=0;reverse<1;++reverse) {
            Word current=w; int sign=1;
            if(reverse) {
                std::reverse(current.begin(),current.end());
                if((w.size()+total*(total-1)/2)%2)sign=-1;
            }
            for(unsigned shift=0;shift<w.size();++shift) {
                if(best.empty()||current<best) { best=current; bestsign=sign; }
                else if(current==best && bestsign!=sign)bestsign=0;
                int first=current.front();
                if(odd[first] && (total-1)%2)sign=-sign;
                std::rotate(current.begin(),current.begin()+1,current.end());
            }
        }
        if(!bestsign)return 0;
        coefficient*=bestsign;w=std::move(best);
    }
    for(unsigned i=0;i<traces.size();++i)if(parity(traces[i]))
        for(unsigned j=i+1;j<traces.size();++j)if(parity(traces[j])) {
            if(traces[i]==traces[j])return 0;
            if(traces[i]>traces[j])coefficient=-coefficient;
        }
    std::sort(traces.begin(),traces.end());return coefficient;
}
std::string encode(const Multi& traces) {
    std::string out;out.reserve(20);
    for(const auto& w:traces) { for(int x:w)out.push_back(char(x));out.push_back(0); }
    return out;
}
