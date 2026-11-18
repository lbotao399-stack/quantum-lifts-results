// Integer 2 Q1 for the pinned symplectic trace convention.
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
        for(int reverse=0;reverse<2;++reverse) {
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
const std::vector<Pattern>& patterns(const Multi& input,const Word& par,int i,int j,int leftodd,int rightodd) {
    std::string key;key.push_back(char(input.size()));
    for(const auto& w:input)key.push_back(char(w.size()));
    for(int x:par)key.push_back(char(x));
    for(int x:{i,j,leftodd,rightodd})key.push_back(char(x));
    auto found=cache.find(key);if(found!=cache.end())return found->second;
    Multi ids;std::vector<std::pair<int,int>> location;int k=0;
    for(unsigned t=0;t<input.size();++t) {
        Word w;for(unsigned u=0;u<input[t].size();++u){w.push_back(k++);location.push_back({t,u});}
        ids.push_back(std::move(w));
    }
    Word newpar=par;newpar[i]=leftodd;newpar[j]=rightodd;
    Word position(k);position[i]=0;position[j]=1;int next=2;
    for(int x=0;x<k;++x)if(x!=i && x!=j)position[x]=next++;
    int front=0;
    for(int x=0;x<i;++x)front+=par[i]*par[x];
    for(int x=0;x<j;++x)if(x!=i)front+=par[j]*par[x];
    std::vector<Pattern> result;
    auto add=[&](Multi traces,int c) {
        Word written;
        for(const auto& w:traces) {
            if(w.size()==1)return;
            written.insert(written.end(),w.begin(),w.end());
        }
        int exponent=front;
        for(unsigned x=0;x<written.size();++x)
            for(unsigned y=x+1;y<written.size();++y)
                if(position[written[x]]>position[written[y]])exponent+=newpar[written[x]]*newpar[written[y]];
        if(exponent%2)c=-c;
        result.push_back({c,std::move(traces)});
    };
    auto slice=[](const Word& w,int start,int end){return Word(w.begin()+start,w.begin()+end);};
    auto join=[](std::initializer_list<Word> parts){Word out;for(const auto& w:parts)out.insert(out.end(),w.begin(),w.end());return out;};
    auto [p,u]=location[i];auto [q,v]=location[j];
    for(int cx=0;cx<2;++cx)for(int cy=0;cy<2;++cy) {
        int sign=(cx+cy)%2?-1:1,g1=u+cx,g2=v+cy;
        Multi rest;
        for(unsigned t=0;t<ids.size();++t)if(int(t)!=p && int(t)!=q)rest.push_back(ids[t]);
        if(p!=q) {
            auto A=slice(ids[p],0,g1),B=slice(ids[p],g1,ids[p].size());
            auto D=slice(ids[q],0,g2),E=slice(ids[q],g2,ids[q].size());
            auto first=rest;first.push_back(join({A,E,D,B}));add(std::move(first),sign);
            std::reverse(D.begin(),D.end());std::reverse(E.begin(),E.end());
            rest.push_back(join({A,D,E,B}));add(std::move(rest),-sign*((D.size()+E.size())%2?-1:1));
        } else {
            int lo=std::min(g1,g2),hi=std::max(g1,g2);
            auto A=slice(ids[p],0,lo),B=slice(ids[p],lo,hi),C=slice(ids[p],hi,ids[p].size());
            auto first=rest;int factor=1;
            if(B.empty())factor*=N;else first.push_back(B);
            auto AC=join({A,C});if(AC.empty())factor*=N;else first.push_back(AC);
            add(std::move(first),sign*factor);
            std::reverse(B.begin(),B.end());auto second=join({A,B,C});
            if(second.empty())add(rest,sign*(B.size()%2?-1:1)*N);
            else {rest.push_back(second);add(rest,sign*(B.size()%2?-1:1));}
        }
    }
    return cache.emplace(std::move(key),std::move(result)).first->second;
}
void action(const Multi& input,std::int64_t coefficient,std::unordered_map<std::string,std::int64_t>& output) {
    Word letters,par;
    for(const auto& w:input)for(int x:w){letters.push_back(x);par.push_back(odd[x]);}
    for(unsigned i=0;i<letters.size();++i)for(unsigned j=i+1;j<letters.size();++j)
        for(const auto& rule:rules[letters[i]][letters[j]])
            for(const auto& pattern:patterns(input,par,i,j,odd[rule.l],odd[rule.r])) {
                Multi target=pattern.traces;
                for(auto& w:target)for(int& x:w)x=x==int(i)?rule.l:x==int(j)?rule.r:letters[x];
                int sign=canonical(target);if(!sign)continue;
                auto key=encode(target);auto& value=output[key];
                value=norm(value+coefficient*rule.c*pattern.c*sign);
                if(!value)output.erase(key);
            }
}
template<class T> void put(std::ofstream& out,const T& value){out.write(reinterpret_cast<const char*>(&value),sizeof(T));}
template<class T> void get(std::ifstream& in,T& value){in.read(reinterpret_cast<char*>(&value),sizeof(T));}
void save(const std::string& path,std::uint64_t cursor,const std::unordered_map<std::string,std::int64_t>& output) {
    std::ofstream out(path+".tmp",std::ios::binary);put(out,cursor);std::uint64_t count=output.size();put(out,count);
    for(const auto& [word,value]:output) { std::uint32_t length=word.size(),c=value;put(out,length);out.write(word.data(),length);put(out,c); }
    out.close();if(!out)throw std::runtime_error("checkpoint write failed");
    if(std::rename((path+".tmp").c_str(),path.c_str()))throw std::runtime_error("checkpoint replace failed");
}
int main(int argc,char** argv) {
    if(argc<3)throw std::runtime_error("q1_native input checkpoint [seconds]");
    double limit=argc>3?std::stod(argv[3]):3600;
    std::ifstream in(argv[1]);int alphabet,rulecount;std::uint64_t total;
    in>>prime>>N>>alphabet;odd.resize(alphabet+1);
    for(int i=1;i<=alphabet;++i)in>>odd[i];
    rules.resize(alphabet+1,std::vector<std::vector<Rule>>(alphabet+1));
    in>>rulecount;
    for(int i=0;i<rulecount;++i){int left,right,c,l,r;in>>left>>right>>c>>l>>r;rules[left][right].push_back({c,l,r});}
    in>>total;std::uint64_t cursor=0;
    std::unordered_map<std::string,std::int64_t> output;
    std::ifstream old(argv[2],std::ios::binary);
    if(old) {
        std::uint64_t count;get(old,cursor);get(old,count);
        for(std::uint64_t i=0;i<count;++i){std::uint32_t length,c;get(old,length);std::string key(length,'\0');old.read(key.data(),length);get(old,c);output[key]=c;}
        if(!old)throw std::runtime_error("truncated checkpoint");
    }
    auto start=std::chrono::steady_clock::now(),last=start;
    for(std::uint64_t index=0;index<total;++index) {
        std::int64_t c;int number;in>>c>>number;Multi word(number);
        for(auto& w:word){int length;in>>length;w.resize(length);for(int& x:w)in>>x;}
        if(!in)throw std::runtime_error("truncated input");
        if(index<cursor)continue;
        action(word,c,output);
        if((index+1)%1024==0 || index+1==total) {
            auto now=std::chrono::steady_clock::now();double elapsed=std::chrono::duration<double>(now-start).count();
            if(std::chrono::duration<double>(now-last).count()>30 || elapsed>limit || index+1==total) {
                save(argv[2],index+1,output);
                std::cout<<"{\"processed\":"<<index+1<<",\"total\":"<<total<<",\"output_terms\":"<<output.size()<<",\"seconds\":"<<elapsed<<",\"cache_patterns\":"<<cache.size()<<"}"<<std::endl;
                last=now;
            }
            if(elapsed>limit)break;
        }
    }
}
