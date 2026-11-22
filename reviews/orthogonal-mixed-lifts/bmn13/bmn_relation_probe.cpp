// Finite orthogonal trace relations and their Q0 images in a three-scalar highest weight.
// Compile with c++ -std=c++17 -O3. All printed entries lift from rationals regular at the prime.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>
using std::array; using std::map; using std::pair; using std::string; using std::vector;
int prime=1000003;
int mod(int64_t n){int r=n%prime;return r<0?r+prime:r;}
int power(int n,int k){int r=1;for(;k;k/=2,n=mod((int64_t)n*n))if(k%2)r=mod((int64_t)r*n);return r;}
int add(int a,int b){int c=a+b;return c>=prime?c-prime:c;}
int mul(int a,int b){return mod((int64_t)a*b);}
struct Mat{int a=0,b=0,c=0,d=0;};
Mat plus(Mat x,Mat y){return {add(x.a,y.a),add(x.b,y.b),add(x.c,y.c),add(x.d,y.d)};}
Mat times(Mat x,Mat y){return {mod((int64_t)x.a*y.a+(int64_t)x.b*y.c),mod((int64_t)x.a*y.b+(int64_t)x.b*y.d),
                              mod((int64_t)x.c*y.a+(int64_t)x.d*y.c),mod((int64_t)x.c*y.b+(int64_t)x.d*y.d)};}
Mat scale(Mat x,int c){return {mul(x.a,c),mul(x.b,c),mul(x.c,c),mul(x.d,c)};}
Mat adj(Mat x){return {x.d,mod(-x.b),mod(-x.c),x.a};}
using Poly=array<int,16>;
struct RNG{
    std::mt19937_64 gen;
    explicit RNG(uint64_t s):gen(s){}
    int pick(int bound){return gen()%bound;}
    int entry(){return mod(pick(9)-4);}
    template<class T> void shuffle(vector<T>& v){for(int i=int(v.size())-1;i>0;i--)std::swap(v[i],v[pick(i+1)]);}
};
struct ProbeBank{
    int total,epsilons,merged=-1;
    bool onePsi;
    map<string,pair<int,int>> canonical;
    map<string,int> traceIndex;
    vector<string> traces;
    map<vector<int>,int> termIndex;
    vector<vector<int>> terms;
    vector<vector<pair<int,int>>> probes;
    vector<vector<int>> shapes;
    explicit ProbeBank(int boxes,int eps,bool psi=false):total(boxes),epsilons(eps),onePsi(psi){}
    pair<int,int> trace(const string& word){
        auto found=canonical.find(word);
        if(found!=canonical.end())return found->second;
        const int len=word.size(),rev=len%2?-1:1;
        map<string,int> orbit;
        for(string t:{word,string(word.rbegin(),word.rend())}){
            int s=t==word?1:rev;
            if(t==word && t==string(word.rbegin(),word.rend()) && rev<0)s=0;
            for(int i=0;i<len;i++){
                string rot=t.substr(i)+t.substr(0,i);
                if(orbit.count(rot)&&orbit[rot]!=s)orbit[rot]=0;
                else if(!orbit.count(rot))orbit[rot]=s;
            }
        }
        if(std::any_of(orbit.begin(),orbit.end(),[](const auto& kv){return kv.second==0;})){
            for(const auto& kv:orbit)canonical[kv.first]={0,-1};
            return {0,-1};
        }
        string standard=orbit.begin()->first;
        int idx;
        if(traceIndex.count(standard))idx=traceIndex.at(standard);
        else{idx=traces.size();traceIndex[standard]=idx;traces.push_back(standard);}
        int s=orbit.begin()->second;
        for(const auto& kv:orbit)canonical[kv.first]={s*kv.second,idx};
        return canonical.at(word);
    }
    void collect(const vector<int>& shape,const vector<int>& slots,int block,int sign,string& word,map<vector<int>,int>& result){
        if(block==epsilons){
            int pos=0;
            string fused;
            if(onePsi){
                int i=word[merged]-'0',j=word[merged+1]-'0';
                if(i==j)return;
                if(i>j){std::swap(i,j);sign=-sign;}
                fused=word.substr(0,merged)+char('2'+i+j)+word.substr(merged+2);
            }
            const string& written=onePsi?fused:word;
            vector<int> ids;
            for(int len:shape){
                auto [s,id]=trace(written.substr(pos,len));pos+=len;
                if(!s)return;
                sign*=s;ids.push_back(id);
            }
            std::sort(ids.begin(),ids.end());
            result[ids]+=sign;return;
        }
        array<char,3> perm={'0','1','2'};
        do{
            int neg=(perm[0]>perm[1])+(perm[0]>perm[2])+(perm[1]>perm[2]);
            for(int i=0;i<3;i++)word[slots[3*block+i]]=perm[i];
            collect(shape,slots,block+1,neg%2?-sign:sign,word,result);
        }while(std::next_permutation(perm.begin(),perm.end()));
    }
    void partitions(int left,int low,vector<int>& v,vector<vector<int>>& all){
        if(left==0){if(v.size()>=2)all.push_back(v);return;}
        for(int i=low;i<=std::min(left,8);i++){
            v.push_back(i);partitions(left-i,i,v,all);v.pop_back();
        }
    }
    void build(int count,uint64_t seed){
        RNG rng(seed);vector<vector<int>> all;vector<int> tmp;partitions(total-int(onePsi),2,tmp,all);
        // Short factors expose different coordinates from long trace words.
        // Repeat each partition in random order. Zero projected probes are discarded.
        int ntry=0;std::set<map<vector<int>,int>> seen;
        while(int(probes.size())<count){
            if(ntry%all.size()==0)rng.shuffle(all);
            auto shape=all[ntry%all.size()];ntry++;
            if(onePsi)merged=rng.pick(total-1);
            vector<int> slots(total);std::iota(slots.begin(),slots.end(),0);rng.shuffle(slots);
            string word(total,'0');map<vector<int>,int> coeffs;
            collect(shape,slots,0,1,word,coeffs);
            for(auto it=coeffs.begin();it!=coeffs.end();)if(!it->second)it=coeffs.erase(it);else ++it;
            if(coeffs.empty()||seen.count(coeffs))continue;
            // Avoid keeping a second copy of large probes merely for de-duplication.
            if(coeffs.size()<40)seen.insert(coeffs);
            vector<pair<int,int>> probe;
            for(auto& [key,val]:coeffs){
                auto [it,fresh]=termIndex.emplace(key,terms.size());
                if(fresh)terms.push_back(key);
                probe.emplace_back(it->second,mod(val));
            }
            shapes.push_back(shape);probes.push_back(std::move(probe));
        }
    }
};
struct Recipe{
    int marks=0,full=0,maxlen=1;
    map<string,vector<pair<int,Mat>>> tiles;
    map<string,int> closed;
    int closedBit=0;
    bool varyScalars=false;
    array<int,3> variationAxis{};
    // Differentiating a closed trace on labelled sites counts every cyclic cut and
    // its reverse. This includes all stabilizers for repeated letters.
    void addClosed(const string& w,int coeff){
        string back(w.rbegin(),w.rend());
        for(size_t i=0;i<w.size();i++){
            string a=w.substr(i)+w.substr(0,i),b=back.substr(i)+back.substr(0,i);
            closed[a]=add(closed[a],coeff);
            closed[b]=add(closed[b],w.size()%2?mod(-coeff):coeff);
        }
    }
    // A marked input path is an ordered product of adjoint matrices. P is replaced by a commutator.
    void addPath(const string& pattern,int bit,RNG& rng){
        map<string,int> words={{"",1}};
        for(char letter:pattern){
            map<string,int> next,local;
            if(letter=='P'){
                for(int i=0;i<3;i++)for(int j=i+1;j<3;j++){
                    int c=rng.entry();local[string{char('0'+i),char('0'+j)}]=c;
                    local[string{char('0'+j),char('0'+i)}]=mod(-c);
                }
            }else for(int i=0;i<3;i++)local[string(1,char('0'+i))]=rng.entry();
            for(auto& [w,a]:words)for(auto& [v,b]:local)next[w+v]=mul(a,b);
            words.swap(next);
        }
        Mat m={rng.entry(),rng.entry(),rng.entry(),rng.entry()};map<string,Mat> contributions;
        for(auto& [w,c]:words){
            contributions[w]=plus(contributions[w],scale(m,c));
            string rev(w.rbegin(),w.rend());
            contributions[rev]=plus(contributions[rev],scale(adj(m),w.size()%2?mod(-c):c));
            maxlen=std::max(maxlen,int(w.size()));
        }
        for(auto& [w,matrix]:contributions)tiles[w].push_back({bit,matrix});
    }
    Recipe()=default;
    void initScalars(RNG& rng){
        for(int i=0;i<3;i++){
            int a=rng.entry();tiles[string(1,char('0'+i))].push_back({0,{a,rng.entry(),rng.entry(),mod(-a)}});
        }
    }
    void addExpanded(const map<string,int>& path,int bit,Mat m){
        map<string,Mat> row;
        for(const auto& [w,c]:path){
            row[w]=plus(row[w],scale(m,c));
            string rev(w.rbegin(),w.rend());
            row[rev]=plus(row[rev],scale(adj(m),w.size()%2?mod(-c):c));
            maxlen=std::max(maxlen,int(w.size()));
        }
        for(auto& [w,value]:row)tiles[w].push_back({bit,value});
    }
    Recipe(int seed,int top,int colours,int mode){
        RNG rng(10000019ULL*seed+324553);
        initScalars(rng);
        if(mode==11){
            assert(top-3>colours);
            array<int,3> u{},v{};
            for(auto& a:u)a=rng.entry();
            for(auto& a:v)a=rng.entry();
            // Q0 tr(P_uv X_u) = tr([X_u,X_v] X_u) = 0 by cyclicity.
            for(int i=0;i<3;i++)for(int j=i+1;j<3;j++){
                int c=mod(int64_t(u[i])*v[j]-int64_t(u[j])*v[i]);
                for(int k=0;k<3;k++)addClosed(string{char('2'+i+j),char('0'+k)},mul(c,u[k]));
            }
            closedBit=1;marks=1;full=1;return;
        }
        if(mode==2||mode==3||mode==4||mode==5){
            assert(top-(mode==4?2:3)>colours);
            // The determinant contains either fourteen singletons (one P), or
            // fourteen scalar singletons. Q0(P) ranges over the commutators.
            array<int,3> u{},v{};
            for(auto& x:u)x=rng.entry();
            for(auto& x:v)x=rng.entry();
            if(mode!=3){
                for(int i=0;i<3;i++)for(int j=0;j<3;j++)
                    addClosed(string{char('0'+i),char('0'+j)},mul(u[i],v[j]));
                if(mode==4){closedBit=1;marks=1;}
                else{
                    addPath(mode==2?"P":"ss",1,rng);closedBit=2;marks=2;
                }
            }else{
                for(int i=0;i<3;i++)for(int j=i+1;j<3;j++){
                    int coeff=rng.entry();
                    for(int k=0;k<3;k++){
                        int c=mul(coeff,v[k]);
                        addClosed(string{char('0'+i),char('0'+j),char('0'+k)},c);
                        addClosed(string{char('0'+j),char('0'+i),char('0'+k)},mod(-c));
                    }
                }
                closedBit=1;marks=1;
            }
            full=(1<<marks)-1;return;
        }
        bool derivative=mode==1;
        int arity=top-int(derivative);
        int excess=arity-colours-1;
        assert(excess>=0 && excess<=3);
        int cost=seed%5<3?excess:rng.pick(excess+1);
        vector<string> shapes;
        while(cost){
            int c=1+rng.pick(cost);shapes.push_back(string(c+1,'s'));cost-=c;
        }
        if(derivative){
            int choices=1;for(auto& shape:shapes)choices+=shape.size();
            int chosen=rng.pick(choices);
            if(chosen==0)shapes.push_back("P");
            else{--chosen;for(auto& shape:shapes){if(chosen<int(shape.size())){shape[chosen]='P';break;}chosen-=shape.size();}}
        }
        marks=shapes.size();assert(marks<=4);full=(1<<marks)-1;
        for(int k=0;k<marks;k++)addPath(shapes[k],1<<k,rng);
    }
    Poly trace(const string& w) const{
        // Sum over circular tilings on labelled sites. Assign a tiling uniquely by
        // the tile containing site 0. This cancels the logarithm's 1/(number of tiles)
        // against the permitted cuts. The omitted 1/(2*length) is constant per probe.
        Poly values{};int len=w.size();string twice=w+w;
        for(int size=1;size<=std::min(len,maxlen);size++)for(int before=0;before<size;before++){
            int start=(len-before)%len;
            string first=twice.substr(start,size);
            auto hit=tiles.find(first);if(hit==tiles.end())continue;
            string rotated=twice.substr(start,len);
            vector<array<Mat,16>> dp(len+1);
            for(auto& [mask,m]:hit->second)dp[size][mask]=plus(dp[size][mask],m);
            for(int pos=size;pos<len;pos++)for(int k=1;k<=std::min(maxlen,len-pos);k++){
                auto it=tiles.find(rotated.substr(pos,k));if(it==tiles.end())continue;
                for(auto& [mark,m]:it->second)for(int mask=0;mask<=full;mask++)if(!(mask&mark)){
                    Mat v=dp[pos][mask];if(!(v.a||v.b||v.c||v.d))continue;
                    dp[pos+k][mask|mark]=plus(dp[pos+k][mask|mark],scale(times(v,m),prime-1));
                }
            }
            for(int mask=0;mask<=full;mask++)values[mask]=add(values[mask],add(dp[len][mask].a,dp[len][mask].d));
        }
        if(closedBit){
            auto it=closed.find(w);
            if(it!=closed.end())values[closedBit]=add(values[closedBit],it->second);
        }
        return values;
    }
    vector<int> evaluate(const ProbeBank& bank) const{
        vector<Poly> tr;tr.reserve(bank.traces.size());
        for(const auto& w:bank.traces){
            auto pos=w.find_first_of("345");
            if(!varyScalars||pos==string::npos){tr.push_back(trace(w));continue;}
            int k=w[pos]-'3';int i=(k==2?1:0),j=(k==0?1:2);
            string replaced=w;replaced[pos]='0'+j;auto u=trace(replaced);
            replaced[pos]='0'+i;auto v=trace(replaced);Poly combined{};
            for(int m=0;m<=full;m++)combined[m]=mod(int64_t(variationAxis[i])*u[m]-int64_t(variationAxis[j])*v[m]);
            tr.push_back(combined);
        }
        vector<int> term;term.reserve(bank.terms.size());
        for(auto& ids:bank.terms){
            Poly v{};v[0]=1;
            for(int id:ids){
                Poly next{};
                for(int mask=0;mask<=full;mask++)if(v[mask]){
                    int avail=full^mask,sub=avail;
                    for(;;){next[mask|sub]=add(next[mask|sub],mul(v[mask],tr[id][sub]));if(!sub)break;sub=(sub-1)&avail;}
                }
                v=next;
            }
            term.push_back(v[full]);
        }
        vector<int> row;
        for(auto& p:bank.probes){int64_t v=0;for(auto [j,c]:p)v+=(int64_t)c*term[j];row.push_back(mod(v));}
        return row;
    }
};
// Relations with one field of three flavour slots, or two fields with two
// flavour slots. In the latter case the two possible variations have opposite
// signs after the two odd polarization parameters have been extracted.
struct LowerRecipe{
    vector<pair<int,Recipe>> variations;
    LowerRecipe(int seed,int top,int colours,int mode){
        RNG rng(7123991ULL*seed+185869);
        Recipe base;base.initScalars(rng);
        const bool inputOne=mode==8||mode==9,three=mode==6;
        int arity=top-(inputOne?1:2),excess=arity-colours-1;
        assert(excess>=0&&excess<=2);
        vector<string> paths;
        int cost=seed%3?excess:rng.pick(excess+1);
        while(cost){int k=1+rng.pick(cost);paths.push_back(string(k+1,'s'));cost-=k;}
        vector<char> inserts=(three?vector<char>{'f'}:(inputOne?vector<char>{'a'}:vector<char>{'a','b'}));
        for(char ins:inserts){
            vector<pair<int,int>> choices={{-1,-1}};
            for(int i=0;i<int(paths.size());i++)for(int j=0;j<int(paths[i].size());j++)
                if(paths[i][j]=='s')choices.push_back({i,j});
            auto [i,j]=choices[rng.pick(choices.size())];
            if(i<0)paths.emplace_back(1,ins);else paths[i][j]=ins;
        }
        assert(paths.size()<=4);
        base.marks=paths.size();base.full=(1<<base.marks)-1;
        map<char,array<int,3>> components;
        for(char l:string("ab"))for(int& c:components[l])c=rng.entry();
        vector<vector<array<int,3>>> scalars(paths.size());
        vector<Mat> endpoints;
        for(int i=0;i<int(paths.size());i++){
            for(char l:paths[i]){(void)l;scalars[i].push_back({rng.entry(),rng.entry(),rng.entry()});}
            endpoints.push_back({rng.entry(),rng.entry(),rng.entry(),rng.entry()});
        }
        for(char differentiated:(inputOne?(mode==8?string("-"):string("a")):(three?string("f"):string("ab")))){
            Recipe result=base;
            for(int n=0;n<int(paths.size());n++){
                map<string,int> expanded={{"",1}};
                for(int t=0;t<int(paths[n].size());t++){
                    char l=paths[n][t];map<string,int> atom;
                    if(l=='s'){
                        for(int i=0;i<3;i++)atom[string(1,char('0'+i))]=scalars[n][t][i];
                    }else if(l=='f'){
                        // Q0 f = -[X0,P12]+[X1,P02]-[X2,P01].
                        for(int i=0;i<3;i++){
                            char a='0'+i,b='5'-i;int c=i==1?1:prime-1;
                            atom[string{a,b}]=c;atom[string{b,a}]=mod(-c);
                        }
                    }else for(int i=0;i<3;i++)for(int j=i+1;j<3;j++){
                        int c=components[l][i+j-1];
                        if(l==differentiated){
                            char a='0'+i,b='0'+j;
                            atom[string{a,b}]=c;atom[string{b,a}]=mod(-c);
                        }else atom[string(1,char('2'+i+j))]=c;
                    }
                    map<string,int> next;
                    for(auto& [w,c]:expanded)for(auto& [v,d]:atom)next[w+v]=add(next[w+v],mul(c,d));
                    expanded.swap(next);
                }
                result.addExpanded(expanded,1<<n,endpoints[n]);
            }
            variations.emplace_back(differentiated=='b'?prime-1:1,std::move(result));
        }
    }
    vector<int> evaluate(const ProbeBank& bank)const{
        vector<int> result(bank.probes.size());
        for(const auto& [sign,recipe]:variations){
            auto r=recipe.evaluate(bank);
            for(int j=0;j<int(r.size());j++)result[j]=add(result[j],mul(sign,r[j]));
        }
        return result;
    }
};
int main(int argc,char** argv){
    if(argc==4&&string(argv[1])=="--product-parameters"){
        int start=std::stoi(argv[2]),count=std::stoi(argv[3]);
        for(int i=start;i<start+count;i++){
            RNG rng(10000019ULL*i+324553);
            std::cout<<i;
            for(int k=0;k<15;k++)std::cout<<' '<<rng.pick(9)-4;
            std::cout<<'\n';
        }
        return 0;
    }
    if(argc!=9){std::cerr<<"usage: binary colours epsilons extra probes first_seed count mode(0=raw,1=Q0,2=closed_ss,3=closed_Ps,4/5=raw_ss) prime\n";return 2;}
    int n=std::stoi(argv[1]),eps=std::stoi(argv[2]),extra=std::stoi(argv[3]),probes=std::stoi(argv[4]);
    int start=std::stoi(argv[5]),count=std::stoi(argv[6]),mode=std::stoi(argv[7]);
    assert(mode>=0&&mode<=12);
    prime=std::stoi(argv[8]);assert(prime>=1000000&&prime<1100000);
    int top=3*eps+extra;ProbeBank bank(top,eps,(mode>=6&&mode<=8)||mode>=10);bank.build(probes,523483+19*top);
    std::cerr<<"probes "<<bank.probes.size()<<" distinct traces "<<bank.traces.size()<<" terms "<<bank.terms.size()<<"\n";
    for(int i=start;i<start+count;i++){
        vector<int> row;
        if(mode==10||mode==12){
            Recipe recipe(i,top-1,n,mode==10?0:4);recipe.varyScalars=true;
            RNG rng(568559ULL+19237ULL*i);
            for(int& a:recipe.variationAxis)a=rng.entry();
            row=recipe.evaluate(bank);
        }else row=(mode>=6&&mode<=9)?LowerRecipe(i,top,n,mode).evaluate(bank):Recipe(i,top,n,mode).evaluate(bank);
        std::cout<<i;for(int x:row)std::cout<<' '<<x;std::cout<<'\n';
    }
}
