// Exact integer coefficients of a BMN invariant of at most twelve letters on three colours.
// One even scalar is diagonal. Up to nine other matrix entries are encoded
// injectively by their exponent multiplicities (a sorted multiset), not a hash.
#define J34_COLOUR_ENGINE_ONLY
#include "colour_polynomial.cpp"
#include <sstream>
struct Wide {
    Plan plan; int scalar, n, length; int64_t weight; int species[7]{}, val[9]{};
    uint64_t choose[66][10]{}; vector<int> colour; unordered_map<uint64_t,int64_t> totals;
    unsigned long long visits=0; unsigned __int128 absolute=0;
    Wide(){ for(int i=0;i<66;i++){choose[i][0]=1;for(int j=1;j<=9&&j<=i;j++)choose[i][j]=i==0?0:choose[i-1][j-1]+choose[i-1][j];}}
    void monomial(int k,int sg,int xfirst,int xthird){
        if(k==n){
            int inv=0;for(int i=0;i<n;i++)if(odd.count(val[i]/8))for(int j=i+1;j<n;j++)if(odd.count(val[j]/8)){
                if(val[i]==val[j])return;
                inv+=(val[i]>val[j]);}
            int a[9];copy(val,val+n,a);sort(a,a+n);uint64_t idx=0;
            for(int j=0;j<n;j++)idx+=choose[a[j]+j][j+1];
            check(idx<choose[56+n-1][n],"rank out of bounds");
            int64_t base=((inv+xthird)%2?-sg:sg)*weight;
            for(int i=0;i<=xthird;i++){int64_t v=base*binom(xthird,i);auto key=(idx<<4)|(xfirst+i);
                totals[key]+=v;visits++;absolute+=uint64_t(v<0?-v:v);}
            return;
        }
        auto [l,a,b]=plan.edges[k];int code=3*colour[a]+colour[b];
        if(code!=8){val[k]=8*l+code;monomial(k+1,sg,xfirst,xthird);}
        else{val[k]=8*l;monomial(k+1,-sg,xfirst,xthird);val[k]=8*l+4;monomial(k+1,-sg,xfirst,xthird);}
    }
    void vertices(int k){
        if(k==plan.vertices){int x=0,z=0;for(int i=0;i<k;i++){x+=(colour[i]==0)*plan.powers[i];z+=(colour[i]==2)*plan.powers[i];}monomial(0,1,x,z);return;}
        for(int i=0;i<3;i++){colour[k]=i;vertices(k+1);}
    }
    void run(const Term&t){
        plan=Plan();for(auto &s:t.traces){int len=s.size()/2;vector<int> labs(len),ids(len,-1);int n0=0;
            for(int i=0;i<len;i++){labs[i]=code(s,2*i);check(0<=labs[i]&&labs[i]<7,"invalid species");if(labs[i]!=scalar){ids[i]=plan.vertices++;n0++;}else plan.diagonal++;}
            if(n0==0){plan.powers.push_back(len);plan.vertices++;continue;}
            for(int i=0;i<len;i++)if(labs[i]!=scalar){int j=(i+1)%len;while(labs[j]==scalar)j=(j+1)%len;
                plan.edges.emplace_back(labs[i],ids[i],ids[j]);plan.order++;plan.powers.push_back(0);}
            for(int i=0;i<len;i++)if(labs[i]!=scalar){int j=(i+1)%len,z=0;while(labs[j]==scalar){j=(j+1)%len;z++;}plan.powers[ids[j]]+=z;}
        }
        check(plan.diagonal+plan.order==length&&plan.order<=9&&int(plan.powers.size())==plan.vertices,"length/order");
        n=plan.order;weight=t.coefficient;colour.resize(plan.vertices);vertices(0);
    }
};
int main(int argc,char**argv){try{
    check(argc==3||argc==4,"source_closure words.txt output.json [length]");ifstream in(argv[1]);ofstream out(argv[2]);int scalar,count;string os;
    const int length=argc==4?stoi(argv[3]):12;check(length>=2&&length<=12,"expected 2 to 12 letters");
    check(bool(in>>scalar>>os>>count)&&count>0&&out.good(),"bad header");for(unsigned i=0;i<os.size();i+=2)odd.insert(code(os,i));check(!odd.count(scalar),"diagonal odd");
    vector<Term> words;unsigned __int128 mass=0;vector<int> content0;
    for(int j=0;j<count;j++){int nt;Term t;check(bool(in>>t.coefficient>>nt)&&nt>0&&t.coefficient!=0,"input");vector<int> c;
        for(int k=0;k<nt;k++){string s;check(bool(in>>s)&&s.size()>=4&&s.size()%2==0,"trace");for(unsigned i=0;i<s.size();i+=2)c.push_back(code(s,i));t.traces.push_back(s);}
        sort(c.begin(),c.end());if(j==0)content0=c;check(int(c.size())==length&&c==content0,"inhomogeneous input");mass+=uint64_t(t.coefficient<0?-t.coefficient:t.coefficient);words.push_back(t);}
    string extra;check(!(in>>extra),"trailing");Wide w;w.scalar=scalar;w.length=length;unsigned __int128 bound=0;int zeros=std::count(content0.begin(),content0.end(),scalar);
    for(auto &t:words){int pure=0;for(auto &s:t.traces){bool all=true;for(unsigned i=0;i<s.size();i+=2)all&=code(s,i)==scalar;pure+=all;}
        unsigned __int128 walk=1;for(int i=0;i<length-zeros+pure;i++)walk*=3;bound+=walk*(1u<<length)*uint64_t(t.coefficient<0?-t.coefficient:t.coefficient);}
    check(bound<uint64_t(INT64_MAX),"integer accumulator bound");for(auto &t:words)w.run(t);
    check(w.absolute<=bound,"realized sum above bound");size_t nz=0;uint64_t first=UINT64_MAX;for(auto [k,v]:w.totals)if(v){nz++;first=min(k,first);}
    out<<"{\"formal_words\":"<<count<<",\"distinct_monomials\":"<<w.totals.size()<<",\"nonzero_monomials\":"<<nz<<",\"first_key\":"<<(nz?first:0)<<",\"first_integer_coefficient\":"<<(nz?w.totals[first]:0)<<",\"contributions\":"<<w.visits<<",\"sum_abs_bounds\":"<<uint64_t(bound)<<",\"sum_abs_actual\":"<<uint64_t(w.absolute)<<"}\n";check(out.good(),"output");return 0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
