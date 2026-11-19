// Exact invariant-polynomial test on sl(3). Letter 00 is even; restrict it to
// diag(x,y,-x-y). All other entries remain general, with one linear constraint
// on each diagonal. Six nonzero-letter occurrences, six copies of 00.
// The output counts polynomial coefficients, not evaluation points.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
using namespace std;
using Int=int64_t; using Key=uint64_t;
void check(bool b,const string &m){if(!b)throw runtime_error(m);}
using Terms=unordered_map<Key,Int>;
using Mono=pair<Key,int>;
struct Trace{int n=0;vector<int> axial; vector<tuple<int,int,int>> entries;};
struct Plan{int vertices=0,order=0,diagonal=0;vector<int> powers; vector<tuple<int,int,int>> edges;};
int code(const string &s,size_t i){return stoi(s.substr(i,2),nullptr,16);}
// Label ordering/parities: sorted(letters_up_to(1,1)), as recorded by the
// Python driver. The odd codes are explicitly supplied by that driver.
set<int> odd;
struct Term{Int coefficient;vector<string> traces;};
int binom(int a,int b){int r=1;for(int i=1;i<=b;++i)r=r*(a-i+1)/i;return r;}
Plan makeplan(const vector<string>&w){Plan p;
 for(auto &s:w){int N=int(s.size()/2);vector<int> a(N),col(N,-1);int count=0;
  for(int i=0;i<N;i++){a[i]=code(s,2*i);if(a[i]){col[i]=p.vertices++;count++;}else p.diagonal++;}
  if(!count){p.powers.push_back(N);p.vertices++;continue;}
  for(int i=0;i<N;i++)if(a[i]){int j=(i+1)%N,gap=0;while(!a[j]){gap++;j=(j+1)%N;}check(col[i]>=0&&col[j]>=0,"edge");p.edges.emplace_back(a[i],col[i],col[j]);p.order++;p.powers.push_back(0);}
  for(int i=0;i<N;i++)if(a[i]){int j=(i+1)%N,gap=0;while(!a[j]){gap++;j=(j+1)%N;}p.powers[col[j]]+=gap;}
 }
 check(int(p.powers.size())==p.vertices && p.order==6 && p.diagonal==6,"content outside supported domain");return p;}
struct Engine{Terms packed;Int leaves=0,bound=0;vector<int> vs;array<int,6> vars{},sites{};Plan p;Int weight;
 static Key key(const array<int,6>&v,int e){Key k=0;for(int i:v)k=(k<<8)|unsigned(i);return (k<<3)|e;}
 void monomial(int i,int sg,int x1,int x2){
  if(i==6){int inv=0;for(int a=0;a<6;a++)for(int b=a+1;b<6;b++)if(odd.count(vars[a]/8)&&odd.count(vars[b]/8)){if(vars[a]==vars[b])return;inv+=vars[a]>vars[b];}
   auto sorted=vars;sort(sorted.begin(),sorted.end());Int c=(inv%2?-sg:sg)*weight;if(x2%2)c=-c;
   for(int e=0;e<=x2;e++){Int v=c*binom(x2,e);packed[key(sorted,x1+e)]+=v;leaves++;bound+=llabs(v);}return;}
  auto [l,a,b]=p.edges[i];int e=vs[a]*3+vs[b];if(e!=8){vars[i]=l*8+e;monomial(i+1,sg,x1,x2);}else{
   vars[i]=l*8;monomial(i+1,-sg,x1,x2);vars[i]=l*8+4;monomial(i+1,-sg,x1,x2);}}
 void colors(int k){if(k==p.vertices){int x1=0,x2=0;for(int i=0;i<k;i++){if(vs[i]==0)x1+=p.powers[i];if(vs[i]==2)x2+=p.powers[i];}monomial(0,1,x1,x2);return;}
  for(int i=0;i<3;i++){vs[k]=i;colors(k+1);}}
 void apply(const Term &t){p=makeplan(t.traces);weight=t.coefficient;vs.resize(p.vertices);colors(0);}
};
int main(int argc,char**argv){try{
 check(argc==4||argc==5,"usage: colour_polynomial input.txt odd_hex output.json [control]");string oddhex=argv[2];for(size_t i=0;i<oddhex.size();i+=2)odd.insert(code(oddhex,i));check(!odd.count(0),"diagonal scalar must be even");
 ifstream in(argv[1]);int n;check(bool(in>>n)&&n>0,"input");vector<Term> all;Int abs_total=0;__int128 worst=0;
 for(int i=0;i<n;i++){Term t;int nt;check(bool(in>>t.coefficient>>nt)&&nt>0,"trace count");for(int j=0;j<nt;j++){string s;check(bool(in>>s) && s.size()>=4 && !(s.size()%2),"trace");t.traces.push_back(s);for(size_t k=0;k<s.size();k+=2)check(code(s,k)<31,"letter id");}auto p=makeplan(t.traces);Int walks=1;for(int k=0;k<p.vertices;k++)walks*=3;
 // At most 2^6 choices for 33 entries; the l1 norm of the last diagonal factor <=2^6.
 worst+=(__int128)llabs(t.coefficient)*walks*4096;all.push_back(t);abs_total+=llabs(t.coefficient);}
 string extra;check(!(in>>extra)&&worst<numeric_limits<Int>::max(),"unproved int64 bound");Engine eng;for(auto &t:all)eng.apply(t);
 int nz=0;Int maxabs=0;Key first=numeric_limits<Key>::max();for(auto [k,v]:eng.packed)if(v){nz++;maxabs=max(maxabs,Int(llabs(v)));first=min(first,k);}check(eng.bound<=worst,"bound");
 Engine ctrl;if(argc==5){ctrl.apply(all[0]);}int cn=0;for(auto [k,v]:ctrl.packed)cn+=v!=0;
 ofstream out(argv[3]);check(bool(out),"output");out<<"{\"formal_terms\":"<<n<<",\"numerator_abs_sum\":"<<abs_total<<",\"distinct_collected_coordinates\":"<<eng.packed.size()<<",\"nonzero_coordinates\":"<<nz<<",\"maximum_abs_residual\":"<<maxabs<<",\"contributions\":"<<eng.leaves<<",\"accumulated_abs_contributions\":"<<eng.bound<<",\"a_priori_int64_bound\":"<<Int(worst)<<",\"first_residual_key\":"<<(nz?first:0)<<",\"first_residual\":"<<(nz?eng.packed.at(first):0)<<",\"first_word_control_nonzero_coordinates\":"<<cn<<"}\n";check(bool(out),"write");return 0;
 }catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
