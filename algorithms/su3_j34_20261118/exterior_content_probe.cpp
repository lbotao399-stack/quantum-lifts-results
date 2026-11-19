// Exact specialization modulo 1000003 of Tr on traceless three-colour matrices.
// Each odd species takes values in the exterior algebra of seven generators.
// Its degree-one coefficient vectors are independent deterministic integers.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;
namespace{
constexpr int p=1000003,d=7,B=1<<d;
int grade[B],sg[B][B];vector<int> masks[d+1];
void check(bool yes,const char* text){if(!yes)throw runtime_error(text);}
struct Poly{int n=0;array<int,B> a{};};
int reduce(long long x){x%=p;return x<0?int(x)+p:int(x);}
Poly times(const Poly& a,const Poly& b){Poly r;r.n=a.n+b.n;check(r.n<=d,"exterior degree");for(int i:masks[r.n]){long long sum=0;for(int j:masks[a.n])if((j&i)==j)sum+=(long long)sg[j][i^j]*a.a[j]*b.a[i^j];r.a[i]=reduce(sum);}return r;}
uint32_t draw(uint32_t &u){u^=u<<13;u^=u>>17;u^=u<<5;return u;}
bool odd[31];Poly M[31][3][3];unordered_map<string,Poly> traces;
void init(uint32_t seed){for(int i=0;i<B;i++){grade[i]=__builtin_popcount((unsigned)i);masks[grade[i]].push_back(i);for(int j=0;j<B;j++){int k=0;for(int q=0;q<d;q++)if((i>>q)&1)k+=__builtin_popcount((unsigned)(j&((1<<q)-1)));sg[i][j]=((i&j)?0:(k%2?-1:1));}}
 for(int l=0;l<31;l++)for(int i=0;i<3;i++)for(int j=0;j<3;j++){auto &a=M[l][i][j];a.n=odd[l];for(int k:masks[a.n])a.a[k]=(i==2&&j==2)?reduce(-M[l][0][0].a[k]-M[l][1][1].a[k]):int(draw(seed)%23)-11;}}
const Poly& tr(const string &key){auto it=traces.find(key);if(it!=traces.end())return it->second;array<Poly,9> acc{};for(int i=0;i<3;i++)acc[4*i].a[0]=1;check(!key.empty()&&key.size()%2==0,"trace key");for(unsigned h=0;h<key.size();h+=2){int l=stoi(key.substr(h,2),nullptr,16);check(0<=l&&l<31,"letter");array<Poly,9> next{};int n=acc[0].n;for(auto &x:next)x.n=n+odd[l];check(n+odd[l]<=d,"trace degree");for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k:masks[n+odd[l]]){long long sum=0;if(!odd[l]){for(int u=0;u<3;u++)sum+=(long long)acc[3*i+u].a[k]*M[l][u][j].a[0];}else for(int q=0;q<d;q++)if((k>>q)&1){int m=k^(1<<q);for(int u=0;u<3;u++)sum+=(long long)sg[m][1<<q]*acc[3*i+u].a[m]*M[l][u][j].a[1<<q];}next[3*i+j].a[k]=reduce(sum);}acc=next;}
 Poly out;out.n=acc[0].n;for(int m:masks[out.n])out.a[m]=reduce(acc[0].a[m]+acc[4].a[m]+acc[8].a[m]);return traces.emplace(key,out).first->second;}
}
int main(int argc,char**argv){try{ios::sync_with_stdio(false);cin.tie(nullptr);check(argc==3,"use: exterior_content_probe seed odd_hex < words");string os=argv[2];for(unsigned h=0;h<os.size();h+=2)odd[stoi(os.substr(h,2),nullptr,16)]=true;init(stoul(argv[1]));string g;int n,c,degree;while(cin>>g>>degree>>n){check(0<=degree&&degree<=d&&n>0,"header");traces.clear();Poly sum;sum.n=degree;for(int j=0;j<n;j++){int m;check(bool(cin>>c>>m)&&m>0&&m<7,"input");Poly v;v.a[0]=1;for(int k=0;k<m;k++){string t;check(bool(cin>>t),"word");v=times(v,tr(t));}check(v.n==degree,"bad parity");for(int k:masks[degree])sum.a[k]=reduce((long long)sum.a[k]+(long long)c*v.a[k]);}
 cout<<g<<' '<<degree;for(int k:masks[degree])cout<<' '<<sum.a[k];cout<<'\n';}check(cin.eof(),"incomplete input");return 0;}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
