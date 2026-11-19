// Contract the 8-slot determinant with ten ordinary skew-colour matrices with homogeneous Grassmann entries.
// Input layout: two ordered sets of eight distinct slots, followed by two fixed metric pairs.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using Key=array<int,5>;
struct Block{int code,len,oid[2];};
int wt[5]={1,8,64,512,4096};
int main(int argc,char**argv){
 assert(argc==2);ifstream table(argv[1]);int n;table>>n;assert(n==224);map<Key,int> ids;vector<array<long long,6>> cov(n);
 for(int i=0;i<n;i++){int l;table>>l;Key k{};for(int j=0;j<l;j++)table>>k[j];sort(k.begin(),k.end());assert(!ids.count(k));ids[k]=i;for(auto &v:cov[i])table>>v;}
 int labs[10],sl[20],mate[20];array<int,8> perm{0,1,2,3,4,5,6,7};vector<pair<array<int,8>,int>> permutations;
 do{int v=0;for(int i=0;i<8;i++)for(int j=0;j<i;j++)v+=perm[i]<perm[j];permutations.push_back({perm,1-2*(v%2)});}while(next_permutation(perm.begin(),perm.end()));
 while(cin>>labs[0]){for(int i=1;i<10;i++)cin>>labs[i];for(int &v:sl)cin>>v;int idodd[10],nb=0;for(int i=0;i<10;i++)idodd[i]=labs[i]<2?nb++:-1;assert(nb==4);
 vector<long long> raw(n);mate[sl[16]]=sl[17];mate[sl[17]]=sl[16];mate[sl[18]]=sl[19];mate[sl[19]]=sl[18];
 for(auto &pp:permutations){for(int i=0;i<8;i++){int a=sl[i],b=sl[8+pp.first[i]];mate[a]=b;mate[b]=a;}
 int sign=pp.second,visited=0,blocks=0;Block bl[5];bool ok=true;
 for(int i=0;i<10&&ok;i++){if(visited&(1<<i))continue;int k=2*i;Block b{0,0,{-1,-1}};
 do{int l=k/2,s=labs[l];assert(!(visited&(1<<l)));visited|=1<<l;if(k%2)sign=-sign;b.len++;b.code+=wt[s];if(s<2){if(b.oid[s]!=-1){ok=false;break;}b.oid[s]=idodd[l];}k=mate[k^1];}while(k!=2*i);
 if(b.len%2){ok=false;break;}if(ok){assert(blocks<5);bl[blocks++]=b;}}
 if(!ok)continue;
 // Python's generator order is by degree, then the count vector lexicographically (odd types first).
 auto order=[](const Block&a,const Block&b){if(a.len!=b.len)return a.len<b.len;for(int i=0;i<5;i++){int x=(a.code>>(3*i))&7,y=(b.code>>(3*i))&7;if(x!=y)return x<y;}return false;};sort(bl,bl+blocks,order);
 Key ky{};int odd[4],idx=0;
 for(int b=0;b<blocks;b++){ky[b]=bl[b].code;const auto&o=bl[b].oid;if(b&&ky[b]==ky[b-1]&&((o[0]>=0)^(o[1]>=0))){ok=false;break;}for(int j=0;j<2;j++)if(o[j]>=0)odd[idx++]=o[j];}
 if(!ok)continue;assert(idx==4);for(int i=0;i<4;i++)for(int j=0;j<i;j++)if(odd[i]<odd[j])sign=-sign;
 sort(ky.begin(),ky.end());auto it=ids.find(ky);assert(it!=ids.end());raw[it->second]+=sign;
 }
 array<long long,6> out{};for(int i=0;i<n;i++){cout<<raw[i]<<' ';for(int j=0;j<6;j++)out[j]+=raw[i]*cov[i][j];}for(auto v:out)cout<<v<<' ';cout<<'\n';}
}
