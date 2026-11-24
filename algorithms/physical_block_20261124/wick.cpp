// Exact USp(2r) Wick contractions, with tr(T_A T_B)=delta_AB.
// Elementary covariance, with J=[0,I;-I,0]:
// sum_A T_A(i,j) conjugate(T_A(k,l))
// = (delta_ik delta_jl + J_il J_kj)/2.
// The returned integer is 2^length times the colour/Fermi contraction.
#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>

namespace {
constexpr int MAX=12;
using Array=std::array<int,MAX>;
struct DSU {
  int parent[2*MAX],parity[2*MAX];
  DSU(int n){for(int i=0;i<n;++i){parent[i]=i;parity[i]=0;}}
  std::pair<int,int> find(int x){
    int p=0;while(parent[x]!=x){p^=parity[x];x=parent[x];}return {x,p};
  }
  bool unite(int x,int y,int p){
    auto a=find(x),b=find(y);
    if(a.first==b.first)return (a.second^b.second)==p;
    parent[a.first]=b.first;parity[a.first]=a.second^b.second^p;return true;
  }
};
std::unordered_map<std::string,int64_t> cache;
int64_t colour(const Array& left,const Array& right,int n,int size){
  std::string key;key.reserve(2*n+2);key.push_back(n);key.push_back(size);
  for(int i=0;i<n;++i)key.push_back(left[i]);
  for(int i=0;i<n;++i)key.push_back(right[i]);
  auto got=cache.find(key);if(got!=cache.end())return got->second;
  int64_t total=0;
  for(unsigned mask=0;mask<(1u<<n);++mask){
    DSU uf(2*n);int exponent=__builtin_popcount(mask);bool valid=true;
    for(int i=0;i<n&&valid;++i){
      if(mask&(1u<<i)){
        valid=uf.unite(i,n+right[i],1)&&uf.unite(left[i],n+i,1);
      }else valid=uf.unite(i,n+i,0)&&uf.unite(left[i],n+right[i],0);
    }
    if(!valid)continue;
    int powers[2*MAX]={};
    for(int i=0;i<n;++i)if(mask&(1u<<i)){
      auto a=uf.find(i),b=uf.find(left[i]);
      powers[a.first]^=1;powers[b.first]^=1;exponent+=a.second+b.second;
    }
    int64_t contribution=(exponent&1)?-1:1;
    for(int i=0;i<2*n;++i)if(uf.parent[i]==i){
      if(powers[i]){contribution=0;break;}contribution*=size;
    }
    total+=contribution;
  }
  cache.emplace(std::move(key),total);return total;
}

int64_t pair_value(const int* letters_a,const int* succ_a,const int* map_a,
                   const int* letters_b,const int* succ_b,const int* parity,
                   int n,int size){
  Array matching{},inverse{},left{},right{};
  for(int i=0;i<n;++i)left[map_a[i]]=map_a[succ_a[i]];
  int64_t total=0;
  auto visit=[&](auto&& self,int pos,unsigned used)->void{
    if(pos==n){
      for(int i=0;i<n;++i)inverse[matching[i]]=i;
      for(int i=0;i<n;++i)right[map_a[i]]=map_a[inverse[succ_b[matching[i]]]];
      int sign=0;
      for(int i=0;i<n;++i)if(parity[letters_a[i]])
        for(int j=i+1;j<n;++j)if(parity[letters_a[j]]&&matching[i]>matching[j])++sign;
      int64_t value=colour(left,right,n,size);
      total+=(sign&1)?-value:value;return;
    }
    for(int j=0;j<n;++j)if(!(used&(1u<<j))&&letters_a[pos]==letters_b[j]){
      matching[pos]=j;self(self,pos+1,used|(1u<<j));
    }
  };
  visit(visit,0,0);return total;
}
}

extern "C" int wick_gram(int n,int size,int count,const int* letters,
                          const int* successors,const int* canonical_maps,
                          const int* parities,int64_t* out){
  if(n>MAX||n<1||size<2||size%2)return -1;
  try{
    for(int i=0;i<count;++i)for(int j=0;j<=i;++j){
      int64_t value=pair_value(letters+i*n,successors+i*n,canonical_maps+i*n,
                              letters+j*n,successors+j*n,parities,n,size);
      out[i*count+j]=out[j*count+i]=value;
    }
  }catch(...){return -2;}
  return 0;
}
extern "C" uint64_t wick_cache_size(){return cache.size();}
extern "C" void wick_clear_cache(){cache.clear();}
