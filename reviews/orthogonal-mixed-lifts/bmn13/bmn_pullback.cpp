// Pull a trace covector through the flavour projector and six odd substitutions.
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
int fieldPrime;
int rem(long long a){a%=fieldPrime;return a<0?a+fieldPrime:a;}
int inv(int a){int k=fieldPrime-2,r=1;while(k){if(k%2)r=rem(1LL*r*a);a=rem(1LL*a*a);k/=2;}return r;}
int fac(int n){int a=1;for(int i=2;i<=n;i++)a=rem(1LL*a*i);return a;}
bool odd(char c){return c=='3'||c=='4';}
struct Tr{int sign=0,stabilizer=0,parity=0;string canonical;};
unordered_map<string,Tr> traces;
Tr one(const string& w){
 auto it=traces.find(w);if(it!=traces.end())return it->second;
 int n=w.size(),o=count_if(w.begin(),w.end(),odd);map<string,int> orbit;bool zero=false;
 for(int r=0;r<=1;r++){
  string s=r?string(w.rbegin(),w.rend()):w;int a=r && (n+o*(o-1)/2)%2?-1:1;
  for(int k=0;k<n;k++){
   string t=s.substr(k)+s.substr(0,k);
   auto [it,fresh]=orbit.emplace(t,a);if(!fresh&&it->second!=a)zero=true;
   if(odd(s[k])&&(o-1)%2)a=-a;
  }
 }
 int orient=orbit.begin()->second;string key=orbit.begin()->first;
 for(auto& [t,v]:orbit)traces[t]={zero?0:v*orient,int(2*n/orbit.size()),o%2,key};
 return traces.at(w);
}
vector<string> split(const string& w){vector<string> out(1);for(char c:w){if(c=='/')out.emplace_back();else out.back()+=c;}return out;}
struct Word{string key;int sign=1,metric=1;};
Word canon(const vector<string>& w){
 Word out;vector<Tr> list;
 for(const string& t:w){
  Tr u=one(t);if(!u.sign){out.sign=0;return out;}
  for(auto& prev:list)if(u.parity&&prev.parity){
    if(u.canonical==prev.canonical){out.sign=0;return out;}
    if(prev.canonical>u.canonical)out.sign=-out.sign;
  }
  list.push_back(u);out.sign*=u.sign;out.metric=rem(1LL*out.metric*u.stabilizer);
 }
 sort(list.begin(),list.end(),[](const auto& a,const auto& b){return a.canonical<b.canonical;});
 int repeat=0;string last;
 for(auto& u:list){
  if(u.canonical!=last){last=u.canonical;repeat=0;}
  out.metric=rem(1LL*out.metric*(++repeat));
  if(!out.key.empty())out.key+='/';out.key+=u.canonical;
 }
 return out;
}
using Table=unordered_map<string,int>;
void accumulate(Table& table,const Word& w,long long coef){if(w.sign){int& v=table[w.key];v=rem(v+w.sign*coef);}}
void clean(Table& table){for(auto it=table.begin();it!=table.end();){if(!it->second)it=table.erase(it);else ++it;}}
int sg(const array<int,3>& p){return ((p[0]>p[1])+(p[0]>p[2])+(p[1]>p[2]))%2?-1:1;}
int cell(int i,int j){return 1<<(3*(3*i+j));}
map<int,int> projector(){
 map<int,int> dual;array<int,3> a={0,1,2};
 do{array<int,3> b={0,1,2};do{
  int m=2;
  for(int i=0;i<3;i++)m+=cell(i,a[i])+cell(i,b[i]);
  dual[m]+=sg(a)*sg(b);
 }while(next_permutation(b.begin(),b.end()));}while(next_permutation(a.begin(),a.end()));
 int norm=0;
 for(auto& [m,c]:dual){int v=c;for(int k=0;k<9;k++)v*=fac((m>>(3*k))&7);norm+=c*v;c=v;}
 int scale=inv(norm);
 for(auto& [m,v]:dual)v=rem(1LL*v*scale);
 return dual;
}
Table highest(const Table& initial){
 auto proj=projector();Table result;
 for(const auto& [w,c]:initial){
  auto ts=split(w);vector<pair<int,int>> slots;vector<int> old;
  for(int i=0;i<int(ts.size());i++)for(int j=0;j<int(ts[i].size());j++)
    if(ts[i][j]<'3'){slots.push_back({i,j});old.push_back(ts[i][j]-'0');}
  assert(slots.size()==8);vector<int> assignment={0,0,0,0,1,1,2,2};
  do{
   int m=0;for(int j=0;j<8;j++)m+=cell(old[j],assignment[j]);
   auto it=proj.find(m);if(it==proj.end()||!it->second)continue;
   for(int j=0;j<8;j++)ts[slots[j].first][slots[j].second]='0'+assignment[j];
   accumulate(result,canon(ts),1LL*c*it->second);
  }while(next_permutation(assignment.begin(),assignment.end()));
 }
 clean(result);return result;
}
Table lower(const Table& input,int spin,int flav){
 Table result;
 for(auto& [w,c]:input){
  auto ts=split(w);int sign=1;
  for(auto& tr:ts)for(auto& l:tr){
   if(l==spin){l=flav;accumulate(result,canon(ts),sign*c);l=spin;}
   if(odd(l))sign=-sign;
  }
 }
 clean(result);return result;
}
int main(int argc,char** argv){
 if(argc!=4){cerr<<"usage: bmn_pullback prime input.dat output.dat\n";return 2;}
 fieldPrime=stoi(argv[1]);ifstream in(argv[2]);assert(in);Table t;int c;string w;
 while(in>>c>>w){auto u=canon(split(w));assert(u.sign);accumulate(t,u,1LL*rem(c)*inv(u.metric));}
 auto p=highest(t);assert(highest(p)==p);t=std::move(p);cerr<<"highest weights "<<t.size()<<'\n';
 for(int k=5;k>=0;k--){t=lower(t,'3'+k/3,'0'+k%3);cerr<<"reverse "<<k<<" support "<<t.size()<<'\n';}
 ofstream out(argv[3]);assert(out);vector<pair<string,int>> sorted(t.begin(),t.end());sort(sorted.begin(),sorted.end());
 for(auto& [s,v]:sorted)out<<rem(1LL*v*canon(split(s)).metric)<<' '<<s<<'\n';
 return 0;
}
