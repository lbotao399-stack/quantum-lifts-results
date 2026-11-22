// Contract Q1 of the two explicit BMN product primitives with an ordered trace covector.
#define main bmn_pullback_driver
#include "bmn_pullback.cpp"
#undef main
#include <unordered_set>
using Mat=array<int,4>;
Mat mm(const Mat& x,const Mat& y){return {rem(1LL*x[0]*y[0]+1LL*x[1]*y[2]),rem(1LL*x[0]*y[1]+1LL*x[1]*y[3]),rem(1LL*x[2]*y[0]+1LL*x[3]*y[2]),rem(1LL*x[2]*y[1]+1LL*x[3]*y[3])};}
string rev(string s){reverse(s.begin(),s.end());return s;}
set<string> orb(const string& t){
 set<string> res;for(const string& w:{t,rev(t)})for(int i=0;i<int(w.size());i++)res.insert(w.substr(i)+w.substr(0,i));
 return res;
}
struct Source{
 array<Mat,3> m;array<int,3> u,v;
 unordered_map<string,int> logarithm,homotopy,aux;
 int au(const string& w){
  auto it=aux.find(w);if(it!=aux.end())return it->second;
  Mat a={1,0,0,1};for(char l:w)a=mm(a,m[l-'0']);
  int val=rem(a[0]+a[3]);if(w.size()<=8)aux[w]=val;return val;
 }
 int logcoef(const string& w){
  auto it=logarithm.find(w);if(it!=logarithm.end())return it->second;
  int value=rem(1LL*(w.size()%2?1:-1)*au(w)*inv(one(w).stabilizer));
  return logarithm[w]=value;
 }
 int h(const string& t){
  auto it=homotopy.find(t);if(it!=homotopy.end())return it->second;
  if(t.size()%2)return 0;
  long long c=0;
  for(const auto& w:orb(t)){
   auto pos=w.find_first_of("567");assert(pos!=string::npos);
   int k=w[pos]-'5';char l=k==2?'1':'0',r=k==0?'1':'2';
   string mid=w.substr(0,pos),tail=w.substr(pos+1);
   for(size_t j=0;j<=tail.size();j++){
    string prefix=rev(tail.substr(j)),end=tail.substr(0,j);
    c+=au(prefix+l+mid+r+end)-au(prefix+r+mid+l+end);
   }
  }
  return homotopy[t]=rem(1LL*rem(c)*inv(4*(1+t.size())));
 }
 int coefficient(const vector<string>& ts){
  string active;map<string,int> others;int odds=1;for(const string& t:ts){
   if(t.find_first_of("567")!=string::npos)active=t;else{others[t]++;odds+=t.size()%2;}
  }
  assert(!active.empty());int hp=h(active);if(!hp)return 0;
  int ans=0;
  for(auto& [b,ct]:others){
   if(b.size()!=2)continue;int i=b[0]-'0',j=b[1]-'0';
   int coef=u[i]*v[j]+(i!=j?u[j]*v[i]:0);if(!coef)continue;
   for(const auto& [w,n]:others){int power=n-(w==b);if(!power)continue;
    int l=logcoef(w);coef=rem(1LL*coef*inv(fac(power)));
    for(int m=0;m<power;m++)coef=rem(1LL*coef*l);
   }
   ans=rem(ans+coef);
  }
  return rem(1LL*ans*rem(1LL*hp*inv(odds)));
 }
};
void emit_parent(const vector<string>& ts,unordered_set<string>& parents){
 if(none_of(ts.begin(),ts.end(),[](const string& t){return t.size()==2;}))return;
 for(int a=0;a<3;a++)for(int choice=0;choice<2;choice++){
  auto all=ts;bool bilinear=false,ok=true;
  for(auto& t:all){
   if(t.size()<2){ok=false;break;}
   for(char& l:t)if(l=='3'||l=='4')l=(l=='3'+choice)?char('7'-a):char('0'+a);
   if(t.find_first_of("567")!=string::npos){if(t.size()%2){ok=false;break;}}
   else if(t.size()==2)bilinear=true;
  }
  if(ok&&bilinear){auto u=canon(all);if(u.sign)parents.insert(u.key);}
 }
}
void candidates(const string& w,unordered_set<string>& parents){
 auto ts=split(w);vector<int> marked;
 for(int i=0;i<int(ts.size());i++)if(ts[i].find_first_of("34")!=string::npos)marked.push_back(i);
 if(marked.size()==1){
  int i=marked[0];auto rest=ts;rest.erase(rest.begin()+i);
  for(const string& t:orb(ts[i])){
   for(int k=0;k<=int(t.size());k++){
    string x=t.substr(0,k),y=t.substr(k);
    auto m=rest;m.push_back(rev(x)+y);emit_parent(m,parents);
    if(x.size()>=2&&y.size()>=2&&x.find_first_of("34")!=string::npos&&y.find_first_of("34")!=string::npos){
     m=rest;m.push_back(x);m.push_back(y);emit_parent(m,parents);m.back()=rev(y);emit_parent(m,parents);
    }
   }
  }
 }
 for(int j=0;j<int(ts.size());j++)for(int k=j+1;k<int(ts.size());k++){
  if(find(marked.begin(),marked.end(),j)==marked.end()&&find(marked.begin(),marked.end(),k)==marked.end())continue;
  if(marked.size()==2 && (j!=marked[0]||k!=marked[1]))continue;
  auto m=ts;m.erase(m.begin()+k);m.erase(m.begin()+j);m.emplace_back();
  for(const string& u:orb(ts[j]))for(const string& v:orb(ts[k])){m.back()=u+v;emit_parent(m,parents);}
 }
}
int pull(const vector<string>& ts,const Table& dual,int n){
 int psi=-1,ip=-1,flav=-1;for(int i=0;i<int(ts.size());i++){
  auto j=ts[i].find_first_of("567");if(j!=string::npos){psi=i;ip=j;flav='7'-ts[i][j];}
 }
 assert(psi>=0);int total=0;
 auto add_variant=[&](vector<string> traces,int coeff){
  bool seenA=false;int orient=1;
  for(auto& tr:traces){if(tr.size()==1)return;for(char l:tr){if(l=='8')seenA=true;if(l=='9'&&!seenA)orient=-1;}}
  for(int switchSpins=0;switchSpins<2;switchSpins++){
   for(auto& tr:traces)for(char& l:tr){
    if(l=='8')l='3';else if(l=='9')l='4';else if(l=='3')l='4';else if(l=='4')l='3';
   }
   auto u=canon(traces);if(!u.sign)continue;
   auto it=dual.find(u.key);
   if(it!=dual.end())total=rem(total+1LL*coeff*orient*(switchSpins?-1:1)*u.sign*it->second);
  }
 };
 for(int q=0;q<int(ts.size());q++)for(int j=0;j<int(ts[q].size());j++){
  if(ts[q][j]!='0'+flav)continue;
  auto raw=ts;raw[psi][ip]='8';raw[q][j]='9';
  for(int a=0;a<2;a++)for(int b=0;b<2;b++){
   int s=((a+b)%2?-1:1)*(flav==1?1:-1);
   if(psi!=q){
    string A=raw[psi].substr(0,ip+a),B=raw[psi].substr(ip+a),D=raw[q].substr(0,j+b),E=raw[q].substr(j+b);
    vector<string> rest;for(int k=0;k<int(raw.size());k++)if(k!=psi&&k!=q)rest.push_back(raw[k]);
    rest.push_back(A+E+D+B);add_variant(rest,s);rest.back()=A+rev(D)+rev(E)+B;
    add_variant(rest,ts[q].size()%2?s:-s);
   }else{
    int lo=min(ip+a,j+b),hi=max(ip+a,j+b);
    string A=raw[q].substr(0,lo),B=raw[q].substr(lo,hi-lo),C=raw[q].substr(hi);
    auto rest=raw;rest.erase(rest.begin()+q);auto sp=rest;int factor=1;
    if(B.empty())factor*=n;else sp.push_back(B);
    if(A.empty()&&C.empty())factor*=n;else sp.push_back(A+C);
    add_variant(sp,s*factor);rest.push_back(A+rev(B)+C);add_variant(rest,B.size()%2?s:-s);
   }
  }
 }
 return total;
}
int main(int argc,char** argv){
 if(argc!=6&&argc!=7){cerr<<"usage: bmn_source_pair prime dual.dat parameters.dat number outputs.dat [--all]\n";return 2;}
 const bool all=argc==7&&string(argv[6])=="--all";assert(argc==6||all);
 fieldPrime=stoi(argv[1]);ifstream din(argv[2]),params(argv[3]);assert(din&&params);Table dual;string w;long long c;
 while(din>>c>>w){auto v=canon(split(w));assert(v.sign);accumulate(dual,v,rem(c));}
 clean(dual);vector<Source> src(stoi(argv[4]));
 for(auto& a:src){int idx;params>>idx;for(auto& m:a.m){params>>m[0]>>m[1]>>m[2];m[3]=-m[0];}
  for(int& x:a.u)params>>x;for(int& x:a.v)params>>x;}
 unordered_set<string> parents;int ix=0;
 for(auto& [s,v]:dual){candidates(s,parents);if((++ix)%1500==0)cerr<<"outputs "<<ix<<" candidate parents "<<parents.size()<<'\n';}
 cerr<<"full parent support "<<parents.size()<<'\n';
 vector<int> scores(src.size());int nonzero=0,seen=0;ofstream out(argv[5]);assert(out);
 vector<string> ordered(parents.begin(),parents.end());sort(ordered.begin(),ordered.end());
 for(const auto& p:ordered){
  auto t=split(p);vector<int> vals;bool any=false;for(auto& x:src){int v=x.coefficient(t);vals.push_back(v);any|=v!=0;}
  if(any||all){int c=pull(t,dual,13);if(c&&any)nonzero++;if((c&&any)||all){out<<p<<' '<<c;for(int j=0;j<int(src.size());j++){scores[j]=rem(scores[j]+1LL*c*vals[j]);out<<' '<<vals[j];}out<<'\n';}}
  if((++seen)%100000==0)cerr<<"paired "<<seen<<'\n';
 }
 cerr<<"contributing "<<nonzero<<" scores";for(int x:scores)cerr<<' '<<x;cerr<<'\n';
}
