// One target polynomial; many exact matrix substitutions on seven Grassmann generators.
#define J34_EXTERIOR_LIBRARY
#include "exterior_content_probe.cpp"
#include <fstream>
struct Term {int c; vector<unsigned> ids;};
int main(int argc,char** argv){try{
 ios::sync_with_stdio(false);cin.tie(nullptr);check(argc==4,"use: centralizer_scan words assignments odd_hex");
 string os=argv[3];for(size_t h=0;h<os.size();h+=2)odd[stoi(os.substr(h,2),nullptr,16)]=true;
 init(238121); ifstream w(argv[1]), a(argv[2]); int degree,count;check(bool(w>>degree>>count)&&degree>=0&&degree<=d&&count>0,"target");
 vector<Term> terms;unordered_map<string,unsigned> index;vector<string> names;int c,nt;
 for(int i=0;i<count;i++){check(bool(w>>c>>nt)&&0<=c&&c<p&&0<nt&&nt<7,"word");Term term{c,{}};
  for(int j=0;j<nt;j++){string s;check(bool(w>>s)&&s.size()>=4&&s.size()%2==0,"trace");auto v=index.emplace(s,names.size());if(v.second)names.push_back(s);term.ids.push_back(v.first->second);}terms.push_back(move(term));}
 string key,tail;check(!(w>>tail),"extra words");int nm;while(a>>key>>nm){check(nm>0,"matrix count");for(int h=0;h<nm;h++){int l,ix;check(bool(a>>l>>ix)&&l>=0&&l<31&&((!odd[l]&&ix==-1)||(odd[l]&&0<=ix&&ix<d)),"matrix key");
   long long trc=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++){int x;check(bool(a>>x),"matrix entry");if(i==j)trc+=x;M[l][i][j].a[ix<0?0:(1<<ix)]=reduce(x);}check(!trc,"traceless");}
  traces.clear();vector<Poly> pol;pol.reserve(names.size());for(const auto& name:names)pol.push_back(tr(name));
  Poly sum;sum.n=degree;for(const auto& t:terms){auto v=pol[t.ids[0]];for(unsigned i=1;i<t.ids.size();i++)v=times(v,pol[t.ids[i]]);check(v.n==degree,"degree");for(int k:masks[degree])sum.a[k]=reduce(sum.a[k]+(long long)t.c*v.a[k]);}
  // Try another assignment only if the previous one vanishes. The input is a single graph.
  for(int k:masks[degree])if(sum.a[k]){cout<<key<<' '<<k<<' '<<sum.a[k]<<'\n';return 0;}
 }check(a.eof(),"bad assignment");cout<<"none\n";return 0;
}catch(const exception& e){cerr<<e.what()<<'\n';return 1;}}
