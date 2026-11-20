// Export each trace polynomial separately. Text rows: number of nonzero
// monomials followed by (packed monomial, signed integer coefficient) pairs.
#define J34_COLOUR_ENGINE_ONLY
#include "colour_polynomial.cpp"
int main(int argc,char **argv){try{
 check(argc==4||argc==5,"usage: colour_columns input.txt odd_hex output.txt [input10]");bool parent=argc==5;check(!parent||string(argv[4])=="input10","mode");
 string os=argv[2];for(size_t i=0;i<os.size();i+=2)odd.insert(code(os,i));check(!odd.count(0),"odd diagonal");
 ifstream in(argv[1]);ofstream out(argv[3]);int count;check(bool(in>>count)&&count>0&&out.good(),"input or output");
 for(int i=0;i<count;++i){int nt;Term t;t.coefficient=1;check(bool(in>>nt)&&nt>0&&nt<=6,"nt");
  for(int j=0;j<nt;j++){string s;check(bool(in>>s)&&s.size()>=4&&s.size()%2==0,"trace");t.traces.push_back(s);for(size_t k=0;k<s.size();k+=2){int l=code(s,k);check(l>=0&&l<31,"letter");}}
  auto p=makeplan(t.traces);check(p.order+p.diagonal==(parent?10:11) && p.order<=7,"unexpected length or nondiagonal count");
  // At most 3^11 walks. Replacing traceless entries and expanding all diagonal powers gives l1 norm at most 2^11 per walk.
  // Nonzero-letter entry codes are >=8; zero padding identifies their number. Length fixes the total diagonal degree.
  Engine e;e.apply(t);check(e.bound<=Int(177147)*2048,"integer bound");vector<pair<Key,Int>> monomials;for(auto [k,v]:e.packed)if(v)monomials.emplace_back(k,v);sort(monomials.begin(),monomials.end());out<<monomials.size();for(auto [k,v]:monomials)out<<' '<<k<<' '<<v;out<<'\n';
 }
 string s;check(!(in>>s)&&out.good(),"trailing input or output failure");return 0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
