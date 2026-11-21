// Reuse the exact cyclic canonicalization and extremal projector of the parent.
#include "su_full_projection_native_v1.cpp"
#include <sstream>
#include <limits>
// Codes 0,1,2: S; 3,4: T_{+,0}, T_{+,1}. Store the source Taylor convention.
std::array<std::vector<std::pair<int,int>>,5> deriv;
std::vector<int> swapped;
void require(bool ok,const char* s){if(!ok)throw std::runtime_error(s);}
void add(Expression& e,const std::string& key,int64_t c){if(!(c=norm(c)))return;auto& v=e[key];v=norm(v+c);if(!v)e.erase(key);}
void plus(Expression& a,const Expression& b,int64_t c){for(const auto& [k,v]:b)add(a,k,c*v);}
void setup(const std::vector<int>& table){int sector[5]={};const uint64_t offset=0;const unsigned char blank=0;int64_t zero=0;require(full_projection(1000033,table.size()/3,table.data(),sector,0,&offset,&blank,&zero)==0,"alphabet");
    std::map<std::array<int,3>,int> ids;for(unsigned i=1;i<projection_letters.size();++i)ids.emplace(projection_letters[i],i);
    for(auto& v:deriv)v.assign(projection_letters.size(),{0,0});swapped.assign(projection_letters.size(),0);
    for(const auto& [l,id]:ids){auto [u,v,mask]=l;swapped[id]=ids.at({v,u,mask});int n=__builtin_popcount((unsigned)mask);
        for(int op=0;op<5;++op){int t=op<3?op:op-3, bit=1<<t;
            int sign=(n+__builtin_popcount((unsigned)(mask&(bit-1))))%2 ? -1:1;
            if(op<3){int f=u+v+n-1;if((mask&bit)&&f)deriv[op][id]={ids.at({u,v,mask^bit}),f*sign};}
            else if(u>0 && !(mask&bit))deriv[op][id]={ids.at({u-1,v,mask|bit}),sign};
        }
    }
}
Expression apply_odd(const Expression& in,int op){Expression out;for(const auto& [k,c]:in){auto w=decode_projection_word(k);int before=1;
    for(unsigned i=0;i<w.size();++i)for(unsigned j=0;j<w[i].size();++j){int l=w[i][j];auto [to,coef]=deriv[op][l];if(coef){auto m=w;m[i][j]=to;int s=canonical(m);if(s)add(out,encode(m),s*before*coef*c);}if(odd[l])before=-before;}}
    return out;}
Expression spin_exchange(const Expression& in){Expression out;for(const auto& [k,c]:in){auto w=decode_projection_word(k);for(auto& tr:w)for(int& l:tr)l=swapped[l];int s=canonical(w);require(s!=0,"spin exchange lost word");add(out,encode(w),s*c);}return out;}
void highest(const Expression& x){for(int op:{0,2,4,6})require(transformed(x,op).empty(),"raising residual");}
Expression project_parent(Expression x){x=root_project(x,0,1,0,1);x=root_project(x,2,3,3,2);x=root_project(x,4,5,3,1);x=root_project(x,6,7,0,1);highest(x);return x;}
std::pair<Expression,Expression> images(const Expression& x,bool check){const int64_t inv=inverse_projection(3);require(inv==666689,"1/3");auto f1=transformed(x,5);auto s=apply_odd(x,0);plus(s,apply_odd(f1,1),inv);plus(s,apply_odd(transformed(f1,1),2),inv);
    auto t=apply_odd(x,4);plus(t,apply_odd(f1,3),-inv);t=spin_exchange(t);
    // Each positive root vanishes, so the final extremal projection is the identity.
    if(check){highest(s);highest(t);}return {std::move(s),std::move(t)};}
int main(int argc,char** argv){try{require(argc==5,"lower alphabet input output-prefix mode(parent|raw)");std::ifstream alph(argv[1]);std::vector<int> tab;int value;while(alph>>value)tab.push_back(value);require(tab.size()==119*3,"incomplete alphabet");setup(tab);
    std::ifstream input(argv[2],std::ios::binary);require(input.good(),"input");std::map<std::vector<int>,Expression> shapes;while(true){int n=input.get();if(n<0){require(input.eof(),"read");break;}require(n>=3&&n<=48,"key length");unsigned char b[4];input.read((char*)b,4);uint32_t c=b[0]|(b[1]<<8)|(b[2]<<16)|(uint32_t(b[3])<<24);std::string k(n,0);input.read(k.data(),n);require(input.good() && c>0 && c<(unsigned)prime && k.back()==0,"record");auto m=decode_projection_word(k);auto canonicalm=m;require(canonical(canonicalm)==1 && canonicalm==m,"uncanonical input");std::vector<int> shp;for(const auto& w:m)shp.push_back(w.size());std::sort(shp.begin(),shp.end());add(shapes[shp],k,c);}
    std::string prefix=argv[3],mode=argv[4];require(mode=="parent"||mode=="raw","mode");std::ofstream so(prefix+".S20.bin",std::ios::binary),tp(prefix+".Tp21.bin",std::ios::binary);require(so.good()&&tp.good(),"open output");
    uint64_t ns=0,nt=0,nx=0;for(auto& [shp,inputx]:shapes){auto x=mode=="parent"?project_parent(std::move(inputx)):std::move(inputx);nx+=x.size();auto [s,t]=images(x,mode=="parent");ns+=s.size();nt+=t.size();
        for(auto p:{std::make_pair(&s,&so),std::make_pair(&t,&tp)}){std::map<std::string,int64_t> ordered(p.first->begin(),p.first->end());for(const auto& [k,c]:ordered){p.second->put((char)k.size());for(int j=0;j<4;++j)p.second->put(char((c>>(8*j))&255));p.second->write(k.data(),k.size());}require(p.second->good(),"output");}}
    so.close();tp.close();require(so.good()&&tp.good(),"incomplete output");std::cout<<nx<<" "<<ns<<" "<<nt<<"\n";return 0;}catch(std::exception& e){std::cerr<<e.what()<<"\n";return 1;}}
