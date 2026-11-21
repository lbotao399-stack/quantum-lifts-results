// Read a pinned lowering once. Split by (theta degree, total jet order), without reprojecting.
#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
void need(bool ok,const char* what){if(!ok)throw runtime_error(what);}
int main(int argc,char** argv){try{
 need(argc==5,"split_profiles alphabet binary outdir expected-mu-as-five-digits");string mu=argv[4];need(mu.size()==5,"weight");
 ifstream alphabet(argv[1]);array<array<int,3>,120> letters{};int a,b,t,n=0;while(alphabet>>a>>b>>t){need(n<119,"alphabet size");letters[++n]={a,b,t};}need(n==119&&alphabet.eof(),"alphabet");
 filesystem::path out=argv[3];filesystem::create_directories(out);ifstream input(argv[2],ios::binary);need(input.good(),"input");
 map<array<unsigned char,9>,ofstream> files;map<array<unsigned char,9>,uint64_t> counts;uint64_t total=0;
 while(true){int len=input.get();if(len<0){need(input.eof()&&!input.bad(),"header");break;}need(11<=len&&len<=13 || len==10,"record length");unsigned char data[32];input.read((char*)data,4+len);need(input.good()&&data[3+len]==0,"truncated");uint32_t c=data[0]|(data[1]<<8)|(data[2]<<16)|(uint32_t(data[3])<<24);need(c&&c<1000033,"coefficient");array<unsigned char,9> profile{};array<int,5> weight{};int pos=0,word=0;for(int i=0;i<len;i++){int l=data[4+i];if(!l){need(word>=2,"trace length");word=0;continue;}need(l<=119&&pos<9,"letter");word++;auto [x,y,mask]=letters[l];weight[0]+=x;weight[1]+=y;for(int j=0;j<3;j++)weight[2+j]+=(mask>>j)&1;profile[pos++]=5*__builtin_popcount(unsigned(mask))+x+y;}
 need(pos==9&&word==0,"nine letters");for(int j=0;j<5;j++)need(weight[j]==mu[j]-'0',"unexpected weight");sort(profile.begin(),profile.end());auto &count=counts[profile];if(!count){string name="p";for(int x:profile)name+="_"+to_string(x);auto &o=files[profile];o.open(out/(name+".bin"),ios::binary);need(o.good(),"open output");}auto &o=files.at(profile);o.put(char(len));o.write((char*)data,4+len);need(o.good(),"write");count++;total++;
 }
 for(auto &[p,f]:files){f.close();need(f.good(),"close output");}for(auto &[p,c]:counts){for(auto x:p)cout<<int(x)<<' ';cout<<c<<'\n';}cerr<<total<<" verified records in "<<counts.size()<<" profiles\n";return 0;
 }catch(const exception& e){cerr<<e.what()<<'\n';return 1;}}
