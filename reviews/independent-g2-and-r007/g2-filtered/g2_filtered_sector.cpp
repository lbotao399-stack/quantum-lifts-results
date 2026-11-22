// Export classical maps in a fixed generator/evaluation basis for one G2 charge family.
#define main classical_reference_main
#include "g2_dotted_chain.cpp"
#undef main
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <sstream>

namespace fs = std::filesystem;
using Weight = std::array<int,5>;

static void store_matrix(const fs::path &path, const Mtx &rows, size_t width) {
  std::ofstream f(path.string()+".part", std::ios::binary);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  u32 header[4] = {0x47325130, (u32)P, (u32)rows.size(), (u32)width};
  f.write(reinterpret_cast<const char *>(header),sizeof(header));
  for (const auto &r : rows) {
    if (r.size()!=width) throw std::logic_error("matrix width");
    f.write(reinterpret_cast<const char *>(r.data()),r.size()*sizeof(u32));
  }
  f.close();
  fs::rename(path.string()+".part",path);
}

static void initialize_exceptional(const std::vector<Sector> &family) {
  NF=7; NN=49; DIMG=14; RANK=2; GNAME="G2";
  set_mm(); set_bmm(); set_g2_basis();
  TB=g2_basis;
  std::vector<std::vector<u64>> gram(DIMG,std::vector<u64>(DIMG));
  for (int i=0;i<DIMG;i++) for (int j=0;j<DIMG;j++) gram[i][j]=mtrmul(TB[i],TB[j]);
  auto inverse=inv_matrix(gram);
  for (int i=0;i<DIMG;i++) {
    Mat dual{};
    for (int j=0;j<DIMG;j++) {
      Mat term=TB[j]; mscale(term,inverse[i][j]); maddto(dual,term);
    }
    TD.push_back(dual);
    for (int j=0;j<DIMG;j++) if (mtrmul(dual,TB[j])!=u64(i==j))
      throw std::logic_error("trace dual");
  }
  POSR={{1,0,0},{0,1,0},{1,1,0},{1,-1,0},{2,1,0},{1,2,0}};
  ADJW={{0,0,0},{0,0,0}};
  for (auto r:POSR) { ADJW.push_back(r); ADJW.push_back({-r[0],-r[1],0}); }
  int max_y=0;
  std::fill(BOX,BOX+5,0);
  for (const auto &s:family) {
    max_y=std::max(max_y,s.Ymax);
    for (int j=0;j<5;j++) BOX[j]=std::max(BOX[j],s.WS[j]);
  }
  build_char_grid(max_y+1); build_letters(); build_q(); build_raise();
  KEEPRAW=true; HWCHECK=true;
}

static void fixed_level(const Sector &s, int length, int base, Lev &l) {
  RS=base+97*length;
  for(int c:s.WS) RS=11*RS+c;
  build_lev(s,l,length);
}

static Mtx evaluation(const Lev &l) {
  Mtx out(l.dim,std::vector<u32>(l.dim));
  size_t start=0;
  for (const auto &s:l.sk) {
    for (size_t i=0;i<s.piv.size();i++) for (size_t j=0;j<s.gens.size();j++)
      out[start+i][start+j]=s.raw[j][s.piv[i]];
    start+=s.piv.size();
  }
  if (start!=l.dim) throw std::logic_error("basis dimension");
  return out;
}

static std::vector<Sector> charge_family() {
  Weight source{5,2,2,1,1};
  if (const char *setting=std::getenv("G2_SOURCE_WEIGHT")) {
    std::istringstream input(setting);
    for (int &v:source) if (!(input>>v)) throw std::invalid_argument("five source degrees required");
    std::string extra;
    if (input>>extra || source[0]<source[1] || source[1]<source[2] || source[2]<2 ||
        source[3]<source[4] || source[4]<0) throw std::invalid_argument("invalid source degrees");
  }
  std::vector<Sector> v;
  for(int order=-1;order<=2;order++) {
    // There is no incoming positive-order sector when the second dotted degree is zero.
    int shift=(order==-1 && source[4]==0) ? 0 : order;
    int w[]={source[0]-shift,source[1]-shift,source[2]-shift,source[3]+shift,source[4]+shift};
    v.push_back(make_sector(w));
  }
  return v;
}

#ifndef FILTERED_QUANTUM
int main(int argc,char **argv) {
  if(argc!=9) return 2;
  fs::path dir=argv[1];
  int grade=std::stoi(argv[2]),lo=std::stoi(argv[3]),hi=std::stoi(argv[4]);
  P=std::stoul(argv[5]); MERS=P==2147483647;
  if(!MERS && P!=1000033 && P!=1000003) return 2;
  int base=std::stoi(argv[6]); NTHREADS=std::stoi(argv[7]);
  bool ends_only=std::stoi(argv[8]);
  auto family=charge_family();
  initialize_exceptional(family);
  const auto &s=family.at(grade+1);
  fs::create_directories(dir);
  Lev prev;
  for(int y=lo;y<=hi;y++) {
    Lev curr;
    fixed_level(s,y,base,curr);
    auto name=std::to_string(grade)+"_"+std::to_string(y)+".u32";
    if(!ends_only || (lo<y && y<hi)) store_matrix(dir/("eval_"+name),evaluation(curr),curr.dim);
    if(y>lo) store_matrix(dir/("in_"+name),q0_rows(s,prev,curr),prev.dim);
    prev=std::move(curr);
  }
}
#endif
