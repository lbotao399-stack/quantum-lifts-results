// Homology representatives in generators, duals in uncontracted evaluation coordinates.
#define FILTERED_QUANTUM
#include "g2_filtered_sector.cpp"
#include <iostream>

static Mtx read_matrix(const fs::path &p) {
  std::ifstream in(p,std::ios::binary);
  in.exceptions(std::ios::failbit|std::ios::badbit);
  u32 h[4]; in.read(reinterpret_cast<char *>(h),sizeof(h));
  if(h[0]!=0x47325130 || (P && P!=h[1])) throw std::logic_error("matrix header or prime");
  P=h[1]; MERS=P==2147483647;
  Mtx rows(h[2],std::vector<u32>(h[3]));
  for(auto &r:rows) in.read(reinterpret_cast<char *>(r.data()),r.size()*sizeof(u32));
  return rows;
}

static std::vector<size_t> sparse_order(const Mtx &rows) {
  std::vector<size_t> order(rows.size()),size(rows.size());
  std::iota(order.begin(),order.end(),0);
  for(size_t i:order) size[i]=std::count_if(rows[i].begin(),rows[i].end(),[](u32 c){return c;});
  std::stable_sort(order.begin(),order.end(),[&](size_t i,size_t j){return size[i]<size[j];});
  return order;
}

static std::vector<u64> wide(const std::vector<u32> &r) {return {r.begin(),r.end()};}

int main(int argc,char **argv) {
  if(argc!=5) return 2;
  fs::path dir=argv[1]; auto grade=std::string(argv[2]); int y=std::stoi(argv[3]);
  NTHREADS=std::stoi(argv[4]); P=0;
  auto name=[&](const std::string &prefix,int yy) {
    return dir/(prefix+"_"+grade+"_"+std::to_string(yy)+".u32");
  };
  Mtx in=read_matrix(name("in",y)),eval=read_matrix(name("eval",y)),out=read_matrix(name("in",y+1));
  size_t n=eval.size();
  Mtx gen=transpose(eval,n),acted=transpose(out,n);
  if (in.size()!=n || gen.size()!=n || acted.size()!=n) throw std::logic_error("sizes");
  IncEch e; e.N=n;
  for(const auto &v:transpose(in,in.at(0).size())) e.insert(wide(v));
  size_t boundaries=e.R.size();
  std::vector<size_t> quotient;
  for(auto j:sparse_order(acted)) if(e.insert(wide(gen[j]))) quotient.push_back(j);
  if(e.R.size()!=n) throw std::logic_error("singular evaluation");
  Mtx compact(out.size(),std::vector<u32>(quotient.size()));
  for(size_t i=0;i<out.size();i++) for(size_t j=0;j<quotient.size();j++)
    compact[i][j]=out[i][quotient[j]];
  auto cycles=kernel_of(compact,quotient.size());
  auto piv=rref(compact,quotient.size());
  Mtx h;
  for(auto &v:cycles) {
    std::vector<u32> full(n);
    for(size_t j=0;j<quotient.size();j++) full[quotient[j]]=v[j];
    h.push_back(std::move(full));
  }
  auto hv=mul_abt(eval,h,n);
  for(size_t i=0;i<h.size();i++)
    std::cerr<<"cycle "<<i<<" support "<<std::count_if(h[i].begin(),h[i].end(),[](u32 x){return x;})<<"\n";
  if(std::any_of(h.begin(),h.end(),[&](const auto &row){
       auto v=mul_abt(out,Mtx{row},n);return std::any_of(v.begin(),v.end(),[](const auto &a){return a[0];});}))
    throw std::logic_error("cycle");
  store_matrix(name("cycle",y),h,n);
  std::ofstream cut(name("pivots",y));
  for(auto c:piv) cut<<quotient[c]<<"\n";
  if(grade=="2") {
    auto candidates=kernel_of(transpose(in,in.at(0).size()),n);
    auto values=transpose(hv,h.size());
    IncEch test; test.N=h.size();
    Mtx dual;
    for(auto j:sparse_order(candidates)) {
      auto row=transpose(mul_abt(values,Mtx{candidates[j]},n),1)[0];
      if(test.insert(wide(row))) dual.push_back(candidates[j]);
    }
    if(dual.size()!=h.size()) throw std::logic_error("dual dimension");
    store_matrix(name("dual",y),dual,n);
  }
  std::cout<<"{\"grade\":"<<grade<<",\"length\":"<<y<<",\"dimension\":"<<n
           <<",\"incoming\":"<<boundaries<<",\"outgoing\":"<<piv.size()
           <<",\"cohomology\":"<<h.size()<<"}\n";
}
