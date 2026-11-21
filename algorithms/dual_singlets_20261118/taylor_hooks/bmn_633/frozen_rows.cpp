// Re-evaluate the signed-minor recipe directly on its archived integer SO(7) matrices.
#define N4_PARITY_HELPERS_ONLY
#include "../unitary_parity.cpp"
#include <fstream>

static std::vector<int> read_vector() {
  size_t n; std::cin >> n;
  if (n>16) throw std::runtime_error("invalid generator indices");
  std::vector<int> values(n);
  for (auto &value : values) std::cin >> value;
  return values;
}

static std::vector<std::vector<int>> read_lists() {
  size_t n; std::cin >> n;
  if (n>10) throw std::runtime_error("invalid number of generator slots");
  std::vector<std::vector<int>> values(n);
  for (auto &v : values) v = read_vector();
  return values;
}

static size_t skeleton(Lev &level, const std::string &key, int length) {
  for (size_t i=0; i<level.sk.size(); ++i) if (level.sk[i].key==key) return i;
  level.sk.push_back(parse_skel(key,length,0));
  return level.sk.size()-1;
}

int main(int argc,char **argv) {
  if (argc!=3 || !parse_gauge("SO7")) return 1;
  P=std::stoull(argv[1]);
  if (P!=2147483647 && P!=1000033 && P!=1000003) return 1;
  MERS=P==2147483647; RS=317083;
  if (const char *threads=std::getenv("PAR_THREADS")) NTHREADS=std::max(1,std::stoi(threads));
  set_mm(); set_bmm();
  int w[5]={6,3,3,0,0}; Sector sec=make_sector(w); init_all({sec});
  std::cin.exceptions(std::ios::failbit | std::ios::badbit);
  size_t nr,nc; std::cin>>nr>>nc;
  if (nr!=685 || nc!=686) return 2;
  Lev src,next;
  std::vector<std::pair<size_t,size_t>> columns,rows;
  for (size_t j=0; j<nc; ++j) {
    std::string key; std::cin>>key;
    size_t i=skeleton(src,key,10);
    auto &s=src.sk[i];
    Gen g; g.cols=read_lists(); g.edges=read_lists(); g.mu=read_vector(); g.tr=read_lists();
    if (g.cols.size()!=10 || g.edges.size()!=10 || g.mu.size()!=10) return 2;
    columns.emplace_back(i,s.gens.size()); s.gens.push_back(std::move(g)); ++src.dim;
  }
  Mtx values(nr,std::vector<u32>(nc));
  for (size_t k=0; k<nr; ++k) {
    std::string kind,key; size_t n; std::cin>>kind>>key>>n;
    bool classical=kind=="Q0";
    if (!classical && kind!="Cartan") return 2;
    size_t index=skeleton(classical ? next : src,key,classical ? 11 : 10);
    auto &s=classical ? next.sk[index] : src.sk[index];
    std::vector<std::vector<Mat>> point(LT.size());
    for (size_t i=0; i<n; ++i) {
      int type,flav,r,t; size_t copies; std::cin>>type>>flav>>r>>t>>copies;
      int l=ltid(type,flav,r,t);
      if (l<0 || !point[l].empty()) return 2;
      for (size_t v=0; v<copies; ++v) {
        Mat m; mzero(m);
        for (int j=0; j<NN; ++j) {
          long long x; std::cin>>x;
          m[j]=static_cast<u32>((x%static_cast<long long>(P)+P)%P);
          if (!classical && j/NF!=j%NF && m[j]) return 2;
        }
        if (!in_alg(m)) return 2;
        point[l].push_back(m);
      }
    }
    if (classical) {
      rows.emplace_back(index,s.pr.size());
      s.piv.push_back(static_cast<int>(s.pr.size())); s.pr.push_back(std::move(point)); ++next.dim;
    } else {
      rows.emplace_back(SIZE_MAX,0);
      Probes single{std::move(point)};
      par_for(nc,[&](size_t j) {
        auto [idx,pos]=columns[j];
        if (idx==index) values[k][j]=static_cast<u32>(eval_at(sec,s,s.gens[pos],single,0));
      });
    }
  }
  auto image=q0_rows(sec,src,next);
  size_t dummy;
  auto col_offset=col_offsets(src,dummy);
  std::vector<size_t> row_offset;
  size_t offset=0;
  for (const auto &s:next.sk) { row_offset.push_back(offset); offset+=s.piv.size(); }
  for (size_t k=0; k<nr; ++k) {
    auto [i,r]=rows[k];
    if (i==SIZE_MAX) continue;
    for (size_t j=0; j<nc; ++j) {
      auto [h,c]=columns[j];
      values[k][j]=image[row_offset[i]+r][col_offset[h]+c];
    }
  }
  std::ofstream os(argv[2],std::ios::binary);
  os.exceptions(std::ios::failbit | std::ios::badbit);
  for (const auto &v:values) os.write(reinterpret_cast<const char *>(v.data()),v.size()*sizeof(u32));
}
