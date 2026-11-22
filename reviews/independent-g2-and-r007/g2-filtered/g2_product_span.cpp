// Match rational graviton products and all incoming Q0 images in identical G2 coordinates.
#define main native_classical_main
#include "g2_dotted_chain.cpp"
#undef main
#include <fstream>

struct ProductInput {
  std::vector<std::vector<int>> traces, words;
  std::map<std::string, std::vector<int>> words_by_skeleton;
  std::vector<std::vector<std::pair<int,u64>>> columns;
  int number_of_products = 0;
};

static u64 read_residue(std::istream &in) {
  std::string value;
  in >> value;
  u64 result = 0;
  for (char digit : value) {
    if (digit == '-') continue;
    if (digit < '0' || digit > '9') throw std::logic_error("invalid integer");
    result = red(10 * result + static_cast<unsigned>(digit - '0'));
  }
  return value.at(0) == '-' ? negm(result) : result;
}

static ProductInput read_products(const char *path) {
  std::ifstream in(path);
  in.exceptions(std::ios::failbit | std::ios::badbit);
  int nt, nw, np;
  in >> nt >> nw >> np;
  ProductInput packet;
  packet.traces.resize(nt);
  packet.words.resize(nw);
  packet.columns.resize(nw);
  packet.number_of_products = np;
  for (auto &trace : packet.traces) {
    int n;
    in >> n;
    for (int i = 0; i < n; ++i) {
      int t, f, r, s;
      in >> t >> f >> r >> s;
      int l = ltid(t, f, r, s);
      if (l < 0) throw std::logic_error("letter outside sector");
      trace.push_back(l);
    }
  }
  for (int w = 0; w < nw; ++w) {
    int n;
    in >> n;
    std::vector<std::pair<int,int>> content;
    for (int i = 0; i < n; ++i) {
      int t;
      in >> t;
      packet.words[w].push_back(t);
      for (int l : packet.traces.at(t)) content.emplace_back(l,1);
    }
    packet.words_by_skeleton[skel_of(content)].push_back(w);
  }
  for (int j = 0; j < np; ++j) {
    int n;
    in >> n;
    for (int i = 0; i < n; ++i) {
      int w;
      in >> w;
      u64 a = read_residue(in), b = read_residue(in);
      if (!b) throw std::logic_error("denominator vanishes at this prime");
      if (a) packet.columns.at(w).emplace_back(j,mulm(a,invm(b)));
    }
  }
  return packet;
}

static GSc trace_value(const std::vector<int> &trace, const std::vector<GMat> &letters) {
  GMat partial = letters.at(trace.at(0)), next;
  GSc out;
  if (trace.size() == 1) {
    gtr(partial,out);
  } else {
    for (size_t i = 1; i+1 < trace.size() && !partial.empty(); ++i) {
      gmul(partial,letters[trace[i]],next);
      partial.swap(next);
    }
    gtrmul(partial,letters[trace.back()],out);
  }
  return out;
}

static void append_products(Mtx &rows, const Lev &level, const ProductInput &packet) {
  size_t start = 0, previous = rows.at(0).size();
  for (const Skel &s : level.sk) {
    auto found = packet.words_by_skeleton.find(s.key);
    par_for(s.piv.size(),[&](size_t k) {
      auto &row = rows.at(start+k);
      row.resize(previous+packet.number_of_products,0);
      if (found == packet.words_by_skeleton.end()) return;
      std::vector<GMat> letters(LT.size());
      for (size_t l = 0; l < LT.size(); ++l) letter_at(s,s.pr,s.piv[k],l,letters[l]);
      std::vector<GSc> traces(packet.traces.size());
      std::vector<bool> ready(traces.size());
      for (int w : found->second) {
        GSc partial = {{0,1}}, next;
        for (int t : packet.words[w]) {
          if (!ready[t]) {
            traces[t]=trace_value(packet.traces[t],letters);
            ready[t]=true;
          }
          gsmul(partial,traces[t],next);
          partial.swap(next);
          if (partial.empty()) break;
        }
        for (auto [mask,v] : partial) if (mask == s.top_te)
          for (auto [column,c] : packet.columns[w])
            row[previous+column] = addm(row[previous+column],mulm(c,v));
      }
    });
    start += s.piv.size();
  }
  if (start != rows.size()) throw std::logic_error("probe rows do not match");
}

int main(int argc, char **argv) {
  if (argc != 11 && argc != 12) return 1;
  int w[5];
  for (int i = 0; i < 5; ++i) w[i]=std::stoi(argv[i+1]);
  int y=std::stoi(argv[6]);
  RS=std::stoull(argv[7]);
  u64 seed=RS;
  NTHREADS=std::stoi(argv[8]);
  P=std::stoull(argv[9]);
  if ((P!=1000003 && P!=1000033 && P!=2147483647) || NTHREADS<1) return 2;
  MERS=P==2147483647;
  if (argc==12) {
    if (!parse_gauge(argv[11])) return 2;
    set_mm(); set_bmm(); build_algebra();
  } else {
    NF=7; NN=49; RANK=2; DIMG=14; GNAME="G2";
    set_mm(); set_bmm(); set_g2_basis();
    POSR={{1,0,0},{0,1,0},{1,1,0},{1,-1,0},{2,1,0},{1,2,0}};
    ADJW={{0,0,0},{0,0,0}};
    for (auto r:POSR) { ADJW.push_back(r); ADJW.push_back({-r[0],-r[1],0}); }
  }
  Sector sec=make_sector(w);
  std::copy(w,w+5,BOX);
  build_char_grid(sec.Ymax+1); build_letters(); build_q(); build_raise();
  HWCHECK=true;
  ProductInput packet=read_products(argv[10]);
  Lev previous,current;
  build_lev(sec,previous,y-1);
  build_lev(sec,current,y);
  Mtx rows=q0_rows(sec,previous,current);
  append_products(rows,current,packet);
  auto pivot=rref(rows,previous.dim+packet.number_of_products);
  std::vector<int> products;
  for(int column:pivot) if (column>=static_cast<int>(previous.dim))
    products.push_back(column-static_cast<int>(previous.dim));
  printf("{\"seed\":%llu,\"prime\":%llu,\"weight\":[%d,%d,%d,%d,%d],"
         "\"length\":%d,\"dimension\":%zu,\"boundary\":%zu,\"augmented\":%zu,"
         "\"products\":%d,\"selected\":[",
         (unsigned long long)seed,(unsigned long long)P,w[3],w[4],w[0],w[1],w[2],
         y,current.dim,pivot.size()-products.size(),pivot.size(),packet.number_of_products);
  for(size_t i=0;i<products.size();++i) printf("%s%d",i?",":"",products[i]);
  puts("]}");
}
