// Produce the three coordinate maps used to construct the SO(7) rational Cartan-kernel cycle.
#define N4_PARITY_HELPERS_ONLY
#include "../unitary_parity.cpp"
#include <filesystem>
#include <fstream>

static void write_matrix(const std::filesystem::path &file, const Mtx &rows) {
  std::ofstream stream(file, std::ios::binary);
  stream.exceptions(std::ios::failbit | std::ios::badbit);
  const uint32_t sizes[] = {static_cast<uint32_t>(rows.size()), static_cast<uint32_t>(rows.at(0).size())};
  stream.write(reinterpret_cast<const char *>(sizes), sizeof sizes);
  for (const auto &row : rows) {
    if (row.size() != sizes[1]) throw std::runtime_error("ragged matrix");
    stream.write(reinterpret_cast<const char *>(row.data()), row.size()*sizeof(u32));
  }
}

template<class Values>
static void numbers(std::ostream &o, const Values &v) {
  o << '[';
  for (size_t j = 0; j < v.size(); ++j) { if (j) o << ','; o << v[j]; }
  o << ']';
}

static void lists(std::ostream &o, const std::vector<std::vector<int>> &v) {
  o << '[';
  for (size_t j = 0; j < v.size(); ++j) { if (j) o << ','; numbers(o, v[j]); }
  o << ']';
}

static void probes(std::ostream &o, const Skel &s, const Probes &points,
                   const std::vector<int> &ids, bool diagonal) {
  for (int k : ids) {
    o << "{\"kind\":\"" << (diagonal ? "Cartan" : "Q0") << "\",\"skeleton\":\"" << s.key << "\",\"letters\":[";
    bool first = true;
    for (size_t i = 0; i < points[k].size(); ++i) {
      if (points[k][i].empty()) continue;
      if (!first) o << ',';
      first = false;
      o << '[' << LT[i].type << ',' << LT[i].flav << ',' << LT[i].r << ',' << LT[i].s << ",[";
      for (size_t j = 0; j < points[k][i].size(); ++j) {
        if (j) o << ',';
        const Mat &v = points[k][i][j];
        o << '[';
        for (int t = 0; t < NN; ++t) {
          if (t) o << ',';
          long long x = v[t];
          if (x > static_cast<long long>(P/2)) x -= P;
          o << x;
        }
        o << ']';
      }
      o << "]]";
    }
    o << "]}\n";
  }
}

int main(int argc, char **argv) {
  if (argc != 2 || !parse_gauge("SO7")) return 1;
  const std::filesystem::path out(argv[1]);
  std::filesystem::create_directories(out);
  P = 2147483647; MERS = true; RS = 317083;
  if (const char *threads = std::getenv("PAR_THREADS")) NTHREADS = std::max(1, std::stoi(threads));
  set_mm(); set_bmm();
  int weights[5] = {6,3,3,0,0};
  const Sector sector = make_sector(weights);
  init_all({sector});
  Lev previous, next;
  Cohom target;
  target.sec = sector;
  build_lev(sector, previous, 9);
  KEEPRAW = true; build_lev(sector, target.L0, 10); KEEPRAW = false;
  build_lev(sector, next, 11);
  target.c0 = col_offsets(target.L0, target.nc);
  if (previous.dim!=1863 || target.nc!=1837 || next.dim!=772) return 2;
  write_matrix(out/"incoming.bin", q0_rows(sector, previous, target.L0));
  Mtx stacked = q0_rows(sector, target.L0, next);
  Mtx evaluation(target.nc, std::vector<u32>(target.nc));
  std::ofstream words(out/"generators.jsonl"), points(out/"probes.jsonl");
  words.exceptions(std::ios::badbit | std::ios::failbit);
  points.exceptions(std::ios::badbit | std::ios::failbit);
  for (const auto &s : next.sk) probes(points,s,s.pr,s.piv,false);
  compute_res(target);
  for (const auto &block : target.res) {
    const auto &s = target.L0.sk[block.si];
    size_t offset = target.c0[block.si];
    for (size_t j = 0; j < s.gens.size(); ++j) {
      const auto &g = s.gens[j];
      words << "{\"skeleton\":\"" << s.key << "\",\"cols\":"; lists(words,g.cols);
      words << ",\"edges\":"; lists(words,g.edges);
      words << ",\"mu\":"; numbers(words,g.mu);
      words << ",\"tr\":"; lists(words,g.tr); words << "}\n";
      for (size_t k = 0; k < s.piv.size(); ++k) evaluation[offset+k][offset+j] = s.raw[j][s.piv[k]];
    }
    for (const auto &evaluation_row : block.val) {
      std::vector<u32> row(target.nc);
      std::copy(evaluation_row.begin(), evaluation_row.end(), row.begin()+offset);
      stacked.push_back(std::move(row));
    }
    std::vector<int> ids(block.val.size());
    std::iota(ids.begin(), ids.end(), 0);
    probes(points,s,s.dg,ids,true);
  }
  write_matrix(out/"cochains.bin",evaluation);
  write_matrix(out/"closure_Cartan.bin",stacked);
}
