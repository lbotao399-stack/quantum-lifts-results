// Bound the classical cohomology killed by restriction to the Cartan.
#define N4_PARITY_HELPERS_ONLY
#include "unitary_parity.cpp"

int main(int argc, char **argv) {
  if ((argc != 9 && argc != 10) || !parse_gauge(argc == 10 ? argv[9] : "SU2")) return 1;
  int w[5];
  for (int i = 0; i < 5; ++i) w[i] = std::stoi(argv[i+1]);
  const int y = std::stoi(argv[6]);
  P = std::stoull(argv[7]);
  const u64 seed = std::stoull(argv[8]);
  if (P != 1000003 && P != 1000033 && P != 2147483647) return 1;
  if (w[0] < w[1] || w[1] < w[2] || w[3] < w[4] || w[2] < 0 || w[4] < 0) return 1;
  Sector sector = make_sector(w);
  if (y <= 0 || y > sector.Ymax) return 1;
  RS = seed;
  MERS = P == 2147483647;
  VERB = 0;
  if (const char *threads = std::getenv("PAR_THREADS")) NTHREADS = std::max(1, std::stoi(threads));
  set_mm();
  set_bmm();
  init_all({sector});
  Cohom target;
  target.sec = sector;
  target.Y = y;
  build_lev(sector, target.L0, y);
  target.c0 = col_offsets(target.L0, target.nc);
  print_head(argc == 10 ? "classical_Cartan_kernel" : "SU2_classical_Cartan_kernel", w, sector, seed);
  std::cout << ",\"length\":" << y << ",\"cochains\":" << target.nc;
  if (target.nc) {
    build_lev(sector, target.Lm, y-1);
    build_lev(sector, target.Lp, y+1);
    auto incoming = q0_rows(sector, target.Lm, target.L0);
    auto outgoing = q0_rows(sector, target.L0, target.Lp);
    std::cout << ",\"parent\":" << target.Lm.dim << ",\"next\":" << target.Lp.dim;
    std::cout << ",\"incoming\":";
    print_minor(std::move(incoming), target.Lm.dim, true);
    std::cout << ",\"outgoing\":";
    print_minor(outgoing, target.nc, true);
    compute_res(target);
    for (const auto &block : target.res) {
      for (const auto &evaluation : block.val) {
        std::vector<u32> row(target.nc, 0);
        std::copy(evaluation.begin(), evaluation.end(), row.begin() + target.c0[block.si]);
        outgoing.push_back(std::move(row));
      }
    }
    std::cout << ",\"stacked_rows\":" << outgoing.size() << ",\"outgoing_and_Cartan\":";
    print_minor(std::move(outgoing), target.nc, true);
  }
  std::cout << "}\n";
}
