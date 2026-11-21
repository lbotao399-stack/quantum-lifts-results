// A single classical image between complete character-sized highest-weight frames.
#define N4_PARITY_HELPERS_ONLY
#include "unitary_parity.cpp"

int main(int argc, char **argv) {
  if (argc != 10 || !parse_gauge(argv[1])) return 1;
  int w[5];
  for (int i = 0; i < 5; ++i) w[i] = std::stoi(argv[2+i]);
  const int length = std::stoi(argv[7]);
  P = std::stoull(argv[8]);
  if (P != 2147483647 && P != 1000033 && P != 1000003) return 1;
  MERS = P == 2147483647;
  const u64 seed = std::stoull(argv[9]);
  RS = seed;
  if (const char *threads = std::getenv("PAR_THREADS")) NTHREADS = std::max(1, std::stoi(threads));
  set_mm();
  set_bmm();
  Sector sector = make_sector(w);
  if (length < 1 || length >= sector.Ymax) return 1;
  init_all({sector});
  Lev source, target;
  build_lev(sector, source, length);
  build_lev(sector, target, length+1);
  Mtx matrix = q0_rows(sector, source, target);
  print_head("complete_classical_image", w, sector, seed);
  std::cout << ",\"length\":" << length << ",\"dimensions\":[" << source.dim << ',' << target.dim << "],\"minor\":";
  print_minor(std::move(matrix), source.dim, true);
  std::cout << "}\n";
}
