// Evaluate Q0 at a complete target frame and its charge-conjugate frame.
#define main archived_engine_main
#include "engine.cpp"
#undef main
#include <iostream>

static void conjugate_probes(Lev &level) {
  for (Skel &skeleton : level.sk) {
    for (auto &probe : skeleton.pr) {
      for (auto &letter : probe) {
        for (Mat &matrix : letter) {
          Mat transposed{};
          mtranspose(matrix, transposed);
          mneg(transposed);
          matrix = transposed;
        }
      }
    }
  }
}

static std::pair<size_t, size_t> image_ranks(const Sector &sector,
                                            const Lev &source, Lev &target) {
  Mtx original = q0_rows(sector, source, target);
  conjugate_probes(target);
  Mtx conjugated = q0_rows(sector, source, target);
  conjugate_probes(target);
  par_for(original.size(), [&](size_t i) {
    for (size_t j = 0; j < source.dim; ++j) {
      u32 value = original[i][j];
      u32 mirror = conjugated[i][j];
      original[i][j] = static_cast<u32>(addm(value, mirror));
      conjugated[i][j] = static_cast<u32>(subm(value, mirror));
    }
  });
  size_t even = rank_of(std::move(original), source.dim);
  size_t odd = rank_of(std::move(conjugated), source.dim);
  return {even, odd};
}

int main(int argc, char **argv) {
  if (argc != 10) {
    std::cerr << "usage: unitary_parity SU2|SU3 theta1 theta2 theta3 plus minus length seed prime\n";
    return 1;
  }
  if (!parse_gauge(argv[1]) || GT != G_SU) return 1;
  int weights[5];
  for (int i = 0; i < 5; ++i) weights[i] = std::stoi(argv[2 + i]);
  if (weights[0] < weights[1] || weights[1] < weights[2] || weights[2] < 0 ||
      weights[3] < weights[4] || weights[4] < 0) return 1;
  int length = std::stoi(argv[7]);
  u64 seed = std::stoull(argv[8]);
  P = std::stoull(argv[9]);
  if (P != 2147483647 && P != 1000003 && P != 1000033) return 1;
  MERS = (P == 2147483647);
  if (seed == 0) return 1;
  RS = seed;
  if (const char *threads = std::getenv("PAR_THREADS"))
    NTHREADS = std::max(1, std::stoi(threads));
  set_mm();
  set_bmm();
  Sector sector = make_sector(weights);
  if (length < 1 || length > sector.Ymax) return 1;
  init_all({sector});
  Lev previous, centre, following;
  build_lev(sector, previous, length - 1);
  build_lev(sector, centre, length);
  build_lev(sector, following, length + 1);
  auto incoming = image_ranks(sector, previous, centre);
  auto outgoing = image_ranks(sector, centre, following);
  print_head("charge_conjugation", weights, sector, seed);
  std::cout << ",\"length\":" << length << ",\"cochains\":" << centre.dim
            << ",\"even_in\":" << incoming.first << ",\"odd_in\":" << incoming.second
            << ",\"even_out\":" << outgoing.first << ",\"odd_out\":" << outgoing.second
            << "}\n";
  return 0;
}
