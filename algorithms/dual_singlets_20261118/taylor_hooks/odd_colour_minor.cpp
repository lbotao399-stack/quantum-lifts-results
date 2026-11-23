// Q0 on C-odd SU(3) cochains. Select the odd basis before forming the differential.
#define N4_PARITY_HELPERS_ONLY
#include "unitary_parity.cpp"
#include <stdexcept>

static Lev odd_level(const Sector &sector, int length) {
  Lev level;
  level.Y = length;
  for (const auto &entry : hw_by_skel(sector, length)) {
    Skel skeleton = parse_skel(entry.first, length, entry.second);
    skeleton.K = static_cast<int>(skeleton.hw) + MARGIN;
    add_probes(skeleton, skeleton.pr, skeleton.K, false);
    Lev mirror;
    mirror.sk.push_back(skeleton);
    conjugate_probes(mirror);
    IncEch even, odd;
    even.N = odd.N = skeleton.K;
    size_t stalled = 0;
    while (even.R.size() + odd.R.size() < static_cast<size_t>(skeleton.hw)) {
      if (++stalled > std::max<size_t>(50000, 200 * skeleton.hw))
        throw std::runtime_error("Parity basis stalled in " + skeleton.key);
      Gen gen;
      if (!random_gen(sector, skeleton, gen)) continue;
      std::vector<u64> direct, reversed;
      eval_all(sector, skeleton, gen, skeleton.pr, skeleton.K, direct);
      eval_all(sector, skeleton, gen, mirror.sk.front().pr, skeleton.K, reversed);
      for (int i = 0; i < skeleton.K; ++i) {
        u64 x = direct[i], y = reversed[i];
        direct[i] = addm(x, y);
        reversed[i] = subm(x, y);
      }
      bool new_even = even.insert(std::move(direct));
      bool new_odd = odd.insert(std::move(reversed));
      if (new_odd) skeleton.gens.push_back(std::move(gen));
      if (new_even || new_odd) stalled = 0;
    }
    if (even.R.size() + odd.R.size() != static_cast<size_t>(skeleton.hw))
      throw std::runtime_error("Parity ranks exceed the character dimension");
    level.exact += skeleton.hw;
    level.dim += odd.R.size();
    skeleton.piv = std::move(odd.piv);
    level.sk.push_back(std::move(skeleton));
  }
  std::cerr << "Y=" << length << " C-odd=" << level.dim << " full=" << level.exact << '\n';
  return level;
}

static Mtx odd_image(const Sector &sector, const Lev &source, Lev &target) {
  Mtx rows;
  rows.reserve(target.dim);
  for (auto &skeleton : target.sk) {
    if (skeleton.piv.empty()) continue;
    Lev block;
    block.Y = target.Y;
    block.sk.push_back(std::move(skeleton));
    Mtx direct = q0_rows(sector, source, block);
    conjugate_probes(block);
    Mtx reversed = q0_rows(sector, source, block);
    for (size_t i = 0; i < direct.size(); ++i) {
      for (size_t j = 0; j < source.dim; ++j)
        direct[i][j] = static_cast<u32>(subm(direct[i][j], reversed[i][j]));
      rows.push_back(std::move(direct[i]));
    }
    skeleton = std::move(block.sk.front());
  }
  return rows;
}

int main(int argc, char **argv) {
  if (argc != 9 || !parse_gauge("SU3")) return 1;
  int weights[5];
  for (int i = 0; i < 5; ++i) weights[i] = std::stoi(argv[i + 1]);
  if (weights[0] < weights[1] || weights[1] < weights[2] || weights[2] < 0 ||
      weights[3] < weights[4] || weights[4] < 0) return 1;
  int length = std::stoi(argv[6]);
  P = std::stoull(argv[7]);
  const u64 seed = std::stoull(argv[8]);
  if (!seed || (P != 1000003 && P != 1000033 && P != 2147483647)) return 1;
  RS = seed;
  MERS = P == 2147483647;
  if (const char *threads = std::getenv("PAR_THREADS")) NTHREADS = std::max(1, std::stoi(threads));
  set_mm();
  set_bmm();
  Sector sector = make_sector(weights);
  if (length < 1 || length >= sector.Ymax) return 1;
  init_all({sector});
  const Lev source = odd_level(sector, length);
  Lev target = odd_level(sector, length + 1);
  Mtx image = odd_image(sector, source, target);
  if (image.size() != target.dim) throw std::runtime_error("Incomplete target frame");
  print_head("C_odd_projected_classical_image", weights, sector, seed);
  std::cout << ",\"length\":" << length << ",\"odd_source\":" << source.dim
            << ",\"odd_target\":" << target.dim << ",\"full_source\":" << source.exact
            << ",\"full_target\":" << target.exact << ",\"minor\":";
  print_minor(std::move(image), source.dim, true);
  std::cout << "}\n";
}
