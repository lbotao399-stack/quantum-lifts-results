// Charge-conjugation-even BMN cochains in the long-root A2 of G2.
#define N4_PARITY_HELPERS_ONLY
#include "unitary_parity.cpp"
#include "bmn_sampler.h"
#include <fstream>

using Dimensions = std::map<std::pair<int, int>, int>;

static Lev even_level(const Sector &sector, int length, const Dimensions &dimensions,
                      bool probes_only) {
  Lev level;
  level.Y = length;
  for (const auto &entry : hw_by_skel(sector, length)) {
    Skel block = parse_skel(entry.first, length, entry.second);
    int fermions = 0;
    for (const auto &group : block.g)
      if (group.type == BET) fermions += group.k;
    auto found = dimensions.find({length, fermions});
    int count = found == dimensions.end() ? 0 : found->second;
    if (count > block.hw) throw std::runtime_error("even dimension exceeds both parities");
    if (!count) continue;
    int probes = count + MARGIN;
    block.K = probes;
    add_probes(block, block.pr, probes, false);
    if (probes_only) {
      block.piv.resize(probes);
      std::iota(block.piv.begin(), block.piv.end(), 0);
    } else {
      Probes point_pairs = block.pr;
      for (const auto &ordinary : block.pr) {
        auto reflected = ordinary;
        for (auto &letter : reflected)
          for (Mat &matrix : letter) {
            Mat transposed{};
            mtranspose(matrix, transposed);
            mneg(transposed);
            matrix = transposed;
          }
        point_pairs.push_back(std::move(reflected));
      }
      IncEch basis;
      basis.N = probes;
      size_t attempts = 0;
      while (static_cast<int>(basis.R.size()) < count) {
        if (++attempts > 100000u + 2000u * count)
          throw std::runtime_error("BMN sampling stalled before the even character dimension");
        Gen recipe;
        if (!bmn_gen(sector, block, recipe)) continue;
        std::vector<u64> values;
        eval_all(sector, block, recipe, point_pairs, 2*probes, values);
        for (int j = 0; j < probes; ++j) values[j] = addm(values[j], values[j+probes]);
        values.resize(probes);
        if (basis.insert(values)) block.gens.push_back(std::move(recipe));
      }
      block.piv = basis.piv;
    }
    level.dim += count;
    level.sk.push_back(std::move(block));
    std::cerr << "length " << length << ", fermions " << fermions
              << ", even dimension " << count << (probes_only ? ", probes" : ", basis") << std::endl;
  }
  size_t expected = 0;
  for (const auto &[label, value] : dimensions)
    if (label.first == length) expected += value;
  if (expected != level.dim) throw std::runtime_error("incomplete even character dimensions");
  return level;
}

static Mtx even_image(const Sector &sector, const Lev &source, Lev &target) {
  Mtx original = q0_rows(sector, source, target);
  conjugate_probes(target);
  Mtx reflected = q0_rows(sector, source, target);
  conjugate_probes(target);
  par_for(original.size(), [&](size_t i) {
    for (size_t j = 0; j < source.dim; ++j)
      original[i][j] = static_cast<u32>(addm(original[i][j], reflected[i][j]));
  });
  return original;
}

int main(int argc, char **argv) {
  if (argc != 9) {
    std::cerr << "usage: su3_even_bmn w1 w2 w3 length seed prime threads dimensions.txt\n";
    return 1;
  }
  int weights[5] = {std::stoi(argv[1]), std::stoi(argv[2]), std::stoi(argv[3]), 0, 0};
  int length = std::stoi(argv[4]);
  u64 seed = std::stoull(argv[5]);
  P = std::stoull(argv[6]);
  NTHREADS = std::stoi(argv[7]);
  if (weights[0]<weights[1] || weights[1]<weights[2] || weights[2]<0 ||
      length<2 || !seed || NTHREADS<1 ||
      (P!=2147483647 && P!=1000003 && P!=1000033)) return 2;
  Dimensions dimensions;
  std::ifstream input(argv[8]);
  int y, f, count;
  while (input >> y >> f >> count) {
    if (f<0 || count<=0 || !dimensions.emplace(std::make_pair(y, f),count).second) return 2;
  }
  if (!input.eof() || dimensions.empty() || !parse_gauge("SU3")) return 2;
  RS = seed;
  MERS = P == 2147483647;
  set_mm();
  set_bmm();
  Sector sector = make_sector(weights);
  if (length > sector.Ymax) return 2;
  init_all({sector});
  Lev previous = even_level(sector, length-1, dimensions, false);
  Lev centre = even_level(sector, length, dimensions, false);
  Lev following = even_level(sector, length+1, dimensions, true);
  print_head("even_BMN_minors", weights, sector, seed);
  std::cout << ",\"length\":" << length << ",\"cochains\":" << centre.dim
            << ",\"parent\":" << previous.dim << ",\"target_cochains\":" << following.dim
            << ",\"incoming\":";
  print_minor(even_image(sector, previous, centre), previous.dim, true);
  previous = Lev{};
  std::cout << ",\"outgoing\":";
  print_minor(even_image(sector, centre, following), centre.dim, true);
  std::cout << "}\n";
}
