// Evaluate Q0 at a complete target frame and its charge-conjugate frame.
#define main archived_engine_main
#include "engine.cpp"
#undef main
#include <iostream>
#include <numeric>

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

static std::pair<Mtx, Mtx> parity_images(const Sector &sector,
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
  return {std::move(original), std::move(conjugated)};
}

static std::pair<size_t, size_t> image_ranks(const Sector &sector,
                                            const Lev &source, Lev &target) {
  auto images = parity_images(sector, source, target);
  return {rank_of(std::move(images.first), source.dim),
          rank_of(std::move(images.second), source.dim)};
}

static void print_indices(const std::vector<size_t> &indices) {
  std::cout << '[';
  for (size_t i = 0; i < indices.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << indices[i];
  }
  std::cout << ']';
}

// The determinant refers to the selected rows in pivot order and columns in
// increasing order. Row swaps affect this recorded order, hence no extra sign.
static void print_minor(Mtx matrix, size_t width, bool exploit_zeros = false) {
  std::vector<size_t> ids(matrix.size()), rows, cols;
  std::iota(ids.begin(), ids.end(), 0);
  u64 determinant = 1;
  for (size_t j = 0; j < width && rows.size() < matrix.size(); ++j) {
    size_t i = rows.size(), k = i;
    while (k < matrix.size() && !matrix[k][j]) ++k;
    if (k == matrix.size()) continue;
    std::swap(matrix[k], matrix[i]);
    std::swap(ids[k], ids[i]);
    rows.push_back(ids[i]);
    cols.push_back(j);
    const u64 pivot = matrix[i][j], inverse = invm(pivot);
    determinant = mulm(determinant, pivot);
    std::vector<size_t> nonzero;
    if (exploit_zeros) {
      for (size_t c = j; c < width; ++c) if (matrix[i][c]) nonzero.push_back(c);
    }
    const bool sparse = exploit_zeros && nonzero.size()*2 < width-j;
    par_for(matrix.size() - i - 1, [&](size_t offset) {
      auto &row = matrix[i+1+offset];
      if (!row[j]) return;
      const u64 scale = negm(mulm(row[j], inverse));
      if (sparse) {
        for (size_t c : nonzero) row[c] = static_cast<u32>(red(row[c] + scale * matrix[i][c]));
      } else {
        for (size_t c = j; c < width; ++c) row[c] = static_cast<u32>(red(row[c] + scale * matrix[i][c]));
      }
    });
    std::vector<u32>().swap(matrix[i]);
  }
  std::cout << "{\"rank\":" << rows.size() << ",\"det\":" << determinant << ",\"rows\":";
  print_indices(rows);
  std::cout << ",\"columns\":";
  print_indices(cols);
  std::cout << '}';
}

#ifndef N4_PARITY_HELPERS_ONLY
int main(int argc, char **argv) {
  const bool image_only = argc == 11 && std::string(argv[10]) == "image";
  const bool complete_image = argc == 11 && std::string(argv[10]) == "complete-image";
  if (argc != 10 && !image_only && !complete_image) {
    std::cerr << "usage: unitary_parity SU2|SU3 theta1 theta2 theta3 plus minus length seed prime [image|complete-image]\n";
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
  if (image_only || complete_image) {
    build_lev(sector, centre, length);
    following.Y = length + 1;
    if (complete_image) {
      build_lev(sector, following, length + 1);
    } else if (length < sector.Ymax) {
      for (const auto &entry : hw_by_skel(sector, following.Y)) {
        auto skeleton = parse_skel(entry.first, following.Y, entry.second);
        int count = static_cast<int>(std::min<i64>(entry.second, centre.dim)) + 3;
        add_probes(skeleton, skeleton.pr, count, false);
        skeleton.piv.resize(count);
        std::iota(skeleton.piv.begin(), skeleton.piv.end(), 0);
        following.sk.push_back(std::move(skeleton));
      }
    }
    auto images = parity_images(sector, centre, following);
    print_head(complete_image ? "complete_colour_image" : "colour_image", weights, sector, seed);
    std::cout << ",\"length\":" << length << ",\"cochains\":" << centre.dim
              << ",\"target_evaluations\":" << images.first.size() << ",\"even\":";
    print_minor(std::move(images.first), centre.dim);
    std::cout << ",\"odd\":";
    print_minor(std::move(images.second), centre.dim);
    std::cout << "}\n";
    return 0;
  }
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
#endif
