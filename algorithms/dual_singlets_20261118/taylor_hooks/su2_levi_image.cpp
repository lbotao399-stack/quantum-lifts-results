// Restrict an explicit SU(3) image to traceless 2 x 2 matrices and compare it with all classical boundaries.
#define N4_PARITY_HELPERS_ONLY
#include "unitary_parity.cpp"
#include <fstream>

struct TraceTerm {
  u64 coefficient;
  std::vector<size_t> factors;
};

struct RestrictedImage {
  std::vector<std::vector<int>> traces;
  std::map<std::string, std::vector<TraceTerm>> by_skeleton;
};

static RestrictedImage read_image(const char *path, const Sector &sector, int length) {
  std::ifstream in(path);
  size_t letters, traces, terms;
  if (!(in >> letters >> traces >> terms)) throw std::runtime_error("missing image header");
  std::vector<int> alphabet(letters);
  for (int &l : alphabet) {
    int type, flavour, plus, minus;
    in >> type >> flavour >> plus >> minus;
    l = ltid(type, flavour, plus, minus);
    if (l < 0) throw std::runtime_error("letter outside Taylor sector");
  }
  RestrictedImage image;
  for (size_t i = 0; i < traces; ++i) {
    size_t n;
    in >> n;
    if (n < 2 || n > static_cast<size_t>(length)) throw std::runtime_error("invalid trace length");
    std::vector<int> word(n);
    for (int &l : word) { size_t id; in >> id; l = alphabet.at(id); }
    image.traces.push_back(std::move(word));
  }
  for (size_t i = 0; i < terms; ++i) {
    long long signed_value;
    size_t n;
    in >> signed_value >> n;
    TraceTerm term{static_cast<u64>((signed_value % static_cast<long long>(P) + P) % P), {}};
    std::vector<std::pair<int, int>> entries;
    int degrees[5] = {};
    for (size_t j = 0; j < n; ++j) {
      size_t id;
      in >> id;
      term.factors.push_back(id);
      for (int l : image.traces.at(id)) {
        entries.emplace_back(l, 1);
        for (int k = 0; k < 5; ++k) degrees[k] += LT[l].w[k];
      }
    }
    if (!in || entries.size() != static_cast<size_t>(length) ||
        !std::equal(degrees, degrees+5, sector.WS)) throw std::runtime_error("incorrect term weight");
    image.by_skeleton[skel_of(entries)].push_back(std::move(term));
  }
  std::string extra;
  if (in >> extra) throw std::runtime_error("unexpected data after the image");
  return image;
}

static u64 evaluate_image(const RestrictedImage &image, const Skel &s, int probe) {
  auto found = image.by_skeleton.find(s.key);
  if (found == image.by_skeleton.end()) return 0;
  std::vector<GMat> letters(LT.size());
  for (size_t i = 0; i < letters.size(); ++i) letter_at(s, s.pr, probe, i, letters[i]);
  std::unordered_map<size_t, GSc> cache;
  u64 value = 0;
  for (const TraceTerm &term : found->second) {
    GSc product{{0, term.coefficient}}, next;
    for (size_t id : term.factors) {
      auto inserted = cache.try_emplace(id);
      GSc &trace = inserted.first->second;
      if (inserted.second) {
        const auto &word = image.traces[id];
        GMat left = letters[word.front()], tmp;
        for (size_t k = 1; k+1 < word.size() && !left.empty(); ++k) {
          gmul(left, letters[word[k]], tmp);
          left.swap(tmp);
        }
        gtrmul(left, letters[word.back()], trace);
      }
      gsmul(product, trace, next);
      product.swap(next);
      if (product.empty()) break;
    }
    for (const auto &item : product) if (item.first == s.top_te) value = addm(value, item.second);
  }
  return value;
}

int main(int argc, char **argv) {
  if (argc != 10 || !parse_gauge("SU2")) return 1;
  int weights[5];
  for (int i = 0; i < 5; ++i) weights[i] = std::stoi(argv[2+i]);
  int length = std::stoi(argv[7]);
  P = std::stoull(argv[8]);
  if (P != 2147483647 && P != 1000033 && P != 1000003) return 1;
  MERS = P == 2147483647;
  u64 seed = std::stoull(argv[9]);
  RS = seed;
  if (const char *threads = std::getenv("PAR_THREADS")) NTHREADS = std::max(1, std::stoi(threads));
  set_mm();
  set_bmm();
  Sector sector = make_sector(weights);
  init_all({sector});
  auto input = read_image(argv[1], sector, length);
  Lev parent, next;
  Cohom target;
  target.sec = sector;
  build_lev(sector, parent, length-1);
  build_lev(sector, target.L0, length);
  build_lev(sector, next, length+1);
  target.c0 = col_offsets(target.L0, target.nc);
  auto incoming = q0_rows(sector, parent, target.L0);
  auto outgoing = q0_rows(sector, target.L0, next);
  compute_res(target);
  auto stacked = outgoing;
  for (const auto &block : target.res) {
    for (const auto &evaluations : block.val) {
      std::vector<u32> row(target.nc, 0);
      std::copy(evaluations.begin(), evaluations.end(), row.begin()+target.c0[block.si]);
      stacked.push_back(std::move(row));
    }
  }
  std::vector<std::pair<size_t, int>> probes;
  for (size_t k = 0; k < target.L0.sk.size(); ++k)
    for (int i : target.L0.sk[k].piv) probes.emplace_back(k, i);
  std::vector<u32> values(probes.size());
  par_for(probes.size(), [&](size_t i) {
    const auto point = probes[i];
    values[i] = static_cast<u32>(evaluate_image(input, target.L0.sk[point.first], point.second));
  });
  print_head("SU3_image_on_SU2", weights, sector, seed);
  std::cout << ",\"length\":" << length << ",\"cochains\":[" << parent.dim << ',' << target.nc << ',' << next.dim << ']';
  std::cout << ",\"image_values_scaled\":[";
  for (size_t i = 0; i < values.size(); ++i) { if (i) std::cout << ','; std::cout << values[i]; }
  std::cout << "],\"incoming\":"; print_minor(incoming, parent.dim);
  std::cout << ",\"outgoing\":"; print_minor(std::move(outgoing), target.nc);
  std::cout << ",\"outgoing_and_Cartan\":"; print_minor(std::move(stacked), target.nc);
  for (size_t i = 0; i < probes.size(); ++i) incoming.at(i).push_back(values[i]);
  std::cout << ",\"augmented\":"; print_minor(std::move(incoming), parent.dim+1);
  std::cout << "}\n";
}
