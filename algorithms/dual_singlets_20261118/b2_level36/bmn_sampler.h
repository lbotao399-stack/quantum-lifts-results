// Included immediately after random_gen in trace.h. Keep the native evaluator,
// Young antisymmetrizers, finite-colour character dimensions and differential.
#pragma once

template <class T> static void bmn_shuffle(std::vector<T> &v) {
  for (size_t k = v.size(); k > 1; --k) std::swap(v[k-1], v[rnd() % k]);
}

static void bmn_shapes(int left, int minimum, std::vector<int> &part,
                       std::vector<std::vector<int>> &all) {
  if (left == 0) { all.push_back(part); return; }
  for (int k = minimum; k <= left; ++k) {
    part.push_back(k);
    bmn_shapes(left-k, k, part, all);
    part.pop_back();
  }
}

static bool bmn_gen(const Sector &sec, const Skel &S, Gen &G) {
  if (sec.n) return random_gen(sec, S, G);
  std::vector<int> scalars, fermions, strengths;
  for (int p = 0; p < S.Y; ++p) {
    switch (S.g[S.pg[p]].type) {
    case GAM: scalars.push_back(p); break;
    case BET: fermions.push_back(p); break;
    case BF: strengths.push_back(p); break;
    default: return random_gen(sec, S, G);
    }
  }
  G.cols.assign(S.Y, {}); G.edges.assign(S.Y, {}); G.mu.assign(S.Y, 0);
  int nf = int(fermions.size()), nb = int(strengths.size());
  // Every field strength transforms as the flavour determinant. Assign it a
  // full height-three column before distributing the remaining flavour boxes.
  for (int j = 0; j < nb; ++j) G.cols[strengths[j]] = {j,j,j};
  std::vector<std::array<int,3>> missing;
  for (int i = 0; i <= nf; ++i) for (int j = 0; j <= nf-i; ++j) {
    std::array<int,3> m = {i,j,nf-i-j};
    bool ok = true;
    for (int k = 0; k < 3; ++k) if (sec.WS[k]-nb-nf+m[k] < 0) ok = false;
    if (ok) missing.push_back(m);
  }
  if (missing.empty() || nb > sec.WS[2]) return false;
  const auto &m = missing[rnd() % missing.size()];
  std::vector<int> sf, ff;
  std::array<std::vector<int>,3> columns;
  for (int k = 0; k < 3; ++k) {
    sf.insert(sf.end(), sec.WS[k]-nb-nf+m[k], k);
    ff.insert(ff.end(), m[k], k);
    for (int j = nb; j < sec.WS[k]; ++j) columns[k].push_back(j);
    bmn_shuffle(columns[k]);
  }
  if (sf.size() != scalars.size() || ff.size() != fermions.size()) return false;
  bmn_shuffle(sf); bmn_shuffle(ff);
  auto take = [&](int p, int k) {
    G.cols[p].push_back(columns[k].back());
    columns[k].pop_back();
  };
  for (size_t i = 0; i < sf.size(); ++i) take(scalars[i], sf[i]);
  for (size_t i = 0; i < ff.size(); ++i)
    for (int k = 0; k < 3; ++k) if (k != ff[i]) take(fermions[i], k);

  // Give each admissible multitrace partition a turn, including products
  // of two-letter traces that are rare under uniform random permutations.
  static std::map<int,std::vector<std::vector<int>>> cache;
  auto &shapes = cache[S.Y];
  if (shapes.empty()) {
    std::vector<int> part; bmn_shapes(S.Y, 2, part, shapes);
    if (shapes.empty()) return false;
  }
  std::vector<int> order(S.Y);
  std::iota(order.begin(), order.end(), 0); bmn_shuffle(order);
  G.tr.clear();
  int offset = 0;
  for (int k : shapes[rnd() % shapes.size()]) {
    G.tr.emplace_back(order.begin()+offset, order.begin()+offset+k);
    offset += k;
  }
  return true;
}
