// Compute a lower bound on rank Q0 out of the shortest BMN cochains.
// The codomain needs evaluations of images, with no choice of a target basis.
#define main archived_engine_main
#include "engine.cpp"
#undef main

int main(int argc, char **argv) {
  if (argc != 5 || !parse_gauge("Sp2")) return 2;
  int w[5] = {10,4,4,0,0};
  const int length = 10;
  const int probes = atoi(argv[1]);
  if (probes < 1) return 2;
  RS = std::stoull(argv[2]);
  NTHREADS = std::max(1, atoi(argv[3]));
  const std::string selection(argv[4]);
  set_mm(); set_bmm();
  Sector sec = make_sector(w);
  init_all({sec});
  Lev src, dst;
  src.Y = length; dst.Y = length+1;
  for (const auto &kv : hw_by_skel(sec, src.Y)) {
    if (selection != "all" && selection != kv.first) continue;
    Skel s = parse_skel(kv.first, src.Y, kv.second);
    build_skel_basis(sec,s);
    fprintf(stderr, "source %s: %zu/%lld after %zu generators\n",
            kv.first.c_str(), s.gens.size(), (long long)kv.second, s.tried);
    src.dim += s.gens.size();
    src.sk.push_back(std::move(s));
  }
  if (!src.dim) return 3;
  for (const auto &kv : hw_by_skel(sec, dst.Y)) {
    Skel s = parse_skel(kv.first, dst.Y, kv.second);
    const int count = int(std::min<i64>({kv.second+3, i64(probes), i64(src.dim+3)}));
    bool reachable = false;
    for (const auto &a : src.sk) {
      for (size_t i = 0; i < a.g.size(); ++i) {
        std::vector<std::pair<int,int>> added;
        reachable |= merge_groups(a, int(i), s, added);
      }
    }
    if (!reachable) continue;
    add_probes(s,s.pr,count,false);
    s.piv.resize(count); std::iota(s.piv.begin(),s.piv.end(),0);
    dst.sk.push_back(std::move(s));
  }
  Mtx rows = q0_rows(sec,src,dst);
  size_t nrow = rows.size(), rank = rank_of(rows,src.dim);
  printf("{\"weight\":[10,4,4],\"length\":10,\"skeleton\":\"%s\","
         "\"source_dimension\":%zu,\"image_probes\":%zu,\"rank_lower_bound\":%zu,"
         "\"prime\":%llu,\"seed\":%s,\"source_blocks\":[",
         selection.c_str(),src.dim,nrow,rank,(unsigned long long)P,argv[2]);
  for (size_t i=0; i<src.sk.size(); ++i)
    printf("%s[\"%s\",%zu,%zu]", i ? "," : "",src.sk[i].key.c_str(),
           src.sk[i].gens.size(),src.sk[i].tried);
  printf("]}\n");
  return 0;
}
