// The defining-seven trace metric in both fundamental-index routings.
using SparseMap=std::vector<std::vector<std::pair<int,u64>>>;
static SparseMap internal_colour, external_colour;
static void set_colour_maps() {
  internal_colour.resize(NN); external_colour.resize(NN);
  for(int i=0;i<NF;i++) for(int j=0;j<NF;j++)
    for(int k=0;k<NF;k++) for(int l=0;l<NF;l++) {
      const int a=i*NF+j,b=k*NF+l;
      u64 c=0;
      for(int h=0;h<DIMG;h++) c=addm(c,mulm(TB[h][a],TD[h][b]));
      if(!c) continue;
      internal_colour[i*NF+l].push_back({j*NF+k,c});
      external_colour[l*NF+k].push_back({i*NF+j,c});
    }
}
static void map_colour(const GMat &in,const SparseMap &map,GMat &out) {
  out.clear();
  for(const auto &x:in) {
    GE mapped{}; mapped.m=x.m;
    bool any=false;
    for(int i=0;i<NN;i++) {
      u64 c=0;
      for(auto [j,coef]:map[i]) c=addm(c,mulm(x.a[j],coef));
      mapped.a[i]=c; any=any || c;
    }
    if(any) out.push_back(mapped);
  }
}
static void psi3(const GMat &u,const GMat &v,const GMat &w,GSc &out) {
  GMat projected,uv; map_colour(v,internal_colour,projected);
  gmul(u,projected,uv); gtrmul(uv,w,out);
}
static void phi2(const GMat &a,const GSc &mid,const GMat &b,GSc &out) {
  GMat projected,last;
  map_colour(a,external_colour,projected);
  sgmul(mid,b,last); gtrmul(projected,last,out);
}
