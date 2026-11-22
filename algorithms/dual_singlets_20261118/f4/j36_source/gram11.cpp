// Exact descending-root Ward recurrence for the nine (1,1) endpoint covectors.
// Compile with signed overflow trapping. All large accumulators are 128-bit integers.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using Key = std::uint64_t;
using Int = std::int64_t;
using Wide = __int128_t;
using Weight = std::array<int,4>;
using Values = std::array<Int,9>;
using Sparse = std::vector<std::pair<Key,int>>;

struct Row {
    Values n{};
    Int d = 1;
};
Int narrow(Wide n) {
    assert(n >= std::numeric_limits<Int>::min() && n <= std::numeric_limits<Int>::max());
    return static_cast<Int>(n);
}
Wide gcd_wide(Wide a, Wide b) {
    if (a<0) a=-a;
    while (b) { Wide r=a%b; a=b; b=r; }
    return a;
}
std::string decimal(Wide n) {
    if (n<0) return "-"+decimal(-n);
    if (!n) return "0";
    std::string s;
    while(n) { s.push_back('0'+n%10); n/=10; }
    std::reverse(s.begin(),s.end());
    return s;
}
template<typename T> Key encode(const T& a) {
    Key k=1;
    for (int i:a) k=(k<<6)|i;
    return k;
}
std::vector<int> decode(Key k) {
    std::vector<int> a;
    while(k>1) { a.push_back(k&63); k>>=6; }
    std::reverse(a.begin(),a.end());
    return a;
}
int orientation(const std::vector<int>& a) {
    int s=1;
    for (unsigned i=0;i<a.size();++i) for(unsigned j=0;j<i;++j) {
        if (a[i]==a[j]) return 0;
        if (a[j]>a[i]) s=-s;
    }
    return s;
}

struct Recurrence {
    std::array<Weight,52> weights;
    std::array<Weight,4> twice_metric;
    Weight twice_rho{}, highest_weight{};
    Key lowest;
    int scale;
    std::vector<std::array<int,3>> roots;
    Sparse brackets[52][52], dual[52];
    std::unordered_map<Key,Sparse> actions, rotated;
    std::unordered_map<Key,Row> seeds, states;
    std::unordered_map<Key,int> gaps;
    std::map<Weight,std::vector<Key>> triples;
    std::vector<std::pair<Key,Key>> pairs;

    Weight weight(Key k) const {
        Weight w{};
        while(k>1) {
            for(int i=0;i<4;++i) w[i]+=weights[k&63][i];
            k>>=6;
        }
        return w;
    }
    void read() {
        for(auto& w:weights) for(int& v:w) std::cin>>v;
        for(auto& w:twice_metric) for(int& v:w) std::cin>>v;
        std::array<int,5> h,l;
        for(int& v:h) std::cin>>v;
        for(int& v:l) std::cin>>v;
        highest_weight=weight(encode(h)); lowest=encode(l);
        int count;
        std::cin>>scale>>count;
        roots.resize(count);
        for(auto& r:roots) {
            for(int& v:r) std::cin>>v;
            for(int j=0;j<4;++j) twice_rho[j]+=weights[r[0]][j];
        }
        auto entries=[](Sparse& row,int n) {
            for(int k=0;k<n;++k) { int i,v; std::cin>>i>>v; row.emplace_back(i,v); }
        };
        std::cin>>count;
        for(int k=0;k<count;++k) {
            int i,j,n; std::cin>>i>>j>>n; entries(brackets[i][j],n);
        }
        for(auto& row:dual) { int n; std::cin>>n; entries(row,n); }
        std::cin>>count;
        for(int k=0;k<count;++k) {
            std::array<int,3> t; for(int& v:t) std::cin>>v;
            auto& row=seeds[encode(t)];
            for(Int& v: row.n) std::cin>>v;
        }
        assert(std::cin && seeds.size()==236 && roots.size()==24);
    }
    const Sparse& move(Key k,int root) {
        Key id=(k<<6)|root;
        auto it=actions.find(id);
        if(it!=actions.end()) return it->second;
        auto a=decode(k);
        Sparse terms;
        for(unsigned j=0;j<a.size();++j) for(auto [i,v]:brackets[root][a[j]]) {
            auto b=a; b[j]=static_cast<int>(i);
            int s=orientation(b);
            if(s) { std::sort(b.begin(),b.end()); terms.emplace_back(encode(b),s*v); }
        }
        return actions.emplace(id,std::move(terms)).first->second;
    }
    int gap(Key k) {
        auto it=gaps.find(k);
        if(it!=gaps.end()) return it->second;
        auto w=weight(k);
        int d=180;
        for(int i=0;i<4;++i) for(int j=0;j<4;++j)
            d-=w[i]*twice_metric[i][j]*(w[j]-twice_rho[j]);
        gaps.emplace(k,d);
        return d;
    }
    const Row& low(Key v,Key t) {
        static const Row zero{};
        if(v==lowest) {
            auto it=seeds.find(t);
            return it==seeds.end()?zero:it->second;
        }
        Key k=(v<<19)|t;
        auto it=states.find(k);
        if(it!=states.end()) return it->second;
        int divisor=gap(v);
        if(divisor<=0) return zero;
        std::array<Wide,9> total{};
        Wide denominator=1;
        for(auto [e,f,s]:roots) {
            const auto& left=move(v,f);
            if(left.empty()) continue;
            const auto& right=move(t,e);
            for(auto [a,x]:left) for(auto [b,y]:right) {
                const auto& r=low(a,b);
                Wide m=r.d/gcd_wide(denominator,r.d);
                Wide factor=denominator*m/r.d*s*x*y;
                for(unsigned i=0;i<9;++i) total[i]=total[i]*m+factor*r.n[i];
                denominator*=m;
            }
        }
        denominator*=divisor;
        Wide common=denominator;
        for(auto& c:total) { c*=-2; common=gcd_wide(c,common); }
        Row out;
        out.d=narrow(denominator/common);
        for(unsigned i=0;i<9;++i) out.n[i]=narrow(total[i]/common);
        return states.emplace(k,out).first->second;
    }
    const Sparse& rotation(Key k) {
        auto it=rotated.find(k);
        if(it!=rotated.end()) return it->second;
        auto a=decode(k);
        std::vector<int> b(a.size());
        std::map<Key,int> combined;
        auto expand=[&](auto&& self,unsigned p,int c)->void {
            if(p==a.size()) {
                int s=orientation(b);
                if(s) { auto t=b; std::sort(t.begin(),t.end()); combined[encode(t)]+=s*c; }
                return;
            }
            for(auto [i,d]:dual[a[p]]) { b[p]=static_cast<int>(i); self(self,p+1,c*d); }
        };
        expand(expand,0,1);
        Sparse out;
        for(auto [v,c]:combined) if(c) out.emplace_back(v,c);
        return rotated.emplace(k,std::move(out)).first->second;
    }
    void enumerate() {
        for(int a=0;a<52;++a) for(int b=a+1;b<52;++b) for(int c=b+1;c<52;++c) {
            Key t=encode(std::array<int,3>{a,b,c});
            triples[weight(t)].push_back(t);
        }
        for(int a=0;a<52;++a) for(int b=a+1;b<52;++b) for(int c=b+1;c<52;++c)
            for(int d=c+1;d<52;++d) for(int e=d+1;e<52;++e) {
                Key v=encode(std::array<int,5>{a,b,c,d,e});
                auto w=weight(v);
                for(int j=0;j<4;++j) w[j]=-w[j]-highest_weight[j];
                auto it=triples.find(w);
                if(it!=triples.end()) for(Key t:it->second) pairs.emplace_back(v,t);
            }
        assert(pairs.size()==573459);
    }
    void gram() {
        Int common=1;
        std::size_t nonzero=0;
        for(auto [v,t]:pairs) {
            const auto& x=low(v,t);
            common=narrow(Wide(common/std::gcd(common,x.d))*x.d);
            nonzero+=(x.n!=Values{});
        }
        actions.clear();
        std::array<std::array<Wide,9>,9> sums{};
        for(auto [v,t]:pairs) {
            const auto& x=low(v,t);
            if(x.n==Values{}) continue;
            std::array<Wide,9> y{};
            for(auto [vv,s]:rotation(v)) for(auto [tt,r]:rotation(t)) {
                const auto& z=low(vv,tt);
                assert(common%z.d==0);
                Wide factor=Wide(s)*r*(common/z.d);
                for(unsigned j=0;j<9;++j) y[j]+=factor*z.n[j];
            }
            for(unsigned i=0;i<9;++i) for(unsigned j=0;j<9;++j)
                sums[i][j]+=Wide(x.n[i])*(common/x.d)*y[j];
        }
        std::cout<<pairs.size()<<" "<<nonzero<<" "<<common<<" "<<scale<<"\n";
        for(unsigned i=0;i<9;++i) {
            for(unsigned j=0;j<9;++j) {
                assert(sums[i][j]==sums[j][i]);
                std::cout<<decimal(sums[i][j])<<" ";
            }
            std::cout<<"\n";
        }
    }
};

int main() {
    Recurrence calculation;
    calculation.read();
    calculation.enumerate();
    calculation.gram();
}
