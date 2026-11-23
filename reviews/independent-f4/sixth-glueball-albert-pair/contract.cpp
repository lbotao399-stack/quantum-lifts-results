// The two open words are paired through the 26 in End(26).
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

constexpr int D=26, N=52;
using ColourWord=std::array<int,5>;
struct MatrixEntry {int row,col,value;};
struct ProductCoefficient {int colour,value;};
using Sparse=std::vector<MatrixEntry>;
using Vector=std::array<int,D>;
struct Quartic {Vector evaluated{}; std::vector<ProductCoefficient> dual;};
using Bracket=std::vector<ProductCoefficient>;
int prime;
int mod(int64_t k){int r=k%prime;return r<0?r+prime:r;}

template<class T>
void read_array(const char* filename,T& result){
    std::ifstream in(filename,std::ios::binary);
    in.read(reinterpret_cast<char*>(&result),sizeof(result));
    assert(in && in.peek()==std::char_traits<char>::eof());
}

Sparse product(const Sparse& a,const Sparse& b){
    int dense[D][D]{};
    for(const auto& x:a) for(const auto& y:b)
        if(x.col==y.row) dense[x.row][y.col]+=x.value*y.value;
    Sparse result;
    for(int i=0;i<D;++i)for(int j=0;j<D;++j)
        if(dense[i][j]) result.push_back({i,j,dense[i][j]});
    return result;
}

int rep[N][D][D], probes[D][D][D], metric[D][D], brackets[N][N][N];
Sparse colours[N],pairs[N][N];
Bracket algebra[N][N];
std::vector<Quartic> quartics;
int lookup[N][N][N][N];

Vector moments(const Sparse& m){
    Vector v{};
    for(int t=0;t<D;++t){
        int64_t sum=0;
        for(auto x:m)sum+=int64_t(x.value)*probes[t][x.row][x.col];
        v[t]=mod(sum);
    }
    return v;
}

void prepare(){
    read_array("representation.bin",rep);read_array("quartic.bin",probes);
    read_array("metric.bin",metric);read_array("brackets.bin",brackets);
    for(int i=0;i<N;++i){
        for(int a=0;a<D;++a)for(int b=0;b<D;++b)
            if(rep[i][a][b])colours[i].push_back({a,b,rep[i][a][b]});
        for(int j=0;j<N;++j)for(int k=0;k<N;++k)
            if(brackets[i][j][k])algebra[i][j].push_back({k,brackets[i][j][k]});
    }
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)pairs[i][j]=product(colours[i],colours[j]);
    for(int i=0;i<N;++i)for(int j=0;j<N;++j)for(int k=0;k<N;++k)
        assert(moments(product(pairs[i][j],colours[k]))==Vector{});
    for(int i=0;i<N;++i)for(int j=i;j<N;++j)for(int k=j;k<N;++k)for(int l=k;l<N;++l){
        Quartic q;
        q.evaluated=moments(product(pairs[i][j],pairs[k][l]));
        for(int t=0;t<D;++t){
            int sum=0;
            for(int u=0;u<D;++u)sum+=q.evaluated[u]*metric[u][t];
            if(mod(sum))q.dual.push_back({t,mod(sum)});
        }
        assert(q.dual.size()<=2);
        std::array<int,4> a{i,j,k,l};
        do { lookup[a[0]][a[1]][a[2]][a[3]]=quartics.size(); }
        while(std::next_permutation(a.begin(),a.end()));
        quartics.push_back(std::move(q));
    }
}

struct Family {
    std::array<int,4> left{};
    ColourWord rest{};
    std::vector<std::array<int,10>> orientation;
    std::vector<int> indices;
};

int main(int argc,char** argv){
    assert(argc==2);prime=std::stoi(argv[1]);assert(prime==1009 || prime==1013);
    prepare();
    std::vector<Family> families;
    std::map<std::array<int,4>,int> lookup_family;
    std::ifstream p("patterns.txt");
    std::array<int,4> left;ColourWord order;
    int number=0;
    while(p>>left[0]>>left[1]>>left[2]>>left[3]){
        for(auto& i:order)p>>i;
        assert(p && std::is_sorted(left.begin(),left.end()));
        if(!lookup_family.count(left)){
            lookup_family[left]=families.size();
            Family f; f.left=left;f.rest=order;std::sort(f.rest.begin(),f.rest.end());
            families.push_back(f);
        }
        auto& f=families[lookup_family[left]];
        std::array<int,10> orient{};int pair=0;
        for(int i=0;i<5;++i)for(int j=i+1;j<5;++j){
            auto a=std::find(order.begin(),order.end(),f.rest[i]);
            auto b=std::find(order.begin(),order.end(),f.rest[j]);
            assert(a!=order.end() && b!=order.end());
            orient[pair++]=a<b?1:-1;
        }
        f.indices.push_back(number++);f.orientation.push_back(orient);
    }
    std::vector<int64_t> totals(number);
    std::ifstream query("parents.bin",std::ios::binary);
    std::array<unsigned char,11> key;
    int queries=0,controls=0;
    while(query.read(reinterpret_cast<char*>(key.data()),11)){
        int scalar=key[9]+256*key[10];
        for(const auto& f:families){
            auto a=f.left;
            const auto& dual=quartics[lookup[key[a[0]]][key[a[1]]][key[a[2]]][key[a[3]]]].dual;
            if(dual.empty())continue;
            ColourWord w;
            for(int i=0;i<5;++i)w[i]=key[f.rest[i]];
            std::array<int,10> entries{};int pair=0;
            for(int i=0;i<5;++i)for(int j=i+1;j<5;++j){
                int rest[3],n=0;
                for(int t=0;t<5;++t)if(t!=i && t!=j)rest[n++]=w[t];
                int64_t v=0;
                for(auto b:algebra[w[i]][w[j]]){
                    const auto& q=quartics[lookup[b.colour][rest[0]][rest[1]][rest[2]]].evaluated;
                    for(auto x:dual)v+=int64_t(b.value)*x.value*q[x.colour];
                }
                entries[pair++]=mod(v);
            }
            if(controls<25000 && queries%13==0){
                auto m=product(product(pairs[w[0]][w[1]],pairs[w[2]][w[3]]),colours[w[4]]);
                auto v=moments(m);int got=0,want=0;
                for(auto x:dual)got+=2*x.value*v[x.colour];
                for(auto k:entries)want+=k;
                assert(mod(got-want)==0);++controls;
            }
            for(size_t h=0;h<f.indices.size();++h){
                int sum=0;
                for(int t=0;t<10;++t)sum+=f.orientation[h][t]*entries[t];
                totals[f.indices[h]]+=int64_t(scalar)*sum;
            }
        }
        ++queries;
    }
    assert(query.eof() && query.gcount()==0 && queries==2747440 && controls==25000);
    for(auto v:totals)std::cout<<mod(int64_t(mod(v))*((prime+1)/2))<<"\n";
    std::cerr<<queries<<" colour assignments, "<<controls<<" direct five-letter checks, "<<families.size()<<" left sets\n";
}
