// Contract graded trace templates after the signs have been extracted in Python.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using Wide = __int128_t;
struct Entry { int i, j, value; };
using Matrix = std::vector<Entry>;
using Word = std::vector<int>;
constexpr int dim = 26;
constexpr int colours = 52;
int dense[colours][dim][dim]{};
int64_t weight[colours]{};

Matrix multiply(const Matrix& left,const Matrix& right) {
    int entries[dim][dim]{};
    for (auto a:left) for (auto b:right) if (a.j==b.i)
        entries[a.i][b.j]+=a.value*b.value;
    Matrix out;
    for (int i=0;i<dim;i++) for (int j=0;j<dim;j++)
        if (entries[i][j]) out.push_back({i,j,entries[i][j]});
    return out;
}

struct Tables {
    std::array<std::vector<Matrix>,4> words;
    Tables() {
        std::ifstream in("rep");
        words[1].resize(colours);
        for (int i=0;i<colours;i++) {
            int count; in>>count;
            while (count--) {
                Entry e; in>>e.i>>e.j>>e.value;
                words[1][i].push_back(e); dense[i][e.i][e.j]=e.value;
            }
        }
        for (auto& x:weight) in>>x;
        assert(in);
        for (int k=2;k<=3;k++)
            for (const auto& a:words[k-1]) for (const auto& b:words[1])
                words[k].push_back(multiply(a,b));
    }
    const Matrix& half(const Word& slots,const unsigned char* key,int start,int n) const {
        int ix=0;
        for (int i=start;i<start+n;i++) ix=ix*colours+key[slots[i]];
        return words[n][ix];
    }
    int64_t trace(const Word& slots,const unsigned char* key) const {
        int n=static_cast<int>(slots.size());
        assert(2<=n && n<=7);
        int a=n/2, b=n==7 ? 3 : n-a;
        const auto& left=half(slots,key,0,a);
        if (left.empty()) return 0;
        const auto& right=half(slots,key,n-b,b);
        int64_t result=0;
        if (n==7) {
            const auto& middle=dense[key[slots[3]]];
            for (auto x:left) for (auto y:right) if (y.j==x.i)
                result+=static_cast<int64_t>(x.value)*middle[x.j][y.i]*y.value;
        } else {
            for (auto x:left) for (auto y:right) if (x.j==y.i && y.j==x.i)
                result+=static_cast<int64_t>(x.value)*y.value;
        }
        return result;
    }
};

std::string decimal(Wide x) {
    if (!x) return "0";
    bool negative=x<0;
    if (negative) x=-x;
    std::string s;
    while (x) { s.push_back(static_cast<char>('0'+x%10)); x/=10; }
    if (negative) s.push_back('-');
    std::reverse(s.begin(),s.end());
    return s;
}

int main() {
    Tables tables;
    std::ifstream query("query",std::ios::binary), patterns("patterns");
    int count,width,buckets,number;
    query>>count>>width>>buckets;
    assert(query.get()=='\n' && width<=9);
    patterns>>number;
    std::map<Word,int> word_ids;
    std::vector<Word> words;
    using Term = std::pair<int,Word>;
    std::map<Word,std::vector<Term>> groups;
    for (int j=0;j<number;j++) {
        int factors; patterns>>factors;
        Word masks,ids;
        for (int t=0;t<factors;t++) {
            int length; patterns>>length;
            Word w(length); int mask=0;
            for (auto& v:w) { patterns>>v; mask|=1<<v; }
            auto [it,inserted]=word_ids.emplace(w,static_cast<int>(words.size()));
            if (inserted) words.push_back(w);
            masks.push_back(mask); ids.push_back(it->second);
        }
        std::sort(masks.begin(),masks.end());
        groups[masks].push_back({j,ids});
    }
    assert(patterns);
    std::vector<Wide> sums(static_cast<size_t>(number)*buckets);
    std::vector<int> visited(words.size(),-1);
    std::vector<int64_t> traces(words.size());
    unsigned char key[11]{};
    std::array<int64_t,512> sums_of_weights{};
    for (int i=0;i<count;i++) {
        query.read(reinterpret_cast<char*>(key),width+2);
        int bucket=key[width]+256*key[width+1];
        assert(query && bucket<buckets);
        for (int j=0;j<width;j++) {
            assert(key[j]<colours);
            int bit=1<<j;
            for (int k=0;k<bit;k++) sums_of_weights[bit+k]=sums_of_weights[k]+weight[key[j]];
        }
        assert(sums_of_weights[(1<<width)-1]==0);
        Wide* total=sums.data()+static_cast<size_t>(bucket)*number;
        for (const auto& [masks,terms]:groups) {
            if (std::any_of(masks.begin(),masks.end(),[&](int m){return sums_of_weights[m]!=0;})) continue;
            for (const auto& [index,ids]:terms) {
                int64_t v=1;
                for (int id:ids) {
                    if (visited[id]!=i) {
                        traces[id]=tables.trace(words[id],key);
                        visited[id]=i;
                    }
                    v*=traces[id];
                    if (!v) break;
                }
                total[index]+=v;
            }
        }
    }
    assert(query.peek()==std::char_traits<char>::eof());
    std::ofstream out("values");
    for (int i=0;i<buckets;i++) {
        for (int j=0;j<number;j++) {
            Wide v=sums[static_cast<size_t>(i)*number+j];
            if (v) out<<j<<" "<<decimal(v)<<" ";
        }
        out<<'\n';
    }
    assert(out);
}
