// Integer order-ten Brauer contractions of adjoint letters at rank four.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifndef RANK_FOUR_LETTERS
#define RANK_FOUR_LETTERS 12
#endif
constexpr int letters = RANK_FOUR_LETTERS;
using Slots = std::array<int,2*letters>;
using Key = std::array<int,letters/2>;
struct Block {
    int length=0, code=0;
    std::array<int,2> odd{-1,-1};
};
struct Diagram {
    std::array<std::uint8_t,10> partner;
    int sign;
};

void match(std::vector<Diagram>& out, Diagram d, int used, int sign) {
    if (used==1023) {
        d.sign=sign;
        out.push_back(d);
        return;
    }
    int i=0;
    while (used & (1<<i)) ++i;
    for (int j=i+1;j<10;++j) if (!(used & (1<<j))) {
        d.partner[i]=j;
        d.partner[j]=i;
        match(out,d,used | (1<<i) | (1<<j),sign);
        sign=-sign;
    }
}

bool block_less(const Block& a,const Block& b) {
    if (a.length!=b.length) return a.length<b.length;
    for (int i=0;i<5;++i) {
        int x=(a.code>>(3*i))&7, y=(b.code>>(3*i))&7;
        if (x!=y) return x<y;
    }
    return false;
}

int main(int argc,char** argv) {
    assert(argc==3);
    bool sp=std::string(argv[2])=="Sp";
    assert(sp || std::string(argv[2])=="SO");
    std::ifstream table(argv[1]);
    int width;
    table>>width;
    assert(width>0);
    std::map<Key,int> ids;
    for (int i=0;i<width;++i) {
        int size;
        table>>size;
        assert(size<=letters/2);
        Key key{};
        for (int j=0;j<size;++j) table>>key[j];
        std::sort(key.begin(),key.end());
        assert(ids.emplace(key,i).second);
    }
    std::vector<Diagram> pf;
    if (sp) {match(pf,{},0,1);assert(pf.size()==945);}
    std::array<int,letters> species;
    Slots slots;
    while (std::cin>>species[0]) {
        for (int i=1;i<letters;++i) std::cin>>species[i];
        const int odd_count=std::count_if(species.begin(),species.end(),[](int s){return s<2;});
        for (auto& v:slots) std::cin>>v;
        Slots mate{},edge{};
        edge.fill(1);
        for (int k=sp?10:20;k<2*letters;k+=2) {
            int a=slots[k],b=slots[k+1];
            mate[a]=b;
            mate[b]=a;
            edge[b]=-1;
        }
        std::vector<int> raw(width);
        auto add=[&](int sign) {
            int visited=0,blocks=0;
            std::array<Block,letters/2> terms;
            for (int i=0;i<letters;++i) if (!(visited&(1<<i))) {
                int current=2*i;
                Block b;
                do {
                    int letter=current/2,s=species[letter];
                    visited|=1<<letter;
                    if (sp) sign*=edge[current^1];
                    else if (current%2) sign=-sign;
                    b.code+=1<<(3*s);
                    ++b.length;
                    if (s<2) {
                        if (b.odd[s]>=0) return;
                        b.odd[s]=letter;
                    }
                    current=mate[current^1];
                } while (current!=2*i);
                if (b.length%2) return;
                assert(blocks<letters/2);
                terms[blocks++]=b;
            }
            std::sort(terms.begin(),terms.begin()+blocks,block_less);
            Key key{};
            std::array<int,letters> odd{};
            int n=0;
            for (int i=0;i<blocks;++i) {
                const auto& b=terms[i];
                if (i && b.code==terms[i-1].code && ((b.odd[0]>=0) != (b.odd[1]>=0))) return;
                key[i]=b.code;
                for (int j:b.odd) if (j>=0) odd[n++]=j;
            }
            assert(n==odd_count);
            for (int j=0;j<6;++j) for (int i=0;i<j;++i) if (odd[i]>odd[j]) sign=-sign;
            std::sort(key.begin(),key.end());
            raw[ids.at(key)]+=sign;
        };
        if (sp) {
            for (const auto& d:pf) {
                for (int i=0;i<10;++i) {
                    mate[slots[i]]=slots[d.partner[i]];
                    edge[slots[i]]=i<d.partner[i]?1:-1;
                }
                add(d.sign);
            }
        } else {
            // Fix the diagonal, then use Heap transpositions. Each transposition reverses the sign.
            std::array<int,10> perm{},counts{};
            for (int i=0;i<10;++i) {
                perm[i]=i;
                mate[slots[i]]=slots[i+10];
                mate[slots[i+10]]=slots[i];
            }
            int sign=1, i=1;
            add(sign);
            while (i<10) {
                if (counts[i]==i) {counts[i]=0;++i;continue;}
                int j=(i%2)?counts[i]:0;
                std::swap(perm[i],perm[j]);
                for (int k:{i,j}) {
                    int a=slots[k],b=slots[10+perm[k]];
                    mate[a]=b;
                    mate[b]=a;
                }
                sign=-sign;
                add(sign);
                ++counts[i];
                i=1;
            }
        }
        for (int i=0;i<width;++i) if (raw[i]) std::cout<<i<<' '<<raw[i]<<' ';
        std::cout<<'\n';
    }
}
