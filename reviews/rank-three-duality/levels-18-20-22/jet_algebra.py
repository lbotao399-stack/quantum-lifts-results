"""Divided-jet Taylor generators in the native sign convention."""
from collections import Counter
from functools import cache
from itertools import product

def E(i,j): return ("E",i,j)
def P(i): return ("P",i)
def K(i): return ("K",i)
def parity(op): return sum(i>=2 for i in op[1:])%2
def sign(a,b): return -1 if parity(a)*parity(b) else 1


@cache
def bracket(a,b):
    phase = sign(a,b)
    i,j = a[1],b[1]
    if a[0]==b[0]=="E":
        m,n = a[2],b[2]
        return tuple((c,v) for c,v in ((int(m==j),E(i,n)),(-phase*int(n==i),E(j,m))) if c)
    if a[0]=="E" and b[0]=="P": return ((1,P(i)),) if a[2]==j else ()
    if a[0]=="E" and b[0]=="K": return ((-phase,K(a[2])),) if i==j else ()
    if a[0]=="K" and b[0]=="P":
        return ((phase,E(j,i)),) + (tuple((1,E(h,h)) for h in range(5)) if i==j else ())
    if a[0]==b[0]: return ()
    return tuple((-phase*c,op) for c,op in bracket(b,a))


def raising(op):
    return op[0]=="P" or op[0]=="E" and op[1]<2<=op[2]


def transpose(op):
    if op[0]=="E": return E(op[2],op[1])
    return (K if op[0]=="P" else P)(op[1])


def action(letter, op):
    dots,fs=letter[:2],letter[2]
    degree=sum(dots)+fs.bit_count()
    def derivative(i):
        if i<2:
            return (int(dots[i]>0),(*[n-int(h==i) for h,n in enumerate(dots)],fs))
        bit=1<<(i-2)
        return (int(bool(fs&bit))*(-1)**((fs//(2*bit)).bit_count()),(*dots,fs&~bit))
    def multiply(letter,j):
        p,q,f=letter
        if j<2:
            return ((p,q)[j]+1,(p+int(j==0),q+int(j==1),f))
        bit=1<<(j-2)
        return (int(not f&bit)*(-1)**((f//(2*bit)).bit_count()),(p,q,f|bit))
    if op[0]=="P": return multiply(letter,op[1])
    if op[0]=="K":
        c,x=derivative(op[1])
        return c*(degree-1),x
    c,x=derivative(op[2])
    if not c: return 0,letter
    n,x=multiply(x,op[1])
    return c*n,x


def verify_brackets():
    operators=[E(i,j) for i,j in product(range(5),repeat=2)]+[op(i) for op in (P,K) for i in range(5)]
    letters=[(p,q,f) for p,q,f in product(range(4),range(4),range(8)) if 1<=p+q+f.bit_count()<=4]
    for a,b in product(operators,repeat=2):
        for l in letters:
            r=Counter()
            for first,second,fac in ((b,a,1),(a,b,-sign(a,b))):
                n,x=action(l,first)
                if n:
                    m,y=action(x,second)
                    r[y]+=n*m*fac
            for c,op in bracket(a,b):
                n,x=action(l,op)
                r[x]-=c*n
            assert not any(r.values()),(a,b,l,r)
    return len(operators)**2*len(letters)


