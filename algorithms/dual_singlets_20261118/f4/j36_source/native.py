"""Independent component evaluation in the archive's 52-colour basis."""
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations, combinations_with_replacement, product


def matrix(rows):
    return {(i,j): Fraction(c) for i,j,c in rows}


def multiply(a, b):
    right = defaultdict(list)
    for (i,j),c in b.items():
        right[i].append((j,c))
    row = Counter()
    for (i,j),c in a.items():
        for k,d in right[j]:
            row[i,k] += c*d
    return {k:v for k,v in row.items() if v}


def linear(row, table):
    answer = Counter()
    for i,c in row.items():
        for key,v in table[i].items():
            answer[key] += c*v
    return {k:v for k,v in answer.items() if v}


def transpose(a):
    return {(j,i):v for (i,j),v in a.items()}


def check_colour(tables):
    rep = list(map(matrix,tables['representation26']))
    dual = list(map(matrix,tables['nativeDualTranspose26']))
    metric = matrix(tables['metric'])
    inverse = matrix(tables['inverseMetric'])
    pairing = matrix(tables['representationInvariantForm'])
    identity = {(i,i):1 for i in range(52)}
    assert multiply(metric,inverse) == identity
    dual_rows = [{j:v for (j,k),v in inverse.items() if k == i} for i in range(52)]
    assert dual == [transpose(linear(row,rep)) for row in dual_rows]
    assert transpose(pairing) == pairing
    form = [[pairing.get((i,j),Fraction()) for j in range(26)] for i in range(26)]
    for j in range(26):
        pivot = next(i for i in range(j,26) if form[i][j])
        form[j],form[pivot] = form[pivot],form[j]
        for i in range(j+1,26):
            if form[i][j]:
                scale = form[i][j]/form[j][j]
                form[i] = [u-scale*v for u,v in zip(form[i],form[j])]
    for m in rep:
        left = multiply(transpose(m),pairing)
        assert left == {k:-v for k,v in multiply(pairing,m).items()}
    for i,j in combinations_with_replacement(range(52),2):
        gram = sum((c*rep[j].get((y,x),0) for (x,y),c in rep[i].items()),Fraction())/6
        assert gram == metric.get((i,j),0)
    brackets = {(i,j):dict(row) for i,j,row in tables['bracketPairs']}
    brackets.update({(j,i):{a:-c for a,c in row.items()} for (i,j),row in list(brackets.items())})
    lower = [{a:c for (k,a),c in metric.items() if k == i} for i in range(52)]
    covariant = {}
    for b,c in combinations_with_replacement(range(52),2):
        primal = Counter()
        for i,x in dual_rows[b].items():
            for j,y in dual_rows[c].items():
                for k,z in brackets.get((i,j),{}).items():
                    primal[k] += x*y*z
        native = linear(primal, lower)
        comm = Counter(multiply(dual[b],dual[c]))
        comm.subtract(multiply(dual[c],dual[b]))
        assert {k:-v for k,v in comm.items() if v} == linear(native,dual)
        covariant[b,c] = native
        covariant[c,b] = {k:-v for k,v in native.items()}

    # Check the generating quartic identity on the entire Cartan via polarization.
    assert all(i == j for m in rep[:4] for i,j in m)
    weights = [tuple(m.get((i,i),0) for i in range(26)) for m in rep[:4]]
    def t2(i,j):
        return sum(x*y for x,y in zip(weights[i],weights[j]))
    checks = 0
    for i,j,k,l in combinations_with_replacement(range(4),4):
        t4 = sum(a*b*c*d for a,b,c,d in zip(weights[i],weights[j],weights[k],weights[l]))
        assert 36*t4 == t2(i,j)*t2(k,l)+t2(i,k)*t2(j,l)+t2(i,l)*t2(j,k)
        checks += 1
    return dict(native_bracket_pairs=1378, signed_matrix_commutator=-1,
                cartan_quartic_coefficients=checks), dual, covariant


def odd(jet):
    return 1 - jet[2].bit_count()%2


def emission_sign(u, v):
    left, right = u[2], v[2]
    crossings = sum(i>j for i in range(3) if left & (1<<i)
                    for j in range(3) if right & (1<<j))
    return (-1)**((left|right).bit_count() + (1+left.bit_count())*right.bit_count() + crossings)


def components(primitive, matrices, bracket, pins):
    @cache
    def matrix_word(colours):
        if len(colours) == 1:
            return matrices[colours[0]]
        return multiply(matrix_word(colours[:-1]),matrices[colours[-1]])

    @cache
    def trace(colours):
        return sum(v for (i,j),v in matrix_word(colours).items() if i == j)

    def evaluate(pin):
        wanted = sorted(((p,q,s),c) for p,q,s,colours in pin for c in colours)
        assert len(set(wanted)) == len(wanted)
        full = (1<<len(wanted))-1
        exterior = sum(1<<i for i,(jet,_) in enumerate(wanted) if odd(jet))
        @cache
        def order_sign(a,b):
            exponent = 0
            high = a & exterior
            while high:
                bit = high & -high
                exponent += (b & exterior & (bit-1)).bit_count()
                high ^= bit
            return (-1)**exponent
        singles, pairs = defaultdict(list), defaultdict(list)
        for i,(jet,colour) in enumerate(wanted):
            singles[jet].append((1<<i,colour,1))
        for i,j in combinations(range(len(wanted)),2):
            (u,b),(v,c) = wanted[i],wanted[j]
            if u[2] & v[2]:
                continue
            jet = (u[0]+v[0],u[1]+v[1],u[2]|v[2])
            forwards = emission_sign(u,v)
            backwards = emission_sign(v,u)*(-1)**(odd(u)*odd(v))
            for k in bracket[b,c]:
                value = (forwards*bracket[b,c][k] + backwards*bracket[c,b][k])/2
                if value:
                    pairs[jet].append(((1<<i)|(1<<j),k,value))

        @cache
        def factor(word, differentiated):
            result = Counter()
            for site in (range(len(word)) if differentiated else (-1,)):
                options = [(pairs if i == site else singles)[jet] for i,jet in enumerate(word)]
                before = (-1)**sum(map(odd,word[:site])) if differentiated else 1
                for choices in product(*options):
                    used,sign = 0,before
                    for bits,colour,value in choices:
                        if used & bits:
                            break
                        sign *= order_sign(used,bits)*value
                        used |= bits
                    else:
                        result[used] += sign*trace(tuple(t[1] for t in choices))
            return {k:v for k,v in result.items() if v}

        def combine(old, factor):
            row = Counter()
            for a,x in old.items():
                for b,y in factor.items():
                    if not a & b:
                        row[a|b] += order_sign(a,b)*x*y
            return {k:v for k,v in row.items() if v}

        total = Fraction()
        for c,word in primitive:
            value, derivative = {0:1}, {}
            for t in word:
                left = combine(derivative, factor(t,False))
                right = combine(value, factor(t,True))
                # Every accumulated mask has the parity of the preceding traces.
                before = 1 if not value else (-1)**((next(iter(value)) & exterior).bit_count())
                derivative = Counter(left)
                for k,v in right.items():
                    derivative[k] += before*v
                value = combine(value, factor(t,False))
                if not value and not derivative:
                    break
            total += c*derivative.get(full,0)
        return str(total)
    return {name:evaluate(pin) for name,pin in pins.items()}
