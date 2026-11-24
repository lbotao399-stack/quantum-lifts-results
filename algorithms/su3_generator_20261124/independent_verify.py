"""Independent integer-polynomial and exterior-matrix verification of B05.

This checker imports none of the construction algebra, map, or frame modules.
It reconstructs the entire derivative from the printed 12224-term operator,
checks both dotted raising and lowering, and evaluates its full 39829-term
stable residual at every coordinate of the complete SU(3) target frame.
Polynomial combinations use unbounded Python integers.
"""
import argparse,gzip,hashlib,json,math,resource,time
from collections import Counter,defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations
from pathlib import Path
import numpy as np


def load(path):return json.loads(gzip.decompress(path.read_bytes()) if path.suffix==".gz" else path.read_text())
def sparse(entries):return {j:Fraction(n,d) for j,n,d in entries}


class Algebra:
    def __init__(self,letters):
        self.letters=[tuple(letter) for letter in letters];self.ids={letter:i for i,letter in enumerate(self.letters)}
        assert all(mask in (0,1) for a,b,mask in self.letters)
    def odd(self,j):return self.letters[j][2]==0
    @cache
    def trace(self,word):
        if len(word)<2:return 0,()
        choices={};phase=1;current=word
        for _ in word:
            if current in choices and choices[current]!=phase:return 0,()
            choices[current]=phase
            if self.odd(current[0]) and sum(map(self.odd,current[1:]))%2:phase=-phase
            current=current[1:]+current[:1]
        canonical=min(choices);return choices[canonical],canonical
    @cache
    def multi(self,word):
        phase=1;traces=[]
        for trace in word:
            sign,value=self.trace(trace)
            if not sign:return 0,()
            phase*=sign;traces.append(value)
        odd_traces=[trace for trace in traces if sum(map(self.odd,trace))%2]
        if len(set(odd_traces))!=len(odd_traces):return 0,()
        phase*=(-1)**sum(left>right for i,left in enumerate(odd_traces) for right in odd_traces[i+1:])
        return phase,tuple(sorted(traces))
    @cache
    def letter_derivative(self,j):
        a,b,mask=self.letters[j];result=[]
        for u in range(a+1):
            for v in range(b+1):
                if mask==0:
                    if u+v and a+b-u-v:
                        result.append((1,self.ids[u,v,0],self.ids[a-u,b-v,0]))
                else:
                    if u+v:result.append((1,self.ids[u,v,0],self.ids[a-u,b-v,1]))
                    if a+b-u-v:result.append((-1,self.ids[u,v,1],self.ids[a-u,b-v,0]))
        return result
    @cache
    def derivative_word(self,word):
        result=Counter();phase=1
        for i,trace in enumerate(word):
            for j,letter in enumerate(trace):
                for coefficient,left,right in self.letter_derivative(letter):
                    changed=trace[:j]+(left,right)+trace[j+1:]
                    sign,target=self.multi(word[:i]+(changed,)+word[i+1:])
                    if sign:result[target]+=phase*coefficient*sign
                if self.odd(letter):phase=-phase
        return tuple((w,c) for w,c in result.items() if c)
    def derivative(self,polynomial):
        result=Counter()
        for word,c in polynomial.items():
            for target,value in self.derivative_word(word):result[target]+=c*value
        return {w:c for w,c in result.items() if c}
    def spin(self,polynomial,up):
        result=Counter()
        for word,c in polynomial.items():
            for i,trace in enumerate(word):
                for j,letter in enumerate(trace):
                    a,b,mask=self.letters[letter]
                    if (up and not b) or (not up and not a):continue
                    changed=(a+1,b-1,mask) if up else (a-1,b+1,mask)
                    # A changed letter can lie outside the archived charge box;
                    # retain its actual triple in an independent temporary word.
                    traces=tuple(tuple(self.letters[t] for t in tr) for tr in word)
                    value=traces[:i]+(traces[i][:j]+(changed,)+traces[i][j+1:],)+traces[i+1:]
                    sign,target=canonical_triples(value)
                    if sign:result[target]+=c*sign*((a+1) if up else (b+1))
        return {w:c for w,c in result.items() if c}
    def conjugation(self,polynomial):
        result=Counter()
        for word,c in polynomial.items():
            phase=1
            for trace in word:
                count=sum(map(self.odd,trace));phase*=(-1)**(len(trace)+count*(count-1)//2)
            sign,target=self.multi(tuple(trace[::-1] for trace in word))
            if sign:result[target]+=c*sign*phase
        return {w:c for w,c in result.items() if c}
    def skeleton(self,word):
        return tuple(sorted(Counter((mask,a+b) for trace in word for a,b,mask in (self.letters[j] for j in trace)).items()))


@cache
def canonical_triples(multi):
    phase=1;traces=[]
    for trace in multi:
        if len(trace)<2:return 0,()
        choices={};current=trace;sign=1
        for _ in trace:
            if current in choices and choices[current]!=sign:return 0,()
            choices[current]=sign
            if current[0][2]==0 and sum(t[2]==0 for t in current[1:])%2:sign=-sign
            current=current[1:]+current[:1]
        canonical=min(choices);phase*=choices[canonical];traces.append(canonical)
    odd=[trace for trace in traces if sum(t[2]==0 for t in trace)%2]
    if len(set(odd))!=len(odd):return 0,()
    phase*=(-1)**sum(a>b for i,a in enumerate(odd) for b in odd[i+1:])
    return phase,tuple(sorted(traces))


@cache
def exterior_table(k,left,right):
    subsets={degree:list(combinations(range(k),degree)) for degree in (left,right,left+right)}
    indices={degree:{value:i for i,value in enumerate(values)} for degree,values in subsets.items()}
    first=[];second=[];signs=[];starts=[]
    for output in subsets[left+right]:
        starts.append(len(first))
        for a in combinations(output,left):
            b=tuple(j for j in output if j not in a)
            first.append(indices[left][a]);second.append(indices[right][b]);signs.append((-1)**sum(i>j for i in a for j in b))
    return np.asarray(first),np.asarray(second),np.asarray(signs,dtype=np.int64),np.asarray(starts)


def input_matrices(block,algebra):
    random=np.random.default_rng(block["seed"]);result={};degree=block["degree"]
    chosen_points=sorted({column//(degree+1) for column in block["columns"]})
    for letter in block["alphabet"]:
        odd=int(algebra.odd(letter));number=degree+1 if odd else 1
        values=np.zeros((block["points"],number,3,3),np.int64)
        for matrix in values.reshape(-1,3,3):
            for _ in range(2):
                order=random.permutation(3)
                if random.integers(2):matrix[order,np.roll(order,1)]+=random.integers(-1,2,size=3)
                else:
                    coefficient=int(random.integers(-1,2));matrix[order[0],order[0]]+=coefficient;matrix[order[1],order[1]]-=coefficient
        assert not values.trace(axis1=-2,axis2=-1).any() and abs(values).sum(axis=-1).max(initial=0)<=2
        result[letter]=(odd,values[chosen_points])
    reindex={point:i for i,point in enumerate(chosen_points)}
    columns=[reindex[column//(degree+1)]*(degree+1)+column%(degree+1) for column in block["columns"]]
    return result,columns


def evaluate_words(words,block,algebra):
    letters,columns=input_matrices(block,algebra);k=block["degree"]+1
    necklaces=sorted({trace for word in words for trace in word});values={};stack=[];previous=()
    for trace in necklaces:
        shared=0
        while shared<min(len(previous),len(trace)) and previous[shared]==trace[shared]:shared+=1
        del stack[shared:]
        if not stack:
            odd,array=letters[trace[0]];stack.append((array,odd))
        for letter in trace[len(stack):]:
            left,left_degree=stack[-1];right_degree,right=letters[letter]
            if not right_degree:result=left@right
            else:
                a,b,sign,starts=exterior_table(k,left_degree,right_degree)
                result=np.add.reduceat((left[:,a]@right[:,b])*sign[None,:,None,None],starts,axis=1)
            stack.append((result,left_degree+right_degree))
        array,degree=stack[-1];values[trace]=(array.trace(axis1=-2,axis2=-1),degree);previous=trace
    answer=[]
    for word in words:
        value,degree=values[word[0]]
        for trace in word[1:]:
            other,other_degree=values[trace];a,b,sign,starts=exterior_table(k,degree,other_degree)
            value=np.add.reduceat(value[:,a]*other[:,b]*sign[None,:],starts,axis=1);degree+=other_degree
        assert degree==block["degree"]
        bound=math.factorial(degree)*2**sum(map(len,word))*3**len(word)
        assert abs(value).max(initial=0)<=bound<2**63
        answer.append(value.reshape(-1)[columns])
    return np.asarray(answer,np.int64)


def verify(directory):
    start=time.monotonic();basis=load(directory/"basis.json.gz");generator=load(directory/"generator.json.gz")
    maps=load(directory/"maps.json.gz");frames=load(directory/"frames.json.gz")
    words=[tuple(tuple(trace) for trace in word) for word in basis["words"]];algebra=Algebra(basis["letters"])
    polynomial={words[j]:Fraction(n,d) for j,n,d in generator["full_polynomial"]}
    denominator=math.lcm(*(c.denominator for c in polynomial.values()))
    polynomial={word:int(c*denominator) for word,c in polynomial.items()}
    for word in polynomial:
        letters=[algebra.letters[j] for trace in word for j in trace]
        assert len(letters)==7 and sum(t[0] for t in letters)==4 and sum(t[1] for t in letters)==4 and sum(t[2] for t in letters)==3
    assert not algebra.spin(polynomial,True) and not algebra.spin(polynomial,False)
    assert algebra.conjugation(polynomial)=={word:-c for word,c in polynomial.items()}
    derivative=algebra.derivative(polynomial)
    archived={words[j]:Fraction(n,d)*denominator for j,n,d in generator["full_stable_residual"]}
    assert derivative==archived and derivative
    print(json.dumps(dict(stage="independent_polynomial",operator_terms=len(polynomial),stable_residual_terms=len(derivative),seconds=time.monotonic()-start)),flush=True)
    grouped=defaultdict(dict)
    for word,c in derivative.items():grouped[algebra.skeleton(word)][word]=c
    checked=0;blocks=[]
    for i,block in enumerate(frames["8"]):
        key=tuple((tuple(label),count) for label,count in block["skeleton"])
        row=grouped.pop(key,{})
        if row:
            selected=sorted(row);values=evaluate_words(selected,block,algebra)
            coefficients=np.asarray([row[word] for word in selected],dtype=object)
            result=coefficients@values.astype(object)
            assert all(value==0 for value in result)
        checked+=len(block["columns"]);blocks.append(dict(block=i,polynomial_terms=len(row),coordinates=len(block["columns"]),nonzero_residuals=0))
        print(json.dumps(dict(stage="independent_complete_target",**blocks[-1],seconds=time.monotonic()-start)),flush=True)
    assert not grouped and checked==1588
    dual=sparse(generator["dual"])
    pairing=lambda row:sum(Fraction(n,d)*dual.get(j,0) for j,n,d in row)
    assert pairing(generator["representative"])==1
    assert all(pairing(row)==0 for row in maps["6"]["finite"])
    assert all(pairing(row)==0 for row in generator["finite_proper"]+generator["finite_products"])
    receipt=dict(independent_implementation=True,operator_terms=len(polynomial),stable_residual_terms=len(derivative),
                 integral_operator_denominator=denominator,all_charge_checks=True,dotted_raising_residual_terms=0,dotted_lowering_residual_terms=0,
                 colour_conjugation_residual_terms=0,full_target_coordinates_checked=checked,full_polynomial_residual_nonzero=0,
                 boundary_columns_checked=len(maps["6"]["finite"]),proper_product_columns_checked=len(generator["finite_proper"]),
                 all_stable_graviton_columns_checked=len(generator["finite_products"]),representative_dual_pairing=1,blocks=blocks,
                 coefficient_arithmetic="unbounded Python integers and exact Fraction",seconds=time.monotonic()-start,
                 maximum_resident_set_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (directory/"basis.json.gz",directory/"maps.json.gz",directory/"frames.json.gz",directory/"generator.json.gz")})
    return receipt


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("directory",type=Path,nargs="?",default=Path(__file__).resolve().parent/"certificate")
    parser.add_argument("--receipt",type=Path);args=parser.parse_args();assert __debug__
    result=verify(args.directory)
    if args.receipt:args.receipt.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
