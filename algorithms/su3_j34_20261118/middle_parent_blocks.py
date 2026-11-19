"""Allow the two middle predecessor contents singled out by the integer witness.

Cancel their side contents with all sampled side-only parents. Optionally replace
both distinguished predecessor lists with complete canonical multitrace lists.
Always start from the fixed 978-word top vector; all other equations in the
70-content projection are unchanged by construction, in the formal trace algebra.
"""
import argparse
from collections import Counter
from pathlib import Path
from fractions import Fraction
import scalar_filtration as l
import content_parents
import content_dimensions as dims
s,f,p,t=l.s,l.f,l.p,l.tri
GROUPS=tuple(bytes.fromhex(x) for x in ('00000000000105050712','00000000000105050d10'))

def compile_recipe(root,rhsfile,vector,groups=GROUPS,all_parents=False,frames=()):
    keep={c for c in l.contents(11,0,l.CHARGES) if c.count(0)>=6}|set(t.CONTENTS);rhs=p.rhs_rows(rhsfile,keep);assert len(keep)==70
    b={w:c for w,c in rhs.items() if p.content(w) in t.CONTENTS}
    for w,c in l.load_vector(vector).items():p.axpy(b,{u:k for u,k in f.q(w).items() if p.content(u) in t.CONTENTS},-c)
    assert {p.content(u) for u in b}=={t.MID};par={w for u in set(b)|set(rhs) for w in f.parents(u)};extra={}
    if all_parents:
        for g in groups:
            words=content_parents.all_words(g);extra[g.hex()]=len(words);par.update(words)
    frames=tuple(frames);fr={}
    for path in frames:
        frame=s.read(path);g=bytes.fromhex(frame['content']);assert g in groups and g not in fr and len(frame['words'])==dims.dimension(g)==frame['dimension'] and frame['determinant_mod_prime'];fr[g]=len(frame['words']);par={w for w in par if p.content(w)!=g};par.update(s.unpack_hex(w) for w in frame['words'])
    assert not(fr and all_parents)
    ims={};cs={};side=set()
    for w in sorted(par):
        im={u:k for u,k in f.q(w).items() if p.content(u) in keep};cset={p.content(u) for u in im}
        if not cset:continue
        ims[w]=im;cs[w]=cset
        if p.content(w) in groups:side.update(cset-{t.MID})
    chosen={w:im for w,im in ims.items() if p.content(w) in groups or cs[w]<={t.MID} or cs[w]<=side};assert not (set.union(*(cs[w] for w in chosen)) - side - {t.MID});z=l.recipe(root,chosen,b,side|{t.MID},'Fix the 978-word top vector; both distinguished contents and selected middle-only or side-only words can vary.');hist=Counter(p.content(w).hex() for w in chosen);s.write(root/'selection.json',dict(groups=[k.hex() for k in groups],fully_enumerated_formal_words=extra,full_three_colour_parent_dimensions={g.hex():n for g,n in fr.items()},selected_by_content=hist,selected=len(chosen),trace_words=len(z['word']),contents=z['contents']));return z
if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--work',type=Path,required=True);a.add_argument('--rhs',type=Path,default=s.HERE/'negative_Q1.jsonl.gz');a.add_argument('--top',type=Path,default=s.HERE/'scalar_layers/top_vector.json.gz');a.add_argument('--groups',nargs='*',default=[x.hex() for x in GROUPS]);a.add_argument('--all-parents',action='store_true');a.add_argument('--frames',type=Path,nargs='*',default=[]);b=a.parse_args();r=compile_recipe(b.work,b.rhs,b.top,tuple(map(bytes.fromhex,b.groups)),b.all_parents,b.frames);print(len(r['source']),len(r['word']),r['contents'],flush=True)
