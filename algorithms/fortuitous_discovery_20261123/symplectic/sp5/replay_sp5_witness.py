"""Verify the exact rational Sp(5) closure and trace-coordinate separating dual."""
from pfaffian_relations import *
from fractions import Fraction
import argparse,gzip,hashlib

def main(output):
 start=time.monotonic()
 receipt=json.loads((HERE/'sp5_discovery_certificate.json').read_text())
 exported=json.loads(gzip.decompress((HERE/'sp5_new_fortuitous_operator.json.gz').read_bytes()))
 polynomial={m.decode(w):Fraction(c) for c,w in exported['terms']}
 scale=math.lcm(*(c.denominator for c in polynomial.values()))
 integer={w:int(c*scale) for w,c in polynomial.items()}
 variation=m.q(integer)
 relations,_=pickle.loads((HERE/'target_sp5_snapshot_181.pkl').read_bytes())
 source,_=pickle.loads((HERE/'exact_relations_sp5_L13.pkl').read_bytes())
 assert len(source)==188 and len(relations)==181
 coefficients=[Fraction(c) for c in receipt['defect_identity_coefficients']]
 den=math.lcm(scale,*(c.denominator for c in coefficients));combined=Counter()
 for coefficient,row in zip(coefficients,relations):
  for word,c in row.items():combined[word]+=int(coefficient*den)*c
 assert {w:c for w,c in combined.items() if c}=={w:c*(den//scale) for w,c in variation.items()}
 for raising in (1,2,3):
  raised=Counter()
  for word,c in integer.items():
   for image,n in m.raise_word(word,raising).items():raised[image]+=c*n
  assert not any(raised.values())
 print('All operator coefficients, exact closure, and native highest weights passed.',flush=True)
 dual=[(m.decode(w),Fraction(c)) for c,w in receipt['dual_word_coefficients']]
 pairings=[]
 for i,row in enumerate(source):
  image=m.q(row)
  pairing=sum(c*image.get(w,0) for w,c in dual)
  assert pairing==0,('nonzero boundary pairing',i)
  pairings.append(str(pairing))
  if (i+1)%20==0:print('Exact dual source pairings',i+1,'seconds',time.monotonic()-start,flush=True)
 pairing=sum(c*variation.get(w,0) for w,c in dual)/scale
 assert pairing==1
 result={'passed':True,'group':'Sp(5)','native_weight':[3,3,9,0,0],'length':13,'j':36,
  'operator_terms':len(integer),'operator_denominator':scale,'variation_terms':len(variation),
  'full_exact_closure_verified':True,'native_highest_weights_verified':True,
  'dual_field':'Q','dual_nonzero_coefficients':len(dual),'actual_source_Q0_images_recomputed':188,
  'dual_source_pairings':pairings,'dual_defect_pairing':str(pairing),
  'remaining_proof_input':'The identity recipe replay verifies finite identities and complete source character dimension.',
  'seconds':time.monotonic()-start}
 output.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k!='dual_source_pairings'}),flush=True)

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=Path,default=HERE/'exact_witness_replay.json')
 args=a.parse_args();main(args.out)
