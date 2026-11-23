"""Stable native Q0 complex at weight (3,3,9,0,0), with bounded checkpoints.

All polynomial keys and coefficients are the original mixed_incoming/trace_descent
conventions.  Rank computations are modular lower bounds until the exact Q0^2
identity and the character dimension force saturation.
"""
from pathlib import Path
import os, sys, json, time, signal, resource, pickle, hashlib, argparse
from collections import defaultdict

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT/'upstream/algorithms/dual_singlets_20261118'
RUNTIME = ROOT/'upstream/algorithms/su3_independent_prime_20261118/runtime/work/sp3-native.nosync'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['N4_TRACE_RUNTIME'] = str(RUNTIME)
sys.path.insert(0,str(SOURCE))
import mixed_incoming as m
from flint import nmod_mat

WEIGHT = (3,3,9,0,0)
EXPECTED = {12:38160,13:37085,14:17745}

def dump(path, obj):
    tmp = path.with_suffix(path.suffix+'.part')
    tmp.write_bytes(pickle.dumps(obj,protocol=4))
    tmp.replace(path)

def jdump(path,obj):
    tmp = path.with_suffix(path.suffix+'.part')
    tmp.write_text(json.dumps(obj,indent=2)+'\n')
    tmp.replace(path)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source_manifest():
    return {'weight':WEIGHT,'source_files':{
        str(path.relative_to(ROOT)):sha(path) for path in
        [SOURCE/'mixed_incoming.py', SOURCE/'trace_descent.py',Path(m.fw.__file__)]
    },'native_weight_order':'first two Taylor indices, final three theta occupation indices',
       'basis_convention':'trace_descent.highest_weights flattened values in insertion order',
       'differential_convention':'mixed_incoming.q without scaling or normalization'}

def prepare_basis(length):
    final = HERE/f'formal_339_L{length}.pkl'
    if final.exists():
        rows = pickle.loads(final.read_bytes())
        assert len(rows)==EXPECTED[length]
        print('basis cached',length,len(rows),flush=True)
        return rows
    checkpoint = HERE/f'basis_339_L{length}_checkpoint.pkl'
    start = time.monotonic()
    blocks = m.words_by_shape(length,WEIGHT)
    print('enumerated',length,len(blocks),sum(map(len,blocks.values())),flush=True)
    if checkpoint.exists():
        next_block,by_degree = pickle.loads(checkpoint.read_bytes())
        assert isinstance(by_degree,dict)
    else:
        next_block,by_degree=0,{}
    for block_id,block in enumerate(blocks.values()):
        if block_id<next_block:
            continue
        for row in m.singlets(block):
            degrees={m.nferm(w) for w in row}
            assert len(degrees)==1
            by_degree.setdefault(degrees.pop(),[]).append(row)
        if (block_id+1)%100==0 or block_id+1==len(blocks):
            dump(checkpoint,(block_id+1,by_degree))
            print('basis progress',length,block_id+1,len(blocks),sum(map(len,by_degree.values())),
                  'seconds',time.monotonic()-start,flush=True)
    rows=[row for degree_rows in by_degree.values() for row in degree_rows]
    assert len(rows)==EXPECTED[length],(length,len(rows),EXPECTED[length])
    dump(final,rows)
    jdump(HERE/f'formal_339_L{length}.json',dict(source_manifest(),length=length,
          dimension=len(rows),degree_dimensions={str(k):len(v) for k,v in by_degree.items()},
          source_shapes=len(blocks),source_words=sum(map(len,blocks.values())),
          basis_sha256=sha(final),seconds=time.monotonic()-start))
    print('basis completed',length,len(rows),'seconds',time.monotonic()-start,flush=True)
    return rows

def prepare_map(length):
    rows=prepare_basis(length)
    final=HERE/f'q_339_L{length}.pkl'
    if final.exists():
        images=pickle.loads(final.read_bytes())
        assert len(images)==len(rows)
        print('Q cached',length,len(images),flush=True)
        return images
    checkpoint=HERE/f'q_339_L{length}_checkpoint.pkl'
    images=pickle.loads(checkpoint.read_bytes()) if checkpoint.exists() else []
    start=time.monotonic()
    for row_id in range(len(images),len(rows)):
        images.append(m.q(rows[row_id]))
        if (row_id+1)%1000==0 or row_id+1==len(rows):
            dump(checkpoint,images)
            print('Q progress',length,row_id+1,len(rows),'seconds',time.monotonic()-start,flush=True)
    dump(final,images)
    jdump(HERE/f'q_339_L{length}.json',dict(source_manifest(),length=length,
        source_dimension=len(rows),nonzero_images=sum(bool(r) for r in images),
        image_terms=sum(map(len,images)),q_sha256=sha(final),seconds=time.monotonic()-start))
    print('Q completed',length,len(images),'seconds',time.monotonic()-start,flush=True)
    return images

def check_q2(length):
    images=prepare_map(length)
    reportpath=HERE/f'q2_339_L{length}.json'
    report=json.loads(reportpath.read_text()) if reportpath.exists() else {'checked':0}
    start=time.monotonic()
    for row_id in range(report['checked'],len(images)):
        assert not m.q(images[row_id]),('nonzero exact Q0 squared',length,row_id)
        if (row_id+1)%1000==0 or row_id+1==len(images):
            report.update(length=length,checked=row_id+1,source_dimension=len(images),
                exact_integer=True,complete=row_id+1==len(images),q_sha256=sha(HERE/f'q_339_L{length}.pkl'),
                seconds=time.monotonic()-start)
            jdump(reportpath,report)
            print('Q2 progress',length,row_id+1,len(images),'seconds',time.monotonic()-start,flush=True)
    return report

def components(images):
    parent=list(range(len(images)))
    owners={}
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i
    for i,row in enumerate(images):
        for word in row:
            if word in owners:
                parent[find(i)]=find(owners[word])
            else:
                owners[word]=i
    groups={}
    for i,row in enumerate(images):
        if row:
            groups.setdefault(find(i),[]).append(i)
    return sorted(groups.values(),key=lambda g:(len(g),g[0]))

def pivot_columns(rref,rank):
    pivots=[]
    column=0
    for row in range(rank):
        while not rref[row,column]:
            column+=1
        pivots.append(column)
        column+=1
    return pivots

def basis_pivot_words(length,prime):
    """Choose independent ambient-word coordinates on the native HW basis."""
    path=HERE/f'pivot_words_formal_339_L{length}_p{prime}.pkl'
    basis_path=HERE/f'formal_339_L{length}.pkl'
    if path.exists():
        data=pickle.loads(path.read_bytes())
        assert data['basis_sha256']==sha(basis_path)
        return data['pivot_words']
    rows=prepare_basis(length)
    groups=components(rows)
    result=[]
    records=[]
    determinant_product=1
    for block_id,group in enumerate(groups):
        words=sorted({word for i in group for word in rows[i]})
        ids={w:j for j,w in enumerate(words)}
        A=nmod_mat(len(group),len(words),prime)
        for i,source_id in enumerate(group):
            for w,c in rows[source_id].items():
                A[i,ids[w]]=c%prime
        rr,rank=A.rref(inplace=True)
        assert rank==len(group),(length,block_id,rank,len(group))
        pivots=pivot_columns(rr,rank)
        selected_words=[words[j] for j in pivots]
        minor=nmod_mat([[rows[i].get(w,0)%prime for w in selected_words]
                        for i in group],prime)
        determinant=int(minor.det())
        assert determinant
        determinant_product=determinant_product*determinant%prime
        result.extend(selected_words)
        records.append({'source_rows':group,'pivot_words':selected_words,
                        'minor_determinant':determinant})
    assert len(result)==len(set(result))==EXPECTED[length]
    data={'length':length,'prime':prime,'basis_sha256':sha(basis_path),
          'rank':len(result),'pivot_words':sorted(result),'blocks':records,
          'block_determinant_product_mod_prime':determinant_product,
          'determinant_convention':'rows and selected columns concatenated in listed block order'}
    dump(path,data)
    print('basis pivot coordinates',length,len(result),'blocks',len(groups),flush=True)
    return data['pivot_words']

def rank_map(length,prime,save_reduction=False,restrict_target=False):
    images=prepare_map(length)
    qhash=sha(HERE/f'q_339_L{length}.pkl')
    groups=components(images)
    allowed=set(basis_pivot_words(length+1,prime)) if restrict_target else None
    restriction='certified target HW pivot words' if restrict_target else 'all output words'
    work=HERE/f'rank_339_L{length}_p{prime}'
    work.mkdir(exist_ok=True)
    dump(work/'components.pkl',groups)
    start=time.monotonic()
    report={'weight':WEIGHT,'length':length,'prime':prime,'q_sha256':qhash,
            'rank_kind':'modular lower bound','blocks':[],'source_dimension':len(images),
            'coordinate_restriction':restriction}
    print('components',length,len(groups),'largest',max(map(len,groups),default=0),flush=True)
    for block_id,group in enumerate(groups):
        path=work/f'block_{block_id}.pkl'
        if path.exists():
            data=pickle.loads(path.read_bytes())
            assert data['record']['q_sha256']==qhash
            assert data['source_rows']==group
            assert data['record'].get('coordinate_restriction','all output words')==restriction
            record=data['record']
        else:
            words=sorted({word for i in group for word in images[i]
                          if allowed is None or word in allowed})
            ids={word:j for j,word in enumerate(words)}
            print('rank starting',length,block_id,len(group),len(words),flush=True)
            A=nmod_mat(len(group),len(words),prime)
            for i,source_id in enumerate(group):
                for word,coefficient in images[source_id].items():
                    if word in ids:
                        A[i,ids[word]]=coefficient%prime
            if save_reduction:
                rr,rank=A.rref(inplace=True)
                pivots=pivot_columns(rr,rank)
                del rr,A
                columns=nmod_mat(rank,len(group),prime)
                for i,column in enumerate(pivots):
                    for j in range(len(group)):
                        columns[i,j]=images[group[j]].get(words[column],0)%prime
                cr,rank2=columns.rref(inplace=True)
                assert rank2==rank
                selected=pivot_columns(cr,rank)
                del cr,columns
                selected_ids=[group[i] for i in selected]
                minor=nmod_mat(rank,rank,prime)
                for i,source_id in enumerate(selected_ids):
                    row=images[source_id]
                    for j,column in enumerate(pivots):
                        minor[i,j]=row.get(words[column],0)%prime
                determinant=int(minor.det())
                assert determinant
                del minor
                data={'source_rows':group,'words':words,'pivot_columns':pivots,
                      'selected_source_rows':selected_ids,'minor_determinant':determinant}
            else:
                rr,rank=A.rref(inplace=True)
                pivots=pivot_columns(rr,rank)
                del rr,A
                data={'source_rows':group,'words':words,'pivot_columns':pivots}
            record={'block':block_id,'source_rows':len(group),'words':len(words),
                    'rank':rank,'q_sha256':qhash,'prime':prime,
                    'coordinate_restriction':restriction,
                    'seconds':time.monotonic()-start,
                    'maxrss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
            data['record']=record
            dump(path,data)
        report['blocks'].append(record)
        report.update(rank=sum(r['rank'] for r in report['blocks']),
                      completed_blocks=len(report['blocks']),total_blocks=len(groups),
                      complete=len(report['blocks'])==len(groups),
                      seconds=time.monotonic()-start)
        jdump(work/'progress.json',report)
        print('rank block',length,block_id,record['rank'],'cumulative',report['rank'],
              'seconds',time.monotonic()-start,flush=True)
    jdump(HERE/f'rank_339_L{length}_p{prime}.json',report)
    return report

def summarize(prime):
    summary=source_manifest()
    verified_dimensions={}
    for length in (12,13,14):
        path=HERE/f'formal_339_L{length}.json'
        if path.exists():
            verified_dimensions[length]=json.loads(path.read_text())['dimension']
    summary.update(stable_dimensions=verified_dimensions,finite_Sp5_dimensions_given={12:38093,13:36897,14:17562},
                   prime=prime)
    ranks={}
    for length in (12,13):
        path=HERE/f'rank_339_L{length}_p{prime}.json'
        if path.exists():
            report=json.loads(path.read_text())
            if report['complete']:
                ranks[length]=report['rank']
    checks={}
    hw_checks={}
    for length in (12,13):
        path=HERE/f'q2_339_L{length}.json'
        if path.exists():
            checks[length]=json.loads(path.read_text())['complete']
        hw_path=HERE/f'q_hw_339_L{length}.json'
        if hw_path.exists():
            hw_checks[length]=json.loads(hw_path.read_text())['complete']
    summary.update(modular_rank_lower_bounds=ranks,exact_Q0_squared_checks=checks,
                   exact_target_highest_weight_checks=hw_checks)
    if (set(ranks)=={12,13} and checks.get(12) and checks.get(13)
        and hw_checks.get(12) and hw_checks.get(13) and verified_dimensions==EXPECTED):
        upper=EXPECTED[13]-ranks[12]-ranks[13]
        assert upper>=0
        summary['stable_H13_dimension_upper_bound']=upper
        if upper==0:
            summary.update(stable_H13_dimension_exact=0,rational_ranks_exact=ranks,
                proof='Exact Q0 squared zero gives rank_Q12+rank_Q13 <= dim C13; modular minor lower bounds attain this dimension.')
    jdump(HERE/'stable_complex_result.json',summary)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('phase',choices=['basis','map','q2','rank','summary'])
    parser.add_argument('--length',type=int,choices=[12,13,14],default=13)
    parser.add_argument('--prime',type=int,default=2147483647)
    parser.add_argument('--seconds',type=int,default=300)
    parser.add_argument('--memory-gib',type=float,default=5.5)
    parser.add_argument('--save-reduction',action='store_true')
    parser.add_argument('--restrict-target',action='store_true')
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(int(args.memory_gib*1024**3),)*2)
    signal.signal(signal.SIGALRM,lambda *_:(_ for _ in ()).throw(TimeoutError('bounded stable complex job')))
    signal.alarm(args.seconds)
    if args.phase=='basis': prepare_basis(args.length)
    elif args.phase=='map': prepare_map(args.length)
    elif args.phase=='q2': check_q2(args.length)
    elif args.phase=='rank': rank_map(args.length,args.prime,args.save_reduction,args.restrict_target)
    elif args.phase=='summary': summarize(args.prime)
