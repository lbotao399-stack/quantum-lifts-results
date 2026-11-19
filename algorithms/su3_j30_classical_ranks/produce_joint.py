"""Disjoint intervals, load balanced across processes and Slurm tasks by advisory locks."""
import argparse
import fcntl
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import joint as j


def intervals(plan):
    tasks=[]
    for s in plan['sources']:
        i=s['index'];width=sum(plan['targets'][t]['dimension'] for t in s['targets']); n=s['dimension']
        for lo in range(0,n,64):
            hi=min(n,lo+64); tasks.append(((hi-lo)*width,i,lo,hi))
    return [(i,lo,hi) for _,i,lo,hi in sorted(tasks,reverse=True)]

def work(a):
    plan=j.load_plan(); maps=j.runtime(a.root);cache={}; tasks=intervals(plan);acquired=0;previous=None
    # Start with a different part of the same complete list in each lane. Cheap
    # intervals also return early, with a global scan at the end to pick up work.
    mid=len(tasks)*a.lane//a.total_lanes; order=tasks[mid:]+tasks[:mid]
    for i,lo,hi in order:
        f=a.out/f's{i:02d}'/f'{lo:05d}_{hi:05d}'; receipt=f/'result.json'
        if receipt.exists():continue
        f.mkdir(exist_ok=True,parents=True)
        with (f/'interval.lock').open('a') as lock:
            try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:continue
            if receipt.exists():continue
            j.interval(a.root,a.out,plan,maps,i,lo,hi,a.workers,a.batch,cache);acquired+=1
            if i!=previous or acquired%10==0:print(json.dumps(dict(lane=a.lane,source=i,stop=hi,intervals_done=acquired)),flush=True)
            previous=i
    print(json.dumps(dict(finished_lane=a.lane,intervals_done=acquired)),flush=True)

def launch(a):
    j.inputs(a.root);plan=j.load_plan();tasks=intervals(plan);assert sum(hi-lo for _,lo,hi in tasks)==plan['source_dimension']
    assert 1<=a.procs<=16 and a.procs*a.workers==int(os.environ.get('SLURM_CPUS_PER_TASK',a.procs*a.workers))
    shard=int(os.environ.get('SLURM_ARRAY_TASK_ID',0));nshards=int(os.environ.get('SLURM_ARRAY_TASK_COUNT',1));assert 0<=shard<nshards
    def one(k):
        subprocess.run([sys.executable,'-u',__file__,str(a.root),str(a.out),'--lane',str(shard*a.procs+k),'--total-lanes',str(nshards*a.procs),'--workers',str(a.workers),'--batch',str(a.batch)],check=True)
    with ThreadPoolExecutor(a.procs) as pool:list(pool.map(one,range(a.procs)))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('out',type=Path);p.add_argument('--procs',type=int,default=8);p.add_argument('--workers',type=int,default=8);p.add_argument('--batch',type=int,default=8);p.add_argument('--lane',type=int,default=-1);p.add_argument('--total-lanes',type=int);a=p.parse_args();(launch if a.lane==-1 else work)(a)
