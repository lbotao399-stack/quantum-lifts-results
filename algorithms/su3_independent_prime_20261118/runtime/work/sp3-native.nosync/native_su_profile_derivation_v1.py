"""Exact Q0 profile blocks by selecting the pinned letter replacements."""
import sp_environment
import native_su_derivation_v1 as native
from sp_full_profiles import profile
from collections import Counter
import copy
import numpy as np


def evaluate(expressions,frame,source_profile,target_profile):
    source_profile=tuple((tuple(k),int(n)) for k,n in source_profile)
    target_profile=tuple((tuple(k),int(n)) for k,n in target_profile)
    words=[w for e in expressions for w in e]
    if not words:return np.zeros((len(expressions),frame.ncols),dtype=np.int64)
    assert all(profile(tuple(L for t in w for L in t))==source_profile for w in words)
    if getattr(frame,'_native_su_derivation_configuration',None) is None:
        native.word_images([words[0]],frame)
    original=frame._native_su_derivation_configuration
    caches=getattr(frame,'_native_su_profile_configurations',{})
    key=(source_profile,target_profile)
    if key not in caches:
        ids,parities,offsets,matrices,rules=original
        letters={value:letter for letter,value in ids.items()}
        kind=lambda index:(len(letters[index][2]),sum(letters[index][:2]))
        selected=[]
        for source,coefficient,left,right in rules:
            counts=Counter(dict(source_profile))
            if not counts[kind(source)]:continue
            counts[kind(source)]-=1;counts[kind(left)]+=1;counts[kind(right)]+=1
            changed=tuple(sorted((k,n) for k,n in counts.items() if n))
            if changed==target_profile:selected.append((source,coefficient,left,right))
        caches[key]=ids,parities,offsets,matrices,np.asarray(selected,dtype=np.int32).reshape(-1,4)
        frame._native_su_profile_configurations=caches
    isolated=copy.copy(frame)
    isolated._native_su_derivation_configuration=caches[key]
    return native.evaluate(expressions,isolated)
