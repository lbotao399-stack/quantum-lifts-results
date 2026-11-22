"""Charge and exact-mode selection for the separated F4 fifth-glueball Ward integrand."""
from itertools import product
import json

from fifth_instanton import calculate as endpoint


def twice_charge(occupations):
    matter = sum(occupations[:3])
    return tuple(matter-occupations[3]-2*n for n in occupations[:3])


def source_profiles():
    profiles = []
    for ghosts in range(11):
        # Twice the incoming left spin is one. A field strength costs two,
        # a chiral matter fermion or a covariant derivative costs one.
        for matter in (None,0,1,2):
            derivative = int(matter is None)
            fermion = int(matter is not None)
            twice_scalars = tuple(9-ghosts-fermion+2*int(i==matter) for i in range(3))
            if min(twice_scalars) < 0 or any(n%2 for n in twice_scalars):
                continue
            scalars = tuple(n//2 for n in twice_scalars)
            assert 3*(ghosts+fermion)+2*(sum(scalars)+derivative) == 29
            assert (ghosts+derivative)%2 == 0
            counts = tuple(ghosts+sum(scalars)-scalars[i]+int(matter==i) for i in range(3))+(0,)
            assert counts == (9,9,9,0) and twice_charge(counts) == (9,)*3
            profiles.append(dict(gauginos=ghosts,covariant_derivatives=derivative,
                                 matter_flavour=1+matter if matter is not None else None,scalars=scalars))
    assert len(profiles) == 20
    return profiles


def calculate():
    lifted = len(endpoint()["doublets"])
    exact = 4
    capacity = lifted+exact
    assert capacity == 18
    labels = {9:"source",1:"current",-10:"radial"}
    allowed = {label:[] for label in labels.values()}
    for n in product(range(capacity+1),repeat=4):
        charge = twice_charge(n)
        if charge[0] in labels and charge == (charge[0],)*3:
            allowed[labels[charge[0]]].append(n)
    minima = {k:vals[0] for k,vals in allowed.items()}
    assert minima == dict(source=(9,9,9,0),current=(1,1,1,0),radial=(0,0,0,10))
    for key,n in minima.items():
        assert allowed[key] == [tuple(v+k for v in n) for k in range(capacity-max(n)+1)]
    tests = []
    for radial_exact in (0,2):
        possibilities = []
        for source in allowed["source"]:
            for current in allowed["current"]:
                total = tuple(x+y+z for x,y,z in zip(source,current,minima["radial"]))
                assert len(set(total)) == 1
                vertices = capacity-total[0]
                if vertices < 0 or source[3]+current[3]+radial_exact < exact:
                    continue
                # A vertex uses one lifted mode of each flavour and zero exact modes.
                candidates = [s for s in range(exact+1)
                              if s <= source[3] and 0 <= exact-radial_exact-s <= current[3]]
                if candidates:
                    possibilities.append((source[3]+current[3],vertices))
        minimum = min(n for n,_ in possibilities)
        powers = sorted({v for _,v in possibilities})
        assert minimum == exact-radial_exact and powers == list(range(5+radial_exact))
        tests.append(dict(radial_fourth_flavour_exact=radial_exact,radial_fourth_flavour_lifted=10-radial_exact,
                          minimum_source_plus_current_common_shifts=minimum,
                          fourth_flavour_allowed_action_powers=powers))
    return dict(exact_per_flavour=exact,lifted_per_flavour=lifted,orientation="linear radial gauginos",
                twice_R=[9,1,-10],minimum_occupations=minima,source_letter_profiles=source_profiles(),
                linear_radial_tests=tests,physical_nonperturbative_lift="open")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("enable assertions")
    print(json.dumps(calculate(),indent=2))
