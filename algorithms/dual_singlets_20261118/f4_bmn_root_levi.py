"""Check Cartan restriction for simple Levi factors in two exceptional BMN source sectors."""
import argparse
import gzip
import json
from pathlib import Path

from su3_levi_kernel import check, compute, gzip_json, key


PRIME = 2147483647
CHARACTERS = {"SU2": ((2,"SU+"),), "SU3": ((3,"SU+"),(3,"SU-")),
              "Sp2": ((4,"Sp"),), "SO7": ((7,"SO"),), "Sp3": ((6,"Sp"),)}


def sectors(width):
    # A partition of an active factor sits inside the outer rectangle (width,width,width).
    for first in range(width + 1):
        for second in range(first + 1):
            for third in range(second + 1):
                weight = (first, second, third, 0, 0)
                for excess in (0, 1):
                    length = first + second + third - excess
                    if 0 < length <= 3*width - 1:
                        yield weight, length


def verify(records, group, widths=(4, 5)):
    from trace_descent import dimensions
    output = {}
    for width in widths:
        total = classes = kernel = odd_classes = 0
        blocks = []
        for weight, length in sectors(width):
            label = key(weight, length)
            dim, upper, invisible = check(records[label], weight, length, PRIME,
                                          None if group == "SU2" else group)
            independent = sum(sum(dimensions((*weight[3:], *weight[:3]), length, size, colour).values())
                              for size, colour in CHARACTERS[group])
            assert independent == dim
            total += dim
            classes += upper
            kernel += invisible
            if length % 2:
                odd_classes += upper
            if invisible:
                blocks.append(dict(weight=weight, length=length, cartan_kernel_upper=invisible))
        output[str(width)] = dict(level=6*width, total_length=3*width-1,
                                  active_blocks=len(list(sectors(width))), sum_of_cochain_dimensions=total,
                                  sum_of_classical_upper_bounds=classes, sum_of_kernel_upper_bounds=kernel,
                                  sum_of_odd_length_classical_upper_bounds=odd_classes, remaining=blocks)
    return dict(prime=PRIME, group="F4", active_group=group, sectors=output)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--records",type=Path,action="append",default=[])
    p.add_argument("--group",choices=tuple(CHARACTERS),default="SU2")
    p.add_argument("--width",type=int,action="append",choices=(4,5))
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--engine",type=Path)
    p.add_argument("--raw",type=Path)
    p.add_argument("--workers",type=int,default=2)
    args=p.parse_args()
    if not __debug__ or bool(args.engine) != bool(args.raw) or args.workers<1 or not (args.engine or args.records):
        p.error("supply --raw with --engine or provide --records for replay, with assertions enabled")
    from su3_third_levi_weights import read_records
    previous=read_records(args.records)
    widths=sorted(set(args.width or (4,5)))
    needed={key(w,y):(w,y) for w,y in sectors(max(widths))}
    if args.engine:
        todo=[dict(weight=w,length=y) for label,(w,y) in needed.items() if label not in previous]
        previous.update(compute({"new":dict(active=todo)},args.engine,args.raw,args.workers,PRIME,
                                None if args.group=="SU2" else args.group))
        gzip_json(args.out,{label:previous[label] for label in needed})
    else:
        args.out.write_text(json.dumps(verify(previous,args.group,widths),indent=2)+"\n")
