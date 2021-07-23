#!/usr/bin/env python

import glob
import collections
import re

from argparse import ArgumentParser

parser = ArgumentParser()

# data input
parser.add_argument("--unique", action="store_true", default=False, help="Don't use only the CDR3 region as key, but add V and J best hits")
parser.add_argument("--exclude-j", action="store_true", default=False, help="If --unique is set, we may want to only use the V+CDR3 as key")
parser.add_argument("--cdr3nt-only", action="store_true", default=False, help="If --unique is set, we may want to only use the CDR3 (nucleotides) as key")
parser.add_argument("--cdr3aa-only", action="store_true", default=False, help="If --unique is set, we may want to only use the CDR3 (aminoacids) as key")
parser.add_argument("--type", type=str, default='raw', choices=['raw', 'gene', 'family'])
parser.add_argument("indices", type=int, nargs="+", help="the list of MIDs to use")

args = parser.parse_args()

# FIXME it doesn't account for a LIST of hits
def tname(s, t):
    if t == 'gene':
        s = re.sub('\*[0-9]+$', '', s)
        return s
    elif t == 'family':
        s = re.sub('-.*$', '', s)
        return s
    else:
        return s

d = {}

for i in range(len(args.indices)):
    f = open("mid%d_clones.csv" % args.indices[i], "r")

    # remove header
    f.readline()

    for l in f.readlines():
        fields = l.split("\t")

        # fields:
        # 0: count
        # 1: CDR3-nt
        # 2: CDR3-AA
        # 3: V gene
        # 4: J gene

        # we use CDR3-nt as key and we get the list of counts.
        # then we use the index i to increase the correct MID count

        seq = fields[1].strip()

        if args.unique:
            v_hit = fields[3].split(",")[0].strip()

            if args.exclude_j:
                key = "%s_%s" % (seq, tname(v_hit, args.type))
            elif args.cdr3nt_only:
                key = seq
            elif args.cdr3aa_only:
                seq = fields[2].strip()
                key = seq
            else:
                j_hit = fields[4].split(",")[0].strip()
                key = "%s_%s_%s" % (seq, tname(v_hit, args.type), tname(j_hit, args.type))
        else:
            key = seq

        if key not in d:
            d[key] = [0]*len(args.indices)

        d[key][i] += int(float(fields[0]))

    f.close()

od = collections.OrderedDict(sorted(d.items()))

for i in range(len(args.indices)):
    f = open("mid%d_complete.csv" % args.indices[i], "w")
    for k, v in od.items():
        f.write("%s\t%d\n" % (k, v[i]))
    f.close()

f = open("mids_counts.csv", "w")
f.write("\t")
f.write("\t".join(["MID"+str(i) for i in args.indices]))
f.write("\n")
for k, v in od.items():
    f.write(k)
    for e in v:
        f.write("\t%d" % e)
    f.write("\n")
f.close()
