#!/usr/bin/python3

import csv
import pickle
import math
from collections import defaultdict
from ete3 import Tree, faces, TreeStyle, NodeStyle, TextFace, COLOR_SCHEMES
from Bio import AlignIO

# FIXME: figure out some way to make this more generic across multiple
# experiments; for now it just uses the MID number as colour position,
# which severely limits its usage
colours = COLOR_SCHEMES["paired"]

aln = AlignIO.read("f.orig.fasta", "fasta")
seqs = {}
root = None
with open("idmap.csv", newline="") as f:
    for seq_id, aln_ids in csv.reader(f):
        if not root:
            root = seq_id
        else:
            seqs[seq_id] = defaultdict(list)
            for aln_id in list(set(aln_ids.split(":"))):
                seq = aln[int(aln_id)]
                meta = dict(p.split("=") for p in seq.id.split("|"))
                for k, v in meta.items():
                    if k == "abundance":
                        v = int(v)
                    seqs[seq_id][k].append(v)
            tot = sum(seqs[seq_id]["abundance"])
            seqs[seq_id]["percents"] = [100*i/tot for i in seqs[seq_id]["abundance"]]

f = "gctree.out.inference.1.p"

with open(f, "rb") as fd:
    p = pickle.load(fd)

t = p.tree

def layout(n):
    size = max(1, 10 * math.sqrt(n.abundance))
    if n.abundance > 0:
        # don't add a label for internal nodes
        T = TextFace(n.name)
        faces.add_face_to_node(T, n, 1, position="branch-right")
    elif n.is_root():
        T = TextFace(n.name)
        T.hz_align = 1
        T.rotation = -90
        faces.add_face_to_node(T, n, 0)
    if n.abundance > 1:
        cols = [colours[int(s[3:])] for s in seqs[n.name]["sample"]]
        values = seqs[n.name]["percents"]
        F = faces.PieChartFace(values, colors=cols,
                               width=size * 2, height=size * 2)
        F.border.width = None
        # F.opacity = 0.6
        faces.add_face_to_node(F, n, 0, position="branch-right")
        ns = NodeStyle()
        ns["size"] = 0
        n.set_style(ns)


ts = TreeStyle()
ts.layout_fn = layout
ts.mode = "r"
ts.rotation = 90
ts.show_leaf_name = False
# t.show(tree_style=ts)
t.render("tree.svg", w=1280, tree_style=ts)
