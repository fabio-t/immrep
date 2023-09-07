#!/usr/bin/python3

import matplotlib
import csv
import pickle
import math
from collections import defaultdict
from ete3 import Tree, faces, TreeStyle, NodeStyle, TextFace, SequenceFace, COLOR_SCHEMES
from Bio import AlignIO
import pandas as pd

# FIXME: figure out some way to make this more generic across multiple
# experiments; for now it just uses the MID number as colour position,
# which severely limits its usage

colours = matplotlib.colormaps['tab20b'].colors + matplotlib.colormaps['tab20c'].colors + matplotlib.colormaps['Set1'].colors + matplotlib.colormaps['Set2'].colors + matplotlib.colormaps['Dark2'].colors
colours = [matplotlib.colors.rgb2hex(colour) for colour in colours]

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
            seqs[seq_id]["tot_abundance"] = [tot]
            seqs[seq_id]["percents"] = [100*i/tot for i in seqs[seq_id]["abundance"]]

tmp = pd.DataFrame.from_dict(seqs, orient="index")
tmp2 = tmp.applymap(lambda x: ",".join(map(str,x)))
tmp2.to_csv("tree.csv", index_label="node", sep="\t")

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
        cols = [colours[int(s[3:]) % len(colours)] for s in seqs[n.name]["sample"]]
        values = seqs[n.name]["percents"]
        F = faces.PieChartFace(values, colors=cols,
                               width=size * 2, height=size * 2)
        F.border.width = None
        # F.opacity = 0.6
        faces.add_face_to_node(F, n, 0, position="branch-right")
        ns = NodeStyle()
        ns["size"] = 0
        n.set_style(ns)

def layout2(n):
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

    if n.name in seqs:
        # add CDR3
        # print(n.name)
        S = SequenceFace(seqs[n.name]["cdr3nt"][0], 'nt') #codon=seqs[n.name]["cdr3nt"][0])
        faces.add_face_to_node(S, n, 2, position="aligned")
    if n.abundance > 1:
        cols = [colours[int(s[3:]) % len(colours)] for s in seqs[n.name]["sample"]]
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

for n in t.traverse("postorder"):
  # remove singleton leaves (and orphaned ancestral nodes)
  if n.abundance < 2 and n.is_leaf():
      n.detach()
ts = TreeStyle()
ts.layout_fn = layout2
ts.mode = "r"
# ts.rotation = 90
ts.show_leaf_name = False
# t.show(tree_style=ts)
ts.draw_guiding_lines = True
ts.guiding_lines_type = 2
ts.guiding_lines_color = "#CCCCCC"
t.render("tree2.svg", w=1280, tree_style=ts)
