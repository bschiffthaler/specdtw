import os
import logging
from networkit.graph import Graph
from collections import Counter, defaultdict
import numpy as np
import networkit as nk
from copy import deepcopy
import pickle
from numba import jit
import progressbar


def isclose(a, b, rel_tol=1e-9, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class peak(object):
    """docstring for peak"""

    def __init__(
        self,
        index: int,
        ppm: float,
        value: float,
        snr: float,
        scale: float,
        sample: int,
    ):
        super(peak, self).__init__()
        self.index = index
        self.ppm = ppm
        self.value = value
        self.snr = snr
        self.scale = scale
        self.sample = sample

    def __repr__(self):
        return f"peak({self.index}, {self.ppm}, {self.value}, {self.snr}, "
        f"{self.scale}, {self.sample})"


class plate(object):
    """docstring for plate"""

    def __init__(self, file: str, platenum: int):
        super(plate, self).__init__()
        self.file = file
        self.peak_dict = {}
        self.idx_to_sample = {}
        self.platenum = platenum
        logging.debug(f"Parsing {file}...")
        with open(file, "r") as peak_file:
            for line in peak_file:
                index, ppm, value, snr, scale, sample = line.strip().split()
                p = peak(
                    int(index),
                    float(ppm),
                    float(value),
                    float(snr),
                    float(scale),
                    int(sample),
                )
                self.peak_dict[(platenum, int(sample), int(index))] = p
                self.idx_to_sample[int(index)] = int(sample)
        logging.debug(f"Have {len(self.peak_dict)} peaks.")

    def __getitem__(self, item):
        return self.peak_dict[item]

    def __iter__(self):
        self.iter = iter(self.peak_dict)
        return self

    def __next__(self):
        return next(self.iter)


class match(object):
    """docstring for match"""

    def __init__(
        self, dist: str, rownames: str, colnames: str, plate_r: plate, plate_c: plate
    ):
        super(match, self).__init__()
        self.dist = dist
        self.rownames = []
        self.colnames = []
        self.dists = []
        self.first, self.second = os.path.basename(dist).replace(".d", "").split("_")
        logging.debug(f"Parsing match: {(dist, rownames, colnames)}...")
        with open(rownames, "r") as in_r:
            for line in in_r:
                x = int(line.strip())
                s = plate_r.idx_to_sample[x]
                self.rownames.append((plate_r.platenum, s, x))

        with open(colnames, "r") as in_r:
            for line in in_r:
                x = int(line.strip())
                s = plate_c.idx_to_sample[x]
                self.colnames.append((plate_c.platenum, s, x))

        with open(dist, "r") as in_r:
            for line in in_r:
                vals = line.strip().split()
                r, c, w = (int(vals[0]), int(vals[1]), float(vals[2]))
                u = (plate_r.platenum, plate_r.idx_to_sample[r], r)
                v = (plate_c.platenum, plate_c.idx_to_sample[c], c)
                self.dists.append((u, v, w))
        logging.debug(f"Have {len(self.dists)} matches")

    def __getitem__(self, item):
        return self.dists[item]

    def __iter__(self):
        self.iter = iter(self.dists)
        return self

    def __next__(self):
        return next(self.iter)


def prune(ndict, rdict, graph, node):
    """
    Prune all except strongest between-sample connections
        Given a peak P and its neighbors N, we first group N by sample ID. For all
        groups with > 1 members (i.e. multiuple peaks originating from the same
        sample id) we select the strongest link and delete all others
    """
    nn = [rdict[n] for n in graph.iterNeighbors(node)]
    dd = defaultdict(list)
    for n in nn:
        a = n[0]
        b = n[1]
        # Make order consistent
        if b < a:
            tmp = a
            a = b
            b = tmp
        dd[f"{a}_{b}"].append(n)
    for k in dd:
        v = dd[k]
        # Check for sample multiplicity
        if len(v) > 1:
            weights = []
            edge_ids = []
            for other_node in v:
                e_id = ndict[other_node]
                e_w = gr.weight(node, e_id)
                weights.append(e_w)
                edge_ids.append(e_id)
            argmax = np.argmax(weights)
            edge_ids.pop(argmax)
            for edge_id in edge_ids:
                gr.removeEdge(node, edge_id)


def avg_membership(n: int, group: list, graph: nk.Graph):
    weights = []
    for node in group:
        if graph.hasEdge(n, node):
            weights.append(graph.weight(n, node))
        else:
            weights.append(0.0)
    if weights:
        return np.mean(weights)
    else:
        return 0


def best_group(
    n: tuple,
    nd: dict,
    groups: dict,
    graph: nk.Graph,
    min_score: int = 0,
    ret_weights=False,
):
    best_weights = []
    best_indices = []
    n_id = nd[n]
    for new_group in groups:
        grp_members = [nd[x] for x in groups[new_group]]
        w = avg_membership(n_id, grp_members, graph)
        if w > min_score:
            best_weights.append(w)
            best_indices.append(new_group)
    if best_weights:
        am = np.argmax(best_weights)
        if ret_weights:
            return (best_indices[am], best_weights[am], best_weights)
        return (best_indices[am], best_weights[am])
    if ret_weights:
        return (-1, -1, -1)
    return (-1, -1)


def initial_grouping(plates: dict, node_dict: dict, gr: nk.Graph):
    sample_nodes = defaultdict(list)
    for pl in plates:
        for pk in plates[pl]:
            cur_pk = plates[pl][pk]
            sam = (pl, cur_pk.sample)
            sample_nodes[sam].append(pk)

    max_grp = 0
    groups = defaultdict(list)
    for sm in progressbar.progressbar(sample_nodes):
        for pk in sample_nodes[sm]:
            r = best_group(pk, node_dict, groups, gr)
            if r == (-1, -1):
                groups[max_grp].append(pk)
                max_grp += 1
            else:
                groups[r[0]].append(pk)
    return groups


def remove_dupes(group: list[(int, int, int)], node_dict: dict, gr: nk.Graph):
    sam_table = defaultdict(list)
    duplicated = list()
    g_nodes = [node_dict[g] for g in group]
    drops = list()
    for p in group:
        k = (p[0], p[1])
        sam_table[k].append(p)
    for s in sam_table:
        if len(sam_table[s]) > 1:
            duplicated.append(sam_table[s])
    for dupe in duplicated:
        d_nodes = [node_dict[g] for g in dupe]
        d_weights = [avg_membership(d, g_nodes, gr) for d in d_nodes]
        am = np.argmax(d_weights)
        for i in range(len(dupe)):
            if i != am:
                drops.append(dupe[i])
    for drop in drops:
        group.remove(drop)
    return (group, drops)


def prune_groups(groups: dict, node_dict: dict, gr: nk.Graph):
    all_dropped = list()
    for gi in progressbar.progressbar(groups):
        new_grp, dropped = remove_dupes(groups[gi], node_dict, gr)
        all_dropped.append(dropped)
    all_dropped = [i for x in all_dropped for i in x]
    return all_dropped


def assign_new(dropped: list, groups: dict, node_dict: dict, gr: nk.Graph):
    for d in progressbar.progressbar(dropped):
        ng, w = best_group(d, node_dict, groups, gr, 0)
        if ng == -1:
            m = max([k for k in groups]) + 1
            groups[m].append(d)
        else:
            groups[ng].append(d)


def group_sim(x: int, y: int, groups: dict, node_dict: dict, gr: nk.Graph):
    xi = [node_dict[i] for i in groups[x]]
    yi = [node_dict[i] for i in groups[y]]
    weights = []
    for i in xi:
        for j in yi:
            if gr.hasEdge(i, j):
                weights.append(gr.weight(i, j))
            else:
                weights.append(0.0)
    if weights:
        return np.mean(weights)
    return 0.0


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)

    gr = Graph(n=0, weighted=True, directed=False)
    node_dict = {}
    plates = {}

    """
    Iterate through all the plate peak definitions and insert a disconnected node
    into the graph.
    """
    logging.debug("Setting up graph")
    for root, path, file in os.walk("prep"):
        for f in file:
            if f.startswith("g") and f.endswith("tsv"):
                idx = int(f.replace("g", "").replace(".tsv", ""))
                p = plate(os.path.join(root, f), idx)
                plates[p.platenum] = p
                for node in p:
                    cur_node = gr.addNode()
                    node_dict[node] = cur_node

    # Reverse lookup node id -> node name
    rev_node_dict = dict([(v, k) for (k, v) in node_dict.items()])
    logging.debug("Connecting peaks with specDTW weights")
    for root, path, file in os.walk("out"):
        for f in file:
            if f.endswith(".d"):
                rows = os.path.join(root, f.replace(".d", ".r"))
                cols = os.path.join(root, f.replace(".d", ".c"))
                dist = os.path.join(root, f)

                ir, ic = [int(x) for x in f.replace(".d", "").split("_")]

                m = match(dist, rows, cols, plates[ir], plates[ic])
                for edge in m:
                    src, tgt, weight = (node_dict[edge[0]], node_dict[edge[1]], edge[2])
                    gr.addEdge(src, tgt, weight)
                del m

    # Now we prune the graph and keep only the strongest between-sample
    # connections
    logging.debug("Pruning weak connections")
    i = 0
    for n in gr.iterNodes():
        prune(node_dict, rev_node_dict, gr, n)
        i += 1
        if i % 1000 == 0:
            print(i, gr)

    # Some more cleanup post-hoc
    logging.debug("Pruning self-loops")
    gr.removeSelfLoops()
    logging.debug("Pruning remaining multi edges")
    gr.removeMultiEdges()
    logging.debug("Compacting")
    gr.compactEdges()

    # Save graph
    with open("graph.dicts", "wb") as oo:
        pickle.dump({"nd": node_dict, "rd": rev_node_dict, "pl": plates}, oo)
    nk.graphio.writeGraph(gr, "graph.nkbin", nk.graphio.Format(19))

    groups = initial_grouping(plates, node_dict, gr)
    with open("groups.pck", "wb") as oo:
        pickle.dump(groups, oo)
    with open("groups.pck", "rb") as ii:
        groups = pickle.load(ii)

    for _ in range(10):
        logging.info(f"Refine round {_}")
        dropped = prune_groups(groups, node_dict, gr)
        logging.info(f"Need to reassign {len(dropped)} peaks")
        assign_new(dropped, groups, node_dict, gr)

    with open("initial_grouping.txt", "w") as igrp:
        for i in range(len(groups)):
            for j, x in enumerate(groups[i]):
                igrp.write(str(x))
                if j == (len(groups[i]) - 1):
                    igrp.write("\n")
                else:
                    igrp.write(" ")

    with open("initial_grouping_e.txt", "w") as igrp:
        for i, j, k in gr.iterEdgesWeights():
            igrp.write(f"{i} {j} {k}\n")
