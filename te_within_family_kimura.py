#!/usr/bin/env python3
"""
Within-family pairwise Kimura (K2P / K80) divergence for each annotated TE copy.

For every annotated TE copy this reports the AVERAGE Kimura distance of all
pairwise comparisons between that copy and the other copies of the SAME family
(rather than the copy-vs-consensus distance).

Strategy (consensus-anchored pairwise distances):
  1. For each TE locus, extract its genomic sequence. Composite / nested loci
     (underscore-joined family labels, i.e. "roo_Doc") are split into their
     constituent families and each constituent is treated as one copy of that
     family.
  2. Each copy is locally aligned (Smith-Waterman, both strands) to its family
     consensus and PROJECTED into consensus coordinates: for every consensus
     position we record the aligned copy base (insertions relative to the
     consensus are dropped, deletions/uncovered positions are marked missing).
     This puts all copies of a family into a common coordinate frame and a
     common orientation without a separate multiple-sequence alignment.
  3. Within each family, every pair of copies is compared over the consensus
     positions where both have a base (pairwise deletion). Transitions (P) and
     transversions (Q) are counted and the K2P distance is
         K = -1/2 ln(1 - 2P - Q) - 1/4 ln(1 - 2Q)   (Kimura 1980).
  4. For each copy we report the mean / median / sd of its pairwise K2P
     distances to all other members of the family (that share >= --min-overlap
     aligned consensus positions).

Family grouping scope is controlled by --group-scope:
  assembly (default) : copies grouped within each assembly separately
  pooled             : copies of a family pooled across all assemblies

Only Biopython + numpy are required (no external aligner binaries).

NOTE: the annotation files contain only the subset of TEs used in the study
(e.g. interacting TEs), so "within-family" comparisons are among those copies,
not every copy of the family in the genome.
"""

import argparse
import math
import os
import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

COL_CHROM, COL_START, COL_END, COL_FAMILY, COL_CLASS = 22, 23, 24, 25, 26

# base -> code; anything else (N, gap) -> -1
CODE = np.full(256, -1, dtype=np.int8)
for b, c in zip(b"ACGT", (0, 1, 2, 3)):
    CODE[b] = c


def build_aligner():
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 1.0
    a.mismatch_score = -1.0
    a.open_gap_score = -5.0
    a.extend_gap_score = -0.5
    return a


def load_library(path):
    lib = {}
    for rec in SeqIO.parse(path, "fasta"):
        lib[rec.id.split("#")[0]] = str(rec.seq).upper()
    lower = {k.lower(): k for k in lib}
    return lib, lower


def resolve_family(fam, lib, lower):
    if fam in lib:
        return fam
    if fam.endswith("-int") and fam[:-4] in lib:
        return fam[:-4]
    if fam.lower() in lower:
        return lower[fam.lower()]
    if fam.endswith("-int") and fam[:-4].lower() in lower:
        return lower[fam[:-4].lower()]
    return None


def load_loci(path):
    seen, loci = set(), []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) <= COL_CLASS:
                continue
            chrom = f[COL_CHROM]
            try:
                start, end = int(f[COL_START]), int(f[COL_END])
            except ValueError:
                continue
            key = (chrom, start, end, f[COL_FAMILY])
            if key in seen:
                continue
            seen.add(key)
            loci.append(dict(chrom=chrom, start=start, end=end,
                             family=f[COL_FAMILY], cls=f[COL_CLASS]))
    return loci


def best_alignment(aligner, copy_seq, cons_seq):
    """Return (alignment, strand, score, target_used) for the better strand."""
    fwd = aligner.align(copy_seq, cons_seq)
    rc = str(Seq(copy_seq).reverse_complement())
    rev = aligner.align(rc, cons_seq)
    if rev.score > fwd.score:
        return rev[0], "-", rev.score, rc
    return fwd[0], "+", fwd.score, copy_seq


def project_to_consensus(aln, target_used, cons_len):
    """Project an aligned copy into consensus coordinates as an int8 array
    (length cons_len; -1 = missing/uncovered). target = copy(or rc), query =
    consensus, so aln.aligned blocks advance 1:1 within each block."""
    proj = np.full(cons_len, -1, dtype=np.int8)
    tcodes = CODE[np.frombuffer(target_used.encode(), dtype=np.uint8)]
    tblocks, qblocks = aln.aligned  # target blocks, query(=consensus) blocks
    for (ts, te), (qs, qe) in zip(tblocks, qblocks):
        proj[qs:qe] = tcodes[ts:te]
    return proj


def kimura_pair(a, b):
    """K2P distance between two consensus-projected int8 arrays over shared
    non-missing positions. Returns (n_shared, k2p_or_None, p_distance_or_None)."""
    mask = (a >= 0) & (b >= 0)
    n = int(mask.sum())
    if n == 0:
        return 0, None, None
    aa, bb = a[mask], b[mask]
    diff = aa != bb
    # transitions: A<->G (0,2) and C<->T (1,3)
    s = aa + bb
    d = np.abs(aa.astype(np.int16) - bb.astype(np.int16))
    ts = int(np.count_nonzero(diff & (d == 2)))  # 0<->2 and 1<->3 both differ by 2
    tv = int(np.count_nonzero(diff) - ts)
    P, Q = ts / n, tv / n
    w1, w2 = 1 - 2 * P - Q, 1 - 2 * Q
    k2p = (-0.5 * math.log(w1) - 0.25 * math.log(w2)) if (w1 > 0 and w2 > 0) else None
    return n, k2p, (ts + tv) / n


PER_COPY_HEADER = [
    "assembly", "locus_id", "chrom", "start", "end", "locus_len",
    "locus_family", "locus_class", "n_families", "is_composite",
    "family_index", "family", "family_class", "consensus_name", "consensus_len",
    "strand", "cov_consensus", "kimura_vs_consensus",
    "group_scope", "n_family_copies", "n_pairs_used",
    "mean_pairwise_k2p", "median_pairwise_k2p", "sd_pairwise_k2p", "note",
]


def fmt(x, nd=6):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def collect_copies(assembly, genome_path, annotation_path, lib, lower, aligner,
                   exclude_composite=False):
    print(f"[{assembly}] loading genome {os.path.basename(genome_path)} ...", flush=True)
    genome = {r.id: str(r.seq).upper() for r in SeqIO.parse(genome_path, "fasta")}
    loci = load_loci(annotation_path)
    if exclude_composite:
        n_before = len(loci)
        loci = [L for L in loci if "_" not in L["family"]]
        print(f"[{assembly}] excluded {n_before - len(loci)} nested/composite "
              f"loci; {len(loci)} single-family loci remain", flush=True)
    else:
        print(f"[{assembly}] {len(loci)} unique TE loci", flush=True)
    copies = []
    for i, L in enumerate(loci, 1):
        locus_id = f"{assembly}_{L['chrom']}_{L['start']}_{L['end']}"
        if L["chrom"] not in genome:
            copies.append(dict(assembly=assembly, locus_id=locus_id, L=L,
                               locus_len=None, fi=1, n_fam=1, fam=L["family"],
                               fam_cls=L["cls"], cons_key=None, cons_len=None,
                               strand="NA", cov=None, k_cons=None, proj=None,
                               note="chrom_not_in_assembly"))
            continue
        copy_seq = genome[L["chrom"]][L["start"] - 1:L["end"]]
        locus_len = len(copy_seq)
        fams = L["family"].split("_")
        clss = L["cls"].split("_")
        n_fam = len(fams)
        for fi, fam in enumerate(fams, 1):
            fam_cls = clss[fi - 1] if fi - 1 < len(clss) else "NA"
            cons_key = resolve_family(fam, lib, lower)
            if cons_key is None:
                copies.append(dict(assembly=assembly, locus_id=locus_id, L=L,
                                   locus_len=locus_len, fi=fi, n_fam=n_fam,
                                   fam=fam, fam_cls=fam_cls, cons_key=None,
                                   cons_len=None, strand="NA", cov=None,
                                   k_cons=None, proj=None,
                                   note="family_not_in_library"))
                continue
            cons = lib[cons_key]
            aln, strand, score, tgt = best_alignment(aligner, copy_seq, cons)
            proj = project_to_consensus(aln, tgt, len(cons))
            n_cov = int((proj >= 0).sum())
            # consensus-anchored copy-vs-consensus K2P (for reference column)
            cons_codes = CODE[np.frombuffer(cons.encode(), dtype=np.uint8)]
            _, k_cons, _ = kimura_pair(proj, cons_codes)
            copies.append(dict(assembly=assembly, locus_id=locus_id, L=L,
                               locus_len=locus_len, fi=fi, n_fam=n_fam, fam=fam,
                               fam_cls=fam_cls, cons_key=cons_key,
                               cons_len=len(cons), strand=strand,
                               cov=(n_cov / len(cons) if len(cons) else None),
                               k_cons=k_cons, proj=proj, note=""))
        if i % 50 == 0:
            print(f"[{assembly}] {i}/{len(loci)} loci aligned", flush=True)
    return copies


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--run", action="append", nargs=3,
                    metavar=("NAME", "GENOME_FASTA", "ANNOTATION_TXT"),
                    required=True)
    ap.add_argument("--group-scope", choices=["assembly", "pooled"],
                    default="assembly",
                    help="group family copies within each assembly (default) "
                         "or pooled across all assemblies")
    ap.add_argument("--min-overlap", type=int, default=100,
                    help="min shared consensus positions for a pair to count "
                         "(default 100)")
    ap.add_argument("--exclude-composite", action="store_true",
                    help="drop nested/composite loci (underscore-joined family "
                         "labels, e.g. I-element_roo) entirely, so they are "
                         "neither measured nor used as comparison partners")
    args = ap.parse_args()

    lib, lower = load_library(args.library)
    print(f"library: {len(lib)} family consensus sequences", flush=True)
    aligner = build_aligner()

    all_copies = []
    for name, genome_path, ann_path in args.run:
        all_copies.extend(collect_copies(name, genome_path, ann_path,
                                         lib, lower, aligner,
                                         exclude_composite=args.exclude_composite))

    # group copies by (scope, family) using only alignable copies
    groups = {}
    for idx, c in enumerate(all_copies):
        if c["proj"] is None:
            continue
        scope = "pooled" if args.group_scope == "pooled" else c["assembly"]
        groups.setdefault((scope, c["fam"]), []).append(idx)

    # pairwise K2P within each group
    print(f"computing within-family pairwise distances "
          f"(scope={args.group_scope}, {len(groups)} family groups) ...", flush=True)
    mean_k = [None] * len(all_copies)
    med_k = [None] * len(all_copies)
    sd_k = [None] * len(all_copies)
    npairs = [0] * len(all_copies)
    gsize = [0] * len(all_copies)
    per_copy_vals = {idx: [] for g in groups.values() for idx in g}

    for (scope, fam), idxs in groups.items():
        m = len(idxs)
        for idx in idxs:
            gsize[idx] = m
        for ii in range(m):
            for jj in range(ii + 1, m):
                i, j = idxs[ii], idxs[jj]
                n, k, _ = kimura_pair(all_copies[i]["proj"], all_copies[j]["proj"])
                if n < args.min_overlap or k is None:
                    continue
                per_copy_vals[i].append(k)
                per_copy_vals[j].append(k)
    for idx, vals in per_copy_vals.items():
        npairs[idx] = len(vals)
        if vals:
            mean_k[idx] = float(np.mean(vals))
            med_k[idx] = float(np.median(vals))
            sd_k[idx] = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0

    with open(args.out, "w") as out:
        out.write("\t".join(PER_COPY_HEADER) + "\n")
        for idx, c in enumerate(all_copies):
            L = c["L"]
            scope = "pooled" if args.group_scope == "pooled" else c["assembly"]
            row = [
                c["assembly"], c["locus_id"], L["chrom"], str(L["start"]),
                str(L["end"]), fmt(c["locus_len"], 0) if c["locus_len"] else "NA",
                L["family"], L["cls"], str(c["n_fam"]),
                "1" if c["n_fam"] > 1 else "0", str(c["fi"]), c["fam"],
                c["fam_cls"], c["cons_key"] or "NA",
                str(c["cons_len"]) if c["cons_len"] else "NA", c["strand"],
                fmt(c["cov"]), fmt(c["k_cons"]),
                scope if c["proj"] is not None else "NA",
                str(gsize[idx]) if c["proj"] is not None else "NA",
                str(npairs[idx]) if c["proj"] is not None else "NA",
                fmt(mean_k[idx]), fmt(med_k[idx]), fmt(sd_k[idx]),
                c["note"] or ("singleton_no_pairs"
                              if (c["proj"] is not None and npairs[idx] == 0)
                              else ""),
            ]
            out.write("\t".join(row) + "\n")
    print(f"wrote {len(all_copies)} rows -> {args.out}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
