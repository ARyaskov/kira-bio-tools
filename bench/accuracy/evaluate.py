#!/usr/bin/env python3
"""Site- and genotype-level accuracy of a query VCF against a truth VCF.

Standard library only. Every ALT allele of every record becomes one
(chrom, pos, ref, alt) key after trimming the shared suffix and prefix
(the same canonical core the annotator uses), so `A>AT` at 100 and
`GA>GAT` at 99 match. A variant counts for a sample when that sample's
genotype carries the allele; hom/het is compared on true positives.
Records outside the confident-region BED (or the region) are ignored on
both sides.

    evaluate.py --truth truth.vcf.gz --query calls.vcf.gz \
        [--bed confident.bed] [--region chr20[:1-100]] [--min-qual Q] \
        [--truth-sample S] [--query-sample S] [--label NAME] [--json out.json]

Prints one TSV row per class (SNP, INDEL, ALL) and, with --json, writes
the numbers as JSON.
"""
import argparse
import bisect
import gzip
import json
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def canonical(pos, ref, alt):
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


class Regions:
    def __init__(self, bed=None, region=None):
        self.by_chrom = {}
        self.region = None
        if bed:
            with open(bed) as f:
                for line in f:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    c, s, e = line.split("\t")[:3]
                    self.by_chrom.setdefault(c, []).append((int(s), int(e)))
            for c in self.by_chrom:
                iv = sorted(self.by_chrom[c])
                merged = []
                for s, e in iv:
                    if merged and s <= merged[-1][1]:
                        merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                    else:
                        merged.append((s, e))
                self.by_chrom[c] = ([s for s, _ in merged], [e for _, e in merged])
        if region:
            chrom, _, span = region.partition(":")
            if span:
                a, _, b = span.partition("-")
                self.region = (chrom, int(a), int(b) if b else 1 << 60)
            else:
                self.region = (chrom, 1, 1 << 60)

    def contains(self, chrom, pos1):
        if self.region and (chrom != self.region[0] or not self.region[1] <= pos1 <= self.region[2]):
            return False
        if not self.by_chrom:
            return True
        iv = self.by_chrom.get(chrom)
        if not iv:
            return False
        starts, ends = iv
        i = bisect.bisect_right(starts, pos1 - 1) - 1
        return i >= 0 and pos1 - 1 < ends[i]


def load(path, regions, sample=None, min_qual=None, pass_only=False):
    """{(chrom, pos, ref, alt): dosage} for the chosen sample."""
    out = {}
    col = None
    with opener(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.rstrip("\r\n").split("\t")
                samples = cols[9:]
                if not samples:
                    col = None
                elif sample is None:
                    col = 9
                elif sample in samples:
                    col = 9 + samples.index(sample)
                else:
                    sys.exit(f"{path}: no sample {sample!r} (has {samples})")
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos = cols[0], int(cols[1])
            if not regions.contains(chrom, pos):
                continue
            if pass_only and cols[6] not in (".", "PASS"):
                continue
            if min_qual is not None and cols[5] != "." and float(cols[5]) < min_qual:
                continue
            if cols[4] in (".", ""):
                continue
            alts = cols[4].split(",")
            if col is not None:
                fmt = cols[8].split(":")
                gi = fmt.index("GT") if "GT" in fmt else None
                gt = cols[col].split(":")[gi] if gi is not None else "."
                alleles = [a for a in gt.replace("|", "/").split("/")]
                counts = {}
                for a in alleles:
                    if a not in (".", "", "0"):
                        counts[int(a)] = counts.get(int(a), 0) + 1
            else:
                counts = {i + 1: 1 for i in range(len(alts))}
            for ai, alt in enumerate(alts, start=1):
                if alt in ("*", "<*>", "<NON_REF>") or alt.startswith("<"):
                    continue
                n = counts.get(ai, 0)
                if n == 0:
                    continue
                key = (chrom,) + canonical(pos, cols[3], alt)
                out[key] = max(out.get(key, 0), min(n, 2))
    return out


def is_snp(key):
    return len(key[2]) == 1 and len(key[3]) == 1


def evaluate(truth, query):
    res = {}
    for cls, pred in (("SNP", is_snp), ("INDEL", lambda k: not is_snp(k)), ("ALL", lambda k: True)):
        t = {k: v for k, v in truth.items() if pred(k)}
        q = {k: v for k, v in query.items() if pred(k)}
        tp = [k for k in q if k in t]
        fp = len(q) - len(tp)
        fn = len(t) - len(tp)
        gt_ok = sum(1 for k in tp if t[k] == q[k])
        prec = len(tp) / len(q) if q else 0.0
        rec = len(tp) / len(t) if t else 0.0
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
        res[cls] = {
            "truth": len(t),
            "query": len(q),
            "tp": len(tp),
            "fp": fp,
            "fn": fn,
            "precision": prec,
            "recall": rec,
            "f1": f1,
            "gt_concordance": gt_ok / len(tp) if tp else 0.0,
        }
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--bed")
    ap.add_argument("--region")
    ap.add_argument("--min-qual", type=float)
    ap.add_argument("--pass-only", action="store_true", help="ignore query records with a FILTER")
    ap.add_argument("--truth-sample")
    ap.add_argument("--query-sample")
    ap.add_argument("--label", default="query")
    ap.add_argument("--json")
    ap.add_argument("--no-header", action="store_true")
    a = ap.parse_args()
    regions = Regions(a.bed, a.region)
    truth = load(a.truth, regions, a.truth_sample)
    query = load(a.query, regions, a.query_sample, a.min_qual, a.pass_only)
    res = evaluate(truth, query)
    if not a.no_header:
        print("label\tclass\ttruth\tquery\ttp\tfp\tfn\tprecision\trecall\tf1\tgt_concordance")
    for cls in ("SNP", "INDEL", "ALL"):
        r = res[cls]
        print(
            f"{a.label}\t{cls}\t{r['truth']}\t{r['query']}\t{r['tp']}\t{r['fp']}\t{r['fn']}\t"
            f"{r['precision']:.4f}\t{r['recall']:.4f}\t{r['f1']:.4f}\t{r['gt_concordance']:.4f}"
        )
    if a.json:
        with open(a.json, "w") as f:
            json.dump({"label": a.label, "classes": res}, f, indent=2)


if __name__ == "__main__":
    main()
