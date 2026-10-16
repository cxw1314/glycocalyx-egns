#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generic GEO series-matrix parser.

Reads a *_series_matrix.txt.gz (optionally probe->gene re-annotated with a
GPL *.annot.gz SOFT file) and writes:
    <outdir>/expr.pkl   : DataFrame genes x samples (float32)
    <outdir>/meta.csv   : per-sample metadata (parsed characteristics etc.)

Usage:
    python parse_series_matrix.py --matrix <path> --outdir <dir> [--annot <GPL.annot.gz>]
"""
from __future__ import annotations

import argparse
import gzip
import io
import os
import pickle
import re

import numpy as np
import pandas as pd

TAG = re.compile(r"^(?P<key>[^:]{1,60}):\s*(?P<val>.*)$")


def parse_matrix(path: str) -> tuple[pd.DataFrame, dict]:
    meta: dict[str, list[str]] = {}
    table_rows: list[str] = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!"):
                if line.startswith("!series_matrix_table_begin"):
                    in_table = True
                    continue
                if line.startswith("!series_matrix_table_end"):
                    in_table = False
                    continue
                key, _, val = line.partition("\t")
                meta.setdefault(key.strip(), []).append(val)
            elif in_table and line:
                table_rows.append(line)
    if not table_rows:
        raise RuntimeError(f"no expression table found in {path}")
    df = pd.read_csv(io.StringIO("\n".join(table_rows)), sep="\t", index_col=0)
    df.index.name = "ID_REF"
    return df, meta


def parse_annot(path: str) -> dict[str, str]:
    """ID -> gene symbol from a GEO GPL .annot.gz SOFT file."""
    id_col = sym_col = None
    out: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin"):
                header = fh.readline().rstrip("\n").split("\t")
                id_col = header.index("ID") if "ID" in header else 0
                for cand in ("gene_assignment", "Gene Symbol", "gene_symbol",
                             "Symbol", "GENE_SYMBOL"):
                    if cand in header:
                        sym_col = header.index(cand)
                        break
                continue
            if line.startswith("!platform_table_end"):
                break
            if header is None or not line:
                continue
            parts = line.split("\t")
            if len(parts) <= max(id_col, sym_col if sym_col is not None else 0):
                continue
            pid = parts[id_col]
            raw = parts[sym_col] if sym_col is not None else ""
            sym = _first_symbol(raw)
            if sym:
                out[pid] = sym
    return out


def _first_symbol(raw: str) -> str:
    if not raw or raw in ("---", "null", ""):
        return ""
    for chunk in raw.split("///"):
        fields = [f.strip() for f in chunk.split("//")]
        if len(fields) >= 2 and fields[1] and fields[1] != "---":
            cand = fields[1]
            # HTA-style can be "SYMBOL" or "SYMBOL1 /// SYMBOL2"
            cand = cand.split("//")[0].strip()
            if re.fullmatch(r"[A-Za-z0-9._\-@]{1,40}", cand):
                return cand.upper()
    return ""


def collapse_to_genes(df: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    sym = pd.Series({i: mapping.get(i, "") for i in df.index})
    keep = sym[sym != ""]
    sub = df.loc[keep.index]
    # keep the probe with the highest mean intensity per symbol
    sub = sub.assign(__sym=keep.values, __mean=sub.mean(axis=1).values)
    sub = sub.sort_values("__mean", ascending=False).drop_duplicates("__sym")
    sub = sub.drop(columns=["__mean"]).set_index("__sym")
    sub.index.name = "gene"
    return sub


def _flatten(meta: dict, key: str) -> list[str]:
    """Sample-level SOFT fields hold all samples on ONE tab-separated line."""
    out: list[str] = []
    for line in meta.get(key, []):
        out.extend(line.split("\t"))
    return [v.strip().strip('"').strip("'") for v in out]


def build_meta(meta: dict, samples: list[str]) -> pd.DataFrame:
    n = len(samples)
    cols: dict[str, list] = {}
    for key in ("!Sample_geo_accession", "!Sample_title", "!Sample_source_name_ch1",
                "!Sample_organism_ch1", "!Sample_description", "!Sample_platform_id",
                "!Sample_relation"):
        vals = _flatten(meta, key)
        if vals:
            cols[key.lstrip("!")] = (vals + [""] * n)[:n]

    # each characteristics line is ONE characteristic key across all samples
    char_keys: dict[str, list[str]] = {}
    for line in meta.get("!Sample_characteristics_ch1", []):
        cells = line.split("\t")
        for i, cell in enumerate(cells):
            m = TAG.match(cell.strip())
            if not m:
                continue
            k = m.group("key").strip()
            char_keys.setdefault(k, [""] * n)
            if i < n:
                char_keys[k][i] = m.group("val").strip()
    cols.update({f"ch_{k}": v for k, v in char_keys.items()})

    accs = cols.get("Sample_geo_accession") or samples
    if len(accs) != n:
        accs = samples
    df = pd.DataFrame(cols, index=accs)
    df.index.name = "gsm"
    return df


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--annot", default=None)
    ap.add_argument("--log2", action="store_true", help="apply log2(x+1) if values look linear")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print(f"[read] {args.matrix}")
    raw, meta = parse_matrix(args.matrix)
    print(f"       probes={raw.shape[0]} samples={raw.shape[1]}")

    samples = [str(c) for c in raw.columns]
    mdf = build_meta(meta, samples)
    mdf.to_csv(os.path.join(args.outdir, "meta.csv"), encoding="utf-8-sig")
    print(f"[meta] {mdf.shape[0]} rows x {mdf.shape[1]} cols -> meta.csv")
    if "title" in mdf.columns:
        print("       first titles:", list(mdf["title"][:4]))

    df = raw
    if args.annot:
        print(f"[annot] {args.annot}")
        mp = parse_annot(args.annot)
        print(f"        mapped probe ids: {len(mp)}")
        if mp:
            df = collapse_to_genes(df, mp)
            print(f"        genes x samples = {df.shape}")

    vals = df.to_numpy(dtype=np.float64)
    finite = np.isfinite(vals)
    sample_vals = vals[finite][:200000]
    if sample_vals.size and sample_vals.max() > 60:
        print("[log2] values look linear -> log2(x+1)")
        df = np.log2(df.clip(lower=0) + 1)
    else:
        print("[log2] values already log-scale -> kept as-is")

    df = df.astype(np.float32)
    with open(os.path.join(args.outdir, "expr.pkl"), "wb") as fh:
        pickle.dump(df, fh)
    print(f"[ok] expr.pkl {df.shape} -> {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
