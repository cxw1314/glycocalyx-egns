#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Search NCBI GEO (via eutils, stdlib only) for human cardiovascular disease
expression datasets. Outputs a candidate table for manual/automatic curation.

Usage:
    python geo_search.py --out ../results/tables/geo_candidates.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "glycocalyx_cvd"
EMAIL = "researcher@example.org"

# Disease category -> GEO title-level query
QUERIES = {
    "Atherosclerosis_CAD": (
        '(atherosclerosis[Title] OR atherosclerotic[Title] OR "coronary artery disease"[Title] '
        'OR "coronary heart disease"[Title] OR "myocardial infarction"[Title] '
        'OR "coronary atherosclerosis"[Title] OR "carotid plaque"[Title])'
    ),
    "HeartFailure": (
        '("heart failure"[Title] OR "dilated cardiomyopathy"[Title] OR "ischemic cardiomyopathy"[Title] '
        'OR "cardiac failure"[Title])'
    ),
    "Hypertension_Hypertrophy": (
        '(hypertension[Title] OR hypertensive[Title] OR "cardiac hypertrophy"[Title] '
        'OR "left ventricular hypertrophy"[Title] OR "myocardial hypertrophy"[Title])'
    ),
    "AtrialFibrillation": (
        '("atrial fibrillation"[Title] OR "atrial myocardium"[Title] OR "atrial remodeling"[Title])'
    ),
    "Aortic_Vascular": (
        '("aortic aneurysm"[Title] OR "aortic dissection"[Title] OR "aortic stenosis"[Title] '
        'OR "peripheral artery disease"[Title] OR "vascular injury"[Title] OR "aortic valve"[Title])'
    ),
    "Myocarditis": (
        '(myocarditis[Title] OR "myocardial inflammation"[Title] OR "inflammatory cardiomyopathy"[Title] '
        'OR "viral myocarditis"[Title])'
    ),
}


def _get(url: str, params: dict, retries: int = 4, timeout: int = 45) -> bytes:
    q = urllib.parse.urlencode(params)
    full = f"{url}?{q}"
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(full, headers={"User-Agent": TOOL})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"request failed: {full} :: {last}")


def esearch(term: str, retmax: int = 200) -> list[str]:
    params = {
        "db": "gds",
        "term": term,
        "retmax": str(retmax),
        "retmode": "json",
        "tool": TOOL,
        "email": EMAIL,
    }
    data = json.loads(_get(f"{EUTILS}/esearch.fcgi", params).decode("utf-8", "replace"))
    time.sleep(0.4)
    return data.get("esearchresult", {}).get("idlist", [])


def esummary(ids: list[str], chunk: int = 100) -> dict:
    out: dict = {}
    for i in range(0, len(ids), chunk):
        part = ids[i:i + chunk]
        params = {
            "db": "gds",
            "id": ",".join(part),
            "retmode": "json",
            "tool": TOOL,
            "email": EMAIL,
        }
        data = json.loads(_get(f"{EUTILS}/esummary.fcgi", params).decode("utf-8", "replace"))
        out.update(data.get("result", {}))
        time.sleep(0.5)
    return out


def norm_taxon(taxon) -> str:
    if isinstance(taxon, list):
        return ";".join(str(t) for t in taxon)
    return str(taxon or "")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--retmax", type=int, default=200)
    ap.add_argument("--min-samples", type=int, default=6)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    rows = []
    for cat, core in QUERIES.items():
        term = f'{core} AND "Homo sapiens"[Organism] AND gse[Entry Type]'
        print(f"[search] {cat} ...", flush=True)
        ids = esearch(term, args.retmax)
        print(f"         -> {len(ids)} hits", flush=True)
        if not ids:
            continue
        summ = esummary(ids)
        for gid in ids:
            rec = summ.get(gid)
            if not isinstance(rec, dict):
                continue
            acc = rec.get("accession", "") or ""
            if not acc.startswith("GSE"):
                continue
            try:
                n = int(rec.get("n_samples", 0) or 0)
            except Exception:  # noqa: BLE001
                n = 0
            taxon = norm_taxon(rec.get("taxon"))
            if "Homo sapiens" not in taxon:
                continue
            if n < args.min_samples:
                continue
            pmid = rec.get("pubmedids") or []
            if isinstance(pmid, list):
                pmid = ";".join(str(p) for p in pmid[:3])
            s = rec.get("summary") or ""
            s = " ".join(str(s).split())
            rows.append({
                "category": cat,
                "gse": acc,
                "title": " ".join(str(rec.get("title", "")).split()),
                "n_samples": n,
                "platform_gpl": rec.get("gpl", ""),
                "type": rec.get("gdstype", rec.get("gdsType", "")),
                "pubmed": pmid,
                "summary": s[:400],
            })

    rows.sort(key=lambda r: (r["category"], -r["n_samples"]))
    with open(args.out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                           ["category", "gse", "title", "n_samples", "platform_gpl",
                            "type", "pubmed", "summary"])
        w.writeheader()
        w.writerows(rows)
    print(f"[done] {len(rows)} candidate datasets -> {args.out}", flush=True)
    for cat in QUERIES:
        sub = [r for r in rows if r["category"] == cat]
        print(f"  {cat}: {len(sub)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
