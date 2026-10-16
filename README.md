# glycocalyx-egns

Analysis code for the manuscript:

> **A transcriptionally defined endothelial glycocalyx non-resolution state links
> immune proteolysis to post-infarction fibrosis**

This repository contains every script needed to reproduce the module scores,
statistical tests, spatial analyses, communication inference, perturbation
analyses and figures reported in the manuscript. All input data are public; no
new data were generated.

---

## 1. Environment

* Python 3.12.13
* Pinned dependencies: see `requirements.txt`
* No R environment is required (the scTenifoldKnk-inspired analysis is a
  documented Python re-implementation, see Sec. 5 below)

```bash
python -m venv .venv
. .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
