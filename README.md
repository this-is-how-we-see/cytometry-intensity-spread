# Intensity spread is information

Code, per-run reports and derived results for the preprint *Intensity spread is information: what FlowSOM
loses when cytometry data are scaled, binarized or brightness-normalized* (Howie Seay, HOW WE SEE LLC,
ORCID [0000-0001-7002-0187](https://orcid.org/0000-0001-7002-0187)). bioRxiv DOI: to be added.

The study clusters the same cells with FlowSOM from several representations of the data (arcsinh
intensity, z-scored and percentile-scaled intensity, brightness-normalized intensity, lineage markers as
bits, every marker as a bit, reweighted intensity) on six public benchmark datasets with manual gates and
the seven-center HIPC Lyoplate T cell panel, and scores recovery, stability and between-site spread.
Plain arcsinh intensity wins; the reports say where and why each alternative loses.

## Layout

- `manuscript/` the preprint (`preprint.md`, `preprint.pdf`) and its figures
- `results/` one report per run (`*.md`) with the metrics behind it (`*_metrics.csv`, `*_perpop.csv`, `*_meta.json`)
- `code/pilot.py` the representations, FlowSOM runs and scoring; `--datasets`, `--arms` and `--tag` select a run
- `code/` the data extraction, the HIPC loader and transform check, the FlowJo 9 gate application, the between-site analysis, and the figure, table and PDF builders
- `data/README.md` every source, with checksums; the data themselves are not in this repository

## Reproduce

```
python3 -m venv code/.venv && . code/.venv/bin/activate && pip install -e code flowutils markdown
python3 -I code/mirror_zzph.py data/FR-FCM-ZZPH
python3 -I code/rda_extract.py data/HDCytoData/Nilsson_rare_SE.rda data/HDCytoData/Nilsson_rare.npz
python3 -I code/pilot.py . --tag pilot
```

`data/README.md` lists the HDCytoData files and the ImmuneSpace folders. The HIPC files need a free
ImmuneSpace account. Each report in `results/` names the command that produced it.

## License

Code: MIT (see `LICENSE`). Reports, figures and the manuscript: CC-BY 4.0.
