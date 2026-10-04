# Extended Data Fig. 7 — PI16+ stroma at the tumour–epithelial interface

Code for Extended Data Fig. 7 of *A neurotrophic secretome from adipose-progenitor-like stromal cells promotes
prostate cancer lineage plasticity*. One panel: the spatial distribution of the adipose-progenitor-like stromal
state and of the tumour populations across all available multiplex-imaging slides.

The panel is the end of a whole imaging pipeline, not a standalone plot, so the pipeline is what is published
here: raw multiplex images → channel resolution → segmentation and per-cell measurement → per-cell population
classification → abundance and proximity, which is where the spatial maps are drawn.

## Layout and order

```
R/00_config.R              global parameters, dependency bootstrap, slide/batch table, confounding report
R/00b_io.R                 tile reader (tifffile/zarr by default, BioFormats optional)
R/01_channels.R            resolve which cycle and channel carries which marker, per slide
R/02_segment_measure.R     nuclear segmentation and per-cell intensity measurement, tiled
R/03_classify.R            per-cell population calls: ASPC, EPI_CA, NE_CA, Multi, Other
R/04_abundance_proximity.R abundance and proximity statistics -- AND the spatial maps used as the panel
R/05_visualize.R           per-biopsy QC atlas (not used in this figure)
R/06_assemble_base.R       the title / subtitle / legend layer
python/compose.py          clips the panel and merges it onto that layer, losslessly
```

Scripts 00–05 are run in order. The panel is the `pmap` object near the end of `R/04_abundance_proximity.R`,
written as `04_spatial_maps.pdf`.

```bash
Rscript R/00_config.R && Rscript R/01_channels.R && Rscript R/02_segment_measure.R \
  && Rscript R/03_classify.R && Rscript R/04_abundance_proximity.R
Rscript R/06_assemble_base.R
ED7_PANEL="<page holding the spatial maps>" python python/compose.py   # -> figures/Extended Data Fig. 7.pdf
```

The output page is **612 × 430 pt (8.5 × 6 in)**.

## What changed between the code and the published panel

Two differences, both cosmetic, both worth knowing before you re-run this:

1. **The subtitle differs.** `R/04_abundance_proximity.R` writes *"ASPC/EPI_CA/NE_CA use pixel-level corrected
   calls; blank controls are audit only"*. The panel in the figure reads *"all available slides; each panel
   translated to its own origin; shared mm scale"*. No R source file under the project or manuscript trees
   contains that second wording — checked against a control string from the same script, which does match, so
   the search is known to work. It was therefore changed downstream of the script, in the layout step or in a
   copy that is not on this machine. The search covered `.R` files only, not every file type. The title, facet
   labels (`GROUP | slide_id`), populations and colours are identical.
2. **The legend has five entries in the figure** (ASPC, EPI_CA, NE_CA, Multi, Other) where the script defines
   eight, the extra three being the `*_not_callable` audit levels. Those levels are dropped in the published
   panel.

`python/compose.py` therefore clips from the laid-out artwork by default (`ED7_BOX` defaults to the box in the
612 × 792 pt page used for the manuscript). A freshly generated `04_spatial_maps.pdf` is 12 × 8 in, so re-measure
the box and pass it as `ED7_BOX="x0,y0,x1,y1"`.

## Two things the pipeline is careful about

- **The antibody is PDGFa, the ligand (CST 5241T), not PDGFRA, the receptor.** The ASPC definition uses the
  channel that was actually stained; `PDGF_IS_LIGAND` in `R/00_config.R` records this rather than letting the
  marker name imply the receptor.
- **Marker confidence is measured against matched blank controls** at the same threshold, per slide, and reported
  alongside every per-biopsy figure.

## Paths

`ED7_ROOT` points at the raw multiplex images (default `/Volumes/T9/Case Studies`); `ED7_OUT` is where the
figure is written. `python/compose.py` reads `ED7_REPO`, `ED7_OUT`, `ED7_PANEL` and `ED7_BOX`.

## Inputs that are NOT in this repository

**The images.** The multiplex slides are large and are patient material; point `ED7_ROOT` at your own copy. The
intermediate objects the pipeline writes (`results_ASPC_NEPC/`, `logs/`) are gitignored for the same reason.

## Requirements

R >= 4.4. `R/00_config.R` bootstraps its own dependencies (`AUTO_INSTALL <- TRUE`), including the Bioconductor
imaging stack and `reticulate` for the Python `tifffile`/`zarr` reader. Set `READER <- "bioformats"` to use
RBioFormats instead, which needs Java; the default avoids it.
Python >= 3.10 with `pypdf` for the compose step (`pip install -r python/requirements.txt`), and `tifffile` plus
`zarr` for the image reader.
