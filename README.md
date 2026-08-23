# TMT proteomics analysis workshop

A two-part R/Quarto workshop that teaches TMT proteomics analysis end to end:
import and quality control, internal reference scaling (IRS) across plexes,
unsupervised clustering, signature and pathway analysis, and a confounder-aware
group comparison. The dataset is synthetic, built to reproduce the statistical
structure of a real TMT experiment (batch effects, blockwise missingness,
plex-subtype confounding) without containing any real patient, disease or
protein data.

Everything needed to run it, including the data, is in this repository.

## Quick start

```bash
git clone https://github.com/YOUR_USERNAME/pedPMS-proteomics-workshop.git
cd pedPMS-proteomics-workshop
```

Open `pedPMS-proteomics-workshop.Rproj` in RStudio, install the packages listed
in `WORKSHOP_GUIDE.md`, then render both scripts in order:

```r
quarto::quarto_render("pedPMS_1_Import_QC_IRS.qmd")
quarto::quarto_render("pedPMS_2_Exploration.qmd")
```

**Read [`WORKSHOP_GUIDE.md`](WORKSHOP_GUIDE.md) first.** It covers setup,
the input files, the biology and TMT concepts behind the analysis, what each
script does step by step, expected outputs, troubleshooting, and how the
synthetic dataset was constructed.

## Contents

```
pedPMS-proteomics-workshop/
├── pedPMS_1_Import_QC_IRS.qmd   # import, QC, reference-channel IRS
├── pedPMS_2_Exploration.qmd     # clustering, signatures, enrichment, plex-7 comparison
├── WORKSHOP_GUIDE.md            # full written guide (read this first)
├── data_raw/                    # synthetic TMT dataset (quantitative + metadata)
└── results/                     # created when you render; not tracked by git
```

## License

The analysis code is provided for teaching purposes. The dataset is fully
synthetic (see section 11 of `WORKSHOP_GUIDE.md` for how it was derived and
anonymised) and contains no real patient, disease, or protein data.
