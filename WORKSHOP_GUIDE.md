# TMT proteomics analysis workshop: guide

This guide explains everything needed to run the two analysis scripts in this
workshop and to understand what they do. It is written for someone new to
proteomics data analysis. Read it once before you start, then keep it open next
to the scripts.

The workshop uses a synthetic dataset built to look and behave like a real TMT
proteomics experiment. Nothing in it identifies a real patient, disease or
protein, but its numerical structure, its quirks and its analytical challenges
are all realistic, so the skills transfer directly to real data.

## 1. What you will do

You will take a raw table of protein intensities and a couple of sample
spreadsheets, turn them into a clean dataset, check its quality, correct a batch
effect between measurement runs, group the samples by their protein profiles, and
finish with a focused comparison of healthy against disease samples. Two scripts
do this in order:

Script 1, `pedPMS_1_Import_QC_IRS.qmd`, reads the files, builds the sample table,
runs quality control, performs internal reference scaling (IRS), and saves one
object.

Script 2, `pedPMS_2_Exploration.qmd`, reads that object, clusters the samples,
finds and interprets the proteins that mark each cluster, and runs a
plex-effect-free healthy-versus-disease comparison.

Each script is a Quarto document (`.qmd`): a mix of explanatory text and R code
that renders to a single self-contained HTML report.

## 2. Before you start: software and setup

You need R (4.2 or newer), RStudio, and Quarto. RStudio bundles Quarto, so
installing RStudio is usually enough. Check with `quarto --version` in a terminal.

Clone the workshop repository. In a terminal:

```bash
git clone https://github.com/YOUR_USERNAME/pedPMS-proteomics-workshop.git
```

This gives you everything you need in one folder, including the data:

```
pedPMS-proteomics-workshop/
  pedPMS-proteomics-workshop.Rproj
  pedPMS_1_Import_QC_IRS.qmd
  pedPMS_2_Exploration.qmd
  WORKSHOP_GUIDE.md
  data_raw/
    pedPMS_data_preIRS.txt
    pedPMS_data_postIRS.txt
    pedPMS_data.txt
    pedPMS_Samples.csv
    pedPMS_Protein_samples_Batch_1.csv
  results/                        (created automatically when you render)
    tables/                       (CSV outputs, created automatically)
    figures/                      (PNG outputs, created automatically)
```

Open `pedPMS-proteomics-workshop.Rproj` in RStudio. This sets the working
directory to the repo root and is what lets the scripts find their inputs and
outputs without any path editing, on Windows, Mac or Linux alike (see the note
on `here()` below).

Install the R packages the scripts use. Run this once in the R console:

```r
install.packages(c("tidyverse", "readxl", "matrixStats", "cowplot", "ggExtra",
                   "ggpubr", "ggrepel", "rstatix", "DT", "conflicted",
                   "quarto", "circlize", "scales", "here"))

install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "pcaMethods", "limma", "bluster",
                       "igraph", "ComplexHeatmap", "clusterProfiler",
                       "org.Hs.eg.db", "ReactomePA"))
```

The scripts create `results/`, `results/tables/` and `results/figures/`
themselves, so you do not need to make them. Tables are written as CSV (including
the significant proteins with their p-values and FDR), and the main figures are
saved as PNG. `results/` is not tracked by git; it is regenerated every time you
render, so nothing is lost by deleting it and starting over.

### A note on paths and `here()`

Both scripts locate the project with `base_dir <- here::here()` instead of a
hard-coded folder. The `here` package walks up from the current working
directory until it finds a marker of the project root, in this case either the
`.Rproj` file or the `.git` folder created by cloning, and returns that folder
regardless of which machine, operating system, or subfolder you are working
from. This is why nothing needs to be edited after `git clone`: `data_raw/` and
`results/` always resolve relative to the repo root.

This only works if the working directory is somewhere inside the repo when the
script runs. Opening the `.Rproj` file guarantees that. If you instead open a
`.qmd` file directly without opening the project first, set the working
directory manually before running any chunk: RStudio menu, Session > Set
Working Directory > To Project Directory.

## 3. The input files

Three quantitative files and two metadata files ship in `data_raw/` already, so
there is nothing to download or copy in separately.

`pedPMS_data.txt` family is the quantitative matrix: one row per protein group,
one column per measurement. It is tab-separated and encoded UTF-16LE with a
two-row header (the first row is a banner, the second holds the column names).
Sample columns are named `channel Rrun_plex`, for example `3 R1_598` means TMT
channel 3, injection run 1, plex 598. There are three versions of this file:

| File | Reference channel | Plexes aligned? | Use it to |
|---|---|---|---|
| `pedPMS_data_preIRS.txt` | present | no (batch effect present) | learn IRS: the script corrects it |
| `pedPMS_data_postIRS.txt` | present | yes (already corrected) | verify IRS: the script confirms it |
| `pedPMS_data.txt` | absent | yes | a plain version without the reference channel |

For the IRS part of the workshop, point script 1 at `pedPMS_data_preIRS.txt`.
The line to change is `data_txt` near the top of script 1.

`pedPMS_Samples.csv` is the TMT labelling layout: which sample went into which
channel of which plex (exported from the original "Listen" spreadsheet sheet).
It has no header row; the first three rows are reagent notes and the sample rows
begin at row 4. The columns the script reads are the labelled sample id, the TMT
tag, the plex label, the plex id and the tube id.

`pedPMS_Protein_samples_Batch_1.csv` is the sorting and preparation record
(exported from the original "LIST" spreadsheet sheet), keyed by tube id,
carrying the cell number, the lysis volume and the disease subtype for each
sample.

The scripts join these two spreadsheets on tube id to build one sample table.

## 4. The shape of the dataset

| Property | Value |
|---|---|
| Protein groups (rows) | 4,893 |
| Biological samples | 66 |
| Injections per sample | 2 (technical re-injections, R1 and R2) |
| Biological columns | 132 (66 samples x 2) |
| Reference columns | 14 (1 pooled reference per plex x 2 injections) |
| Plexes | 7 (keys 541, 551, 561, 571, 581, 591, 598) |
| Plex sizes | 10, 10, 10, 10, 10, 7, 9 |
| Subtypes | 13, including `control` |
| Healthy controls | 6, all in plex 7 |
| Missing values | about 14%, in whole plex blocks |
| Reference channel | channel 11 in every plex |
| Value scale | linear reporter intensities, strictly positive |

Two structural facts drive the analysis. First, samples were assigned to plexes in
tube order, so subtype and plex are strongly associated (Cramér's V about 0.7).
This confounding means you cannot simply "adjust for plex", because removing plex
would remove biology with it. Second, all six healthy controls sit in plex 7
alongside three disease samples, which makes plex 7 the one place a
healthy-versus-disease comparison is free of any plex effect.

## 5. Concepts you need, in plain terms

TMT (tandem mass tags) let you label up to ten (here eleven) samples with
different chemical tags, mix them, and measure them together in one mass
spectrometry run. Each tag reports its sample's abundance as a "reporter
intensity". A set of samples measured together is a plex.

Channel is the position of a tag within a plex (1 to 10, plus the reference at
11). Plex is one mixed set of tagged samples. Because only so many samples fit in
one plex, a study with 66 samples needs several plexes, and differences between
plexes (a batch effect) have to be corrected.

Injection (R1, R2) means the same plex was run through the instrument twice. The
two runs are technical replicates; the scripts confirm they agree, then average
them.

Reference channel is a pooled sample (a mix of a bit of everything) placed in one
channel of every plex. Because it is the same material in every plex, any
difference in its measured value between plexes must be batch effect, which makes
it the tool for correcting plexes.

Sample loading (SL) normalisation makes the total signal per sample equal, correc-
ting for slightly different amounts of material loaded. Internal reference scaling
(IRS) makes the plexes comparable to each other, using the reference channel. The
usual order is SL then IRS. In this dataset SL is already done (the script checks
it), and IRS is the step the workshop performs.

log2 scale is used for most calculations. Raw intensities span a huge range and
multiply; taking log2 compresses them and turns fold-changes into simple
differences (a difference of 1 on log2 means two-fold).

Missing values in TMT are usually structural, not random: a protein is measured in
a whole plex or in none of it. That is why the scripts never impute (fill in)
missing values; inventing them would mean fabricating a whole plex of data.

Confounding means two variables move together so their effects cannot be
separated. Here subtype moves with plex. The scripts quantify it and work around
it rather than pretending it is not there.

PCA (principal component analysis) summarises thousands of proteins into a few
axes that capture the main variation, so samples can be plotted and inspected.

Clustering (community detection) groups samples that have similar protein
profiles. Signature proteins are the ones that distinguish one group from the
rest. ORA (over-representation analysis) asks whether those proteins fall into
known biological pathways more often than expected.

## 6. What script 1 does, step by step

It loads packages and settings, including a block that resolves function-name
clashes (many Bioconductor packages define their own `select`, `filter` and so on,
which collide with the tidyverse; the script forces the tidyverse versions).

It reads the two spreadsheets and joins them into one sample table with 66 rows,
checking that the tube ids match exactly and that every channel is resolved. This
is where beginners most often hit errors, so the script fails loudly with a clear
message if anything does not line up.

It reads the intensity matrix and splits it into three parts: the 132 biological
columns, the 14 reference columns, and the per-protein annotation. Separating the
reference is essential, because it is a technical channel with no patient behind
it and must not be treated as a sample.

It verifies sample loading by checking that column sums are near-equal, then
performs internal reference scaling. For each protein it computes the level its
reference should have (the average across plexes), works out the per-plex scaling
that brings the reference to that level, and applies the same scaling to the real
samples. A before-and-after plot shows the plexes coming into alignment. On the
pre-IRS file this correction is visible; on the post-IRS file it is nearly
nothing, which is how the report demonstrates that IRS was already done.

It puts the data on log2 scale, averages the two injections of each sample,
measures per-protein missingness and applies a 50% cutoff, confirms the two
injections agree, and shows that missingness is blockwise (hence no imputation).

It runs PCA and quantifies the subtype-plex confounding, and lists the plex-7
samples that enable the later plex-effect-free comparison.

Finally it bundles the cleaned matrix, the sample table and the annotation into a
single object (a SummarizedExperiment) and saves it to `results/` with today's
date in the file name.

## 7. What script 2 does, step by step

It finds and loads the object script 1 saved (by matching the file name pattern,
so the date never has to be typed in), and builds a feature matrix from the 1,000
most variable complete proteins.

It reduces those to ten PCA axes, connects each sample to its ten nearest
neighbours to form a graph, and detects communities in that graph. It then checks
the clustering is stable by repeating it on random subsets of the proteins.

For each community it finds signature proteins with a limma model (the standard
tool for group comparisons, which borrows information across proteins to stay
stable when groups are small), draws a heatmap of the strongest ones with subtype
and plex annotated, and tests the up-proteins for GO pathway enrichment against
the correct background (the tested proteins, not the whole genome).

It then runs the plex-7 comparison described below.

It closes with a summary that states plainly what is exploratory and what is not.

## 8. The new functionality added for this workshop

Two capabilities were added on top of the original pipeline.

Reference-channel IRS. The original real study had no reference channel, so its
IRS could only ever be checked, not demonstrated. This workshop dataset includes a
pooled reference channel in every plex, and script 1 computes IRS from it
directly. The pre-IRS and post-IRS files let you see the same code correct a batch
effect and, on the already-corrected file, confirm that no further correction is
needed.

Plex-7 healthy-versus-disease analysis. Because all six healthy controls share
plex 7 with three disease samples, script 2 runs a comparison confined to plex 7.
Every sample in it belongs to the same plex, so the result carries no plex effect,
unlike a whole-cohort comparison where subtype and plex are tangled. It is a small
comparison (nine samples), so it is framed as a sensitivity check: its top proteins
are candidates for an independent experiment, not a finished disease signature.
This analysis is the reason the three extra controls were placed in plex 7.

## 9. How to run it

Render script 1 first, because script 2 reads what it saves. In the R console:

```r
quarto::quarto_render("pedPMS_1_Import_QC_IRS.qmd")
quarto::quarto_render("pedPMS_2_Exploration.qmd")
```

Or open each file in RStudio and click Render. Each produces a single HTML file
next to the script. If you re-run script 1 (for example after changing the input
file), re-run script 2 afterwards so it picks up the new object.

To switch between learning IRS and verifying it, change one line in script 1:

```r
data_txt = "data_raw/pedPMS_data_preIRS.txt"    # learn IRS (correction visible)
data_txt = "data_raw/pedPMS_data_postIRS.txt"   # verify IRS (already aligned)
data_txt = "data_raw/pedPMS_data.txt"           # plain file, no reference channel
```

Script 1 adapts to the input. With a reference channel (the pre/post files) it
computes IRS. With the plain file, which has no reference channel, it cannot
recompute IRS and instead verifies plex alignment diagnostically, leaving the
matrix unchanged. You do not need to edit anything for this; the script detects
the reference channel automatically.

What the scripts write to disk, besides the two HTML reports:

Into `results/`: the SummarizedExperiment object (script 1).

Into `results/tables/` (CSV): per-sample replicate agreement (script 1); all
signature proteins and the significant subset with p-values and FDR, GO enrichment
per cluster, and the full and significant plex-7 healthy-versus-disease results
with p-values and FDR (script 2).

Into `results/figures/` (PNG): sample-loading and IRS-alignment and PCA plots
(script 1); the signature heatmap and the plex-7 volcano (script 2). Figures are
numbered so they sort in reading order.

## 10. Common problems and fixes

An error like `unable to find an inherited method for function 'select'` means a
Bioconductor package's function is masking the tidyverse one. The scripts prevent
this with the `conflicted` block in the setup chunk. If you add code and hit it,
restart R (Session, Restart R) so the conflict rules take effect, then run from the
top.

Garbled column names on import usually means the file encoding was not honoured.
The read step specifies UTF-16LE for a reason; keep it.

A "file not found" at the import step almost always means `here::here()` is not
resolving to the repo root, usually because the `.qmd` was opened without opening
the `.Rproj` first (see the note on `here()` in section 2). A tube-id mismatch at
the metadata step means the sample count changed. The script expects 66 samples;
if you change that number, update the row range and the count check in the
metadata chunk (they are commented in the code).

An error at the start of script 2 about no RDS found means script 1 has not been
run yet, or did not finish. Run script 1 to completion first and confirm a file
appeared in `results/`.

## 11. How this dataset was prepared

For transparency, since the workshop teaches good data practice, here is what was
done to make the data safe to share.

The dataset was derived from a real TMT proteomics experiment, then fully
anonymised. The disease, the subtype names, the sample ids and every protein and
gene name were replaced with fictional ones. Every intensity had noise added and
was rescaled, so no original measurement survives unchanged, while the statistical
structure that makes the data a good teaching case was preserved on purpose:
the plex-subtype confounding, the blockwise missingness, the technical
re-injections, and the soft, continuous clustering.

A pooled reference channel was then added to every plex so that reference-channel
IRS could be taught rather than only verified, and a per-plex batch effect was
injected to create the pre-IRS file. Because that batch effect is centered across
plexes, IRS recovers the aligned signal almost exactly, which is the clean
before-and-after used in the lesson.

Finally, three additional healthy controls were placed in plex 7, bringing the
total to six, to give the plex-effect-free healthy-versus-disease comparison
enough samples to be worth running. These three were modelled specifically on the
two monocyte controls already in the data, so they are monocyte-like: they
resemble the existing monocyte samples more closely than the megakaryocyte
control. The control group is therefore five monocyte-like samples and one
megakaryocyte.

None of this changes how you analyse the data. It is described here only so the
workshop is honest about where the numbers come from.
