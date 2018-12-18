Code used for the analysis of the 16S dataset of sarcoidosis patients and patients with other interstitial lung disease.

# Requirements

## Python
```bash
# create virt. env.
virtualenv -p /usr/bin/python3 venv
# activate it
source venv/bin/activate
# install requirements
pip install -r requirements.txt
```

## R
See file `requirements_R.txt`.

## Tools
- FastQC: `v0.11.7`
- KronaTools: `2.5`
- Optional: (*see note below*: "Note: precomputed output")
    - LotuS: `1.565`
    - BLAST: `2.6.0+`
        - Requires the NCBI nr/nt database
    - EDirect: `9.50`

You may need to adjust the paths to these tools in `pipeline.json`.
Currently all paths are linking to `tools/`.

## Data
Raw FASTQ files are expected to be in `data/fastq` (otherwise adjust the path in `pipeline.json`) and their naming scheme should be
`{sid}_{aid}_{rid}.fastq` where `sid` is the sample ID (same as in the column `SampleID` in `data/meta_01.tsv`),
`aid` the amplicon ID (`V1V2` or `V3V4`) and `rid` the read ID (`R1` or `R2`).

# Note: precomputed output
## LotuS output
Rerunning LotuS should produce the same OTUs but their IDs can change.
Thus, the output files are provided in `data/precomputed/lotus` and will be used by the script.
If you want to rerun LotuS comment out the currently used rule `rule run_lotus_placeholder` and remove the comment flags for `run_lotus` (see comments in `pipeline.snake`).

## BLASTn hits
Since the results depend on the used database version the output files are provided in `data/precomputed/blastn_nt`
and will be used by the script.
The original rules `run_blastn_asn`, `run_blastn_tab`, and `run_blastn_tax` were commented out and
a placeholder for `run_blastn_tax` was added (see comments in `pipeline.snake`).

# Data
- Raw reads (FASTQ files): [TODO](XXX)
- Sample meta information: files `meta_01.tsv` and `meta_02.tsv` in `data/`

# Pipeline
Run
```bash
# virt. env. must be activated
snakemake -s pipeline.snake
```

## Input data
- Samples
    - `data/meta_01.tsv`: Meta data including IDs and `Indication`
    - `data/ID_map.tsv`: Mapping of sample IDs to new IDs for the manuscript:
        - `NTC` for NTC sample
        - `S<i>` for sarcoidosis samples
        - `C<i>` for "control" samples
    - `data/meta_02.tsv`: Additional sample information, e.g. smoker status etc.
- NGS data
    - `data/fastq_raw`
- Other
    - `data/primer.txt`: Primer sequences (needed for LotuS; see "Note: precomputed output")
    - `data/blastn_nt/`: BLASTn hits (see "Note: precomputed output")

## Steps/resulting files
- Creating quality reports:
    - FASTQC on eahc FASTQ file: `data/fastqc/`
    - MultiQC: Summary over all files: `data/multiqc/`
        - `multiqc.html`: MultiQC report
        - `multiqc_plots.pdf`: Plots based on MultiQC report
            - Read counts per sample, group and amplicon
        - `multiqc_data/multiqc_general_stats.txt`: General statistics of FASTQ files (e.g. read count)
- LotuS analysis to create OTUs: `data/lotus/` (see "Note: precomputed output")
- OTU filtering: `data/filtered_otus/`
    - `<amplicon>.tsv`: Table with OTU taxonomy (LotuS, BLAST), filtering flag, and sample counts
    - `<amplicon>.pdf`: Plots
        - Percentage of counts covered by removed/kept OTUs
        - Phylogenetic tree of OTUs where removed OTUs are highlighted
- Krona plots (taxonomic composition): `data/krona/krona.html`
    - Per amplicon, over all samples and per indication group
- Analysis using `ampvis2`: `data/ampvis2`
    - `<amplicon>.pdf`:
        - Rarefaction curves
        - Alpha-diversity
    - `<amplicon>_alphadiv.tsv`: Alpha-diversity values per sample
- CLR-transformation of OTU counts: `data/clr`
- Analysis of transforomed values: `data/clr_analysis/`
    - `<amplicon>_trans_pca.pdf`: PCA plots
    - `<amplicon>_trans_pheatmap.pdf`: Heatmap plots
- Looking at specific taxa using transformed values: `data/clr_taxa/`
    - `<amplicon>.pdf`: Heatmaps using OTUs from specified taxa
        - All OTUs
        - Per taxon
- Diff. analysis using `ALDEx2`: `aldex2_analysis/`
    - `<amplicon>_aldex2.tsv`: Table with results, important columns are:
        - OTU ID (`OTU`)
        - Difference within/between groups (`diff.win`, `diff.btw`)
        - Effect size (`effect`)
        - P-value of the Wilcoxon test (`wi.ep`)
    - `<amplicon>_aldex2_ov.pdf`: Rresults summary plot
        - Within diff. vs. between diff. + boundaries where both are equal
        - CLR-value vs. between diff.
        - Between diff. vs. log(p-value) + p-value boundary
        - Effect size vs. log(p-value) + cutoff boundaries
