# ONT-SV-CNV

Snakemake pipelines for **structural variant (SV)** and **copy number variation (CNV)** analysis from Oxford Nanopore long-read sequencing data (hg19).

Developed for multimodal sarcoma diagnostics combining nanopore methylation profiling with genomic structural analysis.

## Pipeline overview

**SV pipeline** — Delly (long-read mode) → filtering (BND/TRA/INV/DUP) → GENCODE v19 annotation → circos plot

**CNV pipeline** — mosdepth (coverage) → Spectre (CNV calling) → GENCODE v19 annotation → genome-wide & per-chromosome plots

## Repository structure

```
ONT-SV-CNV/
├── environment_linux.yml          # Conda environment (Linux)
├── environment_macos.yml          # Conda environment (macOS)
├── data/
│   ├── example_sample/            # Example data
│   │   ├── example_sample.bam
│   │   └── example_sample.bam.bai
│   └── reference/                 # User-provided reference files
│       ├── hg19.fasta
│       ├── hg19.fasta.fai
│       ├── gencode.v19.annotation.gtf
│       ├── hg19_blacklist.bed     # Optional
│       └── hg19.mdr               # Optional (auto-generated)
├── workflow/
│   ├── Analysis_ONT_SV            # SV Snakefile
│   ├── Analysis_ONT_CNV           # CNV Snakefile
│   └── scripts/
│       ├── annotate_BND.R
│       ├── annotate_CNV.R
│       ├── plot_circos.R
│       ├── plot_CNV.R
│       └── plot_CNV_chromosome.R
└── Output/                        # Generated automatically
```

---

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/AntoineQB/ONT-SV-CNV.git
cd ONT-SV-CNV
```

### 2. Create the conda environment

**Linux:**
```bash
conda env create -f environment_linux.yml
conda activate ont_sv_cnv
```

**macOS (Intel):**
```bash
conda env create -f environment_macos.yml
conda activate ont_sv_cnv
```

**macOS (Apple Silicon M1/M2/M3/M4):**
```bash
CONDA_SUBDIR=osx-64 conda env create -f environment_macos.yml
conda activate ont_sv_cnv
conda config --env --set subdir osx-64
```

### 3. Download reference files

```bash
mkdir -p data/reference

# hg19 reference genome
wget -O data/reference/hg19.fasta.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip data/reference/hg19.fasta.gz
samtools faidx data/reference/hg19.fasta

# GENCODE v19 annotation
wget -O data/reference/gencode.v19.annotation.gtf.gz \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip data/reference/gencode.v19.annotation.gtf.gz
```

### 4. Prepare your data

Place your BAM file (aligned to hg19, sorted, indexed) as the example file:

```
data/example/example.bam
data/example/example.bam.bai
```

### 5. Run the pipelines

```bash
conda activate ont_sv_cnv

# SV analysis
snakemake -s workflow/Analysis_ONT_SV --cores 8 \
    --config sample="YOUR_SAMPLE"

# CNV analysis
snakemake -s workflow/Analysis_ONT_CNV --cores 8 \
    --config sample="YOUR_SAMPLE"
```

---

## Parameters

### SV pipeline

| Parameter | Description | Default |
|---|---|---|
| `sample` | **Required.** Sample name (must match folder name in `data/`) | — |
| `highlight_chr1` | First chromosome to highlight in circos | `NA` |
| `highlight_chr2` | Second chromosome to highlight in circos | `NA` |
| `sample_label` | Display label for the circos plot | sample name |

```bash
# Example with highlighting
snakemake -s workflow/Analysis_ONT_SV --cores 8 \
    --config sample="MY_SAMPLE" highlight_chr1="18" highlight_chr2="X" sample_label="Patient_01"
```

### CNV pipeline

| Parameter | Description | Default |
|---|---|---|
| `sample` | **Required.** Sample name (must match folder name in `data/`) | — |
| `highlight_chr` | Chromosome to highlight on genome-wide plot | `NA` |
| `bin_size` | Mosdepth coverage bin size (bp) | `1000` |
| `min_cnv_len` | Minimum CNV length for Spectre (bp) | `100000` |
| `cancer` | Enable Spectre cancer mode | `false` |

```bash
# Example with cancer mode
snakemake -s workflow/Analysis_ONT_CNV --cores 8 \
    --config sample="MY_SAMPLE" highlight_chr="12" cancer=true
```

---

## Output

### SV results (`Output/SV_Results/{sample}/`)

| File | Description |
|---|---|
| `{sample}_delly.bcf` | Raw Delly SV calls |
| `{sample}_delly.vcf` | Converted VCF |
| `{sample}_delly_filtered.vcf` | Filtered VCF (BND, TRA, INV, DUP only) |
| `{sample}_BND_annotated.csv` | Annotated breakpoints with gene names |
| `figures/{sample}_circos.pdf/png` | Circos plot |

### CNV results (`Output/CNV_Results/{sample}/`)

| File | Description |
|---|---|
| `{sample}.regions.bed.gz` | Mosdepth coverage bins |
| `{sample}.vcf.gz` | Spectre CNV calls |
| `{sample}_CNV_annotated.csv` | Annotated CNVs with overlapping genes |
| `figures/{sample}_CNV_genomewide.pdf/png` | Genome-wide log2 ratio plot |
| `figures/{sample}_CNV_{chr}.pdf/png` | Per-chromosome plots (chr1–chrX) |

---

## Troubleshooting

- **`command not found`** — Make sure the environment is activated: `conda activate ont_sv_cnv`
- **Apple Silicon install fails** — Use the `CONDA_SUBDIR=osx-64` prefix (see step 2)
- **BAM not found** — Check that your file is at `data/SAMPLE/SAMPLE.bam` with its `.bam.bai` index
- **Re-run after failure** — Snakemake resumes automatically; use `--forceall` for a clean re-run
- **Dry run** — Preview without executing: add `-n -p` to any snakemake command

---

## Citation

If you use this pipeline, please cite:

- **Delly**: Rausch T, et al. *Bioinformatics.* 2012;28(18):i333–i339.
- **Spectre**: Smolka M, et al. [github.com/fritzsedlazeck/Spectre](https://github.com/fritzsedlazeck/Spectre)
- **mosdepth**: Pedersen BS, Quinlan AR. *Bioinformatics.* 2018;34(5):867–868.

## License

This project is licensed under the [MIT License](LICENSE).
