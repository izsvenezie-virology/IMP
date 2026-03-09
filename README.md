# imp: Illumina Multipurpose Pipeline

`imp` is a Nextflow pipeline that processes paired-end Illumina reads.
The main workflow is in [`src/main.nf`](./src/main.nf).

## Feature keys

1. Adapters trimming (`cutadapt`).
1. Amplification primers trimming (`cutadapt`).
1. Automatic best reference search (`blast`).
1. Coverage plots (`tacos`).
1. Degenerated and non-degenerated consensus (`consenser`).

## Requirements

- Nextflow
- Conda/micromamba environment with tools in [`environment.yml`](./environment.yml)
- Input paired-end FASTQ files
- Sample metadata TSV
- Reference database FASTA (only required if any sample needs BLAST-based reference discovery)

## Installation / environment

Create an environment from the provided file:

```bash
mamba create -n imp -f environment.yml
mamba activate imp
```

## Running the pipeline

```bash
nextflow run src/main.nf --raw_reads_folder "/path/to/fastq_dir" --samples_metadata "/path/to/samples_metadata.tsv" -resume
```

Or with a config file:

```bash
nextflow run src/main.nf -c <config_file> -resume
```

## Parameters

### Required parameters

- `raw_reads_folder`: path to raw reads folder.
- `samples_metadata`: path to `TSV` with samples metadata file.

### Optional parameters

- `references_database`: path to `Fasta` for reference search, required for automatic best reference search.
- `adapters`: path to `Fasta` with adapter trimming sequences, by default trims Illumina Nextera adapters.
- `raw_reads_pattern`: pattern for raw reads files detection, by default is `"*_R{1,2}*fastq.gz"`
- `run_name`: default group name used when the `Group` column is absent from the metadata (default `"run"`).
- `virus`: virus abbreviation for extra analyses (default `"Others"`).
  - `AIV`: FluMut and Genin2; subtype search (only for sample using BLAST reference discovery)
  - `SIV`: subtype search (only for samples using BLAST reference discovery)

### Passing parameters

Parameters can be passed via command line prefixing the parameter with `--` (eg. `--raw_reads_folder "/path/to/fastq_dir"`)
or can be passed in a config file formatting the parameters as `params.PARAMETER = VALUE`:

```
params.raw_reads_folder = "/path/to/fastq_dir"
params.samples_metadata = "/path/to/samples_metadata.tsv"
```

## Samples metadata

Tab-separated file with header.

### Required columns

- `Sample`: sample ID matching FASTQ-derived sample names (filename up to `_S##`)
- `Name`: FASTA header text used in consensus sequences (`CHROM` will be replaced with the reference header)

Optional columns:

- `Group`: group label for consensus concatenation. If absent, falls back to `params.run_name`.
- `Reference`: path to a sample-specific reference FASTA. If provided and the file exists, BLAST reference discovery is skipped for that sample.
- `Subset`: subsampling fraction of reads used by BLAST reference discovery (default `1` when `Reference` is not provided).
- `Primers`: path to a primers TSV. If absent, primer trimming is skipped and only adapter trimming is performed.
- `MinimumCoverage`: minimum depth for coverage stats, plots, and consensus masking (default `20`).
- `NoConsensusFilter`: if present and non-empty, disables lofreq default filters for that sample.

Lines starting with `#` are ignored.

Example:

```tsv
Sample	Name	Group	Reference	Subset	Primers	MinimumCoverage	NoConsensusFilter
S01	A/chicken/IT/01/2026	H5
S02	A/chicken/IT/02/2026	H5	/path/to/S02_ref.fa		/path/to/primers.tsv	30	1
```

### Primers TSV

One primer per line, tab-separated:

```
<primer_id>\t<sequence>
```

Example:

```tsv
CommonA-Uni12G	GCCGGAGCTCTGCAGATATCAGCGAAAGCAGG
CommonA-Uni12	GCCAGAGCTCTGCAGATATCAGCAAAAGCAGG
```

### Reference database FASTA

Used by BLAST when `Reference` is not provided for a sample.

- `Fasta` format
- Sequence on 1 line
- Header format: `>ID|TYPE|SUBTYPE|SEGMENT`
  - `ID`: sequence unique ID, required
  - `TYPE`: Virus type (eg. influenza A, B or C), optional
  - `SUBTYPE`: Virus subtype (eg. HA-5 or NA-1), optional
  - `SEGMENT`: Segment name in case of segmented viruses (eg. PB2), required

## Outputs

Main outputs (written relative to `outputDir`, default `.`):

| Path                             | Content                                                          |
| -------------------------------- | ---------------------------------------------------------------- |
| `alignments/`                    | BWA-MEM alignment                                                |
| `alignments/references/`         | Sanitized reference used for alignment                           |
| `alignments/references/headers/` | Reference sequence headers                                       |
| `cleaned_reads/`                 | Adapter/primer-trimmed reads (symlinked)                         |
| `consensus/`                     | Degenerated consensus                                            |
| `consensus/unambiguous/`         | Non-degenerated consensus                                        |
| `coverage/`                      | Coverage plot                                                    |
| `coverage/raw/`                  | Per-base coverage table                                          |
| `reads_quality/raw/`             | FastQC reports on raw reads                                      |
| `reads_quality/clean/`           | FastQC reports on cleaned reads                                  |
| `vcfs/`                          | lofreq variant calls                                             |
| `results/<group>_consensus.fa`   | Grouped consensus sequences                                      |
| `results/<group>_flumut.xlsm`    | FluMut results (`AIV` virus only)                                |
| `results/<group>_genin2.tsv`     | Genin2 results (`AIV` virus only)                                |
| `results/subtypes.tsv`           | Subtypes detected for each sample (`AIV` and `SIV` viruses only) |
| `results/statistics.tsv`         | Aggregated run statistics                                        |
| `warnings`                       | Warnings triggered by the pipeline                               |
