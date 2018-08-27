![alt tag](https://img.shields.io/github/issues/kevchn/quagmir.svg)
![alt tag](https://img.shields.io/wercker/ci/wercker/docs.svg)
![alt tag](https://img.shields.io/dub/l/vibe-d.svg)

# QuagmiR
A python-based miRNA sequencing pipeline for isomiR quantification and analysis.

**This repository is no longer supported.**

**Please see https://github.com/Gu-Lab-RBL-NCI/QuagmiR for the newest version of this tool.**

![alt tag](http://g.recordit.co/GbfgMIq28L.gif)

## Dependencies
* Make sure that you have [Python 3.4+](https://www.python.org/downloads/) installed (type `python --version` in the console)
* Make sure you have the latest version of pip: `pip3 install -U pip`
* Make sure you have [Miniconda](http://conda.pydata.org/docs/install/quick.html) installed

## Installation
1. Download repository: `git clone https://github.com/kevchn/quagmir`
2. Go into local quagmir folder: `cd quagmir`
3. Install Python dependencies: `conda env create -f environment.yml`

## Quickstart

1. Add your **.fastq** samples into the **data** folder (a sample has been provided for testing):
  <pre>
  ├── LICENSE
  ├── README.md
  ├── config.yaml
  ├── Snakefile
  ├── environment.yml
  ├── data/
  │   ├── sample.fastq_ready
  │   └── <b>YOUR_FILE_HERE.fastq_ready</b>
  |   ├── collapsed/
  ├── motif-consensus.fa
  └── results/
  │   └── tabular/
  </pre>

2. Edit the **motif-consensus.fa** file to insert your miRNA information with the following format:
  ```
  >miRNA_name miRNA_motif
  miRNA_consensus_sequence

  >passenger-shRNA-mir21-ORF59-5p-1 ACACCCTGGCCGGGT
  CCGACACCCTGGCCGGGTTGT
  ```

3. Edit the **config.yaml** file to change configuration options if needed (default values fine in most use cases):
  ```
  # DISPLAY
  min_ratio: .001
  min_read: 9

  display_summary: True
  display_sequence_info: True
  display_nucleotide_dist: True

  # FUNCTION
  destructive_motif_pull: False

  # INPUT
  motif_consensus_file: 'motif-consensus.fa'
  ```

4. Run pipeline: `bash run.sh` or activate conda env and run `snakemake`

## Update
Run the commands ```git reset --hard``` and ```git pull```.

## Additional information
Quagmir skips any reads that it considers not to be an isomir of the miRNA, even
if the reads are pulled in by the motif of the miRNA. These skipped reads are
reported in the run log in the logs/ folder.

Currently, Quagmir by default skips reads that:

(1) Have none-templating trimming and tailing combinations on the 5' end

**Biologically, 5' modifications are very rare/non-existent. This removes reads
that supposedly have these modifications. Keeps isomirs that are a result of
imprecise cleavage. False positive: Also keeps potentially 'fake' isomirs that
have non-templating modifications beyond the extent of the provided consensus
mirna.**

(2) Have a mutation 3nt from the consensus end (shares rest of end)

**Most likely to be sequencing errors rather than trimming/tailing if 3nt or
more are tailed that match the canonical sequence (1/4^3 = 1.5%).

Quagmir also has an optional flag (off by default) to skip reads that:

(3) Were reported to be an isomir of a previously reported miRNA

**Biologically, a read can only belong to one miRNA**

## Notes
The step of collapsing sample files takes the longest time, but once the samples are collapsed, and you need to re-run the pipeline, the pipeline will automatically start from the collapsed files and take a far shorter amount of time.

To save time you may consider to skip distance metric (weighted-levenstein distance) which is calculated on group results.

Output will be a **sample.fastq.results.txt** file for each sample in the **results/** folder and grouped results if multiple files are provided. All resulting files are in tabular (TSV) format.

## Planned features
- Have a provided primary miRNA transcript for each miRNA either in
motif-consensus as well as gen-uniq-subseq, or have it pre-provided in local
files. Then implement features to smartly filter out reads that are not
isomiRs of a miRNA because they have statistically improbable behavior (e.g
5p non-templating addition) or sequencing errors.
