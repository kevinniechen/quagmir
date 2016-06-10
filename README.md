# QuagmiR
A python-based miRNA sequencing pipeline for isomiR quantification and analysis

## Installation
1. Download repository: `git clone https://github.com/kevchn/quagmir`
2. Go into local quagmir folder: `cd quagmir`
3. Install python dependencies: `pip install -r requirements.txt`

## Usage

1. Add your .fastq_ready samples into the data folder:
  ```
  ├── LICENSE
  ├── README.md
  ├── Snakefile
  ├── data
  │   ├── sample.fastq_ready
  │   ├── <b>YOUR_SAMPLE_HERE.fastq_ready</b>
  ├── motif-consensus.fa
  ├── report.txt
  └── requirements.txt
  ```

2. Edit the motif-consensus.fa file to insert your miRNA information with the following format:
  ```
  >miRNA_name miRNA_motif
  miRNA_consensus_sequence

  >passenger-shRNA-mir21-ORF59-5p-1 ACACCCTGGCCGGGT
  CCGACACCCTGGCCGGGTTGT
  ```

3. Run pipeline: `snakemake`

## Notes
The step of collapsing sample files takes the longest time, but once the samples are collapsed, and you need to re-run the pipeline, the pipeline will automatically start from the collapsed files and take a far smaller amount of time

Output will be a report-**sample_name**.txt file for each sample
