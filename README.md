# QuagmiR
A python-based miRNA sequencing pipeline for isomiR quantification and analysis

## Necessary input files (see Usage details)
- Unmapped RNA-seq files  with .fastq_ready extension
- File with miRNA information (name, motif, consensus)

## Installation
1. Download repository: `git clone https://github.com/kevchn/quagmir`
2. Go into local quagmir folder: `cd quagmir`
3. Install python dependencies: `pip install -r requirements.txt`

## Usage
1. Add your .fastq_ready samples into the data folder:
<pre>
├── LICENSE
├── README.md
├── Snakefile
├── data
│   ├── sample.fastq_ready
│   ├── <b>YOUR_SAMPLE_HERE.fastq_ready</b>
├── motif-consensus.fa
├── report.txt
└── requirements.txt
</pre>
2. Edit the motif-consensus.fa file to insert your miRNA information with the following format:
<pre>
\>miRNA_name miRNA_motif
miRNA_consensus_sequence

\>passenger-shRNA-mir21-ORF59-5p-1 ACACCCTGGCCGGGT
CCGACACCCTGGCCGGGTTGT
</pre>
3. Run pipeline: `snakemake`
