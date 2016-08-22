![alt tag](https://img.shields.io/github/issues/kevchn/quagmir.svg)
![alt tag](https://img.shields.io/wercker/ci/wercker/docs.svg)
![alt tag](https://img.shields.io/dub/l/vibe-d.svg)

# QuagmiR
A python-based miRNA sequencing pipeline for isomiR quantification and analysis

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

1. Activate your conda environment: `source activate quagmir`
2. Add your .fastq_ready samples into the data folder (a sample has been provided for testing):
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

3. Edit the **motif-consensus.fa** file to insert your miRNA information with the following format:
  ```
  >miRNA_name miRNA_motif
  miRNA_consensus_sequence

  >passenger-shRNA-mir21-ORF59-5p-1 ACACCCTGGCCGGGT
  CCGACACCCTGGCCGGGTTGT
  ```

4. Edit the **config.yaml** file to change configuration options if needed (default values fine in most use cases):
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

5. Run pipeline: `snakemake`

## Update
Run the commands ```git reset --hard``` and ```git pull```.

## Notes
The step of collapsing sample files takes the longest time, but once the samples are collapsed, and you need to re-run the pipeline, the pipeline will automatically start from the collapsed files and take a far shorter amount of time

Output will be a **sample.fastq_ready.results.txt** file for each sample in the **results/** folder
