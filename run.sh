#!/bin/bash

source activate quagmir
snakemake
#####find ./results/tabular/ -size  0 -print0 |xargs -0 rm
