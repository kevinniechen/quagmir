#!/bin/bash

source activate quagmir
snakemake
find ./results/ -size  0 -print0 |xargs -0 rm
find ./group_results/ -size  0 -print0 |xargs -0 rm
