# MarshMountainSort

A Python pipeline for processing and analyzing neural recordings using MountainSort5.

## Overview

This repository contains scripts for:
- Running MountainSort5 spike sorting on neural recordings
- Processing data from Intan and DataWave recording systems
- Batch processing multiple recordings via SLURM job submission
- Generating plots and other visualizations of sorted units

## Usage

The main processing pipeline can be run through `main.py` with various command line arguments:

```bash
python main.py --id 1234 --datadir rhds --intanport A -aI A B
```

To see all available command line arguments and their descriptions, run `python main.py --help`.

See `sbatcher.py` for a batch processing script that runs `main.py` on multiple recordings.

Currently all figures are manually generated with `marshmakestatsandfigs-class.ipynb`.

