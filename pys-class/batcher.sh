#!/bin/sh
module load Python/3.10.8-GCCcore-12.2.0.lua
cd /mnt/isilon/marsh_single_unit/MarshMountainSort
source .venv/bin/activate
cd /mnt/isilon/marsh_single_unit/MarshMountainSort/pys-class
python -u main.py --basefolder /mnt/isilon/marsh_single_unit/MarshMountainSort/ --id 579-580 --datadir rhds --intanport B -aI A B -o 'overnight 1-17-25'
echo "sbatcher.py subprocess finished"
