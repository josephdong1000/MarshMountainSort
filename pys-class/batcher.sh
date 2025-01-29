#!/bin/sh
module load Python/3.10.8-GCCcore-12.2.0.lua
cd /mnt/isilon/marsh_single_unit/MarshMountainSort
source .venv/bin/activate
cd /mnt/isilon/marsh_single_unit/MarshMountainSort/pys-class
python -u main.py --basefolder /mnt/isilon/marsh_single_unit/MarshMountainSort/ --mode depth --id 1155
echo "sbatcher.py subprocess finished"
