#!/bin/sh

cd /mnt/isilon/marsh_single_unit/MarshMountainSort

source venv/bin/activate


python -u /mnt/isilon/marsh_single_unit/MarshMountainSort/pys/marshfindbestdepths.py

echo "Find Best Depths finished."
