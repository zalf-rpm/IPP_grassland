#!/bin/bash -x
conda activate poetry_py310
sbatch -c2 --mem=1G --array=1-1000%40 poetry run python run.py end_year=2098 rcp=85