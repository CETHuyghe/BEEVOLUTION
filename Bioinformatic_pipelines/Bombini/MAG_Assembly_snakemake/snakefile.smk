"""

### Main snakefile to process raw bumble bee reads into metagenomic assembled genomes (MAGs) ###

# Date: 20/04/2026
# Person: Charlotte E.T. Huyghe
# Location: UNIL Engel Lab

"""

# Import modules

import os
import fnmatch
import pandas as pd
import numpy as np
from glob import glob

# define config file

configfile: "./config.yaml"

# define directories

RAW_DIR = config["rawread_dir"]
TMP_DIR = config["tmp_dir"]
SNK_DIR = config["snakefile_dir"]

OUT_DIR = os.path.join(SNK_DIR,"output/")
ENV_DIR = os.path.join(SNK_DIR,"envs/")
PKG_DIR = os.path.join(SNK_DIR,"packages/")
LOG_DIR = os.path.join(SNK_DIR,"logs/")
BSC_DIR = os.path.join(SNK_DIR,"bash_scripts/")

# rule all

rule all:
    input:
        expand(
            os.path.join(OUT_DIR, "06_MAGsQC/03_gtdbtk/.done")
        )

# rules to include

include: "rules/01_preprocessing.smk"
include: "rules/02_assembly.smk"
include: "rules/03_kmerclust.smk"
include: "rules/04_backmap.smk"
include: "rules/05_binning.smk"
include: "rules/06_MAGsQC.smk"
