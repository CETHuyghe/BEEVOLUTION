"""

### Main snakefile to analyse the assembled bumblebee MAGs ###

# Date: 07/09/2026
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

TMP_DIR = config["tmp_dir"]
SNK_DIR = config["snakefile_dir"]

OUT_DIR = os.path.join(SNK_DIR,"output/")
ENV_DIR = os.path.join(SNK_DIR,"envs/")
PKG_DIR = os.path.join(SNK_DIR,"packages/")
LOG_DIR = os.path.join(SNK_DIR,"logs/")
BSC_DIR = os.path.join(SNK_DIR,"bash_scripts/")

# rule all

#rule all:
#    input:
#        os.path.join(OUT_DIR, "01_taxprofiling/04_instrain/01_scafftobin/scafftobin.stb")

rule all:
    input:
        expand(
            os.path.join(OUT_DIR, "01_taxprofiling/04_instrain/02_profile/{sample}"), sample= config["samples"]
        )


# rules to include

include: "rules/01_taxprofiling.smk"
