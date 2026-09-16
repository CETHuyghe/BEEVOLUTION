### Snakemake MAG analyses for taxonomic profiling ###

# Date: 08/09/2026
# Person: Charlotte E.T. Huyghe
# Location: UNIL Engel Lab

# This snakemake file includes the rules used for the taxonomic profiling of Bombini hindgut samples
# It takes the preprocessed reads per sample and maps it to a generated reference database including the dereplicated MAGs and closely related NCBI genomes

OUT_DIR_01 = os.path.join(OUT_DIR,"01_taxprofiling/")

wildcard_constraints:
    sample="[^_]+"

wildcard_constraints:
    mags="[^_]+"


MAGS = config["mags"]

MAG_FILES = [
    f
    for mag_files in MAGS.values()
    for f in mag_files
]

# 1. Dereplication of MAGs

# dereplication of all generated and quality-checked Bombini MAGs, including almost identical genomes from NCBI, using dRep

rule drep:
    input:
        mags = MAG_FILES
    output:
        done = os.path.join(OUT_DIR_01,"01_drep/.done")
    threads: 10
    resources:
        mem_mb=50000,
        runtime="1d",
        account="pengel_spirit"
    params:
        completeness = 50,
        contamination = 5,
        outdir = os.path.join(OUT_DIR_01,"01_drep/")
    conda:
        os.path.join(ENV_DIR,"drep_env.yaml")
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}"

        dRep dereplicate {params.outdir} \
            -g {input.mags} \
            -comp {params.completeness} \
            -con {params.contamination} \
            -p {threads}

        touch "{output.done}"
        """

# Make a list of all remaining dereplicated genomes

rule list_drepmags:
    input:
        mags = os.path.join(OUT_DIR_01,"01_drep/.done")
    output:
        txt = os.path.join(OUT_DIR_01,"01_drep/drepmag_ID.txt"),
        loc = os.path.join(OUT_DIR_01,"01_drep/drepmag_loc.txt")
    threads: 1
    resources:
        mem_mb=1000,
        runtime="1d",
        account="pengel_spirit"
    params:
        indir = os.path.join(OUT_DIR_01,"01_drep/dereplicated_genomes/")
    shell:
        r"""

        ls "{params.indir}" | cut -d '.' -f1,2 > "{output.txt}"

        ls "{params.indir}" > "{output.loc}"

        """

# 2. Generate reference database including all dereplicated genomes for mapping
# concatenate all genomes into one fasta file used as database

rule concat_drepmags:
    input:
        txt = os.path.join(OUT_DIR_01, "01_drep/drepmag_ID.txt"),
        loc = os.path.join(OUT_DIR_01, "01_drep/drepmag_loc.txt")
    output:
        mag = os.path.join(OUT_DIR_01, "02_db_drepmag/db_drep_mags.fa")
    params:
        outdir = os.path.join(OUT_DIR_01, "02_db_drepmag/"),
        fasdir = os.path.join(OUT_DIR_01, "01_drep/dereplicated_genomes")
    threads: 2
    resources:
        mem_mb=1000,
        runtime="1d",
        account="pengel_spirit"
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}"

        > "{output.mag}"

        while read -r mag
        do
            fas=$(grep "^${{mag}}" "{input.loc}")

            fasta="{params.fasdir}/${{fas}}"

            cat "${{fasta}}" | sed "s/^>/>${{mag}}-/" >> "{output.mag}"

        done < "{input.txt}"
        """

# index the concatenated database for bowtie mapping

rule index_bowtiedb:
    input:
        mag = os.path.join(OUT_DIR_01,"02_db_drepmag/db_drep_mags.fa")
    output:
        db = os.path.join(OUT_DIR_01,"02_db_drepmag/db_drep_mags.1.bt2")
    threads: 8
    resources:
        mem_mb=10000,
        runtime="1d",
        account="pengel_spirit"
    params:
        prefix = os.path.join(OUT_DIR_01,"02_db_drepmag/db_drep_mags")
    conda:
        os.path.join(ENV_DIR,"bowtie2_env.yaml")
    shell:
        r"""
        set -euo pipefail

        bowtie2-build --threads "{threads}" "{input.mag}" "{params.prefix}"
  
        """

# 2. Map sample reads to the database including all dereplicated genomes
# use bowtie to map

rule map_drepmags:
    input:
        db = os.path.join(OUT_DIR_01,"02_db_drepmag/db_drep_mags.1.bt2"),
        fastq1 = os.path.join(config["nohost_fastq_dir"],"{sample}_nohost_R1.fastq.gz"),
        fastq2 = os.path.join(config["nohost_fastq_dir"],"{sample}_nohost_R2.fastq.gz")
    output:
        bam = os.path.join(OUT_DIR_01,"03_mapping/{sample}_mapped.bam")
    threads: 4
    resources:
        mem_mb=10000,
        runtime="1d",
        account="pengel_spirit"
    params:
        db = os.path.join(OUT_DIR_01,"02_db_drepmag/db_drep_mags"),
        tmp = os.path.join(TMP_DIR,"{sample}_drepmag_mapped.sam")
    conda:
        os.path.join(ENV_DIR,"bowtie2_env.yaml")
    shell:
        r"""
        set -euo pipefail

        # map fastp sequences where host was removed to MAG database

        bowtie2 --threads "{threads}" \
            -x "{params.db}" \
            -1 "{input.fastq1}" \
            -2 "{input.fastq2}" \
            | samtools sort -o "{params.tmp}" 


        # compress sam to bam

        samtools view -bS "{params.tmp}" > "{output.bam}"

        """

# 4. Construct taxonomic profile for each sample
# generate scaffold-to-bin file of the genome database needed for inStrain

rule scafftobin_stb:
    input:
        db = os.path.join(OUT_DIR_01,"02_db_drepmag/db_drep_mags.fa")
    output:
        stb = os.path.join(OUT_DIR_01,"04_instrain/01_scafftobin/scafftobin.stb")
    params:
        outdir = os.path.join(OUT_DIR_01,"04_instrain/01_scafftobin")
    shell:
        r"""

        mkdir -p "{params.outdir}"

        awk '/^>/ {{
            header = substr($0, 2)
            split(header, parts, "-")
            print header "\t" parts[1]
        }}' {input.db} > {output.stb}

        """

# use inStrain with mapping info to generate taxonomic profiles

rule instrain_profile:
    input:
        db = os.path.join(OUT_DIR_01,"02_db_drepmag/db_drep_mags.fa"),
        bam = os.path.join(OUT_DIR_01,"03_mapping/{sample}_mapped.bam"),
        stb = os.path.join(OUT_DIR_01,"04_instrain/01_scafftobin/scafftobin.stb")
    output:
        profile = os.path.join(OUT_DIR_01,"04_instrain/02_profile/{sample}")
    threads: 10
    resources:
        mem_mb=15000,
        runtime="1d",
        account="pengel_spirit"
    params:
        outdir = os.path.join(OUT_DIR_01,"04_inStrain/02_profile")
    conda:
        os.path.join(ENV_DIR,"instrain_env.yaml")
    shell:
        r"""
        set -euo pipefail

        mkdir -p "{params.outdir}"

        # map fastp sequences where host was removed to MAG database

        inStrain profile "{input.bam}" \
            "{input.db}" \
            -o "{output.profile}" \
            -p "{threads}" \
            -s "{input.stb}" 

        """
