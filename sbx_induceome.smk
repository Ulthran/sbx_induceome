import csv


INDUCEOME_FP = Cfg["all"]["output_fp"] / "virus" / "induceome"
try:
    SBX_INDUCEOME_VERSION = get_ext_version("sbx_induceome")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_INDUCEOME_VERSION = "0.0.0"


def get_induceome_ref_var(sample: str) -> str:
    # Check if the pattern "T(0-9)+_" is in sample and if so return the number
    if re.match(r"T\d+_", sample):
        # Extract the number and return it
        return re.search(r"\d+", re.match(r"T\d+_", sample).group()).group()
    # Using this as a default value
    # Not sure if there's another default we might want
    return "1"


# Ingest reference mapping
SBX_INDUCEOME_REF_FP = Path(Cfg["sbx_induceome"]["reference_fp"])
with open(Cfg["sbx_induceome"]["mapping_fp"]) as f:
    csv_reader = csv.reader(f, delimiter=",")
    header = next(csv_reader)
    sample_id_idx = header.index("Sample_ID")
    strain_idx = header.index("strain")
    SBX_INDUCEOME_REF_MAP = {
        row[sample_id_idx]: SBX_INDUCEOME_REF_FP / f"{row[strain_idx]}.fasta"
        for row in csv_reader
    }
# Update reference mapping with true sample names
# Needed because we only get "PSP####-1" from metadata sheet but the full sample name is "PSP####-1_S#"
# PSP value SHOULD be unique though
NEW_SBX_INDUCEOME_REF_MAP = {}
for sample_id, ref in SBX_INDUCEOME_REF_MAP.items():
    if [x for x in Samples if sample_id in x]:
        NEW_SBX_INDUCEOME_REF_MAP[[x for x in Samples if sample_id in x][0]] = ref
    else:
        print(f"Sample {sample_id} not found in Samples list")
SBX_INDUCEOME_REF_MAP = NEW_SBX_INDUCEOME_REF_MAP
# Remove any references that don't exist
# Necessary in case some controls don't have references
SBX_INDUCEOME_REF_MAP = {
    sample: ref for sample, ref in SBX_INDUCEOME_REF_MAP.items() if ref.exists()
}
print(SBX_INDUCEOME_REF_MAP)
# Reduce samples to only those with references
SBX_INDUCEOME_SAMPLES = {
    sample: fps
    for sample, fps in Samples.items()
    if sample in SBX_INDUCEOME_REF_MAP.keys()
}


localrules:
    all_induceome,


rule all_induceome:
    input:
        expand(INDUCEOME_FP / "peaks" / "{sample}.png", sample=SBX_INDUCEOME_SAMPLES),
        expand(INDUCEOME_FP / "peaks" / "{sample}.csv", sample=SBX_INDUCEOME_SAMPLES),
        expand(INDUCEOME_FP / "blastx" / "{sample}.btf", sample=SBX_INDUCEOME_SAMPLES),
        expand(
            INDUCEOME_FP / "phold" / "{sample}_plot" / "phold.png",
            sample=SBX_INDUCEOME_SAMPLES,
        ),


rule induceome_bwa_index:
    """Index the reference genome for BWA"""
    input:
        list(set(SBX_INDUCEOME_REF_MAP.values())),
    output:
        [
            f"{fp}.{ext}"
            for ext in ["amb", "ann", "bwt", "pac", "sa"]
            for fp in set(SBX_INDUCEOME_REF_MAP.values())
        ],
    log:
        LOG_FP / "induceome_bwa_index.log",
    benchmark:
        BENCHMARK_FP / "induceome_bwa_index.tsv"
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    shell:
        """
        echo "Indexing reference genomes with BWA" > {log}
        for ref in {input}
        do
            bwa index $ref 2>> {log}
        done
        """


rule induceome_bwa_mem:
    """Align reads to the reference genome using BWA"""
    input:
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
        ref=lambda wildcards: SBX_INDUCEOME_REF_MAP[wildcards.sample],
        indexes=lambda wildcards: expand(
            str(SBX_INDUCEOME_REF_MAP[wildcards.sample]) + ".{ext}",
            ext=["amb", "ann", "bwt", "pac", "sa"],
        ),
    output:
        INDUCEOME_FP / "aligned" / "{sample}.sam",
    log:
        LOG_FP / "induceome_bwa_mem_{sample}.log",
    benchmark:
        BENCHMARK_FP / "induceome_bwa_mem_{sample}.tsv"
    threads: Cfg["sbx_induceome"]["threads"]
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    shell:
        "bwa mem -t {threads} {input.ref} {input.reads} > {output} 2> {log}"


rule induceome_samtools_sort:
    """Sort, index, and produce pileups from the aligned reads"""
    input:
        sam=INDUCEOME_FP / "aligned" / "{sample}.sam",
        ref=lambda wildcards: SBX_INDUCEOME_REF_MAP[wildcards.sample],
    output:
        sorted=temp(INDUCEOME_FP / "aligned" / "{sample}.sam.sorted.bam"),
        pileup=INDUCEOME_FP / "pileups" / "{sample}.pileup",
    log:
        LOG_FP / "induceome_samtools_sort_{sample}.log",
    benchmark:
        BENCHMARK_FP / "induceome_samtools_sort_{sample}.tsv"
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    shell:
        """
        samtools view -bS {input.sam} | samtools sort -o - > {output.sorted} 2> {log}
        samtools index {output.sorted} 2>> {log}
        samtools mpileup -A -a -Q 0 -f {input.ref} {output.sorted} > {output.pileup} 2>> {log}
        """


rule induceome_find_peaks:
    """Find coverage peaks in the pileup"""
    input:
        pileup=INDUCEOME_FP / "pileups" / "{sample}.pileup",
        ref=lambda wildcards: SBX_INDUCEOME_REF_MAP[wildcards.sample],
    output:
        peaks_img=INDUCEOME_FP / "peaks" / "{sample}.png",
        peaks_csv=INDUCEOME_FP / "peaks" / "{sample}.csv",
        peaks_contigs=INDUCEOME_FP / "peaks" / "{sample}_contigs.fasta",
    log:
        LOG_FP / "induceome_find_peaks_{sample}.log",
    benchmark:
        BENCHMARK_FP / "induceome_find_peaks_{sample}.tsv"
    params:
        min_width=Cfg["sbx_induceome"]["min_width"],
        smoothing_factor=Cfg["sbx_induceome"]["smoothing_factor"],
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    script:
        "scripts/induceome_find_peaks.py"


###
# BLAST
###


rule induceome_blastx:
    """Run blastx on untranslated genes against a target db and write to blast tabular format."""
    input:
        peaks_contigs=INDUCEOME_FP / "peaks" / "{sample}_contigs.fasta",
    output:
        INDUCEOME_FP / "blastx" / "{sample}.btf",
    benchmark:
        BENCHMARK_FP / "run_induceome_blastx_{sample}.tsv"
    log:
        LOG_FP / "run_induceome_blastx_{sample}.log",
    params:
        blast_db=Cfg["sbx_induceome"]["blast_db"],
    threads: Cfg["sbx_induceome"]["threads"]
    resources:
        mem_mb=24000,
        runtime=720,
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    shell:
        """
        if [ -s {input} ]; then
            export BLASTDB=$(dirname {params.blast_db})
            blastx \
            -query {input} \
            -db $(basename {params.blast_db}) \
            -outfmt '7 "qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"' \
            -num_threads {threads} \
            -evalue 0.05 \
            -max_target_seqs 100 \
            -out {output} \
            2>&1 | tee {log}
        else
            echo "Caught empty query" >> {log}
            touch {output}
        fi
        """


###
# PHAROKKA/PHOLD/PHYNTENY
###


rule induceome_install_ph_databases:
    output:
        pharokka=Path(Cfg["sbx_induceome"]["pharokka_db"])
        / "pharokka_db"
        / ".installed",
        annotations=Path(Cfg["sbx_induceome"]["phold_db"]) / "phold_annots.tsv",
    benchmark:
        BENCHMARK_FP / "induceome_install_ph_databases.tsv"
    log:
        LOG_FP / "induceome_install_ph_databases.log",
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    shell:
        """
        echo "Installing pharokka database" > {log}
        install_databases.py -o $(dirname {output.pharokka}) >> {log} 2>&1
        touch {output.pharokka}

        echo "Installing phold annotations" >> {log}
        phold install --database $(dirname {output.annotations}) >> {log} 2>&1
        """


rule induceome_pharokka:
    input:
        contigs=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
        pharokka_db=Path(Cfg["sbx_induceome"]["pharokka_db"]) / "pharokka_db",
    output:
        pharokka_gbk=INDUCEOME_FP / "pharokka" / "{sample}_pharokka.gbk",
    benchmark:
        BENCHMARK_FP / "induceome_pharokka_{sample}.tsv"
    log:
        LOG_FP / "induceome_pharokka_{sample}.log",
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    shell:
        """
        if [ ! -s {input.contigs} ]; then
            touch {output.pharokka_gbk}
        else
            pharokka run -i {input.contigs} -o $(dirname {output.pharokka_gbk}) --database {input.pharokka_db} --force
        fi
        """


rule induceome_phold_predict:
    input:
        contigs=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
        annotations=Path(Cfg["sbx_induceome"]["phold_db"]) / "phold_annots.tsv",
    output:
        _3di=INDUCEOME_FP / "phold" / "{sample}_predict" / "phold_3di.fasta",
    benchmark:
        BENCHMARK_FP / "induceome_phold_predict_{sample}.tsv"
    log:
        LOG_FP / "induceome_phold_predict_{sample}.log",
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    shell:
        """
        if [ ! -s {input.contigs} ]; then
            touch {output._3di}
        else
            phold predict -i {input.contigs} -o $(dirname {output._3di}) --database $(dirname {input.annotations}) --cpu --force
        fi
        """


rule induceome_phold_compare:
    input:
        contigs=ASSEMBLY_FP / "megahit" / "{sample}_asm" / "final.contigs.fa",
        _3di=INDUCEOME_FP / "phold" / "{sample}_predict" / "phold_3di.fasta",
        annotations=Path(Cfg["sbx_induceome"]["phold_db"]) / "phold_annots.tsv",
    output:
        gbk=INDUCEOME_FP / "phold" / "{sample}_compare" / "phold.gbk",
    benchmark:
        BENCHMARK_FP / "induceome_phold_compare_{sample}.tsv"
    log:
        LOG_FP / "induceome_phold_compare_{sample}.log",
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    threads: 8
    shell:
        """
        if [ ! -s {input._3di} ]; then
            touch {output.gbk}
        else
            phold compare -i {input.contigs} --predictions_dir $(dirname {input._3di}) -o $(dirname {output.gbk}) --database $(dirname {input.annotations}) -t 8 --force
        fi
        """


rule induceome_phold_plot:
    input:
        gbk=INDUCEOME_FP / "phold" / "{sample}_compare" / "phold.gbk",
    output:
        png=INDUCEOME_FP / "phold" / "{sample}_plot" / "phold.png",
    benchmark:
        BENCHMARK_FP / "induceome_phold_plot_{sample}.tsv"
    log:
        LOG_FP / "induceome_phold_plot_{sample}.log",
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    shell:
        """
        phold plot -i {input.gbk} -o $(dirname {output.png}) --force
        """
