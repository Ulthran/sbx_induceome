def get_induceome_path() -> Path:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_induceome":
            return Path(fp)
    raise Error(
        "Filepath for sbx_induceome not found, are you sure it's installed under extensions/sbx_induceome?"
    )


def get_induceome_ref_var(sample: str) -> str:
    # Check if the pattern "T(0-9)+_" is in sample and if so return the number
    if re.match(r"T\d+_", sample):
        # Extract the number and return it
        return re.search(r"\d+", re.match(r"T\d+_", sample).group()).group()
    # Using this as a default value
    # Not sure if there's another default we might want
    return "1"


INDUCEOME_FP = Cfg["all"]["output_fp"] / "virus" / "induceome"
SBX_INDUCEOME_VERSION = open(get_induceome_path() / "VERSION").read().strip()
SBX_INDUCEOME_REF_FP = Path(Cfg["sbx_induceome"]["reference_fp"])
# There are a couple assumptions we're making in this mapping:
# 1. The reference file is named cd{number}.fas
# 2. The sample name is in the format T{number}_....fastq.gz
SBX_INDUCEOME_REF_MAP = {
    sample: SBX_INDUCEOME_REF_FP / f"cd{get_induceome_ref_var(sample)}.fas"
    for sample in Samples
}


try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


localrules:
    all_induceome,


rule all_induceome:
    input:
        expand(INDUCEOME_FP / "pileups" / "{sample}.pileup", sample=Samples),


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
