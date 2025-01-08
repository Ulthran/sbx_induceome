import csv


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
        expand(
            INDUCEOME_FP / "pileups" / "{sample}.pileup", sample=SBX_INDUCEOME_SAMPLES
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
    output:
        peaks_img=INDUCEOME_FP / "peaks" / "{sample}.png",
        peaks_csv=INDUCEOME_FP / "peaks" / "{sample}.csv",
    log:
        LOG_FP / "induceome_find_peaks_{sample}.log",
    benchmark:
        BENCHMARK_FP / "induceome_find_peaks_{sample}.tsv"
    params:
        prominence=Cfg["sbx_induceome"]["prominence"],
        distance=Cfg["sbx_induceome"]["distance"],
    conda:
        "envs/sbx_induceome_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_induceome:{SBX_INDUCEOME_VERSION}"
    script:
        "scripts/induceome_find_peaks.py"