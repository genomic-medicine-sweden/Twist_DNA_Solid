__author__ = "Jonas Almlöf"
__copyright__ = "Copyright 2023, Jonas Almlöf"
__email__ = "jonas.almlof@scilifelab.uu.se"
__license__ = "GPL-3"


rule fuseq_wes:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        ref_json=config.get("fuseq_wes", {}).get("ref_json", ""),
        gtfSqlite=config.get("fuseq_wes", {}).get("gtfSqlite", ""),
        fusiondbFn=config.get("fuseq_wes", {}).get("fusiondbFn", ""),
        paralogdb=config.get("fuseq_wes", {}).get("paralogdb", ""),
    output:
        final_fusions=temp("fusions/fuseq_wes/{sample}_{type}/FuSeq_WES_FusionFinal.txt"),
        fusion_reads=temp("fusions/fuseq_wes/{sample}_{type}/feq_ALL.txt"),
        fusion_split_read_info=temp("fusions/fuseq_wes/{sample}_{type}/splitReadInfo.txt"),
        mate_pair1=temp("fusions/fuseq_wes/{sample}_{type}/FuSeq_WES_MR_fge.txt"),
        mate_pair2=temp("fusions/fuseq_wes/{sample}_{type}/FuSeq_WES_MR_fge_fdb.txt"),
        output_dir=temp(directory("{sample}_{type}")),
        split_read1=temp("fusions/fuseq_wes/{sample}_{type}/FuSeq_WES_SR_fge.txt"),
        split_read2=temp("fusions/fuseq_wes/{sample}_{type}/FuSeq_WES_SR_fge_fdb.txt"),
    log:
        "fusions/fuseq_wes/{sample}_{type}/FuSeq_WES_FusionFinal.txt.log",
    benchmark:
        repeat(
            "fusions/fuseq_wes/{sample}_{type}/FuSeq_WES_FusionFinal.txt.benchmark.tsv",
            config.get("fuseq_wes", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fuseq_wes", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fuseq_wes", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fuseq_wes", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fuseq_wes", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fuseq_wes", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fuseq_wes", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fuseq_wes", {}).get("container", config["default_container"])
    conda:
        "../envs/samtools.yaml"
    message:
        "{rule}: call dna fusion into {output.final_fusions} using bam file {input.bam}"
    shell:
        "sh -c \" . /venv/bin/activate && "
        "python3 fuseq_wes.py "
        "--bam {input.bam} "
        "--gtf {input.ref_json} "
        "--mapq-filter "
        "--outdir {output.output_dir} && "
        "process_fuseq_wes.R "
        "in={output.output_dir} "
        "sqlite={input.gtfSqlite} "
        "fusiondbFn={input.fusiondbFn} "
        "paralogdb={input.paralogdb} "
        "out={output.output_dir}\" &> {log}"
