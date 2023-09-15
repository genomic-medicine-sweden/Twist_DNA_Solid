# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


include: "rules/common_references.smk"


rule all:
    input:
        unpack(compile_output_list),


# purecn will be fetching normaldb and intervals from the config
# and we want to use the newly generated files. There fore we overide
# the config point to files created by references pipeline.
config["purecn"]["normaldb"] = "references/purecn_normal_db/output/normalDB_hg19.rds"
config["purecn"]["intervals"] = "references/purecn_interval_file/targets_intervals.txt"
config['reference']['design_intervals'] = "references/preprocess_intervals/design.preprocessed.interval_list"
config['reference']['design_intervals_gatk_cnv'] = "references/preprocess_intervals/design.preprocessed.interval_list"


module pipeline:
    snakefile:
        "Snakefile"
    config:
        config


# The actual pipeline is imported to create files need to generate
# input files for PoN, Background, Artifacts, etc
use rule * from pipeline


use rule gatk_denoise_read_counts from pipeline with:
    input:
        hdf5PoN="references/create_read_count_panel_of_normals/gatk_cnv_panel_of_normal.hdf5",
        hdf5Tumor="cnv_sv/gatk_collect_read_counts/{sample}_{type}.counts.hdf5",


use rule cnvkit_batch from pipeline with:
    input:
        bam="alignment/samtools_merge_bam/{sample}_{type}.bam",
        bai="alignment/samtools_merge_bam/{sample}_{type}.bam.bai",
        cnv_reference="references/cnvkit_build_normal_reference/cnvkit.PoN.cnn",


module misc:
    snakefile:
        get_module_snakefile(config, "hydra-genetics/misc", path="workflow/Snakefile", tag="v0.2.0")
    config:
        config


module references:
    snakefile:
        github("hydra-genetics/references", path="workflow/Snakefile", tag="788c904")
    config:
        config


use rule * from references as references_*


# Ovveride input to use files generate by imported pipeline
# instead of fetching file paths from input files:
# The following files will be create by the pipeline:
# - alignment/samtools_merge_bam/{sample}_{type}.bam
# - annotation/background_annotation/{sample}_{type}.background_annotation.vcf.gz
# - cnv_sv/svdb_query/{sample}_{type}.pathology_purecn.svdb_query.vcf"
# - qc/add_mosdepth_coverage_to_gvcf/{sample}_{type}.mosdepth.g.vcf.gz"
# - references/collect_read_counts/{sample}_{type}.counts.hdf5
# - references/purecn_interval_file/targets_intervals.txt


####################################################
#              svdb
####################################################
# use vcf create by pipeline. ??????????Shoud we override svdb annotation??????????
use rule svdb_build from references as references_svdb_build with:
    input:
        cnv_vcfs=get_cnv_vcfs(units, "svdb"),


####################################################
#              artifact_panel
####################################################
# use vcf created by pipeline
use rule create_artifact_file from references as references_create_artifact_file with:
    input:
        vcfs=get_vcfs(units, "artifact"),


####################################################
#              gatk pon
####################################################
# use hdf5 files created by reference pipeline references/collect_read_counts/%s_%s.counts.hdf5
use rule create_read_count_panel_of_normals from references as references_create_read_count_panel_of_normals with:
    input:
        bams=get_hdf5(units, "gatk_pon"),


# Use bam files created by pipeline: alignment/samtools_merge_bam/{sample}_{type}.bam
use rule collect_read_counts from references as references_collect_read_counts with:
    input:
        bam=lambda wildcards: "alignment/samtools_merge_bam/%s_%s.bam" % (wildcards.sample, wildcards.type),
        bai=lambda wildcards: "alignment/samtools_merge_bam/%s_%s.bam.bai" % (wildcards.sample, wildcards.type),
        interval="references/preprocess_intervals/design.preprocessed.interval_list",


use rule preprocess_intervals from references as references_preprocess_intervals with:
    output:
        temp("references/preprocess_intervals/design.preprocessed.interval_list"),


####################################################
#              cnvkit pon input override
####################################################
# use bam files create by pipeline: alignment/samtools_merge_bam/{sample}_{type}.bam
use rule cnvkit_build_normal_reference from references as references_cnvkit_build_normal_reference with:
    input:
        bams=get_bams(units, "cnvkit_pon"),
        target="references/cnvkit_create_targets/cnvkit_manifest.target.bed",
        antitarget="references/cnvkit_create_anti_targets/cnvkit_manifest.antitarget.bed",
        ref=config.get("reference", {}).get("fasta", ""),
        mappability=config.get("reference", {}).get("mappability", ""),


####################################################
#              purecn normal input over
####################################################
# Use bam files created by pipeline: alignment/samtools_merge_bam/{sample}_{type}.bam
use rule purecn_bam_list from references as references_purecn_bam_list with:
    input:
        bam_list=get_bams(units, "purecn_mapping_bias"),


use rule bcftools_merge from references as references_bcftools_merge with:
    input:
        vcfs=get_vcfs(units, "purcen_mapping_bias"),
        vcfs_tabix=expand("{dataset}.{ext}", dataset=get_vcfs(units, "purecn_mapping_bias"), ext=["tbi"]),


# Uses to create background
use rule create_background_file from references as references_create_background_file with:
    input:
        gvcfs=get_gvcfs(units, "background"),


# Make use of new interval file created by the pipeline, and wait for it to be created
use rule purecn_coverage from references as references_purecn_coverage with:
    input:
        bam_list_file="references/purecn_bam_list/bam_files.list",
        intervals="references/purecn_interval_file/targets_intervals.txt",
    params:
        intervals="references/purecn_interval_file/targets_intervals.txt",
        extra=config.get("purecn_coverage", {}).get("extra", ""),
