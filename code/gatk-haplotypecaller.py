# Copied from https://snakemake-wrappers.readthedocs.io/en/0.51.3/wrappers/gatk/haplotypecaller.html
# to fix an issue for Joel Erberich with GATK ERC setting assuming GVCF.

import os

from snakemake.shell import shell

known = snakemake.input.get("known", "")
if known:
    known = "--dbsnp " + known

extra = snakemake.params.get("extra", "")
java_opts = snakemake.params.get("java_opts", "")
bams = snakemake.input.bam
if isinstance(bams, str):
    bams = [bams]
bams = list(map("-I {}".format, bams))

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "gatk --java-options '{java_opts}' HaplotypeCaller {extra} "
    "-R {snakemake.input.ref} {bams} "
    "-ERC BP_RESOLUTION "
    "-O {snakemake.output.gvcf} {known} {log}"
)
