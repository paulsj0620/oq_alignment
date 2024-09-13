import os
from bin.config import Utils
from snakemake.utils import min_version
min_version('6.1.1')

utils = Utils(config)
utils.add_config(config)
utils.get_paths(config)
utils.get_ref(config) 
utils.print_config(config) 


targets = list()
targets.append('bwa')
targets.append('gatk')
#targets.append('bedtools')
#targets.append('collectwgsmetrics')

rule all:
    input: utils.get_targets(targets)

include: 'modules/bwa.snakefile'
include: 'modules/gatk.snakefile'
include: 'modules/stats.snakefile'
#include: 'modules/bedtools.snakefile'
include: 'modules/picard.snakefile'
