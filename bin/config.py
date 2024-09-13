

import os
import sys
import yaml
import json
import subprocess

import pandas as pd
from collections import defaultdict


class Utils:
    def __init__(self, config):
        self.config = config

    def get_targets(self, targets):
        ls = list()
        for sample, sample_dic in self.config['samples'].items():
            ls.append(f"analysis/bwa/{sample}/{sample}.sam")
            ls.append(f"analysis/bwa/{sample}/{sample}.s.bam")
            ls.append(f"analysis/bwa/{sample}/{sample}.s.bam.bai")
            ls.append(f"analysis/bwa/{sample}/{sample}.s.md.bam")
            ls.append(f"analysis/bwa/{sample}/{sample}.metrics.txt")
            ls.append(f"analysis/bwa/{sample}/{sample}.s.md.bam.bai")
            ls.append(f"analysis/gatk/{sample}/{sample}.s.md.bam.intervals")
            ls.append(f"analysis/gatk/{sample}/{sample}.s.md.ir.bam")
            ls.append(f"analysis/gatk/{sample}/{sample}.s.md.ir.bam.bai")
            ls.append(f"analysis/gatk/{sample}/{sample}.vcf")
            ls.append(f"analysis/gatk/{sample}/{sample}.vcf.gz")
            ls.append(f"analysis/gatk/{sample}/{sample}.vcf.gz.tbi")
        return ls

    def get_ref(self, config):
        infh = open(config['ref'])
        ref_dic = yaml.safe_load(infh)
        for key, value in ref_dic.items():
            key = f'ref_{key}'
            config[key] = value
        return config

    def get_paths(self, config):
        conda_root = subprocess.check_output('conda info --root', shell=True).decode('utf-8').strip()
        if not "conda_bin_path" in config or not config["conda_bin_path"]:
            config["conda_bin_path"] = os.path.join(conda_root, 'envs', 'snakemake', 'bin')
        return config

    def add_config(self, config):
        config["od_smp_s"] = [sample for sample in config["samples"]]
        return config

    def print_config(self, config):
        print("="*100)
        for key, value in sorted(config.items()):
            if not value:
                pass
                #continue
            if key in ['samples']:
                for subkey, subvalue in sorted(value.items()):
                    for subsubkey, subsubvalue in sorted(subvalue.items()):
                        print(f"{key} : {subkey} : {subsubkey} : {subsubvalue}")
                #print(f"{key} : {value}")
            else:
                print(f"{key} : {value}")
        print("="*100)
        return config

