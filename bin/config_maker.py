import glob
import os, sys

cwd = os.getcwd()
rawdata_path = cwd + '/rawdata'

input_sample_list = open("sample.list.txt", 'r')
output_config_file = open("dev.config.yaml", 'w')
output_config_file.write("samples:" + "\n")

for line in input_sample_list :
    split_line = line.rstrip("\r\n").split("\t")
    TBI_ID = split_line[0]
    CST_ID = split_line[1]

    sample_rawdata_path = rawdata_path + "/" + TBI_ID
    output_config_file.write("  Sample_{0}:".format(CST_ID) + "\n")

    file_list = glob.glob(sample_rawdata_path + "/*gz")
    sorted_file_list = sorted(file_list)

    lane_count = 0

    for item in sorted_file_list :
        if "R1" in item or "_1" in item :
            lane_count += 1
            output_config_file.write("    lane{0}:".format(lane_count) + "\n")
            fr_seq = item
        elif "R2" in item or "_2" in item :
            rr_seq = item
            output_config_file.write("      fq1: {0}".format(fr_seq) + "\n")
            output_config_file.write("      fq2: {0}".format(rr_seq) + "\n")
output_config_file.close()

output_cmd_file = open("dev.cmd.sh", 'w')
output_cmd_file.write('\n')
output_cmd_file.write('wkdir="{0}"\n'.format(cwd))
output_cmd_file.write('conffn="{0}/dev.config.yaml"\n'.format(cwd))
output_cmd_file.write('\n')
output_cmd_file.write("mkdir -p ${wkdir}/logs\n")
output_cmd_file.write('cmd=\"snakemake \\\n')
output_cmd_file.write("    --cores all \\\n")
output_cmd_file.write("    --printshellcmds \\\n")
output_cmd_file.write("    --snakefile /data/SJ/pipeline/onequeue_alignment/Snakefile \\\n")
output_cmd_file.write("    --cluster-config /data/SJ/pipeline/onequeue_alignment/envs/cluster.json \\\n")
output_cmd_file.write("    --max-jobs-per-second 5 \\\n")
output_cmd_file.write("    --jobs 192\n")
output_cmd_file.write("    --configfile ${conffn} \\\n")
output_cmd_file.write("    --directory ${wkdir} \\\n")
output_cmd_file.write("echo $cmd")
