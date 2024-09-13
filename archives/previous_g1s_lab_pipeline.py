import sys
#각 서버별 Reference 경로
#human === /data/ref/human/human_g1k_v37.fasta
#Pig ===== /data/ref/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
#mouse === /data/ref/mouse/GRCm38.fa
#med3 === /data1/ref/human/human_g1k_v37.fasta
#med3_T2T === /data2/new_ref_align/ref/chm13v2.0.re.fa
#240122_T2T === /data/ref/T2T/chm13v2.0.nochr.fa
#Macaque_T2T === /data/ref/crab_eating_macaque/T2T_MFA/T2T_MFA8v_good.fa

#각 서버별 GATK 3.5 경로
#med1 === /data/jaesoon/DB/gunhee/tools/gatk-3.5/GenomeAnalysisTK.jar
#med2 === /data/tools/gatk-3.5/GenomeAnalysisTK.jar
#med3 === /data1/tools/gatk-3.5/GenomeAnalysisTK.jar



### Name_list_file 경로 넣기###
lines_list = open('/content/drive/MyDrive/input/nslist.txt').read().splitlines()
#lines_list = ["Test_amp"]

### 커맨드 넣기 ###
base0 = "bwa mem -M -t 3 -R '@RG\\tID:input\\tLB:1\\tSM:input\\tPL:ILLUMINA' /data/ref/crab_eating_macaque/T2T_MFA/T2T_MFA8v_good.fa input_1.fq.gz input_2.fq.gz > input.sam"
base1 = "samtools view -Shu -@ 3 input.sam | samtools sort -@ 3 - > input.s.bam" #일반 샘플
#base1 = "samtools view -Shu -@ 15 input.sam | samtools sort -m 10G -@ 15 - > input.s.bam" #deep bulk용
base2 = "samtools index input.s.bam"
base3 = "picard MarkDuplicates -Xms8g -Xmx30g -XX:ConcGCThreads=1 REMOVE_DUPLICATES=false REMOVE_SEQUENCING_DUPLICATES=false I= input.s.bam O= input.s.md.bam M= input.matrics.txt VALIDATION_STRINGENCY=LENIENT TMP_DIR= input.markdup.tmp"
base4 = "samtools index input.s.md.bam"
base5 = "java -Xms8g -Xmx30g -jar /data/med3_data/DATA1/tools/gatk-3.5/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /data/ref/crab_eating_macaque/T2T_MFA/T2T_MFA8v_good.fa -I input.s.md.bam -o input.s.md.bam.intervals -nt 3"
base6 = "java -Xms8g -Xmx30g -jar /data/med3_data/DATA1/tools/gatk-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R /data/ref/crab_eating_macaque/T2T_MFA/T2T_MFA8v_good.fa -I input.s.md.bam -targetIntervals input.s.md.bam.intervals -o input.s.md.ir.bam"
base7 = "samtools index input.s.md.ir.bam"
base8 = "gatk HaplotypeCaller -R /data/ref/crab_eating_macaque/T2T_MFA/T2T_MFA8v_good.fa -I input.s.md.ir.bam -O input.vcf"
base9 = "bgzip -c input.vcf > input.vcf.gz"
base10 = "tabix -p vcf input"


n='if [ $? -ne 0 ] ; then\n\t   echo "sorry."\n\t   exit\nfi'

for j in range(len(lines_list)):

  final0 = ''
  final1 = ''
  final2 = ''
  final3 = ''
  final4 = ''
  final5 = ''
  final6 = ''
  final7 = ''
  final8 = ''
  final9 = ''
  final10 = ''


  final0= base0.replace('input', str(lines_list[j]))
  final1=base1.replace('input', str(lines_list[j]))
  final2=base2.replace('input', str(lines_list[j]))
  final3=base3.replace('input', str(lines_list[j]))
  final4=base4.replace('input', str(lines_list[j]))
  final5=base5.replace('input', str(lines_list[j]))
  final6=base6.replace('input', str(lines_list[j]))
  final7=base7.replace('input', str(lines_list[j]))
  final8=base8.replace('input', str(lines_list[j]))
  final9=base9.replace('input', str(lines_list[j]))
  final10=base10.replace('input', str(lines_list[j]))




  with open("/content/drive/MyDrive/output/"+lines_list[j]+".command.txt","w") as f:
    #print("#!/bin/bash",final1,n, final2,n, final3,n, final4,n, final5,n, final6,n, final7,"", sep='\n', file=f)
    print("#!/bin/bash", "#SBATCH -n 3", "#SBATCH -N 1", "#SBATCH -w med3", final0,n, final1,n, final2,n, final3,n, final4,n, final5,n, final6,n, final7,n, final8,n, final9,n, final10,"", sep='\n', file=f)