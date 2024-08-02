#!/bin/bash -l
#SBATCH -A pawsey0012
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=48GB
#SBATCH --partition=long
#SBATCH --time=96:00:00
#SBATCH --mail-user=sarah.beecroft@csiro.au
#SBATCH --mail-type=ALL

alpID=alp016
fastq1=alp016_S2_L004_R1_001.fastq.gz
fastq2=alp016_S2_L004_R2_001.fastq.gz

# Syntax= bash my_name.sh alpID path_to/fastq1 path_to/fastq2
# run from the project's directory in the new server, inside the "working" directory.
# check that a bowtie2 index, a .fai and a .dict file have been made for the reference, and are in the same folder as the reference.
source /software/projects/pawsey0001/sbeecroft/miniforge3/bin/activate alpaca
module load bowtie2/2.4.5--py36hd4290be_0 samtools/1.15--h3843a85_0 gatk4/4.2.5.0--hdfd78af_0

## NOTES
# $alpID must be the seq ID prefix. This script is for one sample at a time, not even a batch
# not sure what the other inputs are yet
# could maybe use NF modules for some of these steps, but I think straight up Sarek might not be the right option here for Kylie

# Makes a directory for gvcf_links which won't be overwritten every time the script is run, also, it creates a directory specific for the animal, and subdirectories for each part of the process
mkdir -p gvcf_links $alpID/{qc_reads,mapped,report,gvcf};

# FastP removes the poly-G tail caused by the NovaSeq 6000 which uses 2-colour chemistry. To use with normal illumina, remove the --trim_poly_g and --poly_g_min_len options.
echo "QC & Trimming started"; date; 
fastp --thread 16 \
    -h $alpID/report/${alpID}_fastp.html \
    -i $fastq1 \
    -I $fastq2 \
    -o $alpID/qc_reads/${alpID}_1.fastq.gz \
    -O $alpID/qc_reads/${alpID}_2.fastq.gz &&
echo "QC & Trimming complete"; date;

# Mapping reads to the reference genome using Bowtie2, with local alignment set to allow identificaiton of e.g. Tn insertions, and with unmapped reads being included in the output file.
echo "Mapping started"; date;
bowtie2 --local -x normalized_VicPac4Dic \
    -p 12 \
    -1 $alpID/qc_reads/${alpID}_1.fastq.gz \
    -2 $alpID/qc_reads/${alpID}_2.fastq.gz | \
    samtools view -@ 16 -bS - | \
    samtools view -@ 16 -bT normalized_VicPac4Dic.fasta - | \
    samtools sort -@ 16 -m 3.6GB T /tmp/$SLURM_JOB_ID/bowtie2 - -o $alpID/mapped/${alpID}_q20.bam
echo "Mapping complete"; date;

echo "Creating index 1"; date;
samtools index $alpID/mapped/${alpID}_q20.bam $alpID/mapped/${alpID}_q20.bai
echo "Index 1 created"; date;

# Marking duplicate reads that are PCR artefacts that are created during library construction, limiting bias at the variant calling step
echo "Marking duplicates started for $alpID"; date;
gatk MarkDuplicatesSpark \
    CREATE_INDEX=TRUE \
    ASSUME_SORT_ORDER=coordinate \
    QUIET=TRUE \
    -I $alpID/mapped/${alpID}_q20_rg.bam \
    -O $alpID/mapped/${alpID}_q20_rg_d.bam \
    -M $alpID/report/${alpID}_q20_rg_d_metrics.txt \
    --tmp-dir /tmp/$SLURM_JOB_ID/markduplicates \
	--conf 'spark.executor.cores=16'
echo "Duplicates marked for $alpID"; date;

# Initial variant calling for the individual animal - output is a genomic variant call format (GVCF) file
echo "Calling variants in $alpID"; date;
gatk HaplotypeCaller \
    -R normalized_VicPac4Dic.fasta \
    -I $alpID/mapped/${alpID}_q20_rg_d.bam \
	-ERC GVCF \
    -O $alpID/gvcf/$alpID.g.vcf \
	--tmp-dir /tmp/$SLURM_JOB_ID/haplotypecaller
echo "Variant call file created for $alpID"; date;

# Generates a hard link of the gvcf in the gvcf_links file for easy merging of all animals gvcfs using CombineGVCFs
ln $alpID/gvcf/$alpID.g.vcf gvcf_links/ &&
echo "$alpID complete" 
