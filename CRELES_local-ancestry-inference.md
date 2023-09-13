---
title: "CRELES: Local Ancestry"
date: "September 2021"
author: "Paola Arguello Pascualli"
output: 
  html_document: 
    keep_md: true
---

# Getting ready
## Accessing server files on cluster


```bash

## Mounting the KoborLab nodes (hpc01/02)
bash mounting_server.sh

## I've also created a symbolic link to my home directory within the Kobor lab server
cd KoborSpace ## /mnt/cifs/username/fs/KoborLab/kobor_space/ppascualli 

```


## Installing packages on cluster

All the programs that we need should be on the /mnt/common/KoborLab space within the cluster. This is to be able to access them without any issues when running on the cluster. Nothing should be called from the long-term storage (KoborSpace: /mnt/cifs/username/fs/KoborLab/kobor_space/). 


```bash

## Anaconda ---------------
cd ## working on home directory 
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh 
## PATH: /home/BCRICWH.LAN/Paola.Arguello/programs/anaconda3/bin


## ShapeIT ---------------
wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz
tar -zxvf shapeit.v2.r900.glibcv2.17.linux.tar.gz
## PATH: /mnt/common/KoborLab/ShapeIT/bin/shapeit


## XGMix -----------------
wget https://github.com/AI-sandbox/XGMix/archive/refs/heads/master.zip
pip install -r requirements.txt

```


# 1. Phasing (ShapeIT)

Version: 
 - *v2.r904*

In order to complete the Local Ancestry pipeline we need to phase the quality controlled CRELES data (CRELES pc7). Humans are diploid organisms which means we have two copies of each allele. Genotyping and sequencing technologies allow us to have information about both copies, however, we cannot know the complete parental haplotypes. Phasing uses information from genetic maps collected from hundreds or thousands of WGS genomes to phase the genotyped data. 


```bash

## Copying the genetic maps from Santigo to the cluster:
cd /mnt/scratch/KoborLab/Paola/Local_Ancestry/genetic_maps
wget https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b37.tar.gz
gzip -d genetic_maps.b37.tar.gz | tar -xv
wget https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz
gzip -d genetic_maps.b38.tar.gz | tar -xv


## Running the pashing algorith for all the chromosomes:
cd /mnt/scratch/KoborLab/Paola/Local_Ancestry/CRELES_pc7/shapeit_results

#### SLURM script -----------------------

#!/bin/bash
#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=paola.arguello@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage: 64 Gb of RAM for the whole job; using 8 CPUs (think 8G RAM/CPU)
#SBATCH --mem=64G

## Using 8 CPUs
#SBATCH --cpus-per-task=8

## Running for a max time of 48 hours
#SBATCH --time=48:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

pwd; hostname; date

for i in {1..23}
do
    pwd; date 
    echo "Starting phasing on Chr $i"
    
    /mnt/common/KoborLab/ShapeIT/bin/shapeit \
    --input-vcf /mnt/scratch/KoborLab/Paola/Local_Ancestry/CRELES_pc7/VCF/CRELES_pc7-updated-chr$i.vcf.gz \
    --input-map /mnt/scratch/KoborLab/Paola/Local_Ancestry/genetic_maps/b37/chr$i.b37.gmap.gz \
    --output-max CRELES_pc7up-chr$i\_phased.haps CRELES_pc7up-chr$i.sample \
    --thread 3 \
    --seed 12345
    
    echo "Converting to vcf format..."
    
    /mnt/common/KoborLab/ShapeIT/bin/shapeit  -convert \
    --input-haps CRELES_pc7up-chr$i\_phased \
    --output-vcf CRELES_pc7up-chr$i\_phased.vcf
    
    echo "Done"; date
    
done

#### ------------------------------------

sbatch phasing_cluster.sh
squeue ## check on the job status

```


# 2. Liftover genome build b37 to b38

The trained models that we will use to infer the LA are in the hg38 build of the human genome, while CRELES is in hg37. Thus, we will convert the genomic coordinates from one build to another. 

## Collapse chromosomes files (bcftools)

Version: 

 - *bcftools 1.9*
 
 - *htslib 1.9*

It is not common but sometimes some positions change coordinates to another chromosome. In order to avoid issues when changing the build is recommended to have all the chromosomes in one file.


```bash

for i in {1..22}
do
  #bgzip CRELES7-chr$i.vcf
  tabix -p vcf CRELES7-chr$i.vcf.gz
done

bcftools concat -o CRELES7_all-chr_b37.vcf \
CRELES7-chr1.vcf.gz CRELES7-chr2.vcf.gz \
CRELES7-chr3.vcf.gz CRELES7-chr4.vcf.gz \
CRELES7-chr5.vcf.gz CRELES7-chr6.vcf.gz \
CRELES7-chr7.vcf.gz CRELES7-chr8.vcf.gz \
CRELES7-chr9.vcf.gz CRELES7-chr10.vcf.gz \
CRELES7-chr11.vcf.gz CRELES7-chr12.vcf.gz \
CRELES7-chr13.vcf.gz CRELES7-chr14.vcf.gz \
CRELES7-chr15.vcf.gz CRELES7-chr16.vcf.gz \
CRELES7-chr17.vcf.gz CRELES7-chr18.vcf.gz \
CRELES7-chr19.vcf.gz CRELES7-chr20.vcf.gz \
CRELES7-chr21.vcf.gz CRELES7-chr22.vcf.gz 

rm *.tbi
```


## Liftover 

Version: 

- Using conda environment located under: 


```bash

cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/2_liftOver_37to38/


# Lift over vcf files to GRCh38
# I use the programs crossmap
# documentation: http://crossmap.sourceforge.net/#convert-vcf-format-files
# Santiago created a conda enviroment: encironments/crossmap.yaml 

## download_GRCh38_genome -------------------------------------
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

## download_ucsc_chain_hg19ToHg38 -----------------------------
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

## lift_over_to_GRCh38 ----------------------------------------
conda activate crossmap
bcftools view ../../1_Phasing_shapeit/CRELES_all_phasing_results_Sb37/CRELES7_all-chr_b37.vcf > CRELES-ALL-GRCh37.vcf
## Crossmap.py file_format, chainfile, input.vcf, refGenome.fa, output_file
CrossMap.py vcf ../../lift_reference_chains/hg19ToHg38.over.chain.gz ./CRELES-ALL-GRCh37.vcf ../../lift_reference_chains/GRCh38.primary_assembly.genome.fa ./CRELES-ALL-GRCh38.vcf > lift-over-to-GRCh38.log 
rm -f CRELES-ALL-GRCh37.vcf

## sort_lifted_and_index ---------------------------------------
# We sort compress and index the vcf
bcftools sort CRELES-ALL-GRCh38.vcf -Oz -o CRELES-ALL-GRCh38-sorted.vcf.gz
bcftools index -t CRELES-ALL-GRCh38-sorted.vcf.gz

## lifted_split_by_chromosome ----------------------------------
for chr in {1..22}
do 
   bcftools view -r $chr CRELES-ALL-GRCh38-sorted.vcf.gz -Oz -o CRELES-chr$chr-GRCh38.vcf.gz
   bcftools index CRELES-chr$chr-GRCh38.vcf.gz --tbi
done

```


# 3. XGMix 

*~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/*

XGMix (https://github.com/AI-sandbox/XGMix) is a program that estimates local ancestry using Machine Learning methods. Santiago Medina shared the 1KGP populations he used to trained his models (outlined below). I will use the same samples but instead of using the Whole Genome Sequencing (WGS) data I will use the overlap of the WGS data and CRELES variants. 

  - EUR: IBR + GBR
  
  - NAT: PEL and MXL that were shown to be representative for NAT ancestry
  
  - AFR: YRI


## 3.0 Modify genetic map


```bash
## XGMix expects a tsv file that has the following columns: chr,pos,pos_cm
## here we reformat the genetic map to be used with XGmix

cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix/genetic_maps

## re-rder and no header
for chrN in {1..22}
do 
    bgzip -d ../genetic_maps/b38/chr$chrN.b38.gmap.gz
    awk -v OFS='\t' '{{print $2, $1, $3}}' ../genetic_maps/b38/chr$chrN.b38.gmap | grep -v '^chr' > chr$chrN.b38.gmap.txt
done

```


## 3.1 1KGP-WGS and CRELES intersection

In order to maximize the best local ancestry inference we will use the intersection of the whole genome sequencing (WGS) data from the 1KGP (*https://www.biorxiv.org/content/10.1101/2021.02.06.430068v1*) and the genotyped data from CRELES (~430K genetic variants) and then we will train XGMix with the intersection. 


```bash
## bcftools version
#bcftools 1.9
#Using htslib 1.9
for i in {3..22}
do
    ##Creating index for reference files to be able to get the intersection 
    bcftools index ../references/OneTGP-references-WGS/vcfs/1TGP-chr$i-snps-vep-mask-HW-GRCh38.vcf.gz 
    ##Extract and write records from A shared by both A and B using exact allele match
    bcftools isec -p chr$i -n=2 -w1 \
    ../references/OneTGP-references-WGS/vcfs/1TGP-chr$i-snps-vep-mask-HW-GRCh38.vcf.gz \
    ../2_liftOver_37to38/CRELES_all-chr_update/CRELES-chr$i-GRCh38.vcf.gz
    ## Get number of overlapping variants:
    grep -v "^#" ./chr$i/0000.vcf | wc -l > chr$i-overlaping_variants.txt
done

## Number of overlapping variants: 354,976
#sum(c(28078,
#29380,
#25265,
#23523,
#21464,
#26190,
#19954,
#17920,
#15560,
#18230,
#17448,
#17264,
#12791,
#11623,
#11164,
#11486,
#10195,
#10747,
#7604,
#8964,
#4977,
#5149))


## Reformatting output to different directories
for i in {3..22}
do
  echo "Starting changing names from Chr $i"
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/chr$i/0000.vcf ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/chr$i/chr$i-1kpg_CR354K.vcf
  
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/chr$i/sites.txt ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/chr$i/sites_chr$i.txt
  
  echo "Moving files"
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/chr$i/chr$i.vcf ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter
  
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/chr$i/sites_chr$i.txt ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter
  
  echo "Removing README.txt"
  rm ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/chr$i/README.txt 
  
done 

```


## 3.2 Create reference and samples map files

For the reference map I only took the information from: 

- European: GBR + IBS (~200)
- African: YRI (~100)
- Native American: MXL+PEL (pre-selected based on their global ancestry (~30) *https://github.com/santiago1234/mxb-genomes/blob/main/resources/1TGP-samples-meta-data/native-american.txt*


```bash
###### Reference ---------------------------------------------
cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_GSA_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter

## From this folder, I took the sample names from the header of Chr8 vcf file
perl -pe "s/ /\n/g" samples.txt > samples1
mv samples1 samples.txt

## Now we get the information from those samples:
cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/genotyping_references/1KGP
cp 1KGP_samp-info.txt ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/2_map_files/needed_files
grep -f 1kgp_samp.txt 1KGP_samp-info.txt | cut -f1,3-4 > reference_map_file.txt
rm samples.txt 1KGP_samp-info.txt 

###### Samples -----------------------------------------------
cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/2_liftOver_37to38/CRELES_all-chr_update

## From this folder, I took the names from the header of the CRELES-ALL-GRCh38-sorted.vcf.gz file
perl -pe "s/\s+/\n/g" samples.txt > samples1
mv samples1 samples.txt
vim samples.txt ## some extra manual work was needed due to the complexity of the names that was not captured by the regular expression

### Getting only relevant information from CRELES pop meta (ID, Pop)
cut -f2,3 CRELES_pop_meta.txt | perl -pe "s/CRELES - Non Nicoyan/Nicoyan/g; s/CRELES - Nicoyan/Non-Nicoyan/g" > CRELES_info.txt

##### Merging files information:(R) ----------------------------
CRELES <- read.table("~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/2_map_files/CRELES_info.txt", header = TRUE)
samples <- read.table("~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/2_map_files/samples.txt", header = FALSE)

samples$IID <- unlist(strsplit(as.character(samples$V1 ), "_"))[seq(2,972,2)]
samples_map <- merge(samples, CRELES, by="IID")
samples_map$Population <- "CRELES"
colnames(samples_map) <- c("IID", "Sample", "Pop", "Population")

write.table(samples_map, "~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/2_map_files/sample_map_file_ALL-info.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
```




## 3.3 Merge reference and CRELES


```bash

## CRELES: -------------------------------------------------------------------------------------------------------------------------------------
## ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/2_liftOver_37to38/CRELES_all-chr_update/CRELES-ALL-GRCh38-sorted.vcf.gz

## 1KGP : -------------------------------------------------------------------------------------------------------------------------------------
# Need to get them all into one file
for i in {1..22}
do
  #bgzip chr$i-1kpg_CR354K.vcf
  tabix -p vcf chr$i-1kpg_CR354K.vcf.gz
done

bcftools concat chr1-1kpg_CR354K.vcf.gz chr2-1kpg_CR354K.vcf.gz chr3-1kpg_CR354K.vcf.gz   chr4-1kpg_CR354K.vcf.gz chr5-1kpg_CR354K.vcf.gz chr6-1kpg_CR354K.vcf.gz chr7-1kpg_CR354K.vcf.gz   chr8-1kpg_CR354K.vcf.gz chr9-1kpg_CR354K.vcf.gz chr10-1kpg_CR354K.vcf.gz   chr11-1kpg_CR354K.vcf.gz   chr12-1kpg_CR354K.vcf.gz   chr13-1kpg_CR354K.vcf.gz chr14-1kpg_CR354K.vcf.gz   chr15-1kpg_CR354K.vcf.gz   chr16-1kpg_CR354K.vcf.gz   chr17-1kpg_CR354K.vcf.gz   chr18-1kpg_CR354K.vcf.gz   chr19-1kpg_CR354K.vcf.gz chr20-1kpg_CR354K.vcf.gz   chr21-1kpg_CR354K.vcf.gz   chr22-1kpg_CR354K.vcf.gz  -Oz -o all-chr_kpg_CR354K.vcf.gz

## ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/all-chr_kpg_CR354K.vcf.gz


###### Merging ----------------------------------
cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/3_Merge_ref_CRELES
cp ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/1_CRELES-1kgp_inter/all-chr_kpg_CR354K.vcf.gz ./
cp ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/2_liftOver_37to38/CRELES_all-chr_update/CRELES-ALL-GRCh38-sorted.vcf.gz ./

tabix -p vcf all-chr_kpg_CR354K.vcf.gz
tabix -p vcf CRELES-ALL-GRCh38-sorted.vcf.gz

bcftools merge all-chr_kpg_CR354K.vcf.gz CRELES-ALL-GRCh38-sorted.vcf.gz -Oz -o merged_1KGP-CR357K_CRELES.vcf.gz

```


## 3.4 Subsetting 

```bash

cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/1kgp-WGS_CRELES_overlap/4_Subsetting

grep -v "Sample" reference_map_file.txt | cut -f1 > references_IDs.txt
grep -v "Sample" sample_map_file.txt | cut -f1 > CRELES_IDs.txt

bcftools view merged_1KGP-CR357K_CRELES.vcf.gz -S CRELES_IDs.txt -Oz -o CRELES_LAI_samples.vcf.gz
bcftools view merged_1KGP-CR357K_CRELES.vcf.gz -S references_IDs.txt -Oz -o 1KGP_LAI_samples.vcf.gz

tabix -p vcf CRELES_LAI_samples.vcf.gz
tabix -p vcf 1KGP_LAI_samples.vcf.gz

### Splitting vcf files again ---------------------------------------
for i in {1..22}
do 
  echo "Starting Chromosome $i"
  
  ## For CRELES
  echo "CRELES subset chr$i"
   bcftools view -r $i ../4_Subsetting/CRELES_LAI_samples.vcf.gz -Oz -o CRELES_LAI_samples_chr$i.vcf.gz
   bcftools index CRELES_LAI_samples_chr$i.vcf.gz --tbi
   
   ## For 1KGP
   echo "1KGP subset chr$i"
   bcftools view -r $i ../4_Subsetting/1KGP_LAI_samples.vcf.gz -Oz -o 1KGP_LAI_samples_chr$i.vcf.gz
   bcftools index 1KGP_LAI_samples_chr$i.vcf.gz --tbi
   
done

```


## 3.5 Run/train XGMix

Santiago Medina, a graduate student in Andres' Moreno lab shared his pre-trained model for XGMix. This model was build using the following reference panels taken primarily from the One-thousand reference project (1KGP) and some samples from the Mexican Biobank (MXB): 

  - EUR: IBR + GBR
  
  - NAT: 50 MXB + 1KGP (PEL and MXL that were shown to be representative for NAT ancestry)
  
  - AFR: YRI



```bash


## Copy files to cluster: /mnt/scratch/KoborLab/Paola/

##### xgmix_CR357K_script.sh -----------------------------------------------------------

#!/bin/bash
#SBATCH --partition=defq

## Change to be your email address
#SBATCH --mail-user=paola.arguello@bcchr.ca
#SBATCH --mail-type=ALL

## CPU Usage: 64 Gb of RAM for the whole job; using 8 CPUs (think 8G RAM/CPU)
#SBATCH --mem=128G

## Using 8 CPUs
#SBATCH --cpus-per-task=8

## Running for a max time of 48 hours
#SBATCH --time=72:00:00

## Using only a single node
#SBATCH --nodes=1

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

pwd; hostname; date
conda init bash
conda activate xgmix

for i in {1..22}
do
    ## The parameters required to run XGMix are outlined below in the same order as they are required:
    ## query_file             = ../5_Splitting_chr/CRELES_LAI_samples_chr$i.vcf.gz
    ## genetic_map_file       = ../0_Modify_genetic_maps/chr$chrN.b38.gmap.txt
    ## output_basename        =  CRELES-CR357K_LAI-3pop_chr$i
    ## chr_nr                 = $chrN
    ## phase                  = True
    ## reference_file         = ../5_Splitting_chr/1KGP_LAI_samples_chr$i.vcf.gz
    ## sample_map_file        = ../2_map_files/reference_map_file.txt
    ## logs                   = chr$i-CRELES_CR357K_LAI_results.log
    
    echo "Starting LAI for chr $chrN"
    python /mnt/common/KoborLab/XGMix/XGMIX.py ../5_Splitting_chr/CRELES_LAI_samples_chr$i.vcf.gz ../0_Modify_genetic_maps/chr$i.b38.gmap.txt CRELES-CR357K_LAI-3pop_chr$i $i True  ../5_Splitting_chr/1KGP_LAI_samples_chr$i.vcf.gz  ../2_map_files/reference_map_file.txt  >  chr$i-CRELES_CR357K_LAI_results.log                                                 
    
done


###### -----------------------------------------------------------------
## sbatch xgmix.sh
## squeue

## Moving outputs to new directories
## First we move from the cluster:
## cp -r 6_XGMix /mnt/cifs/paola.arguello/fs/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results

for i in {1..22}
do 
  ## Results 
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES-CR357K_LAI-3pop_chr$i/CRELES-CR357K_LAI-3pop_chr$i.fb.tsv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES_pop3_xgmix_results
  
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES-CR357K_LAI-3pop_chr$i/CRELES-CR357K_LAI-3pop_chr$i.msp.tsv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES_pop3_xgmix_results
  
  ## Rephased vcfs
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES-CR357K_LAI-3pop_chr$i/query_file_phased.vcf ./CRELES-CR357K_LAI-3pop_chr$i/CRELES_rephased_chr$i.vcf
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES-CR357K_LAI-3pop_chr$i/CRELES_rephased_chr$i.vcf ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/rephased_xgmix_vcfs
  
  ## Models
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES-CR357K_LAI-3pop_chr$i/models/model_chm_$i/model_chm_$i.pkl ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES-CR357K_LAI-3pop_chr$i/models/model_chm_$i/model_CR357K_chr$i.pkl
  
  mv ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES-CR357K_LAI-3pop_chr$i/models/model_chm_$i/model_CR357K_chr$i.pkl ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/models_geno_CR357K
  
done


#### Sanity Check (R) ---------------------------------------------------------

chr.sanity <- read.table("~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_genotyping_training_3pop_Results/CRELES_pop3_xgmix_results/CRELES-CR357K_LAI-3pop_chr10.msp.tsv", header = FALSE)

seq(7,938,1)
seq(8,938,1)

comparing.haps <- function(n1,n2){
  
  p1 <- as.data.frame(table(chr.sanity[,n1]))
  p2 <- as.data.frame(table(chr.sanity[,n2]))
  
  ## For haplotype 1 ---------------
  AFR <- p1[p1$Var1 == 0,2]
  if(length(AFR) == 0){AFR <- 0}
  
  AMR <- p1[p1$Var1 == 1,2]
  if(length(AMR) == 0){AMR <- 0}
  
  EUR <- p1[p1$Var1 == 2,2]
  if(length(EUR) == 0){EUR <- 0}
  
  total <- sum(AFR,AMR,EUR)
  
  AFR.per1 <- (AFR*100/total)
  AMR.per1 <- (AMR*100/total)
  EUR.per1 <- (EUR*100/total)
  
  ### For haplotype 2 ---------------
  AFR <- p2[p2$Var1 == 0,2]
  if(length(AFR) == 0){AFR <- 0}

  AMR <- p2[p2$Var1 == 1,2]
  if(length(AMR) == 0){AMR <- 0}
  
  EUR <- p2[p2$Var1 == 2,2]
  if(length(EUR) == 0){EUR <- 0}
  
  total <- sum(AFR,AMR,EUR)
  
  AFR.per2 <- (AFR*100/total)
  AMR.per2 <- (AMR*100/total)
  EUR.per2 <- (EUR*100/total)
  
  data.frame(hap1=c(AFR.per1,AMR.per1,EUR.per1),
             hap2=c(AFR.per2,AMR.per2,EUR.per2),
             diff=c(AFR.per1,AMR.per1,EUR.per1) - c(AFR.per2,AMR.per2,EUR.per2))
  
  
}

comparing.haps(7,8)
comparing.haps(9,10)
comparing.haps(11,12)
comparing.haps(13,14)
comparing.haps(15,16)
comparing.haps(17,18)
comparing.haps(19,20)
comparing.haps(21,22)
comparing.haps(23,24)
comparing.haps(25,26)
comparing.haps(27,28)

```


# 4. Output modifications
Adapting the snakefile code written by Santiago to be more user friendly. His original code can be found here: **https://github.com/santiago1234/mxb-genomes/blob/main/analysis-doc/210514-ProcessXGMixOutForTracts/Snakefile**

## 4.1 Merge xgmix 


```bash

cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_XGMix_2_Tracts

# Merge the XGMix msp files to have one for all the genome
## python module from Santi's repo copied here:  ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/environments/mxbgenomes

python 1_combine-chromosomes-xgmix.py  "/home/BCRICWH.LAN/paola.arguello/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/3_XGMix_genotyping_training/6_XGMix_results/CRELES_pop3_xgmix_results/msp_files/" "/home/BCRICWH.LAN/paola.arguello/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_XGMix_2_Tracts/1_CRELES3_all_msp.csv"

```

## 4.2 Remove low coverage tracts


```bash

python 2_remove_low_cov_windows.py 1_CRELES3_all_msp.csv 22 2_CRELES3_all_msp_lowCovClean.csv

```

## 4.3 Process admixed samples

```bash

#### Sanity Check: ----------
python 3_process-individuals.py 2_CRELES3_all_msp_lowCovClean.csv ./3_bed_files/ 233_102-R2
python 3_process-individuals.py 2_CRELES3_all_msp_lowCovClean.csv ./3_bed_files/ 338_1031

### Loop --------------------

while read ID; do
  python 3_process-individuals.py 2_CRELES3_all_msp_lowCovClean.csv ./3_bed_files/ $ID
done < 3_CRELES_IDs.txt

```


## 4.4 Tracts input ready


```bash

### Getting the txt file containing all the files: --------
cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_XGMix_2_Tracts/3_bed_files
ls | cut -f5 > ../4_Tracts_input_ready/4_bed-files.txt

### Processing the output: ---------------------------------
while read samplebed; do
  # extract the sample name from the input bed file
  sample_name=$(basename -s .bed $samplebed)
  echo "processing: $sample_name"
  python 4_make_tracts_input.py ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_XGMix_2_Tracts/3_bed_files/$samplebed $sample_name ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_XGMix_2_Tracts/4_Tracts_input_ready
done < 4_bed-files.txt

```


# 5. Karyograms

Located in the following directory: **~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/5_Karyograms**


```r
library(tidyverse)
library(scales)
library(ggthemes)

theme_set(theme_tufte(base_family = "Helvetica"))

save.karyograms <- function(sample){
  
  input_bed <- paste0("~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_XGMix_2_Tracts/3_bed_files/",sample, ".bed")
  output_plot <- paste0("~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/5_Karyograms/karyo", sample,".png")
  bed <- read.table(input_bed, header = TRUE, sep= "\t")
  
  individual <- basename(input_bed) %>%
    str_replace(".bed", "")
  
  
  #### Drawing an indicator variable for Hapl  otypes limits ------------------------
  limit_y_max <- function(haplo) {
    if (haplo == "A") {
      return(Inf)
    }
    if (haplo == "B") {
      return(0.5) #0.49
    }
  }
  
  limit_y_min <- function(haplo) {
    if (haplo == "A") {
      return(0.5)
    }
    if (haplo == "B") {
      return(-Inf)
    }
  }
  
  bed <- bed %>% 
    mutate(ymin_haplo = map_dbl(Haplotype, limit_y_min),
           ymax_haplo = map_dbl(Haplotype, limit_y_max),)
  
  #### Creating the plot: ---------------------------------------
  bed %>% 
    ggplot(aes(x = spos, fill = Ancestry)) +
    geom_rect(aes(xmin = spos, xmax=epos, ymin = ymin_haplo, ymax = ymax_haplo)) +
    geom_hline(yintercept = 0.5, color ="white", size = 0.1) +
    facet_grid(chrn ~.,switch = "y") +
    
    scale_x_continuous(
      expand = c(0,0),
      breaks = c(10, 50, 100, 200) * 1e6,
      labels = label_number(scale = 1e-6, suffix = "Mb")) +
    
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text.y = element_text(angle = 0),
          panel.background = element_blank(),
          legend.position = c(.8, .2)) +
    
    labs(y = "Chromosome") +
    #scale_fill_viridis_d(option = "C") +
    scale_fill_manual(values=c("firebrick3", "darkolivegreen3", "deepskyblue3"))
  labs(x = "Position relative to chromosome start",
       title = individual)
  
  ## Saving the plot: -----------------------------
  ggsave(output_plot, height = 5, width = 4)
  ggsave(str_replace(output_plot, ".png", ".pdf"),
         height = 5, width = 4)
  
}

#### Generate karyograms for all individuals: 
samples <- read.table("~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_XGMix_2_Tracts/3_CRELES_IDs.txt")
for (index in 5:length(samples$V1)) {
  save.karyograms(samples$V1[index])
}
```


# 6. Tracts

Downloaded from: https://github.com/santiago1234/tracts/tree/python3


```python


cd ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/6_Tracts
python3 

####### --------------------------------------------------------------------------

# python packages for data analysis
import numpy as np
import pandas as pd
import os
import sys
sys.path.append("/home/BCRICWH.LAN/Paola.Arguello/KoborLab/kobor_space/ppascualli/Programs/tracts-python3/")
import tracts
from warnings import warn

#### Loading data to Tracts ---------------------------------
path_to_tracts_input = "/home/BCRICWH.LAN/Paola.Arguello/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/4_XGMix_2_Tracts/4_Tracts_input_ready/"
labels = ['AMR', 'EUR', 'AFR']

_files = os.listdir(path_to_tracts_input)
files = [file
        for file in _files
        if file.split('.')[-1] == "bed"]  # only consider bed files

# string between individual label and haploid chromosome id in input file
inter = "_CRELES"
# string at the end of input file. Note that the file should end in ".bed"
end = "_cM.bed"

# here we extract sample name from file name ----------------------------------
names = list(set(file.split('_CRELES_')[0] for file in files))

if len(_files) != len(files):
    warn("some files in the bed directory were ignored, since they do not "
            "end with `.bed`.")
            
# Load the population using the population class's constructor. It
# automatically iterates over individuals and haploid copies (labeled _A"
# and "_B" by default
pop = tracts.population(names=names, fname=(path_to_tracts_input, inter, end))

```


## 5.1 Population global ancestry (from LAI)


```python

## Global ancestry in the population
props = list(pop.get_mean_ancestry_proportions(labels))
props = pd.DataFrame(data={'anc': labels, 'p': props})
props

```

<table class=" lightable-paper lightable-striped" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:left;"> ancestry </th>
   <th style="text-align:right;"> p </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> AMR </td>
   <td style="text-align:right;"> 0.330 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> EUR </td>
   <td style="text-align:right;"> 0.584 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AFR </td>
   <td style="text-align:right;"> 0.085 </td>
  </tr>
</tbody>
</table>


## 5.2 Per individual global ancestry (from LAI)


```python

## Global ancestry per individual
anc_p = pop.get_means(ancestries=labels)
anc_p = pd.DataFrame(anc_p, columns=labels)
anc_p['Samplename'] = names
anc_p.head()

```

<table class=" lightable-paper lightable-striped" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:right;"> EUR </th>
   <th style="text-align:right;"> AMR </th>
   <th style="text-align:right;"> AFR </th>
   <th style="text-align:left;"> Samplename </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.732 </td>
   <td style="text-align:right;"> 0.250 </td>
   <td style="text-align:right;"> 0.017 </td>
   <td style="text-align:left;"> 409_1086-R2 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.523 </td>
   <td style="text-align:right;"> 0.430 </td>
   <td style="text-align:right;"> 0.046 </td>
   <td style="text-align:left;"> 17_3381-R2 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.277 </td>
   <td style="text-align:right;"> 0.602 </td>
   <td style="text-align:right;"> 0.120 </td>
   <td style="text-align:left;"> 438_9119-R2 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.434 </td>
   <td style="text-align:right;"> 0.327 </td>
   <td style="text-align:right;"> 0.237 </td>
   <td style="text-align:left;"> 48_1629 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.731 </td>
   <td style="text-align:right;"> 0.236 </td>
   <td style="text-align:right;"> 0.031 </td>
   <td style="text-align:left;"> 238_5253 </td>
  </tr>
</tbody>
</table>



## 5.3 Distribution of tracts length 

```python
# generate the histogram of tract lengths
(bins, data) = pop.get_global_tractlengths(npts=50)
data = pd.DataFrame(data)
data['CM'] = bins
data.head()

```

<table class=" lightable-paper lightable-striped" style='font-family: "Arial Narrow", arial, helvetica, sans-serif; width: auto !important; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:right;"> EUR </th>
   <th style="text-align:right;"> AMR </th>
   <th style="text-align:right;"> AFR </th>
   <th style="text-align:right;"> CM </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 56669 </td>
   <td style="text-align:right;"> 43599 </td>
   <td style="text-align:right;"> 33147 </td>
   <td style="text-align:right;"> 0.0000 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 12627 </td>
   <td style="text-align:right;"> 13673 </td>
   <td style="text-align:right;"> 4812 </td>
   <td style="text-align:right;"> 0.0568 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 9386 </td>
   <td style="text-align:right;"> 9480 </td>
   <td style="text-align:right;"> 2895 </td>
   <td style="text-align:right;"> 0.1137 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 7558 </td>
   <td style="text-align:right;"> 6664 </td>
   <td style="text-align:right;"> 1774 </td>
   <td style="text-align:right;"> 0.1705 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5843 </td>
   <td style="text-align:right;"> 4616 </td>
   <td style="text-align:right;"> 1057 </td>
   <td style="text-align:right;"> 0.2274 </td>
  </tr>
</tbody>
</table>


## 5.4 Saving everything


```python

# save tables
data.to_csv('tract_length_distribution.csv', index=False)
props.to_csv('global-ancestry-population.csv', index=False)
anc_p.to_csv('global-ancestry-individual.csv', index=False)

```


# Comparing local vs global ancestry prediction

We will contrast this to the Local ancestry results as a sanity check.

  - EUR: IBR + GBR
  
  - NAT: PEL and MXL that were shown to be representative for NAT ancestry
  
  - AFR: YRI



```r
## In server: 
# -  ADMIXTURE: ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Global_Ancestry/ADMIXTURE/merged_801kgp-CR357K_mostInf.3.Q
  
# -  Local Ancestry: ~/KoborLab/kobor_space/ppascualli/Files/CRELES/Local_Ancestry/6_Tracts/global-ancestry-individual.csv


## Global ancestry: ------------------
global <- read.table("/Users/paola.arguello/Documents/UBC/LANGEBIO/ADMIXTURE_files/merged_801kgp-CR357K_mostInf.3.Q")
IDs <- read.table("/Users/paola.arguello/Documents/UBC/LANGEBIO/ADMIXTURE_files/merged_801kgp-CR357K_IDs.txt")
global$V4 <- IDs$V1
colnames(global) <- c("gEUR", "gAFR", "gNAT", "ID")

## Local ancestry: -------------------
local <- read.table("/Users/paola.arguello/Documents/UBC/LANGEBIO/ADMIXTURE_files/global-ancestry-individual.csv", 
                    sep = ",", header = T)

local$Samplename <- unlist(strsplit(as.character(local$Samplename), "_"))[seq(2,932,2)]
colnames(local) <- c("lNAT", "lEUR", "lAFR", "ID")

## Global vs Local: ------------------
global.local <- merge(global, local, by="ID")
#write.table(global.local, file = "~/KoborLab/kobor_space/ppascualli/Files/CRELES/Global_Ancestry/ADMIXTURE/global_local_v2.txt",
#            quote = F, col.names = T, row.names = F)

library(ggplot2)
```

```
## Warning in as.POSIXlt.POSIXct(Sys.time()): unknown timezone 'zone/tz/2021a.2.0/
## zoneinfo/America/Mexico_City'
```

```r
ggplot(global.local, aes(x=lNAT, y=gNAT)) + 
  geom_point(colour="olivedrab3") +
  xlab("Local ancestry") + ylab("Global ancestry") +
  ggtitle("Native American ancestry")
```

![](CRELES_local-ancestry-inference_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

```r
#ggsave("NativeAmerican_GvL.png")

ggplot(global.local, aes(x=lEUR, y=gEUR)) + 
  geom_point(colour="firebrick3") +
  xlab("Local ancestry") + ylab("Global ancestry") +
  ggtitle("African ancestry")
```

![](CRELES_local-ancestry-inference_files/figure-html/unnamed-chunk-25-2.png)<!-- -->

```r
#ggsave("African_GvL.png")

ggplot(global.local, aes(x=lAFR, y=gAFR)) + 
  geom_point(colour="deepskyblue3") +
  xlab("Local ancestry") + ylab("Global ancestry") +
  ggtitle("European ancestry")
```

![](CRELES_local-ancestry-inference_files/figure-html/unnamed-chunk-25-3.png)<!-- -->

```r
#ggsave("European_GvL.png")
```




