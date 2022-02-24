##Mod pipeline post review
module load samtools
module load sratoolkit
module load gcc


##########################################
##### GET SOURCE FILES ###################
##########################################
#This is from Misha, will be deposited on NCBI
wget https://www.dropbox.com/sh/vl2vwz1yc9kyc9u/AAAVRQ9RJgAJhzyXQ2IUuxtXa?dl=0
#rename as a zip file
#then unzip

#my bioproject PRJNA757238
export BioProject=PRJNA473816
/home1/06909/cbscott/edirect/esearch -db sra -query $BioProject | /home1/06909/cbscott/edirect/efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $HOME/edirect/esearch -db sra -query $BioProject | $HOME/edirect/efetch -format runinfo > $BioProject.fullMeta.csv
>gets
for A in `cat $BioProject.SRR`;do
echo "fastq-dump $A" >> gets; done
python2 $HOME/bin/launcher_creator.py -j gets -n gets -a tagmap -e cbscott@utexas.edu -t 20:00:00 -N 12 -w 4 -q skx-normal



#convert source files to fastqs
>tofastq #needs gcc
for file in *.bam; do
echo "$HOME/bam2fastq/bam2fastq $file -o ${file/.bam/}#.fastq" >> tofastq;
done
python2 $HOME/bin/launcher_creator.py -j tofastq -n tofastq -a tagmap -e cbscott@utexas.edu -t 01:00:00 -w 4 -q skx-dev
#clean it up



###### REF GENOME ##################
##Create new "humanfied" genome using A. millepora
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/753/865/GCA_013753865.1_Amil_v2.1/GCA_013753865.1_Amil_v2.1_genomic.fna.gz
#this will be our new mapped genome.
#cat together the chromosomes & rename them
#NOTE! further down the pipeline, the projection-based PCA cannot handle more than 99 chromosomes!
#For this genome, we eventually end up with 16 chromosomes, as it is a chromosome level assembly
#I end up concating the unplaced scaffolds into pseudo chromosomes
#This sample chunk will randomly concat things together... I manually did it separately

$HOME/bin/concatFasta.pl fasta=GCA_013753865.1_Amil_v2.1_genomic.fna num=16
awk '/^>/{print ">" ++i; next}{print}' <  GCA_013753865_cc.fasta > numbered_reduced_Amil.fasta
export GEN_AMIL=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/numbered_reduced_Amil.fasta
#add zoox genomes, add MAG genomes
#these are all in 'mastergen'

#this creates:
export GEN_NEW=/work/06909/cbscott/aDNA_mastergen/mastergenome/Amil_concat/Amil_masterconcat_humanfied.fasta

echo "bowtie2-build $GEN_NEW $GEN_NEW" >btb
python2 $HOME/bin/launcher_creator.py -j btb -n btb -l btbl -t 03:00:00 -a tagmap -e cbscott@utexas.edu -w 1 -q skx-normal
sbatch btbl

samtools faidx $GEN_NEW


##GET MODERN ACROPORID DATA
export BioProject=PRJNA473816
$WORK/software/edirect/esearch -db sra -query $BioProject | $WORK/software/edirect/efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $WORK/software/edirect/esearch -db sra -query $BioProject | $WORK/software/edirect/efetch -format runinfo > $BioProject.fullMeta.csv
>gets
for A in `cat $BioProject.SRR`;do
echo "fastq-dump $A" >> gets; done
python2 $HOME/bin/launcher_creator.py -j gets -n gets -a tagmap -e cbscott@utexas.edu -t 20:00:00 -N 12 -w 4 -q skx-normal


###GET MODERN MICROBE DATA

#PANAMA MICROBIOME
#PRJNA295144 - https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1758-2229.12412
export BioProject=PRJNA295144
/home1/06909/cbscott/edirect/esearch -db sra -query $BioProject | /home1/06909/cbscott/edirect/efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $HOME/edirect/esearch -db sra -query $BioProject | $HOME/edirect/efetch -format runinfo > $BioProject.fullMeta.csv
>panama
for A in `cat $BioProject.SRR`;do
echo "fastq-dump $A" >> panama; done
python2 $HOME/bin/launcher_creator.py -j panama -n panama -a tagmap -e cbscott@utexas.edu -t 10:00:00 -N 1 -w 2 -q skx-normal

#FLORIDA
#https://www.nature.com/articles/s41598-019-54855-y
#PRJNA546259
export BioProject=PRJNA546259
/home1/06909/cbscott/edirect/esearch -db sra -query $BioProject | /home1/06909/cbscott/edirect/efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $HOME/edirect/esearch -db sra -query $BioProject | $HOME/edirect/efetch -format runinfo > $BioProject.fullMeta.csv

>florida
for A in `cat $BioProject.SRR`;do
echo "fastq-dump $A " >> florida; done
python2 $HOME/bin/launcher_creator.py -j florida -n florida -a tagmap -e cbscott@utexas.edu -t 05:00:00 -N 5 -w 2 -q skx-normal

##WESTRICH FLORIDA dataset
#PRJNA299413
export BioProject=PRJNA299413
$WORK/software/edirect/esearch -db sra -query $BioProject | $WORK/software/edirect/efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $WORK/software/edirect/esearch -db sra -query $BioProject | $WORK/software/edirect/efetch -format runinfo > $BioProject.fullMeta.csv

>west
for A in `cat $BioProject.SRR`;do
echo "fastq-dump $A " >> west; done
python2 $HOME/bin/launcher_creator.py -j west -n west -a tagmap -e cbscott@utexas.edu -t 05:00:00 -N 5 -w 10 -q skx-normal


##################################################
########## Ancient Read Mapping ##################
##################################################

#without trimming
>anc_e2e
for F in *.fastq; do
echo "bowtie2 --no-unal --end-to-end -x $GEN_NEW -U $F -S ${F/.Y1.E1.L1_M.fastq}.e2e.sam && \
samtools sort -O bam -o ${F/.Y1.E1.L1_M.fastq}.e2e.sorted.bam ${F/.Y1.E1.L1_M.fastq}.e2e.sam && samtools index -c ${F/.Y1.E1.L1_M.fastq}.e2e.sorted.bam" >> anc_e2e; done
python2 $HOME/bin/launcher_creator.py  -j anc_e2e -n anc_e2e -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 4 -q skx-dev
sbatch anc_e2e.slurm

#with 1bp trimming from each end
>trimmer
for file in *.fastq; do
echo "seqtk trimfq -b 1 -e 1 $file > ${file/.fastq}.trim.fastq" >> trimmer; done
python2 $HOME/bin/launcher_creator.py  -j trimmer -n trimmer -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 4 -q skx-normal
sbatch trimmer.slurm

>anc_trim
for F in *.trim.fastq; do
echo "bowtie2 --no-unal --end-to-end -x $GEN_NEW -U $F -S ${F/.Y1.E1.L1_M.trim.fastq}.e2e.trim.sam && \
samtools sort -O bam -o ${F/.Y1.E1.L1_M.trim.fastq}.e2e.trim.sorted.bam ${F/.Y1.E1.L1_M.trim.fastq}.e2e.trim.sam && samtools index -c ${F/.Y1.E1.L1_M.trim.fastq}.e2e.trim.sorted.bam" >> anc_trim; done
python2 $HOME/bin/launcher_creator.py  -j anc_trim -n anc_trim -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 4 -q skx-normal
sbatch anc_trim.slurm

>anc_local
for F in *.trim.fastq; do
echo "bowtie2 --no-unal --local -x $GEN_NEW -U $F -S ${F/.Y1.E1.L1_M.trim.fastq}.local.trim.sam && \
samtools sort -O bam -o ${F/.Y1.E1.L1_M.trim.fastq}.local.trim.sorted.bam ${F/.Y1.E1.L1_M.trim.fastq}.local.trim.sam && samtools index -c ${F/.Y1.E1.L1_M.trim.fastq}.local.trim.sorted.bam" >> anc_local; done
python2 $HOME/bin/launcher_creator.py  -j anc_local -n anc_local -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 4 -q skx-normal
sbatch anc_local.slurm
##NOW Filter all for mapping quality > 30
>filter_bams
for file in *.sorted.bam;
do echo "samtools view -bSq 30 $file > ${file/.sorted.bam}.filter30.sorted.bam && samtools index -c ${file/.sorted.bam}.filter30.sorted.bam" >> filter_bams; done
python2 $HOME/bin/launcher_creator.py  -j filter_bams -n filter_bams -a tagmap -e cbscott@utexas.edu -t 00:20:00 -N 1 -w 8 -q skx-dev
sbatch filter_bams.slurm

#Create table of mapping by scaffold
for file in *trim.filter30.sorted.bam;
do samtools idxstats $file > ${file}.scaffoldstats.txt; done


##################################################
########## Modern Read Mapping ##################
##################################################
#Also do end to end, even tho map rate was low
module load samtools
>mod_e2e
for F in *.fastq; do
echo "bowtie2 --no-unal --end-to-end -x $GEN_NEW -U $F -S ${F/.fastq}.e2e.sam && \
samtools sort -O bam -o ${F/.fastq}.e2e.sorted.bam ${F/.fastq}.e2e.sam && samtools index -c ${F/.fastq}.e2e.sorted.bam" >> mod_e2e; done
python2 $HOME/bin/launcher_creator.py  -j mod_e2e -n mod_e2e -a tagmap -e cbscott@utexas.edu -t 10:00:00 -N 2 -w 10 -q skx-normal
sbatch mod_e2e.slurm #again almost no local mapping?

>mod_local
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $GEN_NEW -U $F -S ${F/.fastq}.local.sam && \
samtools sort -O bam -o ${F/.fastq}.local.sorted.bam ${F/.fastq}.local.sam && samtools index -c ${F/.fastq}.local.sorted.bam" >> mod_local; done
python2 $HOME/bin/launcher_creator.py  -j mod_local -n mod_local -a tagmap -e cbscott@utexas.edu -t 15:00:00 -N 2 -w 12 -q skx-normal
sbatch mod_local.slurm

>mod_local
for F in `cat todofastqs`; do
echo "bowtie2 --no-unal --local -x $GEN_NEW -U $F -S ${F/.fastq}.local.sam && \
samtools sort -O bam -o ${F/.fastq}.local.sorted.bam ${F/.fastq}.local.sam && samtools index -c ${F/.fastq}.local.sorted.bam" >> mod_local; done
python2 $HOME/bin/launcher_creator.py  -j mod_local -n mod_local -a tagmap -e cbscott@utexas.edu -t 10:00:00 -N 1 -w 12 -q skx-normal
sbatch mod_local.slurm

##################################################
########## Microbe Read Mapping ##################
##################################################
#Also do end to end, even tho map rate was low
>microbe_e2e
for F in *.fastq; do
echo "bowtie2 --no-unal --end-to-end -x $GEN_NEW -U $F -S ${F/.fastq}.e2e.sam && \
samtools sort -O bam -o ${F/.fastq}.e2e.sorted.bam ${F/.fastq}.e2e.sam && samtools index -c ${F/.fastq}.e2e.sorted.bam" >> microbe_e2e; done
python2 $HOME/bin/launcher_creator.py  -j microbe_e2e -n microbe_e2e -a tagmap -e cbscott@utexas.edu -t 08:00:00 -N 2 -w 12 -q skx-normal
sbatch microbe_e2e.slurm

>microbe_local
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $GEN_NEW -U $F -S ${F/.fastq}.local.sam && \
samtools sort -O bam -o ${F/.fastq}.local.sorted.bam ${F/.fastq}.local.sam && samtools index -c ${F/.fastq}.local.sorted.bam" >> microbe_local; done
python2 $HOME/bin/launcher_creator.py  -j microbe_local -n microbe_local -a tagmap -e cbscott@utexas.edu -t 08:00:00 -N 2 -w 12 -q skx-normal
sbatch microbe_local.slurm

>microbe_localquick
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $GEN_NEW -U $F -S ${F/.fastq}.local1.sam && \
samtools sort -O bam -o ${F/.fastq}.local1.sorted.bam ${F/.fastq}.local1.sam && samtools index -c ${F/.fastq}.local1.sorted.bam" >> microbe_localquick; done
python2 $HOME/bin/launcher_creator.py  -j microbe_localquick -n microbe_localquick -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 24 -q skx-dev
sbatch microbe_localquick.slurm


>filter_bams
for file in *.local1.sorted.bam;
do echo "samtools view -bSq 20 $file > ${file/.sorted.bam}.filter20.sorted.bam && samtools index -c ${file/.sorted.bam}.filter20.sorted.bam" >> filter_bams; done
python2 $HOME/bin/launcher_creator.py  -j filter_bams -n filter_bams -a tagmap -e cbscott@utexas.edu -t 00:20:00 -N 1 -w 8 -q skx-dev
sbatch filter_bams.slurm

#Create table of mapping by scaffold
for file in *filter20.sorted.bam;
do samtools idxstats $file > ${file}.scaffoldstats.txt; done


#######################################################
########### VERIFY aDNA DAMAGE ########################
#######################################################

#1. Separate out coral reads
#note the coral scaffolds file is just the numbers 1-27
>coral_sep
for file in *e2e.filter30.sorted.bam; do
echo "cat coralscaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.sorted.bam && samtools index -c ${file/.sorted.bam}.coral.sorted.bam && mv ${file/.sorted.bam}.coral.sorted.bam* ../coral_bams" >> coral_sep; done
python2 $HOME/bin/launcher_creator.py  -j coral_sep -n coral_sep -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 4 -q skx-dev
sbatch coral_sep.slurm

>coral_sep
for file in *e2e.trim.filter30.sorted.bam; do
echo "cat coralscaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.sorted.bam && samtools index -c ${file/.sorted.bam}.coral.sorted.bam && mv ${file/.sorted.bam}.coral.sorted.bam* ../coral_bams" >> coral_sep; done
python2 $HOME/bin/launcher_creator.py  -j coral_sep -n coral_sep -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 4 -q skx-dev
sbatch coral_sep.slurm

###2. MapDamage for each read
>mapdam_coral
for file in ../coral_bams/*.e2e*filter30*sorted.bam;
do echo "mapDamage -i $file -r $GEN_NEW --no-stats --merge-libraries" >> mapdam_coral; done
python2 $HOME/bin/launcher_creator.py  -j mapdam_coral -n mapdam_coral -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 4 -q skx-dev


###ALSO DO THIS FOR moderns bams  - in wd "modern_coral_fastas"
#filter by quality for better plots?
>filter_bams
for file in `cat example_modern_bams`; do
echo "samtools view -bSq 20 $file > ${file/.sorted.bam}.filter20.sorted.bam && samtools index -c ${file/.sorted.bam}.filter20.sorted.bam" >> filter_bams; done
python2 $HOME/bin/launcher_creator.py  -j filter_bams -n filter_bams -a tagmap -e cbscott@utexas.edu -t 00:20:00 -N 1 -w 8 -q skx-dev
sbatch filter_bams.slurm

>mapdam_coral
for file in `cat example_modern_filtered`;
do echo "mapDamage -i $file -r $GEN_NEW --no-stats --merge-libraries" >> mapdam_coral; done
python2 $HOME/bin/launcher_creator.py  -j mapdam_coral -n mapdam_coral -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 4 -q skx-dev



####FOR SOME REASON HAVE WAY MORE ZOOX READS - also must do this, exp for S17466
#note the coral scaffolds file is just the numbers 1-27
>zoox_sep
for zoox in `cat zoox_scaffolds`; do
  for file in *e2e.filter30.sorted.bam; do
    echo "samtools view -bh $file $zoox > ${file/.sorted.bam}.${zoox}.sorted.bam && samtools index -c ${file/.sorted.bam}.${zoox}.sorted.bam" >> zoox_sep; done; done
python2 $HOME/bin/launcher_creator.py  -j zoox_sep -n zoox_sep -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 4 -q skx-dev
sbatch zoox_sep.slurm

>mapdam_zoox
for file in ../zoox_bams/*.e2e*filter30*sorted.bam;
do echo "mapDamage -i $file -r $GEN_NEW --no-stats --merge-libraries" >> mapdam_zoox; done
python2 $HOME/bin/launcher_creator.py  -j mapdam_zoox -n mapdam_zoox -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 4 -q skx-dev


###############################
#SEGREGATE MODERN FILES #######
###############################

>coral_sep
for file in *.sorted.bam; do
echo "cat coralscaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.sorted.bam && samtools index -c ${file/.sorted.bam}.coral.sorted.bam" >> coral_sep; done
python2 $HOME/bin/launcher_creator.py  -j coral_sep -n coral_sep -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 2 -w 16 -q skx-normal
sbatch coral_sep.slurm


python2 $HOME/bin/launcher_creator.py  -j indexer -n indexer -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 16 -q skx-dev



#####################################################
####### DETERMINE PUTATIVELY ANCIENT MICROBES #######
#####################################################

#noncoral/zoox scaffolds
#just make this list by hand
microbe_scaffolds

>microbe_sep
for file in *.e2e.filter30.sorted.bam; do
  for scaff in `cat microbe_scaffolds`; do
 echo "samtools view -bh $file $scaff > ${file/.sorted.bam}.${scaff}.sorted.bam && samtools index -c ${file/.sorted.bam}.${scaff}.sorted.bam && mv ${file/.sorted.bam}.${scaff}.sorted.bam microbe_bams && mv ${file/.sorted.bam}.${scaff}.sorted.bam.csi microbe_bams" >> microbe_sep; done; done
python2 $HOME/bin/launcher_creator.py  -j microbe_sep -n microbe_sep -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 16 -q skx-dev
sbatch microbe_sep.slurm

#run map damage on each assembly
module load intel
module load Rstats
>mapdam_microbes
for file in *.sorted.bam;
do echo "mapDamage -i $file -r $GEN_NEW --no-stats --merge-libraries" >> mapdam_microbes; done
python2 $HOME/bin/launcher_creator.py  -j mapdam_microbes -n mapdam_microbes -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 16 -q skx-dev

###Now use R script
#I recommend writing your own script with the algorithmic details in the paper
#I'm not sure I trust this script, and ended up manually doublechecking most of the files
#anyways
putatively_ancient_microbes.R #(it was local in the AncientDNA folder)




#########################################################
##########   ANGSD PIPELINE #############################
#########################################################

#This time use trimmed e2e files - rather than local mappings
export anc_bam=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/angsd_pipe/ancient_bams
export mod_bam=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/angsd_pipe/modern_bams

angsd -out anc_count -doCounts 1 -dumpCounts 2 -bam $anc_bam #this is pretty speedy

gunzip *.gz
#what we want from regions is stored in the anc_count.pos object


#create an angsd file from all data
cat $anc_bam $mod_bam > all_bams
echo "angsd -bam all_bams -dohaplocall 1 -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out all_sites -nThreads 8" > all_sites
python2 $HOME/bin/launcher_creator.py -j all_sites -n all_sites -a tagmap -e cbscott@utexas.edu -t 10:00:00 -N 1 -w 1 -q skx-normal

#parallelize by scaffold... and let it be threaded
>parallel_short
for i in {1..27}; do
echo "angsd -bam all_bams -r ${i}: -dohaplocall 1 -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out par_short_${i} -nThreads 4" >> parallel_short; done
python2 $HOME/bin/launcher_creator.py -j parallel_angsd -n parallel_angsd -a tagmap -e cbscott@utexas.edu -t 10:00:00 -N 12 -w 1 -q skx-normal

>parallel_redo
for i in 2; do
echo "angsd -bam all_bams -r ${i}: -dohaplocall 1 -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out par_redo_${i} -nThreads 8" >> parallel_redo; done
python2 $HOME/bin/launcher_creator.py -j parallel_redo -n parallel_redo -a tagmap -e cbscott@utexas.edu -t 02:30:00 -N 1 -w 1 -q skx-normal

#we will restrict the resulting files for sites covered in at least one ancient file

echo "angsd -bam all_bams -dohaplocall 1 -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out all_sites_long -nThreads 24" > optimize_mem
python2 $HOME/bin/launcher_creator.py -j optimize_mem -n optimize_mem -a tagmap -e cbscott@utexas.edu -t 20:00:00 -N 1 -w 1 -q skx-normal


for file in par*.beagle.gz;
do gunzip $file; done

#to cat them together, keep one with header
par_redo_1.beagle #keep header
#remove the header from the remaining files
for i in {1..27}; do
sed -i '1d' par_${i}.beagle; done #works

#cat them together
>properorder
for i in {2..27}; do
echo "par_${i}.beagle" >> properorder
done

cat properorder | tr '\n' ' ' | xargs cat > noheader.beagle
cat par_1.beagle noheader.beagle > catted_sites.beagle


#create my region of interest file
awk '{print $1}' anc_count.pos | sed '1d' > anc_chrs
awk '{print $2}' anc_count.pos | sed '1d' > anc_pos
paste -d'_' anc_chrs anc_pos > anc_covered.regions


#subset our beagle file for just our regions of interst
awk 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' anc_covered.regions catted_sites.beagle > subset_sites.beagle
#only 9 million sites in catted_beagle
#only 1145 sites? regions file might be too stringent

gzip subset_sites.beagle #angsd needs a gzipped file
module load gcc
python3 /home1/06909/cbscott/pcangsd/pcangsd.py -beagle subset_sites.beagle.gz

#pull locally and plot?
#YES! You have an rscript locally called pcangsd.R

####NOW ALSO NEED TO DEAL WITH THE ADMIXTURE ANALYSIS
for i in {2..7}; do
NGSadmix -likes subset_sites.beagle.gz -K $i -o amil_admixture.${i}; done
#this is a speedy little function




#########################################################
##########   PROJECTION PIPELINE  #######################
#########################################################
module load samtools
export GEN_NEW=/work/06909/cbscott/aDNA_mastergen/mastergenome/Amil_concat/Amil_masterconcat_humanfied.fasta
export modernbams=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/modern_coral_fastas/modern_bams

>parallelized_calls
for scaff in `cat coralscaffolds`; do
echo "bcftools mpileup -Ou -f $GEN_NEW -r $scaff -b $modernbams | bcftools call -m -v -Oz -o ${scaff}.vcf.gz" >> parallelized_calls; done
python2 $HOME/bin/launcher_creator.py  -j parallelized_calls -n parallelized_calls -a tagmap -e cbscott@utexas.edu -t 20:00:00 -N 2 -w 16 -q skx-normal

>properorder
for i in {1..27}; do
echo "${i}.vcf.gz" >> properorder
done

cat properorder | tr '\n' ' ' | xargs bcftools concat > modern_mil.vcf.gz
#get just the biallelic snps (is this correct?)
bcftools view -m2 -M2 -v snps modern_mil.vcf.gz | sed '/^#/d'| cut -f1,2,4,5 > modern_all_reqcols.txt
wc -l modern_all_reqcols.txt # 10966871

#this is much faster than the old for loop
seq 1 10966871 > snpnums
awk '{ printf "mod_mil"; print }' snpnums > snpids
awk 'BEGIN{while(++i<=10966871){print "0"}}' > genpos

cut -f1 modern_all_reqcols.txt  > snpchr
cut -f2,3,4 modern_all_reqcols.txt > snpposrefalt
paste snpids snpchr genpos snpposrefalt > modern.mil.snp

export GEN_NEW=/work/06909/cbscott/aDNA_mastergen/mastergenome/Amil_concat/Amil_masterconcat_humanfied.fasta
export modernbams=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/modern_coral_fastas/modern_bams
export SNP=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/modern_coral_fastas/modern.mil.snp
export anc_bam=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/angsd_pipe/ancient_bams
export anc_bam=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/angsd_pipe/wo_66
export anc_bam=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/angsd_pipe/66

cut -d'/' -f9 $modernbams | cut -d'.' -f1 > modern_names
cut -d'/' -f9 $anc_bam | cut -d'.' -f1 > anc_names

echo "samtools mpileup -R -B -q10 -Q0 -f $GEN_NEW -b $modernbams | pileupCaller --randomHaploid --sampleNames `cat modern_names | tr '\n' ',' | sed 's/.$//'` -f $SNP -e modern" > modtest #don't use -e modern option.... what does this create?
python2 $HOME/bin/launcher_creator.py -j modtest -n modtest -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 1 -q skx-dev

echo "samtools mpileup -R -B -q10 -Q0 -f $GEN_NEW -b $anc_bam | pileupCaller --randomHaploid --sampleNames `cat anc_names | tr '\n' ',' | sed 's/.$//'` -f $SNP -e ancient_66" > anctest #don't use -e modern option.... what does this create?
python2 $HOME/bin/launcher_creator.py -j anctest -n anctest -a tagmap -e cbscott@utexas.edu -t 00:15:00 -N 1 -w 1 -q skx-dev


#dependency issues? what's going on... do it tomorrow
#ISSUES WITH MERGE IT! MEMEROY CACHE? INSTEAD DO BY HAND
echo "mergeit -p parMergeIt" > mergeit
python2 $HOME/bin/launcher_creator.py -j mergeit -n mergeit -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 1 -q skx-dev


export anc_pre=ancient_S17463 #filtering for just one ancinet sample at at time worked... I suspect S17466 is the issue
export mod_pre=modern
export out_pre=merge_mil

sed -i 's/Unknown/Modern/g' ${mod_pre}.ind.txt
sed -i 's/Unknown/Ancient/g' ${anc_pre}.ind.txt

echo "Modern" > poplist.txt

paste -d'\0' ${anc_pre}.geno.txt ${mod_pre}.geno.txt > ${out_pre}.geno
cat ${anc_pre}.ind.txt ${mod_pre}.ind.txt > ${out_pre}.ind
cp ${mod_pre}.snp.txt ${out_pre}.snp

###NOW CREATE PCA! We will need to project, and I know there's something up with the pop names here
#I think poplist will be the populations you want to use to construct the PCA
#In this case, let's leave it "Unknown"
echo "smartpca -p parSmartPCA" > pca
python2 $HOME/bin/launcher_creator.py -j pca -n pca -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 1 -q skx-dev

sed 1d ${out_pre}.evec | tr -s ' ' | sed -e 's/^[ \t]*//' | awk 'NF{NF-=1};1' | awk -v OFS="\t" '$1=$1' | sed -e '1s/^/Sample_ID\tpc1\tpc2\tpc3\tpc4\tpc5\tpc6\tpc7\tpc8\tpc9\tpc10\n/' > ${out_pre}.evec.txt


###DO PROJECTION
#CREATED R CODE TO ASSIST
#PARALLELIZE
split -l 100000 ancient_S17463.snp.txt snpsub
>parallel_names
for file in `cat subsets`; do
echo "Rscript chr_namer_in.R translate_scaffolds.txt $file ${file}_indexed" >> parallel_names; done
python2 $HOME/bin/launcher_creator.py -j parallel_names -n parallel_names -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 24 -q skx-dev

for line in `cat subsets`;
do echo ${line}_indexed >> order; done

cat order | tr '\n' ' ' | xargs cat > indexed.snp.txt

#now create missing ind files, using new R script: cmd_JackScaffold.R

echo "Rscript cmd_JackScaffold.R indexed.snp.txt modern.geno.txt modern.ind.txt" > mod_missing
python2 $HOME/bin/launcher_creator.py -j mod_missing -n mod_missing -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 1 -q skx-dev

paste -d'\0' ancient.missing.geno modern.missing.geno > merged.missing.geno
cat ancient.missing.ind modern.missing.ind > merged.missing.ind
cp ancient_S17463.snp.txt merged.missing.snp

echo "smartpca -p parSmartPCAMissing" > pca
python2 $HOME/bin/launcher_creator.py -j pca -n pca -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 1 -q skx-dev

sed 1d merged.missing_mill.evec | tr -s ' ' | sed -e 's/^[ \t]*//' | awk 'NF{NF-=1};1' | awk -v OFS="\t" '$1=$1' | sed -e '1s/^/Sample_ID\tpc1\tpc2\tpc3\tpc4\tpc5\tpc6\tpc7\tpc8\tpc9\tpc10\n/' > merged.missing_mill.evec.txt



#########################################################
########   f4STATISTICS - FOR HIGHEST COV SAMPLE  #######
#########################################################
module load sratoolkit

#ah, will need to map some amil reads
export BioProject=PRJNA593014 #these are amil reads

$WORK/software/edirect/esearch -db sra -query $BioProject | $WORK/software/edirect/efetch -format runinfo |cut -d "," -f 1 | grep SRR > $BioProject.SRR && $WORK/software/edirect/esearch -db sra -query $BioProject | $WORK/software/edirect/efetch -format runinfo > $BioProject.fullMeta.csv

#let's just grab the first 15 individuals - don't need the full pop
head -15 $BioProject.SRR > head_${BioProject}.SRR
for A in `cat head_$BioProject.SRR`;do
echo "fastq-dump $A" >> grabber; done
python2 $HOME/bin/launcher_creator.py -j grabber -n grabber -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 8 -q skx-dev

module load samtools
export GEN_NEW=/work/06909/cbscott/aDNA_mastergen/mastergenome/Amil_concat/Amil_masterconcat_humanfied.fasta

>amil_map_local
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $GEN_NEW -U $F -S ${F/.fastq}.local.sam && \
samtools sort -O bam -o ${F/.fastq}.local.sorted.bam ${F/.fastq}.local.sam && samtools index -c ${F/.fastq}.local.sorted.bam" >> amil_map_local; done
python2 $HOME/bin/launcher_creator.py  -j amil_map_local -n amil_map_local -a tagmap -e cbscott@utexas.edu -t 04:00:00 -N 1 -w 4 -q skx-normal
sbatch amil_map_local.slurm


#sort out amils by coral only - don't filter for quality
>coral_sep
for file in *local.sorted.bam; do
echo "cat coral_scaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.sorted.bam && samtools index -c ${file/.sorted.bam}.coral.sorted.bam" >> coral_sep; done
python2 $HOME/bin/launcher_creator.py  -j coral_sep -n coral_sep -a tagmap -e cbscott@utexas.edu -t 01:30:00 -N 1 -w 12 -q skx-dev
sbatch coral_sep.slurm

#create a new list of all modern bams, including amils, then do vcf creation
module load samtools
export GEN_NEW=/work/06909/cbscott/aDNA_mastergen/mastergenome/Amil_concat/Amil_masterconcat_humanfied.fasta
#export modernbams=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/f4stats/modernbams
export modernbams=/scratch/06909/cbscott/aDNA_revamp/fstat/Amil_files/fastq/modern_amil_bams
#export modernbams=/scratch/06909/cbscott/aDNA_revamp/fstat/genoall/modernacrobams

>parallelized_calls
for scaff in `cat coralscaffolds`; do
echo "bcftools mpileup -Ou -f $GEN_NEW -r $scaff -b $modernbams | bcftools call -m -v -Oz -o ${scaff}.vcf.gz" >> parallelized_calls; done
python2 $HOME/bin/launcher_creator.py  -j parallelized_calls -n parallelized_calls -a tagmap -e cbscott@utexas.edu -t 20:00:00 -N 2 -w 16 -q skx-normal


#cat resulting files together
gunzip *.gz
cat properorder | tr '\n' ' ' | xargs cat > amil.vcf
gzip amil.vcf

#do pseudohaplocalls fora ll of the modern (a palm and a mil) individuals
#weird errors?
bcftools view -m2 -M2 -v snps amil.vcf.gz | sed '/^#/d'| cut -f1,2,4,5 > modern_all_reqcols.txt
wc -l modern_all_reqcols.txt #10809704

#this is much faster than the old for loop
seq 1 10809704 > snpnums
awk '{ printf "mod_mil"; print }' snpnums > snpids
awk 'BEGIN{while(++i<=10809704){print "0"}}' > genpos

cut -f1 modern_all_reqcols.txt  > snpchr
cut -f2,3,4 modern_all_reqcols.txt > snpposrefalt
paste snpids snpchr genpos snpposrefalt > amil.snp

export GEN_NEW=/work/06909/cbscott/aDNA_mastergen/mastergenome/Amil_concat/Amil_masterconcat_humanfied.fasta
export modernbams=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/f4stats/modernbams
export SNP=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/f4stats/amil.snp
export anc_bam=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/angsd_pipe/ancient_bams
#export anc_bam=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/angsd_pipe/wo_66
#export anc_bam=/scratch/06909/cbscott/aDNA_revamp/reviews/amil_mapping/angsd_pipe/66

cut -d'/' -f9 $modernbams | cut -d'.' -f1 > modern_names
cut -d'/' -f9 $anc_bam | cut -d'.' -f1 > anc_names

echo "samtools mpileup -R -B -q10 -Q0 -f $GEN_NEW -b $modernbams | pileupCaller --randomHaploid --sampleNames `cat modern_names | tr '\n' ',' | sed 's/.$//'` -f $SNP -e modern" > modtest #don't use -e modern option.... what does this create?
python2 $HOME/bin/launcher_creator.py -j modtest -n modtest -a tagmap -e cbscott@utexas.edu -t 03:00:00 -N 1 -w 1 -q skx-normal

echo "samtools mpileup -R -B -q10 -Q0 -f $GEN_NEW -b $anc_bam | pileupCaller --randomHaploid --sampleNames `cat anc_names | tr '\n' ',' | sed 's/.$//'` -f $SNP -e ancient" > anctest #don't use -e modern option.... what does this create?
python2 $HOME/bin/launcher_creator.py -j anctest -n anctest -a tagmap -e cbscott@utexas.edu -t 00:15:00 -N 1 -w 1 -q skx-dev


#merge results after renaming populations by hand/in R

paste -d'\0' ancient.geno.txt modern.geno.txt > f4all.geno
cat ancient.ind.txt namedmodern.ind.txt > f4all.ind
cp modern.snp.txt f4all.snp


#Try using different arg
genotypename:   f4all.geno
snpname:        f4all.snp
indivname:      f4all.ind
poplistname:    poplist.txt
f4mode: YES

echo "qpDstat -p f4parms" > f4job
python2 $HOME/bin/launcher_creator.py -j f4job_specpops -n f4job_pops -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 1 -q skx-dev
