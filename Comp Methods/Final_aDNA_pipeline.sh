#Pipeline From the top
##########################################
##### SOFTWARE NEEDED ####################
##########################################
bam2fastq
bowtie2
conda
mapDamage
parallel
EIGENSOFT

####################################
######  INSTALLATIONS ##############
####################################
#I used conda to manage most of my installations -
# on our computing system these softwares were
#particularly difficult to wrangle
conda install parallel
conda install SequenceTools
conda install openblas
make * CFLAGS="-I/work2/06909/cbscott/stampede2/software/anaconda3/pkgs/openblas-0.3.3-h9ac9557_1001/include -I/work2/06909/cbscott/stampede2/software/anaconda3/pkgs/gsl-2.6-hf94e986_0/include" LDFLAGS="-L/work2/06909/cbscott/stampede2/software/anaconda3/pkgs/openblas-0.3.3-h9ac9557_1001/lib -L/work2/06909/cbscott/stampede2/software/anaconda3/pkgs/gsl-2.6-hf94e986_0/lib"
##########################################
##### LOAD MODULES    ####################
##########################################
module load samtools
module load python2 #launcher creator does not work with python3
module load gcc
module load Rstats
module load intel/18.0.2
#module load gsl/2.3
module load gsl
module load python3/3.7.0


##########################################
##### GET SOURCE FILES ###################
##########################################
#This is from Misha, will be deposited on NCBI
wget https://www.dropbox.com/sh/vl2vwz1yc9kyc9u/AAAVRQ9RJgAJhzyXQ2IUuxtXa?dl=0
#rename as a zip file
#then unzip

#convert source files to fastqs
>tofastq
for file in *.bam; do
echo "$HOME/bam2fastq/bam2fastq $file -o ${file/.bam/}#.fastq" >> tofastq;
done
python $HOME/bin/launcher_creator.py -j tofastq -n tofastq -a tagmap -e cbscott@utexas.edu -t 01:00:00 -w 4 -q development
#clean it up

##########################################
########## MAP READS  ####################
##########################################
export MASTERGEN=$STOCKYARD/aDNA_mastergen/mastergenome/human-fied_genome/concat_master_genome.fasta

#need default mapping and local mapping (want to make sure we trim damage in our ancient reads)
>anc_e2e
for F in *.fastq; do
echo "bowtie2 --no-unal --end-to-end -x $MASTERGEN -U $F -S ${F/.Y1.E1.L1_M.fastq}.e2e.sam && \
samtools sort -O bam -o ${F/.Y1.E1.L1_M.fastq}.e2e.sorted.bam ${F/.Y1.E1.L1_M.fastq}.e2e.sam && samtools index -c ${F/.Y1.E1.L1_M.fastq}.e2e.sorted.bam" >> anc_e2e; done
python $HOME/bin/launcher_creator.py  -j anc_e2e -n anc_e2e -a tagmap -e cbscott@utexas.edu -t 10:00:00 -N 1 -w 4 -q normal
sbatch anc_e2e.slurm

>anc_local
for F in *.fastq; do
echo "bowtie2 --no-unal --local -x $MASTERGEN -U $F -S ${F/.Y1.E1.L1_M.fastq}.local.sam && \
samtools sort -O bam -o ${F/.Y1.E1.L1_M.fastq}.local.sorted.bam ${F/.Y1.E1.L1_M.fastq}.local.sam && samtools index -c ${F/.Y1.E1.L1_M.fastq}.local.sorted.bam" >> anc_local; done
python $HOME/bin/launcher_creator.py  -j anc_local -n anc_local -a tagmap -e cbscott@utexas.edu -t 10:00:00 -N 1 -w 4 -q normal
sbatch anc_local.slurm

#let's move all of these bams to a new DIRECTORY
#eventually we'll want them in work2
mkdir sorted_bams
cd sorted_bams
#next filter for quality >20 in each file
>filter_bams
for file in *.sorted.bam;
do echo "samtools view -bSq 30 $file > ${file/.sorted.bam}.filter30.sorted.bam && samtools index -c ${file/.sorted.bam}.filter30.sorted.bam" >> filter_bams; done
python $HOME/bin/launcher_creator.py  -j filter_bams -n filter_bams -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 8 -q development
sbatch filter_bams.slurm

#get total mapped read counts
#Local mapping is retaining many many more reads than end-to-end mapping
#Separate bams by species/coral/zooxs/MAGS
mkdir coral_bams
mkdir zoox_bams
mkdir MAG_bams

>coral_sep
for file in *.filter30.sorted.bam; do
echo "cat NumberedCoralScaffolds | tr '\n' ' ' | xargs samtools view -bh $file > ${file/.sorted.bam}.coral.sorted.bam && samtools index -c ${file/.sorted.bam}.coral.sorted.bam && mv ${file/.sorted.bam}.coral.sorted.bam* ../coral_bams" >> coral_sep; done
python $HOME/bin/launcher_creator.py  -j coral_sep -n coral_sep -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 8 -q development
sbatch coral_sep.slurm

>zoox_sep
for zoox in `cat zoox_scaffolds`; do
  for file in *.filter30.sorted.bam; do
    echo "samtools view -bh $file $zoox > ${file/.sorted.bam}.${zoox}.sorted.bam && samtools index -c ${file/.sorted.bam}.${zoox}.sorted.bam && mv ${file/.sorted.bam}.${zoox}.sorted.bam ../zoox_bams && mv ${file/.sorted.bam}.${zoox}.sorted.bam.csi ../zoox_bams" >> zoox_sep; done; done
python $HOME/bin/launcher_creator.py  -j zoox_sep -n zoox_sep -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 8 -q normal
sbatch zoox_sep.slurm


>MAG_sep
for MAG in `cat final_MAG_list`; do
  for file in *.filter30.sorted.bam; do
    echo "samtools view -bh $file $MAG > ${file/.sorted.bam}.${MAG}.sorted.bam && samtools index -c ${file/.sorted.bam}.${MAG}.sorted.bam && mv ${file/.sorted.bam}.${MAG}.sorted.bam ../MAG_bams && mv ${file/.sorted.bam}.${MAG}.sorted.bam.csi ../MAG_bams" >> MAG_sep; done; done
python $HOME/bin/launcher_creator.py  -j MAG_sep -n MAG_sep -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 2 -w 8 -q normal
sbatch MAG_sep.slurm

##### COUNT READS IN EACH FILE ######
>filelist
>readcount
for file in *filter30*.sorted.bam;
do echo ${file/.filter20.coral.sorted.bam} >> filelist && samtools view -c $file >> readcount; done
paste filelist readcount > coral_bams_readcounts.txt

>filelist
>readcount
for file in *filter30*.sorted.bam;
do echo ${file/.sorted.bam} >> filelist && samtools view -c $file >> readcount; done
paste filelist readcount > zoox_bams_readcounts.txt


>filelist
>readcount
for file in *.sorted.bam;
do echo ${file/.sorted.bam} >> filelist && samtools view -c $file >> readcount; done
paste filelist readcount > MAG_bams_readcounts.txt

################################################
############ MAP DAMAGE IN E2E READS ###########
################################################
#coral & zoox
export MASTERGEN=$STOCKYARD/aDNA_mastergen/mastergenome/human-fied_genome/concat_master_genome.fasta
>mapdam_zoox
for file in *.e2e*filter30*sorted.bam;
do echo "mapDamage -i $file -r $MASTERGEN --no-stats --merge-libraries" >> mapdam_zoox; done
python2 $HOME/bin/launcher_creator.py  -j mapdam_zoox -n mapdam_zoox -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 4 -q normal


#####ALGORITHMICALLY SORT MAGS FOR DAMAGE #####
#Run mapDamage on all end to end files
>mapdam_MAG
for file in *.e2e*filter30*.sorted.bam;
do echo "mapDamage -i $file -r $MASTERGEN --no-stats --merge-libraries" >> mapdam_MAG; done
python2 $HOME/bin/launcher_creator.py  -j mapdam_MAG -n mapdam_MAG -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 2 -w 8 -q normal

#avg readlength in local mapping reads
samtools view S17465.e2e.filter20.coral.sorted.bam | head -100 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | head

##########################################
########## METAGENOMICS ##################
##########################################

#move all map damage files locally and pass through PlotDamage.R script to sort anc versus modern mags, create nice plots, etcs
#move resulting file of likely_ancient_MAGS back to TACC
#do for local and e2e
while read -r first second; do
    mv ${first}.local.filter30.${second}.sorted.bam ancient_bams && mv ${first}.local.filter30.${second}.sorted.bam.csi ancient_bams
done < MAG_presence

#get read counts per endo in a table
>filelist
>countlist
for file in *.local.*sorted.bam; do
  echo $file >> filelist && samtools view -c $file >> countlist; done

paste filelist countlist > ancMAG_counts

#now use unique mag list files to subset modern files and do counting!
>MAG_sep
for file in *.sorted.bam; do
echo "cat $SCRATCH/aDNA_revamp/sorted_bams/scaffold_lists/all_anc_MAGs | tr '\n' ' ' | xargs samtools view -bh $file > $SCRATCH/aDNA_revamp/sorted_bams/modern_MAG_bams/${file/.master.sorted.bam}.MAG.sorted.bam && samtools index -c $SCRATCH/aDNA_revamp/sorted_bams/modern_MAG_bams/${file/.master.sorted.bam}.MAG.sorted.bam" >> MAG_sep; done
python $HOME/bin/launcher_creator.py  -j MAG_sep -n MAG_sep -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 2 -w 8 -q normal
sbatch MAG_sep.slurm

>filelist
>countlist
>maglist
for file in *.sorted.bam; do
  echo $file >> filelist && samtools view -c $file >> countlist && samtools view $file | awk '{print $3}' | head -1 >> maglist; done

paste filelist maglist countlist > modMAG_counts
#something is wrong here.... I'm using the wrong metagenome files! needed mapped metagenome files from the project....

##########################################
########## MODERN COMPARISON #############
##########################################

#filter for coral reads
>coral_sep
for file in *.sorted.bam; do
echo "cat $SCRATCH/aDNA_revamp/sorted_bams/scaffold_lists/NumberedCoralScaffolds | tr '\n' ' ' | xargs samtools view -bh $file > $SCRATCH/aDNA_revamp/sorted_bams/modern_coral_bams/${file/.master.sorted.bam}.coral.sorted.bam && samtools index -c $SCRATCH/aDNA_revamp/sorted_bams/modern_coral_bams/${file/.master.sorted.bam}.coral.sorted.bam" >> coral_sep; done
python $HOME/bin/launcher_creator.py  -j coral_sep -n coral_sep -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 2 -w 8 -q normal
sbatch coral_sep.slurm


##########################################
#PIPELINE A - as in Vagheesh et al. 2019 #
##########################################

#DONT filter for quality due to trimming issue from local mapping.
export MASTERGEN=$STOCKYARD/aDNA_mastergen/mastergenome/human-fied_genome/concat_master_genome.fasta
export modernbams=/scratch/06909/cbscott/aDNA_revamp/sorted_bams/modern_coral_bams/modernbams

#modern bams is a list of the names of all of your modern bam files
#parallelized_calls simply runs bcftools for several sets of scaffolds in parallel. Here, each "*.scaff" file has the name of three scaffolds
#This could probably be more efficiently coded on stampede's node system, but oh well.
#Currently, this takes a very long time, and consider parallelizing it further.

>parallelized_calls_normal
echo "parallel -a 1.scaff 'bcftools mpileup -Ou -f $MASTERGEN -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.normal.vcf.gz'" > parallelized_calls_normal
echo "parallel -a 2.scaff 'bcftools mpileup -Ou -f $MASTERGEN -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.normal.vcf.gz'" >> parallelized_calls_normal
echo "parallel -a 3.scaff 'bcftools mpileup -Ou -f $MASTERGEN -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.normal.vcf.gz'" >> parallelized_calls_normal
echo "parallel -a 4.scaff 'bcftools mpileup -Ou -f $MASTERGEN -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.normal.vcf.gz'" >> parallelized_calls_normal
echo "parallel -a 5.scaff 'bcftools mpileup -Ou -f $MASTERGEN -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.normal.vcf.gz'" >> parallelized_calls_normal
echo "parallel -a 6.scaff 'bcftools mpileup -Ou -f $MASTERGEN -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.normal.vcf.gz'" >> parallelized_calls_normal
echo "parallel -a 7.scaff 'bcftools mpileup -Ou -f $MASTERGEN -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.normal.vcf.gz'" >> parallelized_calls_normal
echo "parallel -a 8.scaff 'bcftools mpileup -Ou -f $MASTERGEN -r {} -b $modernbams | bcftools call -m -v -Oz -o {}.normal.vcf.gz'" >> parallelized_calls_normal

python2 $HOME/bin/launcher_creator.py -j parallelized_calls_normal -n parallelized_calls_normal -a tagmap -e cbscott@utexas.edu -t 24:00:00 -N 2 -w 4 -q normal

####CREATE SNP FILE########
###the extensioms here will change based on
gunzip *.gz
>properorder #there are 24 coral scaffolds numbered 1-24.
for i in {1..24}; do
echo "${i}.omit.vcf" >> properorder
done

#I should really make a function out of this at some point. Could definitely be done by  a shell script
cat properorder | tr '\n' ' ' | xargs bcftools concat > modern_omit.vcf
#get just the biallelic snps (is this correct?)
bcftools view -m2 -M2 -v snps modern_omit.vcf | sed '/^#/d'| cut -f1,2,4,5 > modern_all_reqcols.txt
wc -l modern_all_reqcols.txt #write this number down.

#can we try this with more than biallelic snps for better results?
bcftools view modern_all.vcf | sed '/^#/d'| cut -f1,2,4,5 > modern_biplus_reqcols.txt
wc -l modern_biplus_reqcols.txt #1169251

#make this a job
>snpids
>genpos
for i in {1..1169251}; #replace this value with the result of the above 'wc -l' Could be automated in future
do
    echo "biplus_modern${i}" >> snpids;
    echo "0" >> genpos;
done

cut -f1 modern_biplus_reqcols.txt  > snpchr
cut -f2,3,4 modern_biplus_reqcols.txt > snpposrefalt
paste snpids snpchr genpos snpposrefalt > modern.biplus.snp


####CHANGE THIS VARIABLES DEPENDING ON YOUR DATASET-
#use the filtered bams for Q>30
export modernbams=/scratch/06909/cbscott/aDNA_revamp/sorted_bams/modern_coral_bams/modern_bams
export GEN_NEW=$STOCKYARD/aDNA_mastergen/mastergenome/human-fied_genome/concat_master_genome.fasta
export SNP=/scratch/06909/cbscott/aDNA_revamp/sorted_bams/modern_coral_bams/modern.snp
export ancbam_list=/scratch/06909/cbscott/aDNA_revamp/sorted_bams/modern_coral_bams/ancient_bams

#need to get names of samples, did not want to carry extensions through analysis
>apalm_names
for file in `cat $modernbams`; do
  echo ${file/.softclip.coral.sorted.bam} >> modern_names; done

#Here we actually do our pseudohaplo calls
echo "samtools mpileup -R -B -q10 -Q0 -f $GEN_NEW -b $modernbams | pileupCaller --randomHaploid --sampleNames `cat modern_names | tr '\n' ',' | sed 's/.$//'` -f $SNP -e modern" > modtest #don't use -e modern option.... what does this create?
python2 $HOME/bin/launcher_creator.py -j modtest -n modtest -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 1 -q skx-dev

echo "samtools mpileup -R -B -q10 -Q0 -f $GEN_NEW -b $ancbam_list | pileupCaller --randomHaploid --sampleNames `cat anc_names | tr '\n' ',' | sed 's/.$//'` -f $SNP -e ancient_all" > anctest #don't use -e modern option.... what does this create?
python2 $HOME/bin/launcher_creator.py -j anctest -n anctest -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 1 -q skx-dev
#computer is just crazy slow

#Use the Rscript - jackknife.R in order to construct chromosome jack-knifed dataset
#then pipe that data set into mergeit.
#very fast
mergeit -p parMergeIt #this comes from EIGENSOFT

#I think modern poputions should be in poplist... others shouldn't be given?
#otherwise how does it know what to project onto?
echo "Unknown" > poplist.txt  #this could be refined by assigning Ancient vs. Modern pop lists
smartpca -p parSmartPCA #Also from EIGENSOFT?

#grab just the data from the first 10 PCs
sed 1d modern_merged.evec | tr -s ' ' | sed -e 's/^[ \t]*//' | awk 'NF{NF-=1};1' | awk -v OFS="\t" '$1=$1' | sed -e '1s/^/Sample_ID\tpc1\tpc2\tpc3\tpc4\tpc5\tpc6\tpc7\tpc8\tpc9\tpc10\n/' > modern_merged.evec.txt
# Move the resulting file locally to complete analysis



###### ATTEMPT ANGSD HAPLOID CALLS + PCANGSD AS IN STICKLEBACK PAPER #########
###Do this bad boy for ancient and modern files together....####
#workingdir=/scratch/06909/cbscott/aDNA_revamp/angsd_test

#use -rf to specify a region file- let's break this into 4 jobs which can run on 4 nodes for speed.
echo "angsd -bam bamlist -dohaplocall 1 -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -trim 1  -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out angsd_bigrun -nThreads 8" > angsd_bigrun
python2 $HOME/bin/launcher_creator.py -j angsd_bigrun -n angsd_bigrun -a tagmap -e cbscott@utexas.edu -t 10:00:00 -N 1 -w 1 -q skx-normal

#ran out of memory :(

./angsd -GL 1 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam bam.filelist


#lets do a test just for chr11chr11:1-
echo "angsd -bam bamlist -dohaplocall 1 -r 1:1- -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -trim 1  -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out angsdtestsmall -nThreads 24" > angsd_testsmall
python2 $HOME/bin/launcher_creator.py -j angsd_testsmall -n angsd_testsmall -a tagmap -e cbscott@utexas.edu -t 01:00:00 -N 1 -w 1 -q skx-dev
#won't work with intel issues

#only making calls at relevant sites for ancient samples?
#Let's only make calls at sites covered by at least one of the four ancient samples.
#We'll tohaplocall to get this
echo "angsd -bam bamlist -dohaplocall 1 -doCounts 1 -minMinor 1 -trim 1 -out halpocalls_all" > all_haplocalls
python2 $HOME/bin/launcher_creator.py -j all_haplocalls -n all_haplocalls -a tagmap -e cbscott@utexas.edu -t 03:00:00 -N 1 -w 1 -q skx-normal


'misc/haploToPlink angsdput.haplo.gz outputname'

echo "angsd -bam bamlist -rf anc_callable.region -dohaplocall 1 -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -trim 1  -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out angsd_anconly" > angsd_anc
python2 $HOME/bin/launcher_creator.py -j angsd_anc -n angsd_anc -a tagmap -e cbscott@utexas.edu -t 00:10:00 -N 1 -w 1 -q skx-dev
python3 /home1/06909/cbscott/pcangsd/pcangsd.py -beagle angsd_anconly.beagle.gz

#only made calls at 27 sites? will this be enough to tell species apart?
#let's make a region file of chr:site
#SKX NODES like a million times faster than regular!!!

##VERY USEFUL:
#SUBSET BASED ON ROWNAMES
sed 's/:/_/g' anc.regions > formatted.regions #added a header line to formatted.regions
awk 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' formatted.regions angsd_big_sub1.beagle > angsd_big_sub2.beagle

#head -1 angsd_bigtest.beagle > header
#cat header subset_bigtest_nohead > angsd_big_sub.beagle
gzip angsd_big_sub.beagle

python3 /home1/06909/cbscott/pcangsd/pcangsd.py -beagle angsd_big_sub2.beagle.gz

#let's do admixture on this
#first three PCS because something interesting seems to be going on
python3 /home1/06909/cbscott/pcangsd/pcangsd.py -beagle angsd_big_sub2.beagle.gz -e 3 -o big_sub_admix -admix -admix_alpha 50

#created a subsetted beagle file for sites where there was data in the ancient file
python3 /home1/06909/cbscott/pcangsd/pcangsd.py -beagle chr1_subset.beagle.gz


#We're going to use PCangsd, but we're only going to restrict calls to positions that even have data in the ancient samples
#do counts to get counts at each base
angsd -out out -doCounts 1 -dumpCounts 2 -trim 1 -bam ancient_bams
#now we'll restrict our analysis  to just these positions with coverage > 0
#create a region and position file
#this is "ancdata.regions"
echo "angsd -bam bamlist -rf anc_data.regions -dohaplocall 1 -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -trim 1  -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out full_anc_restricted -nThreads 8" > angsd_restrict_full
python2 $HOME/bin/launcher_creator.py -j angsd_restrict_full -n angsd_restrict_full -a tagmap -e cbscott@utexas.edu -t 02:00:00 -N 1 -w 1 -q skx-normal
#why are we using a p-value filter for snps?
#ancient reads really should have been trimmed for 1 bp prior to mapping... this may also increase the number of reads we retain! - a later issue
#in theory could parallelize regions by chromosome
###FOR SOME REASON ANGSD IS WAY SLOWER WITH THESE REGION FILES
for i in {1..24}; do
grep "^${i}:" anc_data.regions > ${i}.regions; done

for i in {1..24}; do
  echo $i && wc -l ${i}.regions; done

>angsd_parallel
for i in {1..24}; do
  echo "angsd -bam bamlist -rf ${i}.regions -dohaplocall 1 -doCounts 1 -minMinor 1 -GL 1 -doGlf 2 -trim 1  -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -out par_anc_restricted -nThreads 4" >> angsd_parallel;
done
python2 $HOME/bin/launcher_creator.py -j angsd_parallel -n angsd_parallel -a tagmap -e cbscott@utexas.edu -t 24:00:00 -N 12 -w 2 -q skx-normal

###ANGSD HAS A SOFTWARE TO CONVERT A HALPOTYPE FILE TO A PLINK FILE!
#then, need to replace all "N's" with zeros using sed
#finally, can use
#first replace calls in haplocalls_all.haplo.gz

#ah! I forgot to subset... can go back and subset .tped file if needed


haploToPlink halpocalls_all.haplo.gz haploall
sed -i 's/N/0/g' haploall.tped # ideally do ths before haploToPlink
plink --tfile haploall --recode --out haploped#to get the binary needed by pcangsd
#then, need
plink --file plink --make-bed --out haploall #to get right format

echo "python3 /home1/06909/cbscott/pcangsd/pcangsd.py -plink haploall" > haploall_pca
python2 $HOME/bin/launcher_creator.py -j haploall_pca -n haploall_pca -a tagmap -e cbscott@utexas.edu -t 00:30:00 -N 1 -w 1 -q skx-dev

#let's do the subset really quick too (different file structure than last subset)
awk 'NR==FNR{a[$1]; next} FNR==1 || $2 in a' formatted.regions haploall.tped > haplosub.tped
plink --tfile haplosub --recode --out haplosub
plink --file haplosub --make-bed --out haplosub
python3 /home1/06909/cbscott/pcangsd/pcangsd.py -plink haplosub -out haplosub


#PCAngsd? #be sure to run with python3 to avoid type errors
#module load intel
module load gcc
python3 /home1/06909/cbscott/pcangsd/pcangsd.py -beagle angsd_long.beagle.gz

####NGSadmix
for i in {2..7}; do
NGSadmix -likes angsd_big_sub2.beagle.gz -K $i -o admix_big_sub.${i}; done


####DO leave one out
#need helper file = JackKnife
mergeit -p parMissingMergeIt
sed 1d anc_mod_missing.evec | tr -s ' ' | sed -e 's/^[ \t]*//' | awk 'NF{NF-=1};1' | awk -v OFS="\t" '$1=$1' | sed -e '1s/^/Sample_ID\tpc1\tpc2\tpc3\tpc4\tpc5\tpc6\tpc7\tpc8\tpc9\tpc10\n/' > anc_mod_missing.evec.txt
