#WORKFLOW SUMMARY FOR CORAL aDNA
A summary of the aDNA_remote_workflow and helpful tips for processing your data.

## Software Needed
Samtools
angsd
EIGENSOFT
sratoolkit
gcc
bam2fastq
python
edirect
bowtie2
seqtk
mapDamage2.0
angsd
sequenceTools


## 1. Grab your data

### Get source files
Download your source files. Our coral files were initially mapped to the human genome to remove potential contamination. The resulting bams were converted back to fastqs using bam2fastq.

### Grab any modern reference files you may need
Use NCBI's edirect to get reads from bioproject. You will need to know the BioProject # only to retrieve both metadata and file names.

Use fastq-dump (from sratoolkit) to download reads.

### Grab your reference genome and format it.
Grab your reference genome(s) from NCBI. For our project, we used the Acropora millepora genome, Symbiodiniaceae genomes, and given metagenome associated genomes. I recommend concatenating all of these together to avoid spurious mappings.

Further, since the projection-based pipeline was originally formulated for human DNA, it can't handle genomes with >99 scaffolds. If you do not have a chromosome-level assembly, you should artificially concatenate your reference together into chromosomes. I recommend just naming these numbers - it will make your life easier later. **Important**: be sure to note the number of chromosomes you created somewhere. Also, create a file where each line is the name of a chromosome in (each) of your references. You will use these to divide your mapped genomes by organism later. *e.g. I have one file named Coral_Scaffolds, one named Zoox_Scaffolds, and one named MAG_names, which allows me to sort coral form symbiont from microbe later.*

## 2. Map and process your data

### Use bowtie2  & samtools to build your reference genome.
This is a computationally "cheap" job. It should only take ~3 hours or less depending on how many references you combined together. Then index the genome using 'samtools faidx'.

### Map ancient reads.
How you map your ancient reads will depend on how the sequencing libraries were constructed. Was the DNA repaired in any way? Our libraries were UDG-half treated, so we know that all ancient DNA damage is restricted to the terminal base pair of each read.

Regardless, you should do one mapping using the default end-to-end settings on bowtie2. You will need to analyze the mismatches in your reads to verify aDNA authenticity using mapDamage later.

In addition to this end to end mapping, map your files in a second way which minimizes terminal base damage. For us, we trimmed the terminal base off of each end of reads and mapped them a second time using bowtie2 end-to-end settings.

### Map modern reads (if you have them)
Map your modern comparison files in an appropriate manner. I used bowtie2's 'local' settings to increase the number of mappings I retained.

### Filter your bams for quality.
You generally will want to keep only "high quality" mapping reads (reads which are good matches to their positions). For my ancient reads, I used a quality of 30. This threshold can be somewhat arbitrary, and I initially tried several different cutoffs to see how it affected my resulting data (e.g., number of reads retained).

### Create a table of mapping by scaffold
Use samtools idxstats to see how many reads mapped to each scaffold. This is valuable in our microbial analysis, because it let us infer relative abundances of different species in our data. If you are only mapping to one species, you can probably skip this step.

## 3. Verify aDNA Damage
To verify aDNA damage, we will look at misincorporation patterns at the ends of reads and read length distribution. Here, you will want to separate your mapped bams back out by species. This is when you will use those named scaffold files from the genome creation in [1].

For each species, run mapDamage2.0 to calculate misincorporation rates at terminal ends of reads. mapDamage will ouput both a pdf for visual inspection, and the actual misincorporation rates as a text file. For few files, visually inspect the .pdf output based on mapDamage's documentation. For many files, I have provided an algorithm called "putatively_ancient_microbes.R" which will calculate how well the misincorporation rates match what you expect from aDNA. If you do not have UDG-half treated libraries DO NOT USE THIS SCRIPT! Even then, the script requires manually checking of pdf files after the fact to verify results. I would not trust it to return results unsupervised, but rather use it as a tool to filter down many results to just the most promising.

## 4. ANGSD pipeline to compare modern and ancient data
The ANGSD pipleline compares modern and ancient files by genotyping all of them together. You will then restrict the analysis to positions covered by at least one ancient file.

### Run ANGSD to genotype files and subset for positions covered in your ancient data.
First you'll create a file that has the coverage at each position of the genome for ancient files. Then, after running angsd to genotype your files, use this list of coverage to subset the resulting genotypes. Some details on how to do this using awk are contained in the remote_pipeline.

### Run NGSadmix to calculate admixture.
Use the subsetted angsd beagle file as an input to NGSadmix to get the resulting files used in "plot_admixture.R". I ran this for 2-7 clusters, which I then visually inspected in R.

## 5. Projection pipeline to compare modern and ancient data.
This pipeline will genotype all of your data (modern and ancient). Then, a PCA will be constructed based on modern variation. Finally, ancient samples are projected on this PCA. This serves to reduce bias from ancient samples when constructing your PCA (e.g., ancient samples may cluster just by virtue of being ancient). On the flip side, this method forces ancient samples to the context of modern day variation.


### Call genotypes
Use bcftools mpileup to call genotypes for each file across your reference genome.

### Do pseudohaplo calls on your data
Each input bam will be 'pseudohaplotyped'. This randomly selects one read from the source file at each genomic position. You will need a ".snp" file containing the snp names, positions, reference and alternative alleles. Details on how to create this from the resulting vcf  of the previous step are given in the remote_pipeline. This will return three files, a .geno file, a .ind file, and a .snp file. Make sure the populations in the .ind file are correctly assigned - **make sure your ancient samples have a different population name than your modern samples!**

### Jackknife chromosomes
Use the given cmd_JackScaffold and chr_namer R scripts to create jackknifed genotype and individual files. Merge your ancient and modern .ind and .geno files together. You can use the .snp file you already have.

Give these files to smartpca in a parameter file. I believe details on these parameter files can be found in the EIGENSOFT documentation.

### Pull down resulting .evec file and visualize in R.

## 6. Calculate f4 statistics for your data.

Using _non-jackknifed_ data (the .ind, .geno, and .snp files returned by pseudohaplocall step) to calculate f4 statistcs. You will need samples from an outgroup, two modern populations, and your ancient samples. Additionally, you will ned a poplist, containing the names of each population you want f4 statistics run on. These names must match the popnames in the .ind file. Each line should contain one population name, and spaces are __not__ permitted in names. Then, simply use EIGENSOFT's "qpDstat" command with a parameter file to calculate f4 statistics for all population combinations. A sample of this parameter file is given at the end of the remote_pipeline. 
