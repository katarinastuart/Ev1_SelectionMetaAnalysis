# Outlier Analysis Workshop

There are lots of interesting patterns that you can extract from genetic marker data. This can include patterns of linkage, balancing selection, or even inbreeding signals. One of the most common ones is to try find sites on the genome that are under divergent selection. The following vignette will take you through the basics of genetic selection analysis. 

The project has been funded by <a href="https://ausevo.com/ECR_grants_2022/">the AES ERC Networking Grant Scheme</a> and <a href="https://genetics.org.au/">GSA</a>.

<h2>Schedule</h2>

**Day 1**<br>
  9:00am Introduction Slides<br>
  9:30am Download data<br>
  10:00am Morning Tea<br>
  10:15am PCAdapt<br>
  12:00pm Lunch<br>
  1:00pm VCFtools<br>
  2:00 Afternoon tea<br>
  2:15pm VCFtool continued and setup for Bayescan and Baypass<br>
  
**Day 2**<br>
  9:00am Bayscan<br>
  10:00am Morning Tea<br>
  10:15am Bayescan continued & Bayepass<br>
  12:00pm Lunch <br>
  1:00pm Baypass<br>
  2:00pm Afternoon tea<br>
  2:15pm Compiling results, group discussiona and metanalysis contribution  

<h2>A <i>fairly</i> brief introduction to Genetic Outlier and Association Analysis</h2>

When we look through a genome to try find loci that are under divergent selection, we often conduct what is called outlier or association analyses. **Outlier analysis** requires just knowledge of the genetics of your samples (plus sample metadata, for example population groupings), and tries to find loci that behave very differently from the underlying patterns across the genome (with the assumption being that the rest of the genome represents patterns of neutral genetic diveristy)<sup>A</sup>. Meanwhile **association analysis** require some sort of covariate data, and tests whether there are any genetic variants statistically associated with this new data (you may have heard the term [GWAS](https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies)). This data can come in the form of phenotype data (e.g. morphology, disease status, physiology measures), or could be spatial (e.g. environmental, climate). Association tests look for sites in the genome where the precense or absence of a variant is highly correlated with the values in the co-variate data, usually through some regression type analysis.

<sup>A</sup> it is very important then to account for any population substructure.

> :round_pushpin: **An important note about additive genetic variance**
>
> It is important to bear in mind how the input genetic data for outlier or association models is being interreted by the model. When dealing with many of these models (and input genotype files) the assumption is that the SNP effects are [additive](https://link.springer.com/referenceworkentry/10.1007/978-3-319-47829-6_5-1). This can be seen from, for example, the way we encode homozygous reference allele, heterozygous, and homozygous alternate allele as "0", "1", and "2" respectively in a BayPass input genofile. For the diploid organism (with two variant copies for each allele) one copy of a variant (i.e. heterozygous) is assumed to have half the effect of having two copies. However, what if the locus in question has dominance effects? This would mean the heterozygous form behaves the same as the homozygous dominant form, and it would be more appropriate to label these instead as "0", "0", "1". But with thousands, if not millions of (most likely) completely unknown variants in a dataset, how can we possibly know? The answer is we cannot. And most models will assume additive effects, because this the simplest assumption. However, by not factoring in dominance effects we could possible be missing many important functional variants, as Reynolds et al. [2021](https://www.nature.com/articles/s41588-021-00872-5) demonstrates. Genomics is full of caveats and pitfalls, which while providing new directions to explore can be a bit overwhelming. Remember, you selection analysis doesn't have to be exhaustive, just make sure it is as fit for purpose within your study design. There is so much going on in just one genome, there is no way you can analyse everything in one go. 

Throughout this vignette I will refer to outlier and association analysis collectively as selection analysis. There are a lot of programs that exist currently. You can find a very long (though not exhaustive) list of them [here](https://bioinformaticshome.com/tools/gwas/gwas.html). 

This vignette still start by covering some very simple outlier analyses:
<ul>
<li><a href="https://bcm-uga.github.io/pcadapt/articles/pcadapt.html"><b>PCAdapt</b></a>, a program that can depect genetic marker outliers without having population<sup>B</sup> designations, using a Principle Component Analysis (PCA) approach.</li>
<li><a href="GBIF"><b>F<sub>ST</sub></b></a> outlier analysis, an approach that uses pairwise comparisons between two populations and the fixation index metric to assess each genetic marker.</li>
</ul>

Next we will conduct some more advanced outlier analysis:
<ul>
<li><a href="https://bcm-uga.github.io/pcadapt/articles/pcadapt.html"><b>Bayescan</b></a>, a program.</li>
<li><a href="GBIF"><b>Baypass</b></a>, a program that elaborates on the bayenv model (another popular association analysis program).</li>
</ul>

<sup>B</sup> I will say population for simplicity throughout this vignette. However, equally we can test for differences between sample sites, subpopulations, and other types of groupings. What counts as one 'group' of organisms will be dependent on your study system or study question.

We'll cover the pre-processing of program specific input files, how to run the programs, how to visualise the output and also in some cases we'll need to take extra steps to map the genetic markers of interest back to the SNP data.  

> :round_pushpin: **An important note about reduced representation verses whole genome sequencing**
>
> Completing outlier analysis is possible and often done on reduced representation data. It is important to remember how your genome coverage (the number of genome variant sites / the genome length <sup>C</sup>) will affect your results and interpretation. Often with WGS data, you will see well resolved 'peaks' with a fairly smooth curve of points leading up to it either side. From this we often infer that the highest point is the genetic variant of interest, and the other sites either side of that exhibit signals of selection because they reside close to, and thus are linked, to the variant of interest. However, consider that even in WGS data, unless we have every single genetic variant represented (which may not be the case, depending on our variant calling and filtering parameters) it is possible that the genetic variant of interest that we have identified is not the main one, but is simply another neighbouring linked SNP to one that is not represented in the data. This problem becomes even more relevant with reduced representation sequencing (RRS), for which the genome coverage may be extremely patchy<sup>C</sup>. Thus with all outlier analysis, but especially so for those using RRS data, remember that your flagged outliers are not exhaustive, and may themselves only be liked to the variant that is truely under selection.

<sup>C</sup> You may also want to consider linkage blocks.


https://onlinelibrary.wiley.com/doi/10.1111/mec.14549

https://pubmed.ncbi.nlm.nih.gov/29486732/

 https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.9176

![ScreenShot](https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/asset/4be56b5b-8593-4116-ab9a-ec7b9e3c9a05/gr1.jpg)

also include  an hour of group discussion on one of the days about research question framing + grant integration






## Define you working directory for this project, and the VCF file location:

```
mkdir ~/outlier_analysis
DIR=~/outlier_analysis
cd $DIR
```

Our data tree will look like:

> outlier_analysis/ <br>
> ├── analysis <br>
> │   ├── bayescan <br>
> │   ├── baypass <br>
> │   ├── pcadapt <br>
> │   ├── summary <br>
> │   └── vcftools_fst <br>
> ├── data  <br>
> ├── programs  <br>
> └── workshop_material <br>


So lets set up our directories to match this

```
mkdir -p {analysis/{bayescan,baypass,pcadapt,summary,vcftools_fst},data,programs,workshop_material}
```

## Project data

The data provided in this workshop contains 5007 SNPs loci for across 39 individuals (13 individuals each from 3 different locations). There is some missingness (i.e. missing SNP calls) within this data.

There is also a metadata file, that contains the individuals unique IDs, their assigned populations, and a wingspan measurement for each individual. 

Let's grab this data from the project's git resository, place the data files into our ``data`` directory, and define the environmental variables ``VCF`` and ``METADATA`` with the locations of the genetic variant and metadata files respectively.

```
cd $DIR/workshop_material
git clone https://github.com/katarinastuart/Ev1_SelectionMetaAnalysis.git
cp $DIR/workshop_material/Ev1_SelectionMetaAnalysis/workshop_files/* $DIR/data
VCF=$DIR/data/starling_3populations.recode.vcf
METADATA=$DIR/data/starling_3populations_metadata.txt
```

> :heavy_exclamation_mark: **Working with your own data** <br> 
> <br>
> Alternatively, you can also use your own data for this workshop. If so, it is a good idea to thin your SNP dataset down to roughly 5,000 SNPs to ensure compute times are not too long. If you have more than 50 individuals you may also want to reduce this too. If you would like to do this, just place your genetic variant and metadata file in the ``data`` directory and define ``VCF`` and ``METADATA`` based on their names. <br>


Across this workshop, we will need the genetic data to be in several different formats. Let's prepare that now. First we convert the VCF to PLINK, and then to BED.

```
cd $DIR/data
module load vcftools/0.1.16
module load plink/1.90b6.7 
vcftools --vcf $VCF --out starling_3populations.plink --plink
plink --file starling_3populations.plink --make-bed --noweb --out starling_3populations
```



## PCAdapt

The PCAdapt manual is available [here](https://bcm-uga.github.io/pcadapt/articles/pcadapt.html).

Brief summary of PCAdapt. [FIX]

Install PCAdapt and set your working directory.

```
module load R/3.5.3
R

install.packages("pcadapt")
library(pcadapt)

setwd("/home/z5188231/outlier_analysis/analysis/pcadapt/")
```

Now let's load in the data - PCAdapt uses bed file types.

```
starling_bed <- "/home/z5188231/outlier_analysis/data/starling_3populations.bed"
starlings_pcadapt <- read.pcadapt(starling_bed, type = "bed")
```

Produce K plot

```
starlings_pcadapt_kplot <- pcadapt(input = starlings_pcadapt, K = 20)
pdf("pcadapt_starlings_kplot.pdf")
plot(starlings_pcadapt_kplot, option = "screeplot")
dev.off()
```

<img src="/images/pcadapt_kplot.PNG" alt="k plot" width="400"/>

K value of 3 is most appropriate, as this is the value of K after which the curve starts to flatten out more.

```
starlings_pcadapt_pca <- pcadapt(starlings_pcadapt, K = 3)
summary(starlings_pcadapt_pca)
```

> :heavy_check_mark: **output** <br>
> &nbsp; &nbsp; &nbsp; &nbsp; Length Class  Mode <br>
> scores &nbsp; 117  -none- numeric <br>
> singular.values &nbsp;     3  -none- numeric <br>
> loadings &nbsp;        15021  -none- numeric <br>
> zscores &nbsp;         15021  -none- numeric <br>
> af &nbsp;               5007  -none- numeric <br>
> maf &nbsp;              5007  -none- numeric <br>
> chi2.stat &nbsp;       5007  -none- numeric <br>
> stat &nbsp;             5007  -none- numeric <br>
> gif &nbsp;                 1  -none- numeric <br>
> pvalues &nbsp;          5007  -none- numeric <br>
> pass &nbsp;             4610  -none- numeric 

Investigate axis projections:

```
poplist.names <- c(rep("Lemon", 13),rep("Warrnambool", 13),rep("Nowra", 13))
print(poplist.names)

pdf("pcadapt_starlings_projection1v2.pdf")
plot(starlings_pcadapt_kplot, option = "scores", i = 1, j = 2, pop = poplist.names)
dev.off()

pdf("pcadapt_starlings_projection5v7.pdf")
plot(starlings_pcadapt_kplot, option = "scores", i = 5, j = 7, pop = poplist.names)
dev.off()
```

Ignore the warning: 

> :heavy_exclamation_mark: Use of `df$Pop` is discouraged. Use `Pop` instead.

<img src="/images/pcadapt_proj1.PNG" alt="projection axis1 axis2" width="300"/> <img src="/images/pcadapt_proj2.PNG" alt="projection axis6 axis7" width="300"/>

Investigate manhattan and Q-Qplot:

> :beginner: **Manhattan plots** are a way to visualise the GWAS (genome-wide association study) p-values (or other statistical values) at each SNP locus along the genome

> :beginner: **Q-Qplots plots** are just a quick way to visually check if your residuals are normally distributed. Check out more information [here](https://data.library.virginia.edu/understanding-q-q-plots/).

```
pdf("pcadapt_starlings_manhattan.pdf")
plot(starlings_pcadapt_pca, option = "manhattan")
dev.off()

pdf("pcadapt_starlings_qqplot.pdf")
plot(starlings_pcadapt_pca, option = "qqplot")
dev.off()
```

<img src="/images/pcadapt_manhattan.PNG" alt="Manhattan" width="400"/> <img src="/images/pcadapt_qq.PNG" alt="Q-Qplot" width="400"/>

Plotting and correcting the pvalues

```
starling_pcadapt_pvalues <- as.data.frame(starlings_pcadapt_pca$pvalues)

library("ggplot2")

pdf("pcadapt_starlings_pvalues.pdf")
hist(starlings_pcadapt_pca$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
dev.off()
```

<img src="/images/pcadapt_pvals.PNG" alt="pvals" width="450"/>

```
starlings_pcadapt_padj <- p.adjust(starlings_pcadapt_pca$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(starlings_pcadapt_padj < alpha)
length(outliers)

write.table(outliers, file="starlings_pcadapt_outliers.txt") 

length(outliers)
```

> :heavy_check_mark: **Output** <br>
> &emsp;
> [1] 3

After this, we will be jumping out of R and back into the command line by using the command: 

```
q()
```

Mapping Outliers: PCAdapt

finding the SNP ID of the outlier variants

```
cd $DIR/analysis
```

The first thing we will do is create list of SNPs in VCF, assign line numbers that can be used to find matching line numbers in outliers (SNP ID is lost in PCadapt & Bayescan, line numbers used as signifiers). 

We create this in the ``analysis`` folder because we will use it for more than just mapping the outlier SNPs for PCAdapt.

```
grep -v "^#" $DIR/data/starling_3populations.recode.vcf | cut -f1-3 | awk '{print $0"\t"NR}' > starling_3populations_SNPs.txt
```

Now let's jump back into the ``pcadapt`` directory to contiue working with our outliers. We grab column 2 of the outlier file using the ``AWK`` command, which contain the number of the outliers

```
cd $DIR/analysis/pcadapt
awk '{print $2}' starlings_pcadapt_outliers.txt > starlings_pcadapt_outliers_numbers.txt
```

We now make a list of outlier SNPS ID's

```
awk 'FNR==NR{a[$1];next} (($4) in a)' starlings_pcadapt_outliers_numbers.txt ../starling_3populations_SNPs.txt   | cut -f3 > pcadapt_outlierSNPIDs.txt
head pcadapt_outlierSNPIDs.txt
```
> :heavy_check_mark: **Output** <br>
> &emsp;
> 230955:72:-
> 238881:46:+
> 286527:46:-



## VCFtools windowed Fst

The VCFTools manual is available [here](https://vcftools.sourceforge.net/man_latest.html).

Fst outliers will allow us to identify SNPs that behave abnormally in pairwise comparisons between populations.

The first things we need to do is use our metadata file (currently defined by the environmental variable ``METADATA``) to make three individual files containing just the list of individuals in each of the populations. We can do this by subseting our sample metadata file, using the command ``grep`` to grab lines that match each population's name, and then using ``awk`` to keep only the first column of metadta, i.e. the sample names.

```
module load vcftools/0.1.16
```

```
cd $DIR/data

grep "Lemon" $METADATA | awk '{print $1}' > individuals_Lemon.txt
grep "War" $METADATA | awk '{print $1}' > individuals_War.txt
grep "Nowra" $METADATA | awk '{print $1}' > individuals_Nowra.txt
```

Now we can pick two populations to compare. Let's work with Lemon (short for Lemon Tree, QLD, AU) and War (short for Warnambool, VIC, AU), and so a SNP-based Fst comparison.

```
cd $DIR/analysis/vcftools_fst

vcftools --vcf $VCF --weir-fst-pop $DIR/data/individuals_Lemon.txt --weir-fst-pop $DIR/data/individuals_War.txt --out lemon_war

head -n 5 lemon_war.weir.fst
```
> :heavy_check_mark: **Output** <br>
> &emsp;
> CHROM   POS     WEIR_AND_COCKERHAM_FST<br>
> starling4       107735  0.160891<br>
> starling4       137462  -0.0805785<br>
> starling4       151332  0.0524489<br>
> starling4       227887  0.0569961<br>

The important column is column 5: the Weighted Fst, from [Weir and Cockerham’s 1984 publication](https://www.jstor.org/stable/2408641). This corrects for INFO NEEDED.

```
wc -l lemon_war.weir.fst 
```

> :heavy_check_mark: **Output** <br>
> &emsp;
> 5008

Notice how there are as many lines as there are SNPs in the data set, plus one for a header. It is always a good idea to check your output, and make sure everything looks as you expect it to!

Next, instead of calculating pairwise populaiton differentiation on a SNP by SNP basis, we will be usiing a sliding window approach. The ``--fst-window-size 50000`` refers to the window size of the genome (in bsae pairs) in which we are calculating one value: all SNPs within this window are used to caluclate Fst. The ``--fst-window-step`` option indicates how many base pairs the window is moving down the genome before calculating Fst for the second window, and then the third, and so on. 

> **Warning**
> &emsp;
> These sliding windows only work on ordered SNPs on the same chromosome/scaffold/contig. If you data is not set up like this (i.e. all your SNPs are on a single pseudo chromosome) then this method is not appropriate for your data, as it will be making an assumption about where the SNPs are located with respect to one another.

```
vcftools --vcf $VCF --fst-window-size 50000 --fst-window-step 10000 --weir-fst-pop $DIR/data/individuals_Lemon.txt --weir-fst-pop $DIR/data/individuals_War.txt --out lemon_war

head lemon_war.windowed.weir.fst
```
> :heavy_check_mark: **Output** <br>
> &emsp;
> CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST <br>
> starling4       60001   110000  1       0.160891        0.160891 <br>
> starling4       70001   120000  1       0.160891        0.160891 <br>
> starling4       80001   130000  1       0.160891        0.160891 <br>
> starling4       90001   140000  2       -0.00374291     0.0401563 <br>


Notice the output is different.

```
wc -l lemon_war.windowed.weir.fst 
```

> :heavy_check_mark: **Output** <br>
> &emsp;
> 10838

Notice the line count is different from the SNP-based Fst comparison; there are more lines in the sliding window based Fst comparison. This is because there are more sliding windows across the chromosome in this data set than there are SNPs. Consider which of these steps is better for your data: in low density SNP datasets, the sliding window approach might not be the best to use.

Now let's plot the Fst across the chromosome. To do this we will add line numbers on our Fst file that will be used to order the Fst measurements across the x-axis of our manhattan plot.

> :beginner: **X-axis values** in the following plot are done y using each ourlier window's line number, as they are in order along the genome. Ourlier windows are equally spces, and so line numbers are sufficient to capture the patterns along the genome. Consider that if you are plotting Fst values for SNPs (rather than windows), they may not be equally spaced along the genome and so SNP position may need to be used to make your manhattan plots.

```
awk '{print $0"\t"NR}' ./lemon_war.windowed.weir.fst  > lemon_war.windowed.weir.fst.edit

module load R/3.5.3
R

library("ggplot2")

setwd("/home/z5188231/outlier_analysis/analysis/vcftools_fst")

windowed_fst <- read.table("lemon_war.windowed.weir.fst.edit", sep="\t", header=TRUE)
str(windowed_fst)

quantile(windowed_fst$WEIGHTED_FST, probs = c(.95, .99, .999))
```

> :heavy_check_mark: **Output** <br>
> &emsp;
>       95%       99%     99.9%<br>
> 0.1948850 0.3501600 0.5741306<br>


Choose the quantile threshold above which SNPs will be classified as outliers. Below, we chose 99% (i.e. the top 1% of SNP windows).

```
pdf("fst_starlings_windowed.pdf", width=10, height=5)
ggplot(windowed_fst, aes(x=X1, y=WEIGHTED_FST)) + 
geom_point() + 
theme_classic() +
geom_hline(yintercept=0.35, linetype="dashed", color = "red")
dev.off()

q()
```

<img src="/images/Fst_Windowed.PNG" alt="Windowed Fst" width="600"/>

Finally, we will generate a list of outier SNP IDs. We do this by grabbing all of the SNPs located in the ourlier windows. 

```
cd $DIR/analysis/vcftools_fst
cat lemon_war.windowed.weir.fst.edit | awk '$5>0.3501600' > lemon_war.windowed.outliers
wc -l lemon_war.windowed.outliers
```
> :heavy_check_mark: **Output** <br>
> &emsp;
> 107 lemon_war_fst.windowed.outliers

```
cut -f1-3 lemon_war.windowed.outliers > lemon_war.windowed.outliers.bed 
vcftools --vcf $VCF --bed lemon_war.windowed.outliers.bed --out fst_outliers --recode
grep -v "#" fst_outliers.recode.vcf | awk '{print $3}' > vcftoolsfst_outlierSNPIDs.txt
wc -l vcftoolsfst_outlierSNPIDs.txt
```
> :heavy_check_mark: **Output** <br>
> &emsp;
> 61

We have a total of 61 outlier SNPs locate across 107 outlier SNP windows.


## Bayescan

The VCFTools manual is available [here](https://github.com/mfoll/BayeScan).

Bayescan identified outlier SNPs based on allele frequencies. More explination about alpha and such.

First, we will need to convert out VCF to the Bayescan format. To do this we will use the genetic file conversion program called [PGDspider](http://www.cmpg.unibe.ch/software/PGDSpider/). 

```
cd $DIR/programs
wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.1.1.5.zip
```

We now run PGDSpider in two steps: first we convert the VCF file to the PGD format, second from PGD format to Bayescan format. To do this we will need to create a SPID file. create a file called *VCF_PGD.spid* using the ``nano`` command. Paste in the below, replacing the location of you metadata file.

include snapshot of SPID:


> \# VCF Parser questions <br>
> PARSER_FORMAT=VCF <br>
> \# Only output SNPs with a phred-scaled quality of at least: <br>
> VCF_PARSER_QUAL_QUESTION= <br>
> \# Select population definition file: <br>
> VCF_PARSER_POP_FILE_QUESTION=/srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/data/3pops_pops.txt <br>
> \# What is the ploidy of the data? <br>
> VCF_PARSER_PLOIDY_QUESTION=DIPLOID <br>
> \# Do you want to include a file with population definitions? <br>
> VCF_PARSER_POP_QUESTION=true <br>
> \# Output genotypes as missing if the phred-scale genotype quality is below: <br>
> VCF_PARSER_GTQUAL_QUESTION= <br>
> \# Do you want to include non-polymorphic SNPs? <br>
> VCF_PARSER_MONOMORPHIC_QUESTION=false <br>
> \# Only output following individuals (ind1, ind2, ind4, ...): <br>
> VCF_PARSER_IND_QUESTION= <br>
> \# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated): <br>
> VCF_PARSER_REGION_QUESTION= <br>
> \# Output genotypes as missing if the read depth of a position for the sample is below: <br>
> VCF_PARSER_READ_QUESTION= <br>
> \# Take most likely genotype if "PL" or "GL" is given in the genotype field? <br>
> VCF_PARSER_PL_QUESTION=false <br>
> \# Do you want to exclude loci with only missing data? <br>
> VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false <br>
>  <br>
> \# PGD Writer questions <br>
> WRITER_FORMAT=PGD <br>


```
head FILE
```


>  au05_men        SOUTH <br>
>  au06_men        SOUTH



Now run the two step convserion.

```
cd $DIR/analysis/bayescan

java -Xmx1024m -Xms512m -jar $DIR/programs/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile $VCF -inputformat VCF -outputfile starling_3populations.pgd -outputformat  PGD -spid VCF_PGD.spid 

java -Xmx1024m -Xms512m -jar $DIR/programs/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputfile starling_3populations.pgd -inputformat PGD -outputfile starling_3populations.bs -outputformat GESTE_BAYE_SCAN
```

Now let's set Bayescan to run.

```
#!/bin/bash
#PBS -N 2021-11-21.bayescan_starling.pbs
#PBS -V
#PBS -l nodes=1:ppn=16
#PBS -l mem=40gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M katarina.stuart@unsw.edu.au
#PBS -m ae

module load bayescan/2.1

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/bayescan

bayescan_2.1 ./starling_3populations.bs -od ./ -threads 16 -n 5000 -thin 10 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 10
```
 
Identify outliers:

```
module load R/3.5.3
R
library(ggplot2)
setwd("/srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/bayescan")
source("/apps/bayescan/2.1/R\ functions/plot_R.r")
outliers.bayescan=plot_bayescan("/srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/bayescan/starling_3population_fst.txt",FDR=0.05)
outliers.bayescan
write.table(outliers.bayescan, file="bayscan_outliers.txt")
```


Mapping Outliers

```
cd $DIR/analysis/bayescan
```

Create list of SNPs in VCF, assign line numbers that can be used to find matching line numbers in outliers (SNP ID is lost in bayescan, line numbers used as signifiers).

```
grep -v "^#" ../../data/starling_3populations.recode.vcf  | cut -f1-3 | awk '{print $0"\t"NR}' > starling_3populations_SNPs.txt

awk '{print $2}' bayscan_outliers.txt > bayscan_outliers_numbers.txt
```

list of outlier SNPS, by matching column 1 of of the outliers list to the fourth column of the whole SNP list data.

```
awk 'FNR==NR{a[$1];next} (($4) in a)' bayscan_outliers_numbers.txt starling_3populations_SNPs.txt   | cut -f3 > bayscan_outliers_SNPs.txt
```

Bayescane Log Plot, colouring the outliers in a different colour.

```
module load R/3.5.3
R
library(ggplot2)
library(dplyr)
setwd("/srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/bayescan")

bayescan.out<- read.table("starling_3population_fst.txt", header=TRUE)
bayescan.out <- bayescan.out %>% mutate(ID = row_number())
bayescan.outiers<- read.table("bayscan_outliers_numbers.txt", header=FALSE)
outliers.plot <- filter(bayescan.out, ID %in% bayescan.outiers[["V1"]])

png("bayescan_outliers.png", width=600, height=350)
ggplot(bayescan.out, aes(x=log10.PO., y=alpha))+
geom_point(size=5,alpha=1)+xlim(-1.3,3.5)+ theme_classic(base_size = 18) + geom_vline(xintercept = 0, linetype="dashed", color = "black", size=3)+
geom_point(aes(x=log10.PO., y=alpha), data=outliers.plot, col="red", fill="red",size=5,alpha=1) + theme(axis.text=element_text(size=18), axis.title=element_text(size=22,face="bold"))
dev.off()
```
<img src="/images/bayescan_outliers.png" alt="Windowed Fst" width="300"/>

## BayPass

The Baypass manual can be found [here](http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.31.pdf).

Baypass requires that the allele frequency data be on a population, not an individual basis. The genotyping data file is simply organized as a matrix with nsnp rows and 2 ∗ npop columns. The row field separator is a space. More precisely, each row corresponds to one marker and the number of columns is twice the number of populations because each pair of numbers corresponds to each allele (or read counts for PoolSeq experiment) counts in one population. To generate this population gene count data we will work with the PLINK file.

```
PLINK=$DIR/data/starling_3populations.plink.ped

cut -f 3- $PLINK > x.delete
paste $DIR/data/3pops_plink.txt x.delete > starling_3populations.plink.ped
rm x.delete 
cp ../../data/starling_3populations.plink.map .
cp ../../data/starling_3populations.plink.log .
```

run the pop based allele frequency calculations

```
plink --file starling_3populations.plink --allow-extra-chr --freq counts --family --out starling_3populations
```

manipulate file so it has baypass format, numbers set for plink output file and pop number for column count

```
tail -n +2 starling_3populations.frq.strat | awk '{ $9 = $8 - $7 } 1' | awk '{print $7,$9}' | tr "\n" " " | sed 's/ /\n/6; P; D' > starling_3populations_baypass.txt
```

Now we can run Baypass.

```
#!/bin/bash
#PBS -N 2022-12-08.baypass_starling.pbs
#PBS -l nodes=1:ppn=16
#PBS -l mem=40gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M katarina.stuart@unsw.edu.au
#PBS -m ae

module load baypass/2.1

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/baypass

g_baypass -npop 3 -gfile ./starling_3populations_baypass.txt -outprefix starling_3populations_baypass -nthreads 16
```

Running in R to make the anapod data

```
module load R/3.6.3
R
setwd("/srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/baypass")
source("/apps/baypass/2.1/utils/baypass_utils.R")
library("ape")
library("corrplot")

omega=as.matrix(read.table("starling_3populations_baypass_mat_omega.out"))
pi.beta.coef=read.table("starling_3populations_baypass_summary_beta_params.out",h=T)$Mean
bta14.data<-geno2YN("starling_3populations_baypass.txt")
simu.bta<-simulate.baypass(omega.mat=omega, nsnp=5000, sample.size=bta14.data$NN, beta.pi=pi.beta.coef,pi.maf=0,suffix="btapods")
```

We now have the simulated geno data.

```
#!/bin/bash
#PBS -N 2022-12-08.baypass_starling2.pbs
#PBS -l nodes=1:ppn=16
#PBS -l mem=40gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M katarina.stuart@unsw.edu.au
#PBS -m ae

module load baypass/2.1

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/baypass

g_baypass -npop 2 -gfile G.btapods  -outprefix G.btapods -nthreads 16 
```

XtX calibration; get the pod XtX

```
module load R/3.6.3
R
setwd("/srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/baypass")
source("/apps/baypass/2.1/utils/baypass_utils.R")
library("ape")
library("corrplot")

pod.xtx=read.table("G.btapods_summary_pi_xtx.out",h=T)$M_XtX
```

We compute the 1% threshold for the simulated neutral data.

```
pod.thresh=quantile(pod.xtx,probs=0.99)
pod.thresh
```

> :heavy_check_mark: **Output** <br>
> &emsp;
> 4.73302

Your values may be slightly different as the simlated data will not be identical.

Next, we filter the data for the outlier snps by identifying those above the threshold.

```
cat starling_3populations_baypass_summary_pi_xtx.out | awk '$6>4.73302 ' > baypass_outliers.txt
```

create list of SNPs in VCF, assign line numbers that can be used to find matching line numbers in outliers (SNP ID is lost in bayescan, line numbers used as signifiers).

```
grep -v "^#" ../../data/starling_3populations.recode.vcf  | cut -f1-3 | awk '{print $0"\t"NR}' > starling_3populations_SNPs.txt
```

List of outlier SNPS

```
awk 'FNR==NR{a[$1];next} (($4) in a)' baypass_outliers.txt starling_3populations_SNPs.txt | cut -f3 > baypass_outliers_SNPlist.txt
wc -l baypass_outliers_SNPlist.txt
```

> :heavy_check_mark: **Output** <br>
> &emsp;
> 275

Low lets find SNPs that are statistically associated with covarate data.

```
#!/bin/bash
#PBS -N 2022-12-08.baypass_starling3.pbs
#PBS -l nodes=1:ppn=16
#PBS -l mem=40gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M katarina.stuart@unsw.edu.au
#PBS -m ae

module load baypass

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/baypass

g_baypass -npop 3 -gfile starling_3populations_baypass.txt -efile baypass_environment_covariate.txt -scalecov -auxmodel -nthreads 16 -omegafile starling_3populations_baypass_mat_omega.out -outprefix starling_3populations_baypass_enviro
```


Next we plot the outliers. We are chosing a <a href="https://www.statology.org/bayes-factor/">BF threshold</a> of 20 dB, which indicates "Strong evidence for alternative hypothesis".

```
module load R/3.6.3
R

library(ggplot2)

setwd("/srv/scratch/z5188231/KStuart.Starling-Aug18/Ev1_SelectionMetaAnalysis/analysis/baypass")

covaux.snp.res.mass=read.table("starling_3populations_baypass_enviro_summary_betai.out",h=T)
covaux.snp.xtx.mass=read.table("starling_3populations_baypass_summary_pi_xtx.out",h=T)$M_XtX

pdf("Baypass_plots.pdf")
layout(matrix(1:3,3,1))
plot(covaux.snp.res.mass$BF.dB.,xlab="Mass",ylab="BFmc (in dB)")
abline(h=20, col="red")
plot(covaux.snp.res.mass$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covaux.snp.xtx.mass, xlab="SNP",ylab="XtX corrected for SMS")
dev.off()
```

<img src="/images/Baypass.PNG" alt="Baypass output" width="500"/>

Finally, lets generate the list of phenotype-associated SNP IDs. 

```
cat starling_3populations_baypass_enviro_summary_betai.out | awk '$6>20' > starling_3populations_baypass_enviro_BF20.txt
starling_3populations_baypass_enviro_BF20.txt
```

> :heavy_check_mark: **Output** <br>
> &emsp;
> 40

Filtering the data sets for SNPS above BFmc threshold and finding which are also divergent amongst populations

```
awk 'FNR==NR{a[$2];next} (($4) in a)' starling_3populations_baypass_enviro_BF20.txt starling_3populations_SNPs.txt | cut -f3 > starling_3populations_baypass_enviro_BF20_SNPlist.txt

comm -12 <(sort starling_3populations_baypass_enviro_BF20_SNPlist.txt) <(sort baypass_outliers_SNPlist.txt) > output.txt
```

> :heavy_check_mark: **Output** <br>
> &emsp;
> 38

## Comparing Outlier Overlap

Now we wil make an upset plot to compare the overlap of outliers detected over our different methods.

KATNOTE: fix names of last 3 files to match first 2 naming scheme (i.e. tool_outlierSNPIDs.txt)

```
cd $DIR/outlier_analysis/summary
cp $DIR/outlier_analysis/pcadapt/pcadapt_outlierSNPIDs.txt .
cp $DIR/outlier_analysis/pcadapt/vcftoolsfst_outlierSNPIDs.txt .
cp $DIR/outlier_analysis/pcadapt/bayscan_outliers_SNPs.txt .
cp $DIR/outlier_analysis/pcadapt/baypass_outliers_SNPlist.txt .
cp $DIR/outlier_analysis/pcadapt/starling_3populations_baypass_enviro_BF20_SNPlist.txt .
```

Now we have a copy of all the SNP IDs for each of out outlier analysis, let's use the R package <a href="https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html">UpSetR</a> to plot the overlap.

```
module load R/4.1.0-gimkl-2020a
R
setwd("/nesi/noackup/uoa02613/...")
pcadapt<-scan("pcadapt_outlierSNPIDs.txt", what = "", quiet=TRUE)
vcftools<-scan("vcftoolsfst_outlierSNPIDs.txt", what = "", quiet=TRUE)
bayescan<-scan("bayscan_outliers_SNPs.txt", what = "", quiet=TRUE)
baypass<-scan("baypass_outliers_SNPlist.txt", what = "", quiet=TRUE)
baypass_env<-scan("starling_3populations_baypass_enviro_BF20_SNPlist.txt", what = "", quiet=TRUE)  #total transcripts

all_outliers <- list(PCAdapt = pcadapt, VCFtools = vcftools, Bayescan = bayescan, Baypass = baypass, BaypassEnv = baypass_env)

#install.packages("UpSetR")
library(UpSetR)

pdf("All_outliers_upsetplot.pdf")
upset(fromList(all_outliers), order.by = "freq", empty.intersections = "on", point.size = 3.5, line.size = 2, mainbar.y.label = "Outlier Count", sets.x.label = "Total Outliers", text.scale = c(1.3, 1.3, 1, 1, 2, 1.3), number.angles = 30 ) 
dev.off() 
```

**THIS IS JUST A DUMMY UPSET PLOT** <p>


<img src="/images/All_outliers_upsetplot.PNG" alt="upset plot of outlier overlaps" width="500"/>


Some comments on the proper overlap results :)

## Outlier Analysis Metanalysis

Flow diagram.



## Funding 
<p align="center">

![ScreenShot](https://storage.corsizio.com/uploads/5cea29e798d9a757e03dba1c/events/6141a39502de2a7ff2f3f120/photo-_3UnQXVrb.jpg)

</p>

Thank you to the <a href="https://ausevo.com/ECR_grants_2022/">AES ERC Networking Grant Scheme</a> for funding this project.

## Project Contributors

https://www.dataschool.io/how-to-contribute-on-github/

