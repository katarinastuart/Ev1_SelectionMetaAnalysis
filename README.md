# Ev1_SelectionMetaAnalysis

There are lots of interesting patterns that you can extract from genetic marker data. This can include patterns of linkage, balancing selection, or even inbreeding signals. One of the most common ones is to try find sites on the genome that are under divergent selection. The following vignette will take you through the basics of genetic selection analysis. 

The project has been funded by the <a href="https://ausevo.com/ECR_grants_2022/">AES ERC Networking Grant Scheme</a>.

<h2>A <i>fairly</i> brief introduction to Genetic Outlier and Association Analysis</h2>

When we look through a genome to try find loci that are under divergent selection, we often conduct what is called outlier or association analyses. **Outlier analysis** requires just knowledge of the genetics of your samples (plus sample metadata, for example population groupings), and tries to find loci that behave very differently from the underlying patterns across the genome (with the assumption being that the rest of the genome represents patterns of neutral genetic diveristy)<sup>A</sup>. Meanwhile **association analysis** require some sort of covariate data, and tests whether there are any genetic variants statistically associated with this new data (you may have heard the term [GWAS](https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies)). This data can come in the form of phenotype data (e.g. morphology, disease status, physiology measures), or could be spatial (e.g. environmental, climate). Association tests look for sites in the genome where the precense or absence of a variant is highly correlated with the values in the co-variate data, usually through some regression type analysis.

<sup>A</sup> it is very important then to account for any population substructure.

> **An important note about additive genetic variance**
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

> **An important note about reduced representation verses whole genome sequencing**
>
> Completing outlier analysis is possible and often done on reduced representation data. It is important to remember how your genome coverage (the number of genome variant sites / the genome length <sup>C</sup>) will affect your results and interpretation. Often with WGS data, you will see well resolved 'peaks' with a fairly smooth curve of points leading up to it either side. From this we often infer that the highest point is the genetic variant of interest, and the other sites either side of that exhibit signals of selection because they reside close to, and thus are linked, to the variant of interest. However, consider that even in WGS data, unless we have every single genetic variant represented (which may not be the case, depending on our variant calling and filtering parameters) it is possible that the genetic variant of interest that we have identified is not the main one, but is simply another neighbouring linked SNP to one that is not represented in the data. This problem becomes even more relevant with reduced representation sequencing (RRS), for which the genome coverage may be extremely patchy<sup>C</sup>. Thus with all outlier analysis, but especially so for those using RRS data, remember that your flagged outliers are not exhaustive, and may themselves only be liked to the variant that is truely under selection.

<sup>C</sup> You may also want to consider linkage blocks.


https://onlinelibrary.wiley.com/doi/10.1111/mec.14549

https://pubmed.ncbi.nlm.nih.gov/29486732/

![ScreenShot](https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/asset/4be56b5b-8593-4116-ab9a-ec7b9e3c9a05/gr1.jpg)





## What data are we starting with

Metadata file, including individual and population names

covariates inc. environmental metadata.

The a genetic file. We will start with a VCF file.

Working with your own data


## PCAdapt

The first outlier program that we will deal with is PCAdapt. The reason I start here is simple: PCAdapt doesn't rely on any metadata. We simply require a data set of SNPs and using an ordination approach (principle component analysis) we identify outlier genetic sites. 

PCAdapt is an R package, therefore we will do all the following analysis in R. R is generally slower than unix command line programs as

<pre class="r"><code>library(dplyr)
library(maps)
library(mapdata)
library(ggplot2)
</code></pre>

You distirbution data will look different depending on where you have pulled it from. The important thing is that each data row is a distinct obervation of an individual, and that two of the supplied column information is the Latitiude and Longitude. In my data I also had a column that contained the country names of obersvation, which allowed me to subset my data to label native and invasive ranges with different colours.


## Funding 
<p align="center">

![ScreenShot](https://storage.corsizio.com/uploads/5cea29e798d9a757e03dba1c/events/6141a39502de2a7ff2f3f120/photo-_3UnQXVrb.jpg)

</p>

Thank you to the <a href="https://ausevo.com/ECR_grants_2022/">AES ERC Networking Grant Scheme</a> for funding this project.

## Project Contributors

https://www.dataschool.io/how-to-contribute-on-github/

