# Ev1_SelectionMetaAnalysis

The following vignette will take you through the basics of genetic selection analysis. 

The project has been funded by the <a href="https://ausevo.com/ECR_grants_2022/">AES ERC Networking Grant Scheme</a>.

<h2><i>Genetic Outlier and Association Analysis</i></h2>

There are lots of interesting patterns that you can extract from geneti marker data. This can include patterns of linkage, ballancing selection, and another example. One of the most common ones is to try find sites on the genome that are under divergent selection.

When we look through a genome to try find loci that are under divergent selection, we often conduct what is called outlier or association analyses. **Outlier analysis** requires just knowledge of the genetics of your samples (along with the usual sample metadata, for example population grouping), and tries to find loci that behave very differently from the rest of the genome (with the assumption being that the rest of the genome represents patterns of neutral genetic diveristy)<sup>A</sup>. Meanwhile **association analysis** require some sort of covariate data, and tests whether there are any genetic variants statistically associated with this new data (you may have heard the term [GWAS](https://www.genome.gov/genetics-glossary/Genome-Wide-Association-Studies). This data can come in the form of phenotype data (e.g. morphology, disease status, physiology measures), or could be spatial (e.g. environmental, spatial). Association tests will look for sites in the genome where the precense or absence of a variant is highly correlated with the values in the co-variate data, usually through some regression type analysis.

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
> Completing outlier analysis 



https://onlinelibrary.wiley.com/doi/10.1111/mec.14549



![ScreenShot](https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/asset/4be56b5b-8593-4116-ab9a-ec7b9e3c9a05/gr1.jpg)


## Where to find data

Great - you're working on a new species and need to make a plot of where it is located for some figure or presentation. It's actually a very quick proccess. First step is to find a data base where you can pull information about the distribution of your species. There are many options here, some will be quite taxa or country specific. Some good starting points are:

<ul>
<li><a href="https://ebird.org/home">eBird</a></li>
<li><a href="GBIF">https://www.gbif.org/</a></li>
<li><a href="https://www.inaturalist.org/">iNaturalist</a></li>
<li><a href="https://mol.org/">Map of Life</a></li>
</ul>

Find you species and retrieve you data.

## Plotting in R

Packages used in the below code.

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

