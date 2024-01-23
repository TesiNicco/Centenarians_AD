# Centenarians_AD
## General information
Script accompanying the manuscript **Cognitively Healthy Centenarians are genetically protected against Alzheimer's Disease**.  
A preprint of the manuscript is available on [medRxiv](https://www.medrxiv.org/content/10.1101/2023.05.16.23290049v1).  
The manuscript is currently under revision at a peer-reviewer journal. A link to the publication will be available upon manuscript acceptance.  
The script is written in **R**.

## What the script does
The script performs all association and downstream analyses presented in the manuscript. The different analyses as presented in the manuscript are:
- single-variant analysis: we compared the frequency of SNP associated with Alzheimer's Disease (AD) between *(i)* AD cases, *(ii)* age-matched controls, and *(iii)* cognitively healthy centenarians
- calculation of PRS and comparison of the PRS between *(i)* AD cases, *(ii)* age-matched controls, and *(iii)* cognitively healthy centenarians
- power analysis, where we compare the statistical power given by the cognitively healthy centenarians as compared to the age-matched controls
- main plots and supplementary plots
- main tables and supplementary tables
- additional analysis arisen during the peer-review process

## How to run the script
To run the script, additional files are needed on top of those provided in the repository:
- the individual genotypes
- the quality control summary done prior to analysis. This QC aims at excluding individuals of non-European ancestry as well as individuals with a family relationship
- the phenotype file of the individuals
*Even when all input files are available, we enourage to execute the script interactively in a R session*

## Provided input files
- *AD_snps_clinicalAD.txt*: A plain-text table containing the summary statistics of the SNP associated with Alzheimer's Disease. The associations are relative to the comparison of clinical AD patients as opposed to healthy age-matched individuals
- *AD_snps_proxyAD.txt*: A plain-text table containing the summary statistics of the SNP associated with Alzheimer's Disease. The associations are relative to the comparison of clinical AD + proxy AD as opposed to healthy age-matched individuals.
Summary statistics are extracted from [Bellenguez et al., 2022](https://www.nature.com/articles/s41588-022-01024-z) manuscript. 

## Glossary
*SNP*: Single Nucleotide Polymorphisms. A type of genetic variation characterized by a single change in the nucleotide sequence  
*AD cases*: Patients with diagnosis of Alzheimer's Disease  
*Age-matched controls*: Cognitively healthy individuals, age-matched with the AD patients  
*Cognitively Healthy Centenarians*: Individuals older than 100 years of age with maintained cognitive health. Individuals are from the [100-plus Study](https://holstegelab.eu/)  
*PRS*: Polygenic Risk Score. A numerical score that quantifies the individual genetic risk to develop a given trait  

## Contact
For further information, questions, or comments, please drop and email to the [the first author](mailto:n.tesi@amsterdamumc.nl).  

