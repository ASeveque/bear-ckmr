# Data and code for "Application of close-kin mark-recapture to an American black bear population using harvest samples".

See [Maya Weissman's repository](https://github.com/mweissman97/MI_bear_pedigree) for the pipeline to analyze nuclear and mitochondrial SNP data for Upper Michigan Black Bears in order to determine Parent-Offspring pairs.

## Abstract
Close-kin mark-recapture (CKMR) offers a promising alternative to traditional mark-recapture methods for estimating wildlife population parameters. CKMR estimates demographic parameters by identifying genetic relationships between individuals, effectively "recapturing" parental genotypes through their offspring’s genetic profiles. One key advantage of CKMR over traditional methods is the ability to estimate abundance without recapture, facilitating the use of lethally collected biological samples. Here, we present an application of CKMR to a terrestrial species using samples collected from harvested individuals. We analyzed genetic data from 1,767 American black bears (Ursus americanus) harvested over two consecutive years in Michigan's Upper Peninsula. Using Genotyping-in-Thousands by sequencing (GT-seq), we identified parent-offspring pairs and developed sex-specific models tailored to black bear reproductive biology and management strategies. Female population estimates aligned with our expectations, based on prior knowledge of the population. In contrast, male population estimates deviated from expectations, potentially due to an unmodeled combination of factors pertaining to male reproductive biology and sampling bias. Our study demonstrates the viability of using harvest-based CKMR for monitoring of terrestrial wildlife populations, but also highlights the need for careful model specification. By providing a novel method to estimate demographic parameters and population size, CKMR using harvest samples can inform management decisions and harvest quotas. This approach leverages existing management practices involving collection of genetic samples, such as mandatory harvest registrations and health monitoring programs. Our findings contribute to the growing body of CKMR applications and expand its use to terrestrial species management.
The manuscript is currently under review. 

## Input data files

* `QCBears_fulldata.csv`: Final set of samples used in the CKMR analysis after quality filtering
  * Sequoia_ID - shortened sample ID (e.g. BB21_1000)
  * Sample - long sample IDs (e.g. BB21_1000_MTU20Oct23_R1)
  * NW_025576331.1_3826608:NW_025578505.1_8669226 - SNP data, where 0 = no copies of variant, 1 = heterozygous, 2 = homozygous with variant, -9 = missing read.
  * Uam_SEXY1:Uam_sry2 - genotypes for sex markers; XX = female, XY = male, 0 = missing
  * Genetic_Sex - most common sex from sex marker columns (Uam_SEXY1:Uam_sry2)
  * sex_agreement_prop - fraction of sex markers (Uam_SEXY1:Uam_sry2) that agree with Genetic_Sex consensus
  * Season_Year - year sample was harvested
  * Species - common name for species (i.e. "Black Bear")
  * Harvest_Date - date sample was harvested (e.g. 8-Sep-21)
  * How_Taken - how sample was harvested (e.g. "Hunting")
  * Harvest_BMU - Bear Management Unit where sample was harvested
  * Harvest_County - County where sample was harvested
  * Registration_Sex - sex identified during registration, may differ from genetic sex
  * BirthYear - estimated year the bear was born based on sample age
  * Latitude - approximate decimal latitude coordinates where sample was harvested
  * Longitude - approximate decimal longitude coordinates where sample was harvested
  * filter_remove - whether the sample passed quality control filters (i.e. "pass" or "fail")  
  
* `BB_POpairs.csv`: Parent-offspring pairs identified by Sequoia (full pedigree + age priors)
  * bear1 - long sample ID corresponding to the younger bear (aka offspring) in the relationship pair (e.g. BB21_1000_MTU20Oct23_R1)
  * bear2 - long sample ID corresponding to the older bear (aka parent) in the relationship pair (e.g. BB21_1000_MTU20Oct23_R1)
  * LLR - log likelihood ratio, or log10 transformed likelihood the pair have the assigned relationship divided by the likehlihood the pair have the next most likely relationship type; higher values of LLR correspond to higher confidence
  * rel_type - relationship type; PO_M = mother-offspring, PO_P = father-offspring
  * pair_id - concatenates bear1 ID and bear2 ID to create a unique identifier for the pair
  * bear1_birthyear - birth year of bear1
  * bear2_birthyear - birth year of bear2
  * age_diff - age of parent at offspring's birth, bear1_birthyear - bear2_birthyear
  * bear1_sex - sex of bear1
  * bear2_sex - sex of bear2
  * bear1_county - harvest county of bear1
  * bear2_county - harvest county of bear2
  * county_difference - spatial relationship between bear1_county and bear2_county; "same county" = harvested in the same county, "neighboring county" = bear1_county borders bear2_county, "further" = bears were not harvested in the same or bordering counties
  * module - sequoia model used; "full ped" = full pedigree model, rather than parent-offspring only pedigree
  * age_prior - whether sequoia model incorporated age priors; "yes ap" = age priors were used
  * LLR_bin - turns LLR into a discrete bin; "negative" = LLR was less than 0 and thus there is low confidence in relationship assignment, ">0" = LLR was greater than 0 and thus there is high confidence in relationship assignment

## R script

* `CKMR_Michigan.R`: R code for running the analysis

