# Data and code for "Application of close-kin mark-recapture to an American black bear population using harvest samples".

See repository mweissman97/MI_bear_pedigree for the pipeline to analyze nuclear and mitochondrial SNP data for Upper Michigan Black Bears in order to determine Parent-Offspring pairs.

* **data:**
  
  * `QCBears_fulldata.csv`: Final set of samples used in the CKMR analysis after quality filtering

  * `BB_POpairs.csv`: Parent-offspring pairs identified by Sequoia (full pedigree + age priors)
  
* **R:**

  * `CKMR_Michigan.R`: R code for running the analysis
  

* **Abstract:** Close-kin mark-recapture (CKMR) offers a promising alternative to traditional mark-recapture methods for estimating wildlife population parameters. CKMR estimates demographic parameters by identifying genetic relationships between individuals, effectively "recapturing" parental genotypes through their offspring’s genetic profiles. One key advantage of CKMR over traditional methods is the ability to estimate abundance without physical recapture events, facilitating the use of lethally collected biological samples. Here, we present the first application of CKMR to a terrestrial species using samples collected from harvested individuals. We analyzed genetic data from 1,767 American black bears (Ursus americanus) harvested over two consecutive years in Michigan's Upper Peninsula. Using Genotyping-in-Thousands by sequencing (GT-seq), we identified parent-offspring pairs and developed sex-specific models tailored to black bear reproductive biology and management strategies. Female population estimates aligned with our expectations, based on prior knowledge of the population. In contrast, male population estimates deviated from expectations, potentially due to an unmodeled combination of factors pertaining to male reproductive biology and sampling bias. Our study demonstrates the viability of using harvest-based CKMR for monitoring of terrestrial wildlife populations, but also highlights the need for careful model specification. By providing a novel method to estimate demographic parameters and population size, CKMR using harvest samples can inform management decisions and harvest quotas. This approach leverages existing management practices involving collection of genetic samples, such as mandatory harvest registrations and health monitoring programs. Our findings contribute to the growing body of CKMR applications and expand its use to terrestrial species management.

The manuscript is currently under review. 

