#### Setting up the stage #### 

set.seed(12345)

# Load packages 
library(readr) # To load csv files
library(nimble) # Run models
library(MCMCvis) # Visualize models
library(ggplot2) # For plots

#  Load samples
Sequoia_bears <- read_csv("./data/QCBears_fulldata.csv") # Full list of bear retained for parentage analysis
POPs <- read_csv("./data/BB_POpairs.csv") # Identified POPs from Sequoia

# Remove Drummond
Sequoia_bears <- subset(Sequoia_bears, Harvest_BMU != "Drummond" | is.na(Harvest_BMU))

# Rename columns
samples <- data.frame(ID = Sequoia_bears$Sample,
                      Sex = Sequoia_bears$Genetic_Sex,
                      BirthYear = Sequoia_bears$BirthYear,
                      SamplingYear = Sequoia_bears$Season_Year,
                      SamplingAge = Sequoia_bears$Age) # Age of the bear when sampled

POPs <- data.frame(bear1 = POPs$bear1, # Potential offspring
                   bear2 = POPs$bear2, # Potential parent
                   rel_type = POPs$rel_type) # Most likely relationship

# Add demographic parameters 

# Maximum age
male.max.age <- max(samples[samples$Sex == "XY", "SamplingAge"], # Oldest male in the sample
                    na.rm = TRUE) 
female.max.age <- max(samples[samples$Sex == "XX", "SamplingAge"], # Oldest female in the sample
                      na.rm = TRUE) 

# Age at maturity
male.first.repro <- 3 # Change to 4 here for alternative scenario
female.first.repro <- c(0.15, 0.47, 0.30, 0.08) # Probabilities of primiparity for 3-6 yo

####  Create pairwise comparisons #### 

#### |- All pairwise comparisons #### 

# Sort data frame by birth year, so that individual 1 is always older than individual 2
samples <- samples[order(samples$BirthYear),]

# Generate all the possible pairwise combinations
pair_df <- data.frame(t(combn(samples$ID, m=2)))

# Rename columns so they can easily be joined
colnames(pair_df) <- c("Ind1", "Ind2") 

# Bring relevant columns for Ind1 and Ind2
Ind1_data <- samples[,c("ID", "Sex", "BirthYear", "SamplingYear")]
Ind2_data <- samples[,c("ID", "Sex", "BirthYear", "SamplingYear")]

# Rename columns (make sure they are in the same order as the previous command)
colnames(Ind1_data) <- c("Ind1", "Ind1_sex", "Ind1_birth", "Ind1_sampling")
colnames(Ind2_data) <- c("Ind2", "Ind2_sex", "Ind2_birth", "Ind2_sampling")

# Combine the two data frames
pair_df_all <- merge(x = pair_df, y = Ind1_data, by = "Ind1")
pair_df_all <- merge(x = pair_df_all, y = Ind2_data, by = "Ind2")

# Reorder columns to make it cleaner (not a requirement)
pair_df_all <- pair_df_all[, c(2, 1, 3, 6, 4, 7, 5, 8)]

#### |-  Assign kinship #### 

# Create a unique identification for Ind1xInd2
pair_df_all$Ind12 <- paste0(pair_df_all$Ind1, "-", pair_df_all$Ind2)
POPs$Ind12 <- paste0(POPs$bear2, "-", POPs$bear1) # Note the reverse order

# Import kinship from POPs to the pairwise comparison data frame
pair_df_all$Kin <- ifelse(
  pair_df_all$Ind12 %in% POPs$Ind12 & pair_df_all$Ind1_sex == "XX", "MO", # If pair is found in POPs and the parent is female -> mother-offspring
  ifelse(pair_df_all$Ind12 %in% POPs$Ind12 & pair_df_all$Ind1_sex == "XY", "FO", # If pair is found in POPs and the parent is male -> father-offspring
         "UR")) # If pair is not found in POPs -> unrelated

fathers_RO <-subset(pair_df_all, Kin == "FO")
min(table(fathers_RO$Ind1))
max(table(fathers_RO$Ind1)) # Maximum number of sampled offspring for one father in the sample is 7

#### |-  Conditional expectation ####

# Assign age of first reproduction to females

# Get sex-specific dataframes 
pair_df_all_males <- pair_df_all[pair_df_all$Ind1_sex == "XY",]
pair_df_all_females <- pair_df_all[pair_df_all$Ind1_sex == "XX",]

# Create separate positive/negative data frame
pair_df_all_females_pos <- pair_df_all_females[pair_df_all_females$Kin == "MO",]
pair_df_all_females_neg <- pair_df_all_females[pair_df_all_females$Kin == "UR",]

# Females that reproduce at 3
# Randomly draw row indices from the entire set
rand_ind_female_3 <- sample(
  nrow(pair_df_all_females_neg),
  size = round(female.first.repro[1] * nrow(pair_df_all_females_neg)),
  replace = FALSE)  
# Extract the drawn rows from the data frame
pair_df_all_females_neg_3 <- pair_df_all_females_neg[rand_ind_female_3,] 
# Remove the drawn rows from the original data frame (so you don't select multiple time)
pair_df_all_females_neg_456 <- pair_df_all_females_neg[-rand_ind_female_3,] 

# Females that reproduce at 4
rand_ind_female_4 <- sample(
  nrow(pair_df_all_females_neg_456),
  size = round(female.first.repro[2] * nrow(pair_df_all_females_neg)),
  replace = FALSE) 
pair_df_all_females_neg_4 <- pair_df_all_females_neg_456[rand_ind_female_4,]
pair_df_all_females_neg_56 <- pair_df_all_females_neg_456[-rand_ind_female_4,]

# females that reproduce at 5
rand_ind_female_5 <- sample(
  nrow(pair_df_all_females_neg_56),
  size = round(female.first.repro[3] * nrow(pair_df_all_females_neg)),
  replace = FALSE) 
pair_df_all_females_neg_5 <- pair_df_all_females_neg_56[rand_ind_female_5, ]

# Females that reproduce at 6
# The number of pairwise comparisons left *should* correspond to the last probability
#in female.first.repro
pair_df_all_females_neg_6 <- pair_df_all_females_neg_56[-rand_ind_female_5,]

# Filter comparisons

# New filter dataframe
pair_df_filtered_males <- pair_df_all_males 

# Remove potential parents that were sampled before reaching age at maturity
pair_df_filtered_males <- subset(pair_df_filtered_males,
                                 Ind1_sampling - Ind1_birth >= male.first.repro)

# Remove potential parents that did not reach age at maturity when offspring was conceived
pair_df_filtered_males <- subset(pair_df_filtered_males,
                                 (Ind2_birth-1) - Ind1_birth >= male.first.repro)

# Remove potential parents that were sampled before offspring was conceived
# (because of lethal sampling)
# Remember that birth happens early in the year, and harvest is in autumn.
# Also remember that males actually reproduce the year prior to Ind2_birth (delayed implantation).
# So if Ind2_birth is January 2021, male can have reproduced in spring 2020 and be harvested in fall 2020.
pair_df_filtered_males <- subset(pair_df_filtered_males,
                                 Ind1_sampling >= (Ind2_birth-1))

# Negative comparisons
# For females, we are stricter than the males: if a female is harvested in fall the same year
# the offspring is born, we exclude the comparison because this is technically not possible to
# sample a cub-of-the-year
pair_df_all_females_neg_3 <- subset(pair_df_all_females_neg_3,
                                    Ind1_sampling - Ind1_birth >= 3 & 
                                      Ind2_birth - Ind1_birth >= 3 &
                                      Ind1_sampling > Ind2_birth)

pair_df_all_females_neg_4 <- subset(pair_df_all_females_neg_4,
                                    Ind1_sampling - Ind1_birth >= 4 & 
                                      Ind2_birth - Ind1_birth >= 4 &
                                      Ind1_sampling > Ind2_birth) 

pair_df_all_females_neg_5 <- subset(pair_df_all_females_neg_5,
                                    Ind1_sampling - Ind1_birth >= 5 & 
                                      Ind2_birth - Ind1_birth >= 5 &
                                      Ind1_sampling > Ind2_birth)

pair_df_all_females_neg_6 <- subset(pair_df_all_females_neg_6,
                                    Ind1_sampling - Ind1_birth >= 6 & 
                                      Ind2_birth - Ind1_birth >= 6 &
                                      Ind1_sampling > Ind2_birth)

# Positive comparisons
# For MO, we must take the least-strict age of primiparity (i.e., 3yo)
pair_df_all_females_pos <-  subset(pair_df_all_females_pos,
                                   Ind1_sampling - Ind1_birth >= 3 & 
                                     Ind2_birth - Ind1_birth >= 3 &
                                     Ind1_sampling > Ind2_birth)

# Then group all possible reproductive ages and positives and negatives comparisons
# in a unique data frame
pair_df_filtered_females <- rbind.data.frame(
  pair_df_all_females_pos,
  pair_df_all_females_neg_3,
  pair_df_all_females_neg_4,
  pair_df_all_females_neg_5,
  pair_df_all_females_neg_6)

#### |- Group comparisons ####

# Extract positive father-offspring pairs
pop_dad_positives <- pair_df_filtered_males[pair_df_filtered_males$Kin == "FO",] 

# Extract unrelated individuals
pop_dad_negatives <- pair_df_filtered_males[pair_df_filtered_males$Kin == "UR",] 

# Group all unique combinations of Ind1 birth and Ind2 birth, and count
pop_dad_positives <- aggregate(pop_dad_positives$Ind1_birth,
                               by = list(pop_dad_positives$Ind1_birth,
                                         pop_dad_positives$Ind2_birth),
                               FUN = length)

pop_dad_negatives <- aggregate(pop_dad_negatives$Ind1_birth,
                               by = list(pop_dad_negatives$Ind1_birth,
                                         pop_dad_negatives$Ind2_birth),
                               FUN = length)

# Rename the columns for clarity 
colnames(pop_dad_positives) <- c("Ind1_birth", "Ind2_birth", "yes") 
colnames(pop_dad_negatives) <- c("Ind1_birth", "Ind2_birth", "no")

# Merge positives and negatives
pop_dad_comps <- merge(pop_dad_positives, pop_dad_negatives, all = TRUE)
pop_dad_comps[is.na(pop_dad_comps)] <- 0
pop_dad_comps$all <- pop_dad_comps$yes + pop_dad_comps$no

# Calculate age of Ind1 (father) at year of birth Ind2 -1
pop_dad_comps$Ind1_age_b2min1 <- (pop_dad_comps$Ind2_birth - 1) -  pop_dad_comps$Ind1_birth

# Extract positive mother-offspring pairs
pop_mom_positives <- pair_df_filtered_females[pair_df_filtered_females$Kin == "MO",]

# Extract unrelated individuals
pop_mom_negatives <- pair_df_filtered_females[pair_df_filtered_females$Kin == "UR",] 

# Group all unique combinations of Ind1 birth and Ind2 birth
pop_mom_positives <- aggregate(pop_mom_positives$Ind1_birth,
                               by = list(pop_mom_positives$Ind1_birth,
                                         pop_mom_positives$Ind2_birth), 
                               FUN = length)

pop_mom_negatives <- aggregate(pop_mom_negatives$Ind1_birth,
                               by = list(pop_mom_negatives$Ind1_birth,
                                         pop_mom_negatives$Ind2_birth),
                               FUN = length)

# Rename the columns for clarity 
colnames(pop_mom_positives) <- c("Ind1_birth", "Ind2_birth", "yes") 
colnames(pop_mom_negatives) <- c("Ind1_birth", "Ind2_birth", "no")

# Merge positives and negatives
pop_mom_comps <- merge(pop_mom_positives, pop_mom_negatives, all = TRUE)
pop_mom_comps[is.na(pop_mom_comps)] <- 0
pop_mom_comps$all <- pop_mom_comps$yes + pop_mom_comps$no

# Calculate age of Ind1 (mother) at year of birth Ind2
pop_mom_comps$Ind1_age_b2 <- pop_mom_comps$Ind2_birth -  pop_mom_comps$Ind1_birth

#### Auxiliary data #### 

# Get sampled adult males
samples.male.adult <- samples[samples$Sex == "XY" &
                                samples$SamplingAge >= male.first.repro,]

# Summary of age distribution
summary.ages.male <- aggregate(samples.male.adult$SamplingAge,
                               by = list(samples.male.adult$SamplingAge),
                               FUN = length)

colnames(summary.ages.male) <- c("age", "N") 

# Fill all possible ages, from age at maturity to maximum age
# Including ages that may not be present in the samples
ages <- data.frame(age = c(male.first.repro:male.max.age)) 
Nages.male <- merge(ages, summary.ages.male, by = "age", all = TRUE)
Nages.male[is.na(Nages.male)] <- 0 # Change NAs to 0's

# Get proportions for each age
Nages.male$prop <- Nages.male$N / sum(Nages.male$N) 

# Add juveniles (but do not count them)
if(male.first.repro == 3) {Nages.male <- rbind(c(1, 0, 0),
                                               c(2, 0, 0),
                                               Nages.male)}

if(male.first.repro == 4) {Nages.male <- rbind(c(1, 0, 0),
                                               c(2, 0, 0),
                                               c(3, 0, 0),
                                               Nages.male)}

age_dist_male <- Nages.male

# Raw values extracted directly from Figure 1 in Costello et al., 2009
costello.age.RO.raw <- c(0, 0, 0.044, 0.065, 0.089, # 1-5
                         0.115, 0.141, 0.162, 0.181, 0.191, # 6-10
                         0.193, 0.187, 0.175, 0.154, 0.132, # 11-15
                         0.105, 0.080, 0.057, # 16-18
                         rep(0.01, times = (male.max.age - 18)))

# Relative values so that sum(costello.age.RO.raw) from 3 to 18yo == 1
costello.age.RO.rel <- c(0, 0, 0.021, 0.031, 0.043,
                         0.056, 0.068, 0.078, 0.087, 0.092,
                         0.093, 0.090, 0.085, 0.074, 0.064,
                         0.051, 0.039, 0.028,
                         rep(0.01, times = (male.max.age - 18)))

age_dist_male$RO <- costello.age.RO.rel[age_dist_male$age]

costello.age.RO <- costello.age.RO.rel

# For the alternative scenario with no age-dependent fecundity, run this
costello.age.RO <- rep(1, times = male.max.age)
age_dist_male$RO <- costello.age.RO[age_dist_male$age]
age_dist_male

#### Create arrays #### 

# All unique b2-1
# Again, a cub born in 2022 relates to male reproductive output from 2021
dad_rownames <- unique(pop_dad_comps$Ind2_birth)-1
# Every age possible
dad_colnames <- c(1:male.max.age) # technically should be 1:male.max.age-1 but it does not matter

# Create empty arrays with the maximum number of dimensions for 
# birth year of offspring -1 and parent's age
FO_array <- pop_dad_all_array <-
  matrix(data = ,
         nrow = length(unique(pop_dad_comps$Ind2_birth)),
         ncol = male.max.age,
         dimnames = list(dad_rownames, dad_colnames))

# Populate the matrices with "yes" (FO) and
#"all" (positive + negative pairwise comparisons)
for (i in (min(unique(pop_dad_comps$Ind2_birth))-1):(max(unique(pop_dad_comps$Ind2_birth))-1)){
  
  # Father-Offspring pairs
  FO_vec <- pop_dad_comps[pop_dad_comps$Ind2_birth == i,c("yes","all","Ind1_age_b2min1")]
  
  # if there is a gap in Ind2_birth in the samples, skip to the next Ind2_birth
  if(nrow(FO_vec)==0){next}
  
  # Populate the arrays for each b2xage
  FO_array[paste0(i), paste0(c(FO_vec$Ind1_age_b2min1))] <- FO_vec$yes 
  pop_dad_all_array[paste0(i), paste0(c(FO_vec$Ind1_age_b2min1))] <- FO_vec$all
}

# Change the NA's into 0's
FO_array[is.na(FO_array)] <- 0
pop_dad_all_array[is.na(pop_dad_all_array)] <- 0

# All unique b2
mom_rownames <- unique(pop_mom_comps$Ind2_birth)
# Every age possible
mom_colnames <- c(1:female.max.age) 

# Create empty array with the maximum number of dimensions for 
# birth year of sibling and parent's age
MO_array <- pop_mom_all_array <-
  matrix(data = , 
         nrow = length(unique(pop_mom_comps$Ind2_birth)),
         ncol = female.max.age,
         dimnames = list(mom_rownames, mom_colnames))

# Populate the matrices with yes (MO) and 
# all (all pairwise comparisons)
for (i in min(unique(pop_mom_comps$Ind2_birth)): max(unique(pop_mom_comps$Ind2_birth))){
  MO_vec <- pop_mom_comps[pop_mom_comps$Ind2_birth == i,c("yes","all","Ind1_age_b2")]
  if(nrow(MO_vec)==0){next} 
  MO_array[paste0(i), paste0(c(MO_vec$Ind1_age_b2))] <- MO_vec$yes
  pop_mom_all_array[paste0(i), paste0(c(MO_vec$Ind1_age_b2))] <- MO_vec$all
}

# Change the NA's into 0's
MO_array[is.na(MO_array)] <- 0 
pop_mom_all_array[is.na(pop_mom_all_array)] <- 0

#### Run CKMR #### 

# We could have males and females in the same nimbleCode
# However changing the male model parameters ever-so-slightly affects the results for females
# because we set a seed for the once at the beginning and then change some of the operations (I think)

#### |- Females ####

POP_model_female <- nimbleCode({
  
  # Priors
  Nf ~ dunif(1, 10000) # Uniform prior for abundance
  rf ~ dunif(-0.5, 0.5) # Uniform prior for population growth

  # Kinship probabilities
    for (b2 in 1:pop_mom_length_b2){ # For each birth year of offspring
    
    for (age in 3:female_max_age) { # For each age of the parent
      
      MO[b2,age] ~ dbinom(
        1/(Nf * exp(rf * (pop_mom_offspring_birthyear[b2] - estimation_year))),
        pop_mom_all_comps[b2,age]) 
      
    } # age 
  } # b2
  
})

# Add constants
my.constants <- list(
  pop_mom_length_b2 = nrow(MO_array),   # Number of cohort comparisons to loop over
  pop_mom_offspring_birthyear = as.numeric(rownames(MO_array)),   # Detailed cohort comparisons to loop over
  female_max_age = female.max.age,
  estimation_year = 2020 # The year to which estimation of N is anchored
)

# Add data
my.data <- list(
  MO = MO_array, # Verified POP
  pop_mom_all_comps = pop_mom_all_array # All comparisons
)

# Parameters
initial.values <- function() list(Nf = rnorm(1, mean = 5000, sd = 1000), 
                                  rf = rnorm(1, mean = 0, sd = 0.05))

parameters.to.save <- c("Nf", "rf")

n.iter <- 100000
n.burnin <- 25000
n.chains <- 4
n.thin <- 10

# Run the model
mcmc.output_female <- nimbleMCMC(code = POP_model_female,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          thin = n.thin,
                          nburnin = n.burnin,
                          nchains = n.chains,
                          progressBar = FALSE)

# Visualize results
MCMCsummary_female <- MCMCsummary(mcmc.output_female, round = 3)
MCMCsummary_female

MCMCtrace(mcmc.output_female, pdf = FALSE)

# Extract the mode (here defined as the maximum density of the posterior distribution, not the actual mode per se)
Nf_output <- c(mcmc.output_female[[1]][,"Nf"],
               mcmc.output_female[[2]][,"Nf"],
               mcmc.output_female[[3]][,"Nf"],
               mcmc.output_female[[4]][,"Nf"])

rf_output <- c(mcmc.output_female[[1]][,"rf"],
               mcmc.output_female[[2]][,"rf"],
               mcmc.output_female[[3]][,"rf"],
               mcmc.output_female[[4]][,"rf"])

max_Nf <- which.max(density(Nf_output)$y)
Nf <- round(density(Nf_output)$x[max_Nf], digits = 0)
Nf

max_rf <- which.max(density(rf_output)$y)
rf <- round(density(rf_output)$x[max_rf], digits = 3)
rf

#### |- Males ####

POP_model_male <- nimbleCode({
  
  # Priors
  Nm ~ dunif(1, 10000)
  rm ~ dunif(-0.5, 0.5)
  
  # Kinship probabilities
  for (b2 in 1:pop_dad_length_b2){ # For each **conception** year of offspring (birth year -1)
    
    for (age in male_first_repro:male_max_age) { # For each age of the parent
      
      FO[b2,age] ~ dbinom(
        RO[age] / sum(TRO[b2,male_first_repro:male_max_age]),
        pop_dad_all_comps[b2,age]) 
      
    } # age
  } # b2
  
  # Separate loop to derive parameter TRO for each year and age of the father
  for (b in 1:pop_dad_length_b2){ # For each **conception** year of offspring (birth year -1)
    
    for (age in male_first_repro:male_max_age){ # For each age of the parent
      
      # Age-weighted TRO for adults.
      TRO[b,age] <- (RO[age] * age_dist_male[age] *
                       Nm * exp(rm * (pop_dad_offspring_birthyear[b] - estimation_year)))
      
    } # age
  } #b
  
  # Derived parameters for ages 1 and 2, in all years, is fixed at 0
  # (potential parents are juveniles, and so their RO is 0)
  TRO[1:pop_dad_length_b2,1:(male_first_repro-1)] <- 0
})

my.constants <- list(
  pop_dad_length_b2 = nrow(FO_array),
  pop_dad_offspring_birthyear = as.numeric(rownames(FO_array)), 
  RO = costello.age.RO, # Reproductive output at age
  age_dist_male = age_dist_male$prop, # Adult age distribution
  male_first_repro = male.first.repro,
  male_max_age = male.max.age,
  estimation_year = 2020 
)

my.data <- list(
  FO = FO_array, # Verified POP
  pop_dad_all_comps = pop_dad_all_array # All comparisons
)

# Parameters
initial.values <- function() list(Nm = rnorm(1, mean = 5000, sd = 1000), 
                                  rm = rnorm(1, mean = 0, sd = 0.05))

parameters.to.save <- c("Nm", "rm")

n.iter <- 100000
n.burnin <- 25000
n.chains <- 4
n.thin <- 10

# Run the model
mcmc.output_male <- nimbleMCMC(code = POP_model_male,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          thin = n.thin,
                          nburnin = n.burnin,
                          nchains = n.chains,
                          progressBar = FALSE)

# Visualize results
MCMCsummary_male <- MCMCsummary(mcmc.output_male, round = 3)
MCMCsummary_male

MCMCtrace(mcmc.output_male, pdf = FALSE)

# Extract the mode 
Nm_output <- c(mcmc.output_male[[1]][,"Nm"],
               mcmc.output_male[[2]][,"Nm"],
               mcmc.output_male[[3]][,"Nm"],
               mcmc.output_male[[4]][,"Nm"])


rm_output <- c(mcmc.output_male[[1]][,"rm"],
               mcmc.output_male[[2]][,"rm"],
               mcmc.output_male[[3]][,"rm"],
               mcmc.output_male[[4]][,"rm"])

max_Nm <- which.max(density(Nm_output)$y)
Nm <- round(density(Nm_output)$x[max_Nm], digits = 0)
Nm

max_rm <- which.max(density(rm_output)$y)
rm <- round(density(rm_output)$x[max_rm], digits = 3)
rm