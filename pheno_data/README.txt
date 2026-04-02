* The phenotypic data include three agronomic traits: 
lodging: rated on a 1-9 scale (with a higher score indicating severe lodging)
plant height: measured from the soil surface to the tip of spike in cm (excluding awns)
heading date (recorded in days from January 1st for winter type and from the sowing date onward for the spring type in days)
and four disease traits: Puccinia hordei, Blumeria graminis hordei, Ramularia collo-cygni, and Rhychnosporium commune: all scored on an ordinal scale from 1 (fully resistant) to 9 (fully susceptible).

* The ‘input_raw_data_spring.tsv’ contains the input raw data of the spring population for the best linear unbiased estimations (BLUEs) and heritability calculation.

* The ‘input_raw_data_winter.tsv’ contains the input raw data of the winter population for the best linear unbiased estimations (BLUEs) and heritability calculation.

* The ‘BLUEs_across_env.tsv’ contains the best linear unbiased estimations (BLUEs) for all traits of the whole population with their respective growth habit and status (Elite/PGR) information, as well as the doi number for each plant genetic resources.
Missing value in the BLUEs file is due to the missing of raw data and the excluding of outliers.

* The R code sample_code.R performs the best linear unbiased estimations (BLUEs) and heritability calculation. 
Slight adjustment has to be done for the model inside the "BLUE" function to analysis other traits if the model cannot converge, and a license for asreml-R is needed.

* The ‘varcomp_spring.tsv’ contains the corresponding heritability and relavent variance component.
Sigma.Gen: genotypic variance.
Sigma_Gen.x.Env: the interaction between genotype and environment.
Sigma_Env: environmental variance.
Sigma_E: the residual variance.
h2.Gen: heritability.
Nr.Rep: the number of replications per genotype.
Nr.Env: the average number of environments in which the genotypes were evaluated.
the abbreviations are the same as those in ‘varcomp_winter.tsv’.