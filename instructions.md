#### Please upload results from your experiment's previous waves as a .csv file in the following format:
**Previous data**: Column labels
* *outcome*:  A binary outcome variable
* *treatment1, treatment2, ...*: Dummy variables for each of the treatments,    
or *treatment*: Categorical variable, with values corresponding to the treatments
* *covar1, covar2, ...*: Covariates (strata will be created from full interactions of controls)  
The app allows for the case that no covariates are provided.

**Replicate draws**: The number of MCMC draws from the posterior. A larger choice yields more precise assignment probabilities, but results in slower computation.  

**Share fully randomized**: Choose the share of observations to be assigned with full randomization, rather than with Thompson sampling.