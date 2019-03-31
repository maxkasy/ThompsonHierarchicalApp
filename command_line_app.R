##############################################
#Invoke this script from the terminal as follows:
#RScript --vanilla command_line_app.R priordata_test_missings.csv
##############################################

source("ReadDataApp.R")
source("MCMC_HierarchicalThompson.R")

datapath=commandArgs(trailingOnly = T)[1]

RR=50000
alpha=0.2
nx=16
key=tibble(Strata=factor(1:nx),
           strata=factor(1:nx))

priordata=ReadDataApp(datapath,key)

Pstar=DtchoiceMCMCProbabilities(priordata$Y,priordata$D,priordata$X, #outcomes, treatments, and covariates thus far
                                priordata$k,priordata$nx, #number of treatments and number of strata
                                RR=RR)
Pactual=(1-alpha) * Pstar + alpha * (1/priordata$k)

filename = paste(Sys.Date(), "_treatmentprobabilities.csv", sep="")
write_csv(Pactual, path=filename)
