# *********************************************************
#
#   Introduction to Phylogenetic Comparative Methods
#
#   Emma Dunne (emma.dunne@fau.de)
#   Late updated: August 2024
#
# _________________________________________________________
#
#   2. Phylogenetic signal 
# 
# *********************************************************

# ******************************************************
#
#   Pterosaur climate niche evolution
#
#   MSc thesis 2024
#
#   Lead: Carolin Grether (carolin.grether@fau.de)
#   Supervisor: Emma Dunne (emma.dunne@fau.de)
# ______________________________________________________
#
#   Phylogenetic signal in flight efficiency data (Venditti et al. 2020)
# 
# ******************************************************

# Prepare for analysis

library(phytools)
library(tidyverse)
library(strap)
library(viridis)
library(ape)

## Test the phlogenetic signal (Pagel's λ and Blomberg's K)
lambdaEff <- phylosig(ptero_tree_eff, Inv_CoT_matrix, method = "lambda", test = TRUE)
lambdaEff 
 # P-Value >0.001 and Lambda 0.9999

KEff <- phylosig(ptero_tree_eff, Inv_CoT_matrix, method = "K", test = TRUE, nsim = 10000)
KEff

## Phylogenetic signal is the pattern where close relatives have more similar 
##    trait values than more distant relatives. 
## There are two popular methods for estimating phylogenetic signal:
##      1. Pagel’s λ (lambda)
##      2. Blomberg’s K

## First, let's check the phylogenetic signal of frog eye size using Pagel’s λ (lambda)

## Estimate Pagel’s λ (lambda) using the function phylosig() in the package 'phytools'
lambdaEye <- phylosig(mytree, logEye, method = "lambda", test = TRUE)
lambdaEye # take a look at the output

## Interpreting the output:
##    P-value is the p value from a likelihood ratio test testing whether  λ 
##    is significantly different from 0 (no phylogenetic signal). 
##    Here λ = 0.814 and P < 0.001. We can interpret this as  λ being significantly 
##    different from 0, i.e. there is significant phylogenetic signal in log eye size


## Next, the same but with Blomberg’s K, also using the phytools package:
KEye <- phylosig(mytree, logEye, method = "K", test = TRUE, nsim = 1000)
## Additionally we add the argument nsim = 1000. This is because we need to use a 
##    randomisation test to determine whether K is significantly different from 0
KEye

## Interpreting the output:
##    The observed value of K is  compared to the randomized values. The p value
##    tells us how many times out of 1000, a randomised value of K is more extreme
##    than the observed value. If this number is low, the p value is low (e.g. 
##    if 5 out of 1000 randomised values of K are more extreme than the observed
##    value p = 5/1000 = 0.005)

##____________________________________________
## Bonus exercise:
##  - What is λ and K for snout vent length?
##____________________________________________