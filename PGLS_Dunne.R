# *********************************************************
#
#   Introduction to Phylogenetic Comparative Methods
#
#   Emma Dunne (emma.dunne@fau.de)
#   Late updated: August 2024
#
# _________________________________________________________
#
#   3. Phylogenetic Generalised Least Squares (PGLS)
# 
# *********************************************************

## Let’s investigate the relationship between eye size in frogs 
##    and their body size. Theory leads us to believe that bigger 
##    animals will have bigger eyes, but is this true in frogs?

## Much of the information below assumes that you are familiar 
##    with fitting and interpreting linear models (regressions, 
##    ANOVA, ANCOVA) in R using the function lm.

## Let's look at the phylogenetic non-independence/pseudoreplication 
##    on the graph by colouring the points by family.
ggplot(mydata, aes(x = log(mass), 
                   y = log(eyesize), 
                   colour = Family)) +
  geom_point() +
  theme_bw()

## ___________________________________________________________
## Interpreting the output:
##    - what patterns do you notice about different families?
##    - is there phylogenetic independence evident?
## ___________________________________________________________



# PGLS --------------------------------------------------------------------

## PGLS is one method that accounts for phylogenetic non-independence 
## Another popular earlier method is independent contrasts (PIC). This method 
##    is really similar to PGLS, in fact it is just a special kind of PGLS 
##    where λ is equal to 1.

## There are multiple functions in several packages to fit a PGLS. 
## Here we will use caper, including its built-in function comparative.data()

mydata2 <- rename(mydata, tiplabel = Binomial)
head(mydata2)

model_PGLS <- pgls(log(eyesize) ~ log(mass), 
                   comparative.data(mytree, mydata2, "tiplabel"), 
                   lambda="ML")


# Check model diagnostic plots --------------------------------------------

## You must always check model diagnostic plots whenever you fit a model in 
##    R to check that your data meet the assumptions of the model. f your 
##    assumptions are not met, then the model is useless and needs to be 
##    modified before you can use it to extract p values etc.


## Set up a plotting window with four panes (there are four diagnostic plots):
par(mfrow = c(2, 2))
# Now plot the model diagnostics
plot(model_PGLS)
# Return the plot window to one pane for plots later
par(mfrow = c(1, 1))

## To interpret these plots, you might find the following resource helpful:
##    http://www.r4all.org/
## It takes practice to know what is “good”, “bad” and “acceptable” with these plots

## 1. Top left: 
##    We should not see any data with a studentized residual >± 3. Any species with such
##    large residuals may overly influence the results of the regression. Often these 
##    are the result of measurement error associated with species pairs joined by very short branches.

## 2. Top right: Normal Q-Q plot
##    This tests for normality of residuals, one of the assumptions of linear models.
##    Points should approximately fall on the line - often look a bit messy

## 3. Bottom left: Fitted versus Residuals plot
##    We should see a fairly random scattering of points often referred to as a 'sky at night distribution'
##    This diagnostic plot is related to the systematic component of the model - any pattern here
##    (e.g. strong positive or negative slopes, or humped or U-shaped patterns) suggests that the model 
##    has not been correctly specified. This might reflect a missing variable or interaction term, or that 
##    you need to transform your variables. Or that you need a different kind of model entirely.

## 4. Bottom tight: Fitted versus Observed plot
##    We should see a correlation among the observed values and fitted values in the model. 
##    This will vary depending on how good the model fit is.


# PGLS output and interpretation ------------------------------------------

## Now let's look at our PGLS output and interpret the results

## Peep at the model outputs using anova()
## This uses sequential sum of squares to tell you whether a model including 
##    your predictor variable(s) is a better fit than a model without your 
##    predictor variable(s).
anova(model_PGLS)

## __________________________________________________________________
## Interpreting the output:
##    - is there a significant effect of log(mass) on log(eyesize)?
## __________________________________________________________________


## Now look at the model coefficients (i.e. the intercept and slope)
summary(model_PGLS)

## ___________________________________________________________________________________
## Reporting the results
##    "There was a significant *negative/positive* relationship between X and Y 
##      (PGLS: slope ± SE = ?? ± ??, t = ?? df = ??, p =/</> ???, λ = ???).
## ___________________________________________________________________________________

## Plot the results:
PGLS_plot <- ggplot(mydata, aes(x = log(mass), 
                                y = log(eyesize))) +
  geom_point() +
  geom_abline(slope = coefficients(model_PGLS)[2], 
              intercept = coefficients(model_PGLS)[1]) +
  theme_minimal()
PGLS_plot