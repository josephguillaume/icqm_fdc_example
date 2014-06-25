## nsga2
if(!require("mco")) install.packages("mco")
## PRIM
if(!require("sdtoolkit")) install.packages("sdtoolkit")
## data
for(x in c("zoo", "latticeExtra", "polynom", "car", "Hmisc","reshape"))
  if(!require(package=x,character.only=TRUE)) install.packages(x)
if(!require("hydromad")) install.packages("hydromad", repos="http://hydromad.catchment.org")
## DREAM
if(!require("dream")) install.packages("dream", repos="http://hydromad.catchment.org")

################################################################################
##_* Load required packages
library(hydromad) #data
library(dream)  #MCMC method (M7)
library(mco) #NSGA2 multi-objective optimisation (M3,M8)
library(sdtoolkit) #Scenario discovery (M2)