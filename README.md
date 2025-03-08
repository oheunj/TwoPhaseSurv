# Leveraging Two-Phase Data for Improved Prediction of Survival Outcomes with Application to Nasopharyngeal Cancer
The manuscript is currently under review:
* __Oh, E. J.,__ Ahn, S., Tham, T., and Qian, M. (under review). Leveraging two-phase data for improved prediction of survival outcomes with application to nasopharyngeal cancer.

# Files in this repository
The source code is currently provided in 'main_sourecode.R'

# Installation
R is a statistical software program, and RStudio is a user interface for R. We recommend that users install both R and R Studio. Both R and RStudio are free and open source.

## Requirements
You may need to install the following dependencies first:
```{r}
library(MASS)
library(survival)
library(survminer)
library(knitr)
library(partykit)
library(MatchIt)
library(expss)
library(glmnet)
library(fastDummies)
library(randomForestSRC)
library(ipred)
library(matrixStats)
library(survAUC)
library(SurvMetrics)
library(pec)
library(RANN)
library(Hmisc)
library(caret)
library(DescTools)
library(mice)
library(smcfcs)
library(survcomp) # this package is now available from Bioconductor only
```

# License
```{r}
Licensed under the GNU General Public License v3.0 (GPL-3.0);
you may not use this file except in compliance with the License.
``
