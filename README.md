
# gesttools Version 1.3.0

<!-- badges: start -->
<!-- badges: end -->

General Purpose g-Estimation for End ofStudy or Time-Varying Outcomes.  
Provides a series of general purpose tools to perform g-estimation using the methods described in [Sjolander and Vansteelandt (2016)](https://www.degruyter.com/document/doi/10.1515/em-2015-0005/html) and [Dukes and Vansteelandt (2018)](https://academic.oup.com/aje/article/187/5/1079/4930819). 
The package allows for g-estimation in a wide variety of circumstances, including an end of study or time-varying outcome, and an exposure that is a binary, continuous, or a categorical variable with three or more categories.
The package also supports g-estimation with time-varying causal effects and effect modification by a confounding variable.


## Installation

You can install gesttools version 1.3.0 from R using:

``` r
install.packages("gesttools")
```
For Further Information. Please refer to the User Manual Provided on GitHub or from CRAN.

## Example

We provide a basic example of gesttools with comments below.

``` r
 #Load the Package
 library(gesttools)
 #Create a simple dataset with exposure A and outcome Y. 
 #Format the data using the FormatData() function and generate history of exposure.
 datas <- dataexamples(n = 1000, seed = 123, Censoring = FALSE)
 data <- datas$datagest
 data <- FormatData(
 data = data, idvar = "id", timevar = "time", An = "A",
 varying = c("Y", "A", "L"), GenerateHistory = TRUE, GenerateHistoryMax = 1)
 #Define Inputs
 idvar <- "id"
 timevar <- "time"
 Yn <- "Y"
 An <- "A"
 Cn <- NA
 #Define the outcome models for each of the three exposure times
 outcomemodels <- list("Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A", "Y~A+L+U+Lag1A")
 #Define the propensity score model of being exposed. Note that time is included 
 propensitymodel <- c("A~L+U+as.factor(time)+Lag1A")
 censoringmodel <- NULL
 EfmVar <- NA
 #Perform g-estimation 
 gestSingle(data, idvar, timevar, Yn, An, Cn, outcomemodels, propensitymodel,
 censoringmodel = NULL, type = 1, EfmVar)

```
