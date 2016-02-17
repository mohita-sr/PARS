# load all the required packages
install.packages("RMySQL")
install.packages("randomForest")
require(RMySQL)
require(lubridate)
require(coefplot)
require(boot)
require(parallel)
require(useful)
require(glmnet)
require(doParallel)
require(reshape2)
require(stringr)
require(useful)
require(randomForest)

myDB = dbConnect(MySQL(), user='***', password='****!', dbname=pars)


