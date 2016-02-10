# load all the required packages
install.packages("RMySQL")
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


myDB = dbConnect(MySQL(), user='***', password='****!', dbname=pars)


