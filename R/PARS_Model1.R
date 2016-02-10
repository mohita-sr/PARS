# load all the required packages
require(lubridate)
require(coefplot)
require(boot)
require(parallel)
require(useful)
require(glmnet)
require(doParallel)
require(reshape2)
require(stringr)


# get data from tables / files
# drugDF <- read.table(file="DataSource/Drug.csv", header = T, stringsAsFactors = F, sep = ",")
# memberDF <- read.table(file="DataSource/MEMBER.csv", header = T, stringsAsFactors = F, sep =",")
# clientDF <- read.table(file="DataSource/CLIENT.csv", header = T, stringsAsFactors = F, sep =",")
# pharmacyDF <- read.table(file="DataSource/PHARMACY.csv", header = T, stringsAsFactors = F, sep =",")
# claimDF1 <- read.table(file="DataSource/Claim-1.csv", header = T, stringsAsFactors = F, sep =",")
# claimDF2 <- read.table(file="DataSource/Claim-2.csv", header = T, stringsAsFactors = F, sep =",")
# claimDF <- rbind(claimDF1, claimDF2)
# memberDetDF <- read.table(file="DataSource/MemberDet.csv", header = T, stringsAsFactors = F, sep =",")
# conditionsDF <- read.table(file="DataSource/Conditions.csv", header = T, stringsAsFactors = F, sep =",")

# head(memberDF)
# head(drugDF)
# head(clientDF)
# head(pharmacyDF)
# head(claimDF)
# head(memberDetDF)
# head(conditionsDF)
# head(parsDF)

######### Create PARS data frame ################

## step1 : Merge member data into a single set and remove duplicate columns

parsDF <- merge(memberDF, memberDetDF)
parsDF <- parsDF[,!duplicated(colnames(parsDF))]

# tempDF$WEIGHT_KG <- (parsDF[parsDF$PGM_TYP_ABBR_NM=="DB", 34]*2)
# head(tempDF[tempDF$PGM_TYP_ABBR_NM=="DB", 34])
# head(parsDF$WEIGHT_KG,20)


## step 2 : Merge the condition details on the common column COND_CD and remove duplicate column names
parsDF <- merge(parsDF, conditionsDF, by= "COND_CD")
parsDF <- parsDF[,!duplicated(colnames(parsDF))]

## step 3 : Merge the claims, drugs, pharmacy details
claimFullDF <- merge(claimDF, drugDF, by="DRUG_PROD_GID")
claimFullDF <- merge(claimFullDF,pharmacyDF, by = "PHMCY_PTY_GID")


## step 4: Full data set for PARS
parsDF <- merge(parsDF, claimFullDF,
                by=c("MBR_ACCT_GID", "LVL1_ACCT_GID","LVL3_ACCT_GID",
                     "QL_BNFCY_ID","EPH_LINK_ID"))


## STEP 5 : Calculate BMI, Frailty, Age
parsDF$BMI <- parsDF$WEIGHT_KG/parsDF$HEIGHT_M
parsDF$FRAILTY <- with(parsDF, BMI>30|BMI<.18)
parsDF$AGE <- (Sys.Date() - as.Date(parsDF$PTNT_BRTH_DT, format('%m/%d/%Y')))/365.25
parsDF$AGE <- round(as.numeric(substr(parsDF$AGE,1,length(tempDF$AGE)-5)),digits=2)


####### predict the conditions  ############
# Predictor variables :
#     Demographic :Age, Gender, Ethnicity
#     Medical History : Prior hospital admission , drugs, frailty(BMI<.18 or > 30),  co-morbidity, mental illness, specific medical diagnosis
#     Other factors : social support (lives alone), recent stressful life event     1/23/1976

tempDF <- parsDF
tempDF$COND_DIAB <- with(tempDF, PGM_TYP_ABBR_NM=='DB')


# Use Elastic net to identify the strongest variables

# Build Predictor Matrix
diabX <- build.x(COND_DIAB ~ AGE + GENDER + ENTHNICITY + BMI + DRUG_CLS_SHRT_NM +
                   PRIOR_HSP_ADMN + MBR_IN_FMLY + FRAILTY,
                 data = tempDF, contrasts=F)

# Create a Response Matrix
diabY <- build.y(COND_DIAB ~ AGE + GENDER + ENTHNICITY + BMI + DRUG_CLS_SHRT_NM +
                   PRIOR_HSP_ADMN + MBR_IN_FMLY  + FRAILTY,
                 data = tempDF)

set.seed(186356140)
diabCV1 <- cv.glmnet(x=diabX, y=diabY, family="binomial", nfold=5)

# Fit the ridge model
set.seed(7162350)
diabCV2 <- cv.glmnet(x=diabX, y=diabY, family="binomial", nfold=5, alpha=0)


### Finding the optimal value of alpha reqquires an additional layer of cross validation

# For a 2 layered cross validation, an observation should fall in the same fold each time
# Also, it is considered better to lean towards lasso than ridge, thus consider alpha > 0.5

# set the seed for repeatability
set.seed(2834673)

# create folds
theFolds <- sample(rep(x=1:5, length.out=nrow(diabX)))

# make sequence of œ values
alphas <- seq(from=0.5, to=1, by=.05)



# set seed for repeatability of random results
set.seed(5127151)

# start a cluster with 2 workers
c1 <- makeCluster(2)

# register the workers
registerDoParallel(c1)

# keep track of timing
before <- Sys.time()

# build foreachloop to run in parallel several arguments
diabDouble <- foreach(i=1:length(alphas), .errorhandling = "pass",
                     .inorder = F, .multicombine = T,
                     .export = c("diabX", "diabY", "alphas", "theFolds"),
                     .packages = "glmnet")   %dopar%
                     {
                       print(alphas[i])
                       cv.glmnet(x=diabX, y=diabY, family="binomial", nfolds=5,
                                 foldid = theFolds, alpha=alphas[i])
                     }

after <- Sys.time()
stopCluster(c1)

after - before


# to find the optimal value of lambda and alpha
# function for extracting info from cv.glmnet
extractGlmnetInfo <- function(object)
{
  #find lambdas
  lambdaMin <- object$lambda.min
  lambda1se <- object$lambda.1se

  #figure out where those lambdas fall in the path
  whichMin <- which(object$lambda == lambdaMin)
  which1se <- which(object$lambda == lambda1se)

  data.frame(lambda.min=lambdaMin, error.min=object$cvm[whichMin],
             lambda.1se=lambda1se, error.1se=object$cvm[which1se])
}


#build a one line data frame with each lambda & its corr error figure
alphaInfo <- plyr::ldply(diabDouble,extractGlmnetInfo)

alphaInfo$Alpha <- alphas


alphaMelt <- melt(alphaInfo, id.vars="Alpha", value.name = "Value", variable.name="Measure")
alphaMelt$Type <- str_extract(string=alphaMelt$Measure, pattern = "(min)|(1se)")

# housekeeping :)
alphaMelt$Measure <- str_replace(string=alphaMelt$Measure, pattern = "\\.(min|1se)",
                                 replacement = "")
alphaCast <- dcast(alphaMelt, Alpha + Type ~ Measure, value.var = "Value")

ggplot(alphaCast, aes(x=Alpha, y=error))+
  geom_line(aes(group=Type)) +
  facet_wrap(~Type, scales="free_y", ncol=1) +
  geom_point(aes(size=lambda))


# we have found the optimal value of œ = 0.85  alphaInfo$Alpha[which.min(alphaInfo$error.1se)]
# we fit the model accordingly

set.seed(5127151)
diabCV3 <- cv.glmnet(x=diabX, y=diabY, family="binomial", nfold=5,
                    alpha = alphaInfo$Alpha[which.min(alphaInfo$error.1se)])
plot(diabCV3)
plot(diabCV3$glmnet.fit, xvar="lambda")
abline(v=log(c(diabCV3$lambda.min, diabCV3$lambda.1se)), lty=2)


#viewing the Coeff plot for glmnet

theCoef <- as.matrix(coef(diabCV3, s="lambda.1se"))
coefDF <- data.frame(Value=theCoef, coefficient=rownames(theCoef))
coefDF <- coefDF[nonzeroCoef(coef(diabCV3, s="lambda.1se")),]

ggplot(coefDF, aes(x=X1, y=reorder(coefficient, X1))) +
  geom_vline(xintercept = 0, color="grey", linetype=2 ) +
  geom_point(color="blue") + labs(x="Value", y="Coefficient", title="Coefficient Plot")


