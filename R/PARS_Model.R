# load all the required packages
require(lubridate)
require(coefplot)
require(boot)

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

#Model 1
beforeTime <- Sys.time()
diab1 <- glm(COND_DIAB ~ AGE + GENDER + ENTHNICITY + BMI + DRUG_CLS_SHRT_NM +
               PRIOR_HSP_ADMN + MBR_IN_FMLY + COND_CD + FRAILTY
             ,data=tempDF, family=binomial(link="logit"))
afterTime <- Sys.time()
summary(diab1)


# model2
beforeTime <- Sys.time()
diab2 <- glm(COND_DIAB ~ AGE + GENDER + ENTHNICITY + BMI + DRUG_CLS_SHRT_NM +
               PRIOR_HSP_ADMN + MBR_IN_FMLY  + FRAILTY
             ,data=tempDF, family=binomial(link="logit"))
afterTime <- Sys.time()
summary(diab2)
sprintf("Time elapsed %s", afterTime-beforeTime)

# model3
beforeTime <- Sys.time()
diab3 <- glm(COND_DIAB ~ AGE + GENDER + ENTHNICITY + BMI + DRUG_CLS_SHRT_NM +
               PRIOR_HSP_ADMN  + FRAILTY
             ,data=tempDF, family=binomial(link="logit"))
afterTime <- Sys.time()
summary(diab3)
sprintf("Time elapsed %s", afterTime-beforeTime)

# model4
beforeTime <- Sys.time()
diab4 <- glm(COND_DIAB ~ AGE + GENDER + ENTHNICITY + BMI + DRUG_CLS_SHRT_NM +
               PRIOR_HSP_ADMN + COND_CD + FRAILTY
             ,data=tempDF, family=binomial(link="logit"))
afterTime <- Sys.time()
summary(diab4)
sprintf("Time elapsed %s", afterTime-beforeTime)


invlogit <- function(x)
{
  1/(1+ exp(-x))
}

invlogit(diab1$coefficients)

anova(diab1, diab2, diab3, diab4)
AIC(diab1, diab2, diab3, diab4)
BIC(diab1, diab2, diab3, diab4)

# Cross Validation
diabCV1 <- cv.glm(tempDF, diab1, K=5)
diabCV2 <- cv.glm(tempDF, diab2, K=5)
diabCV3 <- cv.glm(tempDF, diab3, K=5)
diabCV4 <- cv.glm(tempDF, diab4, K=5)

cvResults <- as.data.frame(rbind(diabCV1$delta, diabCV2$delta, diabCV3$delta, diabCV4$delta ))
names(cvResults) <- c("Error", "Adjusted Error")
cvResults$Model <- sprintf("diab%s",1:4)
cvResults

coefplot(diab2)
