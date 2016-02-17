diabFormula <- COND_DIAB ~ AGE + GENDER + ENTHNICITY + BMI + DRUG_CLS_SHRT_NM + PRIOR_HSP_ADMN + MBR_IN_FMLY  + FRAILTY

diabFormX <- build.x(diabFormula, data = tempDF)
diabFormY <- build.y(diabFormula, data = tempDF)

diabForest <- randomForest(x=diabFormX, y=diabFormY)

diabForest
