#!bin/bash/Rscript
#### Custom scripts for tick pathogen 


stepGAM <- function(dat, predictors, response, constant_pred=NULL, family, ignore.combos) {
  ## For troubleshooting
  # dat <- tck_borrelia_adj
  # predictors = allPred
  # response = "borrPresent"
  # constant_pred <- "offset(log(numberTested))"
  # family=binomial
  # ignore.combos=list(c("domainID","s(domainID, bs='re')"), c("s(logNLtckDensity)","s(logadultDensity)","s(lognymphDensity)")
  #                    , c("s(logNLtckDensity)","s(logadultDensity)"),c("s(logNLtckDensity)","s(lognymphDensity)"))
  n_rows <- sum(sapply(1:length(predictors),FUN=function(x){ncol(combn(predictors,m=x))} ))
  allAIC <- setNames(data.frame(matrix(nrow=n_rows, ncol=(4+length(predictors)))),nm = c("AIC","REML","Dev.expl","formula",predictors))
  i <- 1
  pb <- txtProgressBar(title = "progress bar", min = 0,
                       max = n_rows, style=3)
  for ( p in 1:length(predictors) ) {
    combos <- combn(predictors, m=p)
    for ( c in 1:ncol(combos) ) {
      # skip any combinations that are meant to be ignored
      if (any(sapply(ignore.combos, FUN=function(ic){sum(ic %in% combos[,c]) == length(ic)}))) {
        i = i + 1
        next
      }
      if ( length(constant_pred)>0 ) {
        frml.temp <- paste(c(paste0(response," ~ ",constant_pred),combos[,c]), collapse="+")
      } else {
        frml.temp <- paste0(paste0(response," ~ "),paste(combos[,c], collapse=" + "), sep="")
      }
      # If more coefficients than data
      res <- try(gam(as.formula(frml.temp)
                     , data=dat
                     , method="REML"
                     , family=family))
      if ( inherits(res, "try-error") ) {
        next
      } else {
        gam.temp <- gam(as.formula(frml.temp)
                        , data=dat
                        , method="REML"
                        , family=family)
        
        pred.temp <- (colnames(allAIC) %in% combos[,c])[-c(1,2,3,4)]
        allAIC[i,c("AIC","REML","Dev.expl","formula")]  <- as.vector(c(gam.temp$aic,gam.temp$gcv.ubre[1],summary(gam.temp)$dev.expl,frml.temp))
        allAIC[i, predictors] <- pred.temp
        
      }
      
      i <- i+1
      Sys.sleep(0.1)
      setTxtProgressBar(pb, i, label=paste( round(i/total*100, 0),
                                            "% done"))
    }
  }
  close(pb)
  return(allAIC)
}