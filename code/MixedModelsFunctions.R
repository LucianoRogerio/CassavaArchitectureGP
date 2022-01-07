
if (!any(c("sommer", "lme4")%in% installed.packages())) {
  if(!"sommer" %in% installed.packages()) {install.packages("sommer")
    library(sommer)}
  if(!"lme4" %in% installed.packages()) {install.packages("lme4")
    library(lme4)}
} else {library(sommer); library(lme4)}

#### Mixed models analysis for trials with augmented design
#### e completos - lme4 package
analyzeTrial.lme4 <- function(x){
    modfit <- lmer(Value ~ control + block_number + (1 | accession_name:new), data=x, REML = T)
  return(modfit)
}

analyzeTrialrdMod.lme4 <- function(x, trait){
  modfit <- lm(Value ~ control + block_number, data=x)
  return(modfit)
}

Deviance.MM <- function(x, y) {
  anova(x, y)
}

getVarComp.lme4 <- function(model){
  as.data.frame(VarCorr(model))[,c("grp","vcov")] %>%
    rename(c("grp" = "effect", "vcov" = "VarEstimate")) %>%
    spread(key = effect, value = VarEstimate, fill=NA)
}


getVarOutput.lme4 <- function(modfit){
  var_label <- row.names(modfit)
  if (any(is.na(modfit))){
    varcomp <- data.frame(NA, NA)
  } else{
    varcomp <- data.frame(matrix(unlist(modfit), nrow=nrow(modfit), byrow=T))
  }
  colnames(varcomp) <- c("gen_var", "res_var")
  varcomp <- round(varcomp[,c("gen_var","res_var")], 3)
  var.out <- cbind(var_label, varcomp)
  return(var.out)
}



getBLUPs.lme4 <- function(model){
  as.data.frame(ranef(model)$"accession_name:new" +
                  fixef(model)["(Intercept)"]) %>%
    rename(c("(Intercept)" = names(model)))
     # spread(key = Clone, value = BLUP, fill=NA)
}


##### #### Analise de modelos mistos para ensaio com delineamento de blocos aumentados
#### e completos - sommer package
analyzeTrial.sommer <- function(x){
  if(any(x$studyDesign == "Augmented")){
    #### Augmented block design
    modfit <- mmer(Value ~ check + blockNumber,
                   random = ~ germplasmName:new,
                   data = x,
                   tolparinv = 1e-6,
                   verbose = F,
                   method = "EM",
                   getPEV = T)
  } else if(any(x$studyDesign %in% c("RCBD", "Alpha"))) {
    #### Randomized Complete block design
    modfit <- mmer(Value ~ blockNumber,
                   random = ~ germplasmName,
                   data = x,
                   tolparinv = 1e-6,
                   verbose = F,
                   method = "EM",
                   getPEV = T)
  } else {
    #### Randomized design
    modfit <- mmer(Value ~ 1,
                   random = ~ germplasmName,
                   data = x,
                   tolparinv = 1e-6,
                   verbose = F,
                   method = "EM",
                   getPEV = T)
  }
  return(modfit)
}


getVarComp.sommer <- function(model){
  unlist(model$sigma) %>%
    data.frame(effect = c("Clone", "Residual"), VarEstimate = .) %>%
    spread(key = effect, value = VarEstimate, fill=NA)
}

