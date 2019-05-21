##########
# negbin #
##########

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('dplyr', 'pbapply', 'MASS', 'lme4')

######################################
# Fit Negative Binomial To A Dataset #
######################################

fit.negbin <- function(features, metadata, libSize, ID, transformation, MultTestCorrection){
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a Negative Binomial model. Use NONE.')
  
  paras <- pbapply::pbsapply(1:ncol(features), simplify=FALSE, function(x){
    
    featuresVector <- round(features[, x])
    
    # Scrap All-Zero Features
    if(sum(featuresVector!=0)<1){
      print(paste("Cannot fit model to all zeroes for feature", x, "returning NA"))
      para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
      colnames(para)<-c('coef', 'pval')
      para$metadata<-colnames(metadata)
      para$feature<-colnames(features)[x]
      rownames(para)<-NULL
    }   
    else{
      
      # Fit Model
      dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata, libSize)
      formula<-as.formula(paste("expr ~ ", paste(colnames(metadata), collapse= "+")))
      
      # Library size adjustment
      formula<-update(formula, . ~ . - offset(log(libSize)))
      
      # Random effect adjustment
      if(!length(ID)==length(unique(ID))){
        formula<-update(formula, . ~ . +(1|ID))
          fit <- tryCatch({
            fit1 <- lme4::glmer.nb(formula, data = dat_sub)
          }, error=function(err){
            fit1 <- try({lme4::glmer.nb(formula, data = dat_sub)}) 
            return(fit1)
          })
          
          if (class(fit) != "try-error"){
            para<-as.data.frame(summary(fit)$coefficients)[-1,-c(2:3)]
            colnames(para)<-c('coef', 'pval')
            para$metadata<-colnames(metadata)
            para$feature<-colnames(features)[x]
            rownames(para)<-NULL
          }
          else{
            print(paste("Fitting problem for feature", x, "returning NA"))
            para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
            colnames(para)<-c('coef', 'pval')
            para$metadata<-colnames(metadata)
            para$feature<-colnames(features)[x]
            rownames(para)<-NULL
          }
          } else{
          fit <- tryCatch({
            fit1 <- MASS::glm.nb(formula, data = dat_sub)
          }, error=function(err){
            fit1 <- try({MASS::glm.nb(formula, data = dat_sub)}) 
            return(fit1)
          })
          
          if (class(fit) != "try-error"){
            para<-as.data.frame(summary(fit)$coefficients)[-1,-c(2:3)]
            colnames(para)<-c('coef', 'pval')
            para$metadata<-colnames(metadata)
            para$feature<-colnames(features)[x]
            rownames(para)<-NULL
          }
          else{
            print(paste("Fitting problem for feature", x, "returning NA"))
            para<- as.data.frame(matrix(NA, nrow=ncol(metadata), ncol=2))
            colnames(para)<-c('coef', 'pval')
            para$metadata<-colnames(metadata)
            para$feature<-colnames(features)[x]
            rownames(para)<-NULL
          }
        }
        } 
    return(para)
  })
  
  paras<-do.call(rbind, paras)
  paras$qval<-as.numeric(p.adjust(paras$pval, method = MultTestCorrection))
  paras<-paras[order(paras$qval, decreasing=FALSE),]
  paras<-dplyr::select(paras, c('feature', 'metadata'), everything())
  rownames(paras)<-NULL;
  return(paras) 
}

###############################################
# Fit Negative Binomial To A List of Datasets #
###############################################

list.negbin <-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .packages = c("MASS", "pbapply", "dplyr", "lme4"),
          .export="fit.negbin", 
          .errorhandling = 'remove') %dopar% 
  {
    start.time <- Sys.time()
    features<-physeq$features; 
    metadata<-physeq$metadata;
    libSize <- physeq$libSize;
    ID<-physeq$ID;
    DD<-fit.negbin(features, metadata, libSize, ID, transformation, MultTestCorrection)
    DD$pairwiseAssociation<-paste('pairwiseAssociation', 1:nrow(DD), sep='')
    wh.TP = intersect(grep("[[:print:]]+\\_TP$", DD$metadata), grep("[[:print:]]+\\_TP$", DD$feature))
    newname = paste0(DD$pairwiseAssociation[wh.TP], "_TP")
    DD$pairwiseAssociation[wh.TP] <- newname;
    DD<-dplyr::select(DD, c('pairwiseAssociation', 'feature', 'metadata'), everything())
    stop.time <- Sys.time()
    time<-as.numeric(round(difftime(stop.time, start.time, units="min"),3), units = "mins")
    DD$time<-time
    return(DD)
  }
}

