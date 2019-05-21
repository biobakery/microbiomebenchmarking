#########
# ANCOM #
#########

#####################################
# Install or Load Required Packages #
#####################################

if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('exactRankTests', 'openxlsx', 'DT', 'dplyr', 'coin')
if(! require("ancom.R")) {
  download.file("https://www.niehs.nih.gov/research/resources/software/biostatistics/ancom/ancom_software.zip",destfile = "ANCOM.zip")
  unzip("ANCOM.zip",exdir=getwd())
  install.packages("ancom.R_1.1-3.tar.gz", repos = NULL)
  suppressPackageStartupMessages(library("ancom.R"))
  file.remove('ancom.R_1.1-3.tar.gz')
  file.remove('ancom.R.zip')
  file.remove('ANCOM.zip')
  file.remove('README.First.ANCOM.1.1-3.pdf')
}
library(ancom.R)


##########################
# Fit ANCOM To A Dataset #
##########################

fit.ANCOM = function(features, metadata, libSize, ID, transformation, MultTestCorrection) {
  
  if (transformation!='NONE') stop ('Transformation currently not supported for a default ANCOM model. Use NONE.')
  
  ## ANCOM standard pipeline for DA ##
  # Random Effect Adjustment
  input <- as.data.frame(features)
  input$Group <- metadata
  input$Group <- ifelse(input$Group=='1', 'Group1', 'Group0')
  if(!length(ID)==length(unique(ID))){
    input$ID<-ID
    res <- ANCOM(input, sig=0.05, multcorr = 2, repeated=TRUE) # Same as Weiss et al., 2017
  } else{
    res <- ANCOM(input, sig=0.05, multcorr = 2) # Same as Weiss et al., 2017
    }
  
  # Extract results and enforce meaningful format
  df <- data.frame(coef = res[[1]]) 
  df$pval<-1 # Fake p-values
  df$feature = colnames(features)
  df$metadata = names(metadata)
  df[df$feature %in% res[[2]], 'pval'] <- 0 # Fake p-values
  df$qval<-df$pval # Fake q-values
  df<-df[order(df$qval, decreasing=FALSE),]
  df<-dplyr::select(df, c('feature', 'metadata'), everything())
  rownames(df)<-NULL;
  return(df)
}

###################################
# Fit ANCOM To A List of Datasets #
###################################

list.ANCOM<-function(physeq, transformation='NONE', MultTestCorrection = 'BH'){
  foreach(physeq=physeq, 
          .export="fit.ANCOM",
          .packages=c("dplyr", "exactRankTests", 'openxlsx', "DT", "ancom.R"), 
          .errorhandling = 'remove') %dopar% 
          {
            start.time <- Sys.time()
            features<-physeq$features; 
            metadata<-physeq$metadata;
            libSize <- physeq$libSize;
            ID<-physeq$ID;
            DD<-fit.ANCOM(features, metadata, libSize, ID, transformation, MultTestCorrection)
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



