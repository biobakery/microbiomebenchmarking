##############################
## Synthetic Data Generation #
##############################

# Generate Replicated Simulated Datasets For A Combination of Parameters
trigger_sparseDOSSA_Simulator<-function(noZeroInflate=FALSE,
                                        RandomEffect=FALSE,
                                        metadataType,
                                        nSubjects,
                                        nPerSubject,
                                        nMicrobes,
                                        spikeMicrobes,
                                        nMetadata,
                                        spikeMetadata,
                                        effectSize,
                                        readDepth = 50000,
                                        nIterations = 100,
                                        noParallel = FALSE,
                                        rSeed = 1234,
                                        nCores = 4){
  
  # Create Replicates 
  reps = 1:nIterations
  
  ########################
  # Catch Obvious Errors #
  ########################
  
  # Check Character Values
  if (!metadataType %in% c('UVA', 'UVB', 'MVA', 'MVB'))
    stop('Must be one of the following: UVA, UVB, MVA, or MVB.')
  
  # Check Positive Integer Values
  if (round(nSubjects) != nSubjects || 
      nSubjects<0 ||
      round(nPerSubject) != nPerSubject || 
      nPerSubject<0 ||
      round(nMicrobes) != nMicrobes || 
      nMicrobes<0 ||
      round(nMetadata) != nMetadata ||
      nMetadata<0 ||
      round(readDepth) != readDepth || 
      readDepth<0)
    stop('nSubjects/nPerSubject/nMicrobes/nMetadata/readDepth must be positive integers.')
  
  # Check Proportion Values
  if (spikeMicrobes>1 || spikeMicrobes<=0 || spikeMetadata<=0 || spikeMetadata>1)
    stop('spikeMicrobes/spikeMetadata must be in (0, 1].')
  
  # Check Illegal Combinations 
  if(RandomEffect==TRUE && nPerSubject==1)
    stop('nPerSubject must be greater 1 when RandomEffect is TRUE.')
  
  if(RandomEffect==FALSE && nPerSubject>1)
    stop('nPerSubject must be equal to  1 when RandomEffect is FALSE.')
  
  if(metadataType %in% c('UVA', 'UVB') && (nMetadata!=1 || spikeMetadata!=1))
    stop('Both nMetadata and spikeMetadata must be equal to 1 when metadataType is UVA or UVB.')
  
  if(!metadataType %in% c('UVA', 'UVB') && nMetadata==1)
    stop('nMetadata must be greater than 1 when metadataType is MVA or MVB')
  
  # Define the Simulation Parameters Combinations
  simparams = apply(expand.grid(metadataType,
                                nSubjects,
                                nPerSubject,
                                nMicrobes, 
                                spikeMicrobes, 
                                nMetadata,
                                spikeMetadata,
                                effectSize,
                                readDepth,
                                reps), 1, paste, collapse = '_')
  
  # Define the Labels to Go with Each Element of the Simulation Parameter
  simparamslabels = c("metadataType","nSubjects", "nPerSubject", "nMicrobes", "spikeMicrobes", "nMetadata", "spikeMetadata", "effectSize", "readDepth", "rep")
  
  # Track Start Time
  cat(c("Job started at:",date()), "\n")
  start.time <- Sys.time()
  
  # Set Reproducibility Seed
  set.seed(rSeed) 
  
  # Call Grid Computing Only When Specified
  
  if (noParallel){
    
    # Call SparseDOSSA Wrapper (noParallel)
    simlist <- sparseDOSSA_Wrapper_noParallel(simparams, simparamslabels, noZeroInflate=noZeroInflate)
  }
  else{
    
    # Set Up Clustering Environment
    no_cores <- nCores 
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    
    ####################
    # Data Generation #
    ###################
    
    # Call SparseDOSSA Wrapper 
    simlist <- sparseDOSSA_Wrapper(simparams, simparamslabels, noZeroInflate=noZeroInflate)
    
    # Stop the Cluster 
    stopCluster(cl)
  }
  
  # Set Names
  if (noZeroInflate==TRUE && RandomEffect==TRUE) {
    simnames<- paste('noZeroInflate_RandomEffect', simparams, sep='_')} 
  if (noZeroInflate==TRUE && RandomEffect==FALSE) {
    simnames<- paste('noZeroInflate_noRandomEffect', simparams, sep='_')}  
  if (noZeroInflate==FALSE && RandomEffect==TRUE) {
    simnames<- paste('ZeroInflate_RandomEffect', simparams, sep='_')} 
  if (noZeroInflate==FALSE && RandomEffect==FALSE) {
    simnames<- paste('ZeroInflate_noRandomEffect', simparams, sep='_')} 
  names(simlist) <- simnames
  
  # Track End Time
  cat(c("Job ended at:",date()), "\n")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"),3)
  cat("Computational time:", minutes, "minutes \n")

  # Return
  return(simlist)
}

# Trigger sparseDOSSA 
sparseDOSSA_Wrapper<-function(simparams, simparamslabels, noZeroInflate){
  f<-foreach(i = simparams, .packages = c("sparseDOSSA", "MASS", "stringi"),
          .export = c("generateMetadata"), .errorhandling = 'remove') %dopar% {
            
            # Extract Parameter Strings
            params = strsplit(i, '_')[[1]]
            names(params) <- simparamslabels
            
            # Extract Relevant Parameters
            metadataType = as.character(params["metadataType"]) # Type of Metadata
            nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
            nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
            nSamples<-round(nSubjects*nPerSubject) # Number of Samples
            nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
            spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
            nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
            spikeMetadata<-as.numeric(params["spikeMetadata"])  # Proportion of Spiked-in Metadata
            effectSize<-as.character(params["effectSize"]) # Effect Size
            readDepth<-as.numeric(params["readDepth"]) # Library Size
            
            
            # Initialize
            DD = NULL
            
            # sparseDOSSA Error Control 
            tryAgain = TRUE
            infiniteloopcounter = 1
            while (tryAgain & infiniteloopcounter < 5) {
              
              # Generate Metadata
              FF<-generateMetadata(metadataType=metadataType, 
                                   nSubjects=nSubjects, 
                                   nPerSubject=nPerSubject, 
                                   nMetadata=nMetadata, 
                                   spikeMetadata=spikeMetadata)
              
              # Extract Relevant Information
              UserMetadata<-FF$UserMetadata; 
              Metadatafrozenidx<-FF$Metadatafrozenidx;
              significant_metadata<-FF$significant_metadata; 
              spikeCount<-FF$spikeCount
              
              # Generate sparseDOSSA Synthetic Abundances
              DD<-sparseDOSSA::sparseDOSSA(number_features = nMicrobes,
                                           number_samples = nSamples,
                                           UserMetadata = UserMetadata,
                                           Metadatafrozenidx = Metadatafrozenidx,
                                           datasetCount = 1,
                                           spikeCount = spikeCount,
                                           spikeStrength = effectSize,
                                           read_depth = readDepth,
                                           noZeroInflate = noZeroInflate,
                                           percent_spiked=spikeMicrobes,
                                           write_table = FALSE)
              if (is.null(DD) | inherits(DD, "try-error")) {
                tryAgain = TRUE
                infiniteloopcounter = infiniteloopcounter + 1
              } else {
                tryAgain = FALSE
              }
            }
            if (infiniteloopcounter >= 5) {
              stop("Consistent error found during simulation. Need to investigate cause.")
            }
            
            # Gather sparseDOSSA Outputs
            sparsedossa_results <- as.data.frame(DD$OTU_count)
            rownames(sparsedossa_results)<-sparsedossa_results$X1
            sparsedossa_results<-sparsedossa_results[-1,-1]
            colnames(sparsedossa_results)<-paste('Sample', 1:ncol(sparsedossa_results), sep='')
            data<-as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)),])
            data<-data.matrix(data)
            class(data) <- "numeric"
            truth<-c(unlist(DD$truth))
            truth<-truth[!stringi::stri_detect_fixed(truth,":")]
            truth<-truth[(5+nMetadata):length(truth)]
            truth<-as.data.frame(truth)
            significant_features<-as.vector(truth[seq(1, (as.numeric(spikeCount)+1)*(nMicrobes*spikeMicrobes), (as.numeric(spikeCount)+1)),])
            
            # Separate Metadata and Taxa 
            
            # Extract Metadata
            if (metadataType %in% c('UVA', 'UVB')){
              metadata<-as.data.frame(data[1,])
              colnames(metadata)<-rownames(data)[1]
            } else{
              metadata<-as.data.frame(t(data[(1:nMetadata),]))
            }
            
            # Mark True Positive Metadata - Same Format at Mcmurdie and Holmes (2014)
            which.TP = colnames(metadata) %in% significant_metadata
            meta_newname = paste0(colnames(metadata)[which.TP], "_TP")
            colnames(metadata)[which.TP] <- meta_newname
            
            # Extract Features 
            features<-as.data.frame(t(data[-c(1:nMetadata),]))
            
            # Mark True Positive Features - Same Format at Mcmurdie and Holmes (2014)
            wh.TP = colnames(features) %in% significant_features
            colnames(features)<-paste("Feature", 1:nMicrobes, sep = "")
            newname = paste0(colnames(features)[wh.TP], "_TP")
            colnames(features)[wh.TP] <- newname;
            
            # Add Sample ID 
            ID<-rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)
            
            # Add Library Size (or Sequencing Depth)
            libSize<-rowSums(features)
            
            # Return
            return(list(metadata=metadata, features=features, ID=ID, libSize=libSize))
          }
  return(f)
}

# Trigger sparseDOSSA (noParallel)
sparseDOSSA_Wrapper_noParallel<-function(simparams, simparamslabels, noZeroInflate){
  
  # Intitialize
  pclList<-list()
  
  # Repeated Loop 
  for(i in simparams){
    
    # Extract Parameter Strings
    params = strsplit(i, '_')[[1]]
    names(params) <- simparamslabels
    
    # Extract Relevant Parameters
    metadataType = as.character(params["metadataType"]) # Type of Metadata
    nSubjects <- as.numeric(params["nSubjects"])  # Number of Subjects
    nPerSubject <- as.numeric(params["nPerSubject"])  # Number of Samples Per Subject
    nSamples<-round(nSubjects*nPerSubject) # Number of Samples
    nMicrobes <- as.numeric(params["nMicrobes"])  # Number of Microbes
    spikeMicrobes <- as.numeric(params["spikeMicrobes"]) # Proportion of Spiked-in Microbes
    nMetadata<-as.numeric(params["nMetadata"])  # Number of Metadata
    spikeMetadata<-as.numeric(params["spikeMetadata"])  # Proportion of Spiked-in Metadata
    effectSize<-as.character(params["effectSize"]) # Effect Size
    readDepth<-as.numeric(params["readDepth"]) # Library Size
    
    
    # Initialize
    DD = NULL
    
    # sparseDOSSA Error Control 
    tryAgain = TRUE
    infiniteloopcounter = 1
    while (tryAgain & infiniteloopcounter < 5) {
      
      # Generate Metadata
      FF<-generateMetadata(metadataType=metadataType, 
                           nSubjects=nSubjects, 
                           nPerSubject=nPerSubject, 
                           nMetadata=nMetadata, 
                           spikeMetadata=spikeMetadata)
      
      # Extract Relevant Information
      UserMetadata<-FF$UserMetadata; 
      Metadatafrozenidx<-FF$Metadatafrozenidx;
      significant_metadata<-FF$significant_metadata; 
      spikeCount<-FF$spikeCount
      
      # Generate sparseDOSSA Synthetic Abundances
      DD<-sparseDOSSA::sparseDOSSA(number_features = nMicrobes,
                                   number_samples = nSamples,
                                   UserMetadata = UserMetadata,
                                   Metadatafrozenidx = Metadatafrozenidx,
                                   datasetCount = 1,
                                   spikeCount = spikeCount,
                                   spikeStrength = effectSize,
                                   read_depth = readDepth,
                                   noZeroInflate = noZeroInflate,
                                   percent_spiked = spikeMicrobes,
                                   write_table = FALSE)
      
      if (is.null(DD) | inherits(DD, "try-error")) {
        tryAgain = TRUE
        infiniteloopcounter = infiniteloopcounter + 1
      } else {
        tryAgain = FALSE
      }
    }
    if (infiniteloopcounter >= 5) {
      stop("Consistent error found during simulation. Need to investigate cause.")
    }
    
    # Gather sparseDOSSA Outputs
    sparsedossa_results <- as.data.frame(DD$OTU_count)
    rownames(sparsedossa_results)<-sparsedossa_results$X1
    sparsedossa_results<-sparsedossa_results[-1,-1]
    colnames(sparsedossa_results)<-paste('Sample', 1:ncol(sparsedossa_results), sep='')
    data<-as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)),])
    data<-data.matrix(data)
    class(data) <- "numeric"
    truth<-c(unlist(DD$truth))
    truth<-truth[!stringi::stri_detect_fixed(truth,":")]
    truth<-truth[(5+nMetadata):length(truth)]
    truth<-as.data.frame(truth)
    significant_features<-as.vector(truth[seq(1, (as.numeric(spikeCount)+1)*(nMicrobes*spikeMicrobes), (as.numeric(spikeCount)+1)),])
    
    # Separate Metadata and Taxa 
    
    # Extract Metadata
    if (metadataType %in% c('UVA', 'UVB')){
      metadata<-as.data.frame(data[1,])
      colnames(metadata)<-rownames(data)[1]
    } else{
      metadata<-as.data.frame(t(data[(1:nMetadata),]))
    }
    
    # Mark True Positive Metadata - Same Format at Mcmurdie and Holmes (2014)
    which.TP = colnames(metadata) %in% significant_metadata
    meta_newname = paste0(colnames(metadata)[which.TP], "_TP")
    colnames(metadata)[which.TP] <- meta_newname
    
    # Extract Features 
    features<-as.data.frame(t(data[-c(1:nMetadata),]))
    
    # Mark True Positive Features - Same Format at Mcmurdie and Holmes (2014)
    wh.TP = colnames(features) %in% significant_features
    colnames(features)<-paste("Feature", 1:nMicrobes, sep = "")
    newname = paste0(colnames(features)[wh.TP], "_TP")
    colnames(features)[wh.TP] <- newname;
    
    # Add Sample ID 
    ID<-rep(paste('Subject', 1:nSubjects, sep=''), each = nPerSubject)
    
    # Add Library Size (or Sequencing Depth)
    libSize<-rowSums(features)
    
    # Save
    pclList[[i]]<-list(metadata=metadata, features=features, ID=ID, libSize=libSize)
  }
  
  # Return
  return(pclList)
}



#####################
# Generate Metadata #
#####################

generateMetadata<-function(metadataType, 
                           nSubjects, 
                           nPerSubject, 
                           nMetadata, 
                           spikeMetadata){
  
  # Calculate Number of Samples
  nSamples = round(nSubjects*nPerSubject)
  
  # Create Blocking Variable
  if (nPerSubject==1){  # NO RANDOM EFFECTS 
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=0))
  }
  if (nPerSubject>1){  # SUBJECT-SPECIFIC RANDOM EFFECTS
    subjectRandomEffects <- as.matrix(rnorm(nSubjects,mean=0,sd=1))
  }
  BLOCK <- as.vector(matrix(subjectRandomEffects,nrow=nPerSubject,ncol=length(subjectRandomEffects),byrow=TRUE))
  
  # Specify Mean and Covariance Structure
  mu<-rep(0,nMetadata)
  cov<-diag(1,nMetadata, nMetadata)
  
  if (metadataType == 'MVB'){
    for (i in 1:nMetadata){
      for (j in 1:nMetadata){
        if(i!=j) cov[i,j]=0.5**(abs(i-j)) # AR(1)
      }
    }
  }
  
  # Generate from MVN
  fakeMetadata<-as.matrix(MASS::mvrnorm(n=nSamples, mu,cov))
  
  # Transpose and Add Blocking Structure
  finalMetadata<-apply(fakeMetadata, 2, function(x) x+BLOCK)
  
  #############################
  # Modularize Specific Cases #
  #############################
  
  # Multivariable Scenario - Dichotomize Half of the Features
  if (metadataType %in% c('MVA', 'MVB')){
    t_UserMetadata<-apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0))
    columns_not_to_binarize<-sample(1:nMetadata, nMetadata/2)
    t_UserMetadata[,columns_not_to_binarize]<-finalMetadata[, columns_not_to_binarize]
    UserMetadata<-t(t_UserMetadata)
  }
  
  # Univariate Binary
  else if (metadataType == 'UVB'){
    UserMetadata<-t(apply(finalMetadata, 2, function(x) ifelse(x>median(x), 1, 0)))
  } 
  # Univariate Continuous
  else {
    UserMetadata<-t(finalMetadata)
    } 

  # Collect Relevant Spike-in Information
  spikeCount<- round(nMetadata*spikeMetadata)
  Metadatafrozenidx<-sample(1:nMetadata, spikeCount, replace=FALSE)
  significant_metadata<-paste('Metadata', Metadatafrozenidx, sep='')
  spikeCount<-as.character(spikeCount)
  
  # Return 
  return(list(UserMetadata=UserMetadata, 
              Metadatafrozenidx=Metadatafrozenidx,
              significant_metadata=significant_metadata,
              spikeCount=spikeCount))
}


