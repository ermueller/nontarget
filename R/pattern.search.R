pattern.search  <-  function(
        peaklist,
        iso,
        cutint=min(peaklist[,2]),
        rttol=c(-0.5,0.5),
        mztol=3,
        mzfrac=0.1,
        ppm=TRUE,
        inttol=0.5,
        rules=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),
        deter=FALSE,
        entry=20
    ){
        
        ############################################################################
        # (0) Check Inputs #########################################################
        if(mzfrac>1 || mzfrac<=0){ stop("mzfrac must be >0 and <=1") }
        if(mztol<=0){warning("mztol should be >0!")}
        if(inttol>1 || inttol<0){ stop("inttol must be >0 and <=1") }
        if(length(rttol)!=2){stop("rttol must have a lower and an upper bound!")}
        if(rttol[1]>rttol[2]){stop("minimum > maximum for rttol!")}
        if(length(rules)<11){stop("wrong parameter setting: number of rules < 8!")}
        if(!is.data.frame(peaklist)){stop("peaklist must be a data.frame")}
        if(length(peaklist[1,])>3){stop("peaklist with > 3 columns not allowed")}
        if(!nrow(peaklist)>1){stop("peaklist with one entry - doesn`t make sense ...")}
        if(!is.numeric(peaklist[,1]) || !is.numeric(peaklist[,2]) || !is.numeric(peaklist[,3]) ){stop("peaklist columns not numeric")}
        if(rules[4]==TRUE & any(iso$elements=="C")==FALSE & deter!=TRUE){stop("How is rule #7 supposed to work if carbon is not part of the iso argument? Include carbon or set rules[7] to FALSE.")}
        ############################################################################
        
        
        cat("\n (1) Assemble lists ... ")
        # (1) Define parameters / lists / matrices / ...
        # (1.1) Sort peaklist by rt and mz
        peaklist_order <- order(peaklist[,3],peaklist[,1],decreasing=FALSE)
        peaklist <- peaklist[peaklist_order,]
        
        # Retrieve data from isos
        isos <- iso[[1]]
        isomat<-iso[[2]]
        charge_isos <- abs(iso[[3]])
        n_isos <- iso[[4]]
        elements <- iso[[5]]
        
        # If deter.iso wasn't used
        if(deter==FALSE){
            for(i in 1:length(elements)){
                if(any(isos[,1] == elements[i])!=TRUE){
                    stop(paste("Element ",elements[i]," not found in iso[[5]]",sep=""))
                }
            }
        }else{
            # Set rules to FALSE if deter.iso was used
            rules=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
        }
        
        # (1.2) Objects required for the screening step
        
        # (1.2.1)
        # Variables for mass increments
        n_peaks <- nrow(peaklist)
        ID <-  peaklist_order
        
        iso_name <- rep("none",n_peaks)   # (1) which isotope?
        iso_from <- rep("0",n_peaks)      # (2) from which peak?
        iso_to <- rep("0",n_peaks)      # (3) to which peak?
        iso_tolerance <- rep("0",n_peaks)      # (4) within [1] large or [2] small mass tolerance?
        iso_charge <- rep("0",n_peaks)      # (5) with which charge?
        
        # (1.2.2)
        # Variables for mass increment validation
        # Max relative abundance of all isotopes compared to the monoisotope
        max_dabund <- max(isomat[,3])
        
        # Counts for how many 
        countrem1 <- c(0)
        countrem2 <- c(0)
        countrem3 <- c(0)
        countrem4 <- c(0)
        countrem5 <- c(0)
        countrem6 <- c(0)
        countrem7 <- c(0)
        countrem8 <- c(0)
        countrem9 <- c(0)
        countrem11 <- c(0)
        
        # (1.2.3)
        # Variables for grouping and atom number estimation
        groupcount <- c(1)
        group1 <- rep("0",n_peaks)      # which group? per charge level!
        group2 <- rep("0",n_peaks)      # which tree level?
        groupinfo <- rep()           # how many atoms?
        group3 <- c()                # number of group ...
        group4 <- c()                # ... and ID of peaks in that group!
        group5 <- rep("0",n_peaks)      # store charge level
        group6 <- c()                # store charge level
        
        # (2) Use the C++ "mass"-method from massCpp.cpp to screen for mass increments
        cat("\n (2) Screen for mass increments ... ")
        if(ppm==TRUE){
            ppm2=1
        }else{
            ppm2=2
        }
        
        # Since arrays need to be initialized for C++, we need to this in R beforehand.
        isotope_cpp       <- rep(0,nrow(peaklist)*entry)
        iso_from_cpp      <- rep(0,nrow(peaklist)*entry)
        iso_to_cpp        <- rep(0,nrow(peaklist)*entry)
        tolerance_cpp     <- rep(0,nrow(peaklist)*entry)
        charge_cpp        <- rep(0,nrow(peaklist)*entry)
        maxmass <- max(isomat[,2])
        result <- .C("mass",
                   as.double(peaklist[,1]),
                   as.double(peaklist[,3]),
                   as.integer(n_peaks),
                   as.double(mztol*2),
                   as.double(mzfrac*2),
                   as.double(rttol[1]),
                   as.double(rttol[2]),
                   as.integer(n_isos), # 8
                   as.double(isomat[,2]),
                   as.integer(isomat[,4]),
                   as.double(maxmass),
                   as.integer(isomat[,7]), # 12
                   as.integer(entry),
                   as.integer(ppm2), # 14
                   as.integer(isotope_cpp),
                   as.integer(iso_from_cpp),
                   as.integer(iso_to_cpp),
                   as.integer(tolerance_cpp),
                   as.integer(charge_cpp), # 19
                   PACKAGE="nontarget"
        )
        
        # Generate outputs: 
        isomat[,4] <- result[10]
        
        # Go through all peaks
        for(i in 1:(n_peaks)){
            
            # See if the current peak is any isotope or if it has an isotope
            is_isotope <- any(result[15][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0)
            has_isotope <- any(result[16][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0)
            
            # Concatenate relevant peaks/information/charge/etc with a "/" for each peak
            if(is_isotope){
                
                # Which isotope
                iso_name[i] <- paste(iso_name[i],"/",paste0(result[15][[1]][((i-1)*entry+1):((i-1)*entry+entry)][result[15][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0],collapse="/"),sep="")
                
                # To which peak
                iso_to[i] <- paste(iso_to[i],"/",paste0(result[17][[1]][((i-1)*entry+1):((i-1)*entry+entry)][result[17][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0],collapse="/"),sep="")
                
                # Charge level
                iso_charge[i] <- paste(iso_charge[i],"/",paste0(result[19][[1]][((i-1)*entry+1):((i-1)*entry+entry)][result[19][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0],collapse="/"),sep="")
            }

            if(has_isotope){
                
                # From which peak
                iso_from[i] <- paste(iso_from[i],"/",paste0(result[16][[1]][((i-1)*entry+1):((i-1)*entry+entry)][result[16][[1]][((i-1)*entry+1):((i-1)*entry+entry)]!=0],collapse="/"),sep="")
            }

            # Tolerance: Small or large? 
            for(j in 1:entry){
                if(result[18][[1]][(i-1)*entry+j]==1){
                    iso_tolerance[i] <- paste(iso_tolerance[i],"small",sep="/")
                }
                if(result[18][[1]][(i-1)*entry+j]==2){
                    iso_tolerance[i] <- paste(iso_tolerance[i],"large",sep="/")
                }
            }
        }
        
        if(result[13][[1]]!=entry){cat("WARNING: entry overflow -> links missing! Decrease mztol? Increase entry argument?")}
        rm(result)
        if(deter==FALSE){
            cat("\n (3) Check plausibility ... ")
            
            # This function is because the clean-up after every rule is exactly the same
            e <- environment()
            clean_up <- function(){
                
                # We have to do this with environments since we can't set variables outside a function
                # from within a function otherwise
               
                with(e,{ # Behave as if "e" is actually our environment
                    
                    # Rewrite isomat to change the count of the current isotope
                    if(any(!iso_to_keep)){
                        isomat[as.numeric(iso_name_split[!iso_to_keep]),4] <- isomat[as.numeric(iso_name_split[!iso_to_keep]),4]-1
                    }
                    
                    # Rewrite the entries
                    iso_to[i] <- paste(c("0", iso_to_split[iso_to_keep]),sep="/",collapse="/")
                    iso_name[i] <- paste(c("0",iso_name_split[iso_to_keep]),sep="/",collapse="/")
                    iso_tolerance[i] <- paste(c("0",iso_tolerance_split[iso_to_keep]),sep="/",collapse="/")
                    iso_charge[i] <- paste(c("0",iso_charge_split[iso_to_keep]),sep="/",collapse="/")
                    
                    if(iso_name[i] == "0"){
                        iso_name[i] <- "none"
                    }
                })
            }
            
            # (3) Remove invalid dmass-links based on rules (1) to (3)
            ST_2 <- system.time({

            # Pre-split all entries
            iso_to_split_all        <- lapply(iso_to,function(x)strsplit(x,"/")[[1]][-1])
            iso_name_split_all      <- lapply(iso_name,function(x)strsplit(x,"/")[[1]][-1])
            iso_tolerance_split_all <- lapply(iso_tolerance,function(x)strsplit(x,"/")[[1]][-1])
            iso_charge_split_all    <- lapply(iso_charge,function(x)as.numeric(strsplit(x,"/")[[1]][-1]))
            iso_to_keep <- list()
            
            # Go through all indices where tentatively isotopes have been found
            for(i in .Internal(which(iso_to != "0"))){
                
                    # Split the current values by "/"
                    iso_to_split <- iso_to_split_all[[i]]
                    iso_name_split <- iso_name_split_all[[i]]
                    iso_tolerance_split <- iso_tolerance_split_all[[i]]
                    iso_charge_split <- iso_charge_split_all[[i]]
                    iso_to_keep <- rep(TRUE,length(iso_to_split))
                    
                    # (3.1) Rule 1: intensity ratio check over ALL isotopes at charge level ###
                    if(rules[1]=="TRUE"){
                        
                        # Go through all split elements
                        for(j in 1:length(iso_to_split)){
                            
                            # Exclude peaks from this check with more than one atom distance
                            if(isomat[as.numeric(iso_name_split[j]),5] == 1){
                                
                                # Calculate the minimum possible mass
                                to_peak_expected_int_min      <- peaklist[as.numeric(iso_to_split[j]),2]*(1-inttol)
                                current_peak_expected_int_max <- peaklist[i,2]*(1+inttol)
                                
                                possnumb <- (to_peak_expected_int_min/current_peak_expected_int_max) * max_dabund # Based on isotope with largest relative abundance
                                
                                if(possnumb < 1){
                                    possnumb <- 1
                                }
                                
                                # Largest charge z, minimum mass of hydrogen
                                if((possnumb*3.0078)>(peaklist[i,1]*max(iso_charge_split))){ 
                                    iso_to_keep[j] <- FALSE
                                    countrem1 <- c(countrem1+1)
                                }
                            }
                        }
                        
                        # remove entries in iso_to[i] ("to"), iso_from[->i] ("from") and isomat!
                        iso_to_keep[isomat[as.numeric(iso_name_split),5] != 1] <- TRUE # if double-distanced
                        
                        for(m in 1:length(iso_to_keep)){
                            if(iso_to_keep[m]==FALSE){
                                iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])
                            }
                        }
                        
                        clean_up()
                    } # RULE1
                    
                    #
                    ##########################################################################
                    # (3.2) RULE2: intensity ratio check with LARGE mass tolerance ###########
                    if(any(iso_to_keep)){ # anything left to check?
                        if(rules[2]=="TRUE"){
                            iso_to_split <- c(strsplit(iso_to[i],"/")[[1]][-1])
                            iso_to_keep <- rep(TRUE,length(iso_to_split))
                            iso_name_split <- c(strsplit(iso_name[i],"/")[[1]][-1])
                            that6 <- c(iso_to_split[duplicated(iso_to_split)==FALSE]) 
                            that7 <- c(isomat[as.numeric(iso_name_split),5]==1)
                            #iso_tolerance_split <- c(strsplit(iso_tolerance[i],"/")[[1]][-1])
                            iso_charge_split <- c(as.numeric(strsplit(iso_charge[i],"/")[[1]][-1]))
                            for(j in 1:length(that6)){
                                if(any(that7[iso_to_split==that6[j]])){ # only on single-distanced peaks!
                                    max_dabund2 <- min(1/isomat[as.numeric(iso_name_split[iso_to_split==that6[j]&that7==TRUE]),3])
                                    ratioC <- max(isomat[as.numeric(iso_name_split[iso_to_split==that6[j]&that7==TRUE]),6])
                                    possnumb <- min(((peaklist[as.numeric(iso_to_split[iso_to_split==that6[j]]),2]*(1-inttol))/(peaklist[i,2]*(1+inttol)))*max_dabund2)
                                    if(possnumb<1){possnumb <- c(1)}
                                    possmass <- min(isos[as.logical(match(isos[,1],isos[as.logical(match(isos[,2],isomat[as.numeric(iso_name_split[iso_to_split==that6[j]&that7==TRUE]),1],nomatch=FALSE)),1],nomatch=FALSE)),3])
                                    if(ratioC!=0){
                                        if((possnumb*possmass+((possnumb/ratioC)*12))>(peaklist[i,1]*max(isomat[as.numeric(iso_name_split[iso_to_split==that6[j]&that7==TRUE]),7]))){ # include ratios to C! BEWARE if ration=0->Inf->alway>mass
                                            iso_to_keep[as.numeric(iso_to_split)==as.numeric(that6[j])] <- FALSE
                                            countrem2 <- c(countrem2+1)
                                        } # if ...
                                    }else{
                                        if((possnumb*possmass)>(peaklist[i,1]*max(isomat[as.numeric(iso_name_split[iso_to_split==that6[j]&that7==TRUE]),7]))){ # include ratio to C
                                            iso_to_keep[iso_to_split==that6[j]] <- FALSE
                                            countrem2 <- c(countrem2+1)
                                        } # if ...
                                    }
                                } #
                            } # for
                            # remove entries in iso_from[i] (i.e."to") only  / NOT "from!"
                            # reset: keep non-single distanced peaks!
                            iso_to_keep[isomat[as.numeric(iso_name_split),5]!=1] <- TRUE
                            for(m in 1:length(iso_to_keep)){
                                if(iso_to_keep[m]==FALSE){
                                    iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])
                                }
                            }
                            clean_up()
                        } # if
                    } # RULE2
                    #
                    ##########################################################################
                    # (3.3) RULE3: intensity ratio check with SMALL mass tolerance ###########
                    if(rules[3]=="TRUE"){
                        if(any(iso_to_keep)){# anything left to check?
                            iso_to_split <- c(strsplit(iso_to[i],"/")[[1]][-1])
                            iso_to_keep <- rep(TRUE,length(iso_to_split))
                            iso_name_split <- c(strsplit(iso_name[i],"/")[[1]][-1])
                            that7 <- c(isomat[as.numeric(iso_name_split),5]==1)
                            iso_tolerance_split <- c(strsplit(iso_tolerance[i],"/")[[1]][-1])
                            iso_charge_split <- c(as.numeric(strsplit(iso_charge[i],"/")[[1]][-1]))
                            for(j in 1:length(iso_to_split)){
                                if( iso_tolerance_split[j]=="small" && that7[j]==TRUE ){ # within small mass tolerance? single distanced?
                                    max_dabund2 <- c(1/isomat[as.numeric(iso_name_split[j]),3])
                                    ratioC <- c(isomat[as.numeric(iso_name_split[j]),6])
                                    possnumb <- c(((peaklist[as.numeric(iso_to_split[j]),2]*(1-inttol))/(peaklist[i,2]*(1+inttol)))*max_dabund2)
                                    if(possnumb<1){possnumb <- c(1)}
                                    possmass <- min(isos[as.logical(match(isos[,1],isos[as.logical(match(isos[,2],isomat[as.numeric(iso_name_split[j]),1],nomatch=FALSE)),1],nomatch=FALSE)),3])
                                    if(ratioC!=0){
                                        if((possnumb*possmass+((possnumb/ratioC)*12))>(peaklist[i,1]*c(isomat[as.numeric(iso_name_split[j]),7]))){ # include ratios to C! BEWARE if ration=0->Inf->alway>mass
                                            iso_to_keep[j] <- FALSE
                                        } # if ...
                                    }else{
                                        if((possnumb*possmass)>(peaklist[i,1]*c(isomat[as.numeric(iso_name_split[j]),7]))){ # include ratio to C
                                            iso_to_keep[j] <- FALSE
                                        } # if ...
                                    }
                                }    # if with small mass range
                            }    # for all possible isotopes
                            # remove entries in iso_from[i] (i.e."to") only
                            # reset: keep non-single distanced peaks!
                            iso_to_keep[isomat[as.numeric(iso_name_split),5]!=1] <- TRUE
                            for(m in 1:length(iso_to_keep)){ 
                                if(iso_to_keep[m]==FALSE){
                                    #if(length(iso_to_split[iso_to_split==iso_to_split[m]])==1){ # if several hits: do not remove entry but set to "large"
                                    iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])
                                    countrem3 <- c(countrem3+1)
                                    #}else{
                                    #iso_to_keep[m] <- TRUE
                                    #iso_tolerance_split[m] <- "large"
                                    #countrem3 <- c(countrem3+1)
                                    #}
                                }
                            }
                            clean_up()
                        }# if any in iso_to_keep
                    }# rule 3
                    ##########################################################################
                    # (3.4) RULE 4: minimum intensity for all used isotopes reached? #########
                    if(rules[4]=="TRUE"){
                        if(any(iso_to_keep)){# anything left to check?
                            iso_to_split <- c(strsplit(iso_to[i],"/")[[1]][-1])
                            iso_to_keep <- rep(TRUE,length(iso_to_split))
                            iso_name_split <- c(strsplit(iso_name[i],"/")[[1]][-1])
                            iso_tolerance_split <- c(strsplit(iso_tolerance[i],"/")[[1]][-1])
                            iso_charge_split <- c(as.numeric(strsplit(iso_charge[i],"/")[[1]][-1]))
                            for(j in 1:length(iso_to_split)){
                                if(isomat[as.numeric(iso_name_split[j]),5]==1){ # exclude double-distanced peaks
                                    if(((peaklist[as.numeric(iso_to_split[j]),2]+(inttol*peaklist[as.numeric(iso_to_split[j]),2]))/(peaklist[i,2]-(inttol*peaklist[i,2])))<min(isomat[,3])){
                                        iso_to_keep[j] <- FALSE
                                        countrem4 <- c(countrem4+1)        
                                    }
                                }
                            }
                            # remove entries in iso_to[i] ("to"), iso_from[->i] ("from") and isomat!
                            iso_to_keep[isomat[as.numeric(iso_name_split),5]!=1] <- TRUE # if double-distanced
                            for(m in 1:length(iso_to_keep)){
                                if(iso_to_keep[m]==FALSE){
                                    iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])
                                }
                            }
                            clean_up()
                        }# if any in iso_to_keep
                    } # rule 4
                    ##########################################################################
                    # (3.5) RULE 5: minimum intensity for isotopes within large reached? #####
                    if(rules[5]=="TRUE"){
                        if(any(iso_to_keep)){ # anything left to check?
                            iso_to_split <- c(strsplit(iso_to[i],"/")[[1]][-1])
                            iso_to_keep <- rep(TRUE,length(iso_to_split))
                            iso_name_split <- c(strsplit(iso_name[i],"/")[[1]][-1])
                            iso_tolerance_split <- c(strsplit(iso_tolerance[i],"/")[[1]][-1])
                            iso_charge_split <- c(as.numeric(strsplit(iso_charge[i],"/")[[1]][-1]))
                            for(j in 1:length(iso_to_split)){
                                if(isomat[as.numeric(iso_name_split[j]),5]==1){ # exclude double-distanced peaks
                                    if(
                                        ((peaklist[as.numeric(iso_to_split[j]),2]*(1+inttol))/(peaklist[i,2]*(1-inttol)))
                                        < min(isomat[as.numeric(iso_name_split[iso_to_split==iso_to_split[j]]),3])
                                    ){
                                        iso_to_keep[j] <- FALSE
                                        countrem5 <- c(countrem5+1)        
                                    }
                                }
                            }
                            # remove entries in iso_to[i] ("to"), iso_from[->i] ("from") and isomat!
                            iso_to_keep[isomat[as.numeric(iso_name_split),5]!=1] <- TRUE # if double-distanced
                            for(m in 1:length(iso_to_keep)){
                                if(iso_to_keep[m]==FALSE){
                                    iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])
                                }
                            }
                            clean_up()
                        }# if any in iso_to_keep
                    } # rule 5
                    ##########################################################################
                    # (3.6) RULE 6: minimum intensity for isotopes within small reached? #####
                    if(rules[6]=="TRUE"){
                        if(any(iso_to_keep)){# anything left to check?
                            iso_to_split <- c(strsplit(iso_to[i],"/")[[1]][-1])
                            iso_to_keep <- rep(TRUE,length(iso_to_split))
                            iso_name_split <- c(strsplit(iso_name[i],"/")[[1]][-1])
                            that7 <- c(isomat[as.numeric(iso_name_split),5]==1)
                            iso_tolerance_split <- c(strsplit(iso_tolerance[i],"/")[[1]][-1])
                            iso_charge_split <- c(as.numeric(strsplit(iso_charge[i],"/")[[1]][-1]))
                            for(j in 1:length(iso_to_split)){
                                if( iso_tolerance_split[j]=="small" && that7[j]==TRUE ){
                                    if(
                                        ((peaklist[as.numeric(iso_to_split[j]),2]*(1+inttol))/(peaklist[i,2]*(1-inttol)))
                                        < (isomat[as.numeric(iso_name_split[j]),3])
                                    ){
                                        iso_to_keep[j] <- FALSE
                                        countrem6 <- c(countrem6+1)        
                                    }
                                }
                            }
                            # remove entries in iso_from[i] (i.e."to") only
                            # reset: keep non-single distanced peaks!
                            iso_to_keep[isomat[as.numeric(iso_name_split),5]!=1] <- TRUE
                            for(m in 1:length(iso_to_keep)){ 
                                if(iso_to_keep[m]==FALSE){ 
                                    #if(length(iso_to_split[iso_to_split==iso_to_split[m]])==1){ # if several hits: do not remove entry but set to "large"
                                    iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])
                                    countrem6 <- c(countrem6+1)
                                    #}else{
                                    #iso_to_keep[m] <- TRUE
                                    #iso_tolerance_split[m] <- "large"
                                    #countrem6 <- c(countrem6+1)
                                    #}
                                }
                            }
                            clean_up()
                        } # if any in iso_to_keep
                    } # rule 6
                    ##########################################################################
                    # (3.7) RULE 7: carbon ratio check with SMALL mass tolerance #############
                    if(rules[7]=="TRUE"){
                        if(any(iso_to_keep)){ # anything left to check?
                            iso_to_split <- c(strsplit(iso_to[i],"/")[[1]][-1])
                            iso_to_keep <- rep(TRUE,length(iso_to_split))
                            iso_name_split <- c(strsplit(iso_name[i],"/")[[1]][-1])
                            that7 <- c(isomat[as.numeric(iso_name_split),5]==1)
                            iso_tolerance_split <- c(strsplit(iso_tolerance[i],"/")[[1]][-1])
                            iso_charge_split <- c(as.numeric(strsplit(iso_charge[i],"/")[[1]][-1]))
                            for(j in 1:length(iso_to_split)){
                                if(iso_tolerance_split[j]=="small" && that7[j]==TRUE && isomat[as.numeric(iso_name_split[j]),1]!="13C" && isomat[as.numeric(iso_name_split[j]),6]!=0 ){
                                    max_dabund2 <- c(1/isomat[as.numeric(iso_name_split[j]),3])
                                    ratioC <- c(isomat[as.numeric(iso_name_split[j]),6])
                                    possnumb <- c(((peaklist[as.numeric(iso_to_split[j]),2]*(1-inttol))/(peaklist[i,2]*(1+inttol)))*max_dabund2)
                                    if(possnumb<1){possnumb <- c(1)} # must have at least one atom for a peak signal
                                    numbC <- c(possnumb/ratioC)
                                    if(numbC<1){numbC <- c(1)} # must have at least one carbon atom!
                                    if(peaklist[i,2]*isomat[isomat[,1]=="13C",3][1]*(1/numbC)>cutint){
                                        if(rules[10]==FALSE){ # beware of interfering peaks?
                                            if(any(isomat[as.numeric(iso_name_split),1]=="13C")==FALSE){
                                                iso_to_keep[j] <- FALSE
                                                countrem7 <- c(countrem7+1)
                                            }
                                        }else{
                                            if(
                                                any(
                                                    peaklist[,1]>=(peaklist[i,1]+(1.003355/isomat[as.numeric(iso_name_split[j]),7])-0.5) &
                                                    peaklist[,1]<=(peaklist[i,1]+(1.003355/isomat[as.numeric(iso_name_split[j]),7])+0.5) &
                                                    peaklist[,3]>=(peaklist[i,3]+rttol[1]) &
                                                    peaklist[,3]<=(peaklist[i,3]+rttol[2])
                                                )==FALSE
                                            ){
                                                iso_to_keep[j] <- FALSE
                                                countrem7 <- c(countrem7+1)
                                            }
                                        }
                                    }
                                }
                            }
                            # remove entries in iso_from[i] (i.e."to") only !!!
                            # reset: keep non-single distanced peaks!
                            iso_to_keep[isomat[as.numeric(iso_name_split),5]!=1] <- TRUE
                            for(m in 1:length(iso_to_keep)){if(iso_to_keep[m]==FALSE){iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])}}
                            clean_up()
                        }# if any in iso_to_keep
                    }
                } # for
            })[3]
            #it <- data.frame(ID,getit4,getit2,getit1,getit5,getit6)
            #names(it) <- c("ID","to","from","isotope","tol","charge")
            cat("done.")
        }else{
            cat("\n (3) Plausibility tests skipped. ")
        } # if deter == FALSE
  
        ############################################################################
        cat("\n (4) Group peaks within charge levels: ")
        # (4) group! ###############################################################
        ST_3 <- system.time({
        along <- order(peaklist[,1],decreasing=FALSE)
        for(z in 1:length(charge_isos)){
            group1b <- rep("0",n_peaks)      # which group? per charge level! renew per charge level!
            group2b <- rep("0",n_peaks)      # which interaction level?
            group5b <- rep("0",n_peaks)      # ... and which charge?
            i <- c(1)
            while(i<n_peaks){
                # correct entry, if peak points at itself! ###############################
                these1 <- c(as.numeric(strsplit(iso_to[along[i]],"/")[[1]][-1]))
                if(any(these1==along[i])){ # remove any self-reference
                    these1 <- these1[these1!=along[i]]
                    if(length(these1)==0){
                        iso_name[along[i]] <- "none"
                        iso_to[along[i]] <- "0"
                        iso_tolerance[along[i]] <- "0"
                        iso_charge[along[i]] <- "0"
                    }else{
                        these1 <- strsplit(iso_name[along[i]],"/")[[1]][-1]
                        these4 <- as.numeric(strsplit(iso_to[i],"/")[[1]][-1])
                        these5 <- strsplit(iso_tolerance[along[i]],"/")[[1]][-1]
                        these6 <- strsplit(iso_charge[along[i]],"/")[[1]][-1]
                        these1 <- these1[these4!=along[i]]
                        these5 <- these5[these4!=along[i]]
                        these6 <- these6[these4!=along[i]]
                        these4 <- these4[these4!=along[i]]
                        iso_name[along[i]] <- "none"
                        iso_to[along[i]] <- "0"
                        iso_tolerance[along[i]] <- "0"
                        iso_charge[along[i]] <- "0"
                        # collapes by paste, not loop
                        for(j in 1:length(these1)){iso_name[along[i]] <- paste(iso_name[along[i]],"/",these1[j],sep="")}
                        for(j in 1:length(these4)){iso_to[along[i]] <- paste(iso_to[along[i]],"/",these4[j],sep="")}
                        for(j in 1:length(these5)){iso_tolerance[along[i]] <- paste(iso_tolerance[along[i]],"/",these5[j],sep="")}
                        for(j in 1:length(these6)){iso_charge[along[i]] <- paste(iso_charge[along[i]],"/",these6[j],sep="")}
                    }
                }
                if( (iso_to[along[i]]!="0") & (group1b[along[i]]==0) & (grepl(as.character(charge_isos[z]),iso_charge[along[i]])) ){  # group1b: schon als M+X erfasst
                    ########################################################################
                    these1 <- c(along[i],as.numeric(strsplit(iso_to[along[i]],"/")[[1]][-1]))
                    # these1 <- c(i,as.numeric(strsplit(iso_to[i],"/")[[1]][-1]))
                    these5 <- as.numeric(strsplit(iso_charge[along[i]],"/")[[1]])[-1]
                    these1 <- these1[c(TRUE,these5==charge_isos[z])]
                    these1 <- as.numeric(unique(these1)) # remove double entries
                    group2b[along[i]] <- c("1/0")
                    group2b[these1[these1!=along[i]]] <- paste("2",group2b[these1[these1!=along[i]]],sep="/")
                    allpeaks <- these1
                    if(length(these1)>1){
                        newpeaks1 <- these1[these1!=along[i]]
                    }else{
                        newpeaks1 <- c()
                    }
                    level <- c(3)
                    while(length(newpeaks1)>0){
                        newlevel <- c()
                        newpeaks2 <- c()
                        for(m in 1:length(newpeaks1)){
                            these1 <- c(as.numeric(strsplit(iso_to[newpeaks1[m]],"/")[[1]][-1]))
                            these5 <- as.numeric(strsplit(iso_charge[newpeaks1[m]],"/")[[1]])[-1]
                            these1 <- c(these1[these5==charge_isos[z]])
                            if(length(these1)>0){
                                for(n in 1:length(these1)){
                                    if(any(these1[n]==allpeaks)!=TRUE){
                                        newpeaks2 <- c(newpeaks2,these1[n])
                                        newlevel <- c(newlevel,level)
                                    } # if
                                } # for
                            } # if
                        }
                        # remove double entries in [newpeaks2, newlevel] #####################
                        if(length(newpeaks2)>1){
                            bad <- c()
                            unt <- length(newpeaks2)
                            for(l in 1:(unt-1)){
                                if( any(  (newpeaks2[l]==newpeaks2[(l+1):unt]) & (newlevel[l]==newlevel[(l+1):unt])  ) ){bad <- c(bad,l)}
                            }
                            if(length(bad)>0){
                                newpeaks2 <- newpeaks2[-bad]
                                newlevel <- newlevel[-bad]
                            }
                        }
                        if(length(newpeaks2)>0){
                            for(l in 1:length(newpeaks2)){
                                group2b[newpeaks2[l]] <- paste(newlevel[l],group2b[newpeaks2[l]],sep="/")
                            }
                        }
                        # clean for new round ################################################
                        allpeaks <- c(allpeaks,newpeaks1,newpeaks2)
                        allpeaks <- as.numeric(unique(allpeaks))
                        newpeaks1 <- newpeaks2
                        newpeaks1 <- as.numeric(unique(newpeaks1))
                        level <- c(level+1)
                    } # while
                    these1 <- allpeaks
                    these1 <- these1[these1!=0] # dispensable?
                    ########################################################################
                    # RULE8: calculate feasible mass range -> apply as filter after grouping
                    # shift to these 2, keep these 1 for rules[9]
                    # these2: defines upper bound
                    if(rules[8]=="TRUE"){
                        topint <- peaklist[along[i],2]
                        topmass <- peaklist[along[i],1]
                        topcount <- ceiling(topmass/max(isomat[,6]))
                        topput <- c(max(isomat[,2])/charge_isos[z])
                        max_dabund <- min(1/isomat[,3])
                        toprep <- c(1)
                        while(topint>0 && topcount>0){ #topint>cutint->what if these[1]<cutint? 
                            topmass <- c(topmass+topput)
                            topint <- c(topint*(1/max_dabund)*(topcount/toprep))
                            toprep <- c(toprep+1)
                            topcount <- c(topcount-1)
                        }
                        getit <- c(peaklist[these1,1]<=topmass)
                        if(any(getit==FALSE)){countrem8 <- countrem8+1}
                        getit[c(1,2)] <- TRUE
                        these2 <- these1[getit]
                    }else{
                        these2 <- these1
                    }
                    ########################################################################
                    # RULE6 + RULE7: pattern plausibility ##################################
                    plaus <- TRUE # per se before Rule 6 is evaluated!
                    if(rules[9]=="TRUE"){
                        dat2 <- peaklist[these1,]
                        dat4 <- seq(1:length(these1))
                        these4 <- these1
                        dat4 <- dat4[order(dat2[,1],decreasing=FALSE)]
                        these4 <- these4[order(dat2[,1],decreasing=FALSE)]
                        dat2 <- dat2[order(dat2[,1],decreasing=FALSE),]
                        monomass <- dat2[1,1]
                        monointens <- dat2[1,2]
                        this8 <- as.numeric(strsplit(iso_name[along[i]],"/")[[1]])[-1]                # this isotope in "isomat" ...
                        this8b <- c(strsplit(iso_tolerance[along[i]],"/")[[1]])[-1]
                        this8c <- as.numeric(strsplit(iso_charge[along[i]],"/")[[1]])[-1]
                        this10 <- as.numeric(strsplit(iso_to[along[i]],"/")[[1]])[-1]               # ... for this daughter peak ID
                        this11 <- seq(1:length(this8))
                        this8 <- this8[this8b=="small" & this8c==charge_isos[z]]
                        this10 <- this10[this8b=="small" & this8c==charge_isos[z]]
                        this11 <- this11[this8b=="small" & this8c==charge_isos[z]]
                        this9 <- c()                                                         # ... with this entry in peak group "dat2"
                        for(y in 1:length(this10)){
                            this9 <- c(this9,dat4[these4==this10[y]])
                        }
                        if(length(this8)>0){
                            if(rules[10]==FALSE){                        # problem: e.g. 37Cl shading 2*13C
                                if(ppm==TRUE){
                                    mztol2 <- c((mztol*dat2[1,1]/1e6))
                                }else{
                                    mztol2 <- mztol
                                }
                            }else{
                                mztol2 <- 0.5
                            }
                            plaus <- rep(TRUE,length(this8))
                            for(l in 1:length(this8)){ # over all matches found in "isomat"
                                # get number of atoms per dmass
                                numb <- floor((dat2[this9[l],2]*(1-inttol))/(monointens*(1+inttol)*isomat[this8[l],3]))
                                # more such peaks expected?
                                if(numb>2){
                                    mass3 <- c(dat2[1,1]+isomat[this8[l],2]) # start with M not M+1 mass!
                                    int3 <- c(dat2[this9[l],2])
                                    count <- c(2)
                                    while((int3[length(int3)]*(1-inttol))>=(cutint*(1+inttol)) & count<numb){
                                        mass3 <- c(mass3,(mass3[length(mass3)]+isomat[this8[l],2]))
                                        int3 <- c(int3,int3[length(int3)]*isomat[this8[l],3]*((numb-count+1))/((count)))
                                        count <- c(count+1)
                                    }
                                    if(length(mass3)>1){
                                        mass3 <- mass3[-length(mass3)]
                                        int3 <- int3[-length(int3)]
                                    }
                                    for(j in 1:length(mass3)){ # do not check intensities - only if masses exist
                                        if(any( dat2[,1]<=(mass3[j]+(mztol2)) & dat2[,1]>=(mass3[j]-(mztol2)) )){      # code can be improved: search at generation level + link!
                                            plaus[l] <- TRUE
                                        }else{
                                            plaus[l] <- FALSE
                                        }
                                    } # with mztol*3!
                                } # if still plausible ...
                            } # for ...
                        } # if ...
                    } # on RULE 6
                    # to have more rules inserted, the following is separated from RULE 6:
                    if(all(plaus)){
                        group1b[these2] <- paste(groupcount,group1b[these2],sep="/")
                        group5b[these2] <- paste(charge_isos[z],group5b[these2],sep="/")
                        # estimate number of atoms! ##########################################
                        these16 <- c()
                        if(deter==FALSE){
                            for(b in 1:length(these2)){
                                iso_to_split <- c(strsplit(iso_to[these2[b]],"/")[[1]][-1])
                                iso_name_split <- c(strsplit(iso_name[these2[b]],"/")[[1]][-1])
                                iso_tolerance_split <- c(strsplit(iso_tolerance[these2[b]],"/")[[1]][-1])
                                iso_charge_split <- c(as.numeric(strsplit(iso_charge[these2[b]],"/")[[1]][-1]))
                                if(any(iso_tolerance_split=="small")){
                                    for(d in 1:length(iso_tolerance_split)){
                                        if(iso_tolerance_split[d]=="small" && iso_charge_split[d]==charge_isos[z]){
                                            count <- floor((peaklist[as.numeric(iso_to_split[d]),2]*(1-inttol))/(peaklist[as.numeric(these2[b]),2]*(1+inttol)*isomat[as.numeric(iso_name_split[d]),3]))
                                            if(count<1){
                                                count <- c(1)
                                            } # at least one atom for a peak
                                            these16 <- paste(these16,"/",isos[isos[,2]==isomat[as.numeric(iso_name_split[d]),1],1],":",count,sep="")
                                        }
                                    }
                                }
                            } # for b
                        } # deter
                        groupinfo <- c(groupinfo,paste("minimum atom counts: no information",these16,sep=""))
                        group3 <- c(group3,groupcount)
                        group6 <- c(group6,charge_isos[z])
                        iso_to_split6 <- c(these2[1])
                        for(f in 2:length(these2)){
                            iso_to_split6 <- paste(iso_to_split6,",",these2[f],sep="")
                        }
                        group4 <- c(group4,iso_to_split6)
                        groupcount <- c(groupcount+1)
                        i <- c(i+1)				
                    }else{ # if not everything is plausible ...
                        countrem9 <- c(countrem9+1)
                        # remove affected links:
                        this11 <- this11[plaus==FALSE]       # from rule 6!
                        this8 <- this8[plaus==FALSE]         # from rule 6!
                        if(length(this11)==length(strsplit(iso_name[along[i]],"/")[[1]][-1])){
                            iso_name[along[i]] <- "none"
                            iso_to[along[i]] <- "0"
                            iso_tolerance[along[i]] <- "0"
                            iso_charge[along[i]] <- "0"
                            isomat[this8,4] <- isomat[this8,4]-1
                            group2b[these1] <- "0"
                            i <- c(i+1)
                        }else{
                            split_iso <- strsplit(iso_name[along[i]],"/")[[1]][-1][-this11]
                            iso_name[along[i]] <- "0"
                            for(y in 1:length(split_iso)){iso_name[along[i]] <- paste(iso_name[along[i]],"/",split_iso[y],sep="")}
                            split_to <- strsplit(iso_to[along[i]],"/")[[1]][-1][-this11]
                            iso_to[along[i]] <- "0"
                            for(y in 1:length(split_to)){iso_to[along[i]] <- paste(iso_to[along[i]],"/",split_to[y],sep="")}
                            split_tolerance <- strsplit(iso_tolerance[along[i]],"/")[[1]][-1][-this11]
                            iso_tolerance[along[i]] <- "0"
                            for(y in 1:length(split_tolerance)){iso_tolerance[along[i]] <- paste(iso_tolerance[along[i]],"/",split_tolerance[y],sep="")}
                            split_charge <- strsplit(iso_charge[along[i]],"/")[[1]][-1][-this11]
                            iso_charge[along[i]] <- "0"
                            for(y in 1:length(split_charge)){iso_charge[along[i]] <- paste(iso_charge[along[i]],"/",split_charge[y],sep="")}
                            isomat[this8,4] <- isomat[this8,4]-1
                            group2b[these1] <- "0"
                            if(i>1){ i <- c(i-1)} # still peaks left for grouping? reset to previous peak. 
                        }
                    }
                }else{
                    i <- c(i+1)
                } # if conditions 1-3
                ##########################################################################
            } # while i
            ############################################################################
            # merge results from different charge levels! ##############################
            for(x in 1:n_peaks){ 
                if(group1b[x]!="0"){
                    group1[x] <- paste(group1[x],group1b[x],sep="/")
                    group1[x] <- sub("/0/","/",group1[x])
                    group2[x] <- paste(group2[x],group2b[x],sep="/")
                    group2[x] <- sub("/0/","/",group2[x])
                    group5[x] <- paste(group5[x],group5b[x],sep="/")
                    group5[x] <- sub("/0/","/",group5[x])
                }
            }
            cat(paste(charge_isos[z],"/",sep=""))
        } # for z = charge level ###################################################
        ############################################################################
        for(x in 1:n_peaks){ 
            if(group1[x]!="0"){
                group1[x] <- substr(group1[x],1,nchar(group1[x])-2)
                group2[x] <- substr(group2[x],1,nchar(group2[x])-2)
                group5[x] <- substr(group5[x],1,nchar(group5[x])-2)
            }
        }
        #data.frame(group1,group2,group5)
        ############################################################################
        ############################################################################
        # Rule 11: remove nested groups = merge with larger / equal sized group #####
        if(rules[11]==TRUE && length(charge_isos)>1){
            removals1 <- c()
            removals2 <- c()
            removals3 <- c()
            removals4 <- c()
            countrem11 <- c(0)
            for(i in 1:n_peaks){
                if(group1[i]!="0"){
                    this1 <- as.numeric(strsplit(group1[i],"/")[[1]][-1])
                    if(length(this1)>1){
                        for(n in 1:(length(this1)-1)){
                            for(m in (n+1):length(this1)){
                                if(group6[this1[n]]!=group6[this1[m]]){ # not on same charge level (anyway impossible after grouping)!
                                    this2 <- as.numeric(strsplit(group4[this1[n]],",")[[1]])
                                    this3 <- as.numeric(strsplit(group4[this1[m]],",")[[1]])
                                    if( any(is.na(match(this2,this3)))!=TRUE ){
                                        removals1 <- c(removals1,this1[n]) # keep
                                        removals2 <- c(removals2,this1[m]) # bed into
                                        removals3 <- c(removals3,i)
                                        removals4 <- c(removals4,group6[this1[n]])
                                    }
                                    if( length(this2)!=length(this3) ){ # otherwise done above
                                        if( any(is.na(match(this3,this2)))!=TRUE ){
                                            removals1 <- c(removals1,this1[m])
                                            removals2 <- c(removals2,this1[n])
                                            removals3 <- c(removals3,i)
                                            removals4 <- c(removals4,group6[this1[m]])
                                        }
                                    }
                                }
                            }
                        }
                    }
                } # if
            } # for
            # remove double entries ...
            bad <- c()
            unt <- length(removals1)
            for(i in 1:(length(removals1)-1)){
                if(any(removals1[i]==removals1[(i+1):unt]) & any(removals2[i]==removals2[(i+1):unt])){ bad <- c(bad,i)}
            }
            if(length(bad)>0){
                removals1 <- removals1[-bad]   # group to be merged
                removals2 <- removals2[-bad]   # to be kept
                removals3 <- removals3[-bad]   # i
                removals4 <- removals4[-bad]   # charge level to be merged
                countrem11 <- length(removals4) #
                #data.frame(removals1,removals2,removals4,removals3)
                # correct group1 and group6: merge and delete!
                this1 <- rep(TRUE,length(group4))
                for(i in 1:length(removals3)){
                    group3[removals2[i]] <- paste(group3[removals2[i]],removals1[i],sep="/")
                    group6[removals2[i]] <- paste(group6[removals2[i]],removals4[i],sep="/")
                    this1[removals1[i]] <- FALSE
                }
                group3 <- group3[this1]
                group4 <- group4[this1]
                group6 <- group6[this1]
            }else{
                countrem11 <- c(0)
            } # if bad>0
        } # rule 8
        # data.frame(group3,group4,group6)
        # data.frame(ID,getit4,getit1,getit5,getit6)
        cat("done.")
        ############################################################################
        })[3]
        ############################################################################
        cat("\n (5) Create output... ")
        ############################################################################
        overlap <- rep(0,100)
        for(i in 1:n_peaks){
            if(group1[i]!="0"){
                this11 <- strsplit(group1[i],"/")[[1]]
                this11 <- this11[-length(this11)]
                overlap[length(this11)] <- c(overlap[length(this11)]+1)
            }
        }
        overlap <- overlap[overlap!=0]
        if(length(overlap)>0){
            this11 <- data.frame(seq(1:length(overlap)),overlap,stringsAsFactors=FALSE)
            names(this11) <- c("Number of groups in overlap","Peak counts")
        }else{
            this11 <- "No overlaps detected"
        }
        ############################################################################
        deep <- rep(0,1000)
        for(i in 1:length(group2)){
            if(group2[i]!="0"){deep[as.numeric(strsplit(group2[i],"/")[[1]])] <- c((deep[as.numeric(strsplit(group2[i],"/")[[1]])])+1)
            }}
        deep <- deep[deep!=0]
        if(length(deep)>0){
            deep <- data.frame(seq(1:length(deep)),deep,stringsAsFactors=FALSE)
            names(deep) <- c("interaction level","peak counts")
        }else{
            deep <- "No groups formed"
        }
        ############################################################################
        hits <- data.frame(isomat[,c(1,7,4)],rep(0,length(isomat[,1])),rep("0",length(isomat[,1])),stringsAsFactors=FALSE)
        names(hits) <- c("isotope","charge","peak counts","group counts","element")
        # increment counts
        for(j in 1:length(iso_name)){
            if(iso_name[j]!="none"){
                this1 <- as.numeric(strsplit(iso_name[j],"/")[[1]][-1])
                for(n in 1:length(this1)){
                    hits[this1[n],3] <- c( hits[this1[n],3]+ 1)
                }
            }
        }
        # group counts
        if(length(group4)>0){ # anything found at all?
            for(j in 1:length(group4)){
                hit <- c()
                this <- as.numeric(strsplit(as.character(group4[j]),",")[[1]])
                that <- as.numeric(strsplit(as.character(group6[j]),"/")[[1]])
                for(n in 1:length(this)){
                    this1 <- as.numeric(strsplit(iso_name[this[n]],"/")[[1]][-1])
                    this2 <- as.numeric(strsplit(iso_to[this[n]],"/")[[1]][-1])
                    this3 <- as.numeric(strsplit(iso_charge[this[n]],"/")[[1]][-1])
                    if(length(this1)>0){
                        for(m in 1:length(this1)){
                            if(any(this==this2[m]) & any(that==this3[m])){
                                hit <- c(hit,this1[m])
                            }
                        }
                    }
                }
                hit <- as.numeric(unique(hit))
                hits[hit,4] <- c( hits[hit,4]+ 1)
            }
            hits <- hits[order((hits[,1]),(hits[,2]),decreasing=FALSE),]
            hits[,5] <- as.character(hits[,5])
        }
        if(deter==FALSE){
            hits[,5] <- as.character(hits[,5])
            for(i in 1:length(hits[,1])){
                hits[i,5] <- as.character(iso[[1]][,1][as.character(iso[[1]][,2])==as.character(hits[i,1])])[1]
            }
        }
        ############################################################################
        removals <- data.frame(seq(1:length(rules)),rep(0,length(rules)),stringsAsFactors=FALSE)
        names(removals)=c("Rule","Counts")
        removals[1,2]=countrem1
        removals[2,2]=countrem2
        removals[3,2]=countrem3
        removals[4,2]=countrem4
        removals[5,2]=countrem5
        removals[6,2]=countrem6
        removals[7,2]=countrem7
        removals[8,2]=countrem8
        removals[9,2]=countrem9    
        removals[11,2]=countrem11
        removals <- removals[-c(10),] 
        if(rules[10]==TRUE){
            removals[7,1] <- "7/10"
            removals[9,1] <- "9/10"
        }
        ############################################################################
        # correct entries:
        for(i in 1:n_peaks){
            if(iso_from[i]!="0"){iso_from[i] <- substr(iso_from[i],3,nchar(iso_from[i]))}
            if(iso_to[i]!="0"){
                iso_to[i] <- substr(iso_to[i],3,nchar(iso_to[i]))
                this30 <- as.numeric(strsplit(iso_to[i],"/")[[1]])
                this31 <- as.character(ID[this30[1]])
                if(length(this30)>1){
                    for(j in 2:length(this30)){
                        this31 <- paste(this31,"/",ID[this30[j]],sep="")
                    }
                }
                iso_to[i] <- this31
            }
            if(iso_tolerance[i]!="0"){iso_tolerance[i] <- sub("0/","",iso_tolerance[i])}
            if(iso_charge[i]!="0"){iso_charge[i] <- sub("0/","",iso_charge[i])}
            if(iso_name[i]!="none"){
                this12 <- as.numeric(strsplit(iso_name[i],"/")[[1]][-1])
                if(length(this12)==1){
                    iso_name[i] <- isomat[this12,1]
                }else{
                    iso_name[i] <- isomat[this12[1],1]
                    for(j in 2:length(this12)){iso_name[i] <- paste(iso_name[i],isomat[this12[j],1],sep="/")}
                }
            }
            if(group1[i]!="0"){group1[i] <- substr(group1[i],3,(nchar(group1[i])))}
            if(group5[i]!="0"){group5[i] <- substr(group5[i],3,(nchar(group5[i])))}
            if(group2[i]!="0"){group2[i] <- substr(group2[i],3,(nchar(group2[i])))}
        }
        if(length(groupinfo)>0){
            for(i in 1:length(groupinfo)){
                if(groupinfo[i]!="minimum atom counts: no information"){groupinfo[i] <- sub("no information/","",groupinfo[i])}
            }
        }else{
            groupinfo <- "No groups detected"
        }
        ############################################################################
        groupcount <- data.frame(charge_isos,rep(0,length(charge_isos)),stringsAsFactors=FALSE)
        names(groupcount) <- c("Charge level","Counts")
        seq_charge_isos <- 1:length(charge_isos)
        if(length(group6)>0){
            for(i in 1:length(group6)){
                iso_to_split <- as.numeric(strsplit(as.character(group6[i]),"/")[[1]])
                iso_name_split <- c()
                for(j in 1:length(iso_to_split)){iso_name_split <- c(iso_name_split,seq_charge_isos[charge_isos==iso_to_split[j]])}
                groupcount[iso_name_split,2] <- c(groupcount[iso_name_split,2]+1)
            }
        }else{
            groupcount <- "No groups detected"
        }
        ############################################################################
        pattern_peak_list <- data.frame(peaklist,ID,group1,group2,iso_to,iso_name,iso_tolerance,iso_charge,stringsAsFactors=FALSE)
        pattern_peak_list <- pattern_peak_list[order(ID,decreasing=FALSE),]
        names(pattern_peak_list) <- c(names(peaklist),"peak ID","group ID","interaction level","to ID","isotope(s)","mass tolerance","charge level")
        #
        parameters <- data.frame(rttol[1],rttol[2],mztol,mzfrac,ppm,inttol,cutint,deter,stringsAsFactors=FALSE)
        #
        if(length(group4)>0){
            for(k in 1:length(group3)){
                group3[k] <- paste("/",group3[k],"/",sep="")
            }
            for(k in 1:length(group4)){
                this30 <- as.numeric(strsplit(group4[k],",")[[1]])
                this31 <- as.character(ID[this30[1]])
                if(length(this30)>1){
                    for(j in 2:length(this30)){
                        this31 <- paste(this31,",",ID[this30[j]],sep="")
                    }
                }
                group4[k] <- this31
            }
            grouping <- data.frame(group3,group4,group6,stringsAsFactors=FALSE)
            names(grouping) <- c("group ID","peak IDs","charge level")
        }else{
            grouping <- "no groups assembled"
        }
        #
        pattern <- list(pattern_peak_list,parameters,grouping,groupinfo,groupcount,removals,this11,deep,hits,elements,iso[[3]],rules)
        #
        names(pattern) <- c(
            "Patterns",
            "Parameters",
            "Peaks in pattern groups",
            "Atom counts",
            "Count of pattern groups",
            "Removals by rules",
            "Number of peaks with pattern group overlapping",
            "Number of peaks per within-group interaction levels",
            "Counts of isotopes",
            "Elements",
            "Charges",
            "Rule settings"
        )
        cat("done.\n\n")
        return(list(pattern, ST_2, ST_3))
        #return(pattern)
    }
