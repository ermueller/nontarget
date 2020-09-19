# Rule 1: intensity ratio check over all isotopes at charge level
rule_1_check <- function(peaklist, iso_name, iso_from, iso_to, iso_tolerance, iso_charge, isomat, inttol){
    
    # Split the strings
    iso_to_split <- c(strsplit(iso_to,"/")[[1]][-1])
    iso_name_split <- c(strsplit(iso_name,"/")[[1]][-1])
    iso_tolerance_split <- c(strsplit(iso_tolerance,"/")[[1]][-1])
    iso_charge_split <- c(as.numeric(strsplit(iso_charge,"/")[[1]][-1]))
    
    # Which Isotopes do we keep?
    iso_to_keep  <- rep(TRUE,length(iso_to_split))
    
    # Go through each "to" split element
    for(j in 1:length(iso_to_split)){
        
        # Exclude peaks with more than one atom distance from this check
        if(isomat[as.numeric(iso_name_split[j]),5] == 1){
            
            # Calculate the minimum possible mass
            to_peak_expected_int_min      <- peaklist[as.numeric(iso_to_split[j]),2]*(1-inttol)
            current_peak_expected_int_max <- peaklist[i,2]*(1+inttol)
            
            possnumb <- (to_peak_expected_int_min/current_peak_expected_int_max) * max_dabund # Based on isotope with largest relative abundance
            
            if(possnumb < 1){
                possnumb <- 1
            }
            
            
            if((possnumb*3.0078) > (peaklist[i,1]*max(iso_charge_split))){ # largest charge z, minimum mass of hydrogen
                iso_to_keep[j] <- FALSE
                countrem1 <- countrem1 + 1
            } 
        }
    }
        
    # Exclude peaks with more than one atom distance from this check (Dunno if this line again is necessary?)
    iso_to_keep[isomat[as.numeric(that3),5]!=1] <- TRUE
    
    # If we have to change the isotope information
    if(!all(iso_to_keep)){
        
        # Reduce the count of found isotopes by one for the isotopes that have been remove by the rule
        isomat[as.numeric(iso_name_split[iso_to_keep==FALSE]),4] <- isomat[as.numeric(iso_name_split[iso_to_keep==FALSE]),4]-1
        
        # Go through the list of which ones to keep and remove the ones from iso_from which aren't being kept
        
        for(m in 1:length(iso_to_keep)){
            if(iso_to_keep[m]==FALSE){
                iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])
            }
        }
        
        # We only need to change these if
        
        iso_to_split <- iso_to_split[iso_to_keep]
        iso_name_split <- iso_name_split[iso_to_keep]
        iso_tolerance_split <- iso_tolerance_split[iso_to_keep]
        iso_charge_split <- iso_charge_split[iso_to_keep]
    }
    

    if(any(iso_to_keep)){
        iso_to[i] <- paste(c("0",iso_to_split),sep="/",collapse="/")
        iso_name[i] <- paste(c("0",iso_name_split),sep="/",collapse="/")
        iso_tolerance[i] <- paste(c("0",iso_tolerance_split),sep="/",collapse="/")
        iso_charge[i] <- paste(c("0",iso_tolerance_split),sep="/",collapse="/")
    }else{
        iso_to[i] <- "0"
        iso_name[i] <- "none" # (1) which isotope?
        iso_tolerance[i] <- "0"
        iso_charge[i] <- "0"
    }
} 


# (3.2) RULE2: intensity ratio check with LARGE mass tolerance ###########
rule_2_check <- function(){
    if(any(iso_to_keep)){ # anything left to check?
        if(rules[2]=="TRUE"){
            iso_to_split <- c(strsplit(iso_to[i],"/")[[1]][-1])
            iso_to_keep <- rep(TRUE,length(iso_to_split))
            that3 <- c(strsplit(iso_name[i],"/")[[1]][-1])
            that6 <- c(iso_to_split[duplicated(iso_to_split)==FALSE]) 
            that7 <- c(isomat[as.numeric(that3),5]==1)
            #iso_tolerance_split <- c(strsplit(iso_tolerance[i],"/")[[1]][-1])
            iso_charge_split <- c(as.numeric(strsplit(iso_charge[i],"/")[[1]][-1]))
            for(j in 1:length(that6)){
                if(any(that7[iso_to_split==that6[j]])){ # only on single-distanced peaks!
                    mpoldnew2 <- min(1/isomat[as.numeric(that3[iso_to_split==that6[j]&that7==TRUE]),3])
                    ratioC <- max(isomat[as.numeric(that3[iso_to_split==that6[j]&that7==TRUE]),6])
                    possnumb <- min(((samples[as.numeric(iso_to_split[iso_to_split==that6[j]]),2]*(1-inttol))/(samples[i,2]*(1+inttol)))*mpoldnew2)
                    if(possnumb<1){possnumb <- c(1)}
                    possmass <- min(isos[as.logical(match(isos[,1],isos[as.logical(match(isos[,2],isomat[as.numeric(that3[iso_to_split==that6[j]&that7==TRUE]),1],nomatch=FALSE)),1],nomatch=FALSE)),3])
                    if(ratioC!=0){
                        if((possnumb*possmass+((possnumb/ratioC)*12))>(samples[i,1]*max(isomat[as.numeric(that3[iso_to_split==that6[j]&that7==TRUE]),7]))){ # include ratios to C! BEWARE if ration=0->Inf->alway>mass
                            iso_to_keep[as.numeric(iso_to_split)==as.numeric(that6[j])] <- FALSE
                            countrem2 <- c(countrem2+1)
                        } # if ...
                    }else{
                        if((possnumb*possmass)>(samples[i,1]*max(isomat[as.numeric(that3[iso_to_split==that6[j]&that7==TRUE]),7]))){ # include ratio to C
                            iso_to_keep[iso_to_split==that6[j]] <- FALSE
                            countrem2 <- c(countrem2+1)
                        } # if ...
                    }
                } #
            } # for
            # remove entries in iso_from[i] (i.e."to") only  / NOT "from!"
            # reset: keep non-single distanced peaks!
            iso_to_keep[isomat[as.numeric(that3),5]!=1] <- TRUE
            for(m in 1:length(iso_to_keep)){
                if(iso_to_keep[m]==FALSE){
                    iso_from[as.numeric(iso_to_split[m])] <- sub(paste("/",i,sep=""),"",iso_from[as.numeric(iso_to_split[m])])
                }
            }
            if(any(iso_to_keep==FALSE)){
                isomat[as.numeric(that3[iso_to_keep==FALSE]),4] <- isomat[as.numeric(that3[iso_to_keep==FALSE]),4]-1
            }
            iso_to_split <- iso_to_split[iso_to_keep]
            that3 <- that3[iso_to_keep]
            iso_tolerance_split <- iso_tolerance_split[iso_to_keep]
            iso_charge_split <- iso_charge_split[iso_to_keep]
            if(any(iso_to_keep)){
                that4 <- c()
                that5 <- c()
                iso_tolerance_keep <- c()
                iso_charge_keep <- c()
                for(n in 1:length(iso_to_split)){
                    that4 <- paste(that4,iso_to_split[n],sep="/")
                    that5 <- paste(that5,that3[n],sep="/")
                    iso_tolerance_keep <- paste(iso_tolerance_keep,iso_tolerance_split[n],sep="/")
                    iso_charge_keep <- paste(iso_charge_keep,iso_charge_split[n],sep="/")
                }
                iso_to[i] <- paste("0",that4,sep="")
                iso_name[i] <- paste("0",that5,sep="")
                iso_tolerance[i] <- paste("0",iso_tolerance_keep,sep="")
                iso_charge[i] <- paste("0",iso_charge_keep,sep="")
            }else{
                iso_to[i] <- "0"
                iso_name[i] <- "none" # (1) which isotope?
                iso_tolerance[i] <- "0"
                iso_charge[i] <- "0"
            }
        } # if
    } # RULE2
}