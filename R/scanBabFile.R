readCompressedData <- function(inputFile, path = ".", probeIDs = NULL)
{

    ## sort the probeIDs.  The bab file is sorted in increasing numerical order
    if(!is.null(probeIDs))
        probeIDs <- sort(unique(probeIDs));
  
    ## open connection to the binary file
    con <- file(paste(path, inputFile, sep = .Platform$file.sep), "rb");

    ## read the file header
    header <- readHeader(con);

    outputTmp <- list();

    for(i in 1:header$nProbeIDs) {

        probeID <- readBin(con, integer(), size = 4);
        nbeads <- readBin(con, integer(), size = 4);
        
        if(!is.null(probeIDs))
            probeIDs <- probeIDs[which(probeIDs >= probeID)];
        
        ## are we looking for this probeID?
        if( (is.null(probeIDs)) | (probeID %in% probeIDs) ) {
            outputTmp[[paste(probeID)]] <- matrix(ncol = 3+(4^header$twoChannel), nrow = nbeads);
            
            outputTmp[[paste(probeID)]][,1] <- rep(probeID, nbeads);     
            ## only read the intensities if they are there
            if(probeID) {  
                outputTmp[[paste(probeID)]][,2] <- readIntensities(con, nbead = nbeads);
                if(header$twoChannel)
                    outputTmp[[paste(probeID)]][,5] <- readIntensities(con, nbead = nbeads);
            } ## intensities not stored for probeID 0
            else {
                outputTmp[[paste(probeID)]][,2] <- rep(0, nbeads);
                if(header$twoChannel) 
                    outputTmp[[paste(probeID)]][,5] <- rep(0, nbeads);
            }
                
            coords <- readCoordinates(con = con, nbeads = nbeads, nBytes = header$nBytes, twoChannel = header$twoChannel, offset = header$useOffset, base2 = header$base2)

            outputTmp[[paste(probeID)]][,3:4] <- coords[1:(2*nbeads)];
            if(header$twoChannel)
                outputTmp[[paste(probeID)]][,6:7] <- coords[(2*nbeads+1):length(coords)];
                
            ## skip the locs file index, we don't need it here
            seek(con = con, where = nbeads * 3^header$indexingMethod, origin = "current");
        }
        else { ## We don't want this probe
            ## How many bytes can we skip?
            inten <- as.logical(probeID) * 2^header$twoChannel * ( (2 * nbeads) + ( ( (nbeads - 1)  %/% 4) + 1) )
            coords <- (2 * header$nBytes * nbeads) - (0^(!header$useOffset) * 2 * nbeads);
            index <- 3^header$indexingMethod * nbeads;
            seek(con = con, where = sum(inten, coords, index), origin = "current");
        }
    }
    ## close the bab file
    close(con);

    if(length(outputTmp)) {
        output <- matrix(ncol = ncol(outputTmp[[1]]), nrow = sum(unlist(lapply(outputTmp, nrow))));
        pos <- 1;
        for(i in 1:length(outputTmp)) {
            posEnd <- pos+nrow(outputTmp[[i]])-1;
            output[pos:posEnd,] <- outputTmp[[i]];
            pos <- posEnd+1;
        }
    }
    else { ## tell the user if no probe IDs matched
        message("No matching probe IDs");
        output <- NULL;
    }
    return(output);
}
