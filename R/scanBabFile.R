readCompressedData <- function(inputFile, path = ".", probeIDs = NULL)
{

    ## sort the probeIDs.  The bab file is sorted in increasing numerical order
    if(!is.null(probeIDs))
        probeIDs <- sort(unique(probeIDs));
  
    ## open connection to the binary file
    con <- file(paste(path, inputFile, sep = .Platform$file.sep), "rb");

    ## read the file header
    header <- readHeader(con);

	outputCount <- 0;

	## scan through once to find and count the beads we're interested in
    for(i in 1:header$nProbeIDs) {

		storeTmp <- readBin(con, integer(), size = 4, n = 2);
		probeID <- storeTmp[1];
		nbeads <- storeTmp[2];
        
        ## are we looking for this probeID?
        if( (is.null(probeIDs)) | (probeID %in% probeIDs) ) {
            outputCount <- outputCount + nbeads;
        }
        ## How many bytes can we now skip?
        inten <- as.logical(probeID) * 2^header$twoChannel * ( (2 * nbeads) + ( ( (nbeads - 1)  %/% 4) + 1) )
        coords <- (2 * header$nBytes * nbeads) - (0^(!header$useOffset) * 2 * nbeads);
        index <- 3^header$indexingMethod * nbeads;
        seek(con = con, where = sum(inten, coords, index), origin = "current");
    }

	## reset to the beginning of the file
	seek(con = con, where = 0, origin = "start");
	readHeader(con);

	## create output matrix and position counter
    output <- matrix(NA, ncol = 3 + (4^header$twoChannel), nrow = outputCount);
	pos <- 1;

    for(i in 1:header$nProbeIDs) {

		storeTmp <- readBin(con, integer(), size = 4, n = 2);
		probeID <- storeTmp[1];
		nbeads <- storeTmp[2];

        if(!is.null(probeIDs))
            probeIDs <- probeIDs[which(probeIDs >= probeID)];

		## if we've already got all the requested probes, quit the loop
		if(!is.null(probeIDs) & !length(probeIDs))
			break;
        
        ## are we looking for this probeID?
        if( (is.null(probeIDs)) | (probeID %in% probeIDs) ) {

			posEnd <- pos + nbeads - 1;

            output[pos:posEnd,1] <- rep(probeID, nbeads);     
            ## only read the intensities if they are there
            if(probeID) {  
                output[pos:posEnd,2] <- readIntensities(con, nbead = nbeads);
                if(header$twoChannel)
                    outputTmp[pos:posEnd,5] <- readIntensities(con, nbead = nbeads);
            } ## intensities not stored for probeID 0
            else {
                output[pos:posEnd,2] <- rep(0, nbeads);
                if(header$twoChannel) 
                    outputTmp[pos:posEnd,5] <- rep(0, nbeads);
            }
                
            coords <- readCoordinates(con = con, nbeads = nbeads, nBytes = header$nBytes, twoChannel = header$twoChannel, offset = header$useOffset, base2 = header$base2)

            output[pos:posEnd,3:4] <- coords[1:(2*nbeads)];
            if(header$twoChannel)
                output[pos:posEnd,6:7] <- coords[(2*nbeads+1):length(coords)];
                
            ## skip the locs file index, we don't need it here
            seek(con = con, where = nbeads * 3^header$indexingMethod, origin = "current");

			## increase the position counter
			pos <- pos + nbeads;
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

    if(!nrow(output)) {
		## tell the user if no probe IDs matched
        message("No matching probe IDs");
        output <- NULL;
    }
    return(output);
}
