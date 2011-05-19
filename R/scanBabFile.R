readCompressedData <- function(inputFile, path = ".", probeIDs = NULL)
{

    ## sort the probeIDs.  The bab file is sorted in increasing numerical order
    if(!is.null(probeIDs))
        probeIDs <- sort(unique(probeIDs));
  
    ## open connection to the binary file
    con <- file(paste(path, inputFile, sep = .Platform$file.sep), "rb");

    ## read the file header
    header <- readHeader(con);

    ## scan through once to find and count the beads we're interested in
    ## only need to do this if probeIDs != NULL
    if(!is.null(probeIDs)) {
        outputCount <- 0;
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
    }
    else {
        outputCount <- header$nBeads;
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
                    output[pos:posEnd,5] <- readIntensities(con, nbead = nbeads);
            } ## intensities not stored for probeID 0
            else {
                output[pos:posEnd,2] <- rep(0, nbeads);
                if(header$twoChannel) 
                    output[pos:posEnd,5] <- rep(0, nbeads);
            }
                
            coords <- readCoordinates(con = con, nbeads = nbeads, nBytes = header$nBytes, twoChannel = header$twoChannel, offset = header$useOffset, base2 = header$base2)

            output[pos:posEnd,3:4] <- coords[1:(2*nbeads)];
            if(header$twoChannel)
                output[pos:posEnd,6:7] <- coords[(2*nbeads+1):length(coords)];
                
            ## skip the locs file index, we don't need it here
            seek(con = con, where = nbeads * 3^header$indexingMethod, origin = "current");
            ## we do want it after all!
            #if(header$indexingMethod) {
            #    output[pos:posEnd, ncol(output)] <- readBin(con, integer(), size = 1, n = nbeads, signed = FALSE) * 65536;
            #    output[pos:posEnd, ncol(output)] <- output[pos:posEnd, ncol(output)] + readBin(con, integer(), size = 2, n = nbeads, signed = FALSE)
            #}
            #else {
            #    output[pos:posEnd, ncol(output)] <- readBin(con, integer(), size = 1, n = nbeads, signed = FALSE)
            #}        

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
    
    ## round the coordinates to match the text file
    output[,3:4] <- matrix(.Call("roundLocsFileValues", output[,3:4], PACKAGE = "BeadDataPackR"), ncol = 2);
    if(header$twoChannel)
        output[,6:7] <- matrix(.Call("roundLocsFileValues", output[,6:7], PACKAGE = "BeadDataPackR"), ncol = 2);
    

    if(!nrow(output)) {
	## tell the user if no probe IDs matched
        message("No matching probe IDs");
        output <- NULL;
    }
    else {
        if(header$twoChannel) {
            colnames(output) <- c('ProbeID','Grn','GrnX','GrnY','Red','RedX','RedY')
        } 
        else {
            colnames(output) <- c('ProbeID','Grn','GrnX','GrnY')
        }
    }
    return(output);
}

extractLocsFile <- function(inputFile, path = ".") {
##
## Function will extract only the .locs file from a .bab and return it as a matrix
##
    
    ## open connection to the binary file
    con <- file(paste(path, inputFile, sep = .Platform$file.sep), "rb");

    ## read the file header
    header <- readHeader(con);

    ## create output matrix and position counter
    locs <- matrix(NA, ncol = 2 + (3^header$twoChannel), nrow = header$nBeads);
    pos <- 1
    
    for(i in 1:header$nProbeIDs) {

        ## skip everything other than the locs information
        storeTmp <- readBin(con, integer(), size = 4, n = 2);
        probeID <- storeTmp[1];
        nbeads <- storeTmp[2];
        posEnd <- pos + nbeads - 1;
        
        ## How many bytes can we now skip?
        inten <- as.logical(probeID) * 2^header$twoChannel * ( (2 * nbeads) + ( ( (nbeads - 1)  %/% 4) + 1) );
        seek(con = con, where = inten, origin = "current");
        
        coords <- readCoordinates(con = con, nbeads = nbeads, nBytes = header$nBytes, twoChannel = header$twoChannel, offset = header$useOffset, base2 = header$base2)
        
        locs[pos:posEnd,2:3] <- coords[1:(2*nbeads)];
        if(header$twoChannel)
            locs[pos:posEnd,4:5] <- coords[(2*nbeads+1):length(coords)];
        
        ## read the locs file index
        if(header$indexingMethod) {
            locs[pos:posEnd, 1] <- readBin(con, integer(), size = 1, n = nbeads, signed = FALSE) * 65536;
            locs[pos:posEnd, 1] <- locs[pos:posEnd, 1] + readBin(con, integer(), size = 2, n = nbeads, signed = FALSE)
        }
        else {
            locs[pos:posEnd, 1] <- readBin(con, integer(), size = 1, n = nbeads, signed = FALSE)
        }

        ## increase the position counter
        pos <- pos + nbeads;        

    }
    ## close the bab file
    close(con);
    
    ## if the red channel are just offsets from the green then correct this
    if(header$useOffset)
        locs[,4:5] <- floor(locs[,2:3]) + locs[,4:5];
    
    ## reorder into locs file order
    if(!header$indexingMethod) {
        decoded <- decodeIndices(locs[,1], locs[,2:3], header$nSegs, header$marks, header$coeffs, pb = NULL);
        locs[,2:3] <- reformCoordinates(locs[,2:3], header$nSegs, header$marks);
        locs <- locs[decoded,2:(ncol(locs))]
    }
    else {
        locs <- locs[order(locs[,1]), 2:(ncol(locs))];
    }

    return(locs);   
}
