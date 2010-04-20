compressBeadData <-
function(txtFile, locsGrn, locsRed = NULL, outputFile = NULL, path = NULL, nBytes = 8, base2 = TRUE, fullLocsIndex = FALSE, nrow = NULL, ncol = NULL, progressBar = TRUE) {

    message(paste("\nCompressing", strsplit(txtFile, ".txt")));
    
    if(progressBar) {
        pb <- txtProgressBar(style=3)
        setTxtProgressBar(pb, 0.01)
    }
    else { 
        pb <- NULL
    }
    
    ## set the default name for the output file if one isn't specified
    if(is.null(outputFile))
        outputFile <- paste(strsplit(txtFile, ".txt"), "bab", sep = ".");
    
    ## identify the correct path to each of the file
    if(!is.null(path)) {
        txtFile <- file.path(path, txtFile);
        locsGrn <- file.path(path, locsGrn);
        outputFile <- file.path(path, outputFile);
        if(!is.null(locsRed))
        locsRed <- file.path(path, locsRed)
    }

    ## check we aren't using crazy numbers of bytes for storing the fractional parts
    if(nBytes > 8) { nBytes = 8 }
    else if(nBytes < 0) { nBytes = 0 }

    ## is this two channel?
    if(is.null(locsRed)) {
        twoChannel = FALSE
        ## makes sure we aren't using excessive bytes in the one channel case
        if(nBytes > 4) #message("For single channel data a maximum of 4 bytes can be specified");
            nBytes <- min(4, nBytes); 
    }
    else 
        twoChannel = TRUE
      
    ## read the data
    txt <- readBeadLevelTextFile(txtFile);
    if(progressBar) setTxtProgressBar(pb, 0.05)
    locsGrn <- readLocsFile(locsGrn);
    if(twoChannel)
      locsRed <- readLocsFile(locsRed);

    if(progressBar) setTxtProgressBar(pb, 0.1)
    
    ## combine the files and identify non-decoded beads
    combined <- combineFiles(txt, locsGrn, locsRed, pb = pb);
    
    if(progressBar) setTxtProgressBar(pb, 0.5)
    
    ## if we're using the fitted grid then do it
    if(!fullLocsIndex) {
        res <- createIndices(locsGrn, ncol, nrow, pb = pb);
        ## replace coordinates with shifted ones
        shifts <- res[[3]][seq(1,length(res[[3]]), 3)]
        if(any(as.logical(shifts))) {
            message("DEBUG: applying shifts");
            ## find which segments need to be shifted
            shiftIdx <- which(as.logical(shifts))
            for(i in shiftIdx) {
                ## find the beads in those segments
                segIdx <- which( (combined[,ncol(combined)] > (i*res[[5]][4] + 1)) & (combined[,ncol(combined)] < ((i+1)*res[[5]][4])) );
                ## shift them appropriately
                combined[segIdx,4] <- combined[segIdx,4] + shifts[i];
            }
        }
        indices <- (16 * res[[2]][,1]) + res[[2]][,2];
        if(progressBar) setTxtProgressBar(pb, 0.65)
    }
    
    ## if we're using a full index, reduce its size by one byte per bead
    if(fullLocsIndex) {
      combined <- cbind(combined[, 1:(ncol(combined)-1)], matrix(sapply(combined[,ncol(combined)], reduceIndexSize), ncol = 2, byrow = TRUE))
    }
    ## otherwise order the reduced index and combine them
    else {
      indices <- indices[combined[,ncol(combined)]]
      combined <- cbind(combined, indices);
    }

    ## determine whether we can use offset coords for the red channel
    useOffset <- FALSE
    if(twoChannel) 
      useOffset <- allowOffset(combined[,c(3,4,6,7)]) & (nBytes != 8);
    
    ## open the output file
    con <- file(outputFile, "wb")

    ## write the file header
    writeBabHeader(con = con, version = 1, combined = combined, nBytes = nBytes, twoChannel = twoChannel, useOffset = useOffset, base2 = base2, indexingMethod = fullLocsIndex, res = res);

    ## write the name of the array
    writeArrayName(txtFile, con = con);
  
    if(progressBar) setTxtProgressBar(pb, 0.7)
    
    ## write the body of the file
    writeBabBody(combined, con = con, twoChannel = twoChannel, nBytes = nBytes, useOffset = useOffset, base2 = base2, fullLocsIndex = fullLocsIndex, pb = pb);     
    close(con);
    
    if(progressBar) {
        setTxtProgressBar(pb, 1);
        close(pb);
    }
}

