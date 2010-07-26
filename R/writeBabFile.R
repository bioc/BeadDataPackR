writeBabHeader <- function(con, version, combined, nBytes, twoChannel, useOffset, base2, indexingMethod, res)
## write the file header
{
    writeBin(as.integer(version), con = con, size = 1);
    writeBin(as.integer(nrow(combined)), con = con, size = 4)
    writeBin(as.integer(length(unique(combined[,1]))), con = con, size = 4)
    writeBin(as.integer(c(nBytes, twoChannel, useOffset, base2)), con = con, size = 1);
    writeBin(as.integer(indexingMethod), con = con, size = 1);
    ## if we're using the fitted grid, record the number of segments and the segment information
    if(!indexingMethod) {
        ## number of segments
        writeBin(as.integer(res[[5]][1]), con = con, size = 1);
        writeBin(as.integer(res[[3]]), con = con, size = 2);
        writeBin(as.numeric(res[[4]]), con = con, size = 4);
    }

}

writeArrayName <- function(txtFile, con)
## given the text file to be read, extract the array name by removing ".txt"
## and any path at the front.  Then write the name to the binary file.
{  
    tmp <- strsplit(txtFile, "\\.")[[1]][1];
    tmp <- strsplit(tmp, .Platform$file.sep)[[1]]
    name <- tmp[length(tmp)];
    writeBin(name, con = con);

}

writeCoordinates <- function(coordinates, con, twoChannel, nBytes, useOffset = FALSE, base2 = FALSE, ensureSamePixel = TRUE) {
 
    ## formating to deal with cases where there is only one bead
    if(twoChannel) 
        coords <- matrix(coordinates, ncol = 4)
    else {
        coords <- matrix(coordinates, ncol = 2);
    }

    ## is this lossless?
    if( nBytes == (ncol(coords) * 2) ) {
        writeBin(as.vector(coords[,1:ncol(coords)]), con = con, size = 4);
    }
    ## if not, start compressing
    else if(nBytes) {

        ##get the fractional parts
        frac <- coords[,1:ncol(coords)] - floor(coords[,1:ncol(coords)]);
        ## determine the number of bits required for each coordinate
        nBits <- (2 * 2^(!twoChannel) * nBytes);
        
        ## set a multiplier based on base 2 or base 10
        if(base2)
            mult <- 2^nBits
        else
            mult <- 10^(max(which(2^nBits > 10^(1:5))));           
        
        tmpInts <- round(mult * frac);
        ## if we want to ensure the same pixel is used we fix it here
        if(any(tmpInts == 0) && ensureSamePixel) {
            idx <- which( (tmpInts == 0) & (frac != 0) );
            tmpInts[idx] <- 1;
        }
        ## deal with any cases that have been rounded to the maximal value
        ## we want to increment the integer part in this case
        if(any(tmpInts == mult)) {
            idx <- which(tmpInts == mult)
            coords[idx] = coords[idx] + 1;
            tmpInts[idx] = 0;
        }
        
        ## break those ints into blocks of 8 bits and then back to integers
        bits <- matrix(as.integer(matrix(sapply(tmpInts, FUN = intToBits), ncol = length(frac))[1:nBits,]), nrow = 8);
        ints <- .Call("bitsToInts", bits, PACKAGE = "BeadDataPackR");
          
        ## record the integer parts
        writeBin(as.integer(coords[,1:2]), con = con, size = 2); 
        if(twoChannel) {
            if(useOffset) 
                coords[,3:4] <- floor(coords[,3:4]) - floor(coords[,1:2]);
            writeBin(as.integer(coords[,3:4]), con = con, size = 2^(!useOffset) );
        }       
        
        ## now write the fractional data
        writeBin(as.integer(ints), con = con, size = 1);
    }
    else { ## if we aren't storing a fractional part, we should ceiling the values
        coords <- ceiling(coords);
        writeBin(as.integer(coords[,1:2]), con = con, size = 2);
        
        if(twoChannel) {
            if(useOffset) 
                coords[,3:4] <- coords[,3:4] - floor(coords[,1:2]);
            writeBin(as.integer(coords[,3:4]), con = con, size = 2^(!useOffset) );
        }
    }       
}

writeIntensities <- function(intensities, con) {

  ## write intensities as values between 0 and 65535
  ## if negative or have absolute value greater than 65535
  ## flag them as such

    ## identify those values needing flags and normalize the values to [0, 65535]
    neg <- (intensities < 0)
    intensities <- abs(intensities)
    large <- (intensities > 65535)
    intensities <- intensities %% 65536

    ## write the normalized intensities
    writeBin(as.integer(intensities), con = con, size = 2)
    ## write the flags
    flags <- .Call("composeIntensityFlags", as.integer(neg), as.integer(large), PACKAGE = "BeadDataPackR")
    writeBin(as.integer(flags), con = con, size = 1);
}

writeBabBody <- function(combined, con, twoChannel, nBytes, useOffset, base2, ensureSamePixel, fullLocsIndex, pb) {
           
    ## use the index to create a list, each elemet having probes of one type
    divided <- split(combined, combined[,1])
    
    for(i in 1:length(divided)) {
        
        if(i %/% 1000)
            setTxtProgressBar(pb, 0.7 + (0.3 * i/length(divided)))
        
        current <- matrix(divided[[i]], ncol = 5 + (4^twoChannel))

        writeBin(c(as.integer(current[1,1]), as.integer(nrow(current))), con = con, size = 4)
        if(current[1,1] != 0) 
            writeIntensities(current[,2], con = con)
            
        ## if this is two channel write those intensities
        if(twoChannel & (current[1,1] != 0) )
            writeIntensities(current[,5], con = con);
        
        ## now record the coordinates      
        if(twoChannel) {
            writeCoordinates(current[,c(3,4,6,7)], con = con, twoChannel = twoChannel, nBytes = nBytes, useOffset = useOffset, base2 = base2, ensureSamePixel = ensureSamePixel);
        }
        else {
            writeCoordinates(current[,3:4], con = con, twoChannel = twoChannel, nBytes = nBytes, useOffset = useOffset, base2 = base2, ensureSamePixel = ensureSamePixel);
        }
        
        ## record the index of the locs file
        if(fullLocsIndex) {
            writeBin(as.integer(current[, ncol(current) ]), con = con, size = 1)
            writeBin(as.integer(current[, ncol(current) ]), con = con, size = 2)
        }
        else {
            writeBin(as.integer(current[, ncol(current) ]), con = con, size = 1);
        }
          
    }
}
