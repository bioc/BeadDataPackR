readHeader <- function(con) {

  ## read the header information and return as a list
  header <- list();
  header$version <- readBin(con, integer(), size = 1);
  header$nBeads <- readBin(con, integer(), size = 4);
  header$nProbeIDs <- readBin(con, integer(), size = 4);
  header$nBytes <- readBin(con, integer(), size = 1);
  header$twoChannel <- readBin(con, integer(), size = 1);
  header$useOffset <- readBin(con, integer(), size = 1);
  header$base2 <- readBin(con, integer(), size = 1);
  header$indexingMethod <- readBin(con, integer(), size = 1);
  if(!header$indexingMethod) {
    header$nSegs <- readBin(con, integer(), size = 1);
    header$marks <- readBin(con, integer(), size = 2, n = header$nSegs * 3);
    header$coeffs <- readBin(con, numeric(), size = 4, n = header$nSegs * 6);
  }
  header$arrayName <- readBin(con, character());
  return(header);
}

applyFlags <- function(inten) {

    ## apply the flags to each intensity
  
    inten[1] <- inten[1] + (inten[2] * 65536)
    if(inten[3])
        inten[1] <- -inten[1]

    return(inten[1])
}

readIntensities <- function(con, nbeads)
{
    inten <- readBin(con, integer(), size = 2, signed = FALSE, n = nbeads)
    flags <- readBin(con, integer(), size = 1, signed = TRUE, n = ( ( (nbeads-1) %/% 4 ) + 1 ) )

    ## get the bits we're interested in
    flags2 <- .Call("int2Bits", as.integer(flags), PACKAGE = "BeadDataPackR");
    flags3 <- matrix(c(inten, flags2[1:length(inten),]), ncol = 3)
    
    inten <- apply(flags3, 1, applyFlags)  
    return(inten)
}

readCoordinates <- function(con, nbeads, nBytes, twoChannel, offset = FALSE, base2 = FALSE)
{

    if(nBytes == (4 * (2^twoChannel)) ) {
        coords <- readBin(con, double(), size = 4, n = 2 * (2^twoChannel) * nbeads)
    }
    else {

        
        ## read the integer parts.  Grn first, then Red, which may be using offsets
        coords <- readBin(con, integer(), size = 2, n = 2*nbeads, signed = TRUE);
        if(twoChannel) {
            coords <- c(coords, readBin(con, integer(), size = 2^(!offset), n = 2*nbeads, signed = TRUE));
        }
        
        if(nBytes) {
            
            ## calculate the number of bits required for each coordinate
            nBits <- 2 * 2^(!twoChannel) * nBytes;
            
            ## set a multiplier based on base 2 or base 10
            if(base2)
                div <- 2^nBits - 1
            else
                div <- 10^(min(which(2^nBits < 10^(1:5)))-1);
            
            ## read the stored fractional parts
            tmp <- readBin(con, integer(), size = 1, n = nBytes * nbeads, signed = FALSE);
            ## convert them to bits
            bits <- matrix(as.integer(sapply(tmp, intToBits)[1:8,]), nrow = nBits);
            ## convert back into 2/4 integers
            frac <- apply(bits, 2, reformInt) / div;
            coords <- coords + frac;
        }
    }

  return(coords);
}

