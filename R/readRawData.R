numberOfChannels <-
function(file, sep = "\t") {
    lines <- read.table(file, sep = sep, nrows = 2);
    if(ncol(lines) == 4) ##one channel
        return(1)
    else if(ncol(lines) == 7) ##two channel 
        return(2)
    else ##unexpected number of columns
        return(0)
}

readLocsFile <- function(fileName, verbose = FALSE) {
    if(verbose)
        message(paste("Reading:", fileName));
    con <- file(fileName, "rb")
    
    ##read the first two bytes
    readBin(con, integer(), n = 2, size = 4)
    ##3rd byte tells you how many probes there are
    nprobes <- readBin(con, integer(), n = 1, size = 4)
    
    coords <- matrix(ncol = 2, nrow = nprobes)
    colnames(coords) <- c("X", "Y")
    ##read in the whole file
    tmp <- readBin(con, double(), n = 2*nprobes, size = 4)
    ##store the x and y coords in the two columns
    coords[,1] <- tmp[seq(1, 2*nprobes, 2)]
    coords[,2] <- tmp[seq(2, 2*nprobes, 2)]

    close(con)
    return(coords)
}

readBeadLevelTextFile <-
function(file, sep = "\t", verbose = FALSE) {
    
    channels <- numberOfChannels(file, sep = sep);
    if(verbose)
        message(paste("Reading", file));
    if(channels == 1) 
        data <- matrix(unlist(scan(file, sep = "\t", what = list(integer(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 4)
    else if (channels == 2)
        data <- matrix(unlist(scan(file, sep = "\t", what = list(integer(), integer(), numeric(), numeric(), integer(), numeric(), numeric()), skip = 1, quiet = TRUE)), ncol = 7)
    else
        stop("Unknown input format!\nExpected 4 columns for single channel data or 7 columns for two channel data\n");
    
    return(data);
}
