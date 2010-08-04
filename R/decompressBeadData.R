decompressBeadData <- function(input, inputPath = ".", outputMask = NULL, outputPath = ".", outputNonDecoded = FALSE, roundValues = TRUE, progressBar = TRUE)
{
    
    for(inputFile in input) {
        message(paste("\nDecompressing", inputFile));
        if(progressBar) {
            pb <- txtProgressBar(style=3)
            setTxtProgressBar(pb, 0.00)
        }
        else {
            pb <- NULL;
        }
        
        ## open connection to the binary file and inform the user
        con <- file(paste(inputPath, inputFile, sep = .Platform$file.sep), "rb");
        
        ## read the header
        header <- readHeader(con);
        ## parse the header to see if the user needs to be informed about outdated versions
        parseHeader(header);
        
        if(is.null(outputMask))
            outputMask <- header$arrayName
        
        ## create a matrix to hold the results
        if(!header$twoChannel) {
            txt <- matrix(ncol = 4, nrow = header$nBeads)
            colnames(txt) <- c("Code", "Grn", "GrnX", "GrnY")
            locs <- matrix(ncol = 3, nrow = header$nBeads);
        }
        else {
            txt <- matrix(ncol = 7, nrow = header$nBeads)
            colnames(txt) <- c("Code", "Grn", "GrnX", "GrnY", "Red", "RedX", "RedY")
            locs <- matrix(ncol = 5, nrow = header$nBeads);
        }
        
        if(progressBar) { setTxtProgressBar(pb, 0.02) }
        
        ## create a counter so we know where in the results matrix we stick the next results
        pos <- 1
        #message("Processing beads");
        for(i in 1:header$nProbeIDs) {
        
            ## update the progress bar
            if(progressBar) {
                if(i %/% 1000)
                    setTxtProgressBar(pb, 0.02 + (0.63 * i/header$nProbeIDs))
            }
            ## first 4 bytes are probeID, second are the number of beads of that type
            storeTmp <- readBin(con, integer(), size = 4, n = 2);
            probeID <- storeTmp[1];
            nbeads <- storeTmp[2];
            posEnd <- pos+nbeads-1
            
            ## fill in the probeIDs and intensities
            txt[pos:posEnd,1] <- rep(probeID, nbeads)
            if( probeID == 0 ) {
                txt[pos:posEnd,2] <- rep(0, nbeads)
                if(header$twoChannel)
                    txt[pos:posEnd,5] <- rep(0, nbeads)
            }
            else {
                txt[pos:posEnd,2] <- readIntensities(con, nbead = nbeads);
                if(header$twoChannel)
                    txt[pos:posEnd,5] <- readIntensities(con, nbead = nbeads);
            }
                
            coords <- readCoordinates(con = con, nbeads = nbeads, nBytes = header$nBytes, twoChannel = header$twoChannel, offset = header$useOffset, base2 = header$base2)    

            locs[pos:posEnd, 2:3] <- coords[1:(2*nbeads)];
            if(header$twoChannel) {
                locs[pos:posEnd, 4:5] <- coords[(2*nbeads+1):length(coords)];
            }      

            if(header$indexingMethod) {
                locs[pos:posEnd, 1] <- readBin(con, integer(), size = 1, n = nbeads, signed = FALSE) * 65536;
                locs[pos:posEnd, 1] <- locs[pos:posEnd, 1] + readBin(con, integer(), size = 2, n = nbeads, signed = FALSE)
            }
            else {
                locs[pos:posEnd, 1] <- readBin(con, integer(), size = 1, n = nbeads, signed = FALSE)
            }

            ## update the current position
            pos <- pos + nbeads
        }
        ## close the connection
        close(con)

        ## if the red channel are just offsets from the green then correct this
        if(header$useOffset)
            locs[,4:5] <- floor(locs[,2:3]) + locs[,4:5];

        ## insert the coordinates into the txt file
        if(roundValues) {
            txt[,3:4] <- .Call("roundLocsFileValues", as.numeric(locs[, 2:3]));
            if(header$twoChannel)
                txt[,6:7] <- .Call("roundLocsFileValues", as.numeric(locs[, 4:5])); 
        }
        else {
            txt[,3:4] <- locs[, 2:3];
            if(header$twoChannel)
                txt[,6:7] <- locs[, 4:5]; 
        }
        
        if(!header$indexingMethod) {
            decoded <- decodeIndices(locs[,1], locs[,2:3], header$nSegs, header$marks, header$coeffs, pb = pb);
            locs[,2:3] <- reformCoordinates(locs[,2:3], header$nSegs, header$marks);
            txt[,3:4] <- reformCoordinates(txt[,3:4], header$nSegs, header$marks);
            locs <- locs[decoded,2:(ncol(locs))]
        }
        else {
            locs <- locs[order(locs[,1]), 2:(ncol(locs))];
        }
 
        ## remove the nondecoded beads if desired
        if( (!outputNonDecoded) & (length(which(txt[,1] == 0))) ) 
            txt <- txt[-which(txt[,1] == 0),]
        
        ## write the output files
        write.table(txt, file = paste(outputPath, paste(outputMask, ".txt", sep = ""), sep = .Platform$file.sep), sep = "\t", quote = FALSE, row.names = FALSE)
        
        if(progressBar) { setTxtProgressBar(pb, 0.90) };
        
        writeLocsFile(file = paste(outputPath, paste(outputMask, "_Grn.locs", sep = ""), sep = .Platform$file.sep), t(locs[,1:2]), nBeads = header$nBeads);
        if(header$twoChannel) {
            writeLocsFile(file = paste(outputPath, paste(outputMask, "_Red.locs", sep =""), sep = .Platform$file.sep), t(locs[,3:4]), nBeads = header$nBeads);
        }
        
        if(progressBar) {
            setTxtProgressBar(pb, 1);
            close(pb);
        }
    }
}

