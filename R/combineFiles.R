combineFiles <- function(txt, locsGrn, locsRed = NULL, pb = NULL, verbose = FALSE) {
 
    locsRounded <- matrix(.Call("roundLocsFileValues", locsGrn, PACKAGE = "BeadDataPackR"), ncol = 2);
    txtFileCoords <- txt[,3:4];
   
    locs <- cbind(1:nrow(locsGrn), locsGrn, locsRed)

    ## add an index column to the txt file
    txtFileCoords <- cbind(1:nrow(txtFileCoords), txtFileCoords);
    
    locsKeyMult <- locsRounded[,1] * locsRounded[,2];
    txtKeyMult <- txtFileCoords[,2] * txtFileCoords[,3];
    
    if(!is.null(pb))
        setTxtProgressBar(pb, 0.20);
    
    ## if we find a duplicate entry then it may be due to the multiplication
    if(any(duplicated(locsKeyMult))) {
        locsKeyDiv <- locsRounded[,1] / locsRounded[,2];
        txtKeyDiv <- txtFileCoords[,2] / txtFileCoords[,3];
        idx <- which( (locsKeyDiv %in% txtKeyDiv) & (locsKeyMult %in% txtKeyMult) )
    }
    else {
        idx <- which( locsKeyMult %in% txtKeyMult );   
    }

    ## if there's still some duplicated use string concatonation. It's really slow!!
    if(length(idx) != nrow(txtFileCoords)) {
        if(verbose) message("Using string concatenation");
        locsKey <- paste(locsRounded[,1], locsRounded[,2]);
        txtKey <- paste(txtFileCoords[,2], txtFileCoords[,3]);
        idx <- which(locsKey %in% txtKey);

        ## with iScan we can still get perfect duplicate coordinates
        ## currently we remove them and make them non-decoded
        if(length(idx) > length(txtKey)) {
            duplicateList <- removeDuplicates(locsKey, txtKey);
            ## remove any from the txt file coordinates which have multiple matches in the locs file
            txtKey <- txtKey[ -duplicateList[[2]] ];
            txtFileCoords <- txtFileCoords[ -duplicateList[[2]], ];
            txt <- txt[ -duplicateList[[2]], ];
            idx <- which(locsKey %in% txtKey);
        }
    }
    
    if(!is.null(pb))
        setTxtProgressBar(pb, 0.30);
    
    ## find which are in both file and remove those that aren't from the locsCoords
    nonDeCoords <- locs[-idx,];
    locsFileCoords <- cbind( locs[,1], locsRounded, locs[,2:ncol(locs)] )[idx,];

    if(!is.null(pb))
        setTxtProgressBar(pb, 0.40);
    
    ## order the two files by the rounded x and y coords
    if(verbose) message("Reordering");
    txtFileCoords <- txtFileCoords[order(txtFileCoords[,2], txtFileCoords[,3]),]
    locsFileCoords <- locsFileCoords[order(locsFileCoords[,2], locsFileCoords[,3]),]

    txtFileCoords <- cbind(txtFileCoords, locsFileCoords[,c(1, 4:ncol(locsFileCoords) )]);

    ## combine 
    if(is.null(locsRed)) {
        txtFileCoords <- txtFileCoords[order(txtFileCoords[,1]), c(5,6,4)];
        txtFileCoords <- cbind(txt[,1:2], txtFileCoords);
        undecoded <- cbind(rep(0, nrow(nonDeCoords)), rep(0, nrow(nonDeCoords)), nonDeCoords[,c(2,3,1)])
        result <- rbind(undecoded, as.matrix(txtFileCoords));
        colnames(result) <- c("Code", "Grn", "GrnX", "GrnY", "LocsIdx") 
    }
    else {
        txtFileCoords <- txtFileCoords[order(txtFileCoords[,1]), c(5:8, 4)];
        txtFileCoords <- cbind(txt[,1:2], txtFileCoords[,1:2], txt[,5], txtFileCoords[,3:5])
        undecoded <- cbind(rep(0, nrow(nonDeCoords)), rep(0, nrow(nonDeCoords)), nonDeCoords[,2:3], rep(0, nrow(nonDeCoords)), nonDeCoords[,c(4:5, 1)])
        result <- rbind(undecoded, as.matrix(txtFileCoords));
        colnames(result) <- c("Code", "Grn", "GrnX", "GrnY", "Red", "RedX", "RedY", "LocsIdx")  
    }

    return(result);
}

removeDuplicates <- function(locsKey, txtKey) { 
    locsDup <- c(which(duplicated(locsKey)), which(duplicated(locsKey, fromLast = TRUE)))
    txtDup <- which(txtKey %in% locsKey[locsDup])
    return(list(locsDup, txtDup));
}

