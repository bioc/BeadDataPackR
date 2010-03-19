combineFiles <-
function(txt, locsGrn, locsRed = NULL, pb = NULL) {
 
    locsRounded <- matrix(.Call("roundLocsFileValues", locsGrn), ncol = 2);
    txtFileCoords <- txt[,3:4];
   
    locs <- cbind(1:nrow(locsGrn), locsGrn, locsRed)

    ## add an index column to the txt file
    txtFileCoords <- cbind(1:nrow(txtFileCoords), txtFileCoords);
    
    locsKey <- locsRounded[,1] * locsRounded[,2];
    txtKey <- txtFileCoords[,2] * txtFileCoords[,3];

    if(!is.null(pb))
        setTxtProgressBar(pb, 0.20);
    
    ## if we find a duplicate entry then it may be due to the multiplication, so use string concatonation
    ## this could probably made better by checking the cause of the duplicates!
    if(any(duplicated(locsKey))) {
        dupsLocs <- c(which(duplicated(locsKey)), which(duplicated(locsKey, fromLast = TRUE)));
        dupsTxt <- c(which(duplicated(txtKey)), which(duplicated(txtKey, fromLast = TRUE)));
        locsKey[dupsLocs] <- (locsRounded[dupsLocs,1]) / (locsRounded[dupsLocs,2])
        txtKey[dupsTxt] <- (txtFileCoords[dupsTxt,2]) / (txtFileCoords[dupsTxt,3])
      #locsKey <- paste(locsRounded[,1], locsRounded[,2], sep = "_")
      #txtKey <- paste(txtFileCoords[,2], txtFileCoords[,3], sep = "_") 
    }

    if(any(duplicated(locsKey))) 
        message("Still some duplicates");
    
    if(!is.null(pb))
        setTxtProgressBar(pb, 0.30);
    
    ## find which are in both file and remove those that aren't from the locsCoords
    idx <- which(locsKey %in% txtKey);
    nonDeCoords <- locs[-idx,];
    locsFileCoords <- cbind( locs[,1], locsRounded, locs[,2:ncol(locs)] )[idx,];

    if(!is.null(pb))
        setTxtProgressBar(pb, 0.40);
    
    ## order the two files by the rounded x and y coords
    txtFileCoords <- txtFileCoords[order(txtFileCoords[,2], txtFileCoords[,3]),]
    locsFileCoords <- locsFileCoords[order(locsFileCoords[,2], locsFileCoords[,3]),]

    txtFileCoords <- cbind(txtFileCoords, locsFileCoords[,c(1, 4:ncol(locsFileCoords) )]);

    if(is.null(locsRed)) {
      txtFileCoords <- txtFileCoords[order(txtFileCoords[,1]), c(5,6,4)];
      txtFileCoords <- cbind(txt[,1:2], txtFileCoords);
      undecoded <- cbind(rep(0, nrow(nonDeCoords)), rep(0, nrow(nonDeCoords)), nonDeCoords[,c(2,3,1)])
    }
    else {
      txtFileCoords <- txtFileCoords[order(txtFileCoords[,1]), c(5:8, 4)];
      txtFileCoords <- cbind(txt[,1:2], txtFileCoords[,1:2], txt[,5], txtFileCoords[,3:5])
      undecoded <- cbind(rep(0, nrow(nonDeCoords)), rep(0, nrow(nonDeCoords)), nonDeCoords[,2:3], rep(0, nrow(nonDeCoords)), nonDeCoords[,c(4:5, 1)])
    }
      
    return(rbind(undecoded, as.matrix(txtFileCoords)));
}

