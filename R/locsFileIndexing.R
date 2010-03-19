createIndices <- function(locs, nrow = NULL, ncol = NULL, pb = NULL) 
{
    
    ## used during compression to create an index between 0 & 14 in both
    ## the x&y directions for each bead.  This combined with the predicted
    ## coordinates should allow the grid reconstruction.
    
    ## if we haven't been give the grid size, try to calculate it
    if(is.null(nrow) || is.null(ncol)) {
        segsize <- which((locs[-1,1]-locs[-length(locs[,1]),1])>100)[1]
        seg <- locs[1:segsize,]
        tempy <- which(seg[-length(seg[,2]),2]>seg[-1,2])
        nrow <- as.numeric(names(sort(table(tempy[-1]-tempy[-length(tempy)]),decreasing=T))[1]) 
        ncol <- segsize/nrow
    }

    nSegs <- nrow(locs) / (nrow * ncol);

    res <- list()

    ## store the max ycoord of the previous segment
    prevMaxY <- 0
  
    ## process each segment at a time
    for(i in 1:nSegs ) {
    
        if(!is.null(pb))
            setTxtProgressBar(pb, 0.5 + (0.1 * i / nSegs));
      
        seg <- locs[( (nrow * ncol * (i-1) ) + 1):(nrow * ncol * i),  ]
        grid <- cbind(rep(0:(ncol-1), each = nrow), rep(0:(nrow-1), ncol))

        ## adjust for overlapping segments
        ## store the adjust coords and grid indices
        if( floor(min(seg[,2])) < prevMaxY ) {
            message("Overlapping segment found! Adjusting...");
            shift <- prevMaxY - floor(min(seg[,2])) + 1
            seg[,2] <- seg[,2] + shift
            res[[4*i-1]] <- c(shift, floor(min(seg[,2])), ceiling(max(seg[,2])));
        }
        else {
            res[[4*i-1]] <- c(0, floor(min(seg[,2])), ceiling(max(seg[,2])));
        }
        prevMaxY <- ceiling(max(seg[,2]))
        res[[(4*i)-3]] <- seg
        res[[(4*i)-2]] <- grid %% 15

        ## fit the lm to the adjusted coords and store the coefficients
        lmX <- lm(grid[,1] ~ 1 + seg[,1] + seg[,2])
        lmY <- lm(grid[,2] ~ 1 + seg[,1] + seg[,2])
        res[[4*i]] <- c(lmX$coefficients, lmY$coefficients)
    }

    resList <- list();
    resList[[1]] <- res[[1]];
    resList[[2]] <- res[[2]];
    resList[[3]] <- res[[3]];
    resList[[4]] <- res[[4]];
    if(nSegs > 1) {
        for(i in seq(5, length(res), 4)) {
            resList[[1]] <- rbind(resList[[1]], res[[i]]);
            resList[[2]] <- rbind(resList[[2]], res[[i+1]]);
            resList[[3]] <- c(resList[[3]], res[[i+2]]);
            resList[[4]] <- c(resList[[4]], res[[i+3]]);
        }
    }
    resList[[5]] <- c(nSegs, nrow, ncol, nrow * ncol);
    return(resList);
}



decodeIndices <- function(indices, locs, nSegs, marks, coefficients, pb = NULL) {

    ## used during decompression
    ## decode the 0-14 indicies into one for the whole locs file
    ## indices are 8bit numbers and need to be transformed into 2 4bit numbers
    ## input locs are shifted coords, not the originals

    ## convert the indices to 2 numbers
    indices <- .Call("decodeInd", as.integer(indices), PACKAGE = "BeadDataPackR");

    res <- NULL
    ## loop over each segment
    for(i in 1:nSegs) {
    
        if(!is.null(pb))
            setTxtProgressBar(pb, 0.75 + (0.15 * i / nSegs));
        
        ## extract the beads in this segment
        idx <- which( (locs[,2] >= marks[(3 * i) - 1]) & (locs[,2] <= marks[(3 * i)] ) )
        seg <- cbind(locs[idx,], indices[idx,])

        ## get the model coefficiets for this segment
        coeff <- coefficients[(6*i-5):(6*i)]
        
        ## obtain the predicted coordinates 
        predX <- round(coeff[1] + seg[,1] * coeff[2] + seg[,2] * coeff[3])
        predY <- round(coeff[4] + seg[,1] * coeff[5] + seg[,2] * coeff[6])
        modX <- predX %% 15
        modY <- predY %% 15

        predX2 <- .Call("adjustValues", matrix(as.integer(cbind(predX, modX, seg[,3])), ncol = 3), PACKAGE = "BeadDataPackR");
        predY2 <- .Call("adjustValues", matrix(as.integer(cbind(predY, modY, seg[,4])), ncol = 3), PACKAGE = "BeadDataPackR");
        trueIdx <- .Call("returnTrueIndex", as.integer(predX2), as.integer(predY2), as.integer(max(predY2)+1), PACKAGE = "BeadDataPackR");
        
        res <- c(res, idx[order(trueIdx)] );
    }
    return(res);
}

reduceIndexSize <- function(index) {

    ## transform a 32bit integer into 8bit and 16bit pieces that can be multiplied
    ## this reduces the index size by one byte per bead
  
    return(c((index-1) %/% 65536, (index-1) %% 65536));
}


