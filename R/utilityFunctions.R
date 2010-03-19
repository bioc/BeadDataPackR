
bitsToInt <- function(bits) {  
    return(sum(2^(0:7) * bits))
}

convertTo2Ints <- function(num) {

  res <- NULL;
  res[2] <- num %% 16;
  res[1] <- (num - res[2]) / 16;
  return(res);
}

reformInt <- function(bits) { 
    ## converts a sequence of bits into an integer
    n <- length(bits) - 1;
    return(sum(2^(0:n) * bits));
}


allowOffset <- function(data) {

    ## determine whether the red and green coords are close enough
    ## to store a shift, rather than the red coords themselves.
    x <- max(abs(data[,1] - data[,3])) < 128;
    y <- max(abs(data[,2] - data[,4])) < 128;
  
    return(x & y);
}

reformCoordinates <- function(coords, nSegs, marks) {

  ## take the coordinates that have been shifted and alter them back to how they were originally

  for(i in 1:nSegs) {

    if(marks[(3 * i) - 2]) {
      idx <- which( (coords[,2] >= marks[(3 * i) - 1]) & (coords[,2] <= marks[(3 * i)] ) )
      coords[idx,2] <- coords[idx,2] - marks[(3 * i) - 2];
    }
  
  }
  return(coords);
}

