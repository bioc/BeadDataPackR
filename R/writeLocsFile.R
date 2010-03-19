writeLocsFile <-
function(fileName, coords, nBeads)
{
    con <- file(fileName, "wb");
    writeBin(as.integer(c(1, 0, nBeads)), con = con, size = 4);
    writeBin(as.vector(coords), con = con, size = 4);
    close(con);
}

