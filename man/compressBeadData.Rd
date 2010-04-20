\name{compressBeadData}
\alias{compressBeadData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Write raw bead level data to a compressed format.
}
\description{
    Given raw bead level data, in the form of .txt and .locs file, this function combines the two producing a new file with the data stored in a compressed format.
}
\usage{
compressBeadData(txtFile, locsGrn, locsRed = NULL, outputFile = NULL,
                 path = NULL, nBytes = 8, base2 = TRUE, fullLocsIndex = FALSE,
                 nrow = NULL, ncol = NULL, progressBar = TRUE)
}

\arguments{
\item{txtFile}{
    The name of the .txt file to be read in.
}
\item{locsGrn}{
    The locs file for the green channel.
}
\item{locsRed}{
    The locs file for the red channel.  Only needed for two channel data.
}
\item{outputFile}{
    Name of the file to be created.
}
\item{path}{
    Path to where the input files can be found.  If NULL the current working directory is used.  This is also the directory where the output files will be written.
}
\item{nBytes}{
    Gives the number of bytes that are used to store the fractional parts of the bead coordinates.  For a single channel array the maximum value is 4, whilst it is 8 for a two channel array.  Any number larger than this is automatically set the the maximum value.  If the maximum value is used the coordinates are stored in the .bab file as single precision floating point numbers, as they are in the .locs files.  If a value smaller than the maximum is choosen then the integer parts of each coordinate are stored seperately.  The requested number of bytes are then used to store the fractional parts, with a corresponding loss of precision as the number of bytes decreases.
}
\item{base2}{
     If not using the full precision coordinates, the approximations can be stored as either a binary or decimal fraction.  Using a binary fraction (base2 = TRUE) provides a greater accuracy, but can lead to a meandering number of decimal places in the reconstructed .txt files.  If one wants a consistent number of decimal places, set base2 = FALSE.
}
\item{fullLocsIndex}{
    Default value of 0 uses a linear model fitted to each segment of the array to allow reconstruct the locs file when the file is decompressed.  Using a value of 1 a simple index is used to record the locs file order, but requires more space.
}
\item{ncol}{
    This specifies the number of columns in each grid segment on the array and, if left blank, can normally be infered from the grid coordinates.  However, this can fail for particularly small grids. If one wants or needs to specify them explicitly, these values can be found in the .sdf which accompanies the bead level output from the scanner.  The number of columns per segment can be found within the tag <SizeGridX>
}
\item{nrow}{
    See ncol.  If needed can be found within the <SizeGridY> tag in the .sdf file.
}
\item{progressBar}{
    By default the function uses a \code{\link[utils]{txtProgressBar}} to indicate progress through the compression.  Setting this argument to FALSE supresses the drawing of this progress bar.
}
}

\value{
    Primarily invoked for its side effect, which is to produce a compressed version of the input files. The function returns, invisibly, the total number of beads stored in the compressed file.
}

\details{
    In the future the file names will be determined automatically, rather than requiring manual entry of each.  The path argument may also be amended so there are seperate options for the locations of the input and output files.
}

\author{
    Mike L. Smith
}

\examples{
 
    dataPath <- system.file("data", package = "BeadDataPackR")   
    compressBeadData(txtFile = "example.txt", locsGrn = "example_Grn.locs", outputFile = "example.bab", path = dataPath, nBytes = 4, nrow = 326, ncol = 4)
}

\keyword{ IO }

