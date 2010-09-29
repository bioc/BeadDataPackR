\name{extractLocsFile}
\alias{extractLocsFile}

\title{
    Retrieve only the .locs file information
}
\description{
    Provides a mechanism to extract the information from the original .locs file from a compressed .bab file, without the need to extract the intensity or probe ID values.
}
\usage{
extractLocsFile(inputFile, path = ".")
}

\arguments{
\item{inputFile}{
    The name of the .bab file to be read in.
}
\item{path}{
    Path to where the input file can be found.  Default is the current working directory.
}
}

\value{
    A matrix with two columns (four if two-channel data) containing the X and Y values of the bead centre coordinates supplied in the original .locs file.  For two-channel data the first two columns contain the coordinates from the green channel, with the red channel held in columns three and four.
}


\author{
    Mike L. Smith
}

\examples{
 
    dataPath <- system.file("data", package = "BeadDataPackR")   
    locs <- extractLocsFile(inputFile = "example.bab", path = dataPath)
    locs[1:10,]
}

\keyword{ IO }

