\name{readCompressedData}
\alias{readCompressedData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Extract data for specific bead-types from a compressed file
}
\description{
    Given a list of probeIDs this function can scan a compressed .bab file for matching entries and return the data as a data.frame within R, rather than decompressing the data and generating new files.
}
\usage{
readCompressedData(inputFile, path = ".", probeIDs = NULL)
}

\arguments{
\item{inputFile}{
    The name of the .bab file to be read in.
}
\item{path}{
    Path to where the input file can be found.  Default is the current working directory.
}
\item{probeIDs}{
    List the probe IDs for which data should be obtained.  If left NULL then every probe on the array is returned.
}
}

\value{
    If the requested probe IDs are present the function returns a data.frame with one row per bead.  If the probes are not found in the file then the function returns NULL and informs the user.
}


\author{
    Mike L. Smith
}

\examples{
 
    dataPath <- system.file("data", package = "BeadDataPackR")   
    readCompressedData(inputFile = "example.bab", path = dataPath, probeIDs = c(10008, 10010))
}

\keyword{ IO }

