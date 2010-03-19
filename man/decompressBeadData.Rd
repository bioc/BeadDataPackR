\name{decompressBeadData}
\alias{decompressBeadData}

\title{
    Decompress a file in the beadarray binary format
}
\description{
    Decompressed a file create by BeadDataPackR.  The original files that were compressed will be restored as accurately as possible, depending upon the degree of precision specified during the compression.
}
\usage{
decompressBeadData(inputFile, inputPath = ".", outputMask = NULL, outputPath = ".", outputNonDecoded = FALSE, roundValues = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{inputFile}{
    The name of the file to be read.
}
\item{inputPath}{
    Path where the compress file is located.  The default is to use the current working directory.
}
\item{outputMask}{
    Text specify the names of the output files.  The output files will have ".txt", "_Grn.locs" and (if approriate "_Red.locs") appended to this mask.  If left NULL the original names of the section will be used.
}
\item{outputPath}{
    Path to where the uncompressed version of the files should be written to.  The default is to use the current working directory.
}
\item{outputNonDecoded}{
    If TRUE the undecoded beads will be included in the output .txt file.  They will have ProbeID 0 and intensity 0, but the bead centre coordinates will be included.
}
  \item{roundValues}{
	The original Illumina text files give the bead centre coordinates to 7 significant figures.  When this argument is TRUE decompressed files are also truncated in this manner, whilst FALSE writes them to the full precision they are stored in the compressed file.
}
}

\value{
    Called primarily for its side effect, in which two (or three) files are written to the disk.  These files should be representative of the original files that were compressed.  The function returns, invisibly, the number of lines written in the .txt file.
}

\author{
    Mike L. Smith
}

\examples{ 
        dataPath <- system.file("data", package = "BeadDataPackR")
        decompressBeadData(inputFile = "example.bab", inputPath = dataPath, outputPath = ".")
}

\keyword{ IO }

