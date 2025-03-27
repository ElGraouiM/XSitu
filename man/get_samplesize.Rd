\name{get_samplesize}

\alias{get_samplesize}

  
\title{Get the sample-size associated with a species' range}

\description{
Compute the sample-size (number of zones) associated with a species' range
}

\usage{
get_samplesize(x, fun=NULL, ...)
}

\arguments{
\item{x}{SpatRaster or SpatVector representing an (adjusted) range}
\item{fun}{function that computes the number of zones. The first argument is the total area of the zones in km^2. If \code{NULL}, the default function 
\code{function(A, omega) max(1, round(omega * sqrt(A/pi)))} is used}
\item{...}{additional arguments for \code{fun}}
}

\value{
A list with three items:
 \itemize{
    \code{range} thre range represented as a SpatVector of polygons
	\code{area} the range size in km2
	\code{n} the number of zones
 }
}

\examples{
rnge <- terra::rast(system.file("ex/range.tif", package="XSitu"))
land <- terra::readRDS(system.file("ex/land.rds", package="XSitu"))
sp <- terra::readRDS(system.file("ex/sp.rds", package="XSitu"))

x <- adjust_range(rnge, sp, land, include=50, exclude=250)

k <- get_samplesize(x, omega=1/40)
str(k)
}

\keyword{spatial}

