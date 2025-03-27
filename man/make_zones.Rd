\name{make_zones}

\alias{make_zones}

  
\title{Divide a species' range into zones}

\description{
Divide a species' range into \code{n} zones of approximately equal area
}

\usage{
make_zones(x, range, n, spread=TRUE)
}

\arguments{
\item{x}{SpatRaster with two layers and values represeting longitude (or X) and latitude (or Y) coordinates}
\item{range}{SpatVector with the species' adjusted range}
\item{n}{positive integer. The number of zones to create}
\item{spread}{logical. If \code{TRUE}, seeds are used for the k-means allocated by the size of the subpolygons}
}

\value{
SpatVector 
}

\examples{
rnge <- terra::rast(system.file("ex/range.tif", package="XSitu"))
land <- terra::readRDS(system.file("ex/land.rds", package="XSitu"))
sp <- terra::readRDS(system.file("ex/sp.rds", package="XSitu"))

x <- adjust_range(rnge, sp, land, include=50, exclude=250)
k <- get_samplesize(x, omega=1/40)
xy <- c(terra::init(x, "x"), terra::init(x, "y"))

set.seed(321)
zones <- make_zones(xy, k$range, k$n, spread=TRUE)
terra::plot(zones, col=rainbow(25))
}

\keyword{spatial}

