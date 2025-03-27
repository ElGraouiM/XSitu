\name{adjust_range}

\alias{adjust_range}

  
\title{Adjust the range of a species}

\description{
Adjust the range of a species that was predicted by a species distribution model (SDM). The input is a SpatRaster with values \code{1} where the environment is deemed suitable by the SDM and \code{NA} in other areas. An inclusion buffer is used to include areas near locations with known occurrence, and an exclusion buffer
is used to exclude areas that are further away from these locations. 
}

\usage{
adjust_range(x, sp, land, include=10, exclude=100)
}

\arguments{
\item{x}{SpatRaster with a presence/absence range map of a species}
\item{sp}{SpatVector with known occurrence locations for the species}
\item{land}{SpatVector or SpatRaster identifying areas to always exclude (e.g., not on land)}
\item{include}{positive number. The inclusion distance (km). Areas that are as close to known occurrence locations are included in the range map}
\item{exclude}{positive number. The exclusion distance (km). Areas that are further than this distance to known occurrence locations are excluded from the the range map}
}

\value{
SpatRaster 
}

\examples{
rnge <- terra::rast(system.file("ex/range.tif", package="XSitu"))
land <- terra::readRDS(system.file("ex/land.rds", package="XSitu"))
sp <- terra::readRDS(system.file("ex/sp.rds", package="XSitu"))

x <- adjust_range(rnge, sp, land, include=50, exclude=250)

terra::plot(x)

# more illustration
xr <- sum(c(rnge, x+1), na.rm=TRUE)
bb <- rbind(terra::buffer(sp, 50000), terra::buffer(sp, 250000))
bb$size <- rep(c("include", "exclude"), each=nrow(sp))
bb <- terra::as.lines(terra::aggregate(bb, "size"))
sp$label <- "occurrence"

terra::plot(xr, levels=c("Removed", "Added", "Kept"), type="classes", plg=list(title="Range"), mar=c(1,1,1,6))
terra::lines(land, lwd=4, col="light gray")
terra::plot(bb, "size", col=c("green", "red"), lwd=2, add=TRUE, sort=FALSE, type="classes", 
				plg=list(x=43.1, y=2.46, title="Buffer", lwd=2))
terra::plot(sp, "label", col="blue", add=TRUE, plg=list(x=43.1, y=-2))
}


\keyword{spatial}

