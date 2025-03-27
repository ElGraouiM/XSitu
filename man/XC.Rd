\name{XC}

\alias{XC}

  
\title{Compute the XC score}

\description{
Compute the XC score
}

\usage{
XC(regions, sample, env=NULL, envfun=NULL, minssize=10, maxdist=1500, return_network=FALSE)
}

\arguments{
\item{regions}{SpatVector or polygons}
\item{sample}{SpatVector or points}
\item{env}{SpatRaster with environmental data or NULL}
\item{envfun}{function}
\item{minssize}{numeric}
\item{maxdist}{numeric}
\item{return_network}{logical}
}

\value{
list with a numeric value \code{XC}, the ex-situ conservation score, and a data.frame with the a row for all adjacent zones
 from to adj       xf       yf       xt       yt   w      dst

 \itemize{
    \code{from} from zone ID 
	\code{to}  tp zone ID
	\code{adj} Are the zones adjacent?
	\code{xf} from x coordinate
	\code{yf} from y coordinate
	\code{xt} to x coordinate
	\code{yt} to y coordinate
	\code{w} relative vertex weight. 0 (both zones have a sample), 0.5 (one zone has a sample) and 1 (neither zone has a sample)
	\code{dst} distance  between vertices
 }
}

\examples{
rnge <- terra::rast(system.file("ex/range.tif", package="XSitu"))
land <- terra::readRDS(system.file("ex/land.rds", package="XSitu"))
sp <- terra::readRDS(system.file("ex/sp.rds", package="XSitu"))
mprc <- readRDS(system.file("ex/m_prc.rds", package="XSitu"))
mtmp <- readRDS(system.file("ex/m_tmp.rds", package="XSitu"))
env <- terra::rast(system.file("ex/env.tif", package="XSitu"))

x <- adjust_range(rnge, sp, land, include=50, exclude=250)
k <- get_samplesize(x, omega=1/40)
xy <- c(terra::init(x, "x"), terra::init(x, "y"))

set.seed(321)
zones <- make_zones(xy, k$range, k$n, spread=TRUE)

envdist <- function(x) {
	x$tmp[x$tmp > 13] <- 13 
	x$prc[x$prc > 2000] <- 2000
	p_tmp <- predict(mtmp, x)
	p_pr <- predict(mprc, x)
	rowMeans(cbind(p_tmp, p_pr))
}

x <- XC(zones, sp[sp$sample_type=="seed"], env=NULL, envfun=NULL, minssize=10, maxdist=1500, return_network=FALSE)
x$XC 
head(x$dist)

y <- XC(zones, sp[sp$sample_type=="seed"], env=env, envfun=envdist, minssize=10, maxdist=1500, return_network=FALSE)
y$XC 
head(y$dist)
}

\keyword{spatial}

