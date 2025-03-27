

adjust_range <- function(x, sp, land, include=10, exclude=100) { 

	stopifnot(include <= exclude)

	# km to m
	CAmin <- include * 1000
	CAmax <- exclude * 1000

	if ((CAmax > 0) && (CAmax < Inf)) {
		ca_remove <- terra::buffer(sp, CAmax) 
		x <- terra::mask(x, ca_remove, updatevalue=NA)
	}
	if (CAmin > 0) {
		ca_add <- terra::buffer(sp, CAmin, quadsegs=12) 
		x <- terra::rasterize(ca_add, x, update=TRUE)
	}
	if (!is.null(land)) {
		x <- terra::mask(x, land)
	} 
	# terra::ifel(x==1, 1, NA)
	x
}



EnvDist <- function(env, regions, envfun) {
	e <- terra::extract(env, regions, fun=mean, na.rm=TRUE, ID=FALSE)
	ed1 <- stats::dist(e[,1])
	ed2 <- stats::dist(e[,2])	
	ed <- data.frame(ed1, ed2)
	names(ed) <- names(env)
	# express environmental distance expressed as geographic distance
	envd <- envfun(ed)
	structure(envd, class = 'dist', Size=attr(ed1, "Size"))
}


.n_zones <- function(x, min_area, m=3) {
	round(pmax(1, pmin(x, m*sqrt(x))))
}


make_zones <- function(x, range, n, spread=TRUE) {

	if (spread) { # spread sample across polygons
		drange <- terra::disagg(range)
		rr <- terra::rasterize(drange, terra::rast(x), 1:nrow(drange))
		p <- terra::as.polygons(rr)
		p$area <- terra::expanse(p, "km") / 1000 
		avga <- sum(p$area) / n
		p$n <- round(p$area / avga)
		totn <- sum(p$n)
		while (totn > n) {
			p$error <- p$area - p$n * avga
			i <- which.max(p$error)
			p$n[i] <- p$n[i] - 1
			totn <- sum(p$n)
		}
		while (totn < n) {
			p$error <- p$area - p$n * avga
			i <- which.min(p$error)
			p$n[i] <- p$n[i] + 1
			totn <- sum(p$n)
		}
		p <- p[p$n > 0, ]
		seeds <- lapply(1:nrow(p), function(i) {
			out <- terra::spatSample(p[i,], p$n[i])
			if (nrow(out) < p$n[i]) {
				out <- terra::spatSample(p[i,], p$n[i] * 10)
				out <- out[sample(nrow(out), p$n[i]), ]
			}
			out
		}) 
			
		seeds <- terra::crds(terra::vect(seeds))
		
		km <- terra::k_means(terra::mask(x, rr), seeds, iter.max = 25)
	} else {

		xm <- terra::mask(x, range)
		km <- terra::k_means(xm, n, iter.max = 25)
	}

	terra::as.polygons(km)
}

get_samplesize <- function(x, fun=NULL, ...) {
	if (inherits(x, "SpatRaster")) {
		x <- terra::as.polygons(x)
		x <- x[x[,1,drop=TRUE]==1]
	} else if (!inherits(x, "SpatVector")) {
		stop("x should be a SpatVector")
	}
	a <- sum(terra::expanse(x, unit="km"))
	if (is.null(fun)) {
		fun <- function(A, omega) max(1, round(omega * sqrt(A/pi)))
	}
	n <- fun(a, ...)
	#n <- max(nmin, n)
	#z <- max(1, min(n, round(a/min_area)))
	return(list(range=x, area=a, n=n))
}



small_ssize_penalty <- function(ssize, score, minssize=10) {
# linear, make curvilinear
	ifelse(ssize > minssize, score, score * ssize / minssize)
}



get_network <- function(regions, sample, maxdist=1500) {

	terra::values(regions) <- data.frame(id=1:nrow(regions))
	patches <- terra::disagg(terra::aggregate(regions))
	patches$pid <- 1:nrow(patches) # pid = patch id 

	x <- terra::centroids(regions, inside=TRUE)

	x$pid <- terra::extract(patches, x)$pid
	patches <- patches[unique(x$pid), ]
	np <- nrow(patches)

	x <- terra::round(x, 5)  # to help merge
	xy <- terra::crds(x)

	adj <- terra::adjacent(regions)
	adj <- data.frame(unique(t(apply(adj, 1, sort))))
	if (ncol(adj) == 0) {
		adj <- data.frame(from=NA, to=NA, adj=NA)
		adj <- adj[0,]
	} else {
		colnames(adj) <- c("from", "to")
		adj$adj <- 1
	}
	if (np > 1) {
		up <- sort(unique(x$pid))
		dx <- as.matrix(terra::distance(x))
		diag(dx) <- NA
		i <- match(colnames(dx), x$id)
		pid <- x$pid[i]
		for (p in up) {
			s <- dx[pid != p, pid == p, drop=FALSE]
			rid <- pid[pid != p]
			upp <- sort(unique(rid))
			for (pp in upp) {
				ss <- s[rid == pp, ,drop=FALSE]
				j <- which.min(apply(ss, 1, min))
				k <- which.min(ss[j, ])
				add <- c(sort(as.integer(c(rownames(ss)[j], colnames(ss)[k]))), 0)
				add <- data.frame(from=add[1], to=add[2], adj=add[3])
				adj <- rbind(adj, add)
			}
		}
		adj <- unique(adj)
	}

	colnames(xy) <- c("xf", "yf")
	adj <- cbind(adj, xy[adj$from, , drop=FALSE])
	colnames(xy) <- c("xt", "yt")
	adj <- cbind(adj, xy[adj$to, , drop=FALSE])

	adj <- cbind(adj, w=1, dst=terra::distance(x[adj$from, ], x[adj$to, ], pairwise=TRUE, unit="m")/1000)

	if (np > 2) {
		padj <- adj[adj$adj == 0, ]
		adj <- adj[adj$adj != 0, ]
		nx <- nrow(x)
		padj <- padj[order(padj$dst), ]
		g <- igraph::components( igraph::graph_from_data_frame(adj, directed = FALSE) )
		for (i in 1:nrow(padj)) {
			adj2 <- rbind(adj, padj[i,])
			gg <- igraph::components( igraph::graph_from_data_frame(adj2, directed = FALSE) )
			if ((gg$no < g$no) | (length(gg$membership) > length(g$membership))) {
				g <- gg
				adj <- adj2
			}
		}
	}

#	mxd <- median(adj$dst) * 3
#	adj$dst[adj$dst > mxd]] <- mxd

	adj$dst[adj$dst > maxdist] <- maxdist	
	adj
}


make_spatvect <- function(rr) {
	x <- as.matrix(rr[, c("xf", "yf", "xt", "yt")])
	a <- lapply(1:nrow(x), \(i)  cbind(i, matrix(x[i,], nrow=2, byrow=2)))
	b <- do.call(rbind, a)
	v <- terra::vect(b, "lines", crs="lonlat")
	terra::values(v) <- rr
	v
}


XC <- function(regions, sample, env=NULL, envfun=NULL, minssize=10, maxdist=1500, return_network=FALSE) {

## TODO  RH
# refine the adjust effect such that when you have many observations in one zones
# they can only contribute to their neighbors. Do not increase branch length to avoid that
# one region does not compensate for another

	stopifnot(minssize > 0)

	rr <- get_network(regions, sample, maxdist)
	if (return_network) {
		return(make_spatvect(rr))
	}
	
	if (!is.null(env)) {
		envd <- as.matrix(EnvDist(env, regions, envfun))
		rr$envdst <- envd[as.matrix(rr[,1:2])]
		rr$geodst <- rr$dst
		# sum geo and env dist
		rr$dst <- rr$dst + rr$envdst
	} 
	gg <- igraph::graph_from_data_frame(rr, directed = FALSE)
	igraph::E(gg)$weight <- rr$dst * rr$w
	y <- unique(terra::extract(regions, sample)[,2])
	rr$w <- rowSums(!matrix(as.matrix(rr[,1:2]) %in% y, ncol=2)) / 2
	igraph::E(gg)$weight2 <- rr$dst * rr$w

	if (nrow(sample) <= 0) {
		return(list(XC=0, dist=rr))
	}


	n <- igraph::count_components(gg)
	score <- nodes <- rep(NA, n)
	dg <- igraph::decompose(gg)
	for (k in 1:n) {
		g <- dg[[k]]
		dst <- igraph::distances(g)
		d1 <- sum(dst) 

		igraph::E(g)$weight <- igraph::E(g)$weight2

		if (length(y) > 1) {
			b <- utils::combn(as.character(y), 2)
			nms <- igraph::V(g)$name
			haveb <- apply(matrix(b %in% nms, nrow=2), 2, all)
			b <- b[,haveb,drop=FALSE]
			if (ncol(b) > 0) {
				for (i in 1:ncol(b)) {
					if (!igraph::are_adjacent(g, b[1,i], b[2,i])) {
						g <- igraph::add_edges(g, b[,i], weight=0)
					}
				}
			}
		}
		d2 <- igraph::distances(g)
		score[k] <- 1 - (sum(d2) / d1)
		nodes[k] <- length(g)
	}
	score <- stats::weighted.mean(score, nodes)
	
	score <- small_ssize_penalty(length(sample), score, minssize)
	list(XC=score, dist=rr)
}

