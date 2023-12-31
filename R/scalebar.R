

gcDestination <- function(lon, lat, bearing, dist, dist.units = "km",
    model=NULL, Vincenty=FALSE) {
 # lat, lon : lattitude and longitude in decimal degrees
 # bearing : bearing from 0 to 360 degrees
 # dist : distance travelled
 # dist.units : units of distance "km" (kilometers), "nm" (nautical 
 # miles), "mi" (statute miles)
 # model : choice of ellipsoid model ("WGS84", "GRS80", "Airy", 
 # "International", "Clarke", "GRS67")

    if (!is.numeric(lon)) stop("lon not numeric")
    if (!is.numeric(lat)) stop("lat not numeric")
    if (!is.numeric(bearing)) stop("bearing not numeric")
    if (!is.numeric(dist)) stop("dist not numeric")

    if (length(lon) != length(lat)) stop("lon and lat differ in length")
    if (length(bearing) > 1L && length(lon) > 1L) stop("length mismatch")
    if (length(bearing) > 1L && length(dist) > 1L) stop("length mismatch")

    as.radians <- function(degrees) degrees * pi / 180
    as.degrees <- function(radians) radians * 180 / pi
    as.bearing <- function(radians) (as.degrees(radians) + 360) %% 360

    ellipsoid <- function(model = "WGS84") {
        switch(model,
        WGS84 = c(a = 6378137, b = 6356752.3142, f = 1 / 298.257223563),
        GRS80 = c(a = 6378137, b = 6356752.3141, f = 1 / 298.257222101),
        Airy = c(a = 6377563.396, b = 6356256.909, f = 1 / 299.3249646),
        International = c(a = 6378888, b = 6356911.946, f = 1 / 297),
        Clarke = c(a = 6378249.145, b = 6356514.86955, f = 1 / 293.465),
        GRS67 = c(a = 6378160, b = 6356774.719, f = 1 / 298.25),
        c(a = NA, b = NA, f = NA)
    )}
    
    dist <- switch(dist.units,
        km = dist,
        nm = dist * 1.852,
        mi = dist * 1.609344
    )
    lat <- as.radians(lat)
    lon <- as.radians(lon)
    bearing <- as.radians(bearing)

    if (is.null(model)) {
 # Code adapted from JavaScript by Chris Veness 
 # (scripts@movable-type.co.uk) at
 # http://www.movable-type.co.uk/scripts/latlong.html#ellipsoid
 #   originally from Ed Williams' Aviation Formulary, 
 # http://williams.best.vwh.net/avform.htm
        radius <- 6371
        psi <- dist / radius
        lat2 <- asin(sin(lat) * cos(psi) +  cos(lat) * sin(psi) * cos(bearing))
        lon2 <- lon + atan2(sin(bearing) * sin(psi) * cos(lat), cos(psi) - 
            sin(lat) * sin(lat2))
        if (any(is.nan(lat2)) || any(is.nan(lon2))) warning("Out of range values")
        return(cbind(long=as.degrees(lon2), lat=as.degrees(lat2)))
    }

    ellips <- ellipsoid(model)
    if (is.na(ellips["a"])) stop("no such ellipsoid model")
    if (Vincenty) {
 # Code adapted from JavaScript by Chris Veness 
 # (scripts@movable-type.co.uk) at
 # http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
 # Original reference (http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf):
 #   Vincenty, T. 1975.  Direct and inverse solutions of geodesics on 
 # the ellipsoid with application of nested equations.
 #      Survey Review 22(176):88-93
        dist <- dist * 1000
        sin.alpha1 <- sin(bearing)
        cos.alpha1 <- cos(bearing)
        tan.u1 <- (1 - ellips["f"]) * tan(lat)
        cos.u1 <- 1 / sqrt(1 + (tan.u1 ^ 2))
        sin.u1 <- tan.u1 * cos.u1
        sigma1 <- atan2(tan.u1, cos.alpha1)
        sin.alpha <- cos.u1 * sin.alpha1
        cos.sq.alpha <- 1 - (sin.alpha ^ 2)
        u.sq <- cos.sq.alpha * ((ellips["a"] ^ 2) - (ellips["b"] ^ 2)) /
            (ellips["b"] ^ 2)
        cap.A <- 1 + u.sq / 16384 * (4096 + u.sq * (-768 + u.sq * (320 - 
            175 * u.sq)))
        cap.B <- u.sq / 1024 * (256 + u.sq * (-128 + u.sq * (74 - 47 * u.sq)))

        sigma <- dist / (ellips["b"] * cap.A)
        sigma.p <- 2 * pi
        cos.2.sigma.m <- cos(2 * sigma1 + sigma)
        while(any(abs(sigma - sigma.p) > 1e-12)) {
            cos.2.sigma.m <- cos(2 * sigma1 + sigma)
            sin.sigma <- sin(sigma)
            cos.sigma <- cos(sigma)
            delta.sigma <- cap.B * sin.sigma * (cos.2.sigma.m + cap.B / 4 * 
                (cos.sigma *
                (-1 + 2 * cos.2.sigma.m ^ 2) - cap.B / 6 * cos.2.sigma.m *
                (-3 + 4 * sin.sigma ^ 2) * (-3 + 4 * cos.2.sigma.m ^ 2)))
            sigma.p <- sigma
            sigma <- dist / (ellips["a"] * cap.A) + delta.sigma
        }
        tmp <- sin.u1 * sin.sigma - cos.u1 * cos.sigma * cos.alpha1
        lat2 <- atan2(sin.u1 * cos.sigma + cos.u1 * sin.sigma * cos.alpha1,
            (1 - ellips["f"]) * sqrt(sin.alpha ^ 2 + tmp ^ 2))
        lambda <- atan2(sin.sigma * sin.alpha1, cos.u1 * cos.sigma - sin.u1 * 
            sin.sigma * cos.alpha1)
        cap.C <- ellips["f"] / 16 * cos.sq.alpha * (4 + ellips["f"] * 
            (ellips["f"] - 3 * cos.sq.alpha))
        cap.L <- lambda - (1 - cap.C) * ellips["f"] * sin.alpha *
            (sigma + cap.C * sin.sigma * (cos.2.sigma.m + cap.C * cos.sigma * 
            (-1 + 2 * cos.2.sigma.m ^ 2)))
        lat2 <- as.degrees(lat2)
        lon2 <- as.degrees(lon + cap.L)
    } else {
 # Code adapted from JavaScript by Larry Bogan (larry@go.ednet.ns.ca) 
 # at http://www.go.ednet.ns.ca/~larry/bsc/jslatlng.html
        e <- 0.08181922
        radius <- (ellips["a"] / 1000) * (1 - e^2) / ((1 - e^2 * 
            sin(lat)^2)^1.5)
        psi <- dist / radius
        phi <- pi / 2 - lat
        arc.cos <- cos(psi) * cos(phi) + sin(psi) * sin(phi) * cos(bearing)
        lat2 <- as.degrees((pi / 2) - acos(arc.cos))
        arc.sin <- sin(bearing) * sin(psi) / sin(phi)
        lon2 <- as.degrees(lon + asin(arc.sin))
    }
    return(cbind(long=lon2, lat=lat2))
}

#' scalebar
#' @description Adds a scale bar to maps created with ggplot or ggmap.
#' @param data the same \code{\link{data.frame}} passed to \code{\link{ggplot}} to plot the map. If the \code{class} of \code{data} is not \code{sf}, it must contain columns whose names begin with \code{long} and \code{lat}.
#' @param location string indicating the scale bar's location in the plot. Possible options: "topright" (default), "bottomright", "bottomleft" and "topleft".
#' @param dist distance to represent with each segment of the scale bar.
#' @param dist_unit unit of measurement for \code{dist}. Possbile values: "km" (kilometers) and "m" (meters), "nm" (nautical miles) and "mi" (statue miles).
#' @param transform If TRUE, it is assumed that coordinates are in decimal degrees. If FALSE, it assumed that they are in meters.
#' @param dd2km deprecated. Use transform instead.
#' @param model choice of ellipsoid model ("WGS84", "GRS80", "Airy", "International", "Clarke", or "GRS67") Used when transform is TRUE.
#' @param height number between 0 and 1 to indicate the scale bar's height, as a proportion of the y axis.
#' @param st.dist number between 0 and 1 to indicate the distance between the scale bar and the scale bar's text, as a proportion of the y axis.
#' @param st.bottom logical. If TRUE (default) the scale bar's text is displayed at the bottom of the scale bar, if FALSE, it is displayed at the top.
#' @param st.size number to indicate the scale bar's size. It is passed to \code{size} in \code{\link{annotate}} function.
#' @param st.color color of the scale bar's text. Default is black.
#' @param box.fill fill color of the box. If vector of two colors, the two boxes are filled with a different color. Defaults to black and white.
#' @param box.color color of the box's border. If vector of two colors, the borders of the two boxes are colored differently. Defaults to black.
#' @param border.size number to define the border size.
#' @param anchor named \code{\link{vector}} with coordinates to control the symbol's position. For \code{location = "topright"}, \code{anchor} defines the coordinates of the symbol's topright corner and so forth. The x coordinate must be named as x and the y coordinate as y.
#' @param x.min if \code{data} is not defined, number with the minimum x coordinate. Useful for ggmap.
#' @param x.max if \code{data} is not defined, number with the maximum x coordinate. Useful for ggmap.
#' @param y.min if \code{data} is not defined, number with the minimum y coordinate. Useful for ggmap.
#' @param y.max if \code{data} is not defined, number with the maximum y coordinate. Useful for ggmap.
#' @param facet.var if faceting, character vector of variable names used for faceting. This is useful for placing the scalebar in only one facet and must be used together with \code{facet.lev}.
#' @param facet.lev character vector with the name of one level for each variable in \code{facet.var}. The scale bar will be drawn only in the \code{facet.lev} facet.
#' @param st.inherit logical. Set as FALSE if scalebar has unexpected behavior in animations.
#' @export
#' @examples
#' \dontrun{
#' library(sf)
#' data(domestic_violence)
#'
#' # Map in geographic coordinates
#' ggplot(domestic_violence, aes(fill = Scaled)) +
#'     geom_sf() +
#'     scalebar(domestic_violence, dist = 4, dist_unit = "km",
#'              transform = TRUE, model = "WGS84") +
#'     blank() +
#'     scale_fill_continuous(low = "#fff7ec", high = "#7F0000")
#' 
#' # Map in projected coordinates
#' domestic_violence2 <- st_transform(domestic_violence, 31983)
#' ggplot(domestic_violence2, aes(fill = Scaled)) +
#'     geom_sf() +
#'     scalebar(domestic_violence2, dist = 4, dist_unit = "km",
#'              transform = FALSE, model = "WGS84") +
#'     blank() +
#'     scale_fill_continuous(low = "#fff7ec", high = "#7F0000")
#' }
scalebar <- function(data = NULL, location = "bottomright", dist = NULL, dist_unit = NULL, transform = NULL, dd2km = NULL, model = NULL, height = 0.02, st.dist = 0.02, st.bottom = TRUE, st.size = 5, st.color = "black", box.fill = c("black", "white"), box.color = "black", border.size = 1, x.min = NULL, x.max = NULL, y.min = NULL, y.max = NULL, anchor = NULL, facet.var = NULL, facet.lev = NULL, st.inherit = TRUE){
    if (is.null(data)) {
        if (is.null(x.min) | is.null(x.max) |
            is.null(y.min) | is.null(y.max) ) {
            stop('If data is not defined, x.min, x.max, y.min and y.max must be.')
        }
        data <- data.frame(long = c(x.min, x.max), lat = c(y.min, y.max))
    }
    if (is.null(transform)) {
        stop("transform should be logical.")
    }
    if (any(class(data) %in% "sf")) {
        xmin <- sf::st_bbox(data)["xmin"]
        xmax <- sf::st_bbox(data)["xmax"]
        ymin <- sf::st_bbox(data)["ymin"]
        ymax <- sf::st_bbox(data)["ymax"]
    } else {
        if (any(startsWith(colnames(data), "lat")) & any(startsWith(colnames(data), "long"))) {
            xmin <- min(data$long)
            xmax <- max(data$long)
            ymin <- min(data$lat)
            ymax <- max(data$lat)
        } else {
            stop("'", substitute(data), "' must have columns with names that start with 'lat' and 'long'")
        }
    }
    if (location == 'bottomleft') {
        if (is.null(anchor)) {
            x <- xmin
            y <- ymin
        } else {
            x <- as.numeric(anchor['x'])
            y <- as.numeric(anchor['y'])
        }
        direction <- 1
        
    }
    if (location == 'bottomright') {
        if (is.null(anchor)) {
            x <- xmax
            y <- ymin
        } else {
            x <- as.numeric(anchor['x'])
            y <- as.numeric(anchor['y'])
        }
        direction <- -1
        
    }
    if (location == 'topleft') {
        if (is.null(anchor)) {
            x <- xmin
            y <- ymax
        } else {
            x <- as.numeric(anchor['x'])
            y <- as.numeric(anchor['y'])
        }
        direction <- 1
        
    }
    if (location == 'topright') {
        if (is.null(anchor)) {
            x <- xmax
            y <- ymax
        } else {
            x <- as.numeric(anchor['x'])
            y <- as.numeric(anchor['y'])
        }
        direction <- -1
        
    }
    if (!st.bottom) {
        st.dist <-
            y + (ymax - ymin) * (height + st.dist)
    } else {
        st.dist <- y - (ymax - ymin) * st.dist
    }
    height <- y + (ymax - ymin) * height
    
    if (dist_unit == "m") {
        dist <- dist / 1e3
        dist_unit0 <- "m"
        dist_unit <- "km"
    }
    if (!is.null(dd2km)) {
        if (dd2km) {
            transform <- TRUE
        }
        cat("dd2km is deprecated. Use ggsn::transform instead.")
    }
    if (transform) {
        break1 <- gcDestination(lon = x, lat = y,
                                          bearing = 90 * direction,
                                          dist = dist, dist.units = dist_unit,
                                          model = model)[1, 1]
        break2 <- gcDestination(lon = x, lat = y,
                                          bearing = 90 * direction,
                                          dist = dist*2, dist.units = dist_unit,
                                          model = model)[1, 1]
    } else {
        if (location == 'bottomleft' | location == 'topleft') {
            if (exists("dist_unit0") | (!exists("dist_unit0") & dist_unit == "km")) {
                break1 <- x + dist * 1e3
                break2 <- x + dist * 2e3
            } else if (dist_unit == "nm") {
                break1 <- x + dist * 1852
                break2 <- x + dist * 1852 * 2
            } else if (dist_unit == "mi") {
                break1 <- x + dist * 1609.34
                break2 <- x + dist * 1609.34 * 2
            } else {
                break1 <- x + dist
                break2 <- x + dist
            }
        } else {
            if (exists("dist_unit0") | (!exists("dist_unit0") & dist_unit == "km")) {
                break1 <- x - dist * 1e3
                break2 <- x - dist * 2e3
            } else if (dist_unit == "nm") {
                break1 <- x - dist * 1852
                break2 <- x - dist * 1852 * 2
            } else if (dist_unit == "mi") {
                break1 <- x - dist * 1609.34
                break2 <- x - dist * 1609.34 * 2
            } else {
                break1 <- x - dist
                break2 <- x - dist
            }
        }
        
    }
    
    out_of_range <- function(low, n, high) {
        n < low | n > high 
    }
    
    if (out_of_range(xmin, break1, xmax) | out_of_range(xmin, break2, xmax)) {
        stop("The requested scalebar distance (", 
             substitute(dist), " ", substitute(dist_unit), 
             ") is too large to fit on the map.\n  Try reducing it.")
    }
    
    
    box1 <- data.frame(x = c(x, x, rep(break1, 2), x),
                       y = c(y, height, height, y, y), group = 1)
    box2 <- data.frame(x = c(rep(break1, 2), rep(break2, 2), break1),
                       y=c(y, rep(height, 2), y, y), group = 1)
    if (!is.null(facet.var) & !is.null(facet.lev)) {
        for (i in 1:length(facet.var)){
            if (any(class(data) == "sf")) {
                if (!is.factor(data[ , facet.var[i]][[1]])) {
                    data[ , facet.var[i]] <- factor(data[ , facet.var[i]][[1]])
                }
                box1[ , facet.var[i]] <- factor(facet.lev[i],
                                                levels(data[ , facet.var[i]][[1]]))
                box2[ , facet.var[i]] <- factor(facet.lev[i],
                                                levels(data[ , facet.var[i]][[1]]))
            } else {
                if (!is.factor(data[ , facet.var[i]])) {
                    data[ , facet.var[i]] <- factor(data[ , facet.var[i]])
                }
                box1[ , facet.var[i]] <- factor(facet.lev[i],
                                                levels(data[ , facet.var[i]]))
                box2[ , facet.var[i]] <- factor(facet.lev[i],
                                                levels(data[ , facet.var[i]]))
            }
            
        }
    }
    if (exists("dist_unit0")) {
        legend <- cbind(text = c(0, dist * 1e3, dist * 2e3), row.names = NULL)
    } else {
        legend <- cbind(text = c(0, dist, dist * 2), row.names = NULL)
    }
    gg.box1 <- geom_polygon(data = box1, aes(x, y),
                            fill = utils::tail(box.fill, 1),
                            color = utils::tail(box.color, 1),
                            size = border.size)
    gg.box2 <- geom_polygon(data = box2, aes(x, y), fill = box.fill[1],
                            color = box.color[1],
                            size = border.size)
    x.st.pos <- c(box1[c(1, 3), 1], box2[3, 1])
    if (location == 'bottomright' | location == 'topright') {
        x.st.pos <- rev(x.st.pos)
    }
    label <- NULL
    if (exists("dist_unit0")) {
        legend2 <- cbind(data[1:3, ], x = unname(x.st.pos), y = unname(st.dist),
                         label = paste0(legend[, "text"], c("", "", "m")))
    } else {
        legend2 <- cbind(data[1:3, ], x = unname(x.st.pos), y = unname(st.dist),
                         label = paste0(legend[, "text"], c("", "", dist_unit)))
    }
    if (!is.null(facet.var) & !is.null(facet.lev)) {
        for (i in 1:length(facet.var)){
            if (any(class(data) == "sf")) {
                legend2[ , facet.var[i]] <- factor(facet.lev[i],
                                                   levels(data[ , facet.var[i]][[1]]))
            } else {
                legend2[ , facet.var[i]] <- factor(facet.lev[i],
                                                   levels(data[ , facet.var[i]]))
            }
        }
    } else if (!is.null(facet.var) & is.null(facet.lev)) {
        facet.levels0 <- unique(as.data.frame(data)[, facet.var])
        facet.levels <- unlist(unique(as.data.frame(data)[, facet.var]))
        legend2 <- do.call("rbind", replicate(length(facet.levels),
                                              legend2, simplify = FALSE))
        if (length(facet.var) > 1) {
            facet.levels0 <- expand.grid(facet.levels0)
            legend2[, facet.var] <-
                facet.levels0[rep(row.names(facet.levels0), each = 3), ]
        } else {
            legend2[, facet.var] <- rep(facet.levels0, each = 3)
        }
    }
    if (!st.inherit) {
        legend2 <- legend2[, c("x", "y", "label")]
    }
    gg.legend <- geom_text(data = legend2, aes(x, y, label = label),
                           size = st.size, color = st.color,
                           inherit.aes = st.inherit)
    return(list(gg.box1, gg.box2, gg.legend))
}
