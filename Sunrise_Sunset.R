library(lubridate)

sunrise_time <- function(longitude, latitude, elevation=0, N = yday(Sys.Date())) {

  # Corrected Solar declination (degrees)
  delta <- 23.45 * sin(((360 * (284 + N) / 365) / 360) * 2 * pi)

  # Hour angle at sunrise and sunset (HA in radians)
  tan_lat <- tan(latitude * pi / 180)
  tan_delta <- tan(delta * pi / 180)
  HA <- acos(-tan_lat * tan_delta)  # Keeps in radians

  # Elevation correction (radians)
  h_elev <- 0.032 * sqrt(elevation) * pi / 180

  # Atmospheric refraction correction (radians)
  h_refr <- 0.566 * pi / 180

  # Total horizon correction (radians)
  h_total <- h_elev + h_refr

  # Adjusted HA for elevation and refraction
  HA_adj <- HA + h_total

  # Convert HA to hours (15° per hour)
  HA_hours <- HA_adj * 180 / (15 * pi)

  # Equation of Time (EOT in minutes)
  B <- (360 / 365) * (N - 81) * pi / 180
  EOT <- 9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B)

  # Sunrise time offset from hour angle
  sunrise_utc <- 12 - HA_hours - (longitude / 15) + (EOT / 60)

  # Return
  return(sunrise_utc)
}


sunset_time <- function(longitude, latitude, elevation=0, N = yday(Sys.Date())) {
  # Corrected Solar declination (degrees)
  delta <- 23.45 * sin(((360 * (284 + N) / 365) / 360) * 2 * pi)

  # Hour angle at sunrise and sunset (HA in radians)
  tan_lat <- tan(latitude * pi / 180)
  tan_delta <- tan(delta * pi / 180)
  HA <- acos(-tan_lat * tan_delta)  # Keeps in radians

  # Elevation correction (radians)
  h_elev <- 0.032 * sqrt(elevation) * pi / 180

  # Atmospheric refraction correction (radians)
  h_refr <- 0.566 * pi / 180

  # Total horizon correction (radians)
  h_total <- h_elev + h_refr

  # Adjust HA
  HA_sunset  <- HA + h_total

  # Convert HA to hours (15° per hour)
  HA_sunset_hours  <- HA_sunset * 180 / (15 * pi)

  # Equation of Time (EOT in minutes)
  B <- (360 / 365) * (N - 81) * pi / 180
  EOT <- 9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B)

  # Sunset time offset from hour angle
  sunset_utc  <- 12 + HA_sunset_hours - (longitude / 15) + (EOT / 60)

  # Return
    return(sunset_utc)
}

library(raster)
# load DTM
ire<-raster('~/IrelandGIS/Ireland_GTM_WGS.tif')

library(terra)

sunrise_raster<-function(day=yday(Sys.Date()), npix=1500, dem=ire){
   # Generate X and Y sequences
   x_vals <- seq(extent(dem)[1], extent(dem)[2], length.out = npix)  # Longitude
   y_vals <- seq(extent(dem)[4], extent(dem)[3], length.out = npix)  # Latitude (reversed!)

   # Create a grid of coordinates
   coords <- as.matrix(expand.grid(x_vals, y_vals))

   # Extract elevation values in bulk
   elevations <- extract(dem, coords)

   # Compute sun times in bulk
   hapoints <- matrix(mapply(sunrise_time, coords[, 1], coords[, 2], elevations, MoreArgs = list(day)),
                      nrow = npix, ncol = npix, byrow = TRUE)

   return(raster(extent(dem), nrows=npix, ncols=npix, crs=crs(dem), vals=hapoints))
}

sunset_raster<-function(day=yday(Sys.Date()), npix=1500, dem=ire){
  # Generate X and Y sequences
  x_vals <- seq(extent(dem)[1], extent(dem)[2], length.out = npix)  # Longitude
  y_vals <- seq(extent(dem)[4], extent(dem)[3], length.out = npix)  # Latitude (reversed!)

  # Create a grid of coordinates
  coords <- as.matrix(expand.grid(x_vals, y_vals))

  # Extract elevation values in bulk
  elevations <- extract(dem, coords)

  # Compute sun times in bulk
  hapoints <- matrix(mapply(sunset_time, coords[, 1], coords[, 2], elevations, MoreArgs = list(day)),
                     nrow = npix, ncol = npix, byrow = TRUE)

  return(raster(extent(dem), nrows=npix, ncols=npix, crs=crs(dem), vals=hapoints))
}

sunrise141<-sunrise_raster(day=141)
sunset141<-sunset_raster(day=141)


plot(sunrise141, col=hcl.colors(128,'Inferno',rev=TRUE))
plot(sunset141, col=hcl.colors(128,'Inferno'))

raster::writeRaster(sunrise141, '~/IrelandGIS/sunrise141.tiff', format='GTiff', overwrite=TRUE)
raster::writeRaster(sunset141, '~/IrelandGIS/sunset141.tiff', format='GTiff', overwrite=TRUE)


# High res sunrise for April 26th
sunrise85<-sunrise_raster(day=116,npix=2000)
plot(sunrise85, col=hcl.colors(128,'Inferno',rev=TRUE))

# Load the coastline polygon (e.g., from a shapefile)
coastline <- vect("coastline.shp")  # Make sure it's a polygon, not just lines!

iremask<-function(myraster, coast=coastline) {
  X <- terra::rast(myraster)
  raster_masked <- terra::mask(X, coast)
  raster_cropped <- terra::crop(raster_masked, coast)
  return(raster_cropped)
}

maxpoint<-function(X) {
  max_value <- max(values(X), na.rm = TRUE)
  max_cell <- which.max(values(X))
  max_coords <- xyFromCell(X, max_cell)
  points(max_coords[1], max_coords[2], pch = 16, col = "red", cex = 2)
}

minpoint<-function(X) {
  min_value <- min(values(X), na.rm = TRUE)
  min_cell <- which.min(values(X))
  min_coords <- xyFromCell(X, min_cell)
  points(min_coords[1], min_coords[2], pch = 16, col = 'red', cex = 2)
}

# Plot to check results
plot(iremask(sunrise85), col=hcl.colors(128,'Inferno',rev=TRUE))
maxpoint(sunset_raster_masked)
minpoint(sunset_raster_masked)

# PDF containing one page of sunrise per day of year
pdf('sunrises.pdf')
for(day in 1:365){
  sunrise<-sunrise_raster(day=day,npix=600)
  masked_sunrise<-iremask(sunrise)
  plot(masked_sunrise, col=hcl.colors(128,'Inferno',rev=TRUE))
  minpoint(masked_sunrise)
  # Convert day number to date
  dtxt<-format(as.Date(day, origin = paste0(2025 - 1, "-12-31")), "%B %d")
  title(paste('Sunrise,',dtxt))
}
dev.off()

# PDF containing one page of sunset per day of year
pdf('sunsets.pdf')
for(day in 1:365){
  sunset<-sunset_raster(day=day,npix=1200)
  masked_sunset<-iremask(sunset)
  plot(masked_sunset, col=hcl.colors(128,'Inferno'))
  maxpoint(masked_sunset)
  # Convert day number to date
  dtxt<-format(as.Date(day, origin = paste0(2025 - 1, "-12-31")), "%B %d")
  title(paste('Sunset,',dtxt))
}
dev.off()

