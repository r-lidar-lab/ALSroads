![Lifecycle:Experimental](https://img.shields.io/badge/Lifecycle-Experimental-339999)
# ALSroads

ALSroads provides tools to relocate, measure and estimate the state of forestry roads from an inaccurate map. This project has been made in partnership with the [Ministère des Forêts, de la Faune et des Parcs du Québec (MFFP)](https://mffp.gouv.qc.ca/).  This package is experimental and must be considered as a proof of concept. It contains the implementation of the algorithm described in:

> Roussel, J., Bourdon, J., Morley, I. D., Coops, N. C., & Achim, A. (2022). Correction, update, and enhancement of vectorial forestry road maps using ALS data, a pathfinder, and seven metrics. International Journal of Applied Earth Observation and Geoinformation, 114(September), 103020. https://doi.org/10.1016/j.jag.2022.103020


## Installation

``` r
remotes::install_github("metafor-ulaval/ALSroads")
```

## Example

In the following example we can see a road from reference forestry road maps in red. This road is inaccurately mapped and records neither its class nor its width. The algorithm recomputes the accurate location of the road from a lidar point cloud and estimates its width and its state. Here we have an class 1 road with a width of 8 meters.

```r
library(ALSroads)
library(lidR)
library(sf)
library(raster)
library(mapview)
library(leaflet)

# Load data (LAS tiles, DTM, map)
dir  <- system.file("extdata", "", package="ALSroads")
road <- system.file("extdata", "j5gr_centerline_971487.gpkg", package="ALSroads")
dtm  <- system.file("extdata", "j5gr_dtm.tif", package="ALSroads")
ctg  <- readLAScatalog(dir)
road <- st_read(road, quiet = TRUE)
dtm  <- raster(dtm)
crs(dtm) <- crs(road)

# Relocate the road at the correct location
# Measure the width and estimate its state (operating/decommisionned)
res <- measure_road(ctg, road, dtm)
poly <- sf::st_buffer(res, res$ROADWIDTH/2)

# Display results
url = "https://servicesmatriciels.mern.gouv.qc.ca:443/erdas-iws/ogc/wmts/Inventaire_Ecoforestier/Inventaire_Ecoforestier/default/GoogleMapsCompatibleExt2:epsg:3857/{z}/{y}/{x}.jpg"
m = mapview::mapview(list(road, poly),
  layer.name = c("Inaccurate", "Corrected"),
  color = c("red", "blue"), map.type = "Esri.WorldImagery")
leaflet::addTiles(m@map, url)
```

![](inst/extdata/screenshot.png)

