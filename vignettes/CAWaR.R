## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE, message=FALSE, warning=FALSE----------------------------
# load packages
library(fieldRS)
library(raster)
library(ggplot2)
library(knitr)
library(kableExtra)
library(CAWaR)

## ----message=FALSE-------------------------------------------------------
ndvi.ts <- brick(system.file("extdata", "ndvi.tif", package="fieldRS")) # NDVI raster time series
data(fieldData) # ground truth data
data(referenceProfiles) # target crop types NDVI profiles

## ----echo=FALSE----------------------------------------------------------
ndvi.ts <- extend(ndvi.ts, 30)

## ------------------------------------------------------------------------
sampleTest <- checkSamples(fieldData)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
kable_styling(kable(head(sampleTest, 1), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")

## ------------------------------------------------------------------------
sampleTest <- geCheck(fieldData)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
kable_styling(kable(head(sampleTest, 1), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")

## ------------------------------------------------------------------------
# build polygons
p1 <- Polygons(list(Polygon(data.frame(x=c(1, 5, 10, 2, 1), y=c(10, 9, 8, 7, 10)))), ID=1)
p2 <- Polygons(list(Polygon(data.frame(x=c(2, 6, 5, 4, 2), y=c(10, 9, 7, 4, 10)))), ID=2)
p <- SpatialPolygons(list(p1, p2))

# check overlap
sampleTest <- geCheck(p)

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
kable_styling(kable(head(sampleTest$overlap.df, 1), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")

plot(p)
plot(sampleTest$overlap.shp, col="red", add=TRUE)

## ------------------------------------------------------------------------
sampleCorrect <- labelCheck(fieldData$crop)
sampleCorrect$labels # unique labels

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
kable_styling(kable(head(sampleCorrect$label.count, 3), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")
sampleCorrect$label.count.plot

## ------------------------------------------------------------------------
sampleCorrect <- labelCheck(fieldData$crop, sampleCorrect$labels, c("wheat", "not-wheat", "not-wheat"))
fieldData@data$crop_2 <- sampleCorrect$labels

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
kable_styling(kable(head(sampleCorrect$label.count, 3), format="html", align="c", full_width=TRUE), "stripped", bootstrap_options="responsive")
sampleCorrect$label.count.plot

## ---- out.width="98%", fig.height=5, fig.width=10, dpi=600, fig.align="center", fig.show='hold', echo=FALSE----
include_graphics("percentCover.jpg")

