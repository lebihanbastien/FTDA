#------------------------------------------------
# Load libraries
#------------------------------------------------
library(plyr)
library(ggplot2)
library(reshape2)
library(plotly)
source("multiplot.R")


#------------------------------------------------
# Plotly
#------------------------------------------------
py <- plotly(username="b.le-bihan", key="ado73xmcij")  # open plotly connection

#------------------------------------------------
# Data
#------------------------------------------------
#The orbit
orbit = read.table("orbit.txt", header = F, sep = "")
colnames(orbit) = c("x", "y", "z")

#The events
events = read.table("events.txt", header = F, sep = "")
colnames(events) = c("x", "y", "z")

#------------------------------------------------
# 3D plot
#------------------------------------------------
list1 = list(
  x = orbit$x,
  y = orbit$y,
  z = orbit$z,
  mode = "lines",
  marker = list(
    color = "rgb(127, 127, 127)",
    size = 12,
    symbol = "circle",
    line = list(
      color = "rgb(204, 204, 204)",
      width = 1
    ),
    opacity = 0.9
  ),
  type = "scatter3d"
)

list2 = list(
  x = events$x,
  y = events$y,
  z = events$z,
  mode = "markers",
  marker = list(
    color = "rgb(127, 127, 127)",
    size = 12,
    symbol = "circle",
    line = list(
      color = "rgb(204, 204, 204)",
      width = 1
    ),
    opacity = 0.9
  ),
  type = "scatter3d"
)

data <- list(list1, list2)


layout <- list(margin = list(
  l = 0,
  r = 0,
  b = 0,
  t = 0))

response <- py$plotly(data, kwargs=list(layout=layout, filename="random-walk", fileopt="overwrite"))
url <- response$url