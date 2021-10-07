library(marmap)
library(ggplot2)
library(scico)

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Snail_16S_amplicons/Meta-analysis")

dat <- getNOAA.bathy(50, -165, 30, -30, res=4, keep=TRUE, antimeridian=TRUE)

coord <- data.frame(lon = c(60.53, 144.71, 144.70, 144.72, 144.87, 144.51, 144.71, 155.05, 151.67, 151.67, 151.67, 151.67, 152.10, 152.10, 152.11, 173.92, 182.81, 182.85, 186.43, 185.35, 183.86, 183.82, 183.81, 183.44), lat = c(6.36, 18.21, 18.23, 18.18, 16.96, 15.48, 20.82, -9.80, -3.72, -3.73, -3.73, -3.72, -3.80, -3.79, -3.81, -16.95, -14.76, -14.75, -15.16, -15.41, -20.32, -20.67, -20.76, -21.99))
site <- data.frame(lon = c(60.53, 144.72, 144.71, 155.05, 151.67, 173.92, 182.85, 186.43, 185.35, 183.82), lat = c(6.36, 18.18, 20.82, -9.80, -3.73, -16.95, -14.75, -15.16, -16.00, -20.67), basin = c("Carlsberg Ridge", "Mariana Back-Arc", "Mariana Arc", "Woodlark", "Manus Basin", "North Fiji Basin", "Futuna", "Tonga Arc", "NELSC", "ELSC"))

pdf("Map.pdf")
autoplot(dat, geom="r", coast=FALSE) + scale_fill_scico(palette = 'vik') + geom_point(data = coord, aes(x = lon, y = lat), shape = 21, color = "black", fill = "darkred", size = 2) + labs(x = "Longitude", y = "Latitude") + theme(legend.position = "right") + labs(fill = "Elevation (m)") + geom_text(data = site, aes(x = lon, y = lat, label = basin), size=2.5, fontface="bold", nudge_x = c(1.5, -14, 10, -1.5, 1.5, -6, -7, 1.5, 4, 6), nudge_y = c(-3, 0, 0, -3, 3, -3, 1, 3, -2, -1))
dev.off()

