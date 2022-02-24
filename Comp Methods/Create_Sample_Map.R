###USE MARMAP TO CREATE A BATHYMETRY MAP OF OUR STUDY AREA
#Creates bathymetry map of the florida keys and then adds points for the two reefs
#our ancient cores are sourced from. The rest of this map editing was completed 
#in power
library(marmap)
library(maps)
library(usmap)
library(sf)
library(sp)
library("rnaturalearthhires")
library("rnaturalearth")
library("rnaturalearthdata")
library(RColorBrewer)

ancient_meta <- as.data.frame(matrix(NA, nrow = 4, ncol = 6))
colnames(ancient_meta) <- c("Species", "Region", "Reef", "Latitude", "Longitude", "Study")
ancient_meta$Species <- c("A. palmata", "A. spp", "A. palmata", "A. palmata")
ancient_meta$Region <- "Florida"
ancient_meta$Reef <- c("Looe Key", "Looe Key", "Sombrero Reef", "Sombrero Reef")
ancient_meta$Latitude <- c(rep(24.5418, 2), rep(24.6234, 2))
ancient_meta$Longitude <- c(rep(-81.4039, 2), rep(-81.10696, 2))
ancient_meta$Study <- "aDNA"

world <- ne_countries(scale = 10, returnclass = "sf")
blues <- c("lightblue4", "lightblue3", "lightblue2", "lightblue1",
           "lightblue")
blues <- c("skyblue4", "skyblue3", "skyblue2", "skyblue1")
blues <- rev(brewer.pal(6, "Blues"))[c(-1,-6)]
papoue <- getNOAA.bathy(lon1 = -82, lon2 = -80.2,
                        lat1 = 24.45, lat2 = 25.4, resolution = 1)

dev.off()
png("~/Box Sync/Projects/aDNA_final_docs/FL_keysMap.png", height = 5, width = 7, res = 500, units = "in")
plot(papoue, land = F, image = T, bpal = list(c(0, max(papoue), "grey"),
                                              c(min(papoue),0,blues)), lwd = 0)
plot(world, add = T, col = 'white', border = 'white')
points(ancient_meta$Longitude, ancient_meta$Latitude, pch = 19, cex = 1.5)
dev.off()
