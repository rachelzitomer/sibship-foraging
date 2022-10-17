#!/usr/bin/env Rscript
library(raster)
library(sp)

dir.create("../data")
for (species in c("bomvos", "bomcal"))
{
  Species <- c("bomvos"="Bvos", "bomcal"="Bcal")[species]

  for (region in c("north", "south", "midnorth"))
  {
    floral_cover_at_traps <- as.matrix(read.table(paste0("floral_cover_", region, ".txt"), header=TRUE))
    rownames(floral_cover_at_traps) <- sub(" ", "_", rownames(floral_cover_at_traps))

    colony_count_at_traps <- as.matrix(read.csv(paste0("genotyped_captures_", Species, "_", region, ".csv")))

    capture_count_at_traps <- read.csv(paste0("total_captures_", Species, "_", region, ".csv"))
    rownames(capture_count_at_traps) <- sub(" ", "_", capture_count_at_traps[,2])
    capture_count_at_traps <- as.matrix(capture_count_at_traps[,c("n_F", "n_M", "total")])
    colnames(capture_count_at_traps) <- c("females", "males", "total")

    trap_coordinates <- read.table(paste0(region, "_trap_coords.txt"), header=TRUE)
    rownames(trap_coordinates) <- sub(" ", "_", trap_coordinates[,1])
    trap_coordinates <- as.matrix(trap_coordinates[,c("x", "y")])

    for (resolution in c("30m", "60m", "90m"))
    {
      decorator <- c("30m"="binary", "60m"="proportion", "90m"="proportion")[resolution]
      stand_age <- readAll(raster(paste0(region, "_agemap_", resolution, ".tif")))
      roads <- readAll(raster(paste0(region, "_roads_", decorator, "_", resolution, ".tif")))
      canopy_cover<-readAll(raster(paste0(region,'_canopy_cover_map_', resolution, '.tif')))
      stopifnot(compareCRS(stand_age,canopy_cover))
      stopifnot(compareCRS(stand_age, roads))
      stopifnot(all(rownames(trap_coordinates) == rownames(capture_count_at_traps)))
      stopifnot(all(rownames(trap_coordinates) == colnames(colony_count_at_traps)))
      stopifnot(all(rownames(trap_coordinates) == rownames(floral_cover_at_traps)))
      trap_coordinates <- SpatialPoints(trap_coordinates)
      crs(trap_coordinates) <- crs(stand_age)
      save(
        floral_cover_at_traps,
        stand_age,
        canopy_cover,
        roads,
        trap_coordinates,
        colony_count_at_traps,
        capture_count_at_traps,
        file=paste0("../data/", region, ".", species, ".", resolution, ".RData")
      )
    }
  }
}
