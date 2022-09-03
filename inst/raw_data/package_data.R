library(raster)
library(sp)

dir.create("../data")
for (species in c("bomvos", "bomcal"))
{
  for (region in c("north", "south", "midnorth"))
  {
    floral_cover_at_traps <- as.matrix(read.table(paste0(region, ".floral_cover.txt"), header=TRUE))
    colony_count_at_traps <- as.matrix(read.table(paste0(region, ".colony_count_at_traps.", species, ".txt"), header=TRUE))
    capture_count_at_traps <- as.matrix(read.table(paste0(region, ".capture_count_at_traps.", species, ".txt"), header=TRUE))
    for (resolution in c("30m", "60m", "90m"))
    {
      stand_age <- readAll(raster(paste0(region, ".stand_age.", resolution, ".tif")))
      roads <- readAll(raster(paste0(region, ".roads.", resolution, ".tif")))
      trap_coordinates <- read.table(paste0(region, ".trap_coords.", resolution, ".txt"), header=TRUE)
      stopifnot(all(rownames(trap_coordinates) == rownames(capture_count_at_traps)))
      stopifnot(all(rownames(trap_coordinates) == colnames(colony_count_at_traps)))
      stopifnot(all(rownames(trap_coordinates) == rownames(floral_cover_at_traps)))
      stopifnot(compareCRS(stand_age, roads))
      trap_coordinates <- SpatialPoints(trap_coordinates)
      crs(trap_coordinates) <- crs(stand_age)
      save(
        floral_cover_at_traps,
        stand_age,
        roads,
        trap_coordinates,
        colony_count_at_traps,
        capture_count_at_traps,
        file=paste0("../data/", region, ".", species, ".", resolution, ".RData")
      )
    }
  }
}
