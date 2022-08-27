library(raster)
dir.create("../../data")
for (species in c("bomvos", "bomcal"))
{
  for (region in c("north", "south", "midnorth"))
  {
    floral_cover_at_traps <- read.table(...)
    colony_counts_at_traps <- read.table(...)
    for (resolution c("30m", "60m", "90m"))
    {
      stand_age <- readAll(raster(...))
      roads <- readAll(raster(...))
      trap_coordinates <- read.table(...)
      save(
        floral_cover_at_traps,
        stand_age,
        roads,
        trap_coordinates,
        colony_counts_at_traps,
        file=paste0("../../data/", species, "_", region, "_", resolution, ".RData")
      )
    }
  }
}
