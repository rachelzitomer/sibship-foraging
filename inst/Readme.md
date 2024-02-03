This folder contains data and scripts associated with the manuscript: "Canopy cover has divergent effects on movement of closely related bumble bees in managed conifer forest landscapes" prepared to submit to the journal Landscape Ecology. 

In the example dataset, workers of two bumble bee species (Bombus vosnesenskii and Bombus caliginosus) were sampled and genotyped at 75 sites (1-25) in 3 replicate intensively managed forest landscapes ("south", "midnorth", and "north") in western Oregon to examine effects of forest cover on space use and site-level colony abundance.

raw_data file descriptions:

*colony_abundance.csv includes standardized colony abundance estimates (Nns) along with intermediary calculations for each species x site
*genotyped_captures_species_landscape.csv (6 files)  is a matrix of workers of the specified species genotyped at each site in the specified landscape. Rows represent colonies and columns (named) represent sites.
*total_captures_species_landscape.csv (6 files) is a table containing counts of individual bees of the specified species captured at each site in the specified landscape. Counts are given for each sex separately (n_M/n_F) and together (total).
*floral_cover_landscape.txt (3 files) is a vector of floral density values for each site in the specified landscape
*landscape_latlongs.csv (3 files) is a table of lat/long coordinates (WSG84) for each site in the specified landscape.
*landscape_trap_coords.txt (3 files) is a table  of coordinates (Albers Equal Area Conic Projection) for each site in the specified landscape.
*ocr_age_canopycover.csv is a table of stand age and average canopy cover collected in 60 intensively managed Douglas-fir stands in the central part of the Oregon Coast Range during 3 sampling rounds in the late spring-summer of 2018 and 2019 (data from Zitomer et al. 2023)
*landscape_agemap_##m.tif (18 files) is a raster of forest age (years since stand-replacing disturbance in summer 2021) estimated using the disturbance detection algorithm LandTrendr for the specified landscape. Resolution is specified by ##m.
*landscape_canopy_cover_##m.tif (18 files) is a raster of canopy cover projected from forest age rasters based on stand age and canopy cover data collected in a previous study (see ocr_age_canopycover.csv)
*landscape_roads_binary_##m.tif (18 files) is a binary raster of roads in the specified landscape; resolution is specified by ##m.
*package data.R packages raw_data into files in th data folder

data file descriptions:

*landscape.species.##m.RData (18 files) is packaged data from raw_data ready for import in foraging model analyses; resolution is specified by ##m
*structural_connectivity.RData is packaged data ready for import in structural connectivity analyses

analysis/scripts file descriptions:

*fit_stand_age and fit_stand_age_and_roads establish foraging model fitting and bootstrapping with command line parsing for entering model specs (study landscape, resolution, species, lower/upper bounds and grid search size for theta parameters, # bootstrap simulations and block size) used in fitting for a foraging model containing 1 (stand age) or 2 (stand age and road cover) theta parameters of interest

*GEE_script.txt is code used in Google Earth Engine script to obtain stand age and NCLD rasters

*structural_connectivity.R is code used in structural analysis, including function to create modified Hanski connectivity index

analysis/resolution file descriptions:

*tmp.R is a script demonstrating how to construct the figures in Rplots.pdf, fig/south_constrained_sim_at_null.loglik.png and fig/south_constrained.sim_at_null.landscape_distance_on_capture_rate.png from data/south.bomvos.90m.RData
*tmp2.R is a script demonstrating model fitting and bootstrapping for bomvos.midnorth.90m.RData, produces data objects: bomvos.midnorth.90m.stand_age_and_roads.fitted.RData and bomvos.midnorth.90m.stand_age_and_roads.bootstrap.RData

analysis/stand_age file descriptions:

*landscape.stand_age.R demonstrate resistance and foraging model fitting for effects of stand age on foraging range of B. vosnesenskii in the each study landscape using 90 m resolution data and produces corresponding figures in analysis/stand_age/fig
*midnorth.roads.R demonstrates resistance model and foraging model fitting for effects of road cover on foraging range of B. vosnesenskii in the north study landscape using 90 m resolution data

analysis/singletons file descriptions:
*midnorth.singletons.R demonstrates foraging model fitting with only singleton colonies

