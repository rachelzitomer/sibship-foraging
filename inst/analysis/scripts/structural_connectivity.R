#Script purpose:
#calculate modified Hanski connectivity index
#regress colony abundance ~ connectivity index*species + connectivity index*flower density + (1|landscape)

#I need to choose a focal landscape size based on max estimated daily movement distance for our species
  #Pope and Jha 2018 models make similar assumptions to ours and uses B. vos as the focal species, so seems most suitable:
    #obs. min foraging distance for colonies ranged from 25 to 1500 m
    #estimated foraging distance for individuals ranged from 118 to 3175 m 

max_foraging_distance = 3175

#alpha, on the other hand, should represent (the inverse of) average home range
  #from Pope and Jha 2018;
    #obs. mean +/- min foraging distance for colonies: 357.4 +/- 51.6 m)
    #estimated foraging distance for individuals: 453.7 +/- 15.0 m)

alpha = 1/(453.7/1000)#divide by 1000 to get units in km

#Hanski's connectivity index formula for focal patch i:

#CI[i] = sum for all values i and j: exp(-alpha*d[ij])*A[j]^b

#alpha = 1/average migration distance in km
#d[ij] = distance to patch j from focal patch i in km
#A[j] = area of patch j in m^2
#b = parameter which scales the size of surrounding habitat patches
#moilanen & nieminen suggested that ratio of patch edge to patch size decreases with A^0.5 as patch size increases

b = 0.5

#ea landscape has 25 sites, and for ea patch we want to make a matrix or something with a row for ea neighboring patch

#for (i in 1:length(sites)){
  # for (j in 1:length(neighboring_patches)){
    #CI[j] = exp(-alpha*d[ij])*A[j]^b
  #}
  #sum(CI)
#}

#Some versions of this equation include area of the focal patch (a multiplier before the summation) 
#but we won't always have a true focal patch (e.g. in roadside sampling sites)
#I can use a spatial point (the trap) as the "focal patch" 
#then for the traps that are in patches account for patch area by measuring the distance (0) between the point and the patch it is in. 
#This is, I think, the same as the original IFM equation, although it's not clear to me whether or not the focal patch was included in the measure there. 

#A basic example to make sure the measure changes as expected as connectivity increases/decreases:

example_data<-as.matrix(expand.grid(
  'distance' = seq(0.2,4,length.out = 16),
  'area' = seq(1000, 100000, length.out = 16)
  ))

IFM_example_function<-function(distance, area, alpha, b){
  example_data_vector<-c()
  for(i in 1:length(distance)){
    example_data_vector[i]<-exp(-alpha*distance[i])*area[i]^b
  }
  return(sum(example_data_vector))
}

#same distances, different area:
IFM_example_function(example_data[1:16, 1], example_data[1:16,2], alpha, b)
IFM_example_function(example_data[17:32, 1], example_data[17:32,2], alpha, b)
#IFM increases with area of neighbor patches

#same areas, different distances:
example_data2<-example_data[order(example_data[,1],decreasing = FALSE),]
IFM_example_function(example_data2[1:16,1], example_data2[1:16,2], alpha, b)
IFM_example_function(example_data2[17:32,1],example_data2[17:32,2], alpha, b)
#IFM decreases with distance to neighbor patches

#plot IFM shape across distance range (0 to beyond the buffer range)

plot_ex<-as.matrix(
  cbind('distance' = rep(0,nrow(example_data)), 'IFM' = rep(0,nrow(example_data))),byrow = TRUE)

for(j in 1:nrow(example_data[1:16,])){
  IFM<-IFM_example_function(example_data[j,1], example_data[j,2], alpha, b)
  plot_ex[j,1]<-round(example_data[j,1],2)
  plot_ex[j,2]<-round(IFM,2)
}

plot(plot_ex[,1], plot_ex[,2], xlab = 'distance', ylab = 'IFM')

#IFM is a negative exponential decay function which approaches zero as the distance approaches infinity 
#the slope of the function is determined by alpha and IFM is v small when distance > 4* 1/alpha
#this represents movement (foraging or migration) decreasing in incidence as distance to the neighboring patches increases relative to the movement capacity of the species.
#this does mean that the model is pretty insensitive to the focal landscape distance (estimated max foraging distance)

#steps: 
##for each trap make a raster that is clipped to focal radius
##get distance between focal patch and all other patches within a 3175 m radius and area of all other patches

#load packages-----
#spatial packages
library(raster)
library(sf)
library(rgeos)
#general tidying and visualization packages
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(patchwork)
library(GGally)
library(smoothr)
library(viridis)
#linear modeling packages
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(emmeans)
#autocorrelation packages
library(gstat)
library(ncf)
#data package
library(sibships)

##import data----
load(system.file("data/structural_connectivity.RData", package="sibships"))

#includes:
#colony_abundance.csv (colony_abundance)
#south_early_seral_binary_map.tif (south_binary_map)
#midnorth_early_seral_binary_map.tif (midnorth_binary_map)
#north_early_seral_binary_map.tif (north_binary_map)
#south_latlongs.csv (south_ll)
#midnorth_latlongs.csv (midnorth_ll)
#north_latlongs.csv (north_ll)

#quick summary of sample sizes----
landscape_summary<-colony_abundance %>% 
  group_by(Species, Landscape) %>% 
  summarise(n_total = sum(N),
            n_workers = sum(Ni),
            n_genotyped = sum(Ng))

##check raster res----
res(south_binary_map)==res(midnorth_binary_map)
res(south_binary_map)==res(north_binary_map)

##get crs----
albers<-crs(south_binary_map)

##transform latlong data to sf objects----
south_traps<-south_ll %>% 
  st_as_sf(.,coords=c("Long","Lat"),crs=4326)

south_traps<-st_transform(south_traps,albers)

plot(south_binary_map); plot(south_traps$geometry, pch = 16, col = 'red', add =TRUE)

midnorth_traps<-midnorth_ll %>% 
  st_as_sf(.,coords=c("Long","Lat"),crs=4326)

midnorth_traps<-st_transform(midnorth_traps,albers)

plot(midnorth_binary_map); plot(midnorth_traps$geometry, pch = 16, col = 'red', add =TRUE)

north_traps<-north_ll %>% 
  st_as_sf(.,coords=c("Long","Lat"),crs=4326)

north_traps<-st_transform(north_traps,albers)

plot(north_binary_map); plot(north_traps$geometry, pch = 16, col = 'red', add =TRUE)

##make a function to calculate Hanski connectivity index----

Hanski_ConnIndex<-function(binary_raster, traps_sf, max_foraging_distance, alpha, b){
  traps_sf$conn_index<-NA
  poly<-rasterToPolygons(binary_raster,fun = function(x){x==1}, dissolve = TRUE)
  poly<-as(st_as_sf(poly),'Spatial')
  polys<-raster::disaggregate(poly)
  polys_sf<-st_as_sf(polys)
  #eliminate tiny polygons and gaps:
  area_thresh<-units::set_units(90, m^2)
  polys_filled<-smoothr::fill_holes(polys_sf, area_thresh)
  polys_dropped<-drop_crumbs(polys_filled, area_thresh, drop_empty = TRUE)
  plot(binary_raster); plot(polys_dropped$geometry, add =TRUE);
  polys_dropped$area<-st_area(polys_dropped)#calculate area in m^2 fo each polygon in ls
  for(i in 1:nrow(traps_sf)){
    buffer<-st_buffer(traps_sf[1,], max_foraging_distance)
    site_polys<-polys_dropped[buffer,]
    site_polys_dist<-as.numeric(st_distance(traps_sf[i,],site_polys))/1000#divide by 1000 to get dists in km
    site_CI_list<-c()
    for(j in 1:length(site_polys_dist)){
      site_CI_list[j]<-exp(-alpha*site_polys_dist[j])*as.numeric(site_polys[j,'area'])[1]^b
    }
    traps_sf[i,'conn_index']<-sum(site_CI_list)
  }
  return(traps_sf)
}

##calculate hanski index for each site:----
south_traps_ci<-Hanski_ConnIndex(south_binary_map, south_traps, max_foraging_distance =  max_foraging_distance, alpha =  alpha, b = b)

midnorth_traps_ci<-Hanski_ConnIndex(midnorth_binary_map, midnorth_traps, max_foraging_distance = max_foraging_distance, alpha = alpha, b = b)

north_traps_ci<-Hanski_ConnIndex(north_binary_map, north_traps, max_foraging_distance = max_foraging_distance, alpha = alpha, b = b)

#put all the indices into a single dataframe:---

All_traps_ci<-do.call('rbind', list(south_traps_ci, midnorth_traps_ci, north_traps_ci))

#get x y coords (apart from geom):----
All_traps_ci$x<-st_coordinates(All_traps_ci)[,1]
All_traps_ci$y<-st_coordinates(All_traps_ci)[,2]

#join with colony abundance data----

names(colony_abundance)[1]<-'Site_Name'

All_traps_ci_col<-merge(All_traps_ci[,c('Site_Name','Landscape','Trap_num','conn_index','x','y')], 
                        colony_abundance[,c('Site_Name','Species', 'Nns', 'flower_density')], by = 'Site_Name')

All_traps_ci_col$Species = factor(All_traps_ci_col$Species,
                              levels = c("caliginosus","vosnesenskii"),
                              labels = c("B. caliginosus","B. vosnesenskii"))

All_traps_ci_col<-All_traps_ci_col %>% 
  mutate(Landscape = as.factor(Landscape)) %>% 
  mutate(Landscape = fct_relevel(Landscape, 'South','Midnorth', 'North'))

#summarise and examine variables----

##colony abundance:----

Species_nns_tbl<-summarise_at( group_by(as.data.frame(All_traps_ci_col),Landscape, Species), "Nns",
                              list(min = min, 
                                   median = median, 
                                   mean = mean, 
                                   max = max,
                                   sd = sd) ) %>% arrange(Species) %>% 
  ungroup()


#central tendency statistics for vos are much smaller in South landscape than the other two;
#also smallest in the South for cal, but less dramatic difference
#variance is greatest in the North for both species, but the difference is more extreme for vos. 

##connectivity index----
conn_index_tbl<-summarise_at( group_by(as.data.frame(All_traps_ci_col),Landscape), "conn_index",
                              list(min = min,
                                   median = median,
                                   mean = mean,
                                   max = max,
                                   sd = sd) )

##plot distribution of colony abundance:--------
qplot(x=Nns,facets=~Species,data=All_traps_ci_col,geom="histogram")#similar distribution, but longer tail for vos.

qplot(x=Species,y=Nns,data=All_traps_ci_col,geom="boxplot")#variances looking pretty durn different between species

qplot(x=Nns,facets=~Species+Landscape,data=All_traps_ci_col,geom="histogram")
#similar patterns in each landscape with both species, but there are just way fewer cal colonies

qplot(x=Species,y=Nns,facets = ~Landscape, data=All_traps_ci_col,geom="boxplot")
#differences between species are pretty consistent among landscapes

##plot distribution of conn_index:--------

par(mfrow=c(2,1))

qplot(x=Landscape,y=conn_index,data=All_traps_ci_col,geom="boxplot")
#the ranges of values differs among landscapes but the variance appears pretty consistent
#midnorth has the highest overall connectivity, north has the lowest

qplot(x=conn_index,facets = ~Landscape,data=All_traps_ci_col,geom="histogram")
#tracks pretty well with what we would expect from the boxplot

ggplot(All_traps_ci_col, aes(x = conn_index, y = Nns, group = Species)) +
  geom_point(aes(color = Species)) +
  geom_smooth(method = 'lm', formula = y~x, aes(color = Species), se = FALSE) +
  facet_grid(.~Landscape) +
  theme_classic()

ggplot(All_traps_ci_col, aes(x = conn_index, y = Nns, group = Species)) +
  geom_boxplot(aes(fill = Species)) +
  facet_grid(.~Landscape)

#pattern is pretty consistent among landscapes wrt relationship between species
#however, difference in variance among landscapes is pretty large for vos Nns

#write a dataframe without sf for models:----

All_traps_ci_col_df<-All_traps_ci_col %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry)

#tweedie distributed model:----
##seems like a suitable distribution because we have right-skewed, non-negative, continuous data 

mod1_data<-All_traps_ci_col_df %>% 
  mutate(conn_index = scale(conn_index))

mod1_tweedie<-glmmTMB(Nns~conn_index + flower_density + Species + 
                        conn_index:Species + conn_index:flower_density + 
                        (1|Landscape),
                      data = mod1_data,
                      family = tweedie)

mod1_sims =simulateResiduals(mod1_tweedie)
mod1_data$res<-mod1_sims$scaledResiduals

#residual diagnostics
plot(mod1_sims)
#qqplot seems reasonable, quantile deviation is pretty minor
testQuantiles(mod1_sims)
qplot(x = "res", y = res, data = mod1_data,
      geom = "boxplot",
      main = "Boxplot of standardized residuals")#looks quite symmetrical

#residual dispersion, zero-inflation, outliers
par(mfrow=c(1,3))
testDispersion(mod1_sims)#dispersion = 1.35, this is fine
testZeroInflation(mod1_sims)#looks fine
testOutliers(mod1_sims)

#plot residuals vs. independent variables
plotResiduals(mod1_sims, mod1_data$Species)#non-homogeneneous, but this is okay for Tweedie
plotResiduals(mod1_sims, mod1_data$conn_index)#residuals vs. conn index looks fine
plotResiduals(mod1_sims, mod1_data$Landscape)#non-homogeneneous, but this is okay for Tweedie

#autcorrelation:

site_xy<-aggregate(mod1_data[,5:6],list(mod1_data$Landscape, mod1_data$Site_Name),mean)

mod1_sims_grouped_south = recalculateResiduals(mod1_sims, group = mod1_data$Site_Name, 
                                               sel = mod1_data$Landscape == 'South')

testSpatialAutocorrelation(mod1_sims_grouped_south, 
                           site_xy[site_xy$Group.1=='South',3], 
                           site_xy[site_xy$Group.1=='South',4])#looks fine

#midnorth:
mod1_sims_grouped_midnorth = recalculateResiduals(mod1_sims, group = mod1_data$Site_Name, 
                                               sel = mod1_data$Landscape == 'Midnorth')

testSpatialAutocorrelation(mod1_sims_grouped_midnorth, 
                           site_xy[site_xy$Group.1=='Midnorth',3], 
                           site_xy[site_xy$Group.1=='Midnorth',4])#looks fine

#north:
mod1_sims_grouped_north = recalculateResiduals(mod1_sims, group = mod1_data$Site_Name, 
                                               sel = mod1_data$Landscape == 'North')

testSpatialAutocorrelation(mod1_sims_grouped_north, 
                           site_xy[site_xy$Group.1=='North',3], 
                           site_xy[site_xy$Group.1=='North',4])#looks fine

#type II anova:---
car::Anova(mod1_tweedie,type = 2)

emmeans_species = emmeans(mod1_tweedie, specs = pairwise ~ Species, 
                          type = 'response', ref = 1)
emmeans_species$emmeans
emmeans_species$contrasts

#confidence interval slope of conn_index on colony abundance for each species:----
confint(mod1_tweedie)

slope_confint_by_species = emtrends(mod1_tweedie, ~Species, var = 'conn_index', type = 'response')


slope_confint_by_species = as.data.frame(slope_confint_by_species)

names(slope_confint_by_species)

#plotting:-----

####plot theme:
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 10),
  legend.key = element_blank())

slope_confint_plot<-ggplot(slope_confint_by_species, aes(x = Species, y = conn_index.trend, 
                                                         color = Species)) +
  geom_errorbar(width = 0.1, lwd = 0.8, 
                 aes(ymin = asymp.LCL, ymax = asymp.UCL)) +
  scale_color_viridis(discrete = TRUE, begin = 0, end = 0.5, name ='',option='D') +
  geom_point(size = 1.5) +
  labs(y = 'Slope estimate') +
  geom_hline(yintercept = 0, lty = 2) +
  ylim(-1.5,1.5) +
  # scale_x_continuous(breaks = seq(-2,2, 0.5)) +
  theme_bw(base_size = 15) +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10, face = 'italic'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x.bottom = element_line(color="black", size =0.5),
        # axis.line.y.left=element_blank(),
        axis.title.y.left = element_text(size = 10),
        axis.title.x.bottom = element_text(size=10)) +
  # geom_text(aes(x = conn_index.trend-0.25, label = Species),
  #           nudge_y = 0.2,
  #           size = 4, fontface = 'italic', hjust = 0) +
  scale_x_discrete(limits = rev(slope_confint_by_species$Species))+
  scale_y_continuous(limits = c(-1.0, 1.0), 
                     breaks =c(-1,-0.5,0,0.5, 1), 
                     expand = c(0,0))+
  BioR.theme +
  theme(legend.position = 'none')
  
slope_confint_plot

#linear xy plot:----
range(mod1_data$conn_index)
range(mod1_data$Nns)

slope_by_species = ggeffect(mod1_tweedie, terms = c('conn_index', 'Species'))

slope_by_species

slope_by_species_interaction_plot<-ggplot(slope_by_species, 
                                          aes(x = x, y = predicted, colour = group,lty = group))+
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax =conf.high, fill = group), color = NA, alpha = 0.15)+
  scale_x_continuous(limits = c(-2.5,2.5), expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_color_viridis(discrete = TRUE, begin = 0, end = 0.5 ,labels = c('Bombus caliginosus    ', 'Bombus vosnesenksii'), name ='',option='D') +
  scale_fill_viridis(discrete = TRUE, begin = 0, end = 0.5,labels = c('Bombus caliginosus    ', 'Bombus vosnesenksii'), name = '',option='D') +
  scale_linetype_manual(values = c('solid','dotted'), labels = c('Bombus caliginosus    ', 'Bombus vosnesenksii'), name = '') +
  labs(x = "Standardized early seral forest connectivity" ,
       y = "Mean colony abundance") + # change axis labels
  BioR.theme +
  theme(axis.title.x = element_text(size=10, color="black"),
        axis.title.y = element_text(size = 10, color="black"),
        axis.text = element_text(size = 8, color = 'black'),
        plot.title=element_blank(),
        # plot.margin=unit(c(0,0.25,0,0),"cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x.bottom = element_line(color="black"),
        axis.line.y.left=element_line(color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(face = 'italic'),
        legend.spacing.x = unit(0.2,'cm'))

slope_by_species_interaction_plot

legend<-get_legend(slope_by_species_interaction_plot)

slope_by_species_interaction_plot<-slope_by_species_interaction_plot+
  theme(legend.position = 'none')

plot_content<-(slope_by_species_interaction_plot + slope_confint_plot) +
  plot_layout(widths = c(7,5)) &
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 10, face = 'bold'))

plot_content_legend<-wrap_elements(plot_content)/legend+
  plot_layout(heights = c(20,1))

plot_content_legend