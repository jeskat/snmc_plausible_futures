### This script filters the FIA database to select a set of plots representative
### of Sierra Nevada mixed-conifer forests.


# Dependencies and data ---------------------------------------------------
library(sf)
library(rFIA)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(prospectr)
library(corrplot)
library(terra)
library(tidyterra)

## Use existing PRISM datafile?
prism_exists <- TRUE


## Load Sierra Nevada ecoregion bounds (this data needs to to be downloaded)
sn_bounds <- st_read(here::here('Data/Sierra_Nevada Conservancy_subregions/Sierra_Nevada_Conservancy_Subregions.shp')) 

## Load FIA data for CA (this data needs to be downloaded)
ca <- readFIA(here::here('Data/FIA/'))

## Load PRISM normals for FIA plots (needs to be downloaded)
if(prism_exists){
  prism <- read.csv('Data/prism_normals_1991_2020.csv')
}



# Subset SNMC FIA plots ---------------------------------------------------

## Limit to Sierra Nevada ecoregion
sn <- clipFIA(ca, mostRecent = TRUE, mask = sn_bounds, nCores = 3)

## Identify plots with CA mixed conifer forest types (FORTYPCDCALC = 371)
snmc_structure <- tpa(sn, 
                      areaDomain = FORTYPCDCALC == 371 & PLOT_STATUS_CD == 1, 
                      byPlot = TRUE)

## Subset plots that are entirely forest
snmc_structure <- snmc_structure[snmc_structure$PROP_FOREST==1,]



# Export plot coordinates for PRISM (optional) ---------------------------------
## This only needs to be done once, to get PRISM data. 

if(prism_exists==FALSE){
  ### Get lat/long for PRISM data
  plot_coords <- snmc_structure %>% left_join(sn$PLOT[, c('PLT_CN', 'pltID', 'LAT', 'LON')],
                                        by = c('pltID', 'PLT_CN'))
  plot_coords <- plot_coords[, c('LAT', 'LON', 'pltID')]
  write.csv(plot_coords, file = here::here('Outputs/plots_coordinates_for_prism.csv'), 
            row.names = FALSE,
            col.names = NA)
  }


# Characterize forest structure, composition, and climate -----------------

## Add stand density index (summation method)
calculate_sdi <- function(PLT_CN){
  trees_in_plot <- sn$TREE[sn$TREE$PLT_CN==PLT_CN,]
  live_trees <- trees_in_plot[trees_in_plot$STATUSCD==1,]
  sdi <- sum(live_trees$TPA_UNADJ* (live_trees$DIA/10)**1.6)
  return(sdi)
}

snmc_structure['SDI'] <- sapply(snmc_structure$PLT_CN, calculate_sdi)


## Merge in PRISM normals (1991-2020) 
snmc_structure <- snmc_structure %>% 
  inner_join(prism, by = 'pltID')


## Get basal area and stand density by species for CA mixed conifer forest plots
snmc_spp <- tpa(sn,
                bySpecies = TRUE,
                areaDomain = FORTYPCDCALC == 371 & PLOT_STATUS_CD == 1, 
                byPlot = TRUE)

snmc_spp <- snmc_spp[snmc_spp$PROP_FOREST==1,]

## Add each species' fraction of BA
snmc_spp <- snmc_spp %>%
  group_by(PLT_CN, pltID) %>% #, CONDID
  mutate(frac_ba = BAA / sum(BAA))

## Assign a PFT to each species
species_to_pft = list(
  `white fir` = 'fir',
  `incense-cedar` = 'cedar',
  `sugar pine` = 'pine',
  `ponderosa pine` = 'pine',
  `Douglas-fir` = 'psme',
  `canyon live oak` = 'oak',
  `California red fir` = 'fir',
  `western white pine` = 'other',
  `California black oak` = 'oak',
  `Jeffrey pine` = 'pine',
  `western juniper` = 'other',
  `Pacific dogwood` = 'other',
  `lodgepole pine` = 'other',
  `bigleaf maple` = 'other',
  `Shasta red fir` = 'fir',
  `Oregon white oak` = 'oak',
  `interior live oak` = 'oak'
)

snmc_spp$PFT <- 'other'

for(i in 1:nrow(snmc_spp)){
  if(snmc_spp$COMMON_NAME[i] %in% names(species_to_pft)){
    snmc_spp$PFT[i] = species_to_pft[[snmc_spp$COMMON_NAME[i]]]
  }
}

## Calculate fraction PFT per plot
snmc_pft <- snmc_spp %>% 
  group_by(PLT_CN, pltID, PFT) %>% #CONDID, 
  summarise(frac_pft = sum(frac_ba))

## Convert long dataframe to wide
pft_wide <- snmc_pft %>%
  pivot_wider(names_from = PFT, values_from = frac_pft) %>%
  replace_na(list(cedar = 0, fir = 0, pine = 0, oak = 0, other = 0, psme = 0))


## PFT filters
# Remove plots where >10% BA is "other" PFT (i.e., a PFT we don't model in FATES)
to_filter <- snmc_pft[snmc_pft$PFT == 'other' & snmc_pft$frac_pft>0.1,]
snmc_structure <- snmc_structure[!(snmc_structure$PLT_CN %in% unique(to_filter$PLT_CN)),]

# Remove plots where oak fraction is > 0.5
oak_filter <- snmc_pft[snmc_pft$PFT == 'oak' & snmc_pft$frac_pft>0.5,]
snmc_structure <- snmc_structure[!(snmc_structure$PLT_CN %in% unique(oak_filter$PLT_CN)),]

# Remove plots where Douglas fir fraction is > 0.5
psme_filter <- snmc_pft[snmc_pft$PFT == 'psme' & snmc_pft$frac_pft > 0.5,]
snmc_structure <- snmc_structure[!(snmc_structure$PLT_CN %in% unique(psme_filter$PLT_CN)),]

## "Minimal plots:" Filter plots that are privately owned and/or >10% PSME + other PFT

# Join ownership (and some other useful) information, which are provided at the condition level
plot_conditions <- snmc_structure %>% left_join(sn$COND[, c('PLT_CN', 
                                                            'CONDID', 
                                                            'PHYSCLCD', 
                                                            'SITECLCD', 
                                                            'OWNGRPCD')],
                                                by = c('PLT_CN'))

# Identify privately-owned plots
private <- plot_conditions[plot_conditions$OWNGRPCD==40,]
snmc_structure$Ownership <- 'Not private'
snmc_structure[snmc_structure$PLT_CN %in% private$PLT_CN, 'Ownership'] = "Private"

# Identify plots where fraction Douglas fir is between 0.1 and 0.5 
psme_plots <- pft_wide[pft_wide$psme + pft_wide$other > 0.1 & pft_wide$other <=0.1 & pft_wide$psme <= 0.5,]
snmc_structure$PSME_check = 'Not_psme'
snmc_structure[snmc_structure$PLT_CN %in% psme_plots$PLT_CN, 'PSME_check'] = 'PSME'

# Flag plots that are excluded either for ownership or PSME reasons
snmc_structure$Include <- "No"
snmc_structure[snmc_structure$Ownership == 'Not private' & snmc_structure$PSME_check == 'Not_psme', "Include"] = "Yes"

# Create dataframe of plots included under restrictive criteria
snmc_struct_min <- snmc_structure[snmc_structure$Include=='Yes',]


# PFT composition ordination analysis -------------------------------------

## Reduce PFT dataframe to include the minimal set of plots
pft_wide_min <- pft_wide[pft_wide$pltID %in% unique(snmc_struct_min$pltID),]

## Restructure dataframe for ordination
comp_data <- pft_wide_min[, -1:-2]
rownames(comp_data) <- pft_wide_min$pltID

## NMDS
# k=2 means we want a 2-dimensional solution (a 2D "map")
# Run it directly on the *original* fractional data
nmds_result <- metaMDS(comp_data, k = 2, trymax = 100)

# Check the 'stress' value.
# A stress < 0.2 is considered good for ecological data.
# A stress < 0.1 is great.
print(nmds_result)

# Get the individual plot scores from the NMDS result
plot_scores_nmds <- scores(nmds_result, display = "sites")

# Extract just the NMDS1 scores as the composition variable
nmds1_scores <- plot_scores_nmds[, 1]

# Add this to the original compositional dataframe
pft_wide_min$NMDS1 <- nmds1_scores

## Create NMDS plot
# 'type = "n"' creates an empty plot. We'll add things to it.
plot(nmds_result, type = "n") #main = "NMDS of FIA Plots"

# Add the plots as points
points(nmds_result, display = "sites", cex = 0.8, col = "gray")


# Fit the PFTs as vectors
# This calculates the correlation of each species' percentage 
# with the NMDS1 and NMDS2 axes.
# Use our *original* composition data: 'comp_data'
fit <- envfit(nmds_result, comp_data, na.rm = TRUE)

# Add the fitted vectors to the plot
# This draws the arrows for the PFT
plot(fit, col = "navy", cex = 1)

## Option to save grapic
ggsave('Outputs/Fig_oridination_plot.png', device = 'png',
      width = 4, height = 4, units = 'in')

# See the R-squared and p-values for each PFT (how well it fits the ordination).
print(fit)

## Add NDMS1 the master dataframe
snmc_struct_min <- inner_join(snmc_struct_min, 
                              pft_wide_min[, c('pltID', 'NMDS1')], 
                              by = 'pltID')

### Save selected plots to disk
write.csv(snmc_struct_min, file = here::here('Outputs/candidate_FIA_plots.csv'), 
          row.names = FALSE)



# Kennard-Stone maximin distance algorithm to sample candidate plots --------

# ## Option to start here and read in the plots from disk
# snmc_struct_min <- read.csv(here::here('Outputs/candidate_FIA_plots.csv'))

## Define which variables over which to sample
covs <- c('ppt_mm', 'tmean_degC', 'SDI', 'NMDS1')

## Center and scale variables
scaled_vars <- scale(snmc_struct_min[, covs])

## Select points that are furthest from each other in the 4D space
selection <- kenStone(X = scaled_vars, k = 81)

## Extract the sampled data
sampled_data <- snmc_struct_min[selection$model, ]


## Save the samples plots to disk
write.csv(sampled_data, file = 'Outputs/sampled_fia_plots_81_kenStone.csv', 
          row.names = FALSE)



# Visualizations ----------------------------------------------------------

## Correlation plot ----------------------------------------------------------
corrplot(cor(snmc_struct_min[, covs]), method = 'number', type = 'upper')

## Map variable labels to column names 
corrplot_shorthand <- list(
  'ppt_mm' = 'Precipitation',
  'tmean_degC' = 'Temperature',
  'SDI' = 'SDI',
  'NMDS1' = 'Composition'
)

temp <- sampled_data[, covs]

for(c in colnames(temp)){
  if(c %in% names(corrplot_shorthand)){
    colnames(temp)[colnames(temp)==c] = corrplot_shorthand[c]
  }
}

## Save correlation plot
png(here::here('Outputs/Fig_correlation_plot.png'))
corrplot(cor(temp), 
         method = 'color', type='upper', addCoef.col = 'black', 
         diag = FALSE)
dev.off()


## Plots in climate space ------------------------------------------------------

# Visualize 4D space for candidate and sampled plots
clim_space_p <- ggplot(data = snmc_struct_min, 
                       aes(x = tmean_degC, 
                           y = ppt_mm, 
                           size = SDI, 
                           color = NMDS1)) + 
  geom_point(alpha = 0.75) +
  scale_color_viridis_c() +
  geom_point(data = sampled_data, 
             shape = 1,
             color = 'black',
             stroke = 1.2,
             show.legend = FALSE) +
  xlab('Mean temperature (degrees Celsius)') +
  ylab('Mean annual precipitation (mm)') +
  labs(size= 'SDI \n(trees/acre)', color = "Composition \n(NMDS1)") +
  theme_minimal() +
  theme(axis.text=element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16))

ggsave(here::here('Outputs/Fig_sampled_plots_clim_space.png'),
       plot = clim_space_p, 
       device = 'png',
       width = 6, height = 4, units = 'in')


## Map of plots ------------------------------------------------------------
## Convert lat lon columns to points
sf_points <- st_as_sf(snmc_struct_min, coords = c('LON', 'LAT'), crs = 4326, 
                      remove = FALSE)
test_points <- st_as_sf(sampled_data, coords = c('LON', 'LAT'), crs = 4326, 
                        remove = FALSE)

## Import elevation raster
elev <- rast(here::here('Data/SRTM_clippedCA.tif'))

## Reproject SN bounds
sn_outline <- st_transform(sn_bounds, st_crs(test_points))

# Dissolve interior polygons
sn_outline <- st_union(st_buffer(sn_outline, 0.0001))

## Make map
map_p <- ggplot() +
  geom_spatraster(data = elev, alpha=0.5) +
  scale_fill_wiki_c() +
  geom_sf(data = sn_outline, fill = NA) + 
  geom_sf(data = sf_points, color = 'gray', alpha = 0.75) +
  geom_sf(data = test_points, aes(color = NMDS1), alpha = 0.75) + 
  scale_color_viridis_c() +
  labs(fill= 'Elevation (m)', color = "Composition") +
  scale_x_continuous(breaks = c(-120, -115)) +
  scale_y_continuous(breaks = c(35, 40)) +
  theme(panel.background = element_rect(fill = NA, colour = 'NA')) +
  theme(axis.line = element_line(color='gray'),
        axis.text=element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
map_p


ggsave(here::here('Outputs/Fig_sampled_plots_map.png'), 
       plot = map_p, 
       device = 'png',
       width = 4, height = 4, units = 'in')
