### Clustering analysis for island wide nutrients ####3
### Script by Nyssa Silbiger #####
### Date edited 2021 - 01 - 29 ####
################################

# Load libraries #####
library(tidyverse)
library(here)
library(ggdendro)
library(patchwork)
library(ggmap)
library(viridis)
library(maptools)
library(ggforce)

# Read Data ######
# Water column Data
MayData<-read_csv(here("Data","Lagoon_nutrients_May_2021.csv"))
AugustData<-read_csv(here("Data","Lagoon_nutrients_August_2020.csv"))

## rename the column names to have august and may
MayData<-MayData%>% 
  select(-Entry, Site_Number = Site_number) %>% # make names same across sites
  rename_at(vars(Coral_sampler:Notes_Processing), .funs =  ~ paste(.x,"May", sep = "_"))

AugustData<-AugustData %>%
  select(-c(Lat,Lon, Sample_ID)) %>% # These are not complete for August
  rename_at(vars(Phosphate:Notes_Processing), .funs =  ~ paste(.x,"August", sep = "_"))

# Turb Data
TurbMay<-read_csv(here("Data","Turbinaria_CHN_May_2021_compiled.csv")) %>%
  select(Site_Number = Site_number,  Percent_N_May = Percent_N) # Just pull in the % N Data
  
TurbAugust<-read_csv(here("Data","Turb_CHN_Aug2020_compiled.csv")) %>%
  select(Site_Number,  Percent_N_August = Percent_N) %>%# Just pull in the % N Data
  mutate(Percent_N_August = Percent_N_August,
         Percent_N_August = as.numeric(str_sub(Percent_N_August,end = 4))) # remove the percent sign

# Join the data by sample ID

TurbAll<-left_join(TurbAugust, TurbMay)
WaterColumnAll<-left_join(MayData, AugustData)

# Join all nutrient data
NutrientAll<-left_join(WaterColumnAll, TurbAll) %>%
  select(Site_Number, Lat, Lon,Date_May, Time_May, Date_August, Time_August, Phosphate_May, Phosphate_August, Silicate_May, Silicate_August, Nitrite_plus_Nitrate_May, Nitrite_plus_Nitrate_August,Ammonia_May, Ammonia_August, Percent_N_May, Percent_N_August) %>%
  drop_na() # drop the not complete cases

write_csv(NutrientAll, here("Data","NutrientAll.csv"))

## Organize it
## Scale data and make run a cluster analysis ####

MayData_scaled <-MayData %>%
  select(Phosphate, Silicate, Nitrite_plus_Nitrate,Ammonia) %>%
  mutate_all(.funs = scale)

# Make a distanced based matrix
dist_mat <- dist(MayData_scaled, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')

ggdendrogram(hclust_avg, rotate = FALSE, size = 2)

# make 5 clusters
cut_avg <- cutree(hclust_avg, k = 5)

## k-means
clusters <- kmeans(MayData_scaled, 4) # cluster with 4 groups

MayData$clusters<-factor(clusters$cluster)

# Base Maps
API<-names(read_table(here("Data","API.txt")))
register_google(key = API) ### use your own API


M_coords<-data.frame(lon =	-149.83, lat = -17.55)
M1<-get_map(M_coords, maptype = 'satellite', zoom = 12)

map1<-ggmap(M1)+geom_point(data = MayData, mapping = aes(x=Lon, y=Lat, color = clusters), size = 2, alpha = .60)+
  xlab("")+
  ylab("")


## Run PCA
# Run the PCA
pca_V <- prcomp(MayData_scaled, scale. = FALSE, center = FALSE)

# Extract the scores and loadings
PC_scores <-as_tibble(pca_V$x[,1:2])
PC_loadings<-as_tibble(pca_V$rotation) %>%
  bind_cols(labels = rownames(pca_V$rotation))

# Join PC scores with the original data
MayData <-MayData %>%
  bind_cols(PC_scores)

p1<-MayData %>%
  ggplot(aes(x = PC1, y = PC2, color = clusters, shape = clusters))+
  geom_point() +
  coord_cartesian(xlim = c(-9, 3), ylim = c(-4, 9)) +
 # scale_shape_manual(values = c(22,16))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  ggforce::geom_mark_ellipse(
    aes(fill = clusters, label = paste(clusters), color = clusters), 
    alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# loadings plot 
p2<-PC_loadings %>%
  ggplot(aes(x=PC1, y=PC1, label=labels))+
  geom_segment(aes(x=0,y=0,xend=PC1*5,yend=PC2*5),
               arrow=arrow(length=unit(0.1,"cm")), color = "grey")+
  annotate("text", x = PC_loadings$PC1*5+0.1, y = PC_loadings$PC2*5+.1,
           label = PC_loadings$labels)+
  coord_cartesian(xlim = c(-4, 8), ylim = c(-4, 4)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1+p2