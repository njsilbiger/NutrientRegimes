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
library(ggridges)

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
  rename_at(vars(Phosphate:Notes_Processing), .funs =  ~ paste(.x,"August", sep = "_"))%>%
  group_by(Site_Number) %>% ## August has two sites that were collected twice and have high variablity (99 and 96).  I am going to take the average
  summarise_if(.predicate = is.numeric,.funs = mean)


# Turb Data
TurbMay<-read_csv(here("Data","Turbinaria_CHN_May_2021_compiled.csv")) %>%
  select(Site_Number = Site_number,  Percent_N_May = Percent_N) # Just pull in the % N Data
  
TurbAugust<-read_csv(here("Data","Turb_CHN_Aug2020_compiled.csv")) %>%
  select(Site_Number,  Percent_N_August = Percent_N) %>%# Just pull in the % N Data
  mutate(Percent_N_August = Percent_N_August,
         Percent_N_August = as.numeric(str_sub(Percent_N_August,end = 4))) %>% # remove the percent sign
  group_by(Site_Number) %>% ## August has two sites that were collected twice and have high variablity (99 and 96).  I am going to take the average
  summarise_if(.predicate = is.numeric,.funs = mean)


# Join the data by sample ID

TurbAll<-left_join(TurbAugust, TurbMay)
WaterColumnAll<-left_join(MayData, AugustData)

# Join all nutrient data
NutrientAll<-left_join(WaterColumnAll, TurbAll) %>%
  select(Site_Number, Lat, Lon,Date_May, Time_May, Phosphate_May, Phosphate_August, Silicate_May, Silicate_August, Nitrite_plus_Nitrate_May, Nitrite_plus_Nitrate_August,Ammonia_May, Ammonia_August, Percent_N_May, Percent_N_August) %>%
  drop_na() # drop the not complete cases

# write the cleaned data
write_csv(NutrientAll, here("Data","NutrientAll.csv"))
write_csv(NutrientAll_na, here("Data","NutrientAll_na.csv"))

## Scale the data and select just the numerics
NutrientAll_scaled<-NutrientAll %>%
   select_if(is.numeric) %>% # select just the data
   select(-c(Lat,Lon)) %>% #remove the lat lon data
   mutate_all(.funs = scale) # scale it

# k-means clusters with 1 to 12 clusters
KS_df<-tibble(k = 1:14, SS = NA) # make an empty dataframe to fill in

set.seed(200) # set the seed to get the same one every time
for (i in 1:14){
  clusters <- kmeans(NutrientAll_scaled, i) # cluster with i groups
  KS_df$SS[i]<-clusters$tot.withinss # extract total within SS for all groups
}

# make a skree plot
ggplot(KS_df, aes(x = k, y = SS)) +
  geom_point()+
  geom_line() +
  geom_hline(aes(yintercept = SS[8])) # put a line at 8 groups

### Maybe 8 is the best?

### OK, run the k-means cluster with 8 groups 
clusters <- kmeans(NutrientAll_scaled, 8)

# Bring in the groups into the original dataframe

NutrientAll <-NutrientAll %>%
  mutate(cluster = factor(clusters$cluster))

# Make some distribution plots
# Make the data long for ggridges
NutLong<-NutrientAll %>%
  pivot_longer(cols = Phosphate_May:Percent_N_August, names_to = "NutParam")

ggplot(NutLong,aes(x = value, y = NutParam, fill = NutParam))+
  geom_density_ridges(rel_min_height = 0.01)+
  theme_bw()+
  xlim(0,3.5)+
  theme(legend.position = "null")

# Looks like some silicate outliers
NutrientAll %>%
  filter(Silicate_August>3)

# Base Maps
API<-names(read_table(here("API.txt")))
register_google(key = API) ### use your own API


M_coords<-data.frame(lon =	-149.83, lat = -17.55)
M1<-get_map(M_coords, maptype = 'satellite', zoom = 12)

map1<-ggmap(M1)+geom_point(data = NutrientAll, mapping = aes(x=Lon, y=Lat, color = cluster), size = 2, alpha = .60)+
  xlab("")+
  ylab("")


## Run PCA
# Run the PCA
pca_V <- prcomp(NutrientAll_scaled, scale. = FALSE, center = FALSE)

# Extract the scores and loadings
PC_scores <-as_tibble(pca_V$x[,1:2])
PC_loadings<-as_tibble(pca_V$rotation) %>%
  bind_cols(labels = rownames(pca_V$rotation))

# Join PC scores with the original data
NutrientAll <-NutrientAll %>%
  bind_cols(PC_scores)

p1<-NutrientAll %>%
  ggplot(aes(x = PC1, y = PC2, color = cluster, shape = cluster))+
  geom_point() +
  coord_cartesian(xlim = c(-8, 8), ylim = c(-5, 5)) +
   scale_shape_manual(values = c(16:24))+
  scale_colour_hue(l = 45)+
  scale_fill_hue(l = 45)+
  ggforce::geom_mark_ellipse(
    aes(fill = cluster, label = paste(cluster), color = cluster), 
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
#  coord_cartesian(xlim = c(-4, 8), ylim = c(-4, 4)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p1+p2


