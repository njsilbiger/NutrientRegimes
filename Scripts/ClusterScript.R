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
  rename_at(vars(Coral_sampler:Notes_Processing), .funs =  ~ paste(.x,"May", sep = "_")) %>%
  mutate(Ammonia_May = ifelse(Ammonia_May>0,Ammonia_May,0)) # there are some negative values.  This makes them 0

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

## Read in coral data
CoralData<-read_csv(here("Data","Site_coral_categories.csv"))%>%
  rename(Site_Number = Site)

# Join the data by sample ID

TurbAll<-left_join(TurbAugust, TurbMay)
WaterColumnAll<-left_join(MayData, AugustData)

# Join all nutrient data
NutrientAll<-left_join(WaterColumnAll, TurbAll) %>%
  mutate(HighSi = ifelse(Silicate_May<1.5, "No","Yes")) %>%
#  filter(Silicate_May<1.5)%>% # remove low salinity points
 # filter(Silicate_August<1.5)%>%
  select(Site_Number, Lat, Lon,Date_May, Time_May, Phosphate_May, Phosphate_August, Silicate_May, Silicate_August, Nitrite_plus_Nitrate_May, Nitrite_plus_Nitrate_August,Ammonia_May, Ammonia_August, Percent_N_May, Percent_N_August, HighSi) %>%
  drop_na() # drop the not complete cases

# write the cleaned data
write_csv(NutrientAll, here("Data","NutrientAll.csv"))
#write_csv(NutrientAll_na, here("Data","NutrientAll_na.csv"))

## Scale the data and select just the numerics
NutrientAll_scaled<-NutrientAll %>%
  filter(HighSi == "No") %>%
   select_if(is.numeric) %>% # select just the data
   select(-c(Lat,Lon)) %>% #remove the lat lon data
   mutate_all(.funs = ~log(.x+0.01)) %>% # log transform it first
   mutate_all(.funs = scale) # scale it

# k-means clusters with 1 to 12 clusters
KS_df<-tibble(k = 1:8, SS = NA) # make an empty dataframe to fill in

set.seed(200) # set the seed to get the same one every time
for (i in 1:8){
  clusters <- kmeans(NutrientAll_scaled, i) # cluster with i groups
  KS_df$SS[i]<-clusters$tot.withinss # extract total within SS for all groups
}

# make a skree plot
ggplot(KS_df, aes(x = k, y = SS)) +
  geom_point()+
  geom_line() +
  geom_hline(aes(yintercept = SS[4])) # put a line at 8 groups
ggsave(here("Output","skreeplot.png"))

### Maybe 8 is the best?

### OK, run the k-means cluster with 8 groups 
clusters <- kmeans(NutrientAll_scaled, 4)

# Bring in the groups into the original dataframe

NutrientAll_low <-NutrientAll %>%
  filter(HighSi == "No") %>%
  mutate(cluster = clusters$cluster)

# add in the High silicate values as a new cluster
NutrientAll_high <-NutrientAll %>%
  filter(HighSi == "Yes") %>%
  mutate(cluster = 5)

NutrientAll <-bind_rows(NutrientAll_low, NutrientAll_high) %>%
  mutate(cluster = factor(cluster))

# Make some distribution plots
# Make the data long for ggridges
NutLong<-NutrientAll %>%
  pivot_longer(cols = Phosphate_May:Percent_N_August, names_to = "NutParam")

ggplot(NutLong,aes(x = value, y = NutParam, fill = NutParam))+
  geom_density_ridges(rel_min_height = 0.01)+
  theme_bw()+
  xlim(0,3.5)+
  scale_fill_viridis_d()+
  theme(legend.position = "null")
ggsave(here("Output","ridgeplot.png"))

# Looks like some silicate outliers
NutrientAll %>%
  filter(Silicate_August>3)

# Base Maps
#API<-names(read_table(here("API.txt")))
#register_google(key = API, write=TRUE) ### use your own API


M_coords<-data.frame(lon =	-149.83, lat = -17.55)
M1<-get_map(M_coords, maptype = 'satellite', zoom = 12)

map1<-ggmap(M1)+
  geom_point(data = NutrientAll, mapping = aes(x=Lon, y=Lat, color = cluster), size = 1.5)+
  scale_color_viridis_d(option = "B")+
  xlab("")+
  ylab("")

ggsave(filename =  here("Output","MapClusters.png"), plot =map1)

## Run PCA
# Run the PCA
pca_V <- prcomp(NutrientAll[,6:15], scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores <-as_tibble(pca_V$x[,1:2])
PC_loadings<-as_tibble(pca_V$rotation) %>%
  bind_cols(labels = rownames(pca_V$rotation))

# Join PC scores with the original data
NutrientAll <-NutrientAll %>%
  bind_cols(PC_scores)

p1<-NutrientAll %>%
  ggplot(aes(x = PC1, y = PC2, color = cluster, shape = cluster, fill = cluster))+
  geom_point() +
  coord_cartesian(xlim = c(-8, 8), ylim = c(-5, 5)) +
   scale_shape_manual(values = c(16:24))+
  scale_color_viridis_d(option = "B")+
  scale_fill_viridis_d(option = "B")+
#  scale_colour_hue(l = 45)+
#  scale_fill_hue(l = 45)+
  stat_ellipse(alpha = 0.15, geom = "polygon")+
  # ggforce::geom_mark_ellipse(
  #   aes(fill = cluster, label = paste(cluster), color = cluster), 
  #   alpha = .15, show.legend = FALSE,  label.buffer = unit(1, "mm"))+
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
ggsave(here("Output","PCAclust.png"), width = 8, height = 4)

map1 +(p1/p2)
ggsave(here("Output","PCAclust_map.png"), width = 12, height = 12)

## Make a set of boxplots for each cluster
NutLong %>%
  ggplot(aes(x = NutParam, y = value+.01, color = cluster))+
  geom_boxplot()+
  coord_trans(y="log")+
  facet_wrap(~cluster) +
  scale_color_viridis_d(option = "B")+
  ylab("Concentration")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = -.1)) 
ggsave(here("Output","boxplots_clust.png"), height = 10, width = 8)

# how many samples are in each cluster
clusters$size

# join in the coral data
NutrientAll_coral<-left_join(NutrientAll, CoralData)


# make a boxplot on corals by nutrient regime
NutrientAll_coral %>%
  ggplot(aes(x = cluster, y = log(Poc_count_adj+1)))+
  geom_boxplot()

NutrientAll_coral %>%
  ggplot(aes(x = cluster, y = Acr_count_adj))+
  geom_boxplot()

model<- multinom(Poc_cat~cluster, data = NutrientAll_coral)
tt <- broom::tidy(model,conf.int=TRUE)
#tt <- dplyr::filter(tt, term!="(Intercept)")

ggplot(tt, aes(x = y.level, y = estimate))+
  geom_point()+
  geom_errorbar(aes(ymin = estimate-std.error, ymax = estimate+std.error),width = 0.25 ) +
  coord_flip()+
  geom_hline(aes(yintercept = 0))+
  facet_wrap(~term)
