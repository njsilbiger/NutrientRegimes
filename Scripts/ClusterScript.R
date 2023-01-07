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
MayData<-read_csv(here("Data","Lagoon_nutrients_May_2021.csv")) %>%
  mutate(TimePeriod = "2021-May",
         Year = 2021)%>%
  rename(Site = Site_number)

AugustData<-read_csv(here("Data","Lagoon_nutrients_August_2020.csv"))%>%
mutate(TimePeriod = "2020-August",
       Year = 2020) %>%
  rename(Site = Site_Number)

AprilData<-read_csv(here("Data","Lagoon_nutrients_April_2022.csv"))%>%
  mutate(TimePeriod = "2022-April",
         Year = 2022)


AllWater<-bind_rows(AugustData,MayData,AprilData) %>%
  select(Site,Phosphate, Silicate, Nitrite_plus_Nitrate, Ammonia, Date, Time, TimePeriod, Year)
## Make a long table with all the variables
# Turb Data

TurbMay<-read_csv(here("Data","Turbinaria_CHN_May_2021_compiled.csv"))  %>%
  mutate(TimePeriod = "2021-May",
         Year = 2021)%>%
  rename(Site = Site_number)


TurbAugust<-read_csv(here("Data","Turb_CHN_Aug2020_compiled.csv")) %>%
  mutate(TimePeriod = "2020-August",
         Year = 2020)%>%
  group_by(Site_Number, TimePeriod, Year) %>% ## August has two sites that were collected twice and have high variability (99 and 96).  I am going to take the average
  summarise_if(.predicate = is.numeric,.funs = mean) %>%
  ungroup() %>%
  rename(Site = Site_Number) 

TurbApril<- read_csv(here("Data","Turb_CHN_compiled_April2022.csv"))%>%
  mutate(TimePeriod = "2022-April",
         Year = 2022)%>%
  rename(Site = Site_Number)

Turb2016<- read_csv(here("Data","Turb2016_2019.csv"))

TurbAll <- bind_rows(Turb2016, TurbAugust, TurbMay, TurbApril) %>%
  select(Site, Year, TimePeriod, dN15, Percent_N, Percent_C, C_to_N_ratio)

## Read in the metadata
meta<-read_csv(here("Data","islandwide_sites_allMetadata - islandwide_sites_allMetadata.csv"))

## Join all the data together
AllData<-full_join(TurbAll, AllWater) %>%
  full_join(meta) %>%
  drop_na(Year) %>%
  select(-c("Jan-16":"May-21"))


## Make a ridge plot
AllData %>%
 # filter(TimePeriod != "2020-August") %>%
  pivot_longer(cols = dN15:Ammonia) %>%
  select(TimePeriod, name, value)%>%
  drop_na() %>%
ggplot(aes(x = value, fill = TimePeriod))+
  geom_density(alpha = 0.5)+
  theme_bw()+
 # coord_trans(x= "log")+
  #xlim(0,3.5)+
  scale_fill_viridis_d()+
  facet_wrap(~name, scales = "free")
ggsave(here("Output","ridgeplot.png"), width = 8, height = 6)

# water column nutrient data is off for August 2020... might need to drop

AllData %>%
  filter(TimePeriod != "2020-August") %>%
  drop_na(Habitat, Island_shore,Nitrite_plus_Nitrate) %>%
  droplevels()%>%
 # filter(Habitat != "Bay") %>%
  ggplot(aes(x = Percent_N, y = Nitrite_plus_Nitrate, color = TimePeriod))+
  geom_point()+
  geom_smooth(method = "glm", formula = y~log(x),
              method.args = list(family = gaussian(link = 'log')))+
  coord_trans(y = "log")+
  #geom_smooth(method = "lm")+
  facet_wrap(~TimePeriod, scale = "free")
ggsave(here("Output","NNvspercentN.png"), width = 6, height = 4)

AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  ggplot(aes(x = Distance_to_population_center, y = Percent_N, color = TimePeriod))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Distance to population center (m)",
       y = "%N")+
  facet_wrap(Island_shore~Year, scale = "free", ncol = 5) +
  theme_bw()

ggsave(here("Output","Population_NPercent.pdf"), width = 10, height = 6)

AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  ggplot(aes(x = Distance_to_shore, y = Percent_N, color = TimePeriod))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(x = "Distance to population center (m)",
       y = "%N")+
  facet_wrap(Island_shore~Year, scale = "free", ncol = 5) +
  theme_bw()

## Calculate the slope of the relationship between %N and population center... turn this into Bayesian
Pop_coefs<-AllData %>%
  group_by(Year, Site,Island_shore) %>%
  summarise(Percent_N = mean(Percent_N, na.rm = TRUE), # take average by year for 2016 data
            Distance_to_population_center = mean(Distance_to_population_center, na.rm = TRUE)) %>%
  ungroup()%>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore) %>%
  droplevels()%>%
  group_by(Year, Island_shore) %>%
  nest() %>% 
  mutate(model = 
  map(data,function(df) {
    lm(Percent_N ~ Distance_to_population_center, data = df)})) %>%
  mutate(glance = map(model, broom::tidy)) %>%
  unnest(glance) %>%
  select(Year, Island_shore, term, estimate, std.error, statistic,p.value) %>%
  pivot_wider(values_from = c(estimate:p.value), names_from = term) %>%
  janitor::clean_names() 

Pop_coefs%>%
  ggplot(aes(x = estimate_distance_to_population_center, y = year))+
  geom_vline(xintercept = 0, lty =2)+
  geom_point()+
  geom_errorbar(aes(xmin =estimate_distance_to_population_center-std_error_distance_to_population_center,
                    xmax = estimate_distance_to_population_center+std_error_distance_to_population_center), width = 0)+
  facet_wrap(~island_shore)

Pop_coefs%>%
  ggplot(aes(x = estimate_intercept, y = year))+
  geom_point()+
  geom_errorbar(aes(xmin =estimate_intercept-std_error_intercept,
                    xmax = estimate_intercept+std_error_intercept), width = 0)+
  facet_wrap(~island_shore)


AllData %>%
  ggplot()+
  geom_violin(aes(x = factor(Year), y = Percent_N, fill =factor(Year) ))+
  facet_wrap(~Island_shore)

AllData %>%
  ggplot()+
  geom_boxplot(aes(x = factor(Year), y = Percent_N, fill =factor(Year)), )


p1<-AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>% # the august nitrate data is bad :(
  ggplot(aes(x = Silicate, y = (Nitrite_plus_Nitrate+Ammonia)/Phosphate, color = factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~log(x))+
  labs(y = "N:P")+
  coord_trans(x = "log")+
  facet_wrap(~Island_shore) +
  theme_bw()


p2<-AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>% # the august nitrate data is bad :(
  ggplot(aes(x = Silicate, y = Ammonia, color = factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~log(x))+
  coord_trans(x = "log")+
  facet_wrap(~Island_shore) +
  theme_bw()

p3<-AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>% # the august nitrate data is bad :(
  ggplot(aes(x = Silicate, y = Nitrite_plus_Nitrate, color = factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~log(x))+
  coord_trans(x = "log")+
  facet_wrap(~Island_shore) +
  theme_bw()

# Distance to crest
p4<-AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>% # the august nitrate data is bad :(
  ggplot(aes(x = Distsance.to.crest, y = Nitrite_plus_Nitrate, color = factor(Year)))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x)+
#  coord_trans(x = "log")+
  facet_wrap(~Island_shore, scale = "free_x") +
  theme_bw()


p3/p4 + plot_layout(guides = "collect")
ggsave(here("Output","N_si_crest.png"), width = 6, height = 6)


p1/p2/p3 + plot_layout(guides = "collect")
ggsave(here("Output","N_silicate.pdf"), width = 6, height = 6)

## Read in the turb data from all years and calcualte the mean and variance for each site
TurbAllyears<-read_csv(here("Data","turb_N_allYears.csv"))

AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>% # the august nitrate data is bad :(
  ggplot(aes(x = Distance_to_shore, y = Nitrite_plus_Nitrate))+
  geom_point(aes(color = Silicate))+
  geom_smooth(method = "lm", formula = y~log(x))+
  coord_trans(x = "log")+
  facet_wrap(Year~Island_shore, scales = "free_y") +
  theme_bw()


### plot relationship between distance to shore and N--- look to see if residuals relate to silicate

AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>% # the august nitrate data is bad :(
  ggplot(aes(x = Distance_to_shore, y = Nitrite_plus_Nitrate))+
  geom_point(aes(color = Silicate))+
  geom_smooth(method = "lm", formula = y~x)+
  scale_color_continuous(trans = "log10",type = "viridis")+
#  coord_trans(x = "log")+
  labs(x = "Distance from shore (m)",
       y = "NN umolL-1",
        color = "Silicate")+
  facet_wrap(Year~Island_shore, scale = "free") +
  theme_bw()
ggsave(here("Output","distancevsNN.pdf"), width = 8, height = 5)

AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>% # the august nitrate data is bad :(
  mutate(distblock = ifelse(Distance_to_shore<500,"< 500 m","far"))%>%
  ggplot(aes(x = Silicate, y = Nitrite_plus_Nitrate, color = distblock))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~log(x))+
  #scale_color_continuous(trans = "log10",type = "viridis")+
  coord_trans(x = "log")+
  labs(x = "Silicate",
       y = "NN umolL-1",
       color = "Distance to shore")+
  facet_wrap(Year~Island_shore, scale = "free") +
  theme_bw()
ggsave(here("Output","SivsNN_block.pdf"), width = 8, height = 5)


AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  droplevels()%>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>% # the august nitrate data is bad :(
  ggplot(aes(x = Distance_to_shore, y = Silicate))+
  geom_point()+
  geom_smooth(method = "glm", formula = y~x,
              method.args = list(family = gaussian(link = 'log')))+
  scale_color_continuous(trans = "log10",type = "viridis")+
    coord_trans(y = "log")+
  labs(x = "Distance from shore (m)",
       y = "Si umolL-1",
      # color = "Silicate"
       )+
  facet_wrap(Year~Island_shore, scale = "free") +
  theme_bw()
ggsave(here("Output","DistanceVsSi.pdf"), width = 8, height = 5)

AllData %>%
  #  filter(TimePeriod != "2020-August") %>%
  drop_na(Island_shore,Silicate) %>%
  #filter(Habitat != "Bay") %>%
  filter(Silicate<4, Nitrite_plus_Nitrate<2.5, Phosphate <0.6,# drop a few outliers
         Year != 2020)%>%
  droplevels()%>%
  group_by(Year, Island_shore)%>%
  nest()%>%
  mutate(model = map(data, function(x){lm(Nitrite_plus_Nitrate~Distance_to_shore, data = x)})) %>%
  mutate(resids = map(model, broom::augment)) %>%
  select(Year, Island_shore, data, resids)%>%
  unnest() %>%
  ggplot(aes(x = Silicate, y=.std.resid, color = factor(Year)))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_point()+
  geom_smooth(method = "lm", formula = y~log(x))+
  coord_trans(x = "log")+
  labs(y = "residuals(NN ~ Distance from shore)\n pos values = N is higher than expected given distance to shore",
       color = "Year")+
  facet_wrap(Year~Island_shore, scales = "free") +
  theme_bw()
ggsave(here("Output","residualsNNSi.pdf"), width = 8, height = 5)


########## Old stuff ##################

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

## Read in the turb data from all years and calcualte the mean and variance for each site
TurbAllyears<-read_csv(here("Data","turb_N_allYears.csv"))

TurbAllyears<-TurbAllyears %>%
  group_by(Site)%>%
  summarise(meanN = mean(N, na.rm = TRUE),
            varN = var(N, na.rm = TRUE))%>%
  rename(Site_Number = Site)


## Read in coral data
CoralData<-read_csv(here("Data","Site_coral_categories.csv"))%>%
  rename(Site_Number = Site)

# Join the data by sample ID

TurbAll<-left_join(TurbAugust, TurbMay)
WaterColumnAll<-left_join(MayData, AugustData)

## Create a dataset with the may water columne data and the tirn data

May_turbdata<-left_join(MayData, TurbAllyears) %>%
  mutate(HighSi = ifelse(Silicate_May<1.5, "No","Yes")) %>%
  #  filter(Silicate_May<1.5)%>% # remove low salinity points
  # filter(Silicate_August<1.5)%>%
  select(Site_Number, Lat, Lon,Date_May, Time_May, Phosphate_May,Silicate_May,  Nitrite_plus_Nitrate_May, meanN,varN, HighSi) %>%
  drop_na() # drop the not complete cases

ggplot(May_turbdata, aes(x = meanN, y = Nitrite_plus_Nitrate_May))+
  geom_point()+
  geom_smooth(method = "lm")+
  annotate("text", x = .75, y = .75, label = "p < 0.001")+
  annotate("text", x = .75, y = .85, label = "r = 0.1")+
  labs(y = "N+N (umol L-1) in May",
      x = "Mean Turb %N across all years")+
  theme_bw()
  ggsave(here("Output","NNvsTurb.png"))
  

mod1<-lm(Nitrite_plus_Nitrate_May~meanN, data = May_turbdata )
anova(mod1)
summary(mod1)


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

## Scale the data 
May_turbdata_scaled <- May_turbdata %>%
  select_if(is.numeric) %>% # select just the data
  select(-c(Lat,Lon)) %>% #remove the lat lon data
  mutate_all(.funs = ~log(.x+0.01)) %>% # log transform it first
  mutate_all(.funs = ~scale(.x, center = FALSE)) # scale it
  

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

### Same thing with the may only and turb data ###
for (i in 1:8){
  clusters <- kmeans(May_turbdata_scaled, i) # cluster with i groups
  KS_df$SS[i]<-clusters$tot.withinss # extract total within SS for all groups
}

# make a skree plot
ggplot(KS_df, aes(x = k, y = SS)) +
  geom_point()+
  geom_line() +
  geom_hline(aes(yintercept = SS[5])) # put a line at 8 groups
ggsave(here("Output","skreeplot2.png"))

### Maybe 8 is the best?

### OK, run the k-means cluster with 5 groups 
clusters <- kmeans(May_turbdata_scaled, centers = 5)

# Bring into dataframe
May_turbdatawClusters<-May_turbdata %>%
  mutate(cluster = factor(clusters$cluster))


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
  ylab("")+
  ggtitle("Nyssa")

ggsave(filename =  here("Output","MapClusters.png"), plot =map1, width = 5, height = 5)

map1_may<-ggmap(M1)+
  geom_point(data = May_turbdatawClusters, mapping = aes(x=Lon, y=Lat, color = cluster), size = 2)+
 # scale_color_viridis(discrete = TRUE, option = "plasma")+
  scale_color_manual(values =  c("orange","yellow","pink","green","purple"))
  xlab("")+
  ylab("")

ggsave(filename =  here("Output","MapClusters_MayTurb.png"), plot =map1_may, width = 5, height = 5)

## Run PCA
# Run the PCA
pca_V <- prcomp(NutrientAll[,6:15], scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores <-as_tibble(pca_V$x[,1:2])
PC_loadings<-as_tibble(pca_V$rotation) %>%
  bind_cols(labels = rownames(pca_V$rotation))

# Join PC scores with the original data
NutrientAll<-NutrientAll %>%
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


#### for turb data
# Run the PCA
pca_V <- prcomp(log(May_turbdatawClusters[,6:10]), scale. = TRUE, center = TRUE)

# Extract the scores and loadings
PC_scores <-as_tibble(pca_V$x[,1:2])
PC_loadings<-as_tibble(pca_V$rotation) %>%
  bind_cols(labels = rownames(pca_V$rotation))

# Join PC scores with the original data
May_turbdatawClustersPCA<-May_turbdatawClusters %>%
  bind_cols(PC_scores) 

p1<-May_turbdatawClustersPCA %>%
  ggplot(aes(x = PC1, y = PC2, color = cluster, shape = cluster, fill = cluster))+
  geom_point() +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
  scale_shape_manual(values = c(16:24))+
  scale_color_manual(values =  c("orange","yellow","pink","green","purple"))+
  scale_fill_manual(values =  c("orange","yellow","pink","green","purple"))+
  
#  scale_color_viridis_d(option = "B")+
#  scale_fill_viridis_d(option = "B")+
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
ggsave(here("Output","PCAclust2.png"), width = 8, height = 4)

map1_may +(p1/p2)
ggsave(here("Output","PCAclust_map2.png"), width = 12, height = 12)


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

## Make a set of boxplots for each cluster
May_turbdatawClusters %>%
  pivot_longer(cols = Phosphate_May:varN, names_to = "NutParam")%>%
  ggplot(aes(x = NutParam, y = value+.01, color = cluster))+
  geom_boxplot()+
  coord_trans(y="log")+
  facet_wrap(~cluster) +
  scale_color_manual(values =  c("orange","yellow","pink","green","purple"))+
#  scale_color_viridis_d(option = "B")+
  ylab("Concentration")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = -.1)) 
ggsave(here("Output","boxplots_clust2.png"), height = 10, width = 8)


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

### Make some plots showing relationship between May and August data

ggplot(NutrientAll, aes(x = Nitrite_plus_Nitrate_May, y = Nitrite_plus_Nitrate_August))+
  geom_point()

ggplot(NutrientAll, aes(x = Silicate_May, y = Silicate_August))+
  geom_point()+
  coord_trans("log10")

ggplot(NutrientAll, aes(x = Ammonia_May, y = Ammonia_August))+
  geom_point()
