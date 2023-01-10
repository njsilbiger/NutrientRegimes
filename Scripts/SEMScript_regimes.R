## putting together structural equation model ####
#### By Nyssa Silbiger #####
#### Created on 1/8/2023 #####


library(tidyverse)
library(googlesheets4)
library(here)
library(brms)
library(bayesplot)
library(viridis)
library(patchwork)
library(ggridges)
library(ggtext)
library(broom) 
library(ggdist)


### read in data ####
data<-read_csv(here("Data","DataTurbWater_google.csv"))

data2021<- data %>% 
  filter(Year == 2021,
         Silicate <10)  %>%# remove outlier 
  drop_na(Island_shore)



### plot some relationships #####

## Si to nitrogen
data2021 %>%
  ggplot(aes(x = Silicate, y = Nitrite_plus_Nitrate))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~log(x))+
  coord_trans(x = "log")


## Distance to crest to nitrogen
data2021 %>%
  ggplot(aes(x = Distance_to_crest, y = Nitrite_plus_Nitrate))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x)
 # coord_trans(x = "log")


## microbes to nutrients
data2021 %>%
  ggplot(aes(y = Nitrite_plus_Nitrate, x = Microbial_PCoA_Axis1))+
  geom_point()+
  geom_smooth(method = "lm", formula = y~x)


### NN to fDom
  data2021 %>%
  ggplot(aes(x = Nitrite_plus_Nitrate, y = log(Marine.Humic.like)))+
  geom_point()+
  geom_smooth(method = "lm")

  ## fdom to microbes
  data2021 %>%
   # filter(Marine.Humic.like<0.047) %>%
    ggplot(aes(x = log(Marine.Humic.like), y = Microbial_PCoA_Axis1))+
    geom_point()+
    geom_smooth(method = "lm")

  ### Start the Bayesian Model ####
  
  DataModel<- data2021 %>%
    select(Site, Silicate,Distance_to_crest,Island_shore,Marine.Humic.like,Nitrite_plus_Nitrate,Microbial_PCoA_Axis1, Microbial_Species_Richness ) %>%
    mutate(Silicate = log(Silicate), # log transform the silicate and marine humics
           Marine.Humic.like = log(Marine.Humic.like),
           Microbial_Species_Richness = log(Microbial_Species_Richness)) %>%
    mutate_at(vars(Silicate,Distance_to_crest,Marine.Humic.like,Nitrite_plus_Nitrate,Microbial_PCoA_Axis1,Microbial_Species_Richness), function(x) scale(x, center = TRUE, scale = TRUE)) %>%
    drop_na()
  
  # write models 
NNmod<-bf(Nitrite_plus_Nitrate~Silicate + Distance_to_crest+Island_shore+ Microbial_Species_Richness, family = "student")
fDOMmod<-bf(Marine.Humic.like~Silicate +Nitrite_plus_Nitrate, family = "student")
Microbemod<-bf(Microbial_PCoA_Axis1~Marine.Humic.like, family = "student") # possibly use richness
Microbemod<-bf(Microbial_Species_Richness~Marine.Humic.like, family = "student") # 

# run models
fit_brms <- brm(
  NNmod+
    fDOMmod+
    Microbemod+
    set_rescor(FALSE),
  data=DataModel, cores=4, chains = 3)

fit_loo<-loo(fit_brms, reloo = TRUE) # looks good!

p1<-pp_check(fit_brms, resp="NitriteplusNitrate") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("N+N")

p2<-pp_check(fit_brms, resp="MarineHumiclike") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Humics")

p3<-pp_check(fit_brms, resp="MicrobialPCoAAxis1") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Microbes")

p3<-pp_check(fit_brms, resp="MicrobialSpeciesRichness") +
  scale_color_manual(values=c("red", "black"))+
  ggtitle("Microbes")

p1 + p2 + p3
##### Make Marginal effects plots ####

## NN ~ Si
R<-conditional_effects(fit_brms, "Silicate", resp = "NitriteplusNitrate",  resolution = 1000)
R1<-R$NitriteplusNitrate.NitriteplusNitrate_Silicate %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center"),
         lower = lower__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center"),
         upper = upper__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center")) %>%
  mutate(Silicate = Silicate*attr (DataModel$Silicate,"scaled:scale")+attr(DataModel$Silicate,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = exp(Silicate), y = estimate), linewidth = 1, color = 'grey')+
  geom_ribbon(aes(x = exp(Silicate),ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = data2021, aes(x = Silicate, y = Nitrite_plus_Nitrate, color = Distance_to_crest, shape = Island_shore)) +
  xlab(expression(atop("Silicate", paste("(",mu,"mol L"^-1,")"))))+
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
  #   ggtitle("Model 2")+
  coord_trans(x = "log")+
  # scale_x_continuous(breaks = c(0.2,1,5,10,25))+
  #  scale_y_continuous(breaks = c(0,0.1,1,5,10,30))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

## NN ~ dist to crest
R<-conditional_effects(fit_brms, "Distance_to_crest", resp = "NitriteplusNitrate",  resolution = 1000)
R2<-R$NitriteplusNitrate.NitriteplusNitrate_Distance_to_crest %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center"),
         lower = lower__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center"),
         upper = upper__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center")) %>%
  mutate(Distance_to_crest = Distance_to_crest*attr (DataModel$Distance_to_crest,"scaled:scale")+attr(DataModel$Distance_to_crest,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = Distance_to_crest, y = estimate), linewidth = 1, color = 'grey')+
  geom_ribbon(aes(x = Distance_to_crest,ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = data2021, aes(x = Distance_to_crest, y = Nitrite_plus_Nitrate, color = Silicate, shape = Island_shore)) +
  xlab(expression("Distance to shore (m)"))+ 
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

## Nutrients by side of island
R<-conditional_effects(fit_brms, "Island_shore", resp = "NitriteplusNitrate",  resolution = 1000)
R3<-R$NitriteplusNitrate.NitriteplusNitrate_Island_shore%>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center"),
         lower = lower__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center"),
         upper = upper__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center")) %>%
  ggplot(aes(x = Island_shore, y = estimate))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .2)+
  xlab(expression("Island Shore"))+ 
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

## fDOM ~ NN
R<-conditional_effects(fit_brms, "Nitrite_plus_Nitrate", resp = "MarineHumiclike",  resolution = 1000)
R4<-R$MarineHumiclike.MarineHumiclike_Nitrite_plus_Nitrate %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(DataModel$Marine.Humic.like,"scaled:scale")+attr(DataModel$Marine.Humic.like,"scaled:center"),
         lower = lower__*attr(DataModel$Marine.Humic.like,"scaled:scale")+attr(DataModel$Marine.Humic.like,"scaled:center"),
         upper = upper__*attr(DataModel$Marine.Humic.like,"scaled:scale")+attr(DataModel$Marine.Humic.like,"scaled:center")) %>%
  mutate(Nitrite_plus_Nitrate = Nitrite_plus_Nitrate*attr (DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = Nitrite_plus_Nitrate, y = exp(estimate)), linewidth = 1, color = 'grey')+
  geom_ribbon(aes(x = Nitrite_plus_Nitrate,ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = data2021, aes(x = Nitrite_plus_Nitrate, y = Marine.Humic.like)) +
  ylab(expression(atop("fDOM", "(Ramen units)")))+
  xlab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
  #   ggtitle("Model 2")+
  coord_trans(y = "log")+
   theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

## microbes ~ fDOM
## fDOM ~ NN
R<-conditional_effects(fit_brms, "Marine.Humic.like", resp = "MicrobialSpeciesRichness",  resolution = 1000)
R5<-R$MicrobialSpeciesRichness.MicrobialSpeciesRichness_Marine.Humic.like %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(DataModel$Microbial_Species_Richness,"scaled:scale")+attr(DataModel$Microbial_Species_Richness,"scaled:center"),
         lower = lower__*attr(DataModel$Microbial_Species_Richness,"scaled:scale")+attr(DataModel$Microbial_Species_Richness,"scaled:center"),
         upper = upper__*attr(DataModel$Microbial_Species_Richness,"scaled:scale")+attr(DataModel$Microbial_Species_Richness,"scaled:center")) %>%
  mutate(Marine.Humic.like = Marine.Humic.like*attr (DataModel$Marine.Humic.like,"scaled:scale")+attr(DataModel$Marine.Humic.like,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = exp(Marine.Humic.like), y = exp(estimate)), linewidth = 1, color = 'grey')+
  geom_ribbon(aes(x = exp(Marine.Humic.like),ymin=exp(lower), ymax=exp(upper)), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = data2021, aes(x = Marine.Humic.like, y = Microbial_Species_Richness)) +
  xlab(expression(atop("fDOM", "(Ramen units)")))+
  ylab("Microbial community (Richness)")+
  #   ggtitle("Model 2")+
  coord_trans(x = "log", y = "log")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')

## NN ~ microbes
R<-conditional_effects(fit_brms, "Microbial_Species_Richness", resp = "NitriteplusNitrate",  resolution = 1000)
R5<-R$NitriteplusNitrate.NitriteplusNitrate_Microbial_Species_Richness %>% # back transform the scaled effects for the plot
  mutate(estimate = estimate__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center"),
         lower = lower__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center"),
         upper = upper__*attr(DataModel$Nitrite_plus_Nitrate,"scaled:scale")+attr(DataModel$Nitrite_plus_Nitrate,"scaled:center")) %>%
  mutate(Microbial_Species_Richness = Microbial_Species_Richness*attr (DataModel$Microbial_Species_Richness,"scaled:scale")+attr(DataModel$Microbial_Species_Richness,"scaled:center")
  )%>%
  ggplot()+ # back trasform the log transformed data for better visual
  geom_line(aes(x = exp(Microbial_Species_Richness), y = estimate), linewidth = 1, color = 'grey')+
  geom_ribbon(aes(x = exp(Microbial_Species_Richness),ymin=lower, ymax=upper), linetype=1.5, alpha=0.3, fill = "grey")+
  geom_point(data = data2021, aes(x = Microbial_Species_Richness, y = Nitrite_plus_Nitrate)) +
  ylab(expression(atop("Nitrate + Nitrite", paste("(",mu,"mol L"^-1,")"))))+
  xlab("Microbial community (Richness)")+
  #   ggtitle("Model 2")+
  coord_trans(x = "log")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = 'bottom')


R1 + R2 + R5

## Coefficient plot
post <- as_draws_df(fit_brms)

Cof<-post %>% 
  select(starts_with("b"),-ends_with("Intercept")) %>%
  gather() %>% 
  group_by(key)%>%
  median_hdci()%>%
  mutate(sig = ifelse(sign(.lower)==sign(.upper),'yes','no'))%>%# if not significant make it grey
  separate(col = key,into = c("b", "dependent", "independent"),sep = "_") %>%
  mutate(dependent = factor(dependent, levels = c("NitriteplusNitrate","MarineHumiclike","MicrobialSpeciesRichness")))%>%
  mutate(dependent = recode(dependent,"NitriteplusNitrate" = "NOx","MarineHumiclike" = "fDOM","MicrobialSpeciesRichness" = "Microbial community")) 

Cof%>%
  ggplot(aes(x = value, y = independent, shape = sig)) +  # note how we used `reorder()` to arrange the coefficients
  geom_vline(xintercept = 0, alpha = 1/10, color = 'firebrick4') +
  geom_point(size = 2)+
  geom_errorbarh(aes(xmin = .lower, xmax = .upper), height = 0)+
  #scale_alpha_manual(values = c(0.2,1))+
  scale_shape_manual(values = c(16,1))+
  #scale_color_brewer(palette = "Set2")+
  scale_color_manual(values = c("firebrick4"), name = " ")+
  labs(x = NULL, y = NULL) +
  theme_bw() +
  guides(shape = "none")+
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = alpha("firebrick4", 1/4), linetype = 3),
        axis.text.y  = element_text(hjust = 0),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=14),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold")
  )+
  facet_grid(~dependent, scales = "free_y", space='free')  
