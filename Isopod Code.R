setwd("/Users/Nao/Desktop/Isopod/")

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(lme4)
library(nlme)
library(tidyverse)
library(tidyr)
library(multcomp)
library(forcats)

Isopod <- read.csv('Isopod Data.csv', header = T)
Isopod$Location <- as.factor(Isopod$Location)
Isopod$Quadrat <- as.factor(Isopod$Quadrat)
Isopod$Moon_Rise <- as.factor(Isopod$Moon_Rise)
Isopod$Cloud <- as.factor(Isopod$Cloud)
Isopod$Direct_Moonlight <- as.factor(Isopod$Direct_Moonlight)
Isopod$Direct_ALAN <- as.factor(Isopod$Direct_ALAN)

Isopod$Cloud[Isopod$Cloud == 'Y (Partial)'] <- 'N' #Treat partial cloud cover as no cloud cover

Isopod1 <- Isopod %>%
  mutate(Light_Categories = ifelse(Direct_ALAN == 'Y', 'ALAN', 
                                   ifelse(Cloud == 'Y', 'Clouded',
                                          ifelse(Moon_Rise == 'Y', 'Moonlit', 'Moonless')))) #Add light categories
Isopod1$Light_Categories <- as.factor(Isopod1$Light_Categories)

Isopod1$Total <- Isopod1$In+Isopod1$Out #Add Total Abundance

summary(Isopod1)
#Implausible Temp on row 38
#Implausible Lux on row 41
Isopod2 <- Isopod1[-41,] #Remove implausible Lux value

Isopod3 <- Isopod1[-38,] %>% 
  filter(!is.na(Humidity)) #Remove implausible Temp and NA Humidity

Isopod4 <- Isopod1[-38,] #Removed implausible Temp

########################## LUX & LIGHT CATEGORY ANALYSIS########################
### Data Overview
ggplot(Isopod2, aes(x=Light_Categories, y=log(Lux), fill=Light_Categories))+
  geom_boxplot()
hist(log(Isopod2$Lux), breaks = 15) #gaussian ok?

### Model Lux and Light Category
m1 <- glmer(log(Lux) ~ Light_Categories + (1|Quadrat), data = Isopod2, family = gaussian)
summary(m1)

m2 <- glmer(log(Lux) ~ 1 + (1|Quadrat), data = Isopod2, family = gaussian)
anova(m1, m2)  #Light category is a significant predictor of Lux p<0.001

###LUX ~ Light Category MAM
LuxMAM <- glmer(log(Lux) ~ Light_Categories + (1|Quadrat),
                data = Isopod2, family = gaussian)
summary(LuxMAM)

plot(resid(LuxMAM)) #Residuals pretty evenly distributed
plot(LuxMAM, resid(., scaled=TRUE) ~
       fitted(.) | Quadrat, abline = 0) #Residuals by quadrats
plot(LuxMAM, resid(., scaled=TRUE) ~
       fitted(.) | Location, abline = 0) #Residuals by location

#residual normality
qqnorm(resid(LuxMAM)) 
qqline(resid(LuxMAM)) #Skewed Tails

#confidence interval
exp(confint(LuxMAM)) #Doesn't span 0

###Graph Lux MAM 
Figure1 <- Isopod2 %>%
  mutate(Light_Categories = fct_relevel(Light_Categories,
                                        "ALAN", "Clouded", "Moonlit", "Moonless")) %>%
  ggplot(aes(x = Light_Categories, y = log(Lux), col = Light_Categories)) +
  geom_boxplot(data = Isopod2, aes(x = Light_Categories, y = log(Lux), group = Light_Categories), show.legend = FALSE)+
  theme_bw(base_size = 14) + 
  labs(x = "Light Condition", y = "Log Light Intensity (Lux)") +
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))+
  geom_text(aes(x = 2.5, y = 0.2,
                label = "a                 b                c                 a"),
            stat = "unique",
            size = 5, color = "black")
  
Figure1


###Tukey Test
post_hoc <- glht(LuxMAM, linfct = mcp(Light_Categories = "Tukey"))
summary(post_hoc)
# Clouded - ALAN == 0      -0.8883     0.2102  -4.227   <0.001 ***
# Moonless - ALAN == 0     -2.2239     0.2292  -9.703   <0.001 ***
# Moonlit - ALAN == 0      -0.3757     0.2304  -1.631   0.3553    
# Moonless - Clouded == 0  -1.3355     0.1640  -8.143   <0.001 ***
# Moonlit - Clouded == 0    0.5126     0.1665   3.078   0.0109 *  
# Moonlit - Moonless == 0   1.8481     0.1929   9.583   <0.001 ***



############################### TOTAL ABUNDANCE ###############################
### Total Abundance Overview
#plot abundance at each quadrat
ggplot(Isopod1, aes(x=Quadrat, y=In+Out, fill=Quadrat))+
  geom_boxplot()

hist(Isopod1$Total, breaks = 30) # 0 heavy distribution

#plot abundance by log(lux)
ggplot(Isopod1, aes(x = log(Lux), y = Total, colour = Light_Categories))+
  geom_point(aes(x = log(Lux), y = Total, colour = Light_Categories))

#plot abundance by light category
ggplot(Isopod1, aes(x = Light_Categories, y = Total, colour = Light_Categories))+
  geom_boxplot(aes(x = Light_Categories, y = Total, colour = Light_Categories))+
  geom_point(aes(x = Light_Categories, y = Total, colour = Light_Categories))


### Total Abundance Mixed Effect Model
# Use Isopod 3 => removed implausible temp and NA humidity
m1 <- glmer(Total ~ Light_Categories*Temp*Humidity + (1|Quadrat),
            data = Isopod3, family = poisson)
summary(m1)

#check for over dispersion
#Method 1:
dispersion <- 608/98 #6.2 Overdispersed => use negative binomial
#Method 2: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#overdispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m1) 


m2 <- glmer.nb(Total ~ Light_Categories*Temp*Humidity + (1|Quadrat),
               data = Isopod3)
m3 <- glmer.nb(Total ~ Light_Categories + Temp + Humidity +
                 Light_Categories:Temp + Light_Categories:Humidity + Temp:Humidity +
                 (1|Quadrat),
               data = Isopod3)
anova(m2, m3) #3way interaction not significant p=0.7105

m4 <- update(m3, ~.- Temp:Humidity)
anova(m3, m4) #temp:hum not significant p=0.05326

m5 <- update(m4, ~.- Light_Categories:Humidity)
anova(m4, m5) #light:humidity not significant p=0.1492

m6 <- update(m5, ~.- Light_Categories:Temp)
anova(m5, m6) #light:temp not significant p=0.1603

m7 <- update(m6, ~.- Humidity)
anova(m6, m7) #humidity not significant p=0.7059 => refit model model with NA humididty values (Isopod4)

m1 <- glmer.nb(Total ~ Light_Categories*Temp + (1|Quadrat),
               data = Isopod4)
m2 <- glmer.nb(Total ~ Light_Categories + Temp + (1|Quadrat),
               data = Isopod4)
anova(m1, m2) #temp interaction not significant p=0.8197

m3 <- update(m2, ~.- Temp)
anova(m2, m3) # temp not sig p=0.131 => add row 38 back in

m4 <-glmer.nb(Total ~ Light_Categories + (1|Quadrat),
              data = Isopod1)

m5 <- glmer.nb(Total ~ 1 + (1|Quadrat),
               data = Isopod1)
anova(m4, m5) #light categories significant 0.02301

### Total Abundance MAM
TotalMAM <- glmer.nb(Total ~ Light_Categories + (1|Quadrat),
                     data = Isopod1)
summary(TotalMAM)

#Check Residuals
qqnorm(resid(TotalMAM)) 
qqline(resid(TotalMAM))

#Check confidence interval
exp(confint(TotalMAM, method = 'boot', nsim=1000)) # don't cross 0
#.sig01                   1.4149922 2.192121
#(Intercept)              3.1417797 6.735384
#Light_CategoriesClouded  0.7157853 1.563391
#Light_CategoriesMoonless 0.8834463 1.949250
#Light_CategoriesMoonlit  0.5277924 1.209239


#Tukey Test
post_hoc <- glht(TotalMAM, linfct = mcp(Light_Categories = "Tukey"))
summary(post_hoc) #only moonlit significanly lower than moonless
# Clouded - ALAN == 0      0.05291    0.19671   0.269  0.99294   
# Moonless - ALAN == 0     0.27832    0.20513   1.357  0.51714   
# Moonlit - ALAN == 0     -0.21902    0.21274  -1.030  0.72487   
# Moonless - Clouded == 0  0.22541    0.12666   1.780  0.27465   
# Moonlit - Clouded == 0  -0.27193    0.14135  -1.924  0.21062   
# Moonlit - Moonless == 0 -0.49734    0.15822  -3.143  0.00868 **



### Plot
Figure2 <- ggplot(Isopod1, aes(x = Light_Categories, y = log(Total), col = Light_Categories)) +
  geom_boxplot(aes(x = Light_Categories, y = log(Total), col = Light_Categories), size=0.5, shape = 3, show.legend = F) +
  theme_bw(base_size = 14)+
  labs(x = "Light Condition", y = "Log Total Sea Slater Count") +
  geom_text(aes(x = 3.5, y = 3.6,
                label = "*               *"),
            stat = "unique",
            size = 10, color = "black")+
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))
Figure2


########################## MICROHABITAT PREFRENCE ##############################
### Graph Proportion by Lux. 
Isopod1$OutProportion <- Isopod1$Out/Isopod1$Total
Isopod1$OutProportion[is.na(Isopod1$OutProportion)] <- 0

ggplot(Isopod1, aes(x = Light_Categories, y = OutProportion, fill = Quadrat))+
  geom_point(aes(x = Light_Categories, y = OutProportion, colour = Quadrat))


### Proportion Out of Crevice  Mixed Effect Model
# Use Data without improbable temp and NA humidities first
Isopod3$OutProportion <- Isopod3$Out/Isopod3$Total
Isopod3$OutProportion[is.na(Isopod3$OutProportion)] <- 0

summary(Isopod3)

m1 <- glmer(OutProportion ~ Light_Categories*Temp*Humidity + (1|Quadrat),
            data = Isopod3, family = binomial, weights = Total) #Produces error: Downdated VtV is not positive definite => Over fitted model?

Isopod3$ID <- seq.int(nrow(Isopod3))
Isopod3$ID <- as.factor(Isopod3$ID)

# remove 3 way interaction and include ID to counter overdispersion
m1 <- glmer(OutProportion ~ Light_Categories + Temp + Humidity +
              Light_Categories:Temp + Temp:Humidity + Light_Categories:Humidity + (1|Quadrat) + (1|ID),
            data = Isopod3, family = binomial, weights = Total)
summary(m1)
m2 <-  glmer(OutProportion ~ Light_Categories + Temp + Humidity +
              Light_Categories:Temp + Light_Categories:Humidity + (1|Quadrat) + (1|ID),
              data = Isopod3, family = binomial, weights = Total)
anova(m1, m2) #temp:hum not sig p=0.6182

m3 <- glmer(OutProportion ~ Light_Categories + Temp + Humidity +
              Light_Categories:Temp + (1|Quadrat) + (1|ID),
            data = Isopod3, family = binomial, weights = Total)
anova(m2, m3) #Light Cat: Hum not significant p=0.2731

m4 <- glmer(OutProportion ~ Light_Categories + Temp + Humidity + (1|Quadrat) + (1|ID),
            data = Isopod3, family = binomial, weights = Total)
anova(m3, m4) #Light Category:Hum not significant

m5 <- glmer(OutProportion ~ Light_Categories + Temp + (1|Quadrat) + (1|ID),
            data = Isopod3, family = binomial, weights = Total)
anova(m4, m5) #humidity significant p<0.001

m6<- glmer(OutProportion ~ Light_Categories + Humidity + (1|Quadrat) + (1|ID),
            data = Isopod3, family = binomial, weights = Total)
anova(m4, m6) # temp significant

m7 <- glmer(OutProportion ~ Temp + Humidity + (1|Quadrat) + (1|ID),
            data = Isopod3, family = binomial, weights = Total)
anova(m4, m7) #Light Category not significant p=0.2133

###Prop MAM
PropMAM <- glmer(OutProportion ~ Temp + Humidity + (1|Quadrat) + (1|ID),
                     data = Isopod3, family = binomial, weights = Total)
plot(resid(PropMAM))

#check confidence intervals => Error, bootstrap failed
confint(PropMAM)
confint(PropMAM, method = "boot", nsim = 1000) 


##Plot this
x <- seq(min(Isopod3$Temp),max(Isopod3$Temp),length.out=50)
y <- seq(min(Isopod3$Humidity),max(Isopod3$Humidity),length.out=50)
a <- levels(Isopod3$Quadrat)
b <- levels(as.factor(Isopod3$ID))
simData <-  expand.grid(Temp = x, Humidity = y, Quadrat = a, ID = b) #Create matrix with all combinations
simData <- mutate(simData, OutProportion = predict(PropMAM, data.frame(simData), type = 'response', allow.new.levels = TRUE)) #Add predicted values from model
simData <- mutate(simData, unique_id = group_indices(.data = simData, Temp, Humidity)) #Create an ID for each unique temp and humidity combination

simData2 <- simData %>% 
  group_by(unique_id) %>% 
  mutate(mean = mean(OutProportion)) #Find mean of repeated temp and humidity combinations
keeps <- c('Humidity', 'Temp','unique_id', 'mean', 'OutProportion')
simData2 <- simData2[keeps]
simData2 <- unique(simData2) #Clean up and remove columns and repeats not necessary for plot

Figure3a <- ggplot(simData2, aes(x = Temp, y = Humidity, z = OutProportion))+
  geom_tile(aes(fill = OutProportion))+
  geom_point(data = Isopod3, aes(x = Temp, y = Humidity))+
  scale_fill_continuous(type = 'viridis')+
  labs(x = 'Temperature (\u00B0C)',y = 'Humidity (%)', fill ='Proportion of 
  Sea Slaters
  Outside of 
  Crevices')+
  theme_bw()



### Proportion by Light and Temperature (reanalysis without humidity => bigger dataset)
# Use isopod4 (removed implausible temp value, keep NA humidity)
Isopod4$OutProportion <- Isopod4$Out/Isopod4$Total
Isopod4$OutProportion[is.na(Isopod4$OutProportion)] <- 0

m1 <- glmer(OutProportion ~ Light_Categories*Temp + (1|Quadrat),
            data = Isopod4, family = binomial, weights = Total)
summary(m1)
dispersal <- 316.4/150 #2.109 is overdispersed => use ID as random effect to control for this

Isopod4$ID <- seq.int(nrow(Isopod4))
Isopod4$ID <- as.factor(Isopod4$ID)

m2 <- glmer(OutProportion ~ Light_Categories*Temp + (1|Quadrat) + (1|ID),
            data = Isopod4, family = binomial, weights = Total)
m3 <- glmer(OutProportion ~ Light_Categories+Temp + (1|Quadrat) + (1|ID),
            data = Isopod4, family = binomial, weights = Total)
anova(m2, m3) #interaction between temp and light condition present p=0.001034

### Temp Proportion MAM
TempPropMAM <- glmer(OutProportion ~ Light_Categories*Temp + (1|Quadrat) + (1|ID),
                     data = Isopod4, family = binomial, weights = Total)
summary(TempPropMAM)
 
# Tukey
post_hoc <- glht(TempPropMAM, linfct = mcp(Light_Categories = "Tukey"))
summary(post_hoc) #moonless-ALAN and Moonless-Clouded Significant
# Clouded - ALAN == 0      -0.4308     2.2029  -0.196  0.99730   
# Moonless - ALAN == 0      7.6280     2.7090   2.816  0.02432 * 
# Moonlit - ALAN == 0       4.6101     2.8091   1.641  0.35141   
# Moonless - Clouded == 0   8.0587     2.3664   3.405  0.00371 **
# Moonlit - Clouded == 0    5.0409     2.5223   1.999  0.18557   
# Moonlit - Moonless == 0  -3.0178     2.9115  -1.037  0.72466   

#Check Confidence => error in creating confidence intervals
exp(confint(TempPropMAM, method = 'boot', nsim=1000))
exp(confint(TempPropMAM))

#Compare Slopes => same significant pairings as Tukey
library(emmeans)
emt = emtrends(TempPropMAM, 'Light_Categories', var = 'Temp')
emt    # estimated slopes
pairs(emt)

### Plot MAM
simDF <- data.frame(Temp = rep(seq(min(Isopod4$Temp),max(Isopod4$Temp),length.out=1000),4),
                    Light_Categories = c(rep('ALAN', 1000), rep('Clouded',1000), rep('Moonless', 1000), rep('Moonlit', 1000)),
                    Quadrat = c(rep(levels(Isopod4$Quadrat), 181), 'CA1', 'CA2','CA3','CA4','CA5','CA6','EG1','EG2','EG3','EG4','EG5','EG6','WG1','WG2','WG3','WG4','WG5','WG6'))
simDF$ID <- seq.int(nrow(simDF))
simDF$Pred <- predict(TempPropMAM, data.frame(simDF), type = 'response', allow.new.levels = TRUE)

Figure3b <- ggplot(Isopod4, aes(x = Temp, y = OutProportion, colour = Light_Categories))+
  geom_point(aes(x = Temp, y = OutProportion, colour = Light_Categories))+
  geom_smooth(data = simDF, aes(y=Pred, colour = Light_Categories), level = 0, show.legend = F)+
  theme_bw(base_size = 14)+
  labs(x = "Temperature (\u00B0C)", y = "Proportion of Sea Slaters Outside of Crevices", colour = 'Light Categories') +
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))


###Figure 3
dev.new()
Figure3 <- grid.arrange(Figure3a, Figure3b, nrow = 2)


########################## Swanpool Analysis ##############################
IsopodSP <- Isopod1[Isopod1$Location == 'Swanpool',] #isolate SP data
summary(IsopodSP)
summary(IsopodSP$Light_Categories)
#sample size:
# ALAN  Clouded Moonless  Moonlit 
# 15        8        5        4 



###Total Abundance at different light categories
#plots
plot(IsopodSP$Total~IsopodSP$Light_Categories)

FigureAi<- ggplot(IsopodSP, aes(x=Light_Categories, y=log(Total), colour=Light_Categories))+
  geom_boxplot()+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_blank())+
  labs(x = NULL, y = "log(Total Isopod Count)", colour = 'Light Categories') +
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))+
  geom_text(aes(x = 0.8, y = 3.7,
                label = "Swanpool"),
            stat = "unique",
            size = 5, color = "black")+
  theme(plot.margin=unit(c(0,0,0,1), 'cm'))
FigureAi


###SP Abundance
m1 <- glmer(Total ~ Light_Categories + (1|Quadrat),
               data = IsopodSP, family = poisson)
summary(m1)
dispersion <- 194/27 #overdispersed  => use negative binomial

m1 <- glmer.nb(Total ~ Light_Categories + (1|Quadrat),
          data = IsopodSP)
summary(m1)

m2 <- glmer.nb(Total ~ 1 + (1|Quadrat),
               data = IsopodSP)

anova(m1, m2) #light category not significant p=0.1088


### Proportion
IsopodSP$OutProportion <- IsopodSP$Out/IsopodSP$Total
IsopodSP$OutProportion[is.na(IsopodSP$OutProportion)] <- 0
IsopodSP$ID <- seq.int(nrow(IsopodSP))
IsopodSP$ID <- as.factor(IsopodSP$ID)
summary(IsopodSP)

FigureAii <- ggplot(IsopodSP, aes(x=Light_Categories, y=OutProportion, colour = Light_Categories))+
  geom_boxplot()+
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(angle = 90))+
  labs(x = "Light Categories", y = "Proportion of Sea Slaters Outside of Crevices", colour = 'Light Categories') +
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))+
  theme(plot.margin=unit(c(0,0,0,1), 'cm'))
FigureAii

m1 <- glmer(OutProportion ~ Light_Categories + (1|Quadrat),
          data = IsopodSP, family = binomial, weights = Total)
summary(m1)
dispersion <- 92/24 #overdispersed, include ID 

m2 <- glmer(OutProportion ~ Light_Categories + (1|Quadrat) + (1|ID),
          data = IsopodSP, family = binomial, weights = Total)
summary(m2)

m3 <- glmer(OutProportion ~ 1 + (1|Quadrat) + (1|ID),
          data = IsopodSP, family = binomial, weights = Total)
anova(m2, m3) #light category not significant p=0.09431

#Combined Figures
dev.new()
ggarrange(FigureAi, FigureAii,
                     nrow = 2,
                     labels = c("(i)", "(ii)"),
  common.legend = TRUE,
  legend = "right")
dev.off()


########################## West Gylly Analysis ##############################
IsopodWG <- Isopod1[Isopod1$Location == 'West Gylly',] #isolate WG data
summary(IsopodWG)
summary(IsopodWG$Light_Categories)
# ALAN  Clouded Moonless  Moonlit 
# 0       22        8       10 


###Abundance
m1 <- glmer(Total ~ Light_Categories + (1|Quadrat), 
            data = IsopodWG, family = poisson)
summary(m1)
dispersion <- 225/36 #overdispersed => negative binomial

m2 <- glmer.nb(Total ~ Light_Categories + (1|Quadrat), 
            data = IsopodWG)
summary(m2)

m3 <- glmer.nb(Total ~ 1 + (1|Quadrat),
          data = IsopodWG)
anova(m2, m3) #light category not significant p=0.4578

FigureAiii<- ggplot(IsopodWG, aes(x=Light_Categories, y=log(Total), colour=Light_Categories))+
  geom_boxplot()+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_blank())+
  labs(x = NULL, y = '', colour = 'Light Categories') +
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))+
  geom_text(aes(x = 1, y = 3,
                label = "West Gyllyngvase"),
            stat = "unique",
            size = 5, color = "black")+
  theme(plot.margin=unit(c(0,0,0,1), 'cm'))
FigureAiii

###Proportion
IsopodWG$OutProportion <- IsopodWG$Out/IsopodWG$Total
IsopodWG$OutProportion[is.na(IsopodWG$OutProportion)] <- 0
IsopodWG$ID <- seq.int(nrow(IsopodWG))
IsopodWG$ID <- as.factor(IsopodWG$ID)
summary(IsopodWG)

m1 <- glmer(OutProportion ~ Light_Categories + (1|Quadrat),
          data = IsopodWG, family = binomial, weights = Total)
summary(m1)
dispersion <- 105/36 #overdispersed, include ID

m2 <- glmer(OutProportion ~ Light_Categories + (1|Quadrat) + (1|ID),
          data = IsopodWG, family = binomial, weights = Total)
summary(m2)

m3 <- glmer(OutProportion ~ 1 + (1|Quadrat) + (1|ID),
          data = IsopodWG, family = binomial, weights = Total)
anova(m2, m3) #light category not significant p=0.2213


#Figure Aiv
FigureAiv <- ggplot(IsopodWG, aes(x=Light_Categories, y=OutProportion, colour = Light_Categories))+
  geom_boxplot()+
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(angle = 90))+
  labs(x = "Light Categories", y = '', colour = 'Light Categories') +
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))+
  theme(plot.margin=unit(c(0,0,0,1), 'cm'))
FigureAiv


########################## Castle Analysis ##############################
IsopodCA <- Isopod1[Isopod1$Location == 'Castle',] #isolate Castle data
summary(IsopodCA)
summary(IsopodCA$Light_Categories)
# ALAN  Clouded Moonless  Moonlit 
# 0       20       11        9 

plot(log(IsopodCA$Total)~IsopodCA$Light_Categories)

###Abundance
m1 <- glmer(Total ~ Light_Categories + (1|Quadrat), 
            data = IsopodCA, family = poisson)
summary(m1)
dispersion <- 191/36 #overdispersed => negative binomial

m2 <- glmer.nb(Total ~ Light_Categories + (1|Quadrat), 
               data = IsopodWG)
summary(m2)

m3 <- glmer.nb(Total ~ 1 + (1|Quadrat),
               data = IsopodWG)
anova(m2, m3) #light category not significant p=0.4578

FigureAiii<- ggplot(IsopodWG, aes(x=Light_Categories, y=log(Total), colour=Light_Categories))+
  geom_boxplot()+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_blank())+
  labs(x = NULL, y = '', colour = 'Light Categories') +
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))+
  geom_text(aes(x = 1, y = 3,
                label = "West Gyllyngvase"),
            stat = "unique",
            size = 5, color = "black")+
  theme(plot.margin=unit(c(0,0,0,1), 'cm'))
FigureAiii

###Proportion
IsopodWG$OutProportion <- IsopodWG$Out/IsopodWG$Total
IsopodWG$OutProportion[is.na(IsopodWG$OutProportion)] <- 0
IsopodWG$ID <- seq.int(nrow(IsopodWG))
IsopodWG$ID <- as.factor(IsopodWG$ID)
summary(IsopodWG)

m1 <- glmer(OutProportion ~ Light_Categories + (1|Quadrat),
            data = IsopodWG, family = binomial, weights = Total)
summary(m1)
dispersion <- 105/36 #overdispersed, include ID

m2 <- glmer(OutProportion ~ Light_Categories + (1|Quadrat) + (1|ID),
            data = IsopodWG, family = binomial, weights = Total)
summary(m2)

m3 <- glmer(OutProportion ~ 1 + (1|Quadrat) + (1|ID),
            data = IsopodWG, family = binomial, weights = Total)
anova(m2, m3) #light category not significant p=0.2213


#Figure Aiv
FigureAiv <- ggplot(IsopodWG, aes(x=Light_Categories, y=OutProportion, colour = Light_Categories))+
  geom_boxplot()+
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(angle = 90))+
  labs(x = "Light Categories", y = '', colour = 'Light Categories') +
  scale_color_manual(values = c("ALAN" = "#dc267f", "Moonless" = "#785ef0", "Moonlit" = "#ffb000", "Clouded" = "#648fff"))+
  theme(plot.margin=unit(c(0,0,0,1), 'cm'))
FigureAiv




#Combine Site specific Figures
dev.new()
ggarrange(FigureAi, FigureAiii, FigureAii, FigureAiv,
          nrow = 2,
          ncol = 2,
          labels = c("(i)", "(iii)", '(ii)', '(iv)'),
          common.legend = TRUE,
          legend = "right")
dev.off()
########################## Lunar Rhythm ##############################
install.packages('datawizard')
library(datawizard)
Isopod1$Light_Categories <- fct_collapse(Isopod1$Light_Categories, Natural = c("Moonless","Moonlit")) 
Isopod1$Moon_Phase_Name <- categorize(Isopod1$Moon_Phase, "equal_length", n_groups = 3,
                                 labels = c("New", "Quarter", "Full"))
summary(Isopod1)

ggplot(Isopod1, aes(x=Moon_Phase_Name, y=OutProportion, colour = Light_Categories))+
  geom_boxplot()+
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(angle = 90))+
  scale_color_manual(values = c("ALAN" = "#dc267f", "Natural" = "#ffb000", "Clouded" = "#648fff"))

m1 <- glmer(OutProportion ~ Light_Categories*Moon_Phase_Name + (1|Quadrat),
            data = Isopod1, family = binomial, weights = Total)
summary(m1)
dispersal <- 424/150 # over dispersed 2.827

Isopod1$ID <- seq.int(nrow(Isopod1))
Isopod1$ID <- as.factor(Isopod1$ID)

m1 <- glmer(OutProportion ~ Light_Categories*Moon_Phase_Name + (1|Quadrat)+ (1|ID),
            data = Isopod1, family = binomial, weights = Total) 

m2 <- glmer(OutProportion ~ Light_Categories+Moon_Phase_Name + (1|Quadrat)+ (1|ID),
            data = Isopod1, family = binomial, weights = Total)
anova(m1, m2) #interaction not significant

m3 <- glmer(OutProportion ~ Light_Categories+ (1|Quadrat)+ (1|ID),
            data = Isopod1, family = binomial, weights = Total)
anova(m2, m3) # moon phase is a significant factor p<0.001

m4 <- glmer(OutProportion ~ Moon_Phase_Name + (1|Quadrat)+ (1|ID),
            data = Isopod1, family = binomial, weights = Total)
anova(m2, m4) #light category is significant p=0.01281

summary(m2)





ggplot(Isopod1, aes(x=Moon_Phase_Name, y=Total, colour = Light_Categories))+
  geom_boxplot()+
  theme_bw(base_size = 10)+
  theme(axis.text.y = element_text(angle = 90))+
  scale_color_manual(values = c("ALAN" = "#dc267f", "Natural" = "#ffb000", "Clouded" = "#648fff"))

########################## R Citing ##############################
citation()
version$version.string
packageVersion('dplyr')
packageVersion('forcats')
packageVersion('ggplot2')
packageVersion('ggpubr')
packageVersion('gridExtra')
packageVersion('lme4')
packageVersion('multcomp')
packageVersion('nlme')
packageVersion('tidyr')
packageVersion('tidyverse')





