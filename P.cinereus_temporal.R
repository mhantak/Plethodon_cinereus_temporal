#Rapid phenotypic change in a polymorphic salamander over 43 years

library(plyr)
library(ggplot2)
library(car)
library(tidyverse)
library(lme4)
library(effects)
library(MuMIn)
library(lmerTest)
library(lattice)
library(gridExtra)

#Data
Pc.data1 <- read.csv("P.cinereus_all_data_11_22_20.csv", header = TRUE, stringsAsFactors = FALSE)
str(Pc.data1)

table(Pc.data1$year)

Pc.data1$Site <- factor(Pc.data1$Site)
Pc.data1$season <- factor(Pc.data1$season)
Pc.data1$Morph <- factor(Pc.data1$Morph)
Pc.data1$Clade <- factor(Pc.data1$Clade)

Pc.data1$year <- as.numeric(Pc.data1$year)
Pc.data1$MAP <- as.numeric(Pc.data1$MAP)

hist(Pc.data1$year) 
hist(Pc.data1$MAP) 
hist(Pc.data1$MAT) 
hist(Pc.data1$SVL_mm) 

##Scale and center continuous predictors (year, MAT and MAP)
Pc.data2 <- transform(Pc.data1, year=scale(year), MAT=scale(MAT), MAP=scale(MAP), TempOfCollectionMonth=scale(TempOfCollectionMonth), PrecOfCollectionMonth=scale(PrecOfCollectionMonth), Temp_Seasonality=scale(Temp_Seasonality), Prec_seasonality=scale(Prec_seasonality))
str(Pc.data2)

hist(Pc.data2$year) 
hist(Pc.data2$MAP) 
hist(Pc.data2$MAT) 

library(PerformanceAnalytics)
test_vars <- c("MAT", "MAP", "year")
Pc.data3 <- Pc.data2[test_vars]
chart.Correlation(Pc.data3, histogram=TRUE, pch=19) 

#SVL cut-off 
Pc.data2_noJUV <- Pc.data2 %>% filter(SVL_mm > 33.999) %>% droplevels()

hist(Pc.data2_noJUV$SVL_mm) 

#####################################################################################################################################################################
##SVL model

svl.clade.f <- lm(SVL_mm ~ Clade + Morph + MAT + year + MAP + season + Morph:year + Morph:MAT + Morph:MAP + 
                    Morph:season, data = Pc.data2_noJUV)

#dredge
options(na.action=na.fail)
svl.clade.models <- dredge(svl.clade.f)
svl.clade.models

#Top model
svl.clade.reduced <- lm(SVL_mm ~ Clade + Morph + MAT + year + season + Morph:year + Morph:MAT, data = Pc.data2_noJUV)

summary(svl.clade.reduced)
vif(svl.clade.reduced)
MuMIn::r.squaredGLMM(svl.clade.reduced)
plot(allEffects(svl.clade.reduced))

s.plot1 <- plot(predictorEffect("year", svl.clade.reduced), lines=list(multiline=TRUE, col=c("red", "black")), 
                confint=list(style="auto"), main=F, lattice=list(key.args=list(x=.01,y=.97,corner=c(0,1),
                columns=1, border=FALSE, cex=.9, cex.title=1)), rug=FALSE, ylab="Snout-vent length (mm)", xlab="Year", ylim=c(39, 44.5),
                axes=list(x=list(cex=1)))

s.plot2 <- plot(predictorEffect("MAT", svl.clade.reduced), lines=list(multiline=TRUE, col=c("red", "black")), 
                confint=list(style="auto"), main=F, lattice=list(key.args=list(x=.76,y=.97,corner=c(0,1),
                columns=1, border=FALSE, cex=1, cex.title=1.2)), rug=FALSE, ylab="Snout-vent length (mm)", xlab="Mean Annual Temperature", ylim=c(39, 44.5),
                axes=list(x=list(cex=1)))

s.plot3 <- plot(effects::effect("season", svl.clade.reduced), main=F, ylab="Snout-vent length (mm)", xlab="Season", rug=FALSE, colors="blue")

s.plot4 <- plot(effects::effect("Clade", svl.clade.reduced), main=F, ylab="Snout-vent length (mm)", xlab="Mitochondrial Clade", rug=FALSE, colors="blue")

grid.arrange(s.plot1, s.plot2, s.plot3, s.plot4, nrow = 2)

################################################################################

###Morph models

### Striped = 1; Unstriped = 0
Pc.data2$Morph2[Pc.data2$Morph == "S"] <- "1"
Pc.data2$Morph2[Pc.data2$Morph == "U"] <- "0"
Pc.data2$Morph2 <- factor(Pc.data2$Morph2) #for logistic models

Morph.global <- glm(Morph2 ~ year + MAT + MAP + season + Clade, data = Pc.data2, family=binomial())

#dredge
options(na.action=na.fail)
morph.clade.models <- dredge(Morph.global)
morph.clade.models

Morph.reduced <- glm(Morph2 ~ year + MAT + MAP + Clade, data = Pc.data2, family=binomial())

summary(Morph.reduced)
vif(Morph.reduced)
plot(allEffects(Morph.reduced))

m.plot1 <- plot(effects::effect("year", Morph.reduced), main=F, ylab="Proportion of striped morphs", xlab="Year", rug=FALSE, colors="forestgreen")
m.plot2 <- plot(effects::effect("MAT", Morph.reduced), main=F, ylab="Proportion of striped morphs", xlab="Mean Annual Temperature", rug=FALSE, colors="forestgreen")
m.plot3 <- plot(effects::effect("MAP", Morph.reduced), main=F, ylab="Proportion of striped morphs", xlab="Mean Annual Precipitation", rug=FALSE, colors="forestgreen")
m.plot4 <- plot(effects::effect("Clade", Morph.reduced), main=F, ylab="Proportion of striped morphs", xlab="Mitochondrial Clade", rug=FALSE, colors="forestgreen")

grid.arrange(m.plot1, m.plot2, m.plot3, m.plot4, nrow = 2)

