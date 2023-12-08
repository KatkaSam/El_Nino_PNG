library(lme4)
library(bbmle)
library(colorspace)
library(emmeans)
library(merTools)

## Read data for individual feeding guilds
Data <- read.csv("C:/Users/ksam/Desktop/Central_Brain/Manuscripts/5_El_Nino_Shifts_RESUBMIT_BIOTROPICA/Guilds_analysispoints.csv")
summary(Data)
Data$Nino <- as.factor(Data$Time)
levels(Data$Nino) <- c("not", "yes", "not")
summary(Data)
## Models

###### Models and figures for frugivores  ####
AbMod.Null <- glmer(Fr_Abund ~ 1 + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Time <- glmer(Fr_Abund ~ Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Nino <- glmer(Fr_Abund ~ Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevS <- glmer(Fr_Abund ~ Elevation + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Elev <- glmer(Fr_Abund ~ poly(Elevation, 2) + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Add <- glmer(Fr_Abund ~ poly(Elevation, 2)+Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.NinoAdd <- glmer(Fr_Abund ~ poly(Elevation, 2)+Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Full <- glmer(Fr_Abund ~ poly(Elevation, 2)*Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.NinoFull <- glmer(Fr_Abund ~ poly(Elevation, 2)*Nino + (1|Elevation:Year), data = Data, family = "poisson")

#AIC table
AICctab(AbMod.Full, AbMod.NinoFull, AbMod.NinoAdd, AbMod.Add, AbMod.ElevS, AbMod.Elev, AbMod.Nino, AbMod.Time, AbMod.Null)
summary(AbMod.NinoFull)
confint(AbMod.NinoFull)

library("qpcR")
AIC <- read.delim ("clipboard")
summary (AIC)
aic<-AIC[,2:2]
aic
akaike.weights(aic)

## Richness
RichMod.Null <- glmer(Fr_SP ~ 1 + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Time <- glmer(Fr_SP ~ Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Nino <- glmer(Fr_SP ~ Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Elev <- glmer(Fr_SP ~ poly(Elevation, 2) + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Add <- glmer(Fr_SP ~ poly(Elevation, 2)+Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.NinoAdd <- glmer(Fr_SP ~ poly(Elevation, 2)+Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Full <- glmer(Fr_SP ~ poly(Elevation, 2)*Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.NinoFull <- glmer(Fr_SP ~ poly(Elevation, 2)*Nino + (1|Elevation:Year), data = Data, family = "poisson")
## newly added models
RichMod.ElevLong <- glmer(Fr_SP ~ Elevation*Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevRes <- glmer(Fr_SP ~ Elevation*Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevAndLong <- glmer(Fr_SP ~ Elevation+Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevS <- glmer(Fr_SP ~ Elevation + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevANDRes <- glmer(Fr_SP ~ Elevation+Nino + (1|Elevation:Year), data = Data, family = "poisson")

#AIC table
AICctab(RichMod.Full, RichMod.NinoFull, RichMod.NinoAdd, RichMod.ElevS, RichMod.Add, RichMod.Elev, RichMod.Nino, 
        RichMod.Time, RichMod.Null, RichMod.ElevLong, RichMod.ElevRes, RichMod.ElevAndLong, RichMod.ElevANDRes)
summary(RichMod.NinoFull)
confint(RichMod.NinoFull)

#Generate the fitted lines for each "season" separately
NewDataPred <- data.frame(Elevation = rep(seq.int(1,4000),3),
                          Time = rep(c("t1", "t2", "t3"), each = 4000),
                          Nino = rep(c("not", "yes", "not"), each = 4000),
                          Year = NA)
NewDataPred$Abundance <- predict(AbMod.NinoFull, newdata = NewDataPred, re.form = NA, type = "response")
NewDataPred$AbundIntervals <- predictInterval(AbMod.NinoFull, newdata = NewDataPred, 
                                              which = "fixed", level = 0.4, stat = c("median"),
                                              n.sims = 20000, type = "probability")
NewDataPred$Species <- predict(RichMod.NinoFull, newdata = NewDataPred, re.form = NA, type = "response")
NewDataPred$SpecIntervals <- predictInterval(RichMod.NinoFull, newdata = NewDataPred, 
                                             which = "fixed", stat = c("median"), level = 0.4, 
                                             n.sims = 20000, type = "probability")
summary(NewDataPred)

# save the data
write.csv(NewDataPred, "C:/Users/ksam/Desktop/Central_Brain/Manuscripts/5_El_Nino_Shifts_RESUBMIT_BIOTROPICA/Revision_autum2021/results_frugivores.csv")
linedata <- read.csv("C:/Users/ksam/Desktop/Central_Brain/Manuscripts/5_El_Nino_Shifts_RESUBMIT_BIOTROPICA/Revision_autum2021/results_frugivores.csv")

summary(Data)
Data$ElevFig <- Data$Elevation + (as.numeric(Data$Sample)*50)-125
summary(Data)

text_size <-  20
library(scales)
show_col(plasma()(20))

library(ggplot2)
plot_04c <- ggplot(data = Data, aes(x=ElevFig, y = Fr_Abund, color = Nino))+
  geom_point(alpha = 0.5, pch = 16, size = 3,
             position = position_jitterdodge(
               dodge.width = 0,
               jitter.width = 0)) +
  ylim(0,40) +
  xlim(0,3000) +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.fit, color = Nino), size = 3) +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.upr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.lwr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Abundance")), size = 20) +
  scale_fill_manual(values = c("#151B8D", "#FCC01E"))+
  scale_color_manual(values = c("#151B8D","#FCC01E"))+
  theme(
    text = element_text(size = 20),
    axis.text=element_text(color="black", size = 20), 
    legend.position = "top") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate('text', x = 450, y = 40, label = 'c) Frugivores', size = 6)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"))
plot_04c  

plot_05c <- ggplot(data = Data, aes(x=ElevFig, y = Fr_SP, color = Nino))+
  geom_point(alpha = 0.5, pch = 16, size = 3,
             position = position_jitterdodge(
               dodge.width = 0,
               jitter.width = 0)) +
  ylim(0,12) +
  xlim(0,3000) +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.fit, color = Nino), size = 3) +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.upr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.lwr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Richness")), size = 20) +
  scale_fill_manual(values = c("#151B8D", "#FCC01E"))+
  scale_color_manual(values = c("#151B8D","#FCC01E"))+
  theme(
    text = element_text(size = 20),
    axis.text=element_text(color="black", size = 20), 
    legend.position = "top") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate('text', x = 450, y = 12, label = 'c) Frugivores', size = 6)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"))
plot_05c 

################################################################################################################################
################################################################################################################################
###### Models and figures for frugivore-insectivores  ####
AbMod.Null <- glmer(FrIn_Abund ~ 1 + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Time <- glmer(FrIn_Abund ~ Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Nino <- glmer(FrIn_Abund ~ Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Elev <- glmer(FrIn_Abund ~ poly(Elevation, 2) + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Add <- glmer(FrIn_Abund ~ poly(Elevation, 2)+Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.NinoAdd <- glmer(FrIn_Abund ~ poly(Elevation, 2)+Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Full <- glmer(FrIn_Abund ~ poly(Elevation, 2)*Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.NinoFull <- glmer(FrIn_Abund ~ poly(Elevation, 2)*Nino + (1|Elevation:Year), data = Data, family = "poisson")
## newly added models
AbMod.ElevLong <- glmer(FrIn_Abund ~ Elevation*Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevRes <- glmer(FrIn_Abund ~ Elevation*Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevAndLong <- glmer(FrIn_Abund ~ Elevation+Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevS <- glmer(FrIn_Abund ~ Elevation + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevANDRes <- glmer(FrIn_Abund ~ Elevation+Nino + (1|Elevation:Year), data = Data, family = "poisson")

#AIC table
AICctab(AbMod.Full, AbMod.NinoFull, AbMod.NinoAdd, AbMod.Add, AbMod.Elev, AbMod.Nino, AbMod.Time, AbMod.Null,
        AbMod.ElevLong, AbMod.ElevRes,AbMod.ElevAndLong, AbMod.ElevANDRes, AbMod.ElevS)
summary(AbMod.NinoFull)
confint(AbMod.NinoFull)

library("qpcR")
AIC <- read.delim ("clipboard")
summary (AIC)
aic<-AIC[,2:2]
aic
akaike.weights(aic)

## Richness
RichMod.Null <- glmer(FrIn_SP ~ 1 + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Time <- glmer(FrIn_SP ~ Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Nino <- glmer(FrIn_SP ~ Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Elev <- glmer(FrIn_SP ~ poly(Elevation, 2) + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Add <- glmer(FrIn_SP ~ poly(Elevation, 2)+Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.NinoAdd <- glmer(FrIn_SP ~ poly(Elevation, 2)+Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Full <- glmer(FrIn_SP ~ poly(Elevation, 2)*Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.NinoFull <- glmer(FrIn_SP ~ poly(Elevation, 2)*Nino + (1|Elevation:Year), data = Data, family = "poisson")
## newly added models
RichMod.ElevLong <- glmer(FrIn_SP ~ Elevation*Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevRes <- glmer(FrIn_SP ~ Elevation*Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevAndLong <- glmer(FrIn_SP ~ Elevation+Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevS <- glmer(FrIn_SP ~ Elevation + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevANDRes <- glmer(FrIn_SP ~ Elevation+Nino + (1|Elevation:Year), data = Data, family = "poisson")

#AIC table
AICctab(RichMod.Full, RichMod.NinoFull, RichMod.NinoAdd, RichMod.Add, RichMod.Elev, RichMod.Nino, RichMod.Time, RichMod.Null,
        RichMod.ElevLong, RichMod.ElevRes,RichMod.ElevAndLong, RichMod.ElevANDRes, RichMod.ElevS)
summary(RichMod.ElevRes)
confint(RichMod.ElevRes)

library("qpcR")
AIC <- read.delim ("clipboard")
summary (AIC)
aic<-AIC[,2:2]
aic
akaike.weights(aic)


#Generate the fitted lines for each "season" separately
NewDataPred <- data.frame(Elevation = rep(seq.int(1,4000),3),
                          Time = rep(c("t1", "t2", "t3"), each = 4000),
                          Nino = rep(c("not", "yes", "not"), each = 4000),
                          Year = NA)
NewDataPred$Abundance <- predict(AbMod.NinoFull, newdata = NewDataPred, re.form = NA, type = "response")
NewDataPred$AbundIntervals <- predictInterval(AbMod.NinoFull, newdata = NewDataPred, 
                                              which = "fixed", level = 0.4, stat = c("median"),
                                              n.sims = 20000, type = "probability")
NewDataPred$Species <- predict(RichMod.NinoFull, newdata = NewDataPred, re.form = NA, type = "response")
NewDataPred$SpecIntervals <- predictInterval(RichMod.NinoFull, newdata = NewDataPred, 
                                             which = "fixed", stat = c("median"), level = 0.4, 
                                             n.sims = 20000, type = "probability")

summary(NewDataPred)

# save the data
write.csv(NewDataPred, "C:/Users/ksam/Desktop/Central_Brain/Manuscripts/5_El_Nino_Shifts_RESUBMIT_BIOTROPICA/Revision_autum2021/results_frugo-insectivores.csv")
linedata <- read.csv("C:/Users/ksam/Desktop/Central_Brain/Manuscripts/5_El_Nino_Shifts_RESUBMIT_BIOTROPICA/Revision_autum2021/results_frugo-insectivores.csv")


summary(Data)
plot_04b <- ggplot(data = Data, aes(x=ElevFig, y = FrIn_Abund, color = Nino))+
  geom_point(alpha = 0.5, pch = 16, size = 3,
             position = position_jitterdodge(
               dodge.width = 0,
               jitter.width = 0)) +
  ylim(0,50) +
  xlim(0,3000) +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.fit, color = Nino), size = 3) +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.upr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.lwr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Abundance")), size = 20) +
  scale_fill_manual(values = c("#151B8D", "#FCC01E"))+
  scale_color_manual(values = c("#151B8D","#FCC01E"))+
  theme(
    text = element_text(size = 20),
    axis.text=element_text(color="black", size = 20), 
    legend.position = "top") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate('text', x = 800, y = 50, label = 'b) Frugivore-insectivores', size = 6)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"))
plot_04b  

plot_05b <- ggplot(data = Data, aes(x=ElevFig, y = FrIn_SP, color = Nino))+
  geom_point(alpha = 0.5, pch = 16, size = 3,
             position = position_jitterdodge(
               dodge.width = 0,
               jitter.width = 0)) +
  ylim(0, 14) +
  xlim(0,3000) +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.fit, color = Nino), size = 3) +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.upr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.lwr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Richness")), size = 20) +
  scale_fill_manual(values = c("#151B8D", "#FCC01E"))+
  scale_color_manual(values = c("#151B8D","#FCC01E"))+
  theme(
    text = element_text(size = 20),
    axis.text=element_text(color="black", size = 20), 
    legend.position = "top") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate('text', x = 800, y = 14, label = 'b) Frugivore-insectivores', size = 6)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"))
plot_05b 

################################################################################################################################
################################################################################################################################
## Abundance of insectiovres
AbMod.Null <- glmer(In_Abund ~ 1 + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Time <- glmer(In_Abund ~ Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Nino <- glmer(In_Abund ~ Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Elev <- glmer(In_Abund ~ poly(Elevation, 2) + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Add <- glmer(In_Abund ~ poly(Elevation, 2)+Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.NinoAdd <- glmer(In_Abund ~ poly(Elevation, 2)+Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Full <- glmer(In_Abund ~ poly(Elevation, 2)*Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.NinoFull <- glmer(In_Abund ~ poly(Elevation, 2)*Nino + (1|Elevation:Year), data = Data, family = "poisson")
## newly added models
AbMod.ElevLong <- glmer(In_Abund ~ Elevation*Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevRes <- glmer(In_Abund ~ Elevation*Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevAndLong <- glmer(In_Abund ~ Elevation+Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevS <- glmer(In_Abund ~ Elevation + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevANDRes <- glmer(In_Abund ~ Elevation+Nino + (1|Elevation:Year), data = Data, family = "poisson")

#AIC table
AICctab(AbMod.Full, AbMod.NinoFull, AbMod.NinoAdd, AbMod.Add, AbMod.Elev, AbMod.Nino, AbMod.Time, AbMod.Null,
        AbMod.ElevLong, AbMod.ElevRes,AbMod.ElevAndLong, AbMod.ElevANDRes, AbMod.ElevS)
summary(AbMod.NinoFull)
confint(AbMod.NinoFull)

library("qpcR")
AIC <- read.delim ("clipboard")
summary (AIC)
aic<-AIC[,2:2]
aic
akaike.weights(aic)


#AIC table
AICctab(AbMod.Full, AbMod.NinoFull, AbMod.NinoAdd, AbMod.Add, AbMod.Elev, AbMod.Nino, AbMod.Time, AbMod.Null,
        AbMod.ElevLong, AbMod.ElevRes,AbMod.ElevAndLong, AbMod.ElevANDRes, AbMod.ElevS)
summary(AbMod.NinoFull)
confint(AbMod.NinoFull)

library("qpcR")
AIC <- read.delim ("clipboard")
summary (AIC)
aic<-AIC[,2:2]
aic
akaike.weights(aic)


## Richness
RichMod.Null <- glmer(In_SP ~ 1 + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Time <- glmer(In_SP ~ Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Nino <- glmer(In_SP ~ Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Elev <- glmer(In_SP ~ poly(Elevation, 2) + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Add <- glmer(In_SP ~ poly(Elevation, 2)+Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.NinoAdd <- glmer(In_SP ~ poly(Elevation, 2)+Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Full <- glmer(In_SP ~ poly(Elevation, 2)*Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.NinoFull <- glmer(In_SP ~ poly(Elevation, 2)*Nino + (1|Elevation:Year), data = Data, family = "poisson")
## newly added models
RichMod.ElevLong <- glmer(In_SP ~ Elevation*Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevRes <- glmer(In_SP ~ Elevation*Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevAndLong <- glmer(In_SP ~ Elevation+Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevS <- glmer(In_SP ~ Elevation + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevANDRes <- glmer(In_SP ~ Elevation+Nino + (1|Elevation:Year), data = Data, family = "poisson")

#AIC table
AICctab(RichMod.Full, RichMod.NinoFull, RichMod.NinoAdd, RichMod.Add, RichMod.Elev, RichMod.Nino, RichMod.Time, RichMod.Null,
        RichMod.ElevLong, RichMod.ElevRes,RichMod.ElevAndLong, RichMod.ElevANDRes, RichMod.ElevS)
summary(RichMod.ElevRes)
confint(RichMod.ElevRes)

library("qpcR")
AIC <- read.delim ("clipboard")
summary (AIC)
aic<-AIC[,2:2]
aic
akaike.weights(aic)


#AIC table
AICctab(RichMod.Full, RichMod.NinoFull, RichMod.NinoAdd, RichMod.Add, RichMod.Elev, RichMod.Nino, RichMod.Time, RichMod.Null)

#Generate the fitted lines for each "season" separately
NewDataPred <- data.frame(Elevation = rep(seq.int(1,4000),3),
                          Time = rep(c("t1", "t2", "t3"), each = 4000),
                          Nino = rep(c("not", "yes", "not"), each = 4000),
                          Year = NA)
NewDataPred$Abundance <- predict(AbMod.NinoFull, newdata = NewDataPred, re.form = NA, type = "response")
NewDataPred$AbundIntervals <- predictInterval(AbMod.NinoFull, newdata = NewDataPred, 
                                              which = "fixed", level = 0.4, stat = c("median"),
                                              n.sims = 20000, type = "probability")
NewDataPred$Species <- predict(RichMod.NinoFull, newdata = NewDataPred, re.form = NA, type = "response")
NewDataPred$SpecIntervals <- predictInterval(RichMod.NinoFull, newdata = NewDataPred, 
                                             which = "fixed", stat = c("median"), level = 0.4, 
                                             n.sims = 20000, type = "probability")

summary(NewDataPred)

# save the data
write.csv(NewDataPred, "C:/Users/ksam/Desktop/Central_Brain/Manuscripts/5_El_Nino_Shifts_RESUBMIT_BIOTROPICA/Revision_autum2021/results_insectivores.csv")
linedata <- read.csv("C:/Users/ksam/Desktop/Central_Brain/Manuscripts/5_El_Nino_Shifts_RESUBMIT_BIOTROPICA/Revision_autum2021/results_insectivores.csv")


summary(Data)
plot_04a <- ggplot(data = Data, aes(x=ElevFig, y = In_Abund, color = Nino))+
  geom_point(alpha = 0.5, pch = 16, size = 3,
             position = position_jitterdodge(
               dodge.width = 0,
               jitter.width = 0)) +
  ylim(0,50) +
  xlim(0,3000) +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.fit, color = Nino), size = 3) +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.upr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  geom_smooth(data = linedata, aes(x=Elevation, y= AbundIntervals.lwr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Abundance")), size = 20) +
  scale_fill_manual(values = c("#151B8D", "#FCC01E"))+
  scale_color_manual(values = c("#151B8D","#FCC01E"))+
  theme(
    text = element_text(size = 20),
    axis.text=element_text(color="black", size = 20), 
    legend.position = "top") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate('text', x = 450, y = 50, label = 'a) Insectivores', size = 6)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"))
plot_04a  

plot_05a <- ggplot(data = Data, aes(x=ElevFig, y = In_SP, color = Nino))+
  geom_point(alpha = 0.5, pch = 16, size = 3,
             position = position_jitterdodge(
               dodge.width = 0,
               jitter.width = 0)) +
  ylim(0, 20) +
  xlim(0,3000) +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.fit, color = Nino), size = 3) +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.upr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  geom_smooth(data = linedata, aes(x=Elevation, y= SpecIntervals.lwr, color = Nino), alpha = 0.3, size = 0.5, linetype="dashed") +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Richness")), size = 20) +
  scale_fill_manual(values = c("#151B8D", "#FCC01E"))+
  scale_color_manual(values = c("#151B8D","#FCC01E"))+
  theme(
    text = element_text(size = 20),
    axis.text=element_text(color="black", size = 20), 
    legend.position = "top") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate('text', x = 450, y = 20, label = 'a) Insectivores', size = 6)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"))
plot_05a 

################################################################################################################################
################################################################################################################################
## Abundance of insectiovrenectarivores
summary(Data)
AbMod.Null <- glmer(INNE_Abund ~ 1 + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Time <- glmer(INNE_Abund ~ Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Nino <- glmer(INNE_Abund ~ Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Elev <- glmer(INNE_Abund ~ poly(Elevation, 2) + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Add <- glmer(INNE_Abund ~ poly(Elevation, 2)+Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.NinoAdd <- glmer(INNE_Abund ~ poly(Elevation, 2)+Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.Full <- glmer(INNE_Abund ~ poly(Elevation, 2)*Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.NinoFull <- glmer(INNE_Abund ~ poly(Elevation, 2)*Nino + (1|Elevation:Year), data = Data, family = "poisson")
## newly added models
AbMod.ElevLong <- glmer(INNE_Abund ~ Elevation*Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevRes <- glmer(INNE_Abund ~ Elevation*Nino + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevAndLong <- glmer(INNE_Abund ~ Elevation+Time + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevS <- glmer(INNE_Abund ~ Elevation + (1|Elevation:Year), data = Data, family = "poisson")
AbMod.ElevANDRes <- glmer(INNE_Abund ~ Elevation+Nino + (1|Elevation:Year), data = Data, family = "poisson")

#AIC table
AICctab(AbMod.Full, AbMod.NinoFull, AbMod.NinoAdd, AbMod.Add, AbMod.Elev, AbMod.Nino, AbMod.Time, AbMod.Null,
        AbMod.ElevLong, AbMod.ElevRes,AbMod.ElevAndLong, AbMod.ElevANDRes, AbMod.ElevS)
summary(AbMod.NinoFull)
confint(AbMod.NinoFull)

library("qpcR")
AIC <- read.delim ("clipboard")
summary (AIC)
aic<-AIC[,2:2]
aic
akaike.weights(aic)

#AIC table
AICctab(AbMod.Full, AbMod.NinoFull, AbMod.NinoAdd, AbMod.Add, AbMod.Elev, AbMod.Nino, AbMod.Time, AbMod.Null)

## Richness
RichMod.Null <- glmer(INNE_SP ~ 1 + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Time <- glmer(INNE_SP ~ Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Nino <- glmer(INNE_SP ~ Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Elev <- glmer(INNE_SP ~ poly(Elevation, 2) + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Add <- glmer(INNE_SP ~ poly(Elevation, 2)+Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.NinoAdd <- glmer(INNE_SP ~ poly(Elevation, 2)+Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.Full <- glmer(INNE_SP ~ poly(Elevation, 2)*Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.NinoFull <- glmer(INNE_SP ~ poly(Elevation, 2)*Nino + (1|Elevation:Year), data = Data, family = "poisson")
## newly added models
RichMod.ElevLong <- glmer(INNE_SP ~ Elevation*Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevRes <- glmer(INNE_SP ~ Elevation*Nino + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevAndLong <- glmer(INNE_SP ~ Elevation+Time + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevS <- glmer(INNE_SP ~ Elevation + (1|Elevation:Year), data = Data, family = "poisson")
RichMod.ElevANDRes <- glmer(INNE_SP ~ Elevation+Nino + (1|Elevation:Year), data = Data, family = "poisson")

#AIC table
AICctab(RichMod.Full, RichMod.NinoFull, RichMod.NinoAdd, RichMod.Add, RichMod.Elev, RichMod.Nino, RichMod.Time, RichMod.Null,
        RichMod.ElevLong, RichMod.ElevRes,RichMod.ElevAndLong, RichMod.ElevANDRes, RichMod.ElevS)
summary(RichMod.ElevRes)
confint(RichMod.ElevRes)

library("qpcR")
AIC <- read.delim ("clipboard")
summary (AIC)
aic<-AIC[,2:2]
aic
akaike.weights(aic)


summary(Data)
plot_04d <- ggplot(data = Data, aes(x=ElevFig, y = INNE_Abund, color = Nino))+
  geom_point(alpha = 0.5, pch = 16, size = 3,
             position = position_jitterdodge(
               dodge.width = 0,
               jitter.width = 0)) +
  ylim(0,40) +
  xlim(0,3000) +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Abundance")), size = 20) +
  scale_fill_manual(values = c("#151B8D", "#FCC01E"))+
  scale_color_manual(values = c("#151B8D","#FCC01E"))+
  theme(
    text = element_text(size = 20),
    axis.text=element_text(color="black", size = 20), 
    legend.position = "top") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate('text', x = 850, y = 40, label = 'd) Nectarivore-insectivores', size = 6)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"))
plot_04d  

plot_05d <- ggplot(data = Data, aes(x=ElevFig, y = INNE_SP, color = Nino))+
  geom_point(alpha = 0.5, pch = 16, size = 3,
             position = position_jitterdodge(
               dodge.width = 0,
               jitter.width = 0)) +
  ylim(0, 10) +
  xlim(0,3000) +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Richness")), size = 20) +
  scale_fill_manual(values = c("#151B8D", "#FCC01E"))+
  scale_color_manual(values = c("#151B8D","#FCC01E"))+
  theme(
    text = element_text(size = 20),
    axis.text=element_text(color="black", size = 20), 
    legend.position = "top") +
  theme_bw(base_size = 14) + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate('text', x = 850, y = 10, label = 'd) Nectarivore-insectivores', size = 6)+
  theme(axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(axis.ticks = element_line(colour = "black", size = 1, linetype = "solid"))
plot_05d 

require(gridExtra)
plot_4<-grid.arrange(plot_04a, plot_04b, plot_04c, plot_04d, ncol=2)
plot_5<-grid.arrange(plot_05a, plot_05b, plot_05c, plot_05d, ncol=2)

