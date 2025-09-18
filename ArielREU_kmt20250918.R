# Impact of Ivermectin treatment and host sex 
 # on coccidia oocyst shedding in dark-eyed juncos; KMT July 2025

# Prepare -----------------------------------------------------------------
# clean off workspace
rm(list=ls()) #clear workspace
#dev.off() # close graphical device (if necessary)
cat("\014") # clear console

# read packages
library(dplyr)
library(readxl)
library(tidyverse)
library(lme4)
library(ggplot2)
library(DHARMa)
library(visreg)
library(conflicted)
library(sjPlot)
library(ggpubr)
library(rstatix)
library(MuMIn)
library(lmerTest)
library(gridExtra)
library(MASS)
library(AICcmodavg)
library(FSA)
conflicts_prefer(dplyr::summarize)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::recode)
conflicts_prefer(dplyr::rename)

# set working directory 
setwd("~/Documents/KettlabProjects/ATysver_CoccidiaIvermectin/Final_Files")

# read in data 
data=as_tibble(read_excel("master_juncococcidia.xlsx", sheet=2)) # longform data as changes from baseline
data2=as_tibble(read_excel("master_juncococcidia.xlsx", sheet=1)) # all data, long form

# format data, create dataframes
data$sex <- as.factor(data$sex)
data$treatment <- as.factor(data$treatment)
data$band <- as.factor(data$band)
data$group <- as.factor(data$group)
data$group <- relevel(data$group, 3) # make treatment males the reference level
data$change <- as.factor(data$change)
data$mtreat <- as.factor(data$mtreat)

data <- data %>% mutate (group=recode(group, mtreat="experimental male",
                              ftreat="experimental female",
                              mcontrol="control male"))




baseline <- data2 %>% filter(sample=="base2")
baseline$sex <- as.factor(baseline$sex)
baseline$treatment <- as.factor(baseline$treatment)
baseline$group <- as.factor(baseline$group)
mtreats <- data[c(1,7,8)]
mtreats <- unique(mtreats)
baseline <- merge(baseline, mtreats, by="band")
baseline <- unique(baseline)
baseline$count <- round(baseline$count, 0)

# data subsets
malebase <- baseline %>% filter(sex=="0")
roster <- baseline[c(1:4,9,10)]

pd1 <- data %>% filter(change=="baseline to pd1")
pd1$dcount <- round(pd1$dcount, 0)
pd1$dcount2 <- pd1$dcount+44556
pd1$sex <- as.factor(pd1$sex)
pd1$treatment <- as.factor(pd1$treatment)
pd1$malaria <- as.factor(pd1$malaria)
pd1_males <- pd1 %>% filter(sex=="0")
pd1_treat <- pd1 %>% filter(treatment=="1")

pd2 <- data %>% filter(change=="baseline to pd2")
pd2$dcount <- round(pd2$dcount, 0)
pd2$dcount <- pd2$dcount+1
pd2$dcount2 <- pd2$dcount +7304
pd2$sex <- as.factor(pd2$sex)
pd2$treatment <- as.factor(pd2$treatment)
pd2$malaria <- as.factor(pd2$malaria)
pd2_males <- pd2 %>% filter(sex=="0")
pd2_treat <- pd2 %>% filter(treatment=="1")

# make mass dataframes
sub <- data2[c(1,5,7)]
sub1 <- sub %>% filter(sample=="postdose1")
sub1 <- sub1 %>% rename(mass1=mass)
sub2 <- sub%>% filter(sample=="postdose2")
sub2 <- sub2 %>% rename(mass2=mass)
pdmass <- merge(sub1, sub2, by="band")
basemass <- baseline %>% rename(bmass=mass)
mdata <- merge(basemass, pdmass, by="band")
mdata$massdiff1 <- mdata$mass1-mdata$bmass
mdata$massdiff2 <- mdata$mass2-mdata$bmass
mdata$group <- relevel(mdata$group, ref="mtreat")
pdsub <- data[c(1,5,6)]
pdsub1 <- pdsub %>% filter(change=="baseline to pd1")
pdsub1 <- pdsub1 %>% rename(dcount1=dcount)
pdsub2 <- pdsub %>% filter(change=="baseline to pd2")
pdsub2 <- pdsub2 %>% rename(dcount2=dcount)
pdsub <- merge(pdsub1, pdsub2, by="band")
mdata <- merge(mdata, pdsub, by="band") # wide form mass data
m.mdata <- mdata %>% filter(sex=="0")
t.mdata <- mdata %>% filter(treatment=="1")

# long-form mass data
ch1 <- mdata[c(1,15)]
ch2 <- mdata[c(1,16)]
ch1 <- merge(ch1, roster, by="band")
ch2 <- merge(ch2, roster, by="band")
ch1 <- ch1 %>% rename(massdiff=massdiff1)
ch2 <- ch2 %>% rename(massdiff=massdiff2)
ch1$sample <- "baseline to pd1"
ch2$sample <- "baseline to pd2"
long_msub <- rbind(ch1, ch2)
long_msub.m <- long_msub %>% filter(sex=="0")

pdsub1 <- merge(roster, pdsub1, by="band")
pdsub1 <- pdsub1 %>% rename(dcount=dcount1)
pdsub2 <- merge(roster, pdsub2, by="band")
pdsub2 <- pdsub2 %>% rename(dcount=dcount2)
pdlong <- rbind(pdsub1, pdsub2)
pdlong_inoc <- pdlong %>% filter(treatment=="1")
pdlong.m <- pdlong %>% filter(sex=="0")
long_msub_t <- long_msub%>% filter(treatment=="1")
long_msub_m <- long_msub %>% filter(sex=="0")



# a function for visualizing data and checking distributions
funct_view <- function(a) {
  par(mfrow=c(1,2))
  hist(a)
  boxplot(a)
  shapiro.test(a)}

# a function for checking GLM residuals using package DHARMa 
resids <- function(a) {
  b <- simulateResiduals(a, n=100)
  plot(b)
  testResiduals(b)}

# set color schemes
treat.colors <- c( "#648FFF","#DC267F","#FFB000")
sex.colors=c("#FE6100", "#785EF0")

# Data summaries ----------------------------------------------------------

# note that these are on un-transformed data, without an added 1 to the count value
datasum <- data %>% group_by(group, change) %>% 
  summarize(N=n(), mean_dcount=mean(dcount), sd_dcount=sd(dcount), se_dcount=sd_dcount/sqrt(N))

countsum <-  data2%>% group_by(sample, group) %>% 
  summarize(N=n(), mean=mean(count), sd=sd(count), se=sd/sqrt(N))

# Baseline oocyst counts -------------------------------------------

# all juncos 
wilcox.test(baseline$count~baseline$sex)

# do male baseline oocyst counts happen to vary by assigned Ivermectin treatment?
wilcox.test(malebase$count~malebase$treatment) #W=13, p=0.033
ggplot(malebase, aes(x =treatment, y = count)) +
  geom_jitter(height=0, width=0.14, stat = "identity", size=3, shape=21) 
# yes, so we will focus on changes from baseline

# Baseline-pd1 oocyst count changes -------------------------

# treatment birds only
hist(pd1_treat$dcount)
shapiro.test(pd1_treat$dcount)
wilcox.test(pd1_treat$dcount~pd1_treat$sex) #W = 31, p-value = 0.9578

# male birds only
hist(pd1_males$dcount)
shapiro.test(pd1_males$dcount)
wilcox.test(pd1_males$dcount~pd1_males$treatment) #W = 26, p-value = 0.4173

# by group
pd1_grp <-kruskal_test(dcount~group, data=pd1)
summary(pd1_grp)


# view oocyst changes for all birds, both timeframes 
ggplot(pdlong, aes(x =change, y = dcount, col=treatment)) +
  geom_jitter(width=0.2, size=3,  aes(shape=factor(sex)))+
  labs(shape="Sex", colour="Treatment")+
  scale_color_manual(values=treat.colors, labels=c("Control", "Ivermectin")) +
  scale_shape(labels=c("Male", "Female"), solid=FALSE)+
  labs(       x="Timeframe",
              y="Change in oocysts per g feces")

# Baseline-pd2 oocyst count changes --------------------

# treatment birds
hist(pd2_treat$dcount)
shapiro.test(pd2_treat$dcount)
wilcox.test(pd2_treat$dcount~pd2_treat$sex) #W = 54, p-value = 0.007493

# males only
hist(pd2_males$dcount)
shapiro.test(pd2_males$dcount)
wilcox.test(pd2_males$dcount~pd2_males$treatment) #W = 13, p-value = 0.03301

# by group
kruskal_test(dcount~group, data=pd2)
dunnTest(dcount~group, data=pd2, method="bh")

# Baseline-pd1 mass change -----------------------------

# treatment birds only
t.oodiff <- lm(massdiff1~dcount1, data=t.mdata)
summary(t.oodiff) #p=0.88

t.sex <- lm(massdiff1~sex, data=t.mdata)
summary(t.sex) #p=0.91

t.inx <- lm(massdiff1~dcount1*sex, data=t.mdata)
summary(t.inx) #p=0.92

# visualize
ggplot(t.mdata, aes(x =dcount1, y = massdiff1, color=sex)) +
  geom_point(stat = "identity", size=3, shape=21) +
  labs(title="Baseline to post dose 1 changes in mass and oocysts",
       x="Change in oocyst count from baseline \nto post dose 1 (oocysts/g feces)",
       y="Change in mass from baseline \nto post dose 1 (g)", colour="Sex")+
scale_color_manual(values=sex.colors, labels=c("Male", "Female")) 

# male birds only
m.oodiff <- lm(massdiff1~dcount1, data=m.mdata)
summary(m.oodiff) #p=0.88

m.treat <- lm(massdiff1~treatment, data=m.mdata)
summary(m.treat) #p=0.07

# by group
group_oodiff <- lm(massdiff1~group, data=mdata)
summary(group_oodiff) #p=0.32


# visualize
ggplot(m.mdata, aes(x =dcount1, y = massdiff1, color=treatment)) +
  geom_point(stat = "identity", size=3, shape=21) +
  labs(title="Baseline to post dose 1 changes in \n mass and oocysts (males only)",
       x="Change in oocyst count from baseline \nto post dose 1 (oocysts/g feces)",
       y="Change in mass from baseline \nto post dose 1 (g)")+
  scale_color_manual(values=treat.colors, labels=c("Control", "Ivermectin")) +
  labs(colour="Treatment")

ggplot(m.mdata, aes(x =treatment, y = massdiff1, color=mtreat)) +
  geom_jitter(width=0.1, stat = "identity", size=3, shape=21) +
  scale_color_manual(values=treat.colors, labels=c("Control", "Plasmodium")) +
  labs(title="Baseline to post dose 1 mass change \n by ivermectin and malaria treatment",
       x="Ivermectin treatment", y="Mass change (g)", colour="Plasmodium \ntreatment")+
  scale_x_discrete(labels=c("Control", "Ivermectin"))






# Baseline-pd2 mass change --------------------------------

# treatment birds only
t.oodiff2 <- lm(massdiff2~dcount2, data=t.mdata)
summary(t.oodiff2) #p >0.2

t.sex2 <- lm(massdiff2~sex, data=t.mdata)
summary(t.sex2) #p >0.2

t.inx2 <- lm(massdiff2~dcount1*sex, data=t.mdata)
summary(t.inx2) #p >0.2

# visualize
ggplot(t.mdata, aes(x =dcount2, y = massdiff2, color=sex)) +
  geom_point(stat = "identity", size=3, shape=21) +
  labs(title="Baseline to post dose 1 changes in mass and oocysts",
       x="Change in oocyst count from baseline \nto post dose 2 (oocysts/g feces)",
       y="Change in mass from baseline \nto post dose 2 (g)", colour="Sex")+
  scale_color_manual(values=sex.colors, labels=c("Male", "Female"))



# male birds only
m.oodiff2 <- lm(massdiff2~dcount2, data=m.mdata)
summary(m.oodiff2) #p > 0.2

m.treat2 <- lm(massdiff2~treatment, data=m.mdata)
summary(m.treat2) #p > 0.2

# visualize
ggplot(m.mdata, aes(x =dcount2, y = massdiff2, color=treatment)) +
  geom_point(stat = "identity", size=3, shape=21) +
  labs(title="Baseline to post dose 2 changes in \n mass and oocysts (males only)",
       x="Change in oocyst count from baseline \nto post dose 2 (oocysts/g feces)",
       y="Change in mass from baseline \nto post dose 2 (g)",
       colour="Treatment")+
  scale_color_manual(values=treat.colors, labels=c("Control", "Ivermectin"))

ggplot(m.mdata, aes(x =treatment, y = massdiff2, color=mtreat)) +
  geom_jitter(width=0.1, stat = "identity", size=3, shape=21) +
  labs(title="Baseline to post dose 2 changes in \n mass and oocysts (males only)",
       x="Ivermectin treatment", y="Change in mass (g)", colour="Plasmodium \ntreatment")+
  scale_color_manual(values=treat.colors, labels=c("Control", "Plasmodium"))+
  scale_x_discrete(labels=c("Control", "Ivermectin"))


# by group
oodiff2_group <- lm(massdiff2~group, data=mdata)
summary(oodiff2_group)


# Plasmodium analyses -----------------------------------------------------

# do male baseline oocyst counts vary by Plasmodium treatment?
hist(malebase$count)
wilcox.test(malebase$count~malebase$mtreat) # W=37, p=0.9626

# baseline to post-dose 1 oocyst count changes
wilcox.test(dcount2~mtreat, data=pd1_males)

m.pd1_malariatreat <- lm(dcount2~mtreat*treatment, data=pd1_males)
summary(m.pd1_malariatreat)
Anova(m.pd1_malariatreat)
plot(m.pd1_malariatreat)

# baseline to post-dose 2 oocyst count changes
wilcox.test(dcount2~mtreat, data=pd2_males)

m.pd2_malaria <- glm.nb(dcount2~mtreat, data=pd2_males)
summary(m.pd2_malaria)
resids(m.pd2_malaria)

m.pd2_malariatreat <- glm.nb(dcount2~malaria*treatment, data=pd2_males)
summary(m.pd2_malariatreat)
resids(m.pd2_malariatreat)

# baseline to post-dose 1 mass changes
wilcox.test(massdiff1~mtreat, data=m.mdata)

m.maltreat <- lm(massdiff1~treatment*mtreat, data=m.mdata)
summary(m.maltreat) #p=0.36

m.malaria <- lm(massdiff1~mtreat, data=m.mdata)
summary(m.malaria) #p=0.46

# baseline to post-dose 2 mass changes
wilcox.test(massdiff2~mtreat, data=m.mdata)

m.malaria2 <- lm(massdiff2~mtreat, data=m.mdata)
summary(m.malaria2) #p > 0.2

m.maltreat2<- lm(massdiff2~treatment*mtreat, data=m.mdata)
summary(m.maltreat2) #p > 0.2

# Visuals -----------------------------------------------------------------

# Fig 2: baseline oocyst counts by sex
ggplot(baseline, aes(x =sex, y = count)) +
  geom_jitter(height=0, width=0.07, stat = "identity", size=3, aes(shape=factor(sex))) +
  labs(       x ="Host sex", y = "Oocyts/gram feces")+theme(legend.position = "none")+
  scale_x_discrete(breaks=c("0", "1"),
                   labels=c("Male \n(n=17)", "Female \n(n=6)")) +
  scale_shape( solid=FALSE)


# Fig 3: Oocysts by sampling point, treatment, and sex
oocystsA <- ggplot(pdlong_inoc, aes(x =change, y = dcount)) +
  geom_jitter(width=0.07, size=3,  aes(shape=factor(sex)))+ 
  labs(shape="Sex")+ylim(-50000,110000)+
  scale_shape(labels=c("Male", "Female"), solid=FALSE)+
  labs(       x="Timeframe",
       y="Change in oocysts per g feces")

oocystsB<-ggplot(pdlong.m, aes(x =change, y = dcount, col=treatment)) +
  geom_jitter(width=0.07, size=3, shape=21)+ ylim(-50000,110000)+
  labs( colour="Treatment")+
  scale_color_manual(values=treat.colors, labels=c("Control", "Ivermectin"))+ 
  labs(x="Timeframe",
       y="Change in oocysts per g feces")

  ggarrange (oocystsA, oocystsB, labels=c("A", "B"), ncol=1, nrow=2)


# Fig 4: mass changes by sex, ivermectin treatment, and timeframe
treatmass <- ggplot(long_msub_t, aes(x =sample, y = massdiff)) +
  geom_jitter(width=0.07, size=3,  aes(shape=factor(sex)))+
  labs(shape="Sex")+
  scale_shape(labels=c("Male", "Female"), solid=FALSE)+
  labs(       x="Timeframe",
       y="Change in mass (g)")

malemass <- ggplot(long_msub_m, aes(x =sample, y = massdiff, col=treatment)) +
  geom_jitter(width=0.07, size=3, shape=21)+
  labs(colour="Treatment")+
  scale_color_manual(values=treat.colors, labels=c("Control", "Ivermectin")) +
  labs(      x="Timeframe",
       y="Change in mass (g)")

ggarrange (treatmass, malemass, labels=c("A", "B"), ncol=1, nrow=2)


# Fig S1: baseline oocyst counts by sex, ivermectin treatment, and plasmodium treatment
ggplot(baseline, aes(x =sex, y = count, color=mtreat)) +
  geom_jitter(height=0, width=0.07, stat = "identity", size=3, aes(shape=factor(sex))) +
  labs(title="Baseline coccidia oocyst counts \n by host sex and Plasmodium inoculation",
       x ="Host sex", y = "Oocyts/gram feces", colour="Plasmodium \ntreatment", shape="Sex")+
  scale_x_discrete(breaks=c("0", "1"),
                   labels=c("M (n=17)", "F (n=6)"))+
  scale_color_manual(values=treat.colors, labels=c("No Plasmodium", "Plasmodium", "None")) +
  scale_shape(solid=FALSE, labels=c("Male", "Female"))


# Fig S2: oocyst changes in males by ivermectin and plasmodium treatment
ggplot(pdlong.m, aes(x =change, y = dcount)) +
  geom_jitter(width=0.07, size=3,  aes(shape=factor(mtreat)))+
  labs(shape="Plasmodium \ntreatment", colour="Ivermectin \ntreatment")+
  labs(     x="Timeframe",y="Change in oocysts per g feces") +
  scale_shape_manual(values=c(0,6), labels=c("No Plasmodium", "Plasmodium"))


# Fig S3: mass changes in males by ivermectin treatment and plasmodium treatment
ggplot(long_msub.m, aes(x =sample, y = massdiff)) +
  geom_jitter(width=0.07, size=3,  aes(shape=factor(mtreat)))+
  labs(shape="Plasmodium \ntreatment", 
       x="Timeframe", y="Mass change (g)")+
  scale_shape_manual(values=c(0,6), labels=c("No Plasmodium", "Plasmodium"))
  
