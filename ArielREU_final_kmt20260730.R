# Impact of Ivermectin treatment and host sex 
 # on non-target coccidia oocyst shedding in dark-eyed juncos; KMT 2026

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
library(glmmTMB)
library(bbmle)
library(janitor)
conflicts_prefer(dplyr::summarize)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::recode)
conflicts_prefer(dplyr::rename)


# read in data 
long=as_tibble(read_csv("coccidia_longdata.csv")) 
wide=as_tibble(read_csv("coccidia_widedata.csv"))

long <- long %>% mutate(across(c(plasmodium_exp, plasmodium_treat, sample, sex, treatment, group), as.factor))
wide <- wide %>% mutate(across(c(plasmodium_exp, plasmodium_treat,  sex, treatment, group), as.factor))
wide$dmass <- wide$mass_postdose - wide$mass_base

# subsets
long.m <- long %>% filter(sex=="Male")
wide.m <- wide %>% filter(sex=="Male")
baseline <- long %>% filter(sample=="baseline")
wide.inoc <- wide %>% filter(treatment==1)
long.inoc <- long %>% filter(treatment==1)

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
treat.colors <- c( "#648FFF","#FFB000")
sex.colors=c("#FE6100", "#785EF0")



# Data summaries, preliminaries ----------------------------------------------------------

# View counts
hist(long$count, breaks = 1000)

# Summarize oocyst count data
countsum <-  long%>% group_by(group, sample) %>% 
  summarize(N=n(), mean=mean(count), sd=sd(count), se=sd/sqrt(N))
print(countsum)

# Do male baseline oocyst counts happen to vary by assigned Ivermectin treatment?
wilcox.test(long.m$count~long.m$treatment) #W=13, p=0.006
  # Yes, so we need to account for individual identity 
  # (mixed model with fixed effect of identity) or use changes in count (dcount)

ggplot(wide.m, aes(x =treatment, y = count_base)) +
  geom_jitter(height=0, width=0.14, stat = "identity", size=3, shape=21) 


# Baseline oocyst counts by sex
ggplot(baseline, aes(x =sex, y = count)) +
  geom_jitter(height=0, width=0.07, stat = "identity", size=3, aes(shape=factor(sex))) +
  labs(       x ="Host sex", y = "Oocyts/gram feces",
              title="Baseline oocyst counts by host sex")+theme(legend.position = "none")+
  scale_x_discrete(breaks=c("0", "1"),
                   labels=c("Male \n(n=17)", "Female \n(n=6)")) +
  scale_shape( solid=FALSE)

# Baseline oocyst counts by sex, excluding Plasmodium-treated males
counts_noexp <- wide %>% filter(plasmodium_exp!="exp")

ggplot(counts_noexp, aes(x =sex, y = count_base)) +
  geom_jitter(height=0, width=0.07, stat = "identity", size=3, aes(shape=factor(sex))) +
  labs(       x ="Host sex", y = "Oocyts/gram feces",
              title="Baseline oocyst counts by host sex")+theme(legend.position = "none")+
  scale_x_discrete(breaks=c("0", "1"),
                   labels=c("Male \n(n=17)", "Female \n(n=6)")) +
  scale_shape( solid=FALSE)


# Baseline and postdose counts by sex and treatment
ggplot(long, aes(x = sample, y = count, color=sex, shape=treatment)) +
  geom_jitter(height=0, width=0.1, size=3, aes(color=factor(sex), shape=factor(treatment)), alpha=0.50) +
  labs(x ="Sample period", y = "Oocyts/gram feces",
       title="Oocysts in ivermectin-treated juncos \nby sampling period and sex",
       color="Sex", shape="Treatment") + 
  scale_color_manual(values = sex.colors, labels= c("Male", "Female"))+
  scale_x_discrete(labels=c("Baseline", "Post-dose"))


# Baseline oocyst counts by sex and Plasmodium treatment-------------------------------------------

### Sex
# Wilcoxon rank sum test: do baseline oocyst counts vary by sex?
wilcox.test(baseline$count~baseline$sex) # yes

# Linear model
hist(wide$count_base, breaks=20) # very right skewed; lm may not be a good fit
base_sex <- lm(count_base~sex, data=wide)
summary(base_sex)
resids(base_sex)
confint(base_sex)

# Wilcoxon rank sum test: do baseline oocyst counts vary by sex (non-Plasmodium birds)?
wilcox.test(counts_noexp$count_base~counts_noexp$sex) #borderline



### Plasmodium treatment
# do male baseline oocyst counts vary by Plasmodium treatment?
hist(wide.m$count_base, breaks=20) #right skewed
wilcox.test(wide.m$count_base~wide.m$plasmodium_treat) #no


### Baseline oocyst counts by sex and Plasmodium treatment
ggplot(baseline, aes(x =sex, y = count, color=plasmodium_treat)) +
  geom_jitter(height=0, width=0.07, stat = "identity", size=3, aes(shape=factor(sex))) +
  labs(title="Baseline coccidia oocyst counts \n by host sex and Plasmodium inoculation",
       x ="Host sex", y = "Oocyts/gram feces", colour="Plasmodium \ntreatment", shape="Sex")+
  scale_x_discrete(breaks=c("0", "1"),
                   labels=c("M (n=17)", "F (n=6)"))+
  scale_color_manual(values=treat.colors, labels=c("No Plasmodium", "Plasmodium", "None")) +
  scale_shape(solid=FALSE, labels=c("Male", "Female"))

### Baseline oocyst counts by sex
ggplot(baseline, aes(x =sex, y = count)) +
  geom_jitter(height=0, width=0.07, stat = "identity", size=4, pch=1) +
  labs(title="Baseline coccidia oocyst burden by host sex",
       x ="Host sex", y = "Oocysts/gram feces")+
  scale_x_discrete(labels=c("M (n=17)", "F (n=6)"))


# Baseline vs. post-treatment oocyst shedding in males (GLMM) --------------------

m_null1 <- glmmTMB(count~1 + (1|band), data=long.m,
                   family=nbinom1, ziformula = ~0)

m_sample1_intx <- glmmTMB(count~sample * treatment + (1|band), data=long.m,
                          family=nbinom1, ziformula = ~0)

m_sample1_add <- glmmTMB(count~sample + treatment + (1|band), data=long.m,
                         family=nbinom1, ziformula = ~0)

m_sample1_add_plas <- glmmTMB(count~sample + treatment + plasmodium_treat+(1|band), data=long.m,
                              family=nbinom1, ziformula = ~0)

AICtab(m_null1, m_sample1_intx, m_sample1_add, m_sample1_add_plas)

resids(m_sample1_add)
summary(m_sample1_add)
Anova(m_sample1_add)
confint(m_sample1_add)

resids(m_sample1_add_plas)
summary(m_sample1_add_plas)
Anova(m_sample1_add_plas) # Plasmodium treatment not significant predictor
confint(m_sample1_add_plas) # Plasmodium treatment not significant predictor

# Baseline to post-treatment oocyst counts in males 
ggplot(wide.m, aes(x =treatment, y = dcount)) +
  geom_jitter(width=0.07, size=4, pch=1)+
  labs( x="Ivermectin treatment",y="Change in oocysts per g feces",
        title="Change in coccidia oocyst burden in male juncos treatment regimen")+
  scale_x_discrete(labels=c("Control", "Ivermectin"))+
  theme(legend.position="none",
        plot.title = element_text(size=16),
        axis.title= element_text(size=14),
        axis.text= element_text(size=12))

# Baseline to post-treatment oocyst counts in males by Plasmodium treatment
ggplot(wide.m, aes(x =plasmodium_treat, y = dcount)) +
  geom_jitter(width=0.07, size=3)+
  labs( x="Plasmodium treatment",y="Change in oocysts per g feces",
        title="Change in oocysts in male juncos \nby treatment regimen")+
  scale_color_manual(values=treat.colors, labels=c("Control", "Ivermectin"))+
  scale_x_discrete(labels=c("Control", "Plasmodium"))


# Baseline vs. post-treatment oocyst shedding in experimentals (GLMM) ------------------------------------------------------------

# glmmTMB documentation: https://journal.r-project.org/articles/RJ-2017-066/RJ-2017-066.pdf 

mf_null1 <- glmmTMB(count~1 + (1|band), data=long.inoc,
                    family=nbinom1, ziformula = ~0)

mf_sample1_intx <- glmmTMB(count~sample * sex + (1|band), data=long.inoc,
                           family=nbinom1, ziformula = ~0)

mf_sample1_add <- glmmTMB(count~sample + sex + (1|band), data=long.inoc,
                          family=nbinom1, ziformula = ~0)


AICtab(mf_null1, mf_sample1_intx, mf_sample1_add)

summary(mf_sample1_add)
resids(mf_sample1_add)
Anova(mf_sample1_add)
confint(mf_sample1_add)

summary(mf_sample1_intx)
resids(mf_sample1_intx)
Anova(mf_sample1_intx, type="III") # sex*sample not significant predictor
confint(mf_sample1_intx) # sex*sample not significant predictor

# Baseline to post-treatment oocyst counts in inoculated juncos 
ggplot(wide.inoc, aes(x =sex, y = dcount)) +
  geom_jitter(width=0.07, size=4, pch=1)+
  labs( x="Sex",y="Change in oocysts per g feces",
        title="Change in coccidia oocyst burden in Ivermectin-treated juncos")+
  scale_x_discrete(labels=c("Males (n=10)", "Females (n=6)"))+
theme(legend.position="none",
      plot.title = element_text(size=16),
      axis.title= element_text(size=14),
      axis.text= element_text(size=12))

# Mass --------------------------------------------------------------------

hist(wide$dmass, breaks=100)
shapiro.test(wide$dmass)

## In Ivermectin-treated juncos, does mass change vary by sex?
dmass.sex <- lm(dmass~sex, data=wide.inoc)
summary(dmass.sex)
#plot(dmass.sex)

ggplot(wide.inoc, aes(x =sex, y = dmass)) +
  geom_jitter(width=0.07, size=4, pch=1)+
  labs( x="Sex",y="Baseline to post-dose mass change (g)",
        title="Mass change by sex")+
  scale_x_discrete(labels=c("Males (n=10)", "Females (n=6)"))+
  theme(        plot.title = element_text(size=16),
                axis.title= element_text(size=14),
                axis.text= element_text(size=12))



## In males, does mass change vary by Ivermectin treatment and/or Plasmodium treatment?
dmass.males <- lm(dmass~treatment+plasmodium_treat, data=wide.m)
summary(dmass.males)
#plot(dmass.males)

ggplot(wide.m, aes(x =treatment, y = dmass, color=plasmodium_treat)) +
  geom_jitter(width=0.07, size=4, pch=1)+
  labs( x="Ivermectin treatment",y="Baseline to post-dose mass change (g)",
        title="Male mass change by treatment", color = "Plasmodium \ntreatment")+
  scale_x_discrete(labels=c("Control (n=7)", "Ivermectin (n=10)"))+
  scale_color_manual(values=c("0"="#785EF0", "1"="#FE6100"), labels=c("Control","Plasmodium"))+
  theme(        plot.title = element_text(size=16),
        axis.title= element_text(size=14),
        axis.text= element_text(size=12))


## In males, does baseline mass vary by Plasmodium treatment?
bmass.males <- lm(mass_base~plasmodium_treat, data=wide.m)
summary(bmass.males) #no
#plot(bmass.males)

ggplot(wide.m, aes(x =plasmodium_treat, y = mass_base)) +
  geom_jitter(width=0.07, size=4, pch=1)+
  labs( x="Plasmodium treatment",y="Baseline mass (g)",
        title="Male baseline mass by Plasmodium treatment")+
  scale_x_discrete(labels=c("Control (n=8)", "Plasmodium (n=9)"))+
theme(legend.position="none",
      plot.title = element_text(size=16),
      axis.title= element_text(size=14),
      axis.text= element_text(size=12))
