# Author: Christina Jonander (christina.jonander@bioenv.gu.se)

#####
# description of R script
#####

# The script calculates algae consumption in zooplankton bottles, and
# subsets a selected endpoint for statistical testing.
# It also tests the differences between control and
# treatments using an ANOVA and a Dunnett's post hoc test
# & generates a bar graph

#####
# clean up
#####

# remove all objects 
rm(list=ls()) 
# clear console window
cat("\014") 

#####
# Load packages
#####

# tidyverse needs to be run after multcomp & car, if it's already loaded,
# restart R before running this script

require(multcomp)
require(car)
require(tidyverse)

#####
# Import & check data
#####

# import data 
zooplankton_data <- read.csv2("testdata.csv")

# check what the data looks like (shows first 6 lines of dataset)
head(zooplankton_data)

# check how many and what type of variables there are 
str(zooplankton_data)

#####
# calculate consumed algae
#####

# calculate the mean start algae concentration & put it in a new
# dataset named "algae_start"
algae_start<-zooplankton_data %>% 
  filter(type=="start") %>% 
  summarise(mean=mean(algae_0))

# store the mean under the name "algae_start_mean"
algae_start_mean<-algae_start$mean

# calculate difference in algae content between 0h and 48h
# in all algae and zooplankton bottles and put the result in a new
# variable named "algae_diff"
zooplankton_data<-zooplankton_data %>% 
  filter(type!="start") %>%  
  mutate(algae_diff=algae_start_mean-algae_48)

# calculate the means for each algae control treatment and put them in 
# a datasety named "algae_means"
algae_means<-zooplankton_data %>% 
  filter(type=="algae") %>% 
  group_by(conc) %>% 
  summarise(mean=mean(algae_diff)) %>% 
  rename(algae_mean=mean)

# add the means of the algae bottles to the zooplankton dataset
zooplankton_data<-left_join(zooplankton_data, algae_means, by="conc")

# subtract the means of algae bottle treatment from each zooplankton bottle data point
# and name the new column with these numbers "consumed_cells_ml"
zooplankton_data<-zooplankton_data %>% 
  filter(type=="zooplankton") %>% 
  mutate(consumed_cells_ml=algae_diff-algae_mean)

#####
# subset a specific endpoint for further analysis
#####

# select only data from zooplankton bottles, then select the endpoint you want to analyze
# such as egg production (eggprod). then rename the variable "eggprod" to "effect". 
# store the data in the new dataset "data1"
# the level "eggprod" can be changed to any other endpoint, such as "consumed_cells_ml"
data1 <-zooplankton_data %>%
  filter(type =="zooplankton") %>% 
  select(conc, eggprod) %>% 
  rename(effect=eggprod)

#####
# data visualisation
#####

# you can change effect for any other endpoint,
# such as "mortality"

# visualize data, box-whiskers plot
# range = 0 means minimum and maximum are shown in boxplot, box shows 50% of data, 
# black line is median (DO NOT include this graph in the lab report,
# it is only for quickly getting an overview of the data)
boxplot(effect~conc, data = data1, range = 0)

# change the variable "conc" to a factor (needs to be done for the ANOVA to work)
data1$conc <- as.factor(data1$conc)

# check if concentration is still numeric - gives TRUE/FALSE as output
is.numeric(data1$conc)

#####
# ANOVA + post hoc test
#####

# fit a linear model and store it under the name "ANOVA"
ANOVA <- lm(effect~conc, data1)

# compute analysis of variance table for the fitted model
anova(ANOVA)

# make sure 0 is considered the reference when performing multiple-comparisons
data1$conc <- relevel(data1$conc, ref = "0")

# comparisons between control and the other treatments
# specify that "conc" is what you group the data by and compare with a Dunnett's post hoc test
# store the result under the name "Dunnetts"
Dunnetts <- glht(ANOVA, linfct = mcp(conc = "Dunnett"))

# summary of post hoc test
summary(Dunnetts)

#####
# Assumptions - Graphical model checking & tests
#####

# there are some assumptions that the data has to fulfill in order for the 
# results of an ANOVA to be reliable, and two of those are that the residuals of the data are
# normally distributed and that the variances are homogeneous

# the residual plot shows whether the variances are homogeneous. 
# "funnel-shaped" data is something that should be avoided here

# residual plot
# 1 = residual plot
# y axis shows variation 
plot(ANOVA, which = 1)


# we can also test if the variances are homogeneous using a
# levene's test (variances are homogeneous if p>0.05)
leveneTest(effect~conc, data=data1)


# the second plot is a quantile-quantile plot (Q-Q plot), where we can see if the residuals are
# normally distributed. An ideal situation is with as little deviation from the line as possible

# Q-Q plot 
# 2 = QQ plot
plot(ANOVA, which = 2)

# we can also use a shapiro-wilks test to test this assumption (residuals are normal if p>0.05)
shapiro.test(ANOVA$residuals)

#####
# Plot with ggplot
#####

# calculating variatiation:

# summarise means and standard deviation (sd) for the variable "effect" grouped by the variable "conc"
# this is done using code from the package "dplyr" in tidyverse
df_means <- data1 %>% 
  select(conc, effect) %>% 
  group_by(conc) %>% 
  drop_na(effect) %>% 
  summarise_all(funs(mean, sd))

# bar graph
ggplot(data=df_means, aes(x=conc, y=mean))+
  geom_bar(stat="identity", fill="#68B684", position = "stack", width=0.5)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = .1)+
  labs(subtitle="", 
       y="name of y axis", 
       x="name of x axis", 
       title="name of endpoint and chemical")+
  theme(axis.text.x = element_text(size=12, colour = "black"), 
        axis.title.x= element_text(size=12, colour = "black", vjust=-0.2),
        axis.text.y = element_text(size=12, colour="black"), 
        axis.title.y= element_text(size=12, colour = "black", margin=margin(0,20,0,0)),
        panel.border = element_rect(size=0.7, fill=NA),
        legend.title = element_text(size=12, color="white"),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.direction = "vertical", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        plot.title = element_text(size=14))

