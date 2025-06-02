# Author: Silvia Rigatti
# Inspiration from a script by Christina Jonander

#####
# description of R script
#####

# The script calculates algae consumption in zooplankton bottles, and subsets a selected endpoint for statistical testing.
# It also tests the differences between control and treatments using an ANOVA and a Dunnett's post hoc test
# & generates a bar graph

#####
# clean up
#####

# remove all objects 
rm(list=ls()) 
# clear console window
cat("\014") 
Sys.setenv(LANG = "en") #Let's keep stuff in English

#####
# Load packages
#####

# tidyverse needs to be run after multcomp & car, if it's already loaded,
# restart R before running this script

require(multcomp)
require(car)
require(tidyverse)
require(ggplot2)

#####
# Import & check data
#####
## change the folder if needed

inwd <- "C:/Users/Silvia/OneDrive - University of Gothenburg/THESIS/R/3rd_experiment"
outwd <- "C:/Users/Silvia/OneDrive - University of Gothenburg/THESIS/R/3rd_experiment/output"

setwd(inwd)

# import data 
# make sure your file to import is called "testdata"
zooplankton_data <- read.csv2("testdata_3.csv")

# check what the data looks like (shows first 6 lines of dataset)
head(zooplankton_data)

# check how many and what type of variables there are 
str(zooplankton_data)

# To substitute "?" with "μ" in the unit column
zooplankton_data$unit <- gsub("Î¼", "μ", zooplankton_data$unit)
colnames(zooplankton_data)[1] <- 'vial_number'


#####
# FEEDING RATE
#####

# calculate the mean start algae concentration & put it in a new
# dataset named "algae_a_start"
algae_a_start<-zooplankton_data %>% 
  filter(type=="a_start") %>% 
  summarise(mean=mean(algae_0))

# store the mean under the name "algae_a_start_mean"
algae_a_start_mean<-algae_a_start$mean

#Do the same for the ZP bottles
algae_z_start<-zooplankton_data %>% 
  filter(type=="z_start") %>% 
  summarise(mean=mean(algae_0))

# store the mean under the name "algae_start_mean"
algae_z_start_mean<-algae_z_start$mean

# calculate difference in algae content between 0h and 48h
# in all algae and zooplankton bottles and put the result in a new
# variable named "algae_diff"
zooplankton_data1<-zooplankton_data %>% 
  filter(!(type%in% c("a_start", "z_start"))) %>% 
  mutate(algae_diff = case_when(
    type == "algae" ~ algae_a_start_mean - algae_48,   # for type = a_start (algae)
    type == "zooplankton" ~ algae_z_start_mean - algae_48    # for type = z_start (zooplankton)
  ))

# calculate the means for each algae control treatment and put them in 
# a dataset named "algae_means"
algae_means<-zooplankton_data1 %>% 
  filter(type=="algae") %>% 
  group_by(conc) %>% 
  summarise(mean=mean(algae_diff)) %>% 
  rename(algae_mean=mean)

# add the means of the algae bottles to the zooplankton dataset
zooplankton_data1<-left_join(zooplankton_data1, algae_means, by="conc")

# subtract the means of algae bottle treatment from each zooplankton bottle data point
# and name the new column with these numbers "consumed_cells_ml"
zooplankton_data1<-zooplankton_data1 %>% 
  filter(type=="zooplankton") %>% 
  mutate(consumed_cells_ml=algae_diff-algae_mean)

# Divide FR per number of copepods
zooplankton_data1<-zooplankton_data1 %>% 
  filter(type=="zooplankton") %>% 
  mutate(FR_copepod=consumed_cells_ml/copepod_n)

zooplankton_data1$exp_time_zp <- as.numeric(zooplankton_data1$exp_time_zp)

# Divide FR per exposure time
zooplankton_data1<-zooplankton_data1 %>% 
  filter(type=="zooplankton") %>% 
  mutate(FR_copepod_day=FR_copepod/exp_time_zp)

#####
# subset a specific endpoint for further analysis
#####

# select only data from zooplankton bottles, then select the endpoint you want to analyze
# such as egg production per copepod (eggprod_copepod). then rename the variable "eggprod" to "effect". 
# store the data in the new dataset "data1"
# the level "eggprod" can be changed to any other endpoint, such as "consumed_cells_ml"
data1 <-zooplankton_data1 %>%
  select(conc, FR_copepod_day) %>% 
  rename(effect=FR_copepod_day)

#####
# data visualisation
#####

# visualize data, box-whiskers plot
# range = 0 means minimum and maximum are shown in boxplot, box shows 50% of data, 
# black line is median

data1$conc <- factor(data1$conc, levels = conc_order)

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

# Define the desired order of conc levels
conc_order <- c(0, 0.01, 0.1, 1, 10, 100)

# Convert conc to a factor with the specified order in your data frame
df_means$conc <- factor(df_means$conc, levels = conc_order)

# bar graph
FR_copepod <- ggplot(data=df_means, aes(x=conc, y=mean))+
  geom_bar(stat="identity", fill="#68B684", position = "stack", width=0.5)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = .1)+
  labs(subtitle="", 
       y=expression("Consumed cells ml"^{-1}~day^{-1}), 
       x=expression("Conc [μg L"^{-1}~"]"), 
       title="Feeding rate per copepod per day")+
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

#Display the ggplot

FR_copepod

#save the ggplot
ggsave(plot = FR_copepod,
       filename = paste0(outwd, "/FR_copepod_day.png"), height = 6, width = 6)




#####
# TOTAL EGG PRODUCTION
#####

# Sum eggs and nauplii 
zooplankton_data<-zooplankton_data %>% 
  filter(type=="zooplankton") %>% 
  mutate(eggprod=eggs+nauplii)

# Divide eggprod per number of copepods
zooplankton_data<-zooplankton_data %>% 
  filter(type=="zooplankton") %>% 
  mutate(eggprod_copepod=ZPprod/copepod_n)


zooplankton_data$exp_time_eggs<-as.numeric(zooplankton_data$exp_time_eggs)

zooplankton_data<-zooplankton_data %>% 
  filter(type=="zooplankton") %>% 
  mutate(eggprod_copepod_day=eggprod_copepod/exp_time_eggs)

#####
# subset a specific endpoint for further analysis
#####

# select only data from zooplankton bottles, then select the endpoint you want to analyze
# such as egg production per copepod (eggprod_copepod). then rename the variable "eggprod" to "effect". 
# store the data in the new dataset "data1"
# the level "eggprod" can be changed to any other endpoint, such as "consumed_cells_ml"
data1 <-zooplankton_data %>%
  filter(type =="zooplankton") %>% 
  select(conc, eggprod_copepod_day) %>% 
  rename(effect=eggprod_copepod_day)

#####
# data visualisation
#####

# REPEAT LINES 135-198 #

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
eggprod_copepod <- ggplot(data=df_means, aes(x=conc, y=mean))+
  geom_bar(stat="identity", fill="#68B684", position = "stack", width=0.5)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = .1)+
  labs(subtitle="", 
       y=expression("N. of eggs"~day^{-1}), 
       x=expression("Conc [μg L"^{-1}~"]"), 
       title="Total egg production per copepod per day")+
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

eggprod_copepod


# Save the ggplot to a PNG file
ggsave(plot = eggprod_copepod,
       filename = paste0(outwd, "/eggprod_copepod_day.png"), height = 6, width = 6)




#####
# EGG HATCHING #
#####

zooplankton_data$exp_time_eggs<-as.numeric(zooplankton_data$exp_time_eggs)

# Calculate egg hatching
zooplankton_data<-zooplankton_data %>% 
  filter(type=="zooplankton") %>% 
  mutate(egghatch_day=nauplii/(eggs+nauplii)/exp_time_eggs)

#####
# subset a specific endpoint for further analysis
#####

# select only data from zooplankton bottles, then select the endpoint you want to analyze
# such as egg production per copepod (eggprod_copepod). then rename the variable "eggprod" to "effect". 
# store the data in the new dataset "data1"
# the level "eggprod" can be changed to any other endpoint, such as "consumed_cells_ml"
data1 <-zooplankton_data %>%
  filter(type =="zooplankton") %>% 
  select(conc, egghatch_day) %>% 
  rename(effect=egghatch_day)

#####
# data visualisation
#####

# REPEAT LINES 135-198 #

#####
# Plot with ggplot
#####

# calculating variatiation:
df_means <- data1 %>% 
  select(conc, effect) %>% 
  group_by(conc) %>% 
  drop_na(effect) %>% 
  summarise_all(funs(mean, sd))

# bar graph
egghatch_day <- ggplot(data=df_means, aes(x=conc, y=mean))+
  geom_bar(stat="identity", fill="#68B684", position = "stack", width=0.5)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = .1)+
  labs(subtitle="", 
       y=expression("egg hatching success "~day^{-1}), 
       x = expression(Conc~"["*mu*"g "*L^{-1}*"]"),
       title="Hatching success per day")+
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

egghatch_day

# Save the ggplot to a PNG file
ggsave(plot = egghatch_day,
       filename = paste0(outwd, "/egghatch_day.png"), height = 6, width = 6)





#####
# ALGAE PARAMETERS #
#####

# NOW PLOTS ALGAE PARAMETERS WITH START, ALGAE AND 
# ZOOPLANKTON VALUES TO DISPLAY ALL TOGETHER 
#####

#re-upload the original df

zooplankton_data <- read.csv2("testdata_3.csv")

# Specify the order of the concentrations
zooplankton_data$conc <- factor(zooplankton_data$conc, levels = c("start", 
                                                                  0.000e+00,
                                                                  1.0e-02,
                                                                  1.0e-01,
                                                                  1.0e+00,
                                                                  1.0e+01, 
                                                                  1.0e+02))

zooplankton_data$median_size <- as.numeric(zooplankton_data$median_size)

# Summarize data
summary_data <- zooplankton_data %>%
  group_by(conc, type) %>%
  summarize(
    mean_size = mean(median_size),
    sd_size = sd(median_size),
    .groups = "drop"
  )

# Add a column to differentiate "start" from other concentrations
summary_data <- summary_data %>%
  mutate(conc_label = ifelse(type %in% c("a_start", "z_start"), "start", as.character(conc)))

# Set the order for 'conc_label', placing "start" first
summary_data$conc_label <- factor(summary_data$conc_label, 
                                  levels = c("start", sort(unique(as.numeric(as.character(summary_data$conc[summary_data$type != "start"]))))))


summary_data$type <- factor(summary_data$type, levels = c("a_start", "z_start", "algae", "zooplankton"))

# Create the barplot
median_size <- ggplot(summary_data, aes(x = conc_label, y = mean_size, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
  labs(title = "Median algae cell size",
       x = expression("Concentration [µg L"^{-1}~"]"),
       y = "Median Size [µm]",
       fill = "Type"  ) +
  geom_errorbar(aes(ymin = mean_size - sd_size, ymax = mean_size + sd_size), 
                position = position_dodge(width = 0.8), width = 0.25) +
  theme(axis.text.x = element_text(size=12, colour = "black"), 
        axis.title.x= element_text(size=12, colour = "black", vjust=-0.2),
        axis.text.y = element_text(size=12, colour="black"), 
        axis.title.y= element_text(size=12, colour = "black", margin=margin(0,20,0,0)),
        panel.border = element_rect(size=0.7, fill=NA),
        legend.title = element_text(size=12, color="white"),
        legend.text = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        plot.title = element_text(size=14))+
  scale_fill_manual(
    values = c("a_start"="pink2", "z_start" = "mediumaquamarine", "algae" = "lightcoral", "zooplankton" = "lightskyblue"),
    labels=c("Start algae bottles", "Start ZP bottles", "Algae bottles", "Zooplankton bottles"))+
  coord_cartesian(ylim=c(5.5,7.5))

median_size

# Save the ggplot to a PNG file
ggsave(plot = median_size,
       filename = paste0(outwd, "/median_size.png"), height = 6, width = 8)





#####
# RFU per cell #
#####

str(zooplankton_data)

zooplankton_data$RFU<-as.numeric(zooplankton_data$RFU)

# Divide RFU per number of cells
zooplankton_data<-zooplankton_data %>% 
  mutate(RFU_per_cell = ifelse(!is.na(algae_48), RFU / algae_48, ifelse(!is.na(algae_0), RFU / algae_0, NA)))
 
zooplankton_data

zooplankton_data$conc <- factor(zooplankton_data$conc, levels = c("start",0, 0.01, 0.1, 1, 10, 100))
                     
# Summarize data
summary_data <- zooplankton_data %>%
  group_by(conc, type) %>%
  summarize(
    RFU = mean(RFU_per_cell),
    sd_size = sd(RFU_per_cell),
    .groups = "drop"
  )

# Add a column to differentiate "start" from other concentrations
summary_data <- summary_data %>%
  mutate(conc_label = ifelse(type %in% c("a_start", "z_start"), "start", as.character(conc)))

# Set the order for 'conc_label', placing "start" first
summary_data$conc_label <- factor(summary_data$conc_label, 
                                  levels = c("start", sort(unique(as.numeric(as.character(summary_data$conc[summary_data$type != "start"]))))))

summary_data$type <- factor(summary_data$type, levels = c("a_start", "z_start", "algae", "zooplankton"))

# Create the barplot
RFU <- ggplot(summary_data, aes(x = conc_label, y = RFU, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
  labs(title = "RFU per cell",
       x = expression("Concentration [µg L"^{-1}~"]"),
       y = expression("RFU cell"^{-1}),
       fill = "Type"  ) +
  geom_errorbar(aes(ymin = RFU - sd_size, ymax = RFU + sd_size), 
                position = position_dodge(width = 0.8), width = 0.25) +
  theme(axis.text.x = element_text(size=12, colour = "black"), 
        axis.title.x= element_text(size=12, colour = "black", vjust=-0.2),
        axis.text.y = element_text(size=12, colour="black"), 
        axis.title.y= element_text(size=12, colour = "black", margin=margin(0,20,0,0)),
        panel.border = element_rect(size=0.7, fill=NA),
        legend.title = element_text(size=12, color="white"),
        legend.text = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        plot.title = element_text(size=14))+
  scale_fill_manual(
    values = c("a_start"="pink2", "z_start" = "mediumaquamarine", "algae" = "lightcoral", "zooplankton" = "lightskyblue"),
    labels=c("Start algae bottles", "Start ZP bottles", "Algae bottles", "Zooplankton bottles"))+
  scale_y_continuous(breaks = seq(0, 0.6, by=0.2), limits = c(0, 0.6))

RFU

# Set the filename for the output image
RFU_1 <- file.path(outwd, "RFU.png")

# Save the ggplot to a PNG file
ggsave(plot = RFU,
       filename = paste0(outwd, "/RFU per cell.png"), height = 6, width = 8)



### STAT ANALYSIS ON RFU PER CELL

# select only data from zooplankton bottles, then select the endpoint you want to analyze
data1 <-zooplankton_data %>%
  select(conc, type, RFU_per_cell) %>% 
  rename(effect=RFU_per_cell)

data2 <- data1 %>%
  mutate(
    conc_label = case_when(
      type %in% c("a_start", "z_start") & conc == 0 ~ type,  # Change conc to "a_start" or "z_start" if type matches
      TRUE ~ as.character(conc)  # Keep the original conc value for other cases
    )
  )

data3 <-data2 %>%
  select(conc_label, effect) %>% 
  rename(conc=conc_label)

#####
# data visualisation
#####

# you can change effect for any other endpoint,
# such as "mortality"

# visualize data, box-whiskers plot
# range = 0 means minimum and maximum are shown in boxplot, box shows 50% of data, 
# black line is median (DO NOT include this graph in the lab report,
# it is only for quickly getting an overview of the data)
boxplot(effect~conc, data = data3, range = 0)

# change the variable "conc" to a factor (needs to be done for the ANOVA to work)
data3$conc <- as.factor(data3$conc)

# check if concentration is still numeric - gives TRUE/FALSE as output
is.numeric(data3$conc)

#####
# ANOVA + post hoc test
#####

# fit a linear model and store it under the name "ANOVA"
ANOVA <- lm(effect~conc, data3)

# compute analysis of variance table for the fitted model
anova(ANOVA)

# make sure 0 is considered the reference when performing multiple-comparisons
data3$conc <- relevel(data3$conc, ref = "a_start")

# comparisons between control and the other treatments
# specify that "conc" is what you group the data by and compare with a Dunnett's post hoc test
# store the result under the name "Dunnetts"
Dunnetts <- glht(ANOVA, linfct = mcp(conc = "Dunnett"))

# summary of post hoc test
summary(Dunnetts)




#####
# ALGAE CELL COUNT #
#####

# Specify the order of the concentrations
zooplankton_data$conc <- factor(zooplankton_data$conc, levels = c("start", 
                                                                  0.000e+00,
                                                                  1.000e-02,
                                                                  1.0e-01,
                                                                  1e+00,
                                                                  1.0e+01, 
                                                                  1.000e+02))

str(zooplankton_data)

zooplankton_data<-zooplankton_data %>% 
  mutate(cell_count = ifelse(!is.na(algae_48), algae_48, ifelse(!is.na(algae_0), algae_0, NA)))

summary_data <- zooplankton_data %>%
  group_by(conc, type, cell_count)

# Summarize data
summary_data <- zooplankton_data %>%
  group_by(conc, type) %>%
  summarize(
    cellcount = mean(cell_count),
    sd_size = sd(cell_count),
    .groups = "drop"
  )

# Add a column to differentiate "start" from other concentrations
summary_data <- summary_data %>%
  mutate(conc_label = ifelse(type %in% c("a_start", "z_start"), "start", as.character(conc)))

# Set the order for 'conc_label', placing "start" first
summary_data$conc_label <- factor(summary_data$conc_label, 
                                  levels = c("start", sort(unique(as.numeric(as.character(summary_data$conc[summary_data$type != "start"]))))))

summary_data$type <- factor(summary_data$type, levels = c("a_start", "z_start", "algae", "zooplankton"))

# Create the barplot
cell__count <- ggplot(summary_data, aes(x = conc_label, y = cellcount, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
  labs(title = "Algae cell count (cells/ml)",
       x = "Concentration [µg L"^{-1}~"]",
       y = expression("Cells ml"^{-1}),
       fill = "Type"  ) +
  geom_errorbar(aes(ymin = cellcount - sd_size, ymax = cellcount + sd_size), 
                position = position_dodge(width = 0.8), width = 0.25) +
  theme(axis.text.x = element_text(size=12, colour = "black"), 
        axis.title.x= element_text(size=12, colour = "black", vjust=-0.2),
        axis.text.y = element_text(size=12, colour="black"), 
        axis.title.y= element_text(size=12, colour = "black", margin=margin(0,20,0,0)),
        panel.border = element_rect(size=0.7, fill=NA),
        legend.title = element_text(size=12, color="white"),
        legend.text = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        plot.title = element_text(size=14))+
  scale_fill_manual(
    values = c("a_start"="pink2", "z_start" = "mediumaquamarine", "algae" = "lightcoral", "zooplankton" = "lightskyblue"),
    labels=c("Start algae bottles", "Start ZP bottles", "Algae bottles", "Zooplankton bottles"))+
  scale_y_continuous(breaks = seq(0, 20000, by = 2500), limits = c(0, 21000))

cell__count

# Save the ggplot to a PNG file
ggsave(plot = cell__count,
       filename = paste0(outwd, "/cell_count.png"), height = 6, width = 8)




### algae growth in algae bottles

# Specify the order of the concentrations
zooplankton_data$conc <- factor(zooplankton_data$conc, levels = c("start", 
                                                                  0.000e+00,
                                                                  1.000e-02,
                                                                  1.0e-01,
                                                                  1e+00,
                                                                  1.0e+01, 
                                                                  1.000e+02))

# Calculate mean of algae_48 for type == 'a_start'
a_start_mean <- mean(zooplankton_data$algae_48[zooplankton_data$type == "a_start"], na.rm = TRUE)

# Create a new column 'algae_growth' and populate it only for rows where type == 'algae'
zooplankton_data$algae_growth <- NA  # Initialize the column with NA
zooplankton_data$algae_growth[zooplankton_data$type == "algae"] <- 
  zooplankton_data$algae_48[zooplankton_data$type == "algae"] - algae_a_start_mean

str(zooplankton_data)

# Summarize data
summary_data <- zooplankton_data %>%
  filter(type=="algae") %>%
  group_by(conc, type) %>%
  summarize(
    algaegrowth = mean(algae_growth),
    sd_size = sd(algae_growth),
    .groups = "drop"
  )

# Create the barplot
plot_algae_growth <- ggplot(summary_data, aes(x = conc, y = algaegrowth, fill = type)) +
  geom_bar(stat = "identity", fill="#68B684", position = position_dodge(width = 0.7), width = 0.65) +
  labs(title = "Algae growth in algae bottles",
       x = "Conc [µg L"^{-1}~"]",
       y = expression("Cells ml"^{-1}),
       fill = "Type"  ) +
  geom_errorbar(aes(ymin = algaegrowth - sd_size, ymax = algaegrowth + sd_size), 
                position = position_dodge(width = 0.8), width = 0.25) +
  theme(axis.text.x = element_text(size=12, colour = "black"), 
        axis.title.x= element_text(size=12, colour = "black", vjust=-0.2),
        axis.text.y = element_text(size=12, colour="black"), 
        axis.title.y= element_text(size=12, colour = "black", margin=margin(0,20,0,0)),
        panel.border = element_rect(size=0.7, fill=NA),
        legend.title = element_text(size=12, color="white"),
        legend.text = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        plot.title = element_text(size=14))+
  scale_y_continuous(breaks = seq(-2000, 3500, by = 500), limits = c(-2000, 3500))

plot_algae_growth

# Save the ggplot to a PNG file
ggsave(plot = plot_algae_growth,
       filename = paste0(outwd, "/algae_growth.png"), height = 6, width = 6)
