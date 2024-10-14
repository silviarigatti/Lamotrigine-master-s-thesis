
###########
# Housekeeping
###########

##To run a line in R, press ctrl+enter

rm(list=ls()) #remove ALL objects 
cat("\014") # clear console window prior to new run
Sys.setenv(LANG = "en") #Let's keep stuff in English


##########
#import necessary packages
#########
#install.packages("tidyverse")
#install.packages("minpack.lm")
options("install.lock"=FALSE)
install.packages("nlstools")
#install.packages("multcomp")
#install.packages("readxl")


# require and library are similar functions, are both used to recall the 
# packages we have to use
require(tidyverse)
require(minpack.lm)
require(nlstools)
require(multcomp)
require(readxl)
library(ggplot2)

require(multcomp)
require(car)
require(tidyverse)
require(readxl)


#########
#set working directory
#########

inwd <- "C:/Users/Silvia/OneDrive - University of Gothenburg/ECOLOGY/2_ZOOPLANKTON LAB/R"
outwd <- "C:/Users/Silvia/OneDrive - University of Gothenburg/ECOLOGY/2_ZOOPLANKTON LAB/R"

setwd(inwd) # this is to not specify the full path of things I want to import
# when I want to import (and I don't have to specify the full path)







#########
# Establish rawdata to fit

#import data from excel

input_data_algae <- read_excel("input_data_algae.xlsx")
View(input_data_algae)


effectdata <- input_data_algae




ggplot(effectdata, aes(x = factor(conc), y = `cells/ml/day`, group = conc)) +
  geom_boxplot() +
  labs(x = "Concentration", y = "Cells/ml/day", title = "Boxplot of cells/ml/day by concentration") +
  scale_fill_discrete(name = "Concentration")


anova_result <- aov(`cells/ml/day` ~ conc, data = effectdata)

# Print ANOVA table
print(summary(anova_result))





options("install.lock"=FALSE)
install.packages("ggpubr")
library(ggpubr)

# Plotting count algae
count_per_day <- ggplot(effectdata, aes(x = factor(conc), y = `cells/ml/day`, fill = factor(conc))) +
  geom_boxplot() +
  labs(x = "Concentration (%)", y = "Cells/ml/day", title = "Algae cells/ml/day") +
  stat_compare_means(method = "anova") +
  scale_fill_discrete(name = "Concentration")+
  theme(legend.position = "none") 


count_per_day

ggsave(plot = count_per_day,
       filename = paste0(outwd, "/Algae cell count.png"),
       width = 6, height = 6)


# plotting FU

FU_per_day <- ggplot(effectdata, aes(x = factor(conc), y = `FU/day`, fill = factor(conc))) +
  geom_boxplot() +
  labs(x = "Concentration (%)", y = "FU", title = "Algae FU/day") +
  stat_compare_means(method = "anova") +
  scale_fill_discrete(name = "Concentration (%)")+
  theme(legend.position = "none") 

FU_per_day

ggsave(plot = FU_per_day,
       filename = paste0(outwd, "/FU.png"),
       width = 6, height = 6)


# plotting (cells/ml)/FU
count_per_FU <- ggplot(effectdata, aes(x = factor(conc), y = `count/FU`, fill = factor(conc))) +
  geom_boxplot() +
  labs(x = "Concentration", y = "FU", title = "Algae cell count/FU") +
  stat_compare_means(method = "anova") +
  scale_fill_discrete(name = "Concentration")+
  theme(legend.position = "none") 


ggsave(plot = count_per_FU,
       filename = paste0(outwd, "/algae cell count per FU.png"),
       width = 6, height = 6)




### DOSE RESPONSE CURVE


#define the Weibull function

Weibull<-function(tet1, tet2,x){
  1-exp(-exp(tet1+tet2*log10(x)))
}


#########
##define the inverse of the Weibull function
#########
iWeibull<-function(tet1,tet2,x){
  10^((log(-log(1-x))-tet1)/tet2)
}

#########
#define the Logit function
#########
Logit<-function(tet1, tet2,x){
  1/(1+exp(-tet1-tet2*log10(x)))
}

#########
##define the inverse of the Logit function
#########
iLogit<-function(tet1,tet2,x){
  10^(-(log(1/x-1)+tet1)/tet2)
}

#########
#define the Probit function
#########
Probit<-function(tet1, tet2, x){
  pnorm(tet1+tet2*(log10(x)))
}

#########
##define the inverse of the Probit function
#########
iProbit<-function(tet1,tet2,x){
  10^((qnorm(x)-tet1)/tet2)
}



data1$inhibition <- (1-(data1$`yield (cells/ml/46h)`)/(18360.33))



effectdata_without_controls<-subset(data1,data1$conc>0)

# Extract only the controls
controls<-subset(data1,data1$conc==0)

# Establish fake concentration, just for plotting. It's a way to trick R, because 
# it's not possible to plot something with concentration 0 on a log scale graph.
controls$conc <- 0.003


ggplot()+
  geom_point(data=effectdata_without_controls,aes(x=conc,y=inhibition), size = 5)+
  geom_point(data=controls,aes(x=conc,y=inhibition), size = 5, color="red")+
  scale_x_log10("conc")+
  scale_y_continuous("effect")+
  ggtitle("Plot of Rawdata")+
  theme_bw()+
  theme(axis.text.x = element_text(size=18, colour = "black"), 
        axis.title.x= element_text(size=18, colour = "black",margin=margin(20,0,0,0)),
        axis.text.y = element_text(size=18, colour="black"), 
        axis.title.y= element_text(size=18, colour = "black", margin=margin(0,20,0,0)),
        panel.border = element_blank(),
        axis.line.x = element_line(color = 'black',size=1),
        axis.line.y = element_line(color = 'black',size=1),
        legend.title = element_blank(),
        legend.key = element_blank(),
        plot.margin=unit(c(10,10,10,10),"mm") # oben, rechts, unten, links
  ) 



# Weibull
# we look for the best fit with a start
nlsLM_result_Weibull<-nlsLM(inhibition~Weibull(tet1, tet2,conc),data=effectdata_without_controls, 
                            start=list(tet1=1,tet2=1))
# start list may be changed if the function gives errors


# Logit
nlsLM_result_Logit<-nlsLM(inhibition~Logit(tet1, tet2,conc),data=effectdata_without_controls, 
                          start=list(tet1=1,tet2=1))

# Probit
nlsLM_result_Probit<-nlsLM(inhibition~Probit(tet1, tet2,conc),data=effectdata_without_controls, 
                           start=list(tet1=1,tet2=1))

#check whether fit was successful and print
#parameter estimates for tet1 and tet2
summary(nlsLM_result_Weibull)
summary(nlsLM_result_Logit)
summary(nlsLM_result_Probit)




# Step 1: generate a sequence of conc/effect data based on the fitted model

# Step 1a: generate a dataframe (named "sequence") that contains a sequence of data between 
# minimum and maximum of the concentration range  that was tested

minimum_conc<-min(effectdata_without_controls$conc)
maximum_conc<-max(effectdata_without_controls$conc)
sequence=data.frame(conc=(seq(minimum_conc,maximum_conc,length.out=100)))

# Step 1b: generate a dataframe (named "plotdata_fit") that contains a sequence
# of effectdata, calculated on the basis of the results from the fitting process
# using the dataframe "sequence" as input-concentrations
plotdata_fit_Weibull<-as.data.frame(predict(nlsLM_result_Weibull,newdata=sequence))
plotdata_fit_Logit<-as.data.frame(predict(nlsLM_result_Logit,newdata=sequence))
plotdata_fit_Probit<-as.data.frame(predict(nlsLM_result_Probit,newdata=sequence))

# Step 1c: Merge dataframes "sequence" and "plotdata_fit" together into the 
# dataframe "plotdata_fit". Name the variables in the dataframe 
# apprpriately "conc" and "effect"
plotdata_fit<-cbind(sequence, plotdata_fit_Weibull)
colnames(plotdata_fit)[1]<-"conc"
colnames(plotdata_fit)[2]<-"effect_Weibull"

plotdata_fit<-cbind(plotdata_fit, plotdata_fit_Logit)
colnames(plotdata_fit)[3]<-"effect_Logit"

plotdata_fit<-cbind(sequence, plotdata_fit_Probit)
colnames(plotdata_fit)[2]<-"effect_Probit"


# Step 2: plot the curve (available in the form of 100 conc/effect datapairs in
# the dataframe "plotdata_fit" and the raw data (available as conc/effect datapairs
# in the dataframe "effectdata_without_controls")

ggplot()+
  geom_point(data=effectdata_without_controls,aes(x=conc,y=inhibition), size = 5)+
  geom_point(data=controls,aes(x=conc,y=inhibition), size = 5, color="red")+
  #geom_line(data=plotdata_fit, aes(x=conc, y=effect_Weibull), colour="blue",size=1.5, linetype="solid")+
  #geom_line(data=plotdata_fit, aes(x=conc, y=effect_Logit), colour="green",size=1.5, linetype="solid")+
  geom_line(data=plotdata_fit, aes(x=conc, y=effect_Probit), colour="green4",size=1.5, linetype="solid")+
  scale_x_log10("Concentration (%)")+
  scale_y_continuous("Inhibition")+
  ggtitle("Inhibition Algal Growth, Fit to Probit model")+
  theme_bw()+
  theme(plot.title = element_text(size=20),
        axis.text.x = element_text(size=18, colour = "black"), 
        axis.title.x= element_text(size=18, colour = "black",margin=margin(20,0,0,0)),
        axis.text.y = element_text(size=18, colour="black"), 
        axis.title.y= element_text(size=18, colour = "black", margin=margin(0,20,0,0)),
        panel.border = element_blank(),
        axis.line.x = element_line(color = 'black',size=1),
        axis.line.y = element_line(color = 'black',size=1),
        legend.title = element_blank(),
        legend.key = element_blank(),
        plot.margin=unit(c(10,10,10,10),"mm") # oben, rechts, unten, links
  )




# the following 3 lines are how to produce a Dunnet test

data1$group<-as.factor(data1$conc)
fit <- aov(inhibition ~ group, data1)
Dunnet <- glht(fit, linfct=mcp(group="Dunnett"))

summary(Dunnet)
# look at the summary. on the left, the first column, the values of the concentration
# on the right, next to the lines, when we see for the first time the * (one, or two) 
# means that that concentration is the LOEC (lowest observed effect concentration)
# The concentration value before the LOEC is the NOEC. 