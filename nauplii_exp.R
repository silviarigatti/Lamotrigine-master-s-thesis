### SCRIPT FOR DATA ANLYSIS OF NAUPLII EXPERIMENT: NAUPLII AND COPEPODITES 
# SURVIVAL FRACTION OVER INITIAL POPULATION, AND 
# NAUPLII MOULTING INTO COPEPODITES


##### Housekeeping

rm(list=ls()) 
#remove ALL objects 
Sys.setenv(LANG = "en")#Let's keep stuff in English
Sys.setlocale("LC_ALL","English")
cat("\014") 
# clear console window prior to new run

#install required packages. 

library(dplyr)
library(ggplot2)
library(RColorBrewer)


options("install.lock"=FALSE)

install.packages("patchwork")
library(patchwork)

install.packages("rlang")
install.packages("ggpmisc")  # Run only once
library(ggpmisc)




##### Set working directory

# R needs you to set the working directory so that it knows from where it 
# should import your datafiles. you can set your working directory from:
# "Session>Set working directory>Choose directory..." and select the folder
# where you have your files.
# another option is to press "ctrl+shift+H" and select the folder 
# (should work on PC and MAC)


##### Import data and start analysis

# import the data to the name of your dataset with taxonomy data

inwd <- "C:/Users/Silvia/OneDrive - University of Gothenburg/THESIS/R/nauplii_exp"
outwd <- "C:/Users/Silvia/OneDrive - University of Gothenburg/THESIS/R/nauplii_exp/output"

setwd(inwd)

data <- read.csv2("nauplii_exp_input.csv")

View(data)


colnames(data)[1] <- 'sample'   #rename 1st column
data$conc <- gsub("Î¼g/L", "µg/L", data$conc)



## try analysis on survival fraction


# Remove rows with NA in survival_fraction
# Convert '#N/A' to actual NA, if needed
#data$survival_fraction[data$survival_fraction == "#N/A"] <- NA

# Then convert to numeric (if it's still character)
data$survival_fraction <- as.numeric(data$survival_fraction)

# Now remove the real NAs
df_clean <- data %>%
  filter(!is.na(survival_fraction))

# Step 1: Group and summarize
avg_data <- df_clean %>%
  group_by(day, conc, organism) %>%
  summarise(
    mean_survival = mean(survival_fraction, na.rm = TRUE),
    .groups = 'drop'
  )

# Summarize: mean and standard deviation
summary_df <- df_clean %>%
  group_by(day, conc, organism) %>%
  summarise(
    mean_survival = mean(survival_fraction, na.rm = TRUE),
    sd = sd(survival_fraction, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()


# Define shapes: dots for copepodites (16), triangles for nauplii (17)
organism_shapes <- c("copepodites" = 17, "nauplii" = 16)

# Plot
survival_plot <- ggplot(summary_df, aes(x = as.numeric(day), y = mean_survival)) +
  geom_line(aes(group = interaction(conc, organism), color = as.factor(conc)), size = 1) +
  geom_point(aes(shape = organism, color = as.factor(conc)),
             size = 3) +
  geom_errorbar(aes(ymin = mean_survival - sd, ymax = mean_survival + sd, color = as.factor(conc)),
                width = 0.3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = organism_shapes) +
  scale_color_brewer(palette = "Set1", name = "Concentration") +
  scale_x_continuous(breaks = 0:14, limits = c(0, 14)) +
  labs(
    x = "Days",
    y = "Mean Survival Fraction over Initial Population",
    shape = "Life stage",
    title = "Survival Fraction Over Time (Mean ± SD)"
  ) +
  theme_light()

survival_plot

#Export the  plot to the working directory as a .png file: 
ggsave(plot = survival_plot,
       filename = paste0(outwd, "/survival_fraction.png"), height = 6, width = 8)







#### ANALYSIS OF MOLTING FRACTION
# Filter out NA values
molting_data <- data %>%
  mutate(molting_fraction_alive = as.numeric(molting_fraction_alive)) %>%
  mutate(molting_fraction_all=as.numeric(molting_fraction_all)) %>%
  filter(!is.na(molting_fraction_alive)) %>%
  filter(!day==1) %>%
  filter(!sample==9)%>%
  filter(!sample==4) %>%
  mutate(molting_fraction_alive=molting_fraction_alive*100) %>%
  mutate(molting_fraction_all=molting_fraction_all*100)





# Plot
ggplot(molting_data, aes(x = as.numeric(day), y = molting_fraction_alive, color = as.factor(conc))) +
  geom_point(size = 2, alpha = 0.7, position = position_dodge(width = 0.2)) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +  # add regression lines without confidence bands
  scale_color_brewer(palette = "Set1", name = "Concentration") +
  scale_x_continuous(breaks = 0:13, limits = c(0, 13)) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20)  # or adjust spacing as needed
  ) +
  labs(
    x = "Day",
    y = "Molting Fraction (Alive)",
    title = "Molting Fraction of Alive Copepodites Over Time"
  ) +
  theme_minimal()

# Plot with equations

plot_molting_alive <- ggplot(molting_data, aes(x = as.numeric(day), y = molting_fraction_alive, color = as.factor(conc))) +
  geom_point(size = 2, alpha = 0.7, position = position_dodge(width = 0.2)) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left", label.y = "top",
    size = 4, show.legend = FALSE
  ) +
  scale_color_brewer(palette = "Set1", name = "Concentration") +
  scale_x_continuous(breaks = 0:13, limits = c(0, 13)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(
    x = "Days",
    y = "Percentage [%]",
    title = "Nauplii Molting into Alive Copepodites Over Time"
  ) +
  theme_light()

ggsave(plot = plot_molting_alive,
       filename = paste0(outwd, "/molting_alive.png"), height = 6, width = 8)



plot_molting_all <- ggplot(molting_data, aes(x = as.numeric(day), y = molting_fraction_all, color = as.factor(conc))) +
  geom_point(size = 2, alpha = 0.7, position = position_dodge(width = 0.2)) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x = "left", label.y = "top",
    size = 4, show.legend = FALSE
  ) +
  scale_color_brewer(palette = "Set1", name = "Concentration") +
  scale_x_continuous(breaks = 0:13, limits = c(0, 13)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(
    x = "Time [d]",
    y = "Percentage [%]",
    title = "Nauplii Molting into Copepodites Over Time"
  ) +
  theme_light()

plot_molting_all
ggsave(plot = plot_molting_all,
       filename = paste0(outwd, "/molting_all.png"), height = 6, width = 8)


### Statistical comparison of the regression slopes between concentrations

molting_data$day <- as.numeric(molting_data$day)
molting_data$conc <- as.factor(molting_data$conc)

## create a model with the regression lines
model_alive <- lm(molting_fraction_alive ~ day * conc, data = molting_data)
summary(model_alive)

anova(model_alive)


par(mfrow = c(2, 2))
plot(model_alive)


install.packages("emmeans")
library(emmeans)

emtr_alive <- emtrends(model_alive, specs = "conc", var = "day")
summary(emtr_alive)

pairs(emtr_alive)


#the same but for all copepodites

model_all <- lm(molting_fraction_all ~ day * conc, data = molting_data)
summary(model_all)

anova(model_all)


par(mfrow = c(2, 2))
plot(model_all)


emtr_all <- emtrends(model_all, specs = "conc", var = "day")
summary(emtr_all)

pairs(emtr_all)










###### ANALYSIS WITH CURVE MODELS
molting_data_1 <- data %>%
  mutate(molting_fraction_alive = as.numeric(molting_fraction_alive)) %>%
  mutate(molting_fraction_all=as.numeric(molting_fraction_all)) %>%
  filter(!is.na(molting_fraction_alive)) %>%
  filter(!sample==9)%>%
  filter(!sample==4) 

molting_data_1 <- subset(molting_data_1, select = -c(status:dead_organisms))
molting_data_1 <- subset(molting_data_1, select = -molting_fraction_all)

library(ggplot2)
library(ggpmisc)  # For stat_poly_eq if you still want it

molting_data_fresh <- molting_data_1 

require(nlstools)
options("install.lock"=FALSE)
install.packages("nls.multstart")  # For robust non-linear fitting
require(nls.multstart)
library(broom)
library(ggpubr)


# Make sure 'conc' is a factor
molting_data_fresh$conc <- as.factor(molting_data_fresh$conc)

# Define Weibull-like function
weibull_model <- function(day, a, b) {
  1 - exp(-(day / a)^b)
}

# Plot raw data
ggplot(molting_data_fresh, aes(x = day, y = molting_fraction_alive, color = conc)) +
  geom_point(alpha = 0.6) +
  labs(title = "Molting Fraction Alive Over Time",
       x = "Day", y = "Molting Fraction (Alive)") +
  theme_minimal()

# Fit Weibull model per concentration group
fits <- molting_data_fresh %>%
  group_by(conc) %>%
  group_map(~nls_multstart(
    molting_fraction_alive ~ 1 - exp(-(day / a)^b),
    data = .x,
    start_lower = c(a = 1, b = 0.1),
    start_upper = c(a = 20, b = 10),
    iter = 500,
    supp_errors = "Y"
  ), .keep = TRUE)


# Tidy parameters
params <- lapply(fits, tidy)
names(params) <- levels(molting_data_fresh$conc)


# Create prediction grid
newdata <- expand.grid(
  day = seq(min(molting_data_fresh$day), max(molting_data_fresh$day), length.out = 100),
  conc = levels(molting_data_fresh$conc)
)

# Predict values
newdata$fit <- NA
for (c in levels(molting_data_fresh$conc)) {
  p <- params[[c]]
  a <- p$estimate[p$term == "a"]
  b <- p$estimate[p$term == "b"]
  newdata$fit[newdata$conc == c] <- weibull_model(newdata$day[newdata$conc == c], a, b)
}

# Plot with fitted curves
ggplot(molting_data_fresh, aes(x = day, y = molting_fraction_alive, color = conc)) +
  geom_point(alpha = 0.5) +
  geom_line(data = newdata, aes(y = fit), size = 1.2) +
  labs(title = "Weibull Fits for Molting Fraction Alive by Concentration",
       x = "Day", y = "Molting Fraction Alive") +
  theme_minimal()

# Model summaries
model_summaries <- lapply(fits, glance)
names(model_summaries) <- levels(molting_data_fresh$conc)
model_summaries
# View parameters
do.call(rbind, params)

# Optional: compare AICs across models
sapply(model_summaries, function(x) x$AIC)


# Load additional required library
library(minpack.lm)

# Ensure factor
molting_data_fresh$conc <- as.factor(molting_data_fresh$conc)

global_model <- nlsLM(
  molting_fraction_alive ~ 1 - exp(-(day / a)^b),
  data = molting_data_fresh,
  start = list(a = 10, b = 2)
)

# Create numeric dummy variables
molting_data_fresh$conc_num <- as.numeric(molting_data_fresh$conc)

# Fit group-specific model with unique a and b per concentration
group_model <- nlsLM(
  molting_fraction_alive ~ 1 - exp(-(day / a_conc[conc_num])^b_conc[conc_num]),
  data = molting_data_fresh,
  start = list(
    a_conc = rep(10, length(levels(molting_data_fresh$conc))),
    b_conc = rep(2, length(levels(molting_data_fresh$conc)))
  ),
  control = nls.lm.control(maxiter = 1000)
)







###### WEIBULL FIT TO THE SCATTERPLOT ####
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(minpack.lm)

# Prepare the dataset (assuming your dataset is in a dataframe called moulting_data_fresh)
# Convert conc as a factor for coloring
molting_data_fresh$conc <- as.factor(molting_data_fresh$conc)

# Define the Weibull model function
weibull_model <- function(day, a, b) {
  return(1 - exp(-(day / a)^b))
}

# Fit the model for each concentration
fits <- molting_data_fresh %>%
  group_by(conc) %>%
  do({
    fit <- nlsLM(molting_fraction_alive ~ weibull_model(day, a, b),
                 data = .,
                 start = list(a = 12, b = 2.5))
    data.frame(summary(fit)$coefficients)
  })

# Check the results
print(fits)




#### add the r squared 


# Assume your dataset is already loaded as 'molting_data_fresh'
# Ensure 'conc' is treated as a factor
molting_data_fresh$conc <- as.factor(molting_data_fresh$conc)

# Create an empty list to hold models and predictions
model_list <- list()
prediction_df <- data.frame()

# Function to calculate R-squared
calc_r_squared <- function(model, data) {
  predicted_values <- fitted(model)
  observed_values <- data$molting_fraction_alive
  ss_res <- sum((observed_values - predicted_values)^2)
  ss_tot <- sum((observed_values - mean(observed_values))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  return(r_squared)
}

# Loop through each concentration and fit Weibull model
r_squared_values <- c()  # store R-squared values for each concentration
for (c in unique(molting_data_fresh$conc)) {
  data_sub <- subset(molting_data_fresh, conc == c)
  
  # Fit the Weibull model: molting_fraction_alive ~ 1 - exp(-(day / a)^b)
  model <- tryCatch({
    nlsLM(molting_fraction_alive ~ 1 - exp(-(day / a)^b),
          data = data_sub,
          start = list(a = 10, b = 2),
          control = nls.lm.control(maxiter = 200))
  }, error = function(e) NULL)
  
  if (!is.null(model)) {
    model_list[[as.character(c)]] <- model
    
    # Calculate R-squared for the model
    r_squared_values[c] <- calc_r_squared(model, data_sub)
    
    # Create sequence of days to predict over
    day_seq <- data.frame(day = seq(min(data_sub$day), max(data_sub$day), length.out = 100))
    preds <- predict(model, newdata = day_seq)
    
    # Combine predictions into a single data frame
    prediction_df <- rbind(prediction_df, 
                           data.frame(day = day_seq$day,
                                      molting_fraction_alive = preds,
                                      conc = c))
  }
}

# Plot

# Define the correct order of concentrations
ordered_concs <- c("0 µg/L", "1 µg/L", "100 µg/L")

# Create the label data frame with matching order
label_positions <- data.frame(
  conc = ordered_concs,
  x = 3.5,
  y = c(0.80, 0.74, 0.68),
  r_squared = round(unlist(r_squared_values[ordered_concs]), 3)
)

# Plot with colored and spaced R² annotations
plot_molt_weibull_R2 <- ggplot(molting_data_fresh, aes(x = day, y = molting_fraction_alive, color = conc)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_line(data = prediction_df, aes(x = day, y = molting_fraction_alive, color = conc), size = 1) +
  scale_color_brewer(palette = "Set1", name = "Concentration") +
  labs(x = "Day", y = "Molting Success", title = "Weibull Model Fit to Alive Molted Copepodites") +
  theme_light() +
  scale_x_continuous(breaks = seq(0, 14, by = 2),
                     minor_breaks = NULL,  # removes minor ticks (like every 0.5)
                     limits = c(0, 14)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))+
  geom_text(data = label_positions,
            aes(x = x, y = y, label = paste("R² =", r_squared), color = conc),
            size = 5,
            show.legend = FALSE)

# display plot
plot_molt_weibull_R2

# save and export plot
ggsave(plot = plot_molt_weibull_R2,
       filename = paste0(outwd, "/plot_molt_weibullfit_r2.png"), height = 6, width = 8)

# Plot the data with the fitted curves #### TO SAVE!###
plot_molt_weibull <- ggplot(molting_data_fresh, aes(x = day, y = molting_fraction_alive, color = conc)) +
  geom_point(size = 2, alpha = 0.6) + # scatterplot
  stat_smooth(method = "nls", formula = y ~ 1 - exp(-(x / a)^b),
              method.args = list(start = c(a = 12, b = 2.5)),
              se = FALSE) + # weibull fits
  theme_light() + 
  scale_x_continuous(breaks = 0:14, limits = c(0, 14)) +
  labs(title = "Alive molted copepodites per day",
       subtitle = "Fit with Weibull curve",
       x = "Day",
       y = "Molting Success") +
  scale_color_brewer(palette = "Dark2", name = "Concentration")

# display plot
plot_molt_weibull

# save and export plot
ggsave(plot = plot_molt_weibull,
       filename = paste0(outwd, "/plot_molt_weibullfit.png"), height = 6, width = 8)


nlsLM_result_Weibull<-nlsLM(molting_fraction_alive~weibull_model(day, a, b),data=molting_data_fresh, 
                            start=list(a=12,b=2.5))

summary(nlsLM_result_Weibull)










####
library(ggplot2)
library(minpack.lm)  # for nlsLM
library(dplyr)

# Assume your dataset is already loaded as 'molting_data_fresh'
# Ensure 'conc' is treated as a factor
molting_data_fresh$conc <- as.factor(molting_data_fresh$conc)

# Create an empty list to hold models and predictions
model_list <- list()
prediction_df <- data.frame()

# Loop through each concentration and fit Weibull model
for (c in unique(molting_data_fresh$conc)) {
  data_sub <- subset(molting_data_fresh, conc == c)
  
  # Fit the Weibull model: molting_fraction_alive ~ 1 - exp(-(day / a)^b)
  model <- tryCatch({
    nlsLM(molting_fraction_alive ~ 1 - exp(-(day / a)^b),
          data = data_sub,
          start = list(a = 10, b = 2),
          control = nls.lm.control(maxiter = 200))
  }, error = function(e) NULL)
  
  if (!is.null(model)) {
    model_list[[as.character(c)]] <- model
    
    # Create sequence of days to predict over
    day_seq <- data.frame(day = seq(min(data_sub$day), max(data_sub$day), length.out = 100))
    preds <- predict(model, newdata = day_seq)
    
    # Combine predictions into a single data frame
    prediction_df <- rbind(prediction_df, 
                           data.frame(day = day_seq$day,
                                      molting_fraction_alive = preds,
                                      conc = c))
  }
}

# Plot
ggplot(molting_data_fresh, aes(x = day, y = molting_fraction_alive, color = conc)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_line(data = prediction_df, aes(x = day, y = molting_fraction_alive, color = conc), size = 1) +
  scale_color_brewer(palette = "Set1", name = "Concentration") +
  labs(x = "Day", y = "Molting Fraction (Alive)", title = "Weibull Model Fit per Concentration") +
  theme_minimal()



par(mfrow = c(1, 3))
for (c in names(model_list)) {
  model <- model_list[[c]]
  plot(fitted(model), residuals(model),
       main = paste("Residuals vs Fitted:", c),
       xlab = "Fitted", ylab = "Residuals")
  abline(h = 0, col = "red")
}
# i see a bad fit...


pseudo_r2 <- function(model, data) {
  rss <- sum(residuals(model)^2)
  tss <- sum((data$molting_fraction_alive - mean(data$molting_fraction_alive))^2)
  1 - rss / tss
}

sapply(names(model_list), function(c) {
  data_sub <- subset(molting_data_fresh, conc == c)
  pseudo_r2(model_list[[c]], data_sub)
})



# Likelihood ratio test (ANOVA) for pairwise model comparison
# First, ensure that the models have been fitted and saved in model_list

# Create a list to store AIC values for each concentration
aic_values <- sapply(model_list, AIC)

# Print AIC values
print(aic_values)

# Compare models using ANOVA (Likelihood Ratio Test)
# Assuming that the models are nested and we want to compare them

# You need to ensure the models are properly fitted for each concentration
# The 'anova()' function compares models if they have the same form
anova_result <- anova(model_list[['0']], model_list[['1']], model_list[['100']])

# Print ANOVA results to check if there's a significant difference between the models
print(anova_result)

# Now let's plot the residuals vs. fitted values for each concentration model
# Plot residuals for each model

# Create a data frame with residuals and fitted values for each model
residuals_fitted_df <- data.frame(
  day = rep(NA, 0),
  residuals = rep(NA, 0),
  fitted_values = rep(NA, 0),
  conc = rep(NA, 0)
)

# Loop through the models to get residuals and fitted values
for (c in names(model_list)) {
  model <- model_list[[c]]
  
  # Get residuals and fitted values
  residuals <- residuals(model)
  fitted_values <- fitted(model)
  conc_values <- rep(c, length(residuals))
  day_values <- molting_data_fresh$day[molting_data_fresh$conc == c]
  
  # Combine the data
  residuals_fitted_df <- rbind(residuals_fitted_df, 
                               data.frame(day = day_values, 
                                          residuals = residuals, 
                                          fitted_values = fitted_values,
                                          conc = conc_values))
}

# Plot residuals vs. fitted values
ggplot(residuals_fitted_df, aes(x = fitted_values, y = residuals, color = conc)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # add zero line
  labs(x = "Fitted Values", y = "Residuals", 
       title = "Residuals vs. Fitted Values for Each Concentration") +
  scale_color_brewer(palette = "Set1", name = "Concentration") +
  theme_minimal()





###################################################
##################################################
#DO THE SAME BUT WITH THREE DIFFERENT WEIBULL MODELS, ONE PER EACH CONC
# I AM USING THE PARAMETERS I FOUND EARLIER TO SET THE A AND B VALUES OF THE WEIBULL FUNCTION


# Load necessary libraries
library(ggplot2)
library(dplyr)

# === 1. Load the data ===
molting_data_fresh$conc <- as.numeric(gsub(" µg/L", "", molting_data_fresh$conc))

# === 2. Weibull-like function ===
weibull_like <- function(x, a, b) {
  1 - exp(-((x / a)^b))
}

# === 3. Fixed parameters per concentration ===
params <- list(
  `0` = c(a = 12.47211, b = 2.323561),
  `1` = c(a = 12.25603, b = 3.1662884),
  `100` = c(a = 11.23793, b = 2.261990)
)

# === 4. Add predictions & residuals ===
molting_data_fresh <- molting_data_fresh %>%
  rowwise() %>%
  mutate(
    pred = weibull_like(day, params[[as.character(conc)]]["a"], params[[as.character(conc)]]["b"]),
    residual = molting_fraction_alive - pred
  ) %>%
  ungroup()

# === 5. Compute R² values ===
install.packages("Metrics")
library(Metrics)
r2_values <- molting_data_fresh %>%
  group_by(conc) %>%
  summarise(
    r2 = cor(molting_fraction_alive, pred)^2
  )

print(r2_values)

# === 6. Pairwise t-tests on residuals ===
res_0 <- molting_data_fresh$residual[molting_data_fresh$conc == 0]
res_1 <- molting_data_fresh$residual[molting_data_fresh$conc == 1]
res_100 <- molting_data_fresh$residual[molting_data_fresh$conc == 100]

t_0_1 <- t.test(res_0, res_1)
t_0_100 <- t.test(res_0, res_100)
t_1_100 <- t.test(res_1, res_100)

print(list(`0_vs_1` = t_0_1$p.value, `0_vs_100` = t_0_100$p.value, `1_vs_100` = t_1_100$p.value))

### do ANOVA on residuals
# Predict values and compute residuals
resid_df <- molting_data_fresh %>%
  rowwise() %>%
  mutate(
    pred = weibull_like(day, params[[conc]]['a'], params[[conc]]['b']),
    resid = molting_fraction_alive - pred
  ) %>%
  ungroup()

# Run ANOVA on residuals
anova_model <- aov(resid ~ conc, data = resid_df)
summary(anova_model)



# === 7. Plot data and Weibull curves ===
# Generate a sequence of x values for curve plotting
curve_data <- do.call(rbind, lapply(c(0, 1, 100), function(c) {
  data.frame(
    day = seq(1, 13, length.out = 300),
    conc = c,
    pred = weibull_like(seq(1, 13, length.out = 300), params[[as.character(c)]]["a"], params[[as.character(c)]]["b"])
  )
}))

curve_data$conc <- factor(curve_data$conc)
molting_data_fresh$conc <- factor(molting_data_fresh$conc)



## write the R2 values found earlier
r2_values <- tibble(
  conc = as.factor(c("0", "1", "100")),
  r2 = c(0.749, 0.800, 0.835)
)

label_positions_1 <- tibble(
  conc = as.factor(c("0", "1", "100")),
  x = 3,  # left side of the plot
  y = c(0.9, 0.84, 0.78)  # stacked so labels don't overlap
)

label_positions_with_r2 <- left_join(label_positions_1, r2_values, by = "conc")



# Plot with colored and spaced R² annotations
plot_molt_weibull_R2_final <- ggplot(molting_data_fresh, aes(x = day, y = molting_fraction_alive, color = conc)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_line(data = curve_data, aes(x = day, y = pred, color = conc), size = 1) +
  scale_color_brewer(palette = "Set1", name = "Concentration") +
  labs(x = "Day", y = "Molting Success", title = "Weibull Model Fit to Alive Molted Copepodites",
       color = "Concentration (µg/L)") +
  theme_light() +
  scale_x_continuous(breaks = seq(0, 14, by = 2),
                     minor_breaks = NULL,  # removes minor ticks (like every 0.5)
                     limits = c(0, 14)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0,1))+
  geom_text(data = label_positions_with_r2,
            aes(x = x, y = y, label = paste0("R² = ", round(r2, 3)), color = conc),
            size = 5,
            show.legend = FALSE)

# display plot
plot_molt_weibull_R2_final

ggsave(plot = plot_molt_weibull_R2_final,
        filename = paste0(outwd, "/plot_molt_weibull_r2_final.png"), height = 6, width = 8)
       

