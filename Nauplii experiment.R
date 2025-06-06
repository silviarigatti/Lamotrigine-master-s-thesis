### SCRIPT FOR DATA ANLYSIS OF NAUPLII EXPERIMENT: NAUPLII AND COPEPODITES 
# SURVIVAL FRACTION OVER INITIAL POPULATION, AND 
# NAUPLII MOULTING INTO COPEPODITES

##### Housekeeping
rm(list=ls())  #remove ALL objects 
Sys.setenv(LANG = "en") #Let's keep stuff in English
Sys.setlocale("LC_ALL","English")
cat("\014") # clear console window prior to new run

# install required packages. 
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggpmisc)

##### Set working directory

# R needs you to set the working directory so that it knows from where it 
# should import your datafiles. you can set your working directory from:
# "Session>Set working directory>Choose directory..." and select the folder
# where you have your files.
# another option is to press "ctrl+shift+H" and select the folder 
# (should work on PC and MAC)


##### Import data and start analysis

# import the data to the name of your dataset with data
inwd <- "..."
outwd <- "..."
setwd(inwd)

data <- read.csv2("nauplii_exp_input.csv")

View(data)

colnames(data)[1] <- 'sample'   #rename 1st column
data$conc <- gsub("Î¼g/L", "µg/L", data$conc)


#####
# SURVIVAL FRACTION
#####

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




#####
# MOLTING FRACTION ANALYSIS WITH WEIBULL-LIKE MODEL
#####

molting_data_fresh <- data %>%
  mutate(molting_fraction_alive = as.numeric(molting_fraction_alive),
         molting_fraction_all = as.numeric(molting_fraction_all)  ) %>%
  filter(!is.na(molting_fraction_alive),
         sample != 9,
         sample != 4 ) %>%
  select(-c(status:dead_organisms, molting_fraction_all))


require(nlstools)
require(nls.multstart)
library(broom)
library(ggpubr)
library(ggplot2)
library(minpack.lm)


# Make sure conc is a factor
molting_data_fresh$conc <- as.factor(molting_data_fresh$conc)

# Define the Weibull-like function ===
weibull_model <- function(day, a, b) {
  1 - exp(-((day / a)^b)) 
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


# Write the parameters you just obtained as a and b per each conc
params <- list(
  `0 µg/L` = c(a = 12.47211, b = 2.323561),
  `1 µg/L` = c(a = 12.25603, b = 3.1662884),
  `100 µg/L` = c(a = 11.23793, b = 2.261990)
)

# Add predictions & residuals
molting_data_fresh <- molting_data_fresh %>%
  rowwise() %>%
  mutate(
    pred = weibull_model(day, params[[as.character(conc)]]["a"], params[[as.character(conc)]]["b"]),
    residual = molting_fraction_alive - pred
  ) %>%
  ungroup()

# Calculate R² values
library(Metrics)
r2_values <- molting_data_fresh %>%
  group_by(conc) %>%
  summarise(
    r2 = cor(molting_fraction_alive, pred)^2
  )

print(r2_values)

# Do a pairwise t-tests on residuals
res_0 <- molting_data_fresh$residual[molting_data_fresh$conc == "0 µg/L"]
res_1 <- molting_data_fresh$residual[molting_data_fresh$conc == "1 µg/L"]
res_100 <- molting_data_fresh$residual[molting_data_fresh$conc == "100 µg/L"]

t_0_1 <- t.test(res_0, res_1)
t_0_100 <- t.test(res_0, res_100)
t_1_100 <- t.test(res_1, res_100)

print(list(`0_vs_1` = t_0_1$p.value, `0_vs_100` = t_0_100$p.value, `1_vs_100` = t_1_100$p.value))


str(params)

### do ANOVA on residuals
# Predict values and compute residuals
resid_df <- molting_data_fresh %>%
  rowwise() %>%
  mutate(
    pred = weibull_model(day, params[[conc]]['a'], params[[conc]]['b']),
    resid = molting_fraction_alive - pred
  ) %>%
  ungroup()

# Run ANOVA on residuals
anova_model <- aov(resid ~ conc, data = resid_df)
summary(anova_model)


# Plot data and Weibull curves
# Generate a sequence of x values for curve plotting
curve_data <- do.call(rbind, lapply(c("0 µg/L", "1 µg/L", "100 µg/L"), function(c) {
  data.frame(
    day = seq(1, 13, length.out = 300),
    conc = c,
    pred = weibull_model(seq(1, 13, length.out = 300), params[[as.character(c)]]["a"], params[[as.character(c)]]["b"])
  )
}))


curve_data$conc <- factor(curve_data$conc)
molting_data_fresh$conc <- factor(molting_data_fresh$conc)


## write the R2 values found earlier
r2_values <- tibble(
  conc = as.factor(c("0 µg/L", "1 µg/L", "100 µg/L")),
  r2 = c(0.749, 0.800, 0.835)
)

# Set the position of the labels on the graph
label_positions_1 <- tibble(
  conc = as.factor(c("0 µg/L", "1 µg/L", "100 µg/L")),
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
