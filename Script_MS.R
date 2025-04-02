

# DATA AND PACKAGES ####
#______________________________________________________________________________
#                  Meta data                                               ====
#______________________________________________________________________________
# Script by: Ellinor Jakobsson, 2024-2025
# Script for plotting all figures and doing all analyses in manuscript:
#"Timing of nutrient input pulses during ice cover regulates lake spring plankton dynamics"

#______________________________________________________________________________
#                  Load data and packages                                  ====
#______________________________________________________________________________
library(rstudioapi)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(lubridate)
library(purrr)
library(conflicted)
library(segmented)
library(ggh4x)
library(forestmangr)
library(tidyverse)
library(broom)
library(nlstools)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::select)
conflicts_prefer(lubridate::yday)
conflicts_prefer(dplyr::summarise)


dir <- "C:/Users/ellja901/Documents/PhD/PhD documents/Projects/Mesocosm-experiment/Mesocosm experiment/Mesocosm-experiment/Cleaned_data/"
fig_dir <- "C:/Users/ellja901/Documents/PhD/PhD documents/Projects/Mesocosm-experiment/Mesocosm experiment/Mesocosm-experiment/Figures/"
setwd(dir)
#Load chl-a data
Chla_data_QC <- readxl::read_xlsx(paste0(dir, "Chla_data_QC_new.xlsx"), sheet = 2)

#Load nutrient data
CNP_data_QC <- readxl::read_xlsx(paste0(dir, "CNP_data_QC.xlsx"), sheet = 2)

#Load snow and ice data
Snow_ice_data_QC <- readxl::read_xlsx(paste0(dir, "Snow_ice_data_QC.xlsx"), sheet = 2)

#Load bacterial abundance, complexity and size data
FCM_data_QC <- readxl::read_xlsx(paste0(dir, "FCM_data_QC.xlsx"), sheet = 2) %>%
  filter(ID != "Lake") %>%
  mutate(Abundance_Bacteria_events_mL_V  = Abundance_Bacteria_events_uL_V*1000)

#Load DF data
DF_data_QC <- readxl::read_xlsx(paste0(dir, "DF_data_QC.xlsx"), sheet = 2)

#Load sensor data
Sensor_data <- readxl::read_xlsx(paste0(dir, "Sensor_data_QC.xlsx"), sheet = 2)

#Load biovolume data
Final_biovolume_data <- read_xlsx(paste0(dir,"Phytoplankton_data.xlsx"))

#load zooplankton data
Zooplankton_data <- readxl::read_xlsx(paste0(dir, "Zooplankton_data.xlsx"))



# MODELS ####
#______________________________________________________________________________
#                  Bacterial abundance - GAM Mixed effect global model     ==== 
#______________________________________________________________________________
#Check for linearity
# Load necessary libraries
library(mgcv)
#library(lmtest)
FCM_data_QC$DOY <- lubridate::yday(FCM_data_QC$Date)
# 1. Plot the data to visually inspect non-linear trends
ggplot(FCM_data_QC, aes(x = DOY, y = Abundance_Bacteria_events_mL_V, color = ID)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# 2. Fit linear, polynomial, and non-linear models
# Linear model
lm_model <- lm(Abundance_Bacteria_events_mL_V ~ DOY, data = FCM_data_QC)
summary(lm_model)

# Polynomial model
poly_model <- lm(Abundance_Bacteria_events_mL_V ~ poly(DOY, 2), data = FCM_data_QC)
summary(poly_model)

# Non-linear model using logistic growth
nls_model <- nls(Abundance_Bacteria_events_mL_V ~ SSlogis(DOY, Asym, xmid, scal), data = FCM_data_QC)
summary(nls_model)

# Compare models using AIC
AIC(lm_model, poly_model, nls_model)
FCM_data_QC <- FCM_data_QC %>% drop_na(Abundance_Bacteria_events_mL_V)
# 3. Fit a Generalized Additive Model (GAM)
gam_model <- gam(Abundance_Bacteria_events_mL_V ~ s(DOY, bs = "cs"), data = FCM_data_QC)
summary(gam_model)

# Plot the fitted GAM
plot(gam_model, residuals = TRUE)

# 4. Examine residuals from the linear model
plot(lm_model$residuals ~ FCM_data_QC$DOY)
abline(h = 0, col = "red")

# 5. Statistical tests
FCM_data_QC$ID <- as.factor(FCM_data_QC$ID)
FCM_data_QC$ID <- relevel(FCM_data_QC$ID, ref = "Control")
model_gamm_simple <- gamm(Abundance_Bacteria_events_mL_V ~ ID + s(DOY), 
                          random = list(Mesocosm = ~1), 
                          data = FCM_data_QC)

model_gamm_interaction <- gamm(Abundance_Bacteria_events_mL_V ~ s(DOY, by = ID), 
                               random = list(Mesocosm = ~1), 
                               data = FCM_data_QC)

summary(model_gamm_simple$gam)
summary(model_gamm_interaction$gam)

#Calculate descriptives
FCM_data_QC %>% group_by(ID) %>% reframe(mean(Abundance_Bacteria_events_mL_V),
                                         sd(Abundance_Bacteria_events_mL_V))
#______________________________________________________________________________
#                   Chl-a - GAM Mixed effect global model                  ==== 
#______________________________________________________________________________

#Check for linearity
Chla_data_QC$ID <- as.factor(Chla_data_QC$ID)
Chla_data_QC$DOY <- yday(Chla_data_QC$Date)
Chla_data_QC <- Chla_data_QC %>% filter(ID != "Lake")

# 1. Plot the data to visually inspect non-linear trends
ggplot(Chla_data_QC, aes(x = DOY, y = Chla, color = ID)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# 2. Fit linear, polynomial, and non-linear models
# Linear model
lm_model <- lm(Chla ~ DOY, data = Chla_data_QC)
summary(lm_model)

# Polynomial model
poly_model <- lm(Chla ~ poly(DOY, 2), data = Chla_data_QC)
summary(poly_model)

# Non-linear model using logistic growth
nls_model <- nls(Chla ~ SSlogis(DOY, Asym, xmid, scal), data = Chla_data_QC)
summary(nls_model)

# Compare models using AIC
AIC(lm_model, poly_model, nls_model)

Chla_data_QC <- Chla_data_QC %>% drop_na(Chla)
# 3. Fit a Generalized Additive Model (GAM)
gam_model <- gam(Chla ~ s(DOY, bs = "cs"), data = Chla_data_QC)
summary(gam_model)

# Plot the fitted GAM
plot(gam_model, residuals = TRUE)

# 4. Examine residuals from the linear model
plot(lm_model$residuals ~ Chla_data_QC$DOY)
abline(h = 0, col = "red")

# 5. Statistical tests
Chla_data_QC$ID <- relevel(Chla_data_QC$ID, ref = "Control")

model_gamm_simple <- gamm(Chla ~ ID + s(DOY), 
                          random = list(Mesocosm = ~1), 
                          data = Chla_data_QC)

# Change the reference level of ID
model_gamm_interaction <- gamm(Chla ~ s(DOY, by = ID), 
                               random = list(Mesocosm = ~1), 
                               data = Chla_data_QC)
summary(model_gamm_simple$gam)
summary(model_gamm_interaction$gam)

# Refit the simple model
model_gamm_simple_releveled <- gamm(Chla ~ s(DOY) + ID, 
                                    random = list(Mesocosm = ~1), 
                                    data = Chla_data_QC)
# Check the summary
summary(model_gamm_simple_releveled$gam)

#Calculate descriptives
Chla_data_QC %>% group_by(ID) %>% reframe(mean(Chla),
                                          sd(Chla))

#______________________________________________________________________________
#                  P-I curve (K) - NLS and GAM Mixed effect global model   ====
#______________________________________________________________________________
#Fit Maelchin-Menkins curve using nls
# 1. Calculate the darkrate and corrections
Final_data <- left_join(DF_data_QC, Chla_data_QC, by = c("Mesocosm", "Date", "ID"))
Final_data <- Final_data %>% filter(Chla > 0) #Remove when chla is zero

Final_data$Photons_corrected <- Final_data$Photons - Final_data$Darkrate #Correct for darkrate
Final_data$Photons_normalized <- Final_data$Photons_corrected/Final_data$Chla #normalize with chl-a
Final_data <- Final_data %>% 
  mutate(Photons_normalized = case_when(Photons_normalized < 0 ~ 0, T~Photons_normalized)) %>%
  drop_na(Photons_normalized) %>%
  filter(ID != "Lake") %>% #remove lake samples
  select(ID, Mesocosm, Date, Light, Photons, Darkrate, Photons_corrected, Photons_normalized)

# 2. Caluclate the mean photons per light intensity
Final_data_mean <- Final_data %>% group_by(ID, Mesocosm, Date, Light) %>%
  reframe(Mean_photons_normalzied = mean(Photons_normalized)) %>%
  filter(Light < 60) #remove points where photounhibition occurs 

fit_michaelis_menten <- function(data) {
  tryCatch(
    nls(Mean_photons_normalzied ~ P_max * Light / (K + Light),
        data = data,
        start = list(P_max = quantile(data$Mean_photons_normalzied, 0.95, na.rm = TRUE),
                     K = quantile(data$Light, 0.25, na.rm = TRUE))),
    error = function(e) return(NA)  # Return NA if the model fails to fit
  )
}
# Fit the model for each Mesocosm and date 
models <- Final_data_mean %>%
  group_by(Mesocosm, Date, ID) %>%
  nest() %>%
  mutate(model = map(data, fit_michaelis_menten))

failed_models <- models %>%
  filter(is.na(model)) %>%
  select(Mesocosm, Date, ID)  # Keep relevant identifiers

print(failed_models)

# Extract K and P_max values from the model fitting
DF_K_data <- models %>%
  mutate(P_max = map_dbl(model, ~ coef(.x)["P_max"]),
         K = map_dbl(model, ~ coef(.x)["K"])) %>%
  select(ID, Mesocosm, Date, K, P_max)

fitted_curves <- models %>%
  mutate(predicted = map(model, ~ {
    # Ensure model is valid before proceeding
    if (inherits(.x, "nls")) {
      # Create a sequence of Light values for prediction
      new_light <- seq(min(Final_data_mean$Light, na.rm = TRUE),
                       max(Final_data_mean$Light, na.rm = TRUE), length.out = 100)
      
      # Generate predicted values
      tibble(Light = new_light, 
             fitted = predict(.x, newdata = data.frame(Light = new_light)))
    } else {
      tibble(Light = numeric(0), fitted = numeric(0))  # Return empty tibble for failed models
    }
  })) %>%
  unnest(predicted)  # Convert list column to dataframe

ggplot(Final_data_mean, aes(x = Light, y = Mean_photons_normalzied)) +
  geom_point(aes(color = ID), alpha = 0.7) +  # Raw data points
  geom_line(data = fitted_curves, aes(x = Light, y = fitted, group = interaction(Mesocosm, Date, ID)), 
            color = "blue", size = 1) +  # Fitted curves
  facet_grid(Date ~ Mesocosm) + 
  lims(x = c(NA, 60)) +
  labs(title = "Michaelis-Menten curve",
       x = "Light Intensity umol photons m-2 s-1",
       y = "Photosynthetic Rate photons /chla") +
  theme_minimal()

library(mgcv)
library(lmtest)

DF_K_data$DOY <- yday(DF_K_data$Date)

# 1. Plot the data to visually inspect non-linear trends
ggplot(DF_K_data, aes(x = DOY, y = K, color = ID)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  theme_minimal()

# 2. Fit linear, polynomial, and non-linear models
# Linear model
lm_model <- lm(K ~ DOY , data = DF_K_data)
summary(lm_model)

# Polynomial model
poly_model <- lm(K ~ poly(DOY, 2), data = DF_K_data)
summary(poly_model)

# Compare models using AIC
AIC(lm_model, poly_model)
Big_data_filter <- DF_K_data %>% drop_na(K)
# 3. Fit a Generalized Additive Model (GAM)
gam_model <- gam(K ~ s(DOY, bs = "cs", k = 4), data = DF_K_data)
summary(gam_model)

# Plot the fitted GAM
plot(gam_model, residuals = TRUE)

# 4. Examine residuals from the linear model
plot(lm_model$residuals ~ Big_data_filter$DOY)
abline(h = 0, col = "red")

# 5. Statistical tests
# Ramsey's RESET test for non-linearity
resettest(lm_model)
Big_data_filter$ID <- as.factor(Big_data_filter$ID)
Big_data_filter$ID <- relevel(Big_data_filter$ID, ref = "Early")

model_gamm_simple <- gamm(K ~ ID + s(DOY, k = 4), 
                          random = list(Mesocosm = ~1), 
                          data = Big_data_filter)

# model_gamm_interaction <- gamm(K ~ s(DOY, by = ID, k = 4), 
#                                random = list(Mesocosm = ~1), 
#                                data = Big_data_filter)

summary(model_gamm_simple$gam)
#summary(model_gamm_interaction$gam)

#______________________________________________________________________________
#                  Trophic modes - NLS and GAM Mixed effect global model   ====
#______________________________________________________________________________
Trophic_mode_data <- Final_biovolume_data
Trophic_mode_data$ID <- factor(Trophic_mode_data$ID)

Trophic_mode_data1 <- Trophic_mode_data %>% filter(Trophic_mode == "Mixotroph")
Trophic_mode_data2 <- Trophic_mode_data %>% filter(Trophic_mode == "Autotroph")

Trophic_mode_data1$ID <- relevel(Trophic_mode_data1$ID, ref = "Control")

l_gamm_simple <- gamm(Sum_biovolume ~ ID + s(DOY, k = 4),
                      random = list(Mesocosm = ~1),
                      data = Trophic_mode_data1)

model_gamm_interaction <- gamm(Sum_biovolume ~ s(DOY, by = ID, k = 4), 
                               random = list(Mesocosm = ~1), 
                               data = Trophic_mode_data1)
summary(l_gamm_simple$gam)
summary(model_gamm_interaction$gam)

# Change the reference level of ID
Trophic_mode_data2$ID <- relevel(Trophic_mode_data2$ID, ref = "Late")
#Trophic_mode_data2 <- Trophic_mode_data2 %>% filter(Date > as.Date("2023-02-06")) 
l_gamm_simple <- gamm(Sum_biovolume ~ ID + s(DOY, k = 4),
                      random = list(Mesocosm = ~1),
                      data = Trophic_mode_data2)

model_gamm_interaction <- gamm(Sum_biovolume ~ s(DOY, by = ID, k = 4), 
                               random = list(Mesocosm = ~1), 
                               data = Trophic_mode_data2)
summary(l_gamm_simple$gam)
summary(model_gamm_interaction$gam)

model_gamm_interaction$gam


#______________________________________________________________________________
#                 Rotifer abundance - NLS and GAM Mixed effect global model ====
#______________________________________________________________________________

Zooplankton_data$DOY <- lubridate::yday(Zooplankton_data$Date)
Zooplankton_data$Treatment <- as.factor(Zooplankton_data$Treatment)

Zooplankton_data$Treatment <- relevel(Zooplankton_data$Treatment, ref = "Late")

model_gamm_simple <- gamm(ind_L ~ s(DOY, k = 4) + Treatment, 
                          random = list(Mesocosm = ~1), 
                          data = filter(Zooplankton_data, Order == "Rotifer"))


model_gamm_interaction <- gamm(ind_L ~ s(DOY, by = Treatment, k = 4), 
                               random = list(Mesocosm = ~1), 
                               data = filter(Zooplankton_data, Order == "Rotifer"))

summary(model_gamm_simple$gam)
summary(model_gamm_interaction$gam)

# PLOTTING ####
#______________________________________________________________________________
#                  FIGURE 1; SENSOR DATA                                   ====
#______________________________________________________________________________
Sensor_data <- Sensor_data %>% separate(Date, into = c("Date", "Time"), sep = " ")
Sensor_data
#ALL SENSOR DATA
Sensor_data_mean <- Sensor_data %>% filter(Date <= as.Date("2023-03-13") & Date >= as.Date("2023-02-06")) #Remove dates when teh sensors were out of the water

#Calculate diurnal fluctuations
Sensor_data_diurnal <- Sensor_data %>% filter(Date <= as.Date("2023-03-13") & Date >= as.Date("2023-02-14")) %>% filter(Sensor == "Temp") %>%
  group_by(Date, Mesocosm) %>%
  reframe(Mean_temp = mean(Measurement),
          SD_temp = sd(Measurement)) %>%
  reframe(Mean_SD = sd(SD_temp),
          Mean_mean = mean(Mean_temp))

#Calculate mean 
Sensor_data_mean <- Sensor_data_mean %>% 
  group_by(Date, Mesocosm, Sensor, ID) %>% 
  drop_na(Measurement) %>%
  reframe(Mean_measurement = mean(Measurement),
          Obs = n(),
          Sum_measurement = sum(Measurement)) %>%
  ungroup() %>%
  group_by(Date, Sensor, ID) %>%
  reframe(Mean_measurement_new = mean(Mean_measurement),
          sd_measurement = sd(Mean_measurement),
          Sum_measurement_new = mean(Sum_measurement),
          sd_max = sd(Sum_measurement),
          Obs = n()) %>%
  mutate(se_measurement = sd_measurement/(sqrt(Obs)),
         se_max = sd_max/(sqrt(Obs)))

Sensor_data_mean$Date <- as.Date(Sensor_data_mean$Date)
Sensor_data_mean <- left_join(Sensor_data_mean, Snow_ice_data_QC)
Sensor_data_mean <- Sensor_data_mean %>% filter(Sensor != "DOconc")

PAR_plot <- ggplot(filter(Sensor_data_mean, Sensor == "PAR"), aes(x = Date, y = Sum_measurement_new, group = ID)) + 
  geom_vline(aes(xintercept = as.POSIXct("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-06")), linetype = "dashed") +
  geom_point(mapping = aes(x = Date, y = Snow_depth_mesocosms*10), fill = "black", alpha = 0.2, size = 3) + 
  geom_line(mapping = aes(col = ID)) + 
  geom_point(mapping = aes(col = ID)) + 
  geom_errorbar(mapping = aes(ymin = Sum_measurement_new - sd_max, ymax = Sum_measurement_new + sd_max, col = ID)) + 
  scale_fill_viridis_d(option = "D", end = 0.6) + 
  theme_bw() + 
  theme(legend.title = element_blank()) +
  scale_color_viridis_d(option = "D", end = 0.6) + 
  labs(y = expression(paste("("*Sigma*mu*"mol "*"m"^-2*"s"^-1*")"))) +
  theme(axis.title.x = element_blank(),
        legend.position = "top") +
  scale_y_continuous(name = expression(paste(Sigma*"PAR ("*mu*"mol "*"m"^-2*"s"^-1*")")),
                     sec.axis = sec_axis(trans=~./10, name="Snow depth (cm)"))

DO_plot <- ggplot(filter(Sensor_data_mean, Sensor == "DOsat"), aes(x = Date, y = Mean_measurement_new, group = ID, col = ID)) + 
  geom_vline(aes(xintercept = as.POSIXct("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-06")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-09")), linetype = "dashed", colour = "red") +
  geom_line() +
  geom_point(mapping = aes(col = ID)) + 
  geom_errorbar(mapping = aes(ymin = Mean_measurement_new - sd_measurement, ymax = Mean_measurement_new + sd_measurement)) + 
  scale_fill_viridis_d(option = "D", end = 0.6) + 
  theme_bw() + 
  theme(legend.title = element_blank()) +
  scale_color_viridis_d(option = "D", end = 0.6) + 
  labs(y = expression(paste("DO (%)")))+ 
  theme(axis.title.x = element_blank())

Temp_plot <- ggplot(filter(Sensor_data_mean, Sensor == "Temp"), aes(x = Date, y = Mean_measurement_new, group = ID, col = ID)) +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-06")), linetype = "dashed") +
  geom_line() + 
  geom_point(mapping = aes(col = ID)) + 
  geom_errorbar(mapping = aes(ymin = Mean_measurement_new - sd_measurement, ymax = Mean_measurement_new + sd_measurement)) + 
  scale_fill_viridis_d(option = "D", end = 0.6) + 
  theme_bw() + 
  theme(legend.title = element_blank()) +
  scale_color_viridis_d(option = "D", end = 0.6) + 
  labs(y = expression(paste("Temperature ("*degree*"C)")))+ 
  theme(axis.title.x = element_blank())

library(ggpubr)
Sensor_plot <- ggarrange(PAR_plot, DO_plot, Temp_plot, ncol = 1, align = "hv", common.legend = T)

# write the figure into dir
setwd(fig_dir)
tiff("Figure_S2.tiff", width = 12, height = 13, unit = "cm", res = 1600) #Set new window size and re-plot previous plot
Sensor_plot
dev.off()
setwd(dir)
#______________________________________________________________________________
#                  FIGURE 2; BACTERIAL ABUNDANCE                           ====       
#______________________________________________________________________________

library(scales)
# plotting

# calculate n
FCM_data_QC %>% group_by(ID) %>% drop_na(Abundance_Bacteria_events_mL_V) %>% count()
# plot bacterial abundance per sampling during the experiment
Figure_2 <- ggplot(FCM_data_QC, aes(x = Date, y = Abundance_Bacteria_events_mL_V, fill = ID, group = interaction(Date, ID))) +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-06")), linetype = "dashed") +
  geom_smooth(se=F, aes(group = ID,  col = ID), alpha = 0.5, linetype = "dashed", linewidth = 0.6, span = 0.4) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x = "", y = expression(paste("Bacteria (cells mL"^"-1"*")"))) +
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_fill_viridis_d(option = "D", end = 0.6) + 
  scale_color_viridis_d(option = "D", end = 0.6) +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  scale_y_continuous(labels = label_scientific()) + 
  facet_wrap(~ID)


# write the figure into dir
setwd(fig_dir)
tiff("Figure_2.tiff", width = 20, height = 6.5, unit = "cm", res = 1200) #Set new window size and re-plot previous plot
Figure_2
dev.off()
setwd(dir)

#______________________________________________________________________________
#                  FIGURE 3; CHLA AND K P-I curve                          ====       
#______________________________________________________________________________
conflict_prefer(name = "yday", "lubridate")

Plotting_data <- Big_data %>% group_by(ID, Date) %>%
  drop_na(Chla) %>%
  reframe(Mean_chla = mean(Chla),
          sd_chla = sd(Chla))

# calculate n
Chla_data_QC %>% group_by(ID) %>% drop_na(Chla) %>% count()

DF_K_data %>% group_by(ID) %>% drop_na(K) %>% count()

# 1. Plot the data to visually inspect non-linear trends
Chla_boxplot <- ggplot(filter(Chla_data_QC, ID != "Lake"), aes(x = as.Date(Date), y = Chla, group = interaction(Date, ID), fill = ID)) + 
  geom_smooth(se=F, aes(group = ID,  col = ID), alpha = 0.5, linewidth = 0.6, span = 0.3, linetype = "dashed") +
  geom_vline(aes(xintercept = as.Date("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.Date("2023-02-06")), linetype = "dashed") +
  geom_boxplot(width = 2.5, position = position_dodge()) +
  theme_minimal() + 
  scale_fill_viridis_d(option = "D", end = 0.6) + theme_bw() +  
  scale_color_viridis_d(option = "D", end = 0.6) +
  labs(x = "", y = expression(paste("Chl-"*alpha*" ("*mu*"g L"^"-1"*")"))) +
  theme(legend.position = "none") +  theme(legend.title = element_blank())+ 
  #  lims(x = c(as.Date("2023-02-05"), as.Date("2023-04-18"))) + 
  facet_wrap(~ID)


PI_curve_K <- ggplot(DF_K_data, aes(x = Date, y = K, fill = ID, group = interaction(Date, ID))) + 
  geom_vline(aes(xintercept = as.POSIXct("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-06")), linetype = "dashed") +
  geom_smooth(se=F, aes(group = ID,  col = ID), alpha = 0.5, linetype = "dashed", linewidth = 0.6) + 
  geom_boxplot() + 
  labs(y = expression(paste(K["P-I curve"]~"(photons"%.%"L"%.%"m"^-1%.%mu*"g"^-1*" chl-"*alpha*")")), x = "") +
  theme_bw() + 
  scale_fill_viridis_d(option = "D", end = 0.6) + theme_bw() +  
  scale_color_viridis_d(option = "D", end = 0.6) +
  theme(legend.title = element_blank()) + theme(legend.title = element_blank()) + 
  lims(x = c(as.POSIXct("2023-02-05"), as.POSIXct("2023-04-18"))) +
  facet_wrap(~ID) + 
  theme(strip.background = element_blank(),
        strip.text = element_blank()) 

library(patchwork)
Figure_3 <- Chla_boxplot/PI_curve_K + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "none") &
  plot_annotation(tag_levels = 'a',
                  tag_sep = '',
                  tag_prefix = '',
                  tag_suffix = '') &
  theme(plot.tag.position = "topleft",
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0)) &
  theme(plot.tag = element_text(face = 'bold'))

# write the figure into dir
setwd(fig_dir)
tiff("Figure_3.tiff", width = 20, height = 12, unit = "cm", res = 1600) # Set new window size and re-plot previous plot
Figure_3
dev.off()
setwd(dir)

#______________________________________________________________________________
#                  FIGURE 4; PHYTOPLANKTON TROPHIC MODE                    ====
#______________________________________________________________________________
# calculate n
Final_biovolume_data %>% group_by(ID, Trophic_mode) %>% drop_na(Sum_biovolume) %>% count()

Trophic_mode_boxplot <- ggplot(filter(Final_biovolume_data), aes(x = Date, y = Sum_biovolume, fill = ID, group = interaction(Date, ID))) + 
  geom_vline(aes(xintercept = as.POSIXct("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-06")), linetype = "dashed") +
  geom_smooth(se=F, aes(group = ID,  col = ID), alpha = 0.5, linetype = "dashed", linewidth = 0.6, span = 0.4) + 
  geom_boxplot() + 
  labs(y = expression(paste("Biovolume ("*mu*"m"^3*"L"^-1*")")), x = "") +
  theme_bw() + 
  scale_fill_viridis_d(option = "D", end = 0.6) + theme_bw() +  
  scale_color_viridis_d(option = "D", end = 0.6) +
  lims(y = c(0, NA)) +
  theme(legend.title = element_blank()) + theme(legend.title = element_blank()) +
  facet_wrap(~Trophic_mode, scales = "free") + theme(legend.position = "top")
Trophic_mode_boxplot

# write the figure into dir
setwd(fig_dir)
tiff("Figure_4.tiff", width = 14, height = 8, unit = "cm", res = 1600) #Set new window size and re-plot previous plot
Trophic_mode_boxplot
dev.off()
setwd(dir)

#______________________________________________________________________________
#                  FIGURE 5; ZOOPLANKTON ABUNDANCE                         ====  
#______________________________________________________________________________
Zooplankton_data

Zooplankton_data$Order <- factor(Zooplankton_data$Order, levels = c("Rotifer", "Copepod"))
# calculate n
Zooplankton_data %>% filter(Order == "Rotifer") %>% group_by(Treatment) %>% count()

# plot zooplankton abundance per sampling during the experiment
Figure_5 <- ggplot(filter(Zooplankton_data, Order == "Rotifer"), aes(x = Date, y = ind_L, fill = Treatment, group = interaction(Date, Treatment))) +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-06")), linetype = "dashed") +
  geom_smooth(se = F, aes(group = Treatment, col = Treatment, fill = Treatment), alpha = 0.5, linetype = "dashed", linewidth = 0.6, span = 0.4) + 
  geom_boxplot() + 
  theme_bw() +
  labs(x = "", y = expression(paste("Rotifer abundance (ind. L"^"-1"*")"))) +
  theme(legend.position = "top",
        legend.title = element_blank()) +  # Added hjust and vjust
  scale_fill_viridis_d(option = "D", end = 0.6) +  
  scale_color_viridis_d(option = "D", end = 0.6) +
  facet_wrap(~Treatment)
#filter(Zooplankton_data_order, Order == "Rotifer")
# write the figure into dir
setwd(fig_dir)
tiff("Figure_5.tiff", width = 17, height = 7, unit = "cm", res = 1600) #Set new window size and re-plot previous plot
Figure_5
dev.off()
setwd(dir)


# plot zooplankton abundance per sampling during the experiment
Figure_SX <- ggplot(filter(Zooplankton_data, Order == "Copepod"), aes(x = Date, y = ind_10L, shape = Treatment, col = Treatment)) +
  scale_shape_manual(values = c(15,16,17)) +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-28")), linetype = "dashed") +
  geom_vline(aes(xintercept = as.POSIXct("2023-02-06")), linetype = "dashed") + 
  geom_point(size = 4, alpha = 0.4) + 
  theme_bw() +
  labs(x = "", y = expression(paste("Copepod abundance (ind. 10 L"^"-1"*")"))) +
  theme(legend.position = "top",
        legend.title = element_blank()) +  # Added hjust and vjust
  scale_fill_viridis_d(option = "D", end = 0.6) +  
  scale_color_viridis_d(option = "D", end = 0.6) 

setwd(fig_dir)
tiff("Figure_SX.tiff", width = 14, height = 8, unit = "cm", res = 600) #Set new window size and re-plot previous plot
Figure_SX
dev.off()
setwd(dir)


# SUPPLEMENTS ####
#______________________________________________________________________________
#                  SENSOR DATA - TABLE S2 + DESCRIPTIVES      ====
#______________________________________________________________________________
#CAlculate mean (2 times)

##Model on DO after addition##
Sensor_data <- Sensor_data %>% separate(Date, into = c("Date", "Time"), sep = " ")
Sensor_data_mean <- Sensor_data %>% filter(Date <= as.Date("2023-03-13") & Date >= as.Date("2023-02-06")) #Remove dates when teh sensors were out of the water

#Calculate mean 
Sensor_data_mean <- Sensor_data_mean %>% 
  group_by(Date, Mesocosm, Sensor, ID) %>% 
  drop_na(Measurement) %>%
  reframe(Mean_measurement = mean(Measurement),
          Obs = n(),
          Sum_measurement = sum(Measurement))

DO_data <- Sensor_data_mean %>%
  filter(Sensor == "DOsat") %>%
  filter(Date > as.Date("2023-02-08")) %>%
  mutate(Treatment = case_when(Date >= as.Date("2023-02-28")~"Post", T~"Pre"))

DO_data$Mesocosm <- as.factor(DO_data$Mesocosm)

#DO_data <- relevel(DO_data$ID)
library(lmerTest)
model <- lmerTest::lmer(Mean_measurement ~ ID * Treatment + (1 | Date) + (1|Mesocosm), data = DO_data)
summary(model)  #included p-values for fixed effects



#PAR DAILY MEAN
Sensor_data_mean <- Sensor_data %>% filter(Date <= as.Date("2023-03-13") & Date >= as.Date("2023-02-06")) %>%
  filter(Sensor == "PAR")

Sensor_data_mean %>%
  summarise(min_PAR = min(Measurement),
            max_PAR = max(Measurement),
            median_PAR = median(Measurement))
#Calculate mean 
PAR_data <- Sensor_data_mean %>% 
  group_by(Date, Mesocosm, Sensor) %>% 
  drop_na(Measurement) %>%
  reframe(Obs = n(),
          mean_measurement = mean(Measurement)) %>%
  ungroup() %>%
  group_by(Sensor, Date) %>%
  reframe(Mean_mean_measurement = mean(mean_measurement),
          sd = sd(mean_measurement),
          Obs = n()) %>%
  group_by(Sensor)

PAR_data$Date <- as.Date(PAR_data$Date)
PAR_data <- left_join(PAR_data, Snow_ice_data_QC)

#Snow depth - PAR
ggplot(PAR_data, aes(y = Mean_mean_measurement, x = Snow_depth_mesocosms)) +
  geom_point(size = 3) + theme_bw()

PAR_data %>% filter(Snow_depth_mesocosms > 0) %>%
  reframe(Mean = mean(Mean_mean_measurement),
          sd = sd(Mean_mean_measurement),
          Obs = n())

# # FIT exponential to PAR and snow depth data
# # Your dataset
# data <- tibble(
#   Date = as.Date(c("2023-02-06", "2023-02-10", "2023-02-14", 
#                    "2023-02-28", "2023-03-03", "2023-03-06", 
#                    "2023-03-10", "2023-03-13")),
#   Snow_depth = c(0, 0, 0, 8, 4, 14, 40, 40),  # Snow depth in mesocosms
#   Mean_PAR = c(1.70, 2.78, 2.00, 1.04, 2.36, 0.864, 0.533, 0.566)  # Mean daily PAR
# )
# 
# # Fit an exponential model
# model <- nls(Mean_PAR ~ I0 * exp(-k_s * Snow_depth),
#              data = data,
#              start = list(I0 = 3, k_s = 0.1))  # Initial guesses
# 
# summary(model)  # Get estimated k_s
# 
# # Snow cover reduces PAR 4.7% per cm (Ks = 0.04694)
# 
# Plot the relationship
ggplot(data, aes(x = Snow_depth, y = Mean_PAR)) +
  geom_point(size = 3) +
  geom_line(aes(y = predict(model)), color = "black", linewidth = 1, linetype = "dotted") +
  labs(x = "Snow Depth (cm)", y = "Mean Daily PAR (µmol m⁻² s⁻¹)",
       title = "Light Attenuation by Snow Depth") +
  theme_minimal()



#Calculate the PAR during the experiment
PAR_data %>% 
  group_by(Sensor) %>%
  reframe(Mean_mean_mean_measurement = mean(Mean_mean_measurement),
          sd_daly = sd(Mean_mean_measurement),
          Obs = n())



#DIURNAL TEMP FLUXES
Diurnal_temp_data <- Sensor_data %>% 
  filter(Sensor == "Temp") %>% 
  filter(Date <= as.Date("2023-03-13") & Date >= as.Date("2023-02-14")) %>% #After temperatures were stable from filling
  mutate(Date = as.Date(Date)) %>%
  group_by(Date, Mesocosm, ID) %>%
  reframe(Max_measurement = max(Measurement),
          Min_measurement = min(Measurement)) %>%
  mutate(Range = Max_measurement - Min_measurement) %>%
  reframe(Mean_range = mean(Range))



Diurnal_plot <- ggplot() + 
  geom_point(Diurnal_temp_data, mapping = aes(x = Date, y = Max_measurement)) + 
  geom_point(Diurnal_temp_data, mapping = aes(x = Date, y = Min_measurement)) + 
  theme_bw() + 
  theme(legend.title = element_blank()) +
  scale_color_viridis_d(option = "D", end = 0.6) + 
  labs(y = expression(paste("Temperature ("*degree*"C)")))+ 
  theme(axis.title.x = element_blank())
Diurnal_plot

#CALCULATE MEAN SNOW AND ICE THICKNESS
Snow_ice_data_QC %>% drop_na(Snow_depth_mesocosms) %>% count()
mean(Snow_ice_data_QC$Snow_depth_mesocosms, na.rm = T)
sd(Snow_ice_data_QC$Snow_depth_mesocosms, na.rm = T)


Snow_ice_data_QC %>% drop_na(Ice_thickness_mesocosms) %>% count()
mean(Snow_ice_data_QC$Ice_thickness_mesocosms, na.rm = T)
sd(Snow_ice_data_QC$Ice_thickness_mesocosms, na.rm = T)


#______________________________________________________________________________
#                  NP uptake rates - FIGURE S4 AND FIGURE S1               ====
#______________________________________________________________________________
CNP_data <- CNP_data_QC %>% select(Mesocosm, Date, ID, DOC, TOC, DN, TN, DP, TP) %>% gather(Nutrient, Conc, -Mesocosm, -Date, -ID) %>% filter(ID != "Lake")

#AVERAGE PER TREATMENT
Split_data <- CNP_data  %>% group_by(ID, Date, Nutrient, Mesocosm) %>% 
  mutate(Paste_ID = paste0(ID, Nutrient, Mesocosm)) %>%
  mutate(DOY = yday(Date)) %>%
  mutate(Conc_ug_L = case_when(Nutrient == "TN" | Nutrient == "DN" ~ Conc*1000, T~Conc))

# Loop to create breakpoint and slope for each Mesocosm and each nutrient
Split_data_unique <- as.data.frame(unique(Split_data$Paste_ID))
names(Split_data_unique) <- "Paste_ID"
Breakpoint <- rep(NA, nrow(Split_data_unique))
Breakpoint <- as.data.frame(Breakpoint) %>% mutate(Paste_ID = "NA")
Slope <- rep(NA, nrow(Split_data_unique))
Slope <- as.data.frame(Slope) %>% mutate(Paste_ID = "NA")

#Calculate breakpoints per mesocosm
for (i in 1:nrow(Split_data_unique)) {
  Breakpoint[i,1] <- davies.test(lm(Conc_ug_L ~ DOY, data = filter(Split_data, Paste_ID == Split_data_unique[i,1])), seg.Z = ~Date)[3]$statistic
  Breakpoint[i,2] <- Split_data_unique[i,1]
}

Split_data2 <- left_join(Split_data, Breakpoint, by = "Paste_ID")
Split_data2 <- Split_data2 %>% 
  ungroup() %>%
  group_by(Nutrient) %>%
  mutate(Mean_breakpoint = mean(Breakpoint)) %>%
  mutate(Mean_breakpoint_date = as.Date(as.POSIXct(Mean_breakpoint))) %>%
  filter(Date < as.Date("2023-03-31")) %>%
  filter(Nutrient == "DP" | Nutrient =="DN") %>%
  ungroup() %>%
  mutate(Mean_breakpoint_date = mean(Mean_breakpoint_date))


Split_data2 <- Split_data2 %>% mutate(Slope_range_min = case_when(Nutrient == "DN"~as.Date("2023-03-13"),
                                                                  Nutrient == "DP"~as.Date("2023-03-13"), 
                                                                  T~NA)) %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(Slope_range_max = as.Date("2023-03-31"))

# Filtering the subset correctly
filtered_data <- Split_data2 %>%
  filter(Nutrient %in% c("DN", "DP") & Date >= Slope_range_min & Date <= Slope_range_max)

# Plot on the range which the slopes are calculated from
p_CNP_slopes <- ggplot() +
  theme_bw() +
  labs(y = expression(paste("Nutrient conc. ("*mu*"g L"^-1*")"))) +
  geom_point(filter(Split_data2, Nutrient %in% c("DN", "DP")), mapping = aes(y = Conc_ug_L, x = Date, group = Date, fill = ID, col = ID)) +
  scale_fill_viridis_d(option = "D") +
  geom_vline(data = data.frame(xint = as.Date("2023-02-28"), ID = "Late"), aes(xintercept = xint), linetype = "dotted") +
  geom_vline(data = data.frame(xint = as.Date("2023-02-06"), ID = "Early"), aes(xintercept = xint), linetype = "dotted") +
  geom_vline(filter(Split_data2, Nutrient %in% c("DN", "DP")), mapping = aes(xintercept = Mean_breakpoint_date), linetype = "dotted", col = "red") +
  facet_grid(Nutrient ~ ID, scales = "free_y") +
  lims(y = c(0, NA)) +
  geom_smooth(data = filtered_data, aes(y = Conc_ug_L, x = Date, group = Mesocosm), method = "lm", se = F, col = "darkred", linetype = "dashed", size = 0.5) +
  scale_color_viridis_d(option = "D", end = 0.6) + 
  theme(legend.position = "none",
        axis.title.x = element_blank())

p_CNP
#write figure for supplements
setwd(fig_dir)
tiff("Figure_S1.tiff", width = 15, height = 8, units = "cm", res = 600)
p_CNP_slopes
dev.off()
setwd(dir)

# Rsqr values# Rsqr velement_blank()alues
latdry_table_r <- filtered_data %>% group_by(Nutrient, ID, Mesocosm) %>% lm_table(Conc ~ Date) %>%
  select(Nutrient, ID, Mesocosm, Rsqr_adj, b1) %>% rename("R2" = "Rsqr_adj", "Slope" = "b1") %>%
  mutate(R2 = round(R2, 3))

# p values
latdry_table_p <- filtered_data %>% nest(data = -c(Nutrient, ID, Mesocosm)) %>%
  mutate(model = map(data, ~lm(Conc ~ Date, data = .)), tidied = map(model, tidy)) %>%
  unnest(tidied) %>% filter(term == "Date") %>% select(Nutrient, ID, Mesocosm, p.value) %>%
  rename("p" = "p.value") %>% mutate(p = round(p, 3))

# combine
latdry_table <- left_join(latdry_table_r, latdry_table_p, by = c("Nutrient", "ID", "Mesocosm"))

Slope_table <- latdry_table %>% ungroup() %>%  select(-R2, -p) %>% spread(Nutrient, value = Slope)

# Plot the results
p_uptake <- ggplot() + geom_boxplot(filter(latdry_table, Nutrient != "DOC"), mapping = aes(x = ID, y = abs(Slope), group = ID, fill = ID)) +
  geom_point(filter(latdry_table, Nutrient != "DOC"), mapping = aes(x = ID, y = abs(Slope), group = ID, fill = ID), alpha = 0.5, size = 1) +
  facet_grid2(~Nutrient, scales = "free", independent = T) + 
  theme_bw() +
  scale_fill_viridis_d(option = "D", end = 0.6) + 
  labs(y = expression(paste("Uptake ("*mu*"g (L"%*%"d)"^-1*")"))) +
  theme(legend.title = element_blank(), axis.title.x = element_blank()) +
  lims(y = c(0, NA))

#write figure for supplements
setwd(fig_dir)
tiff("Figure_S3.tiff", width = 15, height = 6, units = "cm", res = 600)
p_uptake
dev.off()
setwd(dir)

# Plot the results
ggplot() + geom_boxplot(filter(Slope_anova_data, Nutrient != "DOC"), mapping = aes(x = ID, y = Slope, group = ID, fill = ID)) +
  facet_wrap(~Nutrient, scales = "free") + theme_bw() + scale_fill_viridis_d(option = "D", end = 0.6)

p_CNP <- ggplot(filter(Split_data2, Nutrient == "DP" | Nutrient == "DN"), aes(y = Conc, x = Date, group = Date, fill = ID)) + theme_bw() +
  labs(y = "Nutrient conc (mg/L)") + geom_point() + scale_fill_viridis_d(option = "D") +
  geom_vline(data = data.frame(xint=as.POSIXct("2023-02-28"), ID = "Late"), aes(xintercept = xint), linetype = "dotted") +
  geom_vline(data = data.frame(xint=as.POSIXct("2023-02-06"), ID = "Early"), aes(xintercept = xint), linetype = "dotted") +
  facet_grid2(ID~Nutrient, scales = "free", independent = T) + 
  lims(y = c(0,NA)) +
  geom_vline(filter(Split_data2, Nutrient == "DP" | Nutrient == "DN"), mapping = aes(xintercept = Mean_breakpoint_date), linetype = "dotted", col = "red")
p_CNP

#Plotting and analyses
p_CNP <- ggplot(filter(CNP_data, ID != "Lake"), aes(y = Conc, x = Date, group = Date, shape = Form, fill = ID)) + theme_bw() +
  labs(y = "Nutrient conc (mg/L)") + geom_point() + scale_fill_viridis_d(option = "D") + 
  geom_vline(data = data.frame(xint=as.Date("2023-02-28"), ID = "Late"), aes(xintercept = xint), linetype = "dotted") + 
  geom_vline(data = data.frame(xint=as.Date("2023-02-06"), ID = "Early"), aes(xintercept = xint), linetype = "dotted") +
  facet_grid2(Nutrient+Form~ID, scales = "free_y") + lims(y = c(0,NA))
p_CNP

setwd(fig_dir)
tiff("CNP_sd.tiff", width = 20, height = 15, units = "cm", res = 600)
p_CNP
dev.off()
setwd(dir)

# Stats
# ANOVA on slope of DP, and DN
library(multcomp)
library(rstatix)
library(ggpubr)

Slope_anova_data <- latdry_table %>% select(Slope, Nutrient, ID)

Slope_anova_data_DP <- Slope_anova_data %>% filter(Nutrient == "DP") #Select DP data
Anova3 <- aov(Slope~ID, data = Slope_anova_data_DP)
shapiro.test(residuals(Anova3)) #Test for normality
ggqqplot(residuals(Anova3)) #Check QQ plot
summary(Anova3) #Check results
TukeyHSD(Anova3) #Run pairwise comparison

Slope_anova_data_DN <- Slope_anova_data %>% filter(Nutrient == "DN") #Select DN data
Anova3 <- aov(Slope~ID, data = Slope_anova_data_DN)
shapiro.test(residuals(Anova3)) #Test for normality
ggqqplot(residuals(Anova3)) #Check QQ plot
summary(Anova3) #Check results
TukeyHSD(Anova3) #Run pairwise comparison

#______________________________________________________________________________
#                  CNP OVER EXPERIMENT - FIGURE S3                        ====
#______________________________________________________________________________
CNP_data <- CNP_data_QC %>%
  gather(Nutrient, Conc, -Date, -ID, -Mesocosm) %>%
  mutate(Form = case_when(Nutrient %in% c("DOC", "DP", "DN") ~ "Dissolved", 
                          TRUE ~ "Total")) %>%
  mutate(Nutrient = case_when(Nutrient %in% c("DOC", "TOC") ~ "C",
                              Nutrient %in% c("DP", "TP") ~ "P",
                              Nutrient %in% c("DN", "TN") ~ "N", 
                              TRUE ~ NA_character_)) %>%
  mutate(Date = as.Date(Date)) %>%
  mutate(Conc = case_when(Nutrient == "N" ~ Conc*1000, T~Conc)) #Make N into ug/L instead of mg /L

library(ggpubr)
# Create the plot
p_CNP <- ggplot(filter(CNP_data, ID != "Lake"), 
                aes(y = as.numeric(Conc), x = Date, group = Date, shape = Form, col = Form)) + 
  theme_bw() +
  labs(y = NULL, x = "Date") +  # Remove general y-axis label
  geom_point() + 
  scale_fill_viridis_d(option = "D") + 
  geom_vline(data = data.frame(xint = as.Date("2023-02-28"), ID = "Late"), 
             aes(xintercept = xint), linetype = "dotted") + 
  geom_vline(data = data.frame(xint = as.Date("2023-02-06"), ID = "Early"), 
             aes(xintercept = xint), linetype = "dotted") +
  facet_grid2(Nutrient ~ ID, scales = "free_y") + 
  #  scale_color_viridis_d(option = "D", end = 0.8) + 
  theme(legend.title = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12)) + 
  scale_color_manual(values=colorRampPalette(RColorBrewer::brewer.pal(n = 3, name = "BrBG")[1:3])(length(unique(CNP_data$Form)))) 
#scale_color_manual(values = c("#180F3EFF", "#F1605DFF"))

final_plot <- annotate_figure(p_CNP,  
                              left = text_grob(expression(paste("      Conc. ("*mu*"g L"^-1*")"*"          Conc. ("*mu*"g L"^-1*")"*"        Conc. (mg L"^-1*")")), 
                                               rot = 90, vjust = 0.5, hjust = 0.5, size = 11))

final_plot

setwd(fig_dir)
tiff("Figure_S3.tiff", width = 22, height = 12, units = "cm", res = 1200)
final_plot
dev.off()
setwd(dir)


#______________________________________________________________________________
#                  Calculate P assimilation?                        ====
#______________________________________________________________________________
# Load necessary libraries
library(dplyr)
library(tidyr)
conflicts_prefer(dplyr::lag)
# Set phosphorus quotas (min and max values)
Qp_phyto_min <- 0.1   # Min P quota for phytoplankton (pg P / ?g chl-a)
Qp_phyto_max <- 2.0   # Max P quota for phytoplankton (pg P / ?g chl-a)
Qp_bact_min  <- 10    # Min P quota for bacteria (fg P / cell)
Qp_bact_max  <- 50    # Max P quota for bacteria (fg P / cell)

# Convert bacterial quotas to pg P
Qp_bact_min_pg <- Qp_bact_min * 1e-3
Qp_bact_max_pg <- Qp_bact_max * 1e-3

# Calculate changes and P assimilation between timepoints
P_assimilation_results <- Big_data_wide %>%
  group_by(Mesocosm, ID) %>%
  arrange(Date) %>%  # Ensure data is ordered by Date
  mutate(
    delta_chla = Chla - lag(Chla),  # Change in chl-a
    delta_bact = Abundance_Bacteria_events_uL_V - lag(Abundance_Bacteria_events_uL_V),  # Change in bacteria
    delta_DP   = DP - lag(DP),  # Change in dissolved phosphorus (DP)
    P_bact_min  = ifelse(!is.na(delta_bact), delta_bact * Qp_bact_min_pg, 0),  # Bacteria P assimilation (min)
    P_bact_max  = ifelse(!is.na(delta_bact), delta_bact * Qp_bact_max_pg, 0)   # Bacteria P assimilation (max)
  ) %>%
  ungroup() %>%
  select(ID, Mesocosm, Date, delta_chla, delta_bact, delta_DP, 
         P_bact_min, P_bact_max)

# Set the cutoff date
cutoff_date <- as.Date("2023-03-27")

# Calculate cumulative phosphorus assimilation and mean uptake
P_cumulative_summary <- P_assimilation_results %>%
  filter(Date <= cutoff_date) %>%  # Filter rows up to the cutoff date
  group_by(ID, Mesocosm) %>%
  summarise(
    Total_P_bact_min = sum(P_bact_min, na.rm = TRUE),
    Total_P_bact_max = sum(P_bact_max, na.rm = TRUE),
    Total_DP         = sum(delta_DP, na.rm = TRUE) # Sum change in DP over time
  ) %>%
  group_by(ID) %>%
  summarise(
    # Mean and SD of bacterial uptake
    Mean_uptakemax = mean(Total_P_bact_max, na.rm = TRUE),
    sd_uptakemax   = sd(Total_P_bact_max, na.rm = TRUE),
    Mean_uptakemin = mean(Total_P_bact_min, na.rm = TRUE),
    sd_uptakemin   = sd(Total_P_bact_min, na.rm = TRUE),
    
    # Total dissolved P for the entire group (sum DP across Mesocosms)
    Total_DP = sum(Total_DP, na.rm = TRUE),
    
    # Calculate the % of DP taken up by bacteria
    Percent_uptake_min = (sum(Total_P_bact_min, na.rm = TRUE) / Total_DP) * 100,
    Percent_uptake_max = (sum(Total_P_bact_max, na.rm = TRUE) / Total_DP) * 100
  ) %>%
  ungroup()

# View the results
print(P_cumulative_summary)
