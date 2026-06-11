#'
#' # Dry Forest Growth Project
#' 
#' ## Script 2: Rainfall metrics calculation
#'
#' ### Author: German Vargas G.
#' 
#' ### Date: 2026-05-18
#' 
# Load packages ------
library(tidyverse)
library(lubridate)
library(ggpubr)
library(viridis)
#install.packages("tseries")
#install.packages("forecast")
library(forecast)
library(zoo)
#install.packages("imputeTS")
library(imputeTS)

# Load functions ----
source("code/functions/rainfall_seasonality_metrics.R")# Function to calculate rainfall season

# Load data -----
sr_rain <- read_csv(file = "data/SantaRosa_Rainfall_Monthly_1980_2025.csv")
pv_rain <- read_csv(file = "data/PaloVerde_Rainfall_Monthly_1999_2025.csv")

# Visual data summary matrices -----
png(filename = "output/sr_rain_monthly_heatmap.png", width = 5, height = 9, units = "in", res = 300)
sr_rain %>%
  # Convert from wide to long format for easier plotting
  pivot_longer(cols = jan:dec, # Assuming month columns are named jan, feb, ..., dec
               names_to = "month", # New column for month names
               values_to = "rain_mm") %>%
  # Ensure months are ordered correctly on the x-axis
  mutate(month = factor(month, levels = tolower(month.abb)),
         Year  = as.numeric(year))  %>%
  # Basic heatmap of rainfall across years and months
  ggplot(aes(x = month, y = factor(Year), fill = rain_mm)) +
  geom_tile(color = "white", linewidth = 0.1) +
  # viridis "mako" or "cividis" works beautifully for precipitation
  scale_fill_viridis_c(option = "mako", 
                       direction = -1, 
                       na.value = "gray90", # Highlights missing data/sensor dropouts clearly
                       name = "Rainfall\n(mm)") +
  theme_minimal(base_size = 12) +
  labs(title = "Santa Rosa Historical Precipitation Matrix",
       subtitle = "Monthly rainfall totals (1980 - 2025)",
       x = "Month",
       y = "Year") +
  theme_bw() + # Start with a black and white theme for a cleaner look
  theme(plot.title = element_text(face = "bold",size=14), # bold title
        plot.subtitle = element_text(size = 12),# subtitle in smaller font
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
        axis.line = element_line(color = "black"), # Add axis lines
        panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
        axis.text = element_text(size = 12), # Adjust axis text size
        axis.title = element_blank())
dev.off()

png(filename = "output/pv_rain_monthly_heatmap.png", width = 5, height = 7, units = "in", res = 350)
pv_rain %>%
  # Convert from wide to long format for easier plotting
  pivot_longer(cols = jan:dec, # Assuming month columns are named jan, feb, ..., dec
               names_to = "month", # New column for month names
               values_to = "rain_mm") %>%
  # Ensure months are ordered correctly on the x-axis
  mutate(month = factor(month, levels = tolower(month.abb)),
         Year  = as.numeric(year))  %>%
  # Basic heatmap of rainfall across years and months
  ggplot(aes(x = month, y = factor(Year), fill = rain_mm)) +
  geom_tile(color = "white", linewidth = 0.1) +
  # viridis "mako" or "cividis" works beautifully for precipitation
  scale_fill_viridis_c(option = "mako", 
                       direction = -1, 
                       na.value = "gray90", # Highlights missing data/sensor dropouts clearly
                       name = "Rainfall\n(mm)") +
  theme_minimal(base_size = 12) +
  labs(title = "Palo Verde Historical Precipitation Matrix",
       subtitle = "Monthly rainfall totals (1999 - 2025)",
       x = "Month",
       y = "Year") +
  theme_bw() + # Start with a black and white theme for a cleaner look
  theme(plot.title = element_text(face = "bold",size=14), # bold title
        plot.subtitle = element_text(size = 12),# subtitle in smaller font
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
        axis.line = element_line(color = "black"), # Add axis lines
        panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
        axis.text = element_text(size = 12), # Adjust axis text size
        axis.title = element_blank())
dev.off()

# Gap filling -----
# The rainfall climatology for Santa Rosa is complete, so we will only fill gaps for Palo Verde
pv <- ts(data = as.vector(t(as.matrix(pv_rain[2:27,3:14]))),start = c(2000,1),end = c(2025,12),frequency = 12)

# Check the missing data pattern and visualize it
statsNA(pv)

# This will automatically select the best ARIMA model for imputation and apply Kalman smoothing to fill in the missing values. 
# The 'smooth = TRUE' argument ensures that the imputed values are smoothed estimates rather than just point forecasts, 
# which can help maintain the overall structure of the time series.
pv_imp <- na_kalman(x = pv,model = "auto.arima",smooth = TRUE)

# Visualize the imputation results
png(file = "output/Fig_S_ImputedRainfallPV_kalmanfilter.png",width = 8, height = 6, units = "in", res = 300)
ggarrange(
  #check for empty data
  ggplot_na_distribution(x = pv)+
    ylab(bquote(Rainfall~"(mm)")),
  #check imputed data
  ggplot_na_imputations(pv, pv_imp)+
    ylab(bquote(Rainfall~"(mm)")),
  ncol = 1)
dev.off()

# Fill the missing values using the forecasted value from the ARIMA model
pv_rain[2:27,3:14] <- matrix(pv_imp, nrow = 26, ncol = 12, byrow = TRUE)

# Calculate rainfall seasonality metrics ----
pv_metrics <- seasonality(data = data.frame(pv_rain[-1,2:14]), Site = "Palo Verde")
pv_metrics
sr_metrics <- seasonality(data = data.frame(sr_rain[row_number(sr_rain)>15,2:14]), Site = "Santa Rosa")
sr_metrics

# Graphs -----
#### plot rainfall climatologies-----
pv_rain %>%
  pivot_longer(cols = jan:dec, names_to = "Month", values_to = "Rainfall") %>%
  group_by(Month) %>%
  summarise(mean_rainfall = mean(Rainfall, na.rm = TRUE)) %>%
  ggplot(aes(x = Month, y = mean_rainfall)) +
  geom_bar(stat = "identity", fill = "gray80",col="#000000") +
  labs(title = "Average Monthly Rainfall at Palo Verde (2000-2025)",
       x = "Month",
       y = "Average Rainfall (mm)") +
  theme_bw() + # Start with a black and white theme for a cleaner look
  theme(plot.title = element_text(face = "bold",size=14), # bold title
        plot.subtitle = element_text(size = 12),# subtitle in smaller font
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
        axis.line = element_line(color = "black"), # Add axis lines
        panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
        axis.text = element_text(size = 12), # Adjust axis text size
        axis.title = element_text(size = 14))+
  scale_x_discrete(limits = tolower(month.abb))
  
sr_rain %>% 
  filter(year >= 1995) %>% # we need to assess the climate in reference to the last 30 years!
  pivot_longer(cols = jan:dec, names_to = "Month", values_to = "Rainfall") %>%
  group_by(Month) %>%
  summarise(mean_rainfall = mean(Rainfall, na.rm = TRUE)) %>%
  ggplot(aes(x = Month, y = mean_rainfall)) +
  geom_bar(stat = "identity", fill = "gray80",col="#000000") +
  labs(title = "Average Monthly Rainfall at Santa Rosa (1995-2025)",
       x = "Month",
       y = "Average Rainfall (mm)") +
  theme_bw() + # Start with a black and white theme for a cleaner look
  theme(plot.title = element_text(face = "bold",size=14), # bold title
        plot.subtitle = element_text(size = 12),# subtitle in smaller font
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
        axis.line = element_line(color = "black"), # Add axis lines
        panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
        axis.text = element_text(size = 12), # Adjust axis text size
        axis.title = element_text(size = 12))+
  scale_x_discrete(limits = tolower(month.abb))


# Assuming 'complete_palo_verde_precipitation' is your master wide matrix
# Calculate historical monthly baselines (e.g., median, 25th, and 75th percentiles)
p1 <-   pv_rain %>%
  pivot_longer(cols = jan:dec, names_to = "month", values_to = "rain_mm") %>%
  mutate(month = factor(month, levels = tolower(month.abb)),
         Year  = as.numeric(year))  %>%
  group_by(month) %>%
  summarise(median_rain = mean(rain_mm, na.rm = TRUE),
            lower_q     = quantile(rain_mm, 0.10, na.rm = TRUE),
            upper_q     = quantile(rain_mm, 0.90, na.rm = TRUE),
            .groups = "drop") %>% 
  ggplot( aes(x = month, group = 1)) +
  ylim(0, 825) + # Set a consistent y-axis limit for better visual comparison between sites
  # 1. Shaded ribbon representing the interquartile range (historical variation)
  geom_ribbon(aes(ymin = lower_q, ymax = upper_q), fill = "steelblue", alpha = 0.2) +
  # 2. Strong central line tracking the historical median path
  geom_line(aes(y = median_rain), color = "steelblue4", linewidth = 1.2) +
  geom_point(aes(y = median_rain), color = "steelblue4", size = 2.5) +
  theme_bw(base_size = 12) +
  labs(title = "Palo Verde Climatology (2000-2025)",
       x = "Month",y = "Cumulative Rainfall (mm)") +
  theme_bw() + # Start with a black and white theme for a cleaner look
  theme(plot.title = element_text(hjust=0.5,face = "bold",size=14), # bold title
        plot.subtitle = element_text(size = 12),# subtitle in smaller font
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
        axis.line = element_line(color = "black"), # Add axis lines
        panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
        axis.text = element_text(size = 12), # Adjust axis text size
        axis.title = element_text(size = 14))+
  scale_x_discrete(limits = tolower(month.abb))

p2 <-     sr_rain %>% 
  filter(year >= 1995) %>% # we need to assess the climate in reference to the last 30 years!
  pivot_longer(cols = jan:dec, names_to = "month", values_to = "rain_mm") %>%
  mutate(month = factor(month, levels = tolower(month.abb)),
         Year  = as.numeric(year))  %>%
  group_by(month) %>%
  summarise(median_rain = mean(rain_mm, na.rm = TRUE),
            lower_q     = quantile(rain_mm, 0.10, na.rm = TRUE),
            upper_q     = quantile(rain_mm, 0.90, na.rm = TRUE),
            .groups = "drop") %>% 
  ggplot(aes(x = month, group = 1)) +
  ylim(0, 825) + # Set a consistent y-axis limit for better visual comparison between sites
  # 1. Shaded ribbon representing the interquartile range (historical variation)
  geom_ribbon(aes(ymin = lower_q, ymax = upper_q), fill = "steelblue", alpha = 0.2) +
  # 2. Strong central line tracking the historical median path
  geom_line(aes(y = median_rain), color = "steelblue4", linewidth = 1.2) +
  geom_point(aes(y = median_rain), color = "steelblue4", size = 2.5) +
  theme_bw(base_size = 12) +
  labs(title = "Santa Rosa Climatology (1995-2025)",
       x = "Month",y = "Cumulative Rainfall (mm)") +
  theme_bw() + # Start with a black and white theme for a cleaner look
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size=14), # bold title
        plot.subtitle = element_text(size = 12),# subtitle in smaller font
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines if not desired
        axis.line = element_line(color = "black"), # Add axis lines
        panel.background = element_rect(fill = "white", color = "black"), # White background, black border for plot area
        axis.text = element_text(size = 12), # Adjust axis text size
        axis.title = element_text(size = 14))+
  scale_x_discrete(limits = tolower(month.abb))

png(file = "output/Fig_S_Climatologies.png",width = 9, height = 4, units = "in", res = 300)
ggarrange(p1,p2,ncol = 2)
dev.off()

#### Cumulative rainfall profiles -----
library(scales) # For clean plot label formatting

# 1. Transform wide matrix to long and compute cumulative rain per year
pv_rain_cumulative <- pv_rain %>% filter(year > 1999) %>%
  pivot_longer(cols = jan:dec, names_to = "month", values_to = "rain_mm") %>%
  mutate(month = factor(month, levels = tolower(month.abb)),
         Year  = as.numeric(year)) %>%
  # Group by year to calculate a running sum as the months progress
  group_by(Year) %>%
  arrange(month) %>%
  mutate(cum_rain_mm = cumsum(coalesce(rain_mm, 0))) %>%
  ungroup()

# 2. Separate datasets for the historical background, the baseline, and the target anomaly
historical_background <- pv_rain_cumulative %>% filter(Year != 2015)
enso_2015             <- pv_rain_cumulative %>% filter(Year == 2015)

# Calculate historical median profile (excluding the anomaly year)
climatology_median <- historical_background %>%
  group_by(month) %>%
  summarise(mean_cum_rain = mean(cum_rain_mm, na.rm = TRUE), .groups = "drop")

# 1. Establish the maximum height dynamically to position the legend safely
# (Or manually hardcode y_pos if you prefer, e.g., y_pos = 2600)
y_pos <- max(historical_background$cum_rain_mm, na.rm = TRUE) * 0.95

p1 <- ggplot() +
  # 1. Background Spaghetti Lines (Historical years in light grey)
  geom_line(data = historical_background, aes(x = month, y = cum_rain_mm, group = Year), 
            color = "grey75", linewidth = 0.75, alpha = 0.45) +
  
  # 2. Historical Baseline (The Median path in solid dark grey or blue)
  geom_line(data = climatology_median, aes(x = month, y = mean_cum_rain, group = 1), 
            color = "#2b5c8f", linewidth = 1.1) +
  
  # 3. Highlighted Anomaly (ENSO 2015 in stark crimson)
  geom_line(data = enso_2015, aes(x = month, y = cum_rain_mm, group = 1), 
            color = "#b32424", linewidth = 1.1) +
  
  # 4. CUSTOM LEGEND: Historical Median (Line + Text)
  # X-axis factor levels are 1 (jan) to 12 (dec). We place these around jan (1) & feb (2).
  annotate("segment", x = 0.7, xend = 1.3, y = y_pos, yend = y_pos, 
           color = "#2b5c8f", linewidth = 1.2) +
  annotate("text", x = 1.4, y = y_pos, label = "Mean cumulative rainfall (2000-2025)", 
           color = "#2b5c8f", fontface = "bold", hjust = 0, size = 4) +
  
  # 5. CUSTOM LEGEND: 2015 ENSO (Line + Text)
  # Positioned slightly below the median label (shifted down by 7% of total Y scale)
  annotate("segment", x = 0.7, xend = 1.3, y = y_pos * 0.93, yend = y_pos * 0.93, 
           color = "#b32424", linewidth = 1.2) +
  annotate("text", x = 1.4, y = y_pos * 0.93, label = "2015 ENSO", 
           color = "#b32424", fontface = "bold", hjust = 0, size = 4) +
  
  # 6. Typography and Labels
  labs(title = "Palo Verde",
       x = "Month", 
       y = "Cumulative Rainfall (mm)") +
  
  # 7. Fine-tuning theme layouts
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  scale_x_discrete(limits = tolower(month.abb))

# 1. Transform wide matrix to long and compute cumulative rain per year
sr_rain_cumulative <- sr_rain %>% filter(year >= 1995) %>%
  pivot_longer(cols = jan:dec, names_to = "month", values_to = "rain_mm") %>%
  mutate(month = factor(month, levels = tolower(month.abb)),
         Year  = as.numeric(year)) %>%
  # Group by year to calculate a running sum as the months progress
  group_by(Year) %>%
  arrange(month) %>%
  mutate(cum_rain_mm = cumsum(coalesce(rain_mm, 0))) %>%
  ungroup()

# 2. Separate datasets for the historical background, the baseline, and the target anomaly
historical_background <- sr_rain_cumulative %>% filter(Year != 2015)
enso_2015             <- sr_rain_cumulative %>% filter(Year == 2015)

# Calculate historical median profile (excluding the anomaly year)
climatology_median <- historical_background %>%
  group_by(month) %>%
  summarise(mean_cum_rain = mean(cum_rain_mm, na.rm = TRUE), .groups = "drop")

# 1. Establish the maximum height dynamically to position the legend safely
# (Or manually hardcode y_pos if you prefer, e.g., y_pos = 2600)
y_pos <- max(historical_background$cum_rain_mm, na.rm = TRUE) * 0.95

p2 <- ggplot() +
  # 1. Background Spaghetti Lines (Historical years in light grey)
  geom_line(data = historical_background, aes(x = month, y = cum_rain_mm, group = Year), 
            color = "grey75", linewidth = 0.75, alpha = 0.45) +
  
  # 2. Historical Baseline (The Median path in solid dark grey or blue)
  geom_line(data = climatology_median, aes(x = month, y = mean_cum_rain, group = 1), 
            color = "#2b5c8f", linewidth = 1.1) +
  
  # 3. Highlighted Anomaly (ENSO 2015 in stark crimson)
  geom_line(data = enso_2015, aes(x = month, y = cum_rain_mm, group = 1), 
            color = "#b32424", linewidth = 1.1) +
  
  # 4. CUSTOM LEGEND: Historical Median (Line + Text)
  # X-axis factor levels are 1 (jan) to 12 (dec). We place these around jan (1) & feb (2).
  annotate("segment", x = 0.7, xend = 1.3, y = y_pos, yend = y_pos, 
           color = "#2b5c8f", linewidth = 1.2) +
  annotate("text", x = 1.4, y = y_pos, label = "Mean cumulative rainfall (1995-2025)", 
           color = "#2b5c8f", fontface = "bold", hjust = 0, size = 4) +
  
  # 5. CUSTOM LEGEND: 2015 ENSO (Line + Text)
  # Positioned slightly below the median label (shifted down by 7% of total Y scale)
  annotate("segment", x = 0.7, xend = 1.3, y = y_pos * 0.93, yend = y_pos * 0.93, 
           color = "#b32424", linewidth = 1.2) +
  annotate("text", x = 1.4, y = y_pos * 0.93, label = "2015 ENSO", 
           color = "#b32424", fontface = "bold", hjust = 0, size = 4) +
  
  # 6. Typography and Labels
  labs(title = "Santa Rosa",
       x = "Month", 
       y = "Cumulative Rainfall (mm)") +
  
  # 7. Fine-tuning theme layouts
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"), 
        panel.background = element_rect(fill = "white", color = "black"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  scale_x_discrete(limits = tolower(month.abb))
p2
png(file = "output/Fig_S_CumulativeRainfallProfiles.png",width = 9, height = 4, units = "in", res = 300)
ggarrange(p1,p2,ncol = 2,align = "hv")
dev.off()
