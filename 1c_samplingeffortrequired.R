# Determine how many times I have to visit a site in order to be 99% sure that if a BTNW 
# is present it will be detected. 
# Method of determining how many recordings I have to process per site. 

# attach packages
library(tidyverse)
library(unmarked)
library(lubridate)
library(dplyr)

# Set working directory 
setwd("~/Documents/Chapter 2")

# load occupancy/wildtrax data 
my_data <- read.csv("Data/Raw/SpeciesRawDownload-10/BU_Hart-BTNW-BU-2022_basic_summary.csv")
head(my_data)

# Filter data for only BTNW detections
btnw_data <- subset(my_data, species_code == "BTNW")

# For each site, calculate the detection probability
site_list <- unique(my_data$location)
detection_probs <- numeric(length(site_list))

for (i in 1:length(site_list)) {
  site <- site_list[i]
  
  total_visits <- nrow(subset(my_data, location == site))
  btnw_visits <- nrow(subset(btnw_data, location == site))
  
  detection_probs[i] <- btnw_visits / total_visits
}

# Calculate the average detection probability across all sites
average_detection_prob <- mean(detection_probs)
print(average_detection_prob)


#Using unmarked #####################
# Convert species_code to binary: 1 for BTNW and 0 for NONE
my_data$detected <- ifelse(my_data$species_code == "BTNW", 1, 0)

# Convert to wide format
wide_data <- with(my_data, tapply(detected, list(location, recording_date), sum))
wide_data[wide_data > 1] <- 1  # Ensure binary data, in case of multiple detections in a single visit

# Convert wide data to unmarked frame
umf <- unmarkedFrameOccu(wide_data)

# Fit the occupancy model 
mod <- occu(~1 ~ 1, data=umf)
summary(mod)

# Use inverse logit transformation to convert detection beta to a 
# probability between 0 and 1 
logit_to_prob <- function(logit_val) {
  return(exp(logit_val) / (1 + exp(logit_val)))
}

# Extract logit value from model
logit_val <- coef(mod)[2]

# Convert logit to probability
p <- logit_to_prob(logit_val)

# Determine number of visits required
n <- ceiling(log(1 - P) / log(1 - p))
print(n)






