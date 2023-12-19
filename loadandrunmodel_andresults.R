load('readytogo.Rdata')
library(rjags)

#### MCMC parameters
nchains<-2
thinning<-90
niter<-27e4

### fill the list

datalist<-list(nyrs=nyears,npatches=npatches,nstations=nstapP,nvisits=nvis,detection=detections,dist_p=d_pipe,dist_s=d_seis,dist_h=d_harv,dist_w=d_well,prop_s=pro_seis,prop_h=pro_harv,prop_p=pro_pipe,prop_w=pro_well,treat=treats,sage=stage,pcon=peco)

keep_track_of<-c("alpha_c","beta_c_pipe_d","beta_c_pipe_p","beta_c_seis_p","beta_c_seis_d","beta_c_harv_p","beta_c_harv_d","beta_c_well_p","beta_c_well_d","beta_c_sage","beta_c_con","beta_c_patch","beta_c_treat","alpha_e","beta_e_seis_p","beta_e_seis_d","beta_e_harv_p","beta_e_harv_d","beta_e_pipe_p","beta_e_pipe_d","beta_e_well_p","beta_e_well_d","beta_e_sage","beta_e_con","beta_e_patch","beta_e_treat","p","pe","pc","pp")


#### declaring the model

print("declaring the model on")
print(Sys.time())
mod<-jags.model(n.chains=nchains,file='slimmer_ecd_model',data=datalist)

print("Starting Burn-in")
print(Sys.time())
#### burn-in

update(mod,niter)

#### MCMC
print("starting posterior draws")
print(Sys.time())
posterior_draws<-coda.samples(model=mod,variable.names=keep_track_of,thin=thinning,n.iter=niter)
print("finished")
print(Sys.time())
smod<-summary(posterior_draws)

save.image(version=2,file="sled.Rdata")

# ---- Plot Results ----
# Look at model results and make plots 
results <- load("ecdmoar.Rdata")
ls()
smod<-summary(posterior_draws)
plot(smod)
str(smod)
head(smod)
traceplot(posterior_draws)

# Load the coda package
library(coda)
# Extract the parameter from each chain and combine
mcmc_chain = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_harv_d"]))

# Confirm it's numeric
print(is.numeric(mcmc_chain))

# Plot the posterior distribution
plot(density(mcmc_chain), main="Posterior Distribution of beta_c_harv_d", xlab="beta_c_harv_d", ylab="Density")

# ---- Human Footprint ----
#create violin plot to show effects of footprints 
install.packages("ggplot2")
library(ggplot2)
summary(posterior_draws)
beta_harv_d = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_harv_d"]))
#beta_road_d = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_road_d"])) 
beta_pipe_d = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_pipe_d"])) 
beta_seis_d = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_seis_d"])) 

beta_harv_p = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_harv_p"]))
#beta_road_p = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_road_p"])) 
beta_pipe_p = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_pipe_p"])) 
beta_seis_p = unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_seis_p"])) 

# Assuming all beta vectors are of the same length
data_for_plot = data.frame(
  value = c(beta_harv_d, beta_pipe_d, beta_seis_d, 
            beta_harv_p, beta_pipe_p, beta_seis_p),
  variable = factor(rep(c("Harv_d", "Pipe_d", "Seis_d", 
                          "Harv_p", "Pipe_p", "Seis_p"), 
                        each = length(beta_harv_d)))
)
ggplot(data_for_plot, aes(x = variable, y = value)) + 
  geom_violin(trim = FALSE) + 
  labs(title = "Violin Plot of Beta Parameters", x = "Parameter", y = "Value")

# now create plot showing posterior distributions of variables of interest like is automatically done in stan 
# Combining the chains into one mcmc.list object
combined_chains = mcmc.list(posterior_draws)
mcmc_summary = summary(combined_chains)

# Extract mean and quantiles for each parameter
beta_stats = sapply(c("beta_c_harv_d", "beta_c_pipe_d", "beta_c_seis_d", "beta_c_harv_p", "beta_c_pipe_p", "beta_c_seis_p"), function(param) {
  c(mean = mcmc_summary$statistics[param, "Mean"], 
    lower = mcmc_summary$quantiles[param, "2.5%"], 
    upper = mcmc_summary$quantiles[param, "97.5%"])
}, simplify = "data.frame")

# Reshape the beta_stats matrix into a dataframe
df = data.frame(
  variable = colnames(beta_stats),
  mean = beta_stats["mean", ], 
  lower = beta_stats["lower", ], 
  upper = beta_stats["upper", ]
)

library(ggplot2)

# Assuming df is your prepared dataframe with variables and their posterior stats
ggplot(df, aes(y = variable, x = mean, xmin = lower, xmax = upper)) +
  geom_point(color = "red", size = 2) +
  geom_errorbarh(aes(height = 0), color = "blue", size = 1, position = position_dodge(0.2)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  theme_minimal() +
  labs(title = "Posterior Distributions of Footprint Variables and Colonization", x = "Posterior values", y = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank())

# Now make plots for both col and ext 
# Define the parameter names for colonization and extinction
params_colonization = c("beta_c_harv_d", "beta_c_pipe_d", "beta_c_seis_d", "beta_c_harv_p", "beta_c_pipe_p", "beta_c_seis_p")
params_extinction = c("beta_e_harv_d", "beta_e_pipe_d", "beta_e_seis_d", "beta_e_harv_p", "beta_e_pipe_p", "beta_e_seis_p")

# Extract mean and quantiles for each parameter
extract_stats = function(param_name) {
  c(mean = mcmc_summary$statistics[param_name, "Mean"], 
    lower = mcmc_summary$quantiles[param_name, "2.5%"], 
    upper = mcmc_summary$quantiles[param_name, "97.5%"])
}

# Apply the function to both sets of parameters and transpose the matrix
beta_stats_colonization = t(sapply(params_colonization, extract_stats))
beta_stats_extinction = t(sapply(params_extinction, extract_stats))

# Create data frames from the transposed matrices
df_colonization = data.frame(variable = rownames(beta_stats_colonization), beta_stats_colonization)
df_extinction = data.frame(variable = rownames(beta_stats_extinction), beta_stats_extinction)

# Since we now include the variable names as the first column, we set the column names accordingly
names(df_colonization)[2:4] = c("mean", "lower", "upper")
names(df_extinction)[2:4] = c("mean", "lower", "upper")

# Add the process column
df_colonization$process = "Colonization"
df_extinction$process = "Extinction"

# Combine the data frames
df_combined = rbind(df_colonization, df_extinction)

# Assuming df_combined is your prepared dataframe with variables, their posterior stats, and process type
hfcolext <- ggplot(df_combined, aes(y = variable, x = mean, xmin = lower, xmax = upper, color = process)) +
  geom_point(size = 2) +
  scale_y_discrete(labels = c("beta_c_harv_d" = "Proximity to Harvest", "beta_c_pipe_d" = "Proximity to Pipeline", "beta_c_seis_d" = "Proximity to Seismic Line",
                              "beta_c_harv_p" = "Proportion Harvest", "beta_c_pipe_p" = "Proportion Pipeline", "beta_c_seis_p" = "Proportion Seismic Line",
                              "beta_e_harv_d" = "Proximity to Harvest", "beta_e_pipe_d" = "Proximity to Pipeline", "beta_e_seis_d" = "Proximity to Seismic Line",
                              "beta_e_harv_p" = "Proportion Harvest", "beta_e_pipe_p" = "Proportion Pipeline", "beta_e_seis_p" = "Proportion Seismic Line")) +
  geom_errorbarh(aes(height = 0), size = 1, position = position_dodge(0.2)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Colonization" = "blue", "Extinction" = "red")) +
  theme_bw() +  # White background with grids
  labs(x = "Posterior Values", y = "") +
  theme(text = element_text(size = 15, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()) +  # Remove minor gridlines
      #  panel.border = element_blank()) +  # Remove panel border
   guides(color = guide_legend(override.aes = list(size = 5)))

print(hfcolext)
# Now save the plot
ggsave("COLEXT model/Results/HFcolext.png", plot = hfcolext, width = 10, height = 8, dpi = 300)


# ---- Treatment ----
# Make plots of the treatment parameters posterior distribution 
# for both col and ext 
# Define the parameter names for colonization and extinction
params_colonization = c("beta_c_treat[1]", "beta_c_treat[2]", "beta_c_treat[3]")
params_extinction = c("beta_e_treat[1]", "beta_e_treat[2]", "beta_e_treat[3]")

# Extract mean and quantiles for each parameter
extract_stats = function(param_name) {
  c(mean = mcmc_summary$statistics[param_name, "Mean"], 
    lower = mcmc_summary$quantiles[param_name, "2.5%"], 
    upper = mcmc_summary$quantiles[param_name, "97.5%"])
}

# Apply the function to both sets of parameters and transpose the matrix
beta_stats_colonization = t(sapply(params_colonization, extract_stats))
beta_stats_extinction = t(sapply(params_extinction, extract_stats))

# Create data frames from the transposed matrices
df_colonization = data.frame(variable = rownames(beta_stats_colonization), beta_stats_colonization)
df_extinction = data.frame(variable = rownames(beta_stats_extinction), beta_stats_extinction)

# Since we now include the variable names as the first column, we set the column names accordingly
names(df_colonization)[2:4] = c("mean", "lower", "upper")
names(df_extinction)[2:4] = c("mean", "lower", "upper")

# Add the process column
df_colonization$process = "Colonization"
df_extinction$process = "Extinction"

# Combine the data frames
df_combined = rbind(df_colonization, df_extinction)

# Assuming df_combined is your prepared dataframe with variables, their posterior stats, and process type
trtcolext <- ggplot(df_combined, aes(y = variable, x = mean, xmin = lower, xmax = upper, color = process)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(height = 0), size = 1, position = position_dodge(0.2)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Colonization" = "blue", "Extinction" = "red")) +
  theme_bw() +  # White background with grids
  labs(title = "Effects on Colonization and Extinction Probabilities", x = "Posterior values", y = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank()) +  # Remove minor gridlines
  #  panel.border = element_blank()) +  # Remove panel border
  guides(color = guide_legend(override.aes = list(size = 5)))

print(trtcolext)
# Now save the plot
ggsave("COLEXT model/Results/trtcolext.png", plot = trtcolext, width = 10, height = 8, dpi = 300)


# Now do individual footprints and different col/ext parameters 
# Extract the beta_c_seis_d parameter from each chain and combine
# Assuming you have the mean of the 'beta_c_seis_d' parameter from your posterior draws
mean_beta_c_seis_d <- mean(unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_seis_d"])))

# Assuming 'alpha_c' is also a parameter in your model, representing the intercept
mean_alpha_c <- mean(unlist(lapply(posterior_draws, function(chain) chain[,"alpha_c"])))

# Define a range of proximity to seismic lines values (you might want to adjust these values based on your data)
proximity_to_seismic_lines <- seq(-3, 3, length.out = 100) # Assuming this is a standardized scale

# Calculate colonization probabilities using the logistic function
colonization_probability <- 1 / (1 + exp(-(mean_alpha_c + mean_beta_c_seis_d * proximity_to_seismic_lines)))

# Create a data frame for plotting
plot_data <- data.frame(
  ProximityToSeismicLines = proximity_to_seismic_lines,
  ColonizationProbability = colonization_probability
)

# Plot
ggplot(plot_data, aes(x = ProximityToSeismicLines, y = ColonizationProbability)) +
  geom_line() +
  labs(title = "Effect of Proximity to Seismic Lines on Colonization Probability",
       x = "Proximity to Seismic Lines (Standardized)",
       y = "Probability of Colonization") +
  theme_minimal()

# Try using a different way 
# Combine chains into one mcmc.list object for easier manipulation
combined_chains <- mcmc.list(posterior_draws)

# Extract the mean of the beta parameter for proximity to seismic lines
mean_beta_c_seis_d <- summary(combined_chains)$statistics["beta_c_seis_d", "Mean"]

# Assuming 'alpha_c' is also a parameter in your model representing the intercept
mean_alpha_c <- summary(combined_chains)$statistics["alpha_c", "Mean"]

# Define a sequence of standardized proximity values based on your data
standardized_proximity <- seq(min(asc$NEAR.DIST.conventional.seismic, na.rm = TRUE),
                              max(asc$NEAR.DIST.conventional.seismic, na.rm = TRUE), length.out = 100)

# Transform back to the original scale if necessary (asc data must be loaded and contain the proximity variable)
original_proximity <- standardized_proximity * max(asc$NEAR.DIST.conventional.seismic, na.rm = TRUE)

# Calculate colonization probabilities using the logistic function
colonization_probability <- 1 / (1 + exp(-(mean_alpha_c + mean_beta_c_seis_d * standardized_proximity)))

# Create a data frame for plotting
plot_data <- data.frame(
  ProximityToSeismicLine = original_proximity,
  ColonizationProbability = colonization_probability
)

# Plot
ggplot(plot_data, aes(x = ProximityToSeismicLine, y = ColonizationProbability)) +
  geom_line() +
  labs(title = "Effect of Proximity to Seismic Lines on Colonization Probability",
       x = "Proximity to Seismic Lines (meters)",
       y = "Probability of Colonization") +
  theme_minimal()





