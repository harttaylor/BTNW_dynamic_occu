# Loading required data and libraries
print("loading data")
load('readytogo.Rdata')
cat(".")
load('readytogo2.Rdata')
cat(".")
library(runjags)
cat(".")
#### MCMC parameters
nchains <- 6 # Use all but five cores for parallel processing
thinning <- 100
niter <- 1e5


### fill the list
datalist <- list(nyrs = nyears, npatches = npatches, nstations = nstapP, nvisits = nvis, detection = detections, dist_p = d_pipe, dist_s = d_seis, dist_h = d_harv, dist_r = d_road, prop_s = pro_seis, prop_h = pro_harv, prop_r = pro_road, prop_p = pro_pipe, treat = treats, sage = stage, pcon = peco, isaroadat=isaroadat)

keep_track_of <- c("alpha_c", "beta_c_pipe_d", "beta_c_pipe_p", "beta_c_seis_p", "beta_c_seis_d", "beta_c_harv_p", "beta_c_harv_d", "beta_c_road_p", "beta_c_road_d", "beta_c_sage", "beta_c_con", "beta_c_treat", "alpha_e", "beta_e_seis_p", "beta_e_seis_d", "beta_e_harv_p", "beta_e_harv_d", "beta_e_pipe_p", "beta_e_pipe_d", "beta_e_road_p", "beta_e_road_d", "beta_e_sage", "beta_e_con", "beta_e_treat", "p", "pe", "pc", "pp")

#### MCMC
cat("############\n")
print(paste("let's begin with",nchains,"chains\n"))
print(Sys.time())
cat("############\n\n")

cat("Running JAGS model\n")
print(Sys.time())

# Define the JAGS model file
model_file <- 'slimmer_ecd_model'

# Running the model using runjags
results <- runjags::run.jags(
  model = model_file,
  monitor = keep_track_of,
  data = datalist,
  n.chains = nchains,
  adapt = niter/10,
  burnin = niter, 
  sample = niter/(thinning), 
  thin = thinning,
  method = "parallel"
)
save.image(version = 2, file = "parecds.Rdata")
cat("Model run completed\n")
print(Sys.time())

# Extracting results

cat("Posterior draws completed\n")
print(Sys.time())

cat("Image Saved\n")


print(Sys.time())
print('all done')

load("parecds.Rdata")
library(coda)
library(rjags)
library(ggplot2)

# Have to convert the mcmc part of results object to an mcmclist 
posterior_draws <- as.mcmc.list(results$mcmc)
str(posterior_draws)
summaries <- summary(posterior_draws)
print(summaries[[2]][1:30,]) #to see something manageable


# Create density plots 
beta_c_pipe_d_samples <- posterior_draws[,"beta_c_pipe_d"]
plot(density(beta_c_pipe_d_samples), main="Density of beta_c_pipe_d", xlab="Value", ylab="Density")
abline(v=0, col="red")  # Adding a line at zero

# Prepare dataframe for plotting
df_combined <- data.frame(
  variable = c("beta_c_pipe_d", "beta_c_pipe_p", "beta_c_seis_p", "beta_c_seis_d", "beta_c_harv_p", "beta_c_harv_d",
               "beta_c_road_p", "beta_c_road_d", "beta_c_sage", "beta_c_con", "beta_c_treat[1]", "beta_c_treat[2]",
               "beta_c_treat[3]", "beta_e_seis_p", "beta_e_seis_d", "beta_e_harv_p", "beta_e_harv_d", "beta_e_pipe_p",
               "beta_e_pipe_d", "beta_e_road_p", "beta_e_road_d", "beta_e_sage", "beta_e_con", "beta_e_treat[1]",
               "beta_e_treat[2]", "beta_e_treat[3]"),
  mean = c(0.38790350, 0.11214000, -0.10804750, -12.44545000, -0.54616400, -0.46380750, -0.35213500, -0.07901020,
           1.42790500, -0.20658150, 0.02346870, -0.25832000, -0.43642150, -0.23963850, 2.04469500, -1.21891500,
           -1.23478500, -0.03710885, 1.51048000, -0.28125600, 0.15509700, 1.83010000, -0.02786150, -0.24178250,
           -0.17263700, -0.66270300),
  lower = c(-7.21624800, -0.12743328, -0.37356315, -15.78730000, -0.86873162, -8.27986900, -0.55273697, -7.88151950,
            0.87764625, -0.27562942, -4.06182025, -4.34683350, -4.49592725, -0.59526307, -2.27484850, -1.77868100,
            -9.06266850, -0.38176120, -5.93966825, -0.59630945, -7.76025550, 1.04894950, -0.13204005, -3.53736850,
            -3.45802050, -3.89691325),
  upper = c(8.02732850, 0.34342722, 0.15000842, -9.10602875, -0.25271668, 7.03490725, -0.15655083, 7.75164650,
            1.98806850, -0.13618098, 3.20188525, 2.90764950, 2.75292425, 0.10549432, 6.29447150, -0.68996248,
            6.89558375, 0.30498310, 9.43651200, 0.02350015, 7.78903450, 2.61832500, 0.07824772, 3.84445100,
            3.93229000, 3.45139625),
  process = rep(c("Colonization", "Extinction"), each=13)
)

hfcolext <- ggplot(df_combined, aes(y = variable, x = mean, xmin = lower, xmax = upper, color = process)) +
  geom_point(size = 2) +
  scale_y_discrete(labels = c("beta_c_pipe_d" = "Proximity to Pipeline (Colonization)", "beta_c_pipe_p" = "Proportion Pipeline (Colonization)", 
                              "beta_c_seis_p" = "Proportion Seismic Line (Colonization)", "beta_c_seis_d" = "Proximity to Seismic Line (Colonization)",
                              "beta_c_harv_p" = "Proportion Harvest (Colonization)", "beta_c_harv_d" = "Proximity to Harvest (Colonization)", 
                              "beta_c_road_p" = "Proportion Road (Colonization)", "beta_c_road_d" = "Proximity to Road (Colonization)",
                              "beta_c_sage" = "Sage (Colonization)", "beta_c_con" = "Con (Colonization)", "beta_c_treat[1]" = "Treatment 1 (Colonization)",
                              "beta_c_treat[2]" = "Treatment 2 (Colonization)", "beta_c_treat[3]" = "Treatment 3 (Colonization)",
                              "beta_e_seis_p" = "Proportion Seismic Line (Extinction)", "beta_e_seis_d" = "Proximity to Seismic Line (Extinction)",
                              "beta_e_harv_p" = "Proportion Harvest (Extinction)", "beta_e_harv_d" = "Proximity to Harvest (Extinction)", 
                              "beta_e_pipe_p" = "Proportion Pipeline (Extinction)", "beta_e_pipe_d" = "Proximity to Pipeline (Extinction)",
                              "beta_e_road_p" = "Proportion Road (Extinction)", "beta_e_road_d" = "Proximity to Road (Extinction)",
                              "beta_e_sage" = "Sage (Extinction)", "beta_e_con" = "Con (Extinction)", "beta_e_treat[1]" = "Treatment 1 (Extinction)",
                              "beta_e_treat[2]" = "Treatment 2 (Extinction)", "beta_e_treat[3]" = "Treatment 3 (Extinction)")) +
  geom_errorbarh(aes(height = 0), linewidth = 1, position = position_dodge(0.2)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Colonization" = "blue", "Extinction" = "red")) +
  theme_bw() +
  labs(x = "Posterior Values", y = "") +
  theme(text = element_text(size = 15, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5)))

print(hfcolext)
ggsave("COLEXT model/Results/spars_colext.png", plot = hfcolext)


