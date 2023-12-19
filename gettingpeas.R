# premise: the model has finished running (or we loaded the R environment post-analyses)
load("COLEXT model/ecdmoar.Rdata")
# we want to plot p over years, regardless of other dimensions (i.e: patch, station, visit)

#so, we make an array of lists "columns per year"
# we do this because an array of lists is an heterogeneous container

cpy<-array(list())

# we need to initialise it to something
for(n in 1:nyears) cpy[[n]]<-0

# we put in all the column names corresponding to year w
# vec contains our posterior column names
vec<-colnames(posterior_draws[[1]])
for(w in 1:nyears) cpy[[w]]<-c(cpy[[w]],grep(paste("^p\\[[^,]*,[^,]*,",w,",[^,]*\\]$",sep=""), vec, value = TRUE))
# remove the initialising valuse
for(w in 1:nyears) cpy[[w]]<-cpy[[w]][-1]


# now we make another container, for the actual values
# "values per year"
vpy<-array(list())
for(w in 1:nyears) vpy[[w]]<-9999

# we loop through the chains
for(q in 1:length(posterior_draws))
	for(w in 1:nyears)
		for(e in 1:length(cpy[[w]]))
			vpy[[w]]<-c(vpy[[w]],posterior_draws[[q]][,which(vec==cpy[[w]][e])])
			
for(w in 1:nyears) vpy[[w]]<-vpy[[w]][-1]

# just to be safe, you can keep it within 2 sd
for(w in 1:nyears) vpy[[w]]<-vpy[[w]][which(abs(vpy[[w]])>(mean(vpy[[w]])+2*sd(vpy[[w]])))]

# Simple plot 
library(ggplot2)
data_for_plotting <- data.frame(
  year = rep(1:nyears, sapply(vpy, length)),
  value = unlist(vpy)
)
treplot <- ggplot(data_for_plotting, aes(x=year, y=value)) + 
  geom_point() + 
  theme_minimal() + 
  labs(x="Year", y="Occupancy Probability")
ggsave("COLEXT model/Results/YearlyOcc.png")

#Average occupancy 
library(dplyr)

data_summary <- data_for_plotting %>%
  group_by(year) %>%
  summarise(
    mean = mean(value),
    lower = mean - qt(0.975, df=n()-1) * sd(value)/sqrt(n()),
    upper = mean + qt(0.975, df=n()-1) * sd(value)/sqrt(n())
  )
treplot <- ggplot(data_summary, aes(x=year, y=mean)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill="blue") +  # 95% CI
  geom_line() +                                         # Trendline
  geom_point() +                                        # Mean points
  theme_minimal() + 
  labs(x="Year (1993-2018)", y="Average Occupancy Probability")

# Try a violin plot 
vioplot <- ggplot(data_for_plotting, aes(factor(year), value)) + 
  geom_violin(trim=FALSE) +
  labs(x="Year (1993-2018)", y="Occupancy Probability") +
  theme_minimal()

# Try adding a trendline to the violin plot 
# Calculating mean for each year
yearly_means <- aggregate(value ~ year, data_for_plotting, mean)

# Adding the trendline to the plot
viotreplot <- ggplot(data_for_plotting, aes(factor(year), value)) + 
  geom_violin(trim=FALSE) +
  geom_point(data = yearly_means, aes(x = factor(year), y = value), color = "blue") +
  geom_smooth(data = yearly_means, aes(x = factor(year), y = value), method = "lm", color = "red") +
  labs(x="Year (1993-2018)", y="Occupancy Probability") +
  theme_minimal()

##----Now do colonization probability----

# Extracting the 'beta_c_seis_d' parameter samples across all chains
beta_c_seis_d_samples <- unlist(lapply(posterior_draws, function(chain) chain[,"beta_c_seis_d"]))

# Assuming the years are sequenced in the order of the data collection
years <- 1:25  # Replace with actual year values if they are different

# Create a dataframe for plotting
df <- data.frame(
  Year = rep(years, each = length(beta_c_seis_d_samples) / length(years)),
  Beta_c_seis_d = beta_c_seis_d_samples
)

# Calculate summary statistics for each year
summary_df <- aggregate(Beta_c_seis_d ~ Year, data = df, function(x) c(Mean = mean(x), SD = sd(x)))

# Reshape for easier plotting with ggplot2
summary_df <- do.call(data.frame, summary_df)
names(summary_df)[names(summary_df) == "Beta_c_seis_d.Mean"] <- "Mean"
names(summary_df)[names(summary_df) == "Beta_c_seis_d.SD"] <- "SD"

# Create the plot
ggplot(summary_df, aes(x = Year, y = Mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2) +
  labs(title = "Effect of Distance to Seismic Lines on Colonization Probability",
       x = "Year",
       y = "Beta_c_seis_d (Mean Effect Size)") +
  theme_minimal()


