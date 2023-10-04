#################################################
#
# Code for a paper on monitoring sites selection bias
# The code below performs the simulations and makes the plots
#
# monitoring_sites_selection_bias_Simulations_and_Figures.R
#
#
#################################################

# clean working space
rm(list = ls())

# load packages
library(Rfast)
library(tidyr)
library(lme4)
library(ggplot2)
library(tidyverse)
library(viridis)
library(ggpubr)


########## Simulations to reproduce the results of the article to which we are replying ##########

# simulation settings
nsites = 100 ; nsurveys = 1 ; nyears = 20 ; mean.beta = 0 ; sd.lam = c(0, 0)
df_simul_settings <- expand.grid(my_rate = c(0.0005, 0.005, 0.05), # distance of abundance outliers
                                 my_temp_sd = c(10, 100, 300)) # severity of yearly fluctuations
nsim = 1000
year <- (1:nyears) - 1
df_simul_results <- data.frame(matrix(NA, nrow = nsim, ncol = nrow(df_simul_settings)))
colnames(df_simul_results) <- paste0("case_", 1:nrow(df_simul_settings))

# store start time
(start.time_total <- Sys.time())

# set random seed
set.seed(20231003)

# run the simulations
for (j in 1:nrow(df_simul_settings)) {
  rate <- df_simul_settings$my_rate[j]
  sd.N <- df_simul_settings$my_temp_sd[j]
  
  for (i in 1:nsim) {
    
    alpha <- rgamma(n = nsites, shape = 0.1, rate = rate) # rate is the parameter to adjust the distance of outliers
    beta <- rnorm(n = nsites, mean = mean.beta, sd = sd.lam[2]) # population trend, set to zero with no variability
    lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears)) # create a matrix with nsites rows and nyears columns for all parameters
    
    # fill the matrix with the expected abundance
    for (t in (0:nyears)) { 
      lam[, t] <- exp(log(alpha) + beta * year[t])
    }
    
    # fill the matrix with the realised abundance
    for (t in (1:nyears)) { 
      N[, t] <- rnorm(n = nsites, mean = lam[, t], sd = sd.N)
    }
    
    # set negative abundances to zero (as did the authors)
    N <- ifelse(N<0, 0, N)
    
    # select the 10 sites with the greatest initial abundance
    vect_10_selected_sites <- c(Rfast::nth(x = N[, 1], k = 10, num.of.nths = 10, descending = TRUE, index.return = TRUE))
    N_selected <- N[vect_10_selected_sites, ]
    
    # convert simulated data from matrix to dataframe format
    df_counts <- as.data.frame.array(N_selected)
    df_counts <- cbind(site = as.factor(c(1:10)), df_counts) # add column with site id
    df_counts <- pivot_longer(df_counts, !site, names_to = "year", values_to = "count", names_prefix = "V") # convert the dataframe to long format
    df_counts$year <- as.integer(df_counts$year)
    
    # fit the linear model
    my_model <- lmer(log(df_counts$count+1) ~ year + (1|site), data = df_counts)
    
    # compute the apparent_pop_change
    estimated_beta <- coef(my_model)$site$year[1]
    apparent_pop_change <- exp(estimated_beta*19)
    df_simul_results[i, j] <- apparent_pop_change
    
  } # i
} # j

# print end time and total execution time
(end.time_total <- Sys.time())
(time.taken_total <- difftime(end.time_total, start.time_total, units = "mins"))
colMeans(df_simul_results) # mean estimated population change across all simulations

# convert the table to long format 
df_all_results_reproduced <- rbind(data.frame(rate = df_simul_settings$my_rate[1], temp_sd = df_simul_settings$my_temp_sd[1], estim_pop_change = df_simul_results$case_1),
                        data.frame(rate = df_simul_settings$my_rate[2], temp_sd = df_simul_settings$my_temp_sd[2], estim_pop_change = df_simul_results$case_2),
                        data.frame(rate = df_simul_settings$my_rate[3], temp_sd = df_simul_settings$my_temp_sd[3], estim_pop_change = df_simul_results$case_3),
                        data.frame(rate = df_simul_settings$my_rate[4], temp_sd = df_simul_settings$my_temp_sd[4], estim_pop_change = df_simul_results$case_4),
                        data.frame(rate = df_simul_settings$my_rate[5], temp_sd = df_simul_settings$my_temp_sd[5], estim_pop_change = df_simul_results$case_5),
                        data.frame(rate = df_simul_settings$my_rate[6], temp_sd = df_simul_settings$my_temp_sd[6], estim_pop_change = df_simul_results$case_6),
                        data.frame(rate = df_simul_settings$my_rate[7], temp_sd = df_simul_settings$my_temp_sd[7], estim_pop_change = df_simul_results$case_7),
                        data.frame(rate = df_simul_settings$my_rate[8], temp_sd = df_simul_settings$my_temp_sd[8], estim_pop_change = df_simul_results$case_8),
                        data.frame(rate = df_simul_settings$my_rate[9], temp_sd = df_simul_settings$my_temp_sd[9], estim_pop_change = df_simul_results$case_9))
df_all_results_reproduced$rate <- as.factor(df_all_results_reproduced$rate)
df_all_results_reproduced$temp_sd <- as.factor(df_all_results_reproduced$temp_sd)

# save simulation output
write.csv2(df_all_results_reproduced, file = "./output/df_all_results_reproducing_the_authors_results_nsim1000.csv", row.names = FALSE)

# plot the results
plot_reproducing_the_authors_results <- ggplot(df_all_results_reproduced, aes(x = rate, y = estim_pop_change)) + 
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1), fill = "lightgrey") +
  geom_boxplot(fill = NA) +
  facet_grid(. ~ temp_sd, labeller = label_both) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Rate Parameter") +
  ylab("Apparent Population Change")

# save the plot
pdf(file = "./output/Fig0_trend_estimates_reproducing_the_authors_results.pdf", width = 10, height = 7)
print(plot_reproducing_the_authors_results)
dev.off()



########## Simulations with temporal variance depending on site abundance ##########

# simulation settings
nsites = 100 ; nsurveys = 1 ; nyears = 20 ; mean.beta = 0 ; sd.lam = c(0, 0)
my_rate = c(0.0005, 0.005, 0.05)
nsim = 1000
year <- (1:nyears) - 1
df_simul_results <- data.frame(matrix(NA, nrow = nsim, ncol = length(my_rate)))
colnames(df_simul_results) <- paste0("case_", 1:length(my_rate))

# store start time
(start.time_total <- Sys.time())

# set random seed
set.seed(20231002)

# run the simulations
for (j in 1:length(my_rate)) {
  rate <- my_rate[j]

  for (i in 1:nsim) {
    
    alpha <- rgamma(n = nsites, shape = 0.1, rate = rate) # rate is the parameter to adjust the distance of outliers
    beta <- rnorm(n = nsites, mean = mean.beta, sd = sd.lam[2]) # population trend, set to zero with no variability
    lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears)) # create a matrix with nsites rows and nyears columns for all parameters
    
    # fill the matrix with the expected abundance
    for (t in (0:nyears)) { 
      lam[, t] <- exp(log(alpha) + beta * year[t])
    }
    
    # fill the matrix with the realised abundance
    for (t in (1:nyears)) { 
      N[, t] <- rnorm(n = nsites, mean = lam[, t], sd = lam[, t]) # temporal SD equal to expected site abundance
    }
    
    # set negative abundances to zero (as did the authors)
    N <- ifelse(N<0, 0, N)
    
    # select the 10 sites with the greatest initial abundance
    vect_10_selected_sites <- c(Rfast::nth(x = N[, 1], k = 10, num.of.nths = 10, descending = TRUE, index.return = TRUE))
    N_selected <- N[vect_10_selected_sites, ]
    
    # convert simulated data from matrix to dataframe format
    df_counts <- as.data.frame.array(N_selected)
    df_counts <- cbind(site = as.factor(c(1:10)), df_counts) # add column with site id
    df_counts <- pivot_longer(df_counts, !site, names_to = "year", values_to = "count", names_prefix = "V") # convert the dataframe to long format
    df_counts$year <- as.integer(df_counts$year)
    
    # fit the linear model
    my_model <- lmer(log(df_counts$count+1) ~ year + (1|site), data = df_counts)
    
    # compute the apparent_pop_change
    estimated_beta <- coef(my_model)$site$year[1]
    apparent_pop_change <- exp(estimated_beta*19)
    df_simul_results[i, j] <- apparent_pop_change
    
  } # i
} # j

# print end time and total execution time
(end.time_total <- Sys.time())
(time.taken_total <- difftime(end.time_total, start.time_total, units = "mins"))
colMeans(df_simul_results) # mean estimated population change across all simulations

# convert the table to long format 
df_all_results_var_dep_abundance <- rbind(data.frame(rate = my_rate[1], estim_pop_change = df_simul_results$case_1),
                        data.frame(rate = my_rate[2], estim_pop_change = df_simul_results$case_2),
                        data.frame(rate = my_rate[3], estim_pop_change = df_simul_results$case_3))
df_all_results_var_dep_abundance$rate <- as.factor(df_all_results_var_dep_abundance$rate)

# save simulation output
write.csv2(df_all_results_var_dep_abundance, file = "./output/df_all_results_variance_depending_on_abundance_nsim1000.csv", row.names = FALSE)

# plot the results
plot_results_variance_depending_on_abundance <- ggplot(df_all_results_var_dep_abundance, aes(x = rate, y = estim_pop_change)) + 
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1), fill = "lightgrey") +
  geom_boxplot(fill = NA) +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Rate Parameter") +
  ylab("Apparent Population Change")

# save the plot
pdf(file = "./output/Fig1_trend_estimates_variance_depending_on_abundance.pdf", width = 5, height = 7)
print(plot_results_variance_depending_on_abundance)
dev.off()



########## Simulations with temporal variance depending on site abundance WITH A dependence GRADIENT ##########

# simulation settings
nsites = 100 ; nsurveys = 1 ; nyears = 20 ; mean.beta = 0 ; sd.lam = c(0, 0)
rate = 0.0005
my_dep = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1) # coefficient between temporal variance and abundance
nsim = 1000
year <- (1:nyears) - 1
df_simul_results <- data.frame(matrix(NA, nrow = nsim, ncol = length(my_dep)))
colnames(df_simul_results) <- paste0("case_", 1:length(my_dep))

# store start time
(start.time_total <- Sys.time())

# set random seed
set.seed(20231001)

# run the simulations
for (j in 1:length(my_dep)) {
  dep_ratio <- my_dep[j]
  
  for (i in 1:nsim) {
    
    alpha <- rgamma(n = nsites, shape = 0.1, rate = rate) # rate is the parameter to adjust the distance of outliers
    beta <- rnorm(n = nsites, mean = mean.beta, sd = sd.lam[2]) # population trend, set to zero with no variability
    lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears)) # create a matrix with nsites rows and nyears columns for all parameters
    
    # fill the matrix with the expected abundance
    for (t in (0:nyears)) { 
      lam[, t] <- exp(log(alpha) + beta * year[t])
    }
    
    # fill the matrix with the realised abundance
    for (t in (1:nyears)) { 
      N[, t] <- rnorm(n = nsites, mean = lam[, t], sd = dep_ratio*lam[, t]) # temporal SD proportional to expected site abundance
    }
    
    # set negative abundances to zero (as did the authors)
    N <- ifelse(N<0, 0, N)
    
    # select the 10 sites with the greatest initial abundance
    vect_10_selected_sites <- c(Rfast::nth(x = N[, 1], k = 10, num.of.nths = 10, descending = TRUE, index.return = TRUE))
    N_selected <- N[vect_10_selected_sites, ]
    
    # convert simulated data from matrix to dataframe format
    df_counts <- as.data.frame.array(N_selected)
    df_counts <- cbind(site = as.factor(c(1:10)), df_counts) # add column with site id
    df_counts <- pivot_longer(df_counts, !site, names_to = "year", values_to = "count", names_prefix = "V") # convert the dataframe to long format
    df_counts$year <- as.integer(df_counts$year)
    
    # fit the linear model
    my_model <- lmer(log(df_counts$count+1) ~ year + (1|site), data = df_counts)
    
    # compute the apparent_pop_change
    estimated_beta <- coef(my_model)$site$year[1]
    apparent_pop_change <- exp(estimated_beta*19)
    df_simul_results[i, j] <- apparent_pop_change
    
  } # i
} # j

# print end time and total execution time
(end.time_total <- Sys.time())
(time.taken_total <- difftime(end.time_total, start.time_total, units = "mins"))
colMeans(df_simul_results) # mean estimated population change across all simulations


# convert the table to long format 
df_all_results_var_dependence_gradient <- rbind(data.frame(var_dependence = my_dep[1], estim_pop_change = df_simul_results$case_1),
                        data.frame(var_dependence = my_dep[2], estim_pop_change = df_simul_results$case_2),
                        data.frame(var_dependence = my_dep[3], estim_pop_change = df_simul_results$case_3),
                        data.frame(var_dependence = my_dep[4], estim_pop_change = df_simul_results$case_4),
                        data.frame(var_dependence = my_dep[5], estim_pop_change = df_simul_results$case_5),
                        data.frame(var_dependence = my_dep[6], estim_pop_change = df_simul_results$case_6),
                        data.frame(var_dependence = my_dep[7], estim_pop_change = df_simul_results$case_7),
                        data.frame(var_dependence = my_dep[8], estim_pop_change = df_simul_results$case_8),
                        data.frame(var_dependence = my_dep[9], estim_pop_change = df_simul_results$case_9),
                        data.frame(var_dependence = my_dep[10], estim_pop_change = df_simul_results$case_10))
df_all_results_var_dependence_gradient$var_dependence <- as.factor(df_all_results_var_dependence_gradient$var_dependence)

# save simulation output
write.csv2(df_all_results_var_dependence_gradient, file = "./output/df_all_results_var_dependence_gradient_nsim1000.csv", row.names = FALSE)

# plot the results
plot_results_var_dependence_gradient <- ggplot(df_all_results_var_dependence_gradient, aes(x = var_dependence, y = estim_pop_change)) + 
  geom_rect(data = NULL, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 1), fill = "lightgrey") +
  geom_boxplot(fill = NA) +
  theme_classic() +
  theme(text = element_text(size = 19)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Coefficient between temporal variance and abundance") +
  ylab("Apparent Population Change")

# save the plot
pdf(file = "./output/Fig2_trend_estimates_variance_dependence_gradient.pdf", width = 10, height = 7)
print(plot_results_var_dependence_gradient)
dev.off()



########## Make plots to illustrate the yearly abundance fluctuations we simulated ##########

#### Figure 1 : temporal variance depending on site abundance ####

# simulate data the same way as we did above but store the realised abundances
nsites = 100 ; nsurveys = 1 ; nyears = 20 ; mean.beta = 0 ; sd.lam = c(0, 0)
my_rate = c(0.0005, 0.005, 0.05)
nsim = 1000
year <- (1:nyears) - 1
saved_N <- list(NA, NA, NA)

# set random seed
set.seed(123)

# run the simulations
for (j in 1:length(my_rate)) {
  rate <- my_rate[j]
  
  alpha <- rgamma(n = nsites, shape = 0.1, rate = rate) # rate is the parameter to adjust the distance of outliers
  beta <- rnorm(n = nsites, mean = mean.beta, sd = sd.lam[2]) # population trend, set to zero with no variability
  lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears)) # create a matrix with nsites rows and nyears columns for all parameters
  
  # fill the matrix with the expected abundance
  for (t in (0:nyears)) { 
    lam[, t] <- exp(log(alpha) + beta * year[t])
  }
  
  # fill the matrix with the realised abundance
  for (t in (1:nyears)) { 
    N[, t] <- rnorm(n = nsites, mean = lam[, t], sd = lam[, t]) # temporal SD equal to expected site abundance
  }
  
  # set negative abundances to zero (as did the authors)
  N <- ifelse(N<0, 0, N)
  
  # save the matrix with the realised abundance for plotting
  saved_N[[j]] <- N
  
} # j


# format data for plotting
format_my_data <- function(N_matrix) {
  df_N <- as.data.frame.matrix(N_matrix)
  df_N <- df_N[order(df_N$V1, decreasing = TRUE), ]
  colnames(df_N) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
  df_N <- cbind(site_ID = c(1:100), site_status = c(rep("selected", times = 10), rep("not_selected", times = 90)), df_N)
  
  # convert to long format
  df_N_long <- df_N %>%
    pivot_longer(cols = -c(site_ID, site_status), names_to = "year", values_to = "abundance")
  
  # set right column types
  df_N_long$site_ID <- as.factor(df_N_long$site_ID)
  df_N_long$site_status <- as.factor(df_N_long$site_status)
  df_N_long$year <- as.integer(df_N_long$year)
  df_N_long$abundance <- as.double(df_N_long$abundance)
  
  return(df_N_long)
}

df_abundances_0005 <- format_my_data(N_matrix = saved_N[[1]])
df_abundances_005 <- format_my_data(N_matrix = saved_N[[2]])
df_abundances_05 <- format_my_data(N_matrix = saved_N[[3]])


# make the plots
my_palette <- viridis(n = 10)

plot_rate_0005 <- ggplot(df_abundances_0005, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(color = "lightgrey", size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "6"), color = my_palette[6], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "7"), color = my_palette[7], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "8"), color = my_palette[8], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "9"), color = my_palette[9], size = 1) +
  geom_line(data = subset(df_abundances_0005, df_abundances_0005$site_ID == "10"), color = my_palette[10], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Rate = 0.0005")

plot_rate_005 <- ggplot(df_abundances_005, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(color = "lightgrey", size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "6"), color = my_palette[6], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "7"), color = my_palette[7], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "8"), color = my_palette[8], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "9"), color = my_palette[9], size = 1) +
  geom_line(data = subset(df_abundances_005, df_abundances_005$site_ID == "10"), color = my_palette[10], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Rate = 0.005")

plot_rate_05 <- ggplot(df_abundances_05, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(color = "lightgrey", size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "6"), color = my_palette[6], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "7"), color = my_palette[7], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "8"), color = my_palette[8], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "9"), color = my_palette[9], size = 1) +
  geom_line(data = subset(df_abundances_05, df_abundances_05$site_ID == "10"), color = my_palette[10], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Rate = 0.05")

# arrange the plots together to make the figure
Figure1 <- ggarrange(plot_rate_0005, plot_rate_005, plot_rate_05,
                     plot_results_variance_depending_on_abundance,
                     labels = c("(a)", "(b)", "(c)", "(d)"),
                     ncol = 2, nrow = 2)

# save the plot
pdf(file = "./output/Figure1.pdf", width = 10, height = 8)
print(Figure1)
dev.off()
png(file = "./output/Figure1.png", width = 10, height = 8, units = "in", res = 600)
print(Figure1)
dev.off()



#### Figure 2 : variance dependence gradient ####

# simulate data the same way as we did above but store the realised abundances
# simulation settings
nsites = 100 ; nsurveys = 1 ; nyears = 20 ; mean.beta = 0 ; sd.lam = c(0, 0)
rate = 0.0005
my_dep = c(0.2, 0.5, 0.8) # three coefficients to illustrate
nsim = 1000
year <- (1:nyears) - 1
saved_N_dependence <- list(NA, NA, NA)

# set random seed
set.seed(12)

# run the simulations
for (j in 1:length(my_dep)) {
  dep_ratio <- my_dep[j]
  
  alpha <- rgamma(n = nsites, shape = 0.1, rate = rate) # rate is the parameter to adjust the distance of outliers
  beta <- rnorm(n = nsites, mean = mean.beta, sd = sd.lam[2]) # population trend, set to zero with no variability
  lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears)) # create a matrix with nsites rows and nyears columns for all parameters
  
  # fill the matrix with the expected abundance
  for (t in (0:nyears)) { 
    lam[, t] <- exp(log(alpha) + beta * year[t])
  }
  
  # fill the matrix with the realised abundance
  for (t in (1:nyears)) { 
    N[, t] <- rnorm(n = nsites, mean = lam[, t], sd = dep_ratio*lam[, t]) # temporal SD proportional to expected site abundance
  }
  
  # set negative abundances to zero (as did the authors)
  N <- ifelse(N<0, 0, N)
  
  # save the matrix with the realised abundance for plotting
  saved_N_dependence[[j]] <- N
  
} # j


# format data for plotting
df_abundances_dep02 <- format_my_data(N_matrix = saved_N_dependence[[1]])
df_abundances_dep05 <- format_my_data(N_matrix = saved_N_dependence[[2]])
df_abundances_dep08 <- format_my_data(N_matrix = saved_N_dependence[[3]])

# make the plots
my_palette <- viridis(n = 10)

plot_dep02 <- ggplot(df_abundances_dep02, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(color = "lightgrey", size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "6"), color = my_palette[6], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "7"), color = my_palette[7], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "8"), color = my_palette[8], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "9"), color = my_palette[9], size = 1) +
  geom_line(data = subset(df_abundances_dep02, df_abundances_dep02$site_ID == "10"), color = my_palette[10], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 19), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Coefficient = 0.2")

plot_dep05 <- ggplot(df_abundances_dep05, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(color = "lightgrey", size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "6"), color = my_palette[6], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "7"), color = my_palette[7], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "8"), color = my_palette[8], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "9"), color = my_palette[9], size = 1) +
  geom_line(data = subset(df_abundances_dep05, df_abundances_dep05$site_ID == "10"), color = my_palette[10], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 19), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Coefficient = 0.5")

plot_dep08 <- ggplot(df_abundances_dep08, aes(x = year, y = abundance, group = site_ID)) +
  geom_line(color = "lightgrey", size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "1"), color = my_palette[1], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "2"), color = my_palette[2], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "3"), color = my_palette[3], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "4"), color = my_palette[4], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "5"), color = my_palette[5], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "6"), color = my_palette[6], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "7"), color = my_palette[7], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "8"), color = my_palette[8], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "9"), color = my_palette[9], size = 1) +
  geom_line(data = subset(df_abundances_dep08, df_abundances_dep08$site_ID == "10"), color = my_palette[10], size = 1) +
  theme_classic() +
  theme(text = element_text(size = 19), plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Abundance") +
  ggtitle("Coefficient = 0.8")


# arrange the plots together to make the figure
Figure2 <- ggarrange(ggarrange(plot_dep02, plot_dep05, plot_dep08, ncol = 3,
                               labels = c("(a)", "(b)", "(c)"), font.label = list(size = 19)),
                     plot_results_var_dependence_gradient,
                     nrow = 2, labels = c(" ", "(d)"), font.label = list(size = 19))

# save the plot
pdf(file = "./output/Figure2.pdf", width = 13, height = 10)
print(Figure2)
dev.off()
png(file = "./output/Figure2.png", width = 13, height = 10, units = "in", res = 600)
print(Figure2)
dev.off()


