# dependencies
library(rethinking)
library(viridis)
library(testthat)
library(tidyverse)

# create folder to contain plots
dir.create("Figures")

# load in data if necessary
# simulated data
d <- read.csv("d_sim.csv", stringsAsFactors = FALSE)

# when working with real data:
# df used for model
#d <- read.csv("s_person_year_df.csv", stringsAsFactors = FALSE)

# df of birth-death lifecourses
#s <- read.csv(file = "w_s.csv", stringsAsFactors = FALSE)

# load in samples if necessary
# load("post_pois.RData")

#-----------------------------------------------------------------------------------------------------------------------
### plotting total residential moves over the lifetime (only possible with HSN data)

# check first reg address is at age 0
# expect_true(nrow(s[which((s$mun_move_category == "first") & (s$age_at_move != 0)),]) == 0)   # first reg address is always at age 0
# 
# ss <- s[which((s$death_y - s$birth_y) >= 20),]
# 
# ss %>%
#   group_by(person_id) %>%
#   summarize(
#     total_moves = length(nmove)
#   ) %>%
#   as.data.frame(stringsAsFactors = FALSE) -> ss_total_moves
# 
# table(ss_total_moves$total_moves) # this shows us the number of registered addresses logged per person
# 
# dist <- table(ss_total_moves$total_moves - 1) # this shows us the number of moves, as for the individuals with just one registered address, they are now logged as having zero moves
# 
# # for range and IQR
# summary(dist)
# 
# # percent of individuals that have more than 30 moves
# sum(dist[31:82])/sum(dist)*100
# 
# # plotting total moves over lifetime
# png("Figures/hist_s_total_moves.png", res = 300, height = 15, width = 20, units = "cm")
# # plot with -1 to correct nmove from listing order of registration to representing moves
# plot(table(ss_total_moves$total_moves -1),
#      main = "Total number of moves recorded over lifetime",
#      xlab = "Lifetime number of moves",
#      ylab = "Number of RPs",
#      xlim = c(0,131),
#      xaxt = "n",
#      ylim = c(0,250),
#      yaxt = "n",
#      bty = "n",
#      col = "grey41",
#      font.main = 1)
# axis(1, at = seq(0, 200, by = 10))
# axis(2, at = seq(0, 2500, by = 50))
# abline(v = median(ss_total_moves$total_moves - 1), col = plasma(20)[16], lty = 2, lwd = 2)
# dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### Simulating and plotting predicted moves |post-stratification|

person_ids <- sort(unique(d$person_id)) 
n_rp <- 100 #36595  #36595 for full set
rp_sub <- sample(person_ids, size = n_rp)
dm_sim <- subset(d, d$person_id %in% rp_sub)

# re index
dm_sim$person_id <- match(dm_sim$person_id, rp_sub)  

uniq_ages <- length(unique(dm_sim$age))

n_individuals <- rep(0, uniq_ages)
real_moves <- rep(0, uniq_ages)
est_moves_mean <- rep(0, uniq_ages)
est_moves_hpdi <- matrix(0, 2, uniq_ages)
est_moves_pi_25 <- matrix(0, 2, uniq_ages)
est_moves_pi_50 <- matrix(0, 2, uniq_ages)
est_moves_pi_75 <- matrix(0, 2, uniq_ages)
age_indices <- seq(from = 1, to = uniq_ages, by = 1)

pb <- txtProgressBar(min(age_indices), max(age_indices), style=3)
for(age_i in age_indices) {
  
  est_moves_sum <- rep(0, 2000)
  
  dm_individuals_at_age <- dm_sim[which(dm_sim$age == age_i - 1),]
  
  if (nrow(dm_individuals_at_age) == 0) next # should never happen on whole dataset
  
  # loop through all persons at given age
  for(ind in 1:nrow(dm_individuals_at_age)) {
    
    # Simply counting the number of individuals per age group
    # to check that it looks reasonable
    n_individuals[age_i] <- n_individuals[age_i] + 1
    
    # get real moves for individual
    n_moves <- dm_individuals_at_age[ind, "n_moves"]
    # add real moves for individual to real_moves[age_i]
    real_moves[age_i] <- real_moves[age_i] + n_moves
    
    # calculate N lambdas
    person_id <- dm_individuals_at_age[ind, "person_id"]
    lambdas <- exp(post_pois$mu + rnorm(2000,0,post_pois$sd_id) + post_pois$beta[, age_i])
    
    # calculate N est_moves for individual
    est_moves <- sapply(lambdas, function(i) rpois(1, i) )
    
    # Sum est_moves across individuals, keeping the sample dimension
    est_moves_sum <- est_moves_sum + est_moves # est_moves is for one person, here we sum all people in age group
    
    # add mean of est moves
    est_moves_mean[age_i] <- est_moves_mean[age_i] + mean(est_moves)
  }
  
  # calculate HPDI est_moves for individual
  hpdi <- HPDI(est_moves_sum, prob=0.89) 
  PI_25 <- PI(est_moves_sum, prob=0.25)
  PI_50 <- PI(est_moves_sum, prob=0.50)
  PI_75 <- PI(est_moves_sum, prob=0.75)
  
  # add HPDI est_moves for individual to est_moves_hpdi
  est_moves_hpdi[1, age_i] <- est_moves_hpdi[1, age_i] + hpdi[1]
  est_moves_hpdi[2, age_i] <- est_moves_hpdi[2, age_i] + hpdi[2]
  
  # add PI est moves for individual to est_moves_pi
  est_moves_pi_25[1, age_i] <- est_moves_pi_25[1, age_i] + PI_25[1]
  est_moves_pi_25[2, age_i] <- est_moves_pi_25[2, age_i] + PI_25[2]
  
  est_moves_pi_50[1, age_i] <- est_moves_pi_50[1, age_i] + PI_50[1]
  est_moves_pi_50[2, age_i] <- est_moves_pi_50[2, age_i] + PI_50[2]
  
  est_moves_pi_75[1, age_i] <- est_moves_pi_75[1, age_i] + PI_75[1]
  est_moves_pi_75[2, age_i] <- est_moves_pi_75[2, age_i] + PI_75[2]
  
  setTxtProgressBar(pb, age_i)
}
close(pb)

# plotting colors
vir <- plasma(20)
vir_int <- plasma(20, alpha = 0.8)

png("Figures/moves_ages_PI.png", res = 300, height = 15, width = 20, units = "cm")
plot(y = real_moves,
     x = age_indices-1,
     ylim = c(0,max(real_moves)), 
     xlim = c(0,max(age_indices-1)),
     main = "Total moves per age", 
     xlab = "Age", 
     ylab = "Total number of moves", 
     bty = "n", 
     col = "grey41",
     pch = 19,
     xaxt = "n", 
     yaxt = "n", 
     font.main = 1)
axis(1, at = seq(0, 95, by = 10))
axis(2, at = seq(0, 10000, by = 1000))

shade(est_moves_pi_75, age_indices - 1, col = vir_int[16])
shade(est_moves_pi_50, age_indices - 1, col = vir_int[15])
shade(est_moves_pi_25, age_indices - 1, col = vir_int[13])

dev.off()


#-----------------------------------------------------------------------------------------------------------------------
### plotting beta estimate, implied number of moves per year from model estimates

# reformatting dm_sim to obtain averages from data
dm_sim %>%
  group_by(age) %>%
  summarize(
    total_moves = sum(n_moves),
    average_moves = mean(n_moves)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> mean_moves

# processing posterior
exp_mu <- sapply( 1:uniq_ages , function(i) exp( post_pois$mu + post_pois$beta[,i] + rnorm(2000,0,post_pois$sd_id) ) )
beta_mod_int <- apply(exp_mu, 2, HPDI)
beta_mod <- apply(exp_mu, 2, mean)


# plotting
png("Figures/beta_estimates.png", res = 300, height = 15, width = 20, units = "cm")
plot(y = 0:(length(beta_mod)-1),
     x = 0:(length(beta_mod)-1), 
     ylim = c(0, max(beta_mod_int)),  
     xlim = c(0,(length(beta_mod)-1)),
     xlab = "Age", 
     ylab = "Estimated number of moves per year", 
     main = "Age-based estimated moves per year", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     type = "n")
axis(1, at = seq(0,90, by = 10))
axis(2, at = seq(0, max(beta_mod_int), by = 0.1))

shade( apply(exp_mu, 2, PI, 0.75 ) , 0:(length(beta_mod)-1), col = vir_int[16])
shade( apply(exp_mu, 2, PI, 0.5 ) , 0:(length(beta_mod)-1), col = vir_int[15])
shade( apply(exp_mu, 2, PI, 0.25 ) , 0:(length(beta_mod)-1), col = vir_int[13])
points(y = mean_moves$average_moves, x = 0:(length(mean_moves$average_moves)-1), col = "grey 80", pch = 19)
lines(x = 0:(length(beta_mod)-1), y = beta_mod, col = vir[5], lty = 2, lwd = 2 )
dev.off()
#-----------------------------------------------------------------------------------------------------------------------
### plotting individual trajectories - how RPs acquire moves over their observation periods

# create person tables and cumulative sum tables
dm_sim %>%
  group_by(person_id) %>%
  summarize(
    moves = sum(n_moves),
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> dm_sim_person # person table with total number of lifetime moves

# cumulative counts of moves
dm_sim %>%
  group_by(person_id) %>%
  mutate(
    moves_cumsum = cumsum(n_moves),
  ) %>%
  ungroup() %>%
  as.data.frame(stringsAsFactors = FALSE) -> dm_sim_cumsum

plot_cols <- plasma(max(dm_sim_cumsum$moves_cumsum) + 1, alpha = 0.5)

for(i in 1:nrow(dm_sim_person)){
  moves <- dm_sim_person$moves[i]
  dm_sim_person$col_real[i] <- plot_cols[moves + 1]
}

# removing hard-to-see-yellow from plot
#dm_sim_person[dm_sim_person$col_real == "#F0F92180", "col_real"] <- "#F0F62580"

# Plot move accumulation
png("Figures/trajectories_nmoves.png", res = 300, height = 14, width = 14, units = "cm")

# Plotting for real trajectories
plot(dm_sim_cumsum$moves_cumsum ~ dm_sim$age ,
     type = "n", 
     xlab  = "Age", 
     ylab = "Number of moves", 
     main = "Cumulative totals of number of moves per RP",
     font.main = 1,
     bty = "n",
     ylim = c(0,200),
     xlim = c(0,100))

for (i in 1:length(unique(dm_sim$person_id))){
  rp_id <- unique(dm_sim$person_id)[i]
  temp <- dm_sim[which(dm_sim$person_id == rp_id), ]
  temp_m <- dm_sim_cumsum[which(dm_sim_cumsum$person_id == rp_id), ]
  lines(temp$age, temp_m$moves_cumsum, lwd = 1, col = dm_sim_person[which(dm_sim_person$person_id == rp_id),"col_real"])
}

#legend
points(x = rep(100, times = 99), y = seq(from = 100, to = 150, length.out = 99), col = plasma(99) )
text(x = 80, y = 100, "low lifetime moves", adj = 0, cex = 0.5)     # update y coordinate
text(x = 80, y = 150, "high lifetime moves", adj = 0, cex = 0.5)    # update y coordinate    
dev.off()


#-----------------------------------------------------------------------------------------------------------------------
### plotting cohort model

#cohorts <- seq(from = 1850, to = 1922 , by = 1) # note: rewrites cohorts from 2_fit_modes.R

years_length <- length(cohorts) + ncol(m_coh_samples[[length(cohorts)]][["beta"]])
coh_label <- seq(from = 1850, to = 1850 + years_length - 1 , by = 1)

# constructing matrix for plotting
coh_mat <- matrix(data = NA, nrow = 100, ncol = years_length+10)

for(x in 1:length(cohorts)){
  
  exp_mu_coh <- sapply( 1:ncol(m_coh_samples[[x]][["beta"]]) , function(i) exp( m_coh_samples[[x]][["mu"]] + m_coh_samples[[x]][["beta"]][,i] + rnorm(2000,0,m_coh_samples[[x]][["sd_id"]]) ) )
  
  m <- apply(exp_mu_coh, 2, mean)
  
  for(y in 1:length(m)){
    if (x+(y-1) < years_length) {
      # matrix construction for years of time on x axis
      coh_mat[y, x+(y-1)] <-  m[y]
    }
  }
}

coh_mat <- t(coh_mat)

grey <- rgb(0,0,0,alpha=0.2)

# note that observation ends at 1945 in HSN sample

# cohort heatmap
png("Figures/heatmap.png", res = 300, height = 20, width = 20, units = "cm")
image(x = 0:(years_length - 1 + 10) , # + 10 for plotting
      y = 0:100, 
      z = coh_mat, 
      col = plasma(1000), 
      ylim = c(0,100), 
      xlab = "Birth year", 
      ylab = "Age", 
      bty = "n",
      xaxt = "n" )
axis(1, at= seq(0 ,years_length + 10, by = 10), labels= seq(1850, 1950, by = 10))
abline(h = c(20, 40, 60, 80), col = grey )
abline(v = seq(10, 100, by = 10), col = grey)
title("Mobility over age and time", font.main = 1)

#legend
points(x = rep(2, times = 99), y = seq(from = 85, to = 95, length.out = 99), col = plasma(99))
text(x = 2, y = 82, "low mobility", adj = 0, cex = 0.8)     
text(x = 2, y = 98, "high mobility", adj = 0, cex = 0.8)  
dev.off()


### Supplementary plots
#-----------------------------------------------------------------------------------------------------------------------
# individual effects

exp_mu_a <- sapply( 1:36595 , function(i) exp( post_pois$a[,i]  ) )
a_mod <- apply(exp_mu_a, 2, mean)
a_mod_int <- apply(exp_mu_a, 2, PI, 0.5)

results <- matrix(0, 3, length(unique(dm_sim$person_id)))
results[1,] <- a_mod
results[2,] <- a_mod_int[1,]
results[3,] <- a_mod_int[2,]

# order by mean
results <- results[,order(results[1,])]

png("Figures/alphas.png", res = 300, height = 15, width = 20, units = "cm")
plot(y = 0:(length(a_mod)-1),
     x = 0:(length(a_mod)-1), 
     ylim = c(0, max(a_mod_int)),  
     xlim = c(0,(length(a_mod)-1)),
     xlab = "RP index", 
     ylab = "Alpha estimate", 
     main = "Alpha estimate for each RP", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     yaxt = "n",
     type = "n")
axis(1, at = seq(0,40000, by = 10000))
axis(2, at = seq(0, 25, by = 5))

shade(results[c(2,3),] , 0:(length(a_mod)-1), col = vir_int[15])
lines(results[1,], col = vir[5])
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### sample age structure 

# png("Figures/age_representation.png", res = 300, height = 15, width = 20, units = "cm")
# plot(table(d$age),
#      main = "Representation of ages", 
#      xlab = "Age", 
#      ylab = "Number of observations", 
#      xlim = c(0,100),
#      xaxt = "n",
#      ylim = c(0,40000), 
#      bty = "n", 
#      col = "grey41",
#      font.main = 1)
# axis(1, at = seq(0, 100, by = 10))
# dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### sample cohort representation 

# d %>%
#   group_by(person_id) %>%
#   summarise(b_y = first(b_y)) %>%
#   as.data.frame(stringsAsFactors = FALSE) -> cohorts
# 
# png("Figures/birth_year_rep.png", res = 300, height = 15, width = 20, units = "cm")
# plot(table(cohorts$b_y), 
#      main = "Birth year representation in the HSN",
#      font.main = 1,
#      ylim = c(0,max(table(cohorts$b_y))),
#      xlim = c(1830,1940),
#      xaxt = "n",
#      xlab = "Year",
#      ylab = "Number of RPs",
#      col = "grey41",
#      bty = "n")
# axis(1, at = seq(1840, 1930, 10), labels = seq(1840, 1930, 10))
# dev.off()

#------------------------------------------------------------------------------
# pois model estimates

png("Figures/model_estimates.png", res = 300, height = 6, width = 8, units = "cm")
plot(precis(stanfit_pois, pars = c("mu", "rho", "eta"), depth = 2), ylab = "Parameter")
dev.off()

#------------------------------------------------------------------------------
# gamma pois model estimates

png("Figures/negbin_estimates.png", res = 300, height = 6, width = 8, units = "cm")
plot(precis(stanfit_negbin, pars = c("mu", "rho", "eta", "phi"), depth = 2), ylab = "Parameter") 
dev.off()

#------------------------------------------------------------------------------
# gamma pois betas

exp_mu_negbin <- sapply( 1:91 , function(i) exp( post_negbin$mu + post_negbin$beta[,i] + rnorm(2000,0,post_negbin$sd_id) ) )
beta_mod_int_negbin <- apply(exp_mu_negbin, 2, HPDI)
beta_mod_negbin <- apply(exp_mu_negbin, 2, mean)

png("Figures/beta_estimates_neg.png", res = 300, height = 15, width = 20, units = "cm")
plot(y = 0:(length(beta_mod_negbin)-1),
     x = 0:(length(beta_mod_negbin)-1), 
     ylim = c(0, max(beta_mod_int_negbin)),  
     xlim = c(0,(length(beta_mod_negbin)-1)),
     xlab = "Age", 
     ylab = "Estimated number of moves per year", 
     main = "Age-based estimated moves per year", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     type = "n")
axis(1, at = seq(0,90, by = 10))
axis(2, at = seq(0, max(beta_mod_int_negbin), by = 0.1))

shade( apply(exp_mu_negbin, 2, PI, 0.75 ) , 0:(length(beta_mod_negbin)-1), col = vir_int[16])
shade( apply(exp_mu_negbin, 2, PI, 0.5 ) , 0:(length(beta_mod_negbin)-1), col = vir_int[15])
shade( apply(exp_mu_negbin, 2, PI, 0.25 ) , 0:(length(beta_mod_negbin)-1), col = vir_int[13])
points(y = mean_moves$average_moves, x = 0:(length(mean_moves$average_moves)-1), col = "grey 80", pch = 19)
lines(x = 0:(length(beta_mod_negbin)-1), y = beta_mod_negbin, col = vir[5], lty = 2, lwd = 2 )

dev.off()

#-----------------------------------------------------------------------------------------------------------------------
# Rhat values against effective samples

png("Figures/rhat_neff.png", res = 300, height = 15, width = 20, units = "cm")
dashboard(stanfit_pois)
dev.off()

