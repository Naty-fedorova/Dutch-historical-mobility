# plotting and processing code

# dependencies
library(rethinking)
library(viridis)
library(testthat)
library(tidyverse)

# create folder to contain plots
dir.create("Figures")

# load in data
# df used for model
# when working with real data:
d <- read.csv("s_person_year_df.csv", stringsAsFactors = FALSE)

# df of birth-death lifecourses
s <- read.csv(file = "w_s.csv", stringsAsFactors = FALSE)

# when working with simulated data
d <- read.csv("s_person_year_sim.csv", stringsAsFactors = FALSE)

# load in samples if necessary
# load("post_pois.RData")

#-----------------------------------------------------------------------------------------------------------------------
### plotting total residential moves over the lifetime (only possible with HSN data)

# check first reg address is at age 0
expect_true(nrow(s[which((s$mun_move_category == "first") & (s$age_at_move != 0)),]) == 0)   # first reg address is always at age 0

ss <- s[which((s$death_y - s$birth_y) >= 20),]

ss %>%
  group_by(person_id) %>%
  summarize(
    total_moves = length(nmove)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> ss_total_moves

table(ss_total_moves$total_moves) # this shows us the number of registered addresses logged per person

dist <- table(ss_total_moves$total_moves - 1) # this shows us the number of moves, as for the individuals with just one registered address, they are now logged as having zero moves

# for range and IQR
summary(dist)

# percent of individuals that have more than 30 moves
sum(dist[31:82])/sum(dist)*100

# plotting total moves over lifetime
png("Figures/hist_s_total_moves.png", res = 300, height = 15, width = 20, units = "cm")
# plot with -1 to correct nmove from listing order of registration to representing moves
plot(table(ss_total_moves$total_moves -1),
     main = "Total number of moves recorded over lifetime",
     xlab = "Total number of moves",
     ylab = "Number of RPs",
     xlim = c(0,131),
     xaxt = "n",
     ylim = c(0,250),
     yaxt = "n",
     bty = "n",
     col = "grey41",
     font.main = 1)
axis(1, at = seq(0, 200, by = 10))
axis(2, at = seq(0, 2500, by = 50))
abline(v = median(ss_total_moves$total_moves - 1), col = plasma(20)[16], lty = 2, lwd = 2)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### Simulating counterfactuals from model

# select a different sample from s_person_year, keep same n_rp
set.seed(2)
person_ids <- sort(unique(d$person_id)) 
rp_sub <- sample(person_ids, size = n_rp)

dm_sim <- subset(d, d$person_id %in% rp_sub)

dm_sim$person_id <- match(dm_sim$person_id, rp_sub)  
#dm_sim$person_id <- coerce_index(dm_sim$person_id)   

# init
dm_sim$lambda_sim <- NA
dm_sim$a_mod <- NA
dm_sim$y_mod <- NA

# get coefficients from model posterior 
beta_mod <- apply(post_pois$beta, 2, mean)
beta_mod_int <- apply(post_pois$beta, 2, HPDI)
a_mod <- apply(post_pois$a, 2, mean)
a_mod_int <- apply(post_pois$a, 2, HPDI)
mu_mod <- mean(post_pois$mu)

p_ids <- unique(dm_sim$person_id)

# simulating lambda from model posterior 
for(i in 1:n_rp){
  p_id <- p_ids[i]
  ind <- dm_sim$person_id == p_id
  dm_sim$a_mod[ind] <- a_mod[p_ids[i]]
  
  for(j in 1:nrow(dm_sim[ind, ])){
    log_lambda <-  mu_mod + dm_sim$a_mod[ind][j] + beta_mod[dm_sim$age[ind] [j] +1] # +1 as offset so that beta is indexed correctly
    dm_sim$lambda_sim[ind][j] <- exp(log_lambda)
  }
}

expect_true(sum(is.na(dm_sim$lambda_sim)) == 0) # if this is false it is because the sample has ages that were not in the sample the model was run on (i.e. older individuals)

# simulating var number of counterfactual ys from counterfactual lambda
dm_sim$y_mod <- as.list(dm_sim$y_mod)

var <- 100

for(i in 1: nrow(dm_sim)){
  dm_sim$y_mod[i] <- list(rpois(var, dm_sim$lambda_sim[i]))
}

# summarizing simulated counterfactuals
age_mat <- matrix(data = NA, ncol = var + 2, nrow = length(unique(dm_sim$age)))
age_mat[,1] <- min(dm_sim$age):max(dm_sim$age)

for (i in 1:nrow(age_mat)){
  age <- age_mat[i,1]
  age_subset <- dm_sim[which(dm_sim$age == age),]
  
  age_mat[i,2] <- sum(age_subset$n_moves)
  
  for(k in 1:var){
    
    mod_mat <- matrix(data = unlist(age_subset$y_mod), nrow = nrow(age_subset), ncol = var, byrow = TRUE)
    
    age_mat[i,k + 2] <- sum(mod_mat[,k])
  }
}

# remove NAs if there are any
complete.cases(age_mat)
age_mat <- na.omit(age_mat) 

mod_mean <- apply(age_mat[,3:var+2], 1, mean)
mod_HPDI <- apply(age_mat[,3:var+2], 1, HPDI)


### Plotting counterfactual predicted moves

# define colors for plotting
vir <- plasma(20)
vir_int <- plasma(20, alpha = 0.8)


png("Figures/moves_ages.png", res = 300, height = 15, width = 20, units = "cm")
plot(y = age_mat[,2],
     x = age_mat[,1],
     ylim = c(0,max(age_mat[,2])), 
     xlim = c(0,90),
     main = "Total moves per age in the sample", 
     xlab = "Age", 
     ylab = "Total number of moves", 
     bty = "n", 
     col = "grey41",
     pch = 19,
     xaxt = "n", 
     yaxt = "n", 
     font.main = 1)
shade(mod_HPDI, age_mat[,1], col = vir_int[16])
lines(age_mat[,1], mod_mean,  col = vir[16])

axis(1, at = seq(0, 95, by = 10))
axis(2, at = seq(0, 10000, by = 2000))
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### plotting beta estimate + mu, implied number of moves per year from model estimates

# reformatting dm_sim to obtain averages from data
dm_sim %>%
  group_by(age) %>%
  summarize(
    total_moves = sum(n_moves),
    average_moves = mean(n_moves)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> mean_moves

# plotting
png("Figures/beta_estimates.png", res = 300, height = 15, width = 20, units = "cm")
plot(y = 0:90,
     x = 0:(length(beta_mod)-1), 
     ylim = c(0, 0.5),
     xlim = c(0,(length(beta_mod)-1)),
     xlab = "Age", 
     ylab = "Estimated number of moves per year", 
     main = "Age-based estimated moves per year", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     type = "n")
axis(1, at = seq(0,100, by = 10))
axis(2, at = seq(0, 0.5, by = 0.1))
shade(exp(beta_mod_int + mu_mod), 0:(length(beta_mod)-1), col = vir_int[16])
points(y = mean_moves$average_moves, x = 0:90, col = "grey 41", pch = 19)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### plotting individual trajectories - how RPs acquire moves over their observation periods

# select one counterfactual run to compare
for(i in 1:nrow(dm_sim)){
  dm_sim$y_mod_sample[i] <- unlist(dm_sim$y_mod[i])[1] # take the first from the set of simulation runs
}

# create person tables and cumulative sum tables
dm_sim %>%
  group_by(person_id) %>%
  summarize(
    moves = sum(n_moves),
    y_mod = sum(y_mod_sample)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> dm_sim_person

# cumulative counts of moves
dm_sim %>%
  group_by(person_id) %>%
  mutate(
    moves_cumsum = cumsum(n_moves),
    y_mod_cumsum = cumsum(y_mod_sample)
  ) %>%
  ungroup() %>%
  as.data.frame(stringsAsFactors = FALSE) -> dm_sim_cumsum

dm_sim_person$y_mod <- as.factor(dm_sim_person$y_mod, ordered = TRUE)
dm_sim_person$moves <- as.factor(dm_sim_person$moves, ordered = TRUE)

max_col <- max(c(dm_sim_cumsum$y_mod_cumsum, dm_sim_cumsum$moves_cumsum), na.rm = TRUE) # pick highest number of moves
plot_cols <- plasma(max_col, alpha = 0.5)

color_sim <- plot_cols[dm_sim_person$y_mod + 1] # add 1 to index age 0
color_real <- plot_cols[dm_sim_person$moves + 1] # add 1 to index age 0
dm_sim_person$col_sim <- color_sim
dm_sim_person$col_real <- color_real

# Plot move accumulation
png("Figures/trajectories_nmoves.png", res = 300, height = 7.5, width = 7.5, units = "cm")
par(mfrow=c(1,1))

# Plotting for real trajectories
plot(dm_sim_cumsum$y_mod_cumsum ~ dm_sim$age ,
     type = "n", 
     xlab  = "Age", 
     ylab = "Number of moves", 
     main = "Cumulative totals of number of moves per RP",
     font.main = 1,
     bty = "n",
     ylim = c(0, 200),
     xlim = c(0,100))

for (i in 1:length(unique(dm_sim$person_id))){
  rp_id <- unique(dm_sim$person_id)[i]
  temp <- dm_sim[which(dm_sim$person_id == rp_id), ]
  temp_m <- dm_sim_cumsum[which(dm_sim_cumsum$person_id == rp_id), ]
  lines(temp$age, temp_m$moves_cumsum, lwd = 1, col = dm_sim_person[which(dm_sim_person$person_id == rp_id),"col_real"])
}

#legend
points(x = rep(0, times = 99), y = seq(from = 150, to = 190, length.out = 99), col = plasma(99) )
text(x = 2, y = 150, "low lifetime moves", adj = 0, cex = 0.5)     # update y coordinate
text(x = 2, y = 190, "high lifetime moves", adj = 0, cex = 0.5)    # update y coordinate    
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### Plotting alpha estimate + mu, implied moves per year per individual

a <- t(rbind(a_mod_int, a_mod))
a_o <- a[order(a[,"a_mod"]), ]
a_int <- t(a_o[,1:2])

# averages per RP
# reformatting dm_sim to obtain averages from data
dm_sim %>%
  group_by(person_id) %>%
  summarize(
    average_moves_rp = mean(n_moves)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> mean_moves_rp


png("Figures/a_est.png", res = 300, height = 15, width = 20, units = "cm")
plot(exp(a_o[,"a_mod"] + mu_mod),
     xlab = "RP Index", 
     ylab = "Estimated number of moves per year", 
     main = "Individual-based estimated moves per year", 
     font.main = 1,
     bty = "n", 
     ylim = c(0,3.5),
     xlim = c(0, 40000),
     type = "n")
shade(exp(a_int + mu_mod), 1:length(a_mod), col = vir_int[16])
points(sort(mean_moves_rp$average_moves_rp), col = "grey41", pch = 20)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### plotting cohort model

years_length <- length(cohorts) + ncol(m_coh_samples[[length(cohorts)]][["beta"]])
coh_label <- seq(from = 1860, to = 1860 + years_length - 1 , by = 1)

# constructing matrix for plotting
coh_mat <- matrix(data = NA, nrow = 100, ncol = years_length)

for(x in 1:length(cohorts)){
  
  m <- exp(apply(m_coh_samples[[x]][["beta"]], 2, mean) + mean(m_coh_samples[[x]][["mu"]])) 
  
  for(y in 1:length(m)){
    
    # need if statement to skip population of end values if removing extreme values
    if (x+(y-1) < years_length) {
      # matrix construction for years of time on x axis
      coh_mat[y, x+(y-1)] <-  m[y]
    }
    
    # matrix construction with cohort on x axis
    #coh_mat[y, x] <-  m[y]
  }
}

coh_mat <- t(coh_mat)

grey <- rgb(0,0,0,alpha=0.2)

# cohort heatmap
png("Figures/heatmap.png", res = 300, height = 20, width = 20, units = "cm")
image(x = 0:years_length, 
      y = 0:100, 
      z = coh_mat, 
      col = plasma(1000), 
      ylim = c(0,100), 
      xlab = "Year", 
      ylab = "Age", 
      bty = "n",
      xaxt = "n" )
axis(1, at= seq(1 ,years_length, by = 10), labels= seq(1860, 1940, by = 10))
abline(h = c(20, 40, 60, 80), col = grey )
abline(v = seq(10, 80, by = 10), col = grey)
title("Mobility over age and time", font.main = 1)

#legend
points(x = rep(80, times = 99), y = seq(from = 85, to = 95, length.out = 99), col = plasma(99) )
text(x = 71, y = 85, "low mobility", adj = 0, cex = 0.7)     
text(x = 71, y = 95, "high mobility", adj = 0, cex = 0.7)  
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### plotting pop structure for supplementary

png("Figures/age_representation.png", res = 300, height = 15, width = 20, units = "cm")
plot(table(d$age),
     main = "Representation of ages", 
     xlab = "Age", 
     ylab = "Number of observations", 
     xlim = c(0,100),
     xaxt = "n",
     ylim = c(0,40000), 
     bty = "n", 
     col = "grey41",
     font.main = 1)
axis(1, at = seq(0, 100, by = 10))
dev.off()

#-----------------------------------------------------------------------------------------------------------------------
### plotting cohort representation in sample of HSN for supplementary

d %>%
  group_by(person_id) %>%
  summarise(b_y = first(b_y)) %>%
  as.data.frame(stringsAsFactors = FALSE) -> cohorts

png("Figures/birth_year_rep.png", res = 300, height = 15, width = 20, units = "cm")

plot(table(cohorts$b_y), 
     main = "Birth year representation in the HSN",
     font.main = 1,
     ylim = c(0,max(table(cohorts$b_y))),
     xlim = c(1830,1940),
     xaxt = "n",
     xlab = "Year",
     ylab = "Number of RPs",
     col = "grey41",
     bty = "n")
axis(1, at = seq(1840, 1930, 10), labels = seq(1840, 1930, 10))
dev.off()

#------------------------------------------------------------------------------

# plot for pois model estimates
png("Figures/model_estimates.png", res = 300, height = 10, width = 8, units = "cm")
plot(precis(m_nc_real, pars = c("mu", "rho", "eta"), depth = 2), ylab = "Parameter")
dev.off()

#------------------------------------------------------------------------------
# plotting for gamma pois model estimates

png("Figures/negbin_estimates.png", res = 300, height = 10, width = 8, units = "cm")
plot(precis(stanfit_negbin, pars = c("mu", "rho", "eta", "phi"), depth = 2), ylab = "Parameter") 
dev.off()

# gamma pois betas
beta_mod_neg <- apply(post_negbin$beta, 2, mean)
beta_mod_int_neg <- apply(post_negbin$beta, 2, HPDI)
mu_mod_neg <- mean(post_negbin$mu)

png("Figures/beta_estimates_negbin.png", res = 300, height = 15, width = 20, units = "cm")
plot(y = 0:90,
     x = 0:(length(beta_mod_neg)-1), 
     ylim = c(0, 0.5),
     xlim = c(0,(length(beta_mod_neg)-1)),
     xlab = "Age", 
     ylab = "Estimated number of moves per year", 
     main = "Age-based estimated moves per year", 
     font.main = 1,
     bty = "n",
     xaxt = "n",
     type = "n")
axis(1, at = seq(0,100, by = 10))
axis(2, at = seq(0, 0.5, by = 0.1))
shade(exp(beta_mod_int_neg + mu_mod_neg), 0:(length(beta_mod_neg)-1), col = vir_int[16])
points(y = mean_moves$average_moves, x = 0:90, col = "grey 41", pch = 19)
dev.off()

# gamma pois alphas
a_mod_neg <- apply(post_negbin$a, 2, mean)
a_mod_int_neg <- apply(post_negbin$a, 2, HPDI)

a_neg <- t(rbind(a_mod_int_neg, a_mod_neg))
a_o_neg <- a_neg[order(a_neg[,"a_mod_neg"]), ]
a_int_neg <- t(a_o_neg[,1:2])

dm %>%
  group_by(person_id) %>%
  summarize(
    average_moves_rp = mean(n_moves)
  ) %>%
  as.data.frame(stringsAsFactors = FALSE) -> mean_moves_rp_neg

png("Figures/a_est_negbin.png", res = 300, height = 15, width = 20, units = "cm")
plot(exp(a_o_neg[,"a_mod_neg"] + mu_mod_neg),
     xlab = "RP Index", 
     ylab = "Estimated number of moves per year", 
     main = "Individual-based estimated moves per year", 
     font.main = 1,
     bty = "n", 
     ylim = c(0,3.5),
     xlim = c(0, 1000),
     type = "n")
shade(exp(a_int_neg + mu_mod_neg), 1:length(a_mod_neg), col = vir_int[16])
points(sort(mean_moves_rp_neg$average_moves_rp), col = "grey41", pch = 20)
dev.off()




