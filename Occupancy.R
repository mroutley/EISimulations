# =======================================================
# = Functions for calculating p & Ψ from MacKenzie 2002
# Matthew Routley <matt@routleynet.org>                 =
# =======================================================
# base <- "/Users/mroutley/Dropbox/WLNP/" # Change as appropriate
# Required libraries
library(reshape)
library(doBy)
library(ggplot)
# Common values
iterations <- 1000
# =================================================================
# = Calculate the values necessary for estimating occupancy rates =
# =================================================================
site_occupancy_values <- function(data) {
	ndot <- sum(by(data$present, list(data$site), max), na.rm = T)
	nt <- as.vector(by(data$present, list(data$visit), sum))
	Ti <- max(data$visit)
	N <- length(levels(as.factor(data$site)))
	list(ndot=ndot, nt=nt, Ti=Ti, N=N)
}
# ====================================================================
# = Calculate the maximum-likelihood estimates of the occupancy rate =
# ====================================================================
site_occupancy <- function(data) { # Returns the parameters for the presence-absence model
	values <- site_occupancy_values(data)
	nt <- values$nt
	ndot <- values$ndot
	N <- values$N
	Ti <- values$Ti
	param_range <- seq(from=0.1, to=0.9, by=0.05)
	param_space <- expand.grid(psi=param_range, pt=param_range)
	param_space$ll <- NA
	site_occupancy_likelihood <- function(x) { # Returns the maximum-likelihood estimators and likelihood
	# FIXME: The commented code is the more efficient approach, but it doesn't work.
	# The optim function should quickly find the correct response, but it insists on blowing up.
	# This is likely due to the chosen tolerance for determining sufficient differences betwen values.
	# Instead, I'm using param_space to explicitly calculate likelihood values for the entire range of psi & pt.
		psi <- x[["psi"]]
		pt <- x[["pt"]]
		# ifelse (psi > 1 | psi < 0 | pt > 1 | pt < 0, NA, log(psi^ndot * prod((pt^nt)*(1-pt)^(ndot-nt)) * (psi * (1-pt)^Ti + (1-psi))^(N-ndot)))
		psi^ndot * prod((pt^nt)*(1-pt)^(ndot-nt)) * (psi * (1-pt)^Ti + (1-psi))^(N-ndot)
	}
	for (i in 1:nrow(param_space)) {
		param_space$ll[i] <- site_occupancy_likelihood(list(psi=param_space$psi[i], pt=param_space$pt[i]))
	}
	# max_values <- optim(list(psi=0.5, pt=0.5), site_occupancy_likelihood, method = "Nelder-Mead")
	max_values <- param_space[param_space$ll == max(param_space$ll), ]
	# list(par=max_values$par, ll=-2*log(max_values$value), values=values)
	list(par=list(psi=max_values$psi, pt=max_values$pt), ll=max_values$ll, values=values)
}
# ===============================================
# = Bootstrap resampling for error calculations =
# ===============================================
sampled_data <- function(data) { # Returns the data frame after sampling with replacement
	wide <- reshape(data, v.names="present", idvar="site", timevar="visit", direction="wide")
	sampled <- wide[sample(nrow(wide), replace=TRUE), ]
	row.names(sampled) <- NULL
	sampled$site <- 1:nrow(wide)
	reshape(sampled)
}
error_distribution <- function(data) { # Returns a data frame of p & Ψ bootstrap values
	params <- NULL
	for (i in 1:iterations) {
		params <- rbind(params, site_occupancy(sampled_data(data))$par)
	}
	params <- unlist(params)
	params <- data.frame(psi=params[1:iterations], pt=params[(iterations+1):(2*iterations)])
	params
}
