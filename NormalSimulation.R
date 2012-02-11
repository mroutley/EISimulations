start <- Sys.time()
library(nlme)
GenerateSimulation <- function(mean=1, sd=0.1, sites=5, duration=5, change=0.2) {
	# Create a simulated dataframe with intitial and final samples
	# Returns a dataframe
	time <- rep(c("initial", "final"), each=sites)
	initial <- NULL; final <- NULL
	initial <- rnorm(sites, mean=mean, sd=sd)
	final <- rnorm(sites, mean=mean*(1+change)^duration, sd=sd)
	data.frame(time=time, plot=rep(1:sites, 2), measure=c(initial, final))
}
TestSimulation <- function(simulated=GenerateSimulation()) {
	# Test for a significant effect of time on detection
	# Returns a 1 for significant test, 0 for non-significant
	simulation.model <- try(lme(measure ~ time, random= ~ 1|plot, data=simulated))
	# simulation.model <- lm(measure ~ time, data=simulated)
	ifelse(anova(simulation.model)$"p-value"[2] < 0.05, 1, 0)
}
RunSimulation <- function(sites=sites, duration=duration, change=change, mean=mean, sd=sd, iterations=5) {
	# Wraps GenerateSimulation() and TestSimulation() in a loop to generate a sequence of statistical tests
	# Returns the proportion of significant tests
	powerList <- NULL
	for (i in seq(iterations)) {
		simulated <- GenerateSimulation(mean=mean, sd=sd, sites=sites, duration=duration, change=change)
		powerList[i] <- TestSimulation(simulated)
	}
	sum(powerList, na.rm=TRUE)/iterations
}
# ==============================
# = Set the initial parameters =
# ==============================
durations <- seq(5, 10, 5)			# Vector of time periods (usually years)
changes <- seq(-0.1, 0.1, 0.025)	# Vector of percent changes per time period
sites <- seq(4, 6, 2)				# Number of sites visited
mean <- 2.1; sd <- 1.8			# Mean and sd of the simulated distribution
iterations <- 500
# =============================================================
# = Create a dataframe across the combinations of parameters  =
# = Call RunSimulation() for each combination and store power =
# =============================================================
powerSimulation <- expand.grid(durations=durations, changes=changes, sites=sites, mean=mean, sd=sd, power=NA)
for (row in 1:dim(powerSimulation)[1]) {
	powerSimulation[row,]$power <- RunSimulation(duration=powerSimulation[row,]$durations, change=powerSimulation[row,]$changes, sites=powerSimulation[row,]$sites, mean=powerSimulation[row,]$mean, sd=powerSimulation[row,]$sd, iterations=iterations)
	print(powerSimulation[row,]) # Just an update on progress
}
save(powerSimulation, file="QABRPAStreamflowGeorgeLakeMin.Rdata")
print(Sys.time() - start)