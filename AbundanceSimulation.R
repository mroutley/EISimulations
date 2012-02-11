start <- Sys.time()
source("SimulationInput.R")
print(params)
library(nlme)
GenerateSimulation <- function(duration, change) {
	dist <- params$mus
	plots <- params$plots
	time <- rep(c("initial", "final"), each=plots)
	initial <- NULL; final <- NULL
	initial <- rpois(plots, lambda=dist)
	final <- rpois(plots, lambda=dist*(1+change)^duration)
	data.frame(time=time, plot=rep(1:plots, 2), count=c(initial, final))
}
TestSimulation <- function(simulated) {
	simulation.model <- lme(count ~ time, random= ~ 1|plot, data=simulated)
	ifelse(anova(simulation.model)$"p-value"[2] < 0.05, 1, 0)
}
RunSimulation <- function() {
	powerSequence <- data.frame(expand.grid(change=params$changes, duration=params$durations))
	powerSequence$power <- NA
	k <- 1 # Tracks position in sequence
	for (duration in params$durations) {
		for (change in params$changes){
			print(c(change, duration))
			powerList <- NULL
			for (i in seq(params$iterations)) {
				simulated <- GenerateSimulation(duration, change)
				powerList[i] <- TestSimulation(simulated)
			}
			powerSequence$power[k] <- sum(powerList, na.rm=TRUE)/params$iterations
			k <- k+1
		}
	}
	powerSequence # Return value
}
powerSequence <- RunSimulation()
powerSequence <- c(title = title, powerTable = powerSequence, params = params, simtime = (Sys.time() - start))
save(powerSequence, file="SimulationOutput.Rdata")