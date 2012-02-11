start <- Sys.time()
GenerateSimulation <- function(p=0.5, sites=5, frequency=2, threshold=0.5) {
	# Create a simulated dataframe of the number of times a event was detected
	# Returns a dataframe
	initial <- rbinom(sites, frequency, p)
	final <- rbinom(sites, frequency, p*(1-threshold))
	detected <- c(initial, final)
	data.frame(detected = detected, not_detected = frequency-detected, time = ordered(rep(c("initial", "final"), each=sites), levels = c("initial", "final")))
}
TestSimulation <- function(simulated=GenerateSimulation()) {
	# Test for a significant effect of time on detection
	# Returns a 1 for significant test, 0 for non-significant
	model <- glm(cbind(detected, not_detected) ~ time, family=binomial, data = simulated)
	null <- glm(cbind(detected, not_detected) ~ 1, family=binomial, data = simulated)
	ifelse(anova(model, null, test="Chisq")$"P(>|Chi|)"[2] < 0.05, 1, 0)
}
RunSimulation <- function(p, sites, frequency, threshold, iterations=5) {
	# Wraps GenerateSimulation() and TestSimulation() in a loop to generate a sequence of statistical tests
	# Returns the proportion of significant tests
	powerList <- NULL
	for (i in seq(iterations)) {
		simulated <- GenerateSimulation(p, sites, frequency, threshold)
		powerList[i] <- TestSimulation(simulated)
	}
	sum(powerList, na.rm=TRUE)/iterations
}
# ==============================
# = Set the initial parameters =
# ==============================
p0 <- seq(0.4, 0.9, 0.1)	# Initial detection probabilities
frequency <- seq(3, 9, 3)	# Number of visits to each site annualy
threshold <- 0.3			# Detection threshold
sites <- 12					# Number of sites visited
iterations <- 10000
# =============================================================
# = Create a dataframe across the combinations of parameters  =
# = Call RunSimulation() for each combination and store power =
# =============================================================
powerSimulation <- expand.grid(p=p0, frequency=frequency, sites=sites, threshold=threshold, power=NA)
for (row in 1:dim(powerSimulation)[1]) {
	powerSimulation[row,]$power <- RunSimulation(p=powerSimulation[row,]$p, sites=powerSimulation[row,]$sites, frequency=powerSimulation[row,]$frequency, threshold=powerSimulation[row,]$threshold, iterations=iterations)
	print(powerSimulation[row,]) # Just an update on progress
}
save(powerSimulation, file="SimulationOutput.Rdata")
print(Sys.time() - start)