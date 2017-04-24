# QCBS MCMC Workshop
install.packages("coda")


# To explore the Metropolis-Hastings, we are going to test a custom algorithm in R based on a simple linear model. First, we are going to generate the data, build a likelihood function, set our priors, and construct the posterior distribution. As a reminder, we provide Bayes rule below. 

## Bayes rule: 

# $$ P(\theta \mid Data) = \frac{P(Data \mid \theta) \, P(\theta)}{P(Data)} $$
# or it can be rewritten as
# $$ P(\theta \mid Data) \propto P(Data \mid \theta) \, P(\theta) $$


# Specifying the data -----------------------------------------------------
# In this exercise, we are going to set the data according to a linear model: y = ax+b + epsilon. 

true_a <- 5       # slope
true_b <- 0       # Intercept
true_Sd <- 10     # error 
sample_Size <- 31 # sample size 

# create independent x-values (between -15 to 15)
x <- (-(sample_Size-1)/2):((sample_Size-1)/2)

# create dependent values according to ax + b + N(0,sd)
y <-  true_a * x + true_b + rnorm(n = sample_Size,
mean = 0,
sd = true_Sd)

# Look at the data 
par(mfrow=c(1,1))
plot(x, y, 
col = "black", 
bg  = "black", 
pch = 21,
main= "Linear relationship between two variables")
abline(h = 0, v = 0, lty = 3)



# To define our Metropolis-Hastings MCMC, we need to define a statistical model. First, let's take a look at the likelihood.


# Likelihood --------------------------------------------------------------
## Likelihood or L(Data|Parameters)
# With the likelihood, we want to construct a function to answer the question: what is the function of parameters, in a statistical model, given the data. 
# $$P(Data \mid \theta_i)$$
  
likelihood <- function(param){
  # Setting parameters
  a = param[1]
  b = param[2]
  sd = param[3]
  
  # Create the predicted y 
  y_pred = a*x + b
  
  # We look how likely it is to have y, if we assume 
  # it comes from a normal distribution with
  # mean y_pred and a standard deviation of sd
  # dnorm outputs a probability (it's going to be the same length as the y vector)
  singlelikelihoods = dnorm(y, # y is the data we created that is dependent on x              
                            mean = y_pred, 
                            sd = sd, 
                            log = TRUE) # This is giving the log of the
  # probability. It is, partly, to
  # accommodate underflow and stuff
  # like that (log(10^-34)=-78.29)
  sumll = sum(singlelikelihoods) # (The logarithm of a product equals
  # the sum of the logarithms): explained
  # here https://onlinecourses.science.psu.edu/stat504/node/27
  return(sumll)   
}

# Example: plot the likelihood profile of the slope a 
# (we are going to test different slope values to see how it's working)
slopevalues <- function(test_a){
  return(likelihood(c(test_a, 
                      true_b, 
                      true_Sd)))
}

# Here we are going to create a bunch of slope values and test them with our likelihood 
sequence = seq(0, 7, by=.05)
slopelikelihoods <- lapply(sequence, 
                           slopevalues)
plot (sequence, 
      slopelikelihoods, 
      type = "l", 
      xlab = "Values of slope parameter a tested", 
      ylab = "Log likelihood of the slope tested", 
      main = "Log likelihood function of slope");abline(v=true_a, col = "red", lty = 2)



# Prior -------------------------------------------------------------------
## Prior distribution or Prior(Parameters)
# $$P(\theta_i)$$
  
# There is one prior distribution per parameter 
prior <- function(param){
  a = param[1]
  b = param[2]
  sd = param[3]
  
  # The choice of the priors might not be the best. 
  # They were chosen for simplicity in this exercice.
  # As a reminder: 
  # A Uniform distribution is defined by a minimum and a maximum value
  # A Normal distribution is defined by a mean and a standard deviation
  aprior = dunif(a, min=0, max=10, log = TRUE)
  bprior = dnorm(b, mean = 0, sd = 5, log = TRUE)
  sdprior = dunif(sd, min=0, max=30, log = TRUE)
  
  # Addition of the prior (since everything in on a 
  # log scale with "log = TRUE")
  return(aprior + bprior + sdprior)
}

# If you want to look at the uniform density distribution: 
plot(density(runif(1000000,min = 1,max = 10)), 
     main = "Uniform distribution")
# If you want to look at the normal density distribution: 
plot(density(rnorm(1000000, mean = 0, sd = 1)), 
     main = "Normal distribution")



# The posterior or P(Parameters|Data) -------------------------------------

# $$ P(\theta_i \mid Data)$$
# The posterior is a density distribution that reflects the probability of observing the data given some parameter θ~i~. It is proportional to $P(\theta_i \mid Data) P(Data)$ or the likelihood times the prior.
# In this case, since we work with the log-scale (log likelihood and prior), we can sum the two instead of multiplying them. This function will be implemented in the MCMC sampler. 

posterior <- function(param){
  return (likelihood(param) + # Likelihood
            prior(param))     # Prior
}





# Metropolis-Hastings algorithm  ------------------------------------------
# Now let's turn to the Markov chain Monte Carlo (MCMC). That's the process that will help us to determine the posterior distribution with the likelihood and the prior we defined earlier.
# Developed in part with Nicholas Constantine Metropolis (born in US, 1915-1999) and Wilfred Keith Hastings (born in Canada, 1930-2016). 
# The objective is to navigate through the parameter space (mixing), pick some values and test them. 
# Here's a resource if you want to get a sense of how a MCMC algorithm is working:  Hartig, F., J. M. Calabrese, B. Reineking, T. Wiegand, and A. Huth. 2011. Statistical inference for stochastic simulation models - theory and application. Ecology Letters 14:816–827.

proposalfunction <- function(param){
return(rnorm(3, # This outputs 3 values 
mean = param, 
sd = c(0.1,0.5,0.3))) # If you're not sure what 
# to put here, start with 
# 1/20 and see how it goes! 
}

run_metropolis_MCMC <- function(startvalue,  # The start value may bias 
                                # the first step in the MCMC, 
                                # but it is discarded later 
                                iterations){ # The number of times it's 
  # going to test values 
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    
    # Getting new proposed values 
    proposal = proposalfunction(chain[i,])
    
    # Probability of acceptance 
    # Here we compare the ratio of the posterior probability with 
    # the proposal and the posterior probability of the last value 
    probab = exp(posterior(proposal) - posterior(chain[i,])) # p1/p2 = exp[log(p1)-log(p2)]
    if (runif(1) < probab){ # If the probability we get is higher 
      # than a random probability (it could 
      # also be coded as a rbinom())
      chain[i+1,] = proposal
    } else {
      chain[i+1,] = chain[i,]
    }
  } # Closing for loop 
  return(chain)
} # End of run_metropolis_MCMC

startvalue = c(4,0,10)
chain = run_metropolis_MCMC(startvalue = startvalue, 
                            iterations = 10000)

# The burnIn will discard the first iterations that 
# might not be meaningful in the MCMC chain (it also discards the initial values)
# Another way of looking at this is that we need to 
# "warm up" the MCMC chain before we interpret any conclusion from it. 
burnIn = 5000

# Below, we look at the acceptance rate. 
# How often the values proposed from the "proposalfunction" 
# were accepted by the algorithm. For this, we look at duplicated
# lines. If a line was copied entirely, it means that the proposal 
# was rejected.
(acceptance.rate = 1-mean(duplicated(chain[-(1:burnIn),])))

### End of MCMC 




# Summary and diagnostic --------------------------------------------------

# In this part, we are going to look at the posterior distribution with an histogram and use traceplots to diagnose how well our chains were exploring the parameter space (mixing). 

par(mfrow = c(2,3))

# Histograms for posterior 
colnames(chain) = c("a","b","sd")
chain.burn = chain[-(1:burnIn),]

hist.mcmc = function(matrix,
                     xlab = "True value = red line",
                     main = "Posterior", 
                     nclass=30,
                     true.value = NA,
                     col = "red") {
  hist(matrix,nclass=nclass,  
       main=main, xlab=xlab)
  abline(v = mean(matrix), col="blue" )
  abline(v = true.value, col=col )
}

hist.mcmc(chain.burn[,"a"], main="Posterior of a")
hist.mcmc(chain.burn[,"b"], main="Posterior of b")
hist.mcmc(chain.burn[,"sd"], main="Posterior of sd")

# Trace plots function for chain mixing
traceplot.mcmc = function(matrix,
                          xlab = "True value = red line",
                          main = "Chain values", 
                          true.value = NA,
                          col = "red") {
  plot(matrix, type = "l", 
       xlab =xlab, 
       main = main )
  abline(h = true.value, col = col)
}

traceplot.mcmc(matrix = chain.burn[,"a"],
               main = "Chain values of a",
               true.value = true_a)
traceplot.mcmc(matrix = chain.burn[,"b"],
               main = "Chain values of b",
               true.value = true_b)
traceplot.mcmc(matrix = chain.burn[,"sd"],
               main = "Chain values of sd",
               true.value = true_Sd)

# for comparison, see the linear model:
summary(lm(y~x))

par(mfrow=c(1,1))

# Show the data again
plot(y~x,col = "black", bg = "black", pch = 21); abline(h=0,v=0,lty = 3)

# Line from the linear model 
abline(lm(y~x), col = "red",lty = 2)

# Extract the parameters from our Bayesian model by computing the mean of each parameter
bayes.linear.parameters = apply(chain[-(1:burnIn),],2,mean)
abline(a = bayes.linear.parameters[2], # a in this case is the intercept
       b = bayes.linear.parameters[1], # b is the slope
       col = "blue", lty = 3, lwd = 2)




# Coda: Convergence Diagnosis and Output Analysis -------------------------

# There are also packages that are designed to diagnose Bayesian models. Here we are going to explore traceplots of multiple chains, and density plots of the posterior.

# Diagnostic of the parameters by converting our matrix in a MCMC object
library(coda)
class(chain)

# Transform to an MCMC object 
chain.mcmc = mcmc(chain)
class(chain.mcmc)
str(chain.mcmc)

summary(chain.mcmc)

# Draw the diagnostic (traceplot and posterior distribution)
plot(chain.mcmc)


# You can then run multiple independent chains 
# This is particularly helpful to monitor convergence of the MCMC chains. 
# That means that independent chains should be roughly looking the same. 
chain1 = run_metropolis_MCMC(startvalue, 10000)
chain2 = run_metropolis_MCMC(startvalue, 10000)
chain3 = run_metropolis_MCMC(startvalue, 10000)

# Removing the first iterations (with the burn-in)
chain1.b = chain1[-(1:burnIn),]
chain2.b = chain2[-(1:burnIn),]
chain3.b = chain3[-(1:burnIn),]

# Converting to mcmc object 
chain1.mcmc = mcmc(chain1.b)
chain2.mcmc = mcmc(chain2.b)
chain3.mcmc = mcmc(chain3.b)

# Combining to a mcmc.list of multiple chains 
combinedchains = mcmc.list(chain1.mcmc, 
                           chain2.mcmc, 
                           chain3.mcmc)

# Diagnostic plots 
plot(combinedchains)


