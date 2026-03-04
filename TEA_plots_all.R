remove(list = ls())
setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")

source("TEA_functions.R")


M <- 100
B <- 1000
lambda <- 9

params <- rep(0.7, 20)
#params <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
#params <- seq(0.6, 0.9, length = 20)   # true means per version
#params <- c(rep(0.6, 19), 0.9)

a_theta <- 1
b_theta <- 1
mean_theta <- 0.5                          # prior mean
sd_theta   <- 0.5                          # prior sd
sigma      <- 0.5                          # measurement noise
p_eps <- 0.5

###1) Trajectories

# Suppose you ran posterior_sim_gaussian()
res <- posterior_sim_gaussian(params, M, B, lambda, mean_theta, sd_theta, sigma, p_eps)
plot_trajectories(res$thetas, res$epsilons, params)

res <- posterior_sim_binomial(params, M, B, lambda, a_theta, b_theta, p_eps)
plot_trajectories(res$thetas, res$epsilons, params)


###2) Power curves

lambdas <- 1:50 * 4
p1_binomial <- c( 0.453, 0.477, 0.492, 0.557, 0.548, 0.576, 0.571, 0.585, 0.562, 0.562, 0.581, 0.564, 0.573, 0.611, 0.578, 0.614, 0.611, 0.614, 0.632, 0.621, 0.623, 0.622, 0.617, 0.616, 0.618, 0.609, 0.632, 0.664, 0.629, 0.619, 0.627, 0.658, 0.637, 0.66, 0.666, 0.655, 0.662, 0.683, 0.67, 0.677, 0.701, 0.69, 0.647, 0.68, 0.689, 0.694, 0.677, 0.688, 0.67, 0.688 )
p2_binomial <- c( 0.3, 0.355, 0.365, 0.44, 0.421, 0.438, 0.478, 0.514, 0.535, 0.557, 0.575, 0.584, 0.573, 0.628, 0.616, 0.615, 0.638, 0.608, 0.63, 0.659, 0.633, 0.649, 0.659, 0.648, 0.685, 0.675, 0.7, 0.682, 0.667, 0.689, 0.69, 0.695, 0.719, 0.722, 0.695, 0.705, 0.726, 0.715, 0.735, 0.712, 0.723, 0.744, 0.716, 0.731, 0.754, 0.741, 0.746, 0.765, 0.738, 0.74 )
p1_gaussian <- c( 0.447, 0.426, 0.442, 0.429, 0.407, 0.457, 0.412, 0.41, 0.42, 0.416, 0.384, 0.398, 0.396, 0.399, 0.371, 0.412, 0.4, 0.4, 0.365, 0.399, 0.4, 0.377, 0.391, 0.365, 0.379, 0.375, 0.357, 0.356, 0.34, 0.346, 0.341, 0.345, 0.346, 0.337, 0.335, 0.337, 0.349, 0.336, 0.338, 0.335, 0.358, 0.321, 0.33, 0.327, 0.342, 0.329, 0.359, 0.324, 0.33, 0.328 )
p2_gaussian <- c( 0.34, 0.339, 0.303, 0.357, 0.347, 0.385, 0.417, 0.376, 0.41, 0.393, 0.42, 0.422, 0.434, 0.452, 0.438, 0.443, 0.435, 0.448, 0.421, 0.477, 0.482, 0.469, 0.445, 0.469, 0.473, 0.444, 0.447, 0.433, 0.454, 0.48, 0.446, 0.43, 0.483, 0.467, 0.443, 0.459, 0.459, 0.442, 0.462, 0.464, 0.44, 0.464, 0.432, 0.451, 0.442, 0.445, 0.45, 0.472, 0.422, 0.477 )

plot_power_curves(lambdas, p1_binomial, p2_binomial, "gradual", "abrupt")
plot_power_curves(lambdas, p1_gaussian, p2_gaussian, "gradual", "abrupt")


###3) Variance, bias, MSE

