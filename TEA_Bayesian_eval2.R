remove(list = ls())


#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("TEA_functions.R")
#source("TEA_scenarios.R")

a_theta <- 1
b_theta <- 1
p_eps <- 0.5

M <- 1e2
B <- 1e3 #number of estimates


#1) Bayesian decision statistic for different values of K (number of versions)

lambda <- 9 #average sample size per version = 10

x <- 2:50

tea_probs <- vector()
naive_probs <- vector()
for(i in x){
  print(i)
  params <- seq(0.5, 0.8, length.out = i)
  #params <- rep(0.5, i)
  tea <- tea_diff_mean(params, lambda, M, B)
  naive <- naive_diff(params, lambda, M, B)

  tea_probs <- c(tea_probs, sum(tea > 0.2)/length(tea))
  naive_probs <- c(naive_probs, sum(naive > 0.2)/length(naive))
}

cols <- c("#4E79A7", "#59A14F", "#9C755F", "#B07AA1", "#76B7B2")



# first plot
plot(x, tea_probs,
     type = "l",
     col = cols[4],
     ylim = c(0,1),
     xlab = "K",
     ylab = expression(P(abs(hat(theta)[K] - hat(theta)[1]) > 0.2)),
     main = ""
)

# add second line
lines(x, naive_probs,
      col = cols[5],
      type = "l")

# optional legend
legend("topright",
       legend = c("TEA", "naive"),
       col = cols[4:5],
       lty = 1)



####
#2) At which version is the stopping decision made?


params_gradual <- seq(0.6, 0.9, length = 20)   # true means per version
params_abrupt <- c(rep(0.6, 19), 0.9)
params_drift <- c(seq(0.6, 0.9, length = 10), seq(0.9, 0.6, length = 10))
params_stable <- rep(0.7, 20)

params <- params_gradual


x <- 1:50
naive_stops <- vector()
tea_stops <- vector()

start_time <- Sys.time()

for(lambda in x){
  print(lambda)
  naive <- naive_stopping_decision(params, lambda, M, B)
  tea <- tea_stopping_decision(params, lambda, M, B)
  naive_stops <- c(naive_stops, naive)
  tea_stops <- c(tea_stops, tea)
}

end_time <- Sys.time()
time_taken <- end_time - start_time

print(paste0("time taken = ", time_taken))
print(paste0("M =", M))
print(paste0("B =", B))
print(paste0("lambdas <- c(", paste(x, collapse = ", "), ")"))
print(paste0("params <- c(", paste(params, collapse = ", "), ")"))
print(paste0("naive_stops <- c(", paste(naive_stops, collapse = ", "), ")"))
print(paste0("tea_stops <- c(", paste(tea_stops, collapse = ", "), ")"))

#("TEA_Bayesian_eval_plots.R" for the plots)

