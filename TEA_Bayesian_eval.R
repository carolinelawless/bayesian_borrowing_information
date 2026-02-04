remove(list = ls())


#setwd("/Users/clawless/Documents/MediTwin/bayesian_borrowing_information")
setwd("/home/clawless/simulations/bayesian_borrowing_information")
source("TEA_functions.R")
#source("TEA_scenarios.R")



M <- 1e2

a_theta <- 1
b_theta <- 1
p_eps <- 0.5

B <- 1e3 #number of estimates

K <- 20
params <- seq(0.5, 0.8, length.out = K)


#params <- rep(0.6, K)

# 
# lambda <- 4
# lambda <- 9
# 
# # tea <- tea_diff_mean(params, lambda, M, B)
# # naive <- naive_diff(params, lambda, M, B)
# # 
# # sum(tea > 0.2)/length(tea)
# # sum(naive > 0.2)/length(naive)
# 
# x <- 2:20
# x <- 1:4*5
# tea_probs <- vector()
# naive_probs <- vector()
# for(i in x){
#   print(i)
#   params <- seq(0.5, 0.8, length.out = i)
#   #params <- rep(0.5, i)
#   tea <- tea_diff_mean(params, lambda, M, B)
#   naive <- naive_diff(params, lambda, M, B)
#   
#   tea_probs <- c(tea_probs, sum(tea > 0.2)/length(tea))
#   naive_probs <- c(naive_probs, sum(naive > 0.2)/length(naive))
# }
# 
# cols <- c("#4E79A7", "#59A14F", "#9C755F", "#B07AA1", "#76B7B2")
# 
# ylim_range <- range(c(tea_probs, naive_probs))
# 
# # first plot
# plot(x, tea_probs,
#      type = "l",
#      col = cols[4],
#      ylim = c(0,1),
#      xlab = "K",
#      ylab = expression(P(abs(hat(theta)[K] - hat(theta)[1]) > 0.2)),
#      main = ""
# )
# 
# # add second line
# lines(x, naive_probs,
#       col = cols[5],
#       type = "l")
# 
# # optional legend
# legend("topright",
#        legend = c("TEA", "naive"),
#        col = cols[4:5],
#        lty = 1)
# 
# 
# 
# ####


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




