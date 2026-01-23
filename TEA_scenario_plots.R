remove(list = ls())


params_listA <- list( #sc1 - sc8
  c(rep(0.1, 4), 0.2),
  c(rep(0.1, 4), 0.3),
  c(rep(0.1, 4), 0.4),
  c(rep(0.1, 4), 0.5),
  c(rep(0.1, 4), 0.6),
  c(rep(0.1, 4), 0.7),
  c(rep(0.1, 4), 0.8),
  c(rep(0.1, 4), 0.9)
)

params_listB <- list( #sc9 - sc16
  c(0.1, 0.1 + 1:4/(4*10/1)),
  c(0.1, 0.1 + 1:4/(4*10/2)),
  c(0.1, 0.1 + 1:4/(4*10/3)),
  c(0.1, 0.1 + 1:4/(4*10/4)),
  c(0.1, 0.1 + 1:4/(4*10/5)),
  c(0.1, 0.1 + 1:4/(4*10/6)),
  c(0.1, 0.1 + 1:4/(4*10/7)),
  c(0.1, 0.1 + 1:4/(4*10/8))
)

params_listC <- list( #sc17 - sc20
  c(0.1, 0.1, 0.9, 0.1, 0.1),
  c(0.1, 0.5, 0.9, 0.5, 0.1),
  c(0.1, 0.1, 0.1, 0.9, 0.1),
  c(0.1, 0.4, 0.7, 0.9, 0.1)
)


params_listD <- list( #sc21 - sc29
  rep(0.1, 5),
  rep(0.2, 5), 
  rep(0.3, 5),
  rep(0.4, 5),
  rep(0.5, 5),
  rep(0.6, 5),
  rep(0.7, 5),
  rep(0.8, 5),
  rep(0.9, 5)
)





################
x <- 1:5
cols <- c("#4E79A7", "#59A14F", "#9C755F", "#B07AA1", "#76B7B2")
params_list <- params_listA
n_scenarios <- length(params_list)

if (n_scenarios > length(cols)) {
  cols <- colorRampPalette(cols)(n_scenarios)
}

plot(x, params_list[[1]], 
     type = "p", pch = 16, cex = 1.3,
     col = cols[1],
     ylim = c(0, 1),
     xlab = "version",
     ylab = "true θ")

for(i in 2:length(params_list)){
  points(x, params_list[[i]], col = cols[i], pch = 16, cex = 1.3)
}

legend("topleft",
       legend = c("θ = (0.1, 0.1, 0.1, 0.1, 0.2)", "θ = (0.1, 0.1, 0.1, 0.1, 0.3)", "θ = (0.1, 0.1, 0.1, 0.1, 0.4)", "θ = (0.1, 0.1, 0.1, 0.1, 0.5)", "θ = (0.1, 0.1, 0.1, 0.1, 0.6)", "θ = (0.1, 0.1, 0.1, 0.1, 0.7)", "θ = (0.1, 0.1, 0.1, 0.1, 0.8)", "θ = (0.1, 0.1, 0.1, 0.1, 0.9)"),
       col = cols, pch = 16, cex = 0.8)


params_list <- params_listB
n_scenarios <- length(params_list)

if (n_scenarios > length(cols)) {
  cols <- colorRampPalette(cols)(n_scenarios)
}

plot(x, params_list[[1]], 
     type = "p", pch = 16, cex = 1.3,
     col = cols[1],
     ylim = c(0, 1),
     xlab = "version",
     ylab = "true θ")

for(i in 2:length(params_list)){
  points(x, params_list[[i]], col = cols[i], pch = 16, cex = 1.3)
}

legend("topleft",
       legend = c("θ = (0.1, 0.125, 0.15, 0.175, 0.2)", "θ = (0.1, 0.15, 0.2, 0.25, 0.3)", "θ = (0.1, 0.175, 0.25, 0.325, 0.4)", "θ = (0.1, 0.2, 0.3, 0.4, 0.5)", "θ = (0.1, 0.225, 0.35, 0.475, 0.6)", "θ = (0.1, 0.25, 0.4, 0.55, 0.7)", "θ = (0.1, 0.275, 0.45, 0.625, 0.8)", "θ = (0.1, 0.3, 0.5, 0.7, 0.9)"),
       col = cols, pch = 16, cex = 0.8)


params_list <- params_listC
n_scenarios <- length(params_list)

if (n_scenarios > length(cols)) {
  cols <- colorRampPalette(cols)(n_scenarios)
}

plot(x, params_list[[1]], 
     type = "p", pch = 16, cex = 1.3,
     col = cols[1],
     ylim = c(0, 1),
     xlab = "version",
     ylab = "true θ")

for(i in 2:length(params_list)){
  points(x, params_list[[i]], col = cols[i], pch = 16, cex = 1.3)
}

legend("topleft",
       legend = c("θ = (0.1, 0.1, 0.9, 0.1, 0.1)", "θ = (0.1, 0.5, 0.9, 0.5, 0.1)", "θ = (0.1, 0.1, 0.1, 0.9, 0.1)", "θ = (0.1, 0.4, 0.7, 0.9, 0.1)"),
       col = cols, pch = 16, cex = 0.8)



params_list <- params_listD
params_list <- params_list[1:5]


n_scenarios <- length(params_list)

if (n_scenarios > length(cols)) {
  cols <- colorRampPalette(cols)(n_scenarios)
}

plot(x, params_list[[1]], 
     type = "p", pch = 16, cex = 1.3,
     col = cols[1],
     ylim = c(0, 1),
     xlab = "version",
     ylab = "true θ")

for(i in 2:length(params_list)){
  points(x, params_list[[i]], col = cols[i], pch = 16, cex = 1.3)
}

legend("topleft",
       legend = c("θ = (0.1, 0.1, 0.1, 0.1, 0.1)", "θ = (0.2, 0.2, 0.2, 0.2, 0.2)", "θ = (0.3, 0.3, 0.3, 0.3, 0.3)", "θ = (0.4, 0.4, 0.4, 0.4, 0.4)", "θ = (0.5, 0.5, 0.5, 0.5, 0.5)"),
       #legend = c("θ = (0.1, 0.1, 0.1, 0.1, 0.1)", "θ = (0.2, 0.2, 0.2, 0.2, 0.2)", "θ = (0.3, 0.3, 0.3, 0.3, 0.3)", "θ = (0.4, 0.4, 0.4, 0.4, 0.4)", "θ = (0.5, 0.5, 0.5, 0.5, 0.5)", "θ = (0.6, 0.6, 0.6, 0.6, 0.6)", "θ = (0.7, 0.7, 0.7, 0.7, 0.7)", "θ = (0.8, 0.8, 0.8, 0.8, 0.8)", "θ = (0.9, 0.9, 0.9, 0.9, 0.9)"),
       col = cols, pch = 16, cex = 0.8)



