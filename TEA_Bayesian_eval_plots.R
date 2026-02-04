remove(list = ls())


cols <- c("#4E79A7", "#59A14F", "#9C755F", "#B07AA1", "#76B7B2")

x <- 1:50
# first plot
plot(x, tea_stops,
     type = "l",
     col = cols[4],
     ylim = c(0,length(params)),
     xlab = "lambda",
     ylab = "new trial",
     main = ""
)

# add second line
lines(x, naive_stops,
      col = cols[5],
      type = "l")

# optional legend
legend("topright",
       legend = c("TEA", "naive"),
       col = cols[4:5],
       lty = 1)

