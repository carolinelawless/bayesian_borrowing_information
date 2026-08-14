plot(lambda_vector, proba_TEA_0, type = "l",
     col = "blue",
     xlab = "λ",
     ylab = "Probability")

lines(lambda_vector, proba_TEA_0.5, col = "red")
lines(lambda_vector, proba_TEA_1, col = "purple")
lines(lambda_vector, proba_TEA_EB_0.5, col = "pink")
lines(lambda_vector, proba_TEA_EB_1, col = "black")

legend("topright",
       legend = c("TEA_0", "TEA_0.5", "TEA_1", "TEA_EB_0.5", "TEA_EB_1"),
       col = c("blue", "red", "purple", "pink", "black"),
       lty = 1)

