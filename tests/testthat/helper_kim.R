Q <- matrix(c(1.0, 0.1,
              0.1, 1.0), 2, 2, byrow = TRUE)
R <- matrix(c(1.0, 0.2, 0.1,
              0.2, 0.8, 0.5,
              0.1, 0.5, 1.2), 3, 3, byrow = TRUE)
H <- matrix(c(1.0, 0.7,
             0.5, 0.7,
             0.8, 0.1), 3, 2, byrow = TRUE)
F <- matrix(c(0.6, 0.2,
              0.1, 0.3), 2, 2, byrow = TRUE)
x0 <- rep(1, 2)
P0 <- diag(1, 2)

data <- matrix(c(1.04, 2.20, 3.12,
                 1.11, 2.33, 3.34,
                 1.23, 2.21, 3.45), 3, 3, byrow = TRUE)

F1 <- array(F, c(2, 2, 2))
H1 <- array(H, c(3, 2, 2))
p <- matrix(rep(0.5, 4), 2, 2)
