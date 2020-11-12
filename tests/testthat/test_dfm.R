library(dynfactoR)

x <- matrix(rnorm(50*10), 50, 10)

fit <- dfm(x, 2, 2, 1)
fit <- dfm(x, 2, 2, 2)
fit <- dfm(x, 2, 1, 1)

x <- matrix(rnorm(200*10), 200, 10)
W <- as.logical(matrix(rbinom(200*10, 1, 0.01), 200, 10))
x[W] <- NA
fit <- dfm(x, 2, 2, 1, lower_set = 1e-3, max_iter = 1000)
fit <- dfm(x, 2, 2, 2, lower_set = 1e-3, max_iter = 1000)
fit <- dfm(x, 2, 1, 1, lower_set = 1e-3, max_iter = 1000)

x <- matrix(rnorm(200*10), 200, 10)
y <- dfmMS(x[,1:5], nf = 2)
