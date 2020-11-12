context("Kalman filtering and smoothing")

test_that('kalman', {

  ## Test that results from Kalman filtering procedure are
  ## equivalent to those that could be obtained from FKF
  ## or pykalman.
  kf <- KalmanFilter(data, H, Q, R, F, x0, P0)

  expect_that(all(kf$xF - matrix(c(1.35635067, 0.77081876,
                                   1.55174759, 0.46023087,
                                   1.68756582, 0.33374625),
                                 3, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kf$xP - matrix(c(1.0000000, 1.0000000,
                                   0.9679742, 0.3668807,
                                   1.0230947, 0.2932440,
                                   1.0792887, 0.2688805),
                                 4, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kf$Pf[,,1] - matrix(c(0.46616058, -0.16948804,
                                        -0.16948804, 0.55648347),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kf$Pf[,,2] - matrix(c(0.47692309, -0.15019802,
                                        -0.15019802, 0.54950433),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kf$Pf[,,3] - matrix(c(0.47780773, -0.1498213,
                                        -0.1498213, 0.54913824),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(abs(kf$loglik - -17.40283081) < 1e-7, is_true())

  ## Test that results from Kalman smoothing procedure are
  ## equivalent to those that could be obtained from pykalman.
  ks <- KalmanSmoother(F, H, R, kf$xF, kf$xP, kf$Pf, kf$Pp)

  expect_that(all(ks$xS - matrix(c(1.51225349, 0.77961369,
                                   1.69965529, 0.46657133,
                                   1.68756582, 0.33374625),
                                 3, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(ks$Ps[,,1] - matrix(c(0.43626195, -0.17509815,
                                        -0.17509815, 0.54651832),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(ks$Ps[,,2] - matrix(c(0.44567712, -0.15752478,
                                        -0.15752478, 0.53925906),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(ks$Ps[,,3] - matrix(c(0.47780773, -0.1498213,
                                        -0.1498213, 0.54913824),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())


})
