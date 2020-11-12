context("Kalman filtering with NAs")

test_that('kalman', {

  data[1, 1] <- NA
  kf <- KalmanFilter(data, H, Q, R, F, x0, P0)

  ## Test that presence of NAs does not make the Kalman filtering procedure in
  ## correct in some way. Reference values are based on FKF package.
  expect_that(all(kf$xF - matrix(c(1.947193543, 1.0170326503,
                                   1.710003692, 0.4556847315,
                                   1.726467898, 0.3234284077),
                                 3, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kf$Pf[,,1] - matrix(c(0.64776479220, -0.09381059693,
                                        -0.09381059693, 0.58801949583),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kf$Pf[,,2] - matrix(c(0.4880801052, -0.1505185232,
                                        -0.1505185232, 0.5495135325),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kf$Pf[,,3] - matrix(c(0.4784520052, -0.1499921822,
                                        -0.1499921822, 0.5491835627),
                                      2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(abs(kf$loglik - -16.04131559) < 1e-7, is_true())

})
