context("Kim filtering and smoothing")

test_that('kim', {

  ## If input matrices and input probabilities are equal for each state, then
  ## Kim filter and smoother should lead to the same results as Kalman filter
  ## and smoother.

  ## Kim filtering
  kimF <- KimFilterCpp(data, R, Q, H1, F1, x0, P0, p)

  expect_that(abs(sum(kimF$result) - -17.40283081) < 1e-7, is_true())

  expect_that(all(t(kimF$xA[,1,]) - matrix(c(1.35635067, 0.77081876,
                                             1.55174759, 0.46023087,
                                             1.68756582, 0.33374625),
                                           3, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(t(kimF$xA[,2,]) - matrix(c(1.35635067, 0.77081876,
                                             1.55174759, 0.46023087,
                                             1.68756582, 0.33374625),
                                           3, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimF$Pa[1,1][[1]][,,1] - matrix(c(0.46616058, -0.16948804,
                                                    -0.16948804, 0.55648347),
                                                  2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimF$Pa[1,1][[1]][,,2] - matrix(c(0.47692309, -0.15019802,
                                                    -0.15019802, 0.54950433),
                                                  2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimF$Pa[1,1][[1]][,,3] - matrix(c(0.47780773, -0.1498213,
                                                    -0.1498213, 0.54913824),
                                                  2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimF$Pa[2,1][[1]][,,1] - matrix(c(0.46616058, -0.16948804,
                                                    -0.16948804, 0.55648347),
                                                  2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimF$Pa[2,1][[1]][,,2] - matrix(c(0.47692309, -0.15019802,
                                                    -0.15019802, 0.54950433),
                                                  2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimF$Pa[2,1][[1]][,,3] - matrix(c(0.47780773, -0.1498213,
                                                    -0.1498213, 0.54913824),
                                                  2, 2, byrow = TRUE) < 1e-7), is_true())

  ## Check that wrapper function for Kim filter returns the same results
  kimFC <- KimFilter(data, R, Q, H1, F1, x0, P0, p)

  expect_identical(generateArray(kimF$Pa), kimFC$Pa)
  expect_identical(generateArray(kimF$x), kimFC$x)
  expect_identical(generateArray(kimF$P), kimFC$P)

  ## Kim smoothing
  kimS <- KimSmoother(F1, kimF$xA, generateArray(kimF$Pa),
                      generateArray(kimF$x),
                      generateArray(kimF$P),
                      kimF$stateP, kimF$stateP_fut, p)

  expect_that(all(kimS$xF - matrix(c(1.51225349, 0.77961369,
                                     1.69965529, 0.46657133,
                                     1.68756582, 0.33374625),
                                   3, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimS$Pf[,,1] - matrix(c(0.43626195, -0.17509815,
                                          -0.17509815, 0.54651832),
                                        2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimS$Pf[,,2] - matrix(c(0.44567712, -0.15752478,
                                          -0.15752478, 0.53925906),
                                        2, 2, byrow = TRUE) < 1e-7), is_true())

  expect_that(all(kimS$Pf[,,3] - matrix(c(0.47780773, -0.1498213,
                                          -0.1498213, 0.54913824),
                                        2, 2, byrow = TRUE) < 1e-7), is_true())


})
