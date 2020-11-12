context("Field to array conversion")

test_that('kim', {

  kimF <- KimFilterCpp(data, R, Q, H1, F1, x0, P0, p)

  array_Pa <- generateArray(kimF$Pa)
  array_x <- generateArray(kimF$x)
  array_P <- generateArray(kimF$P)

  expect_true(inherits(kimF$Pa, "matrix"))
  expect_true(inherits(kimF$x, "matrix"))
  expect_true(inherits(kimF$P, "matrix"))

  expect_true(inherits(kimF$Pa[1,1], "list"))
  expect_true(inherits(kimF$x[1,1], "list"))
  expect_true(inherits(kimF$P[1,1], "list"))

  expect_true(inherits(kimF$Pa[1,1][[1]], "array"))
  expect_true(inherits(kimF$x[1,1][[1]], "matrix"))
  expect_true(inherits(kimF$P[1,][[1]], "array"))

  expect_true(identical(kimF$Pa[1,1][[1]], array_Pa[,,,1,1]))
  expect_true(identical(kimF$Pa[2,1][[1]], array_Pa[,,,2,1]))

  expect_true(identical(kimF$x[1,1][[1]], array_x[,,1,1]))
  expect_true(identical(kimF$x[2,1][[1]], array_x[,,2,1]))
  expect_true(identical(kimF$x[1,2][[1]], array_x[,,1,2]))
  expect_true(identical(kimF$x[2,2][[1]], array_x[,,2,2]))

  expect_true(identical(kimF$P[1,1][[1]], array_P[,,,1,1]))
  expect_true(identical(kimF$P[2,1][[1]], array_P[,,,2,1]))
  expect_true(identical(kimF$P[1,2][[1]], array_P[,,,1,2]))
  expect_true(identical(kimF$P[2,2][[1]], array_P[,,,2,2]))

})

