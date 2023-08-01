
v1 = (-2:20) %>% setdiff(0)
i1 = c(-2, -1, 1, 1, 2)
f1 = function(x, index_omit = NULL){
  x = x[setdiff(seq_along(x), index_omit)]
  prod(x)
}
v2 = -2:2
i2 = -2:2
f2 = function(x, index_omit = NULL){
  x = x[setdiff(seq_along(x), index_omit)]
  -sum(dnorm(x, log = T))
}

test_that("CoordinateDescent works", {
  expect_equal(CoordinateDescent( FUN = f1, v1, i1, verbose = F)$x,  c(20, -2, 20, 20, 20))
  expect_equal(CoordinateDescent( FUN = f2, v2, i2, verbose = F)$x,  rep(0, 5))
})

test_that("ComputeSingleCoordPaths works", {
  expect_equal(colSums(ComputeSingleCoordPaths(FUN = f1, v1, i1)$objective),
               c( -414, -828,  828,  828,  414))
  expect_equal(colSums(ComputeSingleCoordPaths(FUN = f2, v2, i2)$objective),
               c( 42.97346, 50.47346, 52.97346, 50.47346, 42.97346), tol = 0.1)
})



test_that("OmitSingleCoords works", {
  expect_equal(OmitSingleCoords(FUN = f1, init = i1)$objective[1,],
               c( -2, -4,  4,  4,  2))
  expect_equal(OmitSingleCoords(FUN = f2, init = i2)$objective[1,],
               c(  6.675754, 8.175754, 8.675754, 8.175754, 6.675754 ), tol = 0.001)
})

test_that("SampleUniformly works", {
  expect_silent({
  SampleUniformly(FUN = f1, v1, i1, niter = 3)
  SampleUniformly(f2, v2, i2, niter = 3)
  })
})





