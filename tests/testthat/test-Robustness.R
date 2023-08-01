test_that("GetDefaultDistances all work", {
  x = rpois(16, lambda = 10) %>% matrix(nrow = 4) %>% set_rownames(1:4) %>% set_colnames(1:4)
  y = rpois(16, lambda = 10) %>% matrix(nrow = 4) %>% set_rownames(1:4) %>% set_colnames(1:4)
  for( ii in names(GetDefaultDistances())){
    cat(ii, "\n")
    if(ii == "correlation_distance_slow"){
      testthat::expect_warning(GetDefaultDistances()[[ii]](x, y))
    } else {
      testthat::expect_silent(GetDefaultDistances()[[ii]](x, y))

    }
  }
})


test_that("Faster correlation_distance equals slower", {
  x = rpois(16, lambda = 10) %>% matrix(nrow = 4) %>% set_rownames(1:4) %>% set_colnames(1:4)
  y = rpois(16, lambda = 10) %>% matrix(nrow = 4) %>% set_rownames(1:4) %>% set_colnames(1:4)
  expect_warning(
    expect_equal(
      GetDefaultDistances()[["correlation_distance"]](x, y) %>% setNames(NULL),
      GetDefaultDistances()[["correlation_distance_slow"]](x, y)
    )
  )

})
