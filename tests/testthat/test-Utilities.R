
test_that("PlotEquivalentsTSNE works", {
  expect_silent(
    {
      p = PlotEquivalentsTSNE(
        embedding = matrix(runif(200), ncol = 2) %>% set_colnames(1:2),
        group_labels = rep(c("EPI", "DE"), each = 50),
        pairs = data.frame(
          query = c("ESC", "DE"),
          reference = c("EPI", "DE"),
          colors = c("red", "yellow"),
          stringsAsFactors = F)
      )
      print(p)
    }
  )
})


test_that("SelectInformativeGenes works", {
  expect_silent(SelectInformativeGenes(cellmatch::early_endoderm_reference))
})

test_that("AggregateByCluster works", {
  expect_silent(AggregateByCluster(matrix(rnorm(16), ncol = 4), c('a', "b", "a", "b"), method = "average"))
  expect_warning(AggregateByCluster(matrix(rnorm(16), ncol = 4), c(1, 2, 1, 2),        method = "average"))
})

test_that("string + works", {
  expect_equal("2"+"2", "22")
})

test_that("nwm works", {
  maximilian = 1:5 %>% setNames(letters[1:5]) %>% nwm
  expect_equal(maximilian, "e")
})


test_that("ReplaceNA works", {
  expect_equal(
    ReplaceNA(c(1, 2, NA, NA, 5)),
    c(1, 2, 0, 0, 5)
  )
  expect_equal(
    ReplaceNA(c(1, 2, NA, NA, 5), filler = ""),
    c(1, 2, "", "", 5)
  )
})

test_that("PadVector works", {
  expect_equal(
    PadVector(c(1, 2, 5), n = 5),
    c(1, 2, 5, "", "")
  )
  expect_equal(
    PadVector(c(1, 2, 5), n = 5, filler = 0),
    c(1, 2, 5, 0, 0)
  )
})
