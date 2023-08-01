library(magrittr)
test_that("Model selection works", {
  expect_silent({
    SelectModelThoroughly(
      cellmatch::early_endoderm_reference[1:1000, 1:5] %>% as.matrix,
      cellmatch::early_endoderm_reference[1:1000, 1:5] %>% as.matrix,
      verbose = F
    )
  })
})

test_that("alternate model plotting works", {
  expect_silent( {
    neighbors = MutateModel(
      cellmatch::early_endoderm_reference[1:1000, 1:5] %>% as.matrix,
      cellmatch::early_endoderm_reference[1:1000, 1:5] %>% as.matrix,
      init = colnames(cellmatch::early_endoderm_reference[1:5])
    )
    PlotNeighboringModels(neighbors, results_path = "temp")
  })
})

test_that("alternate model plotting works", {
  expect_silent( {
    neighbors = MutateModel(
      cellmatch::early_endoderm_reference[1:1000, 1:5] %>% as.matrix,
      cellmatch::early_endoderm_reference[1:1000, 1:5] %>% as.matrix,
      init = colnames(cellmatch::early_endoderm_reference[1:5])
    )
    PlotNeighboringModels(neighbors, results_path = "temp")
  })
})


test_that("RunCellMatch  works", {
  RunCellMatch(
    cellmatch::early_endoderm_reference[1:1000, 1:5] %>% as.matrix,
    cellmatch::early_endoderm_reference[1:1000, 1:5] %>% as.matrix,
    results_path = "temp", K = 500, num_init = 2, verbose = F
  )
})



test_that("Mopup", {
  expect_silent(
    unlink("temp", recursive=TRUE)
  )
})



