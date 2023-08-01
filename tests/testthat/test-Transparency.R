test_that("PairedHeatmap works", {
  p = expect_silent(
    PairedHeatmap(
      cellmatch::early_endoderm_reference[1:200, 1:5] %>% as.matrix,
      cellmatch::early_endoderm_reference[1:200, 1:5] %>% as.matrix,
      genes = rownames(cellmatch::early_endoderm_reference[1:200,]),
      scores = runif(200)
    )
  )
})


test_that("EvaluateByGene works", {
  X = expect_message(
    EvaluateByGene(
      query = cellmatch::early_endoderm_reference[1:200, 1:5] %>% as.matrix,
      reference = cellmatch::early_endoderm_reference[1:200, ] %>% as.matrix,
      equivalents = colnames(cellmatch::early_endoderm_reference[1:5]),
      results_path = "temp"
    ), regexp = "geom_smooth"
  )
})


test_that("Mopup", {
  expect_silent(
    unlink("temp", recursive=TRUE)
  )
})

test_that("PlotMissingGenes works", {
  X = expect_silent({
    output = PlotMissingGenes(
      query = cellmatch::early_endoderm_reference[1:200, 1:5] %>% as.matrix,
      reference = cellmatch::early_endoderm_reference[1:200, ] %>% as.matrix,
      query_baseline = "DE",
      reference_baseline = "Gut.tube",
      query_mature = "Gut.tube.Thymus",
      reference_mature = "Gut.tube.Thyroid",
      results_path =  "temp"
    )
  })
})

test_that("Mopup", {
  expect_silent(
    unlink("temp", recursive=TRUE)
  )
})


