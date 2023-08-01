test_that("species conversion works for character vectors", {
  expect_equal(
    get_ortholog(c("Foxn1", "Zfp750", "GARBAGE"), from = "mouse", to = "human") %>% setNames(NULL),
               c("FOXN1", "ZNF750", NA)
    )
  expect_equal(
    get_ortholog(c("FOXN1", "ZNF750", "GARBAGE"), from = "human", to = "mouse") %>% setNames(NULL),
    c("Foxn1", "Zfp750", NA)
  )
})

test_that("merge dupes works ", {
  rn = rep(letters[1:4], each = 2) %>% c("f", "g")
  expect_output({
    cellmatch::merge_dupe_rows(X = matrix(rnorm(100), 10) %>% set_rownames(rn),
                               row.names = rn, verbose = T)
  })
})

test_that("species conversion works for expression matrices", {
  X = cellmatch::early_endoderm_reference %>% as.matrix %>% head(100)
  expect_equal(
    cellmatch::convert_species_rownames(expr = X,     from = "mouse", to = "human") %>% rownames,
    X %>% rownames %>% get_ortholog(from = "mouse", to = "human") %>% (function(x) x[!is.na(x)])
  )
})

