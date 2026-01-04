# Input Validation Tests -------------------------------------------------------
test_that("to_tibble errors with non-flowSet input", {
  expect_error(
    to_tibble(data.frame(x = 1:5)),
    class = "rlang_error"
  )

  expect_error(
    to_tibble(list(a = 1, b = 2)),
    class = "rlang_error"
  )

  expect_error(
    to_tibble(NULL),
    class = "rlang_error"
  )

  expect_error(
    to_tibble(),
    class = "missingArgError"
  )
})

test_that("to_tibble errors with invalid keyword input", {
 data(GvHD, package = "flowCore")

  expect_error(
    to_tibble(GvHD, "invalid_keyword"),
    class = "rlang_error"
  )

  expect_error(
    to_tibble(GvHD, "all", "bad"),
    class = "rlang_error"
  )
})

# Default Behavior Tests -------------------------------------------------------
test_that("to_tibble returns alldata by default", {
  data(GvHD, package = "flowCore")

  result <- GvHD |> to_tibble()

  expect_s3_class(result, "tbl_df")
  expect_shape(result, dim = c(35,4))
  expect_named(result, c("sample_names", "exprs", "keywords", "meta_data"))
})

test_that("to_tibble returns alldata with 'all' keyword", {
  data(GvHD, package = "flowCore")
  result <- GvHD |> to_tibble("all")

  expect_s3_class(result, "tbl_df")
  expect_shape(result, dim = c(35,4))
  expect_named(result, c("sample_names", "exprs", "keywords", "meta_data"))
})

test_that("to_tibble returns only sample_names when requested", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "sample_names")

  expect_s3_class(result, "tbl_df")
  expect_named(result, "sample_names")
  expect_shape(result, dim = c(35,1))
  expect_type(result$sample_names, "character")
  expect_identical(result$sample_names, rownames(flowCore::pData(GvHD)))
})

test_that("to_tibble returns only exprs when requested", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "exprs")

  expect_s3_class(result, "tbl_df")
  expect_named(result, "exprs")
  expect_shape(result, dim = c(35,1))
  expect_type(result$exprs, "list")
  expect_length(result$exprs, length(GvHD))

  # Check that each element is a matrix
  expect_true(all(sapply(result$exprs, is.matrix)))
})

test_that("to_tibble returns only keywords when requested", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "keywords")

  expect_s3_class(result, "tbl_df")
  expect_named(result, "keywords")
  expect_shape(result, dim = c(35,1))
  expect_type(result$keywords, "list")
  expect_equal(length(result$keywords), length(GvHD))

  # Check that each element is a list of keywords
  expect_true(all(sapply(result$keywords, is.list)))
})

test_that("to_tibble returns only meta_data when requested", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "meta_data")

  expect_s3_class(result, "tbl_df")
  expect_named(result, "meta_data")
  expect_shape(result, dim = c(35,1))

  # Check that meta_data matches pData
  expect_equal(result[[1]], flowCore::pData(GvHD))
})

test_that("to_tibble handles multiple keywords", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "sample_names", "exprs")

  expect_s3_class(result, "tbl_df")
  expect_shape(result, dim = c(35,2))
  expect_named(result, c("sample_names", "exprs"))
})

test_that("to_tibble handles three keywords", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "sample_names", "exprs", "keywords")

  expect_s3_class(result, "tbl_df")
  expect_shape(result, dim = c(35,3))
  expect_named(result, c("sample_names", "exprs", "keywords"))
})

test_that("to_tibble handles all keywords except one", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "sample_names", "exprs", "keywords")

  expect_s3_class(result, "tbl_df")
  expect_false("meta_data" %in% names(result))
  expect_true(all(c("sample_names", "exprs", "keywords") %in% names(result)))
})

# Output Structure Tests -------------------------------------------------------

test_that("to_tibble output has correct structure", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD)

  # Check it's a tibble
  expect_s3_class(result, "tbl_df")

  # Check number of rows equals number of samples
  expect_equal(nrow(result), length(GvHD))

  # Check each column type
  expect_type(result$sample_names, "character")
  expect_type(result$exprs, "list")
  expect_type(result$keywords, "list")
  expect_type(result$meta_data, "list")
})

test_that("sample_names matches flowSet sample names", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "sample_names")

  expect_equal(
    result$sample_names,
    rownames(flowCore::pData(GvHD))
  )
})

test_that("exprs contains correct matrices", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "exprs")

  # Check number of matrices
  expect_equal(length(result$exprs), length(GvHD))

  # Check each matrix matches the flowFrame exprs
  for (i in seq_along(result$exprs)) {
    expect_equal(
      result$exprs[[i]],
      flowCore::exprs(GvHD[[i]])
    )
  }
})

test_that("keywords contains correct keyword lists", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "keywords")

  # Check number of keyword lists
  expect_equal(length(result$keywords), length(GvHD))

  # Check each keyword list matches the flowFrame keywords
  for (i in seq_along(result$keywords)) {
    expect_equal(
      result$keywords[[i]],
      flowCore::keyword(GvHD[[i]])
    )
  }
})

test_that("meta_data matches pData", {
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "meta_data")

  expect_equal(
    result[[1]],
    flowCore::pData(GvHD)
  )
})

# Edge Cases -------------------------------------------------------------------

test_that("to_tibble works with single sample flowSet", {
  data(GvHD, package = "flowCore")
  single_sample <- GvHD[1]
  result <- to_tibble(single_sample)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1)
  expect_named(result, c("sample_names", "exprs", "keywords", "meta_data"))
})

test_that("to_tibble handles duplicate keyword specifications", {
  data(GvHD, package = "flowCore")
  # Silently ignored
  expect_no_error(
    to_tibble(GvHD, "exprs", "exprs")
  )
})

test_that("to_tibble works with all keyword mixed with specific keywords", {
  # When "all" is present, it should return alldata
  data(GvHD, package = "flowCore")
  result <- to_tibble(GvHD, "all", "sample_names")

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("sample_names", "exprs", "keywords", "meta_data"))
})
