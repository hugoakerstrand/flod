#' Convert a flowSet into a tibble
#'
#' @description `to_tibble()` creates a [tibble::tibble()] from slots of a
#'   [flowCore::flowSet()]; an S4 class that stores .fcs file data.
#'
#' @param data A flowSet object.
#' @param ... What slots to extract? <kw: default = "all"> //  One or more of
#'   `"all"`, `"sample_names"`, `"exprs"`, `"keywords"`, or `"meta_data"`.
#'
#' @returns A tibble containing the requested data components. The function will
#'   error if `data` is not a flowSet or if invalid keywords are provided.
#'
#'   Additional keywords are ignored if combined with `all`.
#'
#' @export
to_tibble <- function(data, ...) {
  input_keywords <- c(...)

  # Input validation
  if(!any(class(data) == "flowSet")) {
    cli::cli_abort("{.obj data} must have class {.cls flowSet}, not {.cls {class(data)[1]}}.")
  }

  valid_keywords <- c("all", "sample_names", "exprs", "keywords", "meta_data")

 if(!missing(...) && !all(input_keywords %in% valid_keywords)) {
    cli::cli_abort("`...` must be a combination of {.val {valid_keywords}}, not {.val {(input_keywords)}}.")
 }

  # Handle default "all" case to access complete data
  if (missing(...) || "all" %in% input_keywords) {
    input_keywords <- setdiff(valid_keywords, "all")
  }

  # Collect data using input_keywords
  result_list <- list()

  for (kw in input_keywords) {
    result_list[[kw]] <- switch(kw,
                                sample_names = get_sample_names(data),
                                exprs = get_exprs(data),
                                keywords = get_keywords(data),
                                meta_data = get_meta_data(data)
                        )
  }

  do.call(tibble::tibble, result_list)
}

get_sample_names <- function(data) {
  flowCore::pData(data) |> rownames()
}

get_exprs <- function(data) {
  lapply(1:length(data), function(i) flowCore::exprs(data[[i]]))
}

get_keywords <- function(data) {
  lapply(1:length(data), function(i) flowCore::keyword(data[[i]]))
}

get_meta_data <- function(data) {
  flowCore::pData(data)
}
