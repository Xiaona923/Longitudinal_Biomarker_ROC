#' Example Dataset
#'
#' A list of two example data frames used for demonstrating package functionality.
#'
#' @format A named list:
#' \describe{
#'   \item{data.long}{a simulation data in long format. Each subject may contribute multiple records
#'    (id: Subject id,
#'     Xt: Values of longitudinal biomarker, 
#'     vtime: Biomarker measurement time)}
#'   \item{data.short}{a simulation data contains subject-level information.
#'    Each subject only has one record.
#'    (id: Subject id,
#'     Y: Observed event time,
#'     delta: Event indicator. 1 = event, 0 = censored,
#'     Z: a binary covariate,
#'     Zcont: a continuous covariate}
#' }
#' @examples
#' data(example_data)
#' str(example_data)
"example_data"
