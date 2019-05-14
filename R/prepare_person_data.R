#' Creates data on the person level 
#' 
#' @param data A data frame in long form
#' @param out_time The target time for which prediction is needed. 
#' The number correspond to one of the knots. 
#' @param outcome -
#' @param time_var - 
#' @param pat_id - 
#' @param list_bs A list of objects, each of class \code{brokenstick} 
#' or class \code{brokenstick_export}. The list components with the 
#' name specified in \code{outcome} is taken as the external 
#' broken stick model.
#' @return A data frame with one row for each \code{pat_id}, with
#' a column named \code{yhat} containing the broken stick estimates
#' for each person at time point specified in \code{time_var}.
#' @examples 
#' timedata <- create_testdata("smocc")
#' persondata <- create_person_data(data = timedata, out_time = 14/12,
#'     outcome = "hdc", time_var = "age", pat_id = "id", 
#'     list_bs = PMMSKNN::smocc_bs)
#'@export
create_person_data <- function(data, 
                               out_time = 90,
                               outcome = "tug",
                               time_var = "time",
                               pat_id = "patient_id",
                               list_bs = NULL) {
  if (is.null(list_bs)) 
    stop("This function does not (yet) handle a NULL for list_bs")
  if (!outcome %in% names(list_bs)) 
    stop("Broken stick model not found for ", outcome)
  vars <- c(outcome, time_var, pat_id)
  found <- vars %in% names(data)
  if (any(!found)) stop("Variables not found: ", vars[!found])
  
  # find knot at out_time
  knots <- get_knots(list_bs[[outcome]])
  kk <- which(round(knots, 4L) %in% round(out_time, 4L))
  if (length(kk) == 0L) 
    stop("Knot at ", round(out_time, 4L), " not found.")
  
  # calculate brokenstick estimate at out_time knot
  # and reduce to person level
  data %>% 
    group_by(.data[[pat_id]]) %>% 
    mutate(yhat = predict(list_bs[[outcome]], 
                          y = .data[[outcome]],
                          x = .data[[time_var]],
                          at = "knots", 
                          output = "vector")[kk]) %>% 
    slice(1L) %>% 
    ungroup()
}
