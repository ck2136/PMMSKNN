#' Function to add baseline variable
#'
#' `baseline` column is created based on the full dataframe
#' provided by the user. If the `time_var` column contains 
#' more than one potential baseline values (i.e. values <= 0)
#' then `baselinemk()` will arrange the times according to 
#' closest value to 0 then assign that as baseline = 1 
#' while the other potential baseline values will be assigne -1. 
#'
#' @param dftotf data frame to transform (i.e. data frame where a new column - baseline column - will be added).
#'   The data frame must contain specified `pat_id` and
#'   `time_var`.
#' @param pat_id string of characters representing the column
#'   of patient id.
#' @param time_var string of characters representing the column
#'   of time.
#'   
#' @return A data frame that contains a new `baseline` column indicating the baseline observation (i.e. row).
#' 
#' @export
baselinemk <- function(dftotf,
                       pat_id = "patient_id",
                       time_var = "time"
                       ){   
    dftotf <- dftotf %>%
        left_join(
                  dftotf %>%
                      dplyr::select(!!dplyr::sym(pat_id), !!dplyr::sym(time_var)) %>%
                      dplyr::filter(!!dplyr::sym(time_var) <= 0) %>%
                      dplyr::arrange(!!dplyr::sym(pat_id), dplyr::desc(!!dplyr::sym(time_var))) %>%
                      dplyr::group_by(!!dplyr::sym(pat_id)) %>%
                      dplyr::mutate(baseline = if_else(row_number() == 1, 1, -1)) ,
                  by = c(pat_id, time_var)
                  ) %>%
    mutate(baseline = if_else(is.na(.data$baseline), 0,.data$baseline))
}