#' Create datasets for testing purposes
#' 
#' @param dataset_name A string. Currently, only 
#' \code{dataset_name = "smocc"} is supported.
#' @return A data frame in the long form
#' @examples 
#' full <- create_testdata("smocc")
#' nrow(full)
#' @export
create_testdata <- function(dataset_name) {
  if (dataset_name == "smocc") {
    
    smocc <- PMMSKNN::smocc
    # include time-invariant covariates into long data
    person <- smocc[["child"]] %>% 
      dplyr::select(.data$id, .data$edu, .data$twin, .data$agem, 
                    .data$smo, .data$hgtm, .data$hgtf, .data$hdc_0)
    full <- smocc[["time"]] %>% 
      dplyr::select(.data$id, .data$rec, .data$age, .data$sex, 
                    .data$ga, .data$bw, .data$hgt, .data$wgt, 
                    .data$hdc, .data$bmi, 
                    .data$hgt.z, .data$wgt.z, .data$hdc.z, .data$bmi.z) %>% 
      left_join(person, by = "id") %>% 
      dplyr::select(.data$id, .data$rec, .data$age, .data$sex, 
                    .data$ga, .data$bw, .data$edu:.data$hgtf, 
                    .data$hgt:.data$bmi.z, .data$hdc_0)
    
    # for testing, make random selection of 50 persons
    set.seed(98881)
    ids <- sample(unique(full$id), 50)
    full <- full %>% 
      filter(.data$id %in% ids)
    
    # Train and Test split for all TKA outcomes: create
    full <- baselinemk(full, pat_id = "id", time_var = "age")
    
    # divide into train and test subsets
    ntrain <- as.integer(length(ids)/ 3 * 2)
    train_test <- c(rep(1L, ntrain), rep(2L, length(ids) - ntrain))
    df <- tibble(id = ids, train_test = train_test)
    full <- left_join(full, df, by = "id")
    
    return(full)
  }
}
