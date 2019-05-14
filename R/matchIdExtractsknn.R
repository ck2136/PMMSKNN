#' Match ID Function: create array of ID values closest matching to training or testing case
#' 
#' @param data          Full dataset. Contains both Training and Testing data. 
#' @param train_test    Column name indicating whether individual belongs to
#' Training or Testing dataset. Train = 1, Test = 2 by default.
#' @param patid         Column name indicating patient id.
#' @param formula       Formula indicating the variables used for matching.
#' (e.g. \code{ ~ var1 + var2 + var3 }).
#' @return              An k x klist$test (or k x klist$train) array of patient id numbers
#' @export
matchIdExtractsknn <- function(
                 data,
                 train_test = "train__test",
                 patid,
                 formula
                 ){

    # - - - - - - - - - - - - - - - - - - - - #
    # DATA SET GENERATION
    # - - - - - - - - - - - - - - - - - - - - #
    full <- data 
    train <- full %>%
        filter(train_test == 1) %>% 
        distinct_(.dots = patid, .keep_all=TRUE) %>%
        dplyr::select(patid, all.vars(formula))
    test <- full %>%
        filter(train_test == 2) %>% 
        distinct_(.dots = patid, .keep_all=TRUE) %>%
        dplyr::select(patid, all.vars(formula))


    # - - - - - - - - - - - - - - - - - - - - #
    # Extract training patient_id closest to test patient
    # - - - - - - - - - - - - - - - - - - - - #
    nnarraytest <- sapply(1:nrow(test), function(y){

                          temp  <- sapply(2:ncol(test), function(x){
                                              (train[,x] - c(test[y,x]))^2
                        })
                          sqrt(apply(data.frame(temp),1, sum))
                          full %>%
                              filter(train_test == 1) %>%
                              distinct_(.dots = patid)  %>%
                              bind_cols(diff = sqrt(apply(data.frame(temp),1, sum))) %>%
                              arrange(.data$diff) %>%
                              #head(n = 5) %>%
                              dplyr::select(patid) %>% unlist %>% as.vector
                 })

    # - - - - - - - - - - - - - - - - - - - - #
    # Extract training patient_id closest to individual training id 
    # - - - - - - - - - - - - - - - - - - - - #
    nnarraytrain <- sapply(1:nrow(train), function(y){

                               temp  <- sapply(2:ncol(train), function(x){
                                                   (train[-y,x] - c(train[y,x]))^2
                        })

                               #sqrt(apply(data.frame(temp),1, sum))
                               temp2 <- sapply(1:length(temp), function(x) {
                                                   append(temp[[x]],0,y-1)
                        })

                               full %>%
                                   filter(train_test == 1) %>%
                                   distinct_(.dots = patid)  %>%
                                   bind_cols(diff = sqrt(apply(data.frame(temp2),1, sum))) %>%
                                   arrange(.data$diff) %>%
                                   .[-1, ] %>%
                                   #head(n = 5) %>%
                                   dplyr::select(patid) %>% unlist %>% as.vector
                 })

    return(list(
                nnarraytrain = nnarraytrain,
                nnarraytest = nnarraytest
                )
    )

}
