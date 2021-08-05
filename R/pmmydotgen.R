#' PMM \eqn{\dot{Y}} Generating Function: Create Y dots through Predictive Mean Matching
#' 
#' @param df Training data
#' @param formula A string with the name of the splitting variables
#' @param pat_id     Name of the variable with patient identification
#' @param m Integer value indicating the number of replications to obtain
#' \eqn{\beta*} and \eqn{\dot{y}}.
#' @param seed Integer value indicating seed number for random generation
#' of \eqn{\sigma*} and \eqn{\beta*}
#' @param dftest  Test data.frame object 
#' 
#' @return A data frame that includes \eqn{\dot{y}}.
#' 
#' @export
pmmydotgen <- function(df,  # df should be with id column and all other covariates that are going to be included in the formula
                       formula, 
                       pat_id,
                       m, 
                       seed=1234,
                       dftest=NULL
                       ){
    if(pat_id != "id"){
        df <- df %>%
            dplyr::rename(id = pat_id)

        if(!is.null(dftest)){
            dftest <- dftest %>%
                dplyr::rename(id = pat_id)
        }

    }

    if(is.null(dftest)){
        temp <- lapply(1:m, function(y){
                           ydots <- sapply(1:nrow(df), function(x) {
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 1. LM
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      lm <- lm(formula, data=df[-x, ])
                                      fs <- summary(lm)

                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 2. Draw \sigma^2*
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      dfree <- length(lm$fitted.values) - length(lm$coefficients)
                                      set.seed(seed+y)
                                      sigmastar <- c(fs$residuals %*% fs$residuals) / rchisq(1,dfree) 

                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 3. Draw beta* from mvn N(\hat{\beta}, \sigma*)
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      X <- model.matrix(lm)
                                      covmat <- sigmastar[1] * solve(crossprod(X))
                                      #covmat <- make.positive.definite(covmat)
                                      set.seed(seed+y)
                                      betastar <- mvrnorm(n=1,mu=lm$coefficients,Sigma=covmat)
                                      #lm$coef # compare with original

                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 4. calculate \hat{y}_obs and \hat{y}_mis
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # \hat{y}obs ez 
                                      Yobs <- X %*% coef(lm)
                                      # \hat{y}mis ez 
                                      Xmis <- model.matrix(formula(paste0("~", paste0(formula[3]))), data=df[x,] 
                                      ) 
                                      Ymis <- Xmis %*% betastar
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 5. Find delta = abs(\hat{y}obs - \hat{y}mis)
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 6. Randomly sample one from 3 smallest
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      df[-x,] %>%
                                          dplyr::select(id, !!sym(paste0(formula[2]))) %>%
                                          mutate("id = as.integer(id)") %>%
                                          filter(.data$id %in% (df[-x,] %>% 
                                                          dplyr::select(id) %>%
                                                          mutate("id = as.integer(id)") %>%
                                                          bind_cols(
                                                                    # find delta
                                                                    diff = as.vector(abs(Yobs - Ymis[1]))
                                                                    ) %>% 
                                                          arrange(diff) %>%
                                                          dplyr::select(id) %>%
                                                          # sample 5
                                                          head(n=5) %>% unlist %>% as.vector %>%
                                                          sample(., 1))
                                          ) %>%
                                          dplyr::select(!!sym(paste0(formula[2]))) %>% unlist %>% as.vector
           })
                       })

        # take the mean of the 10 values
        temp <- data.frame(temp)
        colnames(temp) <- paste0("ydot", 1:m)
        return(
               bind_cols(
                         temp,
                         ydotavg = apply(data.frame(temp), 1, mean)
               )
        )

    # for testing set
    } else {
        temp <- lapply(1:m, function(y){
                           ydots <- sapply(1:nrow(dftest), function(x) {
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 1. LM
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      lm <- lm(formula, data=df)
                                      fs <- summary(lm)

                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 2. Draw \sigma^2*
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      dfree <- length(lm$fitted.values) - length(lm$coefficients)
                                      set.seed(seed+y)
                                      sigmastar <- c(fs$residuals %*% fs$residuals) / rchisq(1,dfree) 

                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 3. Draw beta* from mvn N(\hat{\beta}, \sigma*)
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      X <- model.matrix(lm)
                                      covmat <- sigmastar[1] * solve(crossprod(X))
                                      #covmat <- make.positive.definite(covmat)
                                      set.seed(seed+y)
                                      betastar <- mvrnorm(n=1,mu=lm$coefficients,Sigma=covmat)
                                      #lm$coef # compare with original

                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 4. calculate \hat{y}_obs and \hat{y}_mis
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # \hat{y}obs ez 
                                      Yobs <- X %*% coef(lm)
                                      # \hat{y}mis ez 
                                      Xmis <- model.matrix(formula(paste0("~", paste0(formula[3]))), data=dftest[x,] 
                                      ) 
                                      Ymis <- Xmis %*% betastar
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 5. Find delta = abs(\hat{y}obs - \hat{y}mis)
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      # 6. Randomly sample one from 3 smallest
                                      # - - - - - - - - - - - - - - - - - - - - - #
                                      df %>%
                                          dplyr::select(id, !!sym(paste0(formula[2]))) %>%
                                          mutate("id = as.integer(id)") %>%
                                          filter(.data$id %in% (df %>% 
                                                          dplyr::select(id) %>%
                                                          mutate("id = as.integer(id)") %>%
                                                          bind_cols(
                                                                    # find delta
                                                                    diff = as.vector(abs(Yobs - Ymis[1]))
                                                                    ) %>% 
                                                          arrange(diff) %>%
                                                          dplyr::select(id) %>%
                                                          # sample 3
                                                          head(n=5) %>% unlist %>% as.vector %>%
                                                          sample(., 1))
                                          ) %>%
                                          dplyr::select(!!sym(paste0(formula[2]))) %>% unlist %>% as.vector
           })
                       })

        # take the mean of the 10 values
        temp <- data.frame(temp)
        colnames(temp) <- paste0("ydot", 1:m)
        return(
               bind_cols(
                         temp,
                         ydotavg = apply(data.frame(temp), 1, mean)
               )
        )

    }
}
