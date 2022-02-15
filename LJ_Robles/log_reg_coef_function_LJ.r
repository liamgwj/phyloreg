# log_reg_coef

## requires that columns of incidence correspond to rows of phydist

function (incidence, phydist){
    
    if (nrow(incidence) == 1) {
        
        # get names of genera with '1' compatibility
        l <- list(names(incidence[, apply(incidence, 2,
                                          function(x) (x == 1))]))
        
        names(l) <- rownames(incidence)
        
    } else {
        
        l <- lapply(apply(apply(incidence, 2, function(x) (x == 1)),
                          1, which),
                    names)
        
    }
    # choose a focal host
    focal <- lapply(l, function(x) sample(x, 1))
    
    df <- do.call(rbind,
                  lapply(names(focal),
                         function(x) data.frame(colnames(incidence),
                                                focal[[x]],
                                                phydist[colnames(incidence),
                                                        focal[[x]]],
                                                as.numeric(incidence[x, ]),
                                                x)))
    
    colnames(df) <- c("tohost", "fromhost", "phydist", "suscept", 
                      "incidence")
    
    log_out <- stats::glm(df$suscept ~ df$phydist,
                          family = stats::binomial(link = "logit"))
    
    stats_log_out <- stats::coef(summary(log_out))
    
    conf_int <- stats::confint.default(log_out)
    
    lrcoeffs <- data.frame(log_out$coefficients[1],
                           stats_log_out[1, 2],
                           stats_log_out[1, 3],
                           stats_log_out[1, 4],
                           conf_int[1, 1],
                           conf_int[1, 2],
                           log_out$coefficients[2],
                           stats_log_out[2, 2],
                           stats_log_out[2, 3],
                           stats_log_out[2, 4],
                           conf_int[2, 1],
                           conf_int[2, 2])
    
    colnames(lrcoeffs) <- c("intercept", "Std. Error", "z value", "Pr(>|z|)",
                            "2.5 %", "97.5 %", "slope", "Std. Error",
                            "z value", "Pr(>|z|)", "2.5 %", "97.5 %")
    
    rownames(lrcoeffs) <- NULL
    
    return(lrcoeffs)
    
}
