contrastsTest <- function(formula, phi.formula,
                          contrasts_DA = NULL,
                          contrasts_DV = NULL,
                          data,
                          link = "logit",
                          phi.link = "logit",
                          sample_data = NULL,
                          taxa_are_rows = TRUE,
                          filter_discriminant = TRUE,
                          fdr_cutoff = 0.05,
                          fdr = "fdr",
                          inits = NULL,
                          try_only = NULL,
                          ...) {
  
  if (is.null(contrasts_DA) && is.null(contrasts_DV)) {
    stop("Must include at least one of contrasts_DA or contrasts_DV!")
  }
  
  num_DA <- length(contrasts_DA)
  num_DV <- length(contrasts_DV)
  
  # Record call
  call <- match.call(expand.dots = TRUE)
  # Record mu link
  link <- match.arg(link, choices = "logit")
  # Record phi link
  phi.link <- match.arg(phi.link, choices = c("fishZ", "logit"))
  
  # Convert phyloseq objects
  if ("phyloseq" %in% class(data)) {
    # Set up response
    taxanames <- phyloseq::taxa_names(data)
  } else if (is.matrix(data) || is.data.frame(data)) {
    
    # use phyloseq
    OTU <- phyloseq::otu_table(data, taxa_are_rows = taxa_are_rows)
    
    # Make sample data
    sampledata <- phyloseq::sample_data(data.frame(
      sample_data,
      row.names = phyloseq::sample_names(OTU)
    ))
    
    # Make phyloseq object
    data <- phyloseq::phyloseq(OTU, sampledata)
    # Set up response
    taxanames <- phyloseq::taxa_names(data)
  } else {
    stop("Input must be either data frame, matrix, or phyloseq object!")
  }
  
  # Set up output
  pvals <- perfDisc_DA <- perfDisc_DV <- matrix(NA, nrow = length(taxanames),
                                                ncol = num_DA + num_DV)
  #model_summaries <- rep(list(NA), length(taxanames))
  # check to make sure inits is of the same length
  if (!is.null(inits)) {
    ncol1 <- ncol(stats::model.matrix(object = formula, data = data.frame(sample_data(data))))
    ncol2 <- ncol(stats::model.matrix(object = phi.formula, data = data.frame(sample_data(data))))
    if (length(inits) != ncol1 + ncol2) {
      stop("inits must match number of regression parameters in formula and phi.formula!")
    }
  }
  
  if (is.null(try_only)) {
    try_only <- length(taxanames)
  }
  # Loop through OTU/taxa
  for (i in 1:try_only) {
    
    # Subset data to only select that taxa
    data_i <- convert_phylo(data, select = taxanames[i])
    
    if (sum(data_i$W) == 0) {
      perfDisc_DA[i] <- TRUE
      perfDisc_DV[i] <- TRUE
    } else {
      # Update formula to match
      formula_i <- stats::update(formula, cbind(W, M - W) ~ .)
      
      # Fit unrestricted model
      mod <- suppressWarnings(try(bbdml(formula = formula_i, phi.formula = phi.formula,
                                        data = data_i, link = link, phi.link = phi.link,
                                        inits = inits, ...), silent = TRUE))
      
      if (class(mod) != "try-error") {
        # If both models fit, otherwise keep as NA
        #model_summaries[[i]] <- suppressWarnings(summary(mod))
        
        if (num_DA > 0) {
          for (contr in 1:num_DA) {
            tmp <- try(waldchisq(mod = mod, contrasts_DA = contrasts_DA[[contr]],
                                 contrasts_DV = NULL), silent = TRUE)
            if (class(tmp) != "try-error") {
              pvals[i, contr] <- tmp
            }
            perfDisc_DA[i, contr] <- mod$sep_da
            perfDisc_DV[i, contr] <- mod$sep_dv
          }
        }
        
        if (num_DV > 0) {
          for (contr in 1:num_DV) {
            tmp <- try(waldchisq(mod = mod, contrasts_DA = NULL,
                                 contrasts_DV = contrasts_DV[[contr]]), silent = TRUE)
            if (class(tmp) != "try-error") {
              pvals[i, contr + num_DA] <- tmp
            }
            perfDisc_DA[i, contr + num_DA] <- mod$sep_da
            perfDisc_DV[i, contr + num_DA] <- mod$sep_dv
          }
        }
        
      }
    }
  }
  
  disc_vec_da <- disc_vec_dv <- signif_vec <- rep(list(NA), num_DA + num_DV)
  post_fdr <- matrix(NA, nrow = length(taxanames), ncol = num_DA + num_DV)
  
  for (contr in 1:(num_DA + num_DV)) {
    ind_disc_da <- which(perfDisc_DA[, contr] == TRUE)
    ind_disc_dv <- which(perfDisc_DV[, contr] == TRUE)
    disc_vec_da[[contr]] <- taxanames[ind_disc_da]
    disc_vec_dv[[contr]] <- taxanames[ind_disc_dv]
    
    ind_disc <- union(ind_disc_da, ind_disc_dv)
    
    if (filter_discriminant && length(ind_disc) > 0) {
      # Want to keep same length, rest will ignore NAs
      pvals[ind_disc, contr] <- NA
    }
    
    if (all(is.na(pvals[, contr]))) {
      next
    }
    post_fdr[, contr] <- stats::p.adjust(pvals[, contr], method = fdr)
    names(pvals[, contr]) <- names(post_fdr[, contr]) <- taxanames
    # Record significant taxa
    signif_vec[[contr]] <- taxanames[which(post_fdr[, contr] < fdr_cutoff)]
    #signif_models <- model_summaries[which(post_fdr < fdr_cutoff)]
  }
  
  
  
  
  structure(
    list("p" = pvals, "p_fdr" = post_fdr,
         "significant_taxa" = signif_vec,
         #"significant_models" = signif_models,
         #"all_models" =  model_summaries,
         "contrasts_DA" = contrasts_DA,
         "contrasts_DV" = contrasts_DV,
         "discriminant_taxa_DA" = disc_vec_da,
         "discriminant_taxa_DV" = disc_vec_dv,
         "data" = data),
    class = "contrastsTest"
  )
}








#' Wald-type t test
#'
#' @param mod an object of class \code{bbdml}
#'
#' @return Matrix with wald test statistics and p-values. Only performs univariate tests.
#'
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' soil <- soil_phylo %>%
#' phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
#' tax_glom("Phylum")
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil)
#' waldt(mod1)
#' }
#' @export
waldt <- function(mod) {
  # Covariance matrix
  covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
  if ("try-error" %in% class(covMat)) {
    warning("Singular Hessian! Cannot calculate p-values in this setting.", immediate. = TRUE)
    np <- length(mod$param)
    se <- tvalue <- pvalue <- rep(NA, np)
  } else {
    # Standard errors
    se <- sqrt(diag(covMat))
    # test statistic
    tvalue <- mod$param/se
    # P-value
    pvalue <- 2*stats::pt(-abs(tvalue), mod$df.residual)
  }
  # make table
  coef.table <- cbind(mod$param, se, tvalue, pvalue)
  dimnames(coef.table) <- list(names(mod$param),
                               c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  return(coef.table)
}





#' Wald-type chi-squared test statistic
#'
#' This is a helper function and not intended for users
#'
#' @param mod an object of class \code{bbdml}
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#' @param contrasts_DA List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{formula}. Note that this is only available with \code{"Wald"} value for \code{test}.
#' @param contrasts_DV List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{phi.formula}. Note that this is only available with \code{"Wald"} value for \code{test}.
#'
#' @return Test statistic for Wald test.
#'
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' soil <- soil_phylo %>%
#' phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
#' tax_glom("Phylum")
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil)
#'
#' waldchisq_test(mod = mod1, restrictions = 2)
#' }
waldchisq_test <- function(mod, restrictions = NULL, restrictions.phi = NULL,
                           contrasts_DA = NULL, contrasts_DV = NULL) {
  if (length(restrictions) == 0 && length(restrictions.phi) == 0 &&
      is.null(contrasts_DA) && is.null(contrasts_DV)) {
    stop("No restrictions or contrasts provided!")
  }
  
  if (is.null(contrasts_DA) && is.null(contrasts_DV)) {
    # Covariance matrix - I_n^-1
    covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
    if ("try-error" %in% class(covMat)) {
      stop("Singular Hessian!")
    }
    
    #### Get index of terms to be tested ####
    if (!is.null(restrictions) && !is.numeric(restrictions)) {
      if (is.character(restrictions)) {
        restrictions <- getRestrictionTerms(mod = mod, restrictions = restrictions)$mu
      } else {
        stop("restrictions must be either character vector or integer vector!")
      }
    }
    if (!is.null(restrictions.phi) && !is.numeric(restrictions.phi)) {
      if (is.character(restrictions.phi)) {
        restrictions.phi <- getRestrictionTerms(mod = mod, restrictions.phi = restrictions.phi)$phi
      } else {
        stop("restrictions.phi must be either character vector or integer vector!")
      }
    }
    if (is.null(attr(restrictions.phi, "added"))) {
      restrictions.phi <- restrictions.phi + mod$np.mu
    }
    
    index <- c(restrictions, restrictions.phi)
    
    
    
    cov_test <- covMat[index, index]
    par_test <- mod$param[index]
    chi.val <- c(crossprod(par_test, chol2inv(chol(cov_test))) %*% par_test)
    attr(chi.val, "df") <- length(index)
    return(chi.val)
    #end restrictions if, begin contrasts
  } else {
    covMat <- try(chol2inv(chol(hessian(mod))), silent = TRUE)
    if ("try-error" %in% class(covMat)) {
      stop("Singular Hessian!")
    }
    if (!is.null(contrasts_DA)) {
      contr_vec <- suppressWarnings(limma::makeContrasts(contrasts = contrasts_DA, levels = colnames(mod$X.mu)))
      contr_vec <- c(contr_vec, rep(0, mod$np.phi))
    } else if (!is.null(contrasts_DV)) {
      contr_vec <- suppressWarnings(limma::makeContrasts(contrasts = contrasts_DV, levels = colnames(mod$X.phi)))
      contr_vec <- c(rep(0, mod$np.mu), contr_vec)
    }
    
    par_contr <- crossprod(contr_vec, mod$param)
    cov_contr <- c(crossprod(contr_vec, covMat)) %*% contr_vec
    chi.val <- c(crossprod(par_contr, chol2inv(chol(cov_contr))) %*% par_contr)
    attr(chi.val, "df") <- 1
    return(chi.val)
  }
}


#' Wald-type chi-squared test
#'
#' @param mod an object of class \code{bbdml}
#' @param mod_null Optional. An object of class \code{bbdml}, should be nested within \code{mod}. If not included, need to include \code{restrictions} or \code{restrictions.phi}.
#' @param restrictions Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the abundance to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the abundance.
#' @param restrictions.phi Optional. Defaults to \code{NULL}. Numeric vector indicating the parameters associated with the dispersion to test, or character vector with name of variable to test. Note that \code{1} is the intercept associated with the dispersion.
#' @param contrasts_DA List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{formula}. Note that this is only available with \code{"Wald"} value for \code{test}.
#' @param contrasts_DV List. Optional. Constructs a contrast matrix. List elements should be characters specifying contrasts in the parameters within \code{phi.formula}. Note that this is only available with \code{"Wald"} value for \code{test}.
#'
#' @return P-value from Wald test.
#'
#' @examples
#' \dontrun{
#' data(soil_phylo)
#' soil <- soil_phylo %>%
#' phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
#' phyloseq::tax_glom("Phylum")
#' mod1 <- bbdml(formula = OTU.1 ~ DayAmdmt,
#' phi.formula = ~ DayAmdmt,
#' data = soil)
#'
#' mod2 <- bbdml(formula = OTU.1 ~ 1,
#' phi.formula = ~ 1,
#' data = soil)
#'
#' # Example using mod_null
#' waldchisq(mod = mod1, mod_null = mod2)
#'
#' # Example using restrictions and restrictions.phi
#' waldchisq(mod = mod1, restrictions = 2, restrictions.phi = 2)
#' waldchisq(mod = mod1, restrictions = "DayAmdmt", restrictions.phi = "DayAmdmt")
#' waldchisq(mod = mod1, restrictions = 2, restrictions.phi = "DayAmdmt")
#' }
#' @export
waldchisq <- function(mod, mod_null = NULL, restrictions = NULL,
                      restrictions.phi = NULL,
                      contrasts_DA = NULL, contrasts_DV = NULL) {
  if (!is.null(mod_null)) {
    tmp <- getRestrictionTerms(mod = mod, mod_null = mod_null)
    restrictions <- tmp$mu
    restrictions.phi <- tmp$phi
  }
  chi.val <- try(waldchisq_test(mod, restrictions = restrictions,
                                restrictions.phi = restrictions.phi,
                                contrasts_DA = contrasts_DA,
                                contrasts_DV = contrasts_DV),
                 silent = TRUE)
  if (class(chi.val) == "try-error") {
    return(NA)
  }
  dof.dif <- attr(chi.val, "df")
  return(stats::pchisq(chi.val, dof.dif, lower.tail = FALSE))
}



