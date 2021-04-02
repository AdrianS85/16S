betta <- function (chats, ses, X = NA, initial_est = NULL) 
{
  if (isTRUE(is.na(X))) {
    X <- matrix(rep(1, length(chats)), ncol = 1)
  }
  consider <- !(is.na(chats) | is.na(ses) | apply(is.na(X), 
                                                  1, sum))
  chats_effective <- chats[consider]
  ses_effective <- ses[consider]
  X_effective <- as.matrix(X[consider, ])
  n <- dim(X_effective)[1]
  p <- dim(X_effective)[2]
  likelihood <- function(input) {
    ssq_u <- input[1]
    beta <- input[2:length(input)]
    W <- diag(1/(ssq_u + ses_effective^2))
    -0.5 * (sum(log(ssq_u + ses_effective^2) + (chats_effective - 
                                                  X_effective %*% beta)^2/(ssq_u + ses_effective^2)) + 
              log(det(t(X_effective) %*% W %*% X_effective)))
  }
  if (any(is.null(initial_est))) {
    initial_est <- c(var(chats_effective), solve(t(X_effective) %*% 
                                                   X_effective) %*% t(X_effective) %*% chats_effective)
  }
  output <- try(optim(initial_est, likelihood, hessian = FALSE, 
                      control = list(fnscale = -1), lower = c(0, rep(-Inf, 
                                                                     p)), method = "L-BFGS-B"), silent = TRUE)
  i = 0
  while ("try-error" %in% class(output) & i < 200) {
    i <- i + 1
    perturb <- rnorm(n = length(initial_est), mean = c(0, 
                                                       0), sd = 0.001 * i * abs(initial_est))
    initial_est_perturbed <- pmax(c(0, rep(-Inf, p)), initial_est + 
                                    perturb)
    output <- try(optim(initial_est_perturbed, likelihood, 
                        hessian = FALSE, control = list(fnscale = -1), lower = c(0, 
                                                                                 rep(-Inf, p)), method = "L-BFGS-B"), silent = TRUE)
  }
  if ("try-error" %in% class(output)) {
    stop(paste("The starting value and 200 perturbations were not", 
               "enough to find a maximum likelihood solution.", 
               "Please try again with a new choice of `initial_est`."))
  }
  ssq_u <- output$par[1]
  beta <- output$par[2:length(output$par)]
  W <- diag(1/(ssq_u + ses_effective^2))
  vars <- 1/diag(t(X_effective) %*% W %*% X_effective)
  global <- t(beta) %*% (t(X_effective) %*% W %*% X_effective) %*% 
    beta
  Q <- sum((chats_effective - X_effective %*% beta)^2/ses_effective^2)
  R <- diag(ses_effective^2)
  G <- diag(ssq_u, n)
  mytable <- list()
  mytable$table <- cbind(Estimates = beta, `Standard Errors` = sqrt(vars), 
                         `p-values` = round(2 * (1 - pnorm(abs(beta/sqrt(vars)))), 
                                            3))
  mytable$cov <- solve(t(X_effective) %*% W %*% X_effective)
  mytable$ssq_u <- ssq_u
  mytable$homogeneity <- c(Q, 1 - pchisq(Q, n - p))
  mytable$global <- c(global, 1 - pchisq(global, p - 1))
  us <- c(ssq_u * W %*% (chats_effective - X_effective %*% 
                           beta))
  blups <- rep(NA, length(chats))
  blups[consider] <- c(X_effective %*% beta + us)
  mytable$blups <- blups
  var_matrix <- matrix(NA, nrow = (n + p), ncol = (n + p))
  var_matrix[1:p, 1:p] <- t(X_effective) %*% solve(R) %*% 
    X_effective
  var_matrix[(p + 1):(n + p), (p + 1):(n + p)] <- solve(R) + 
    MASS::ginv(G)
  var_matrix[1:p, (p + 1):(n + p)] <- t(X_effective) %*% solve(R)
  var_matrix[(p + 1):(n + p), 1:p] <- solve(R) %*% X_effective
  var_matrix_inv <- MASS::ginv(var_matrix)
  blupvars <- rep(NA, length(chats))
  blupvars[consider] <- (cbind(X_effective, diag(1, n)) %*% 
                           var_matrix_inv %*% t(cbind(X_effective, diag(1, n)))) %>% 
    diag %>% sqrt %>% c
  mytable$blupses <- blupvars
  logLhat <- -0.5 * (n * log(2 * pi) + sum(log(ssq_u + ses_effective^2) + 
                                             (chats_effective - X_effective %*% beta)^2/(ssq_u + 
                                                                                           ses_effective^2)))
  mytable$loglikelihood <- logLhat
  mytable$aic <- -2 * logLhat + 2 * (1 + p)
  mytable$aicc <- mytable$aic + (2 * (1 + p)^2 + 2 * (1 + 
                                                        p))/(n - (1 + p) - 1)
  mytable$r_squared_wls <- 1 - sum((chats_effective - X_effective %*% 
                                      beta)^2)/(sum(chats_effective^2) - n * (mean(chats_effective))^2)
  return(mytable)
}