#' @useDynLib optimist
#' @importFrom Rcpp sourceCpp
NULL

#' @title plot_xts
#' @description Plots time series from xts
#' @param x xts
plot_xts <- function(x, ...){
  x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("date") %>%
    mutate(date = as.Date(date)) %>%
    gather(key, value, -date) %>%
    ggplot(aes(x = date, y = value, color = key)) +
    geom_line() +
    labs(x = "", y = "", color = "", ...)
}

#' @title get_prices_yahoo
#' @description Builds an xts containing prices downloaded from Yahoo! Finance.
#' @param yahoo_id character array containining yahoo ids.
#' @param from start date. Retrieve data no earlier than this date. (2007-01-01).
#' @param to end date. (Sys.Date())
#' @param column indicating which column should be returned: "Open", "High", "Low", "Close", "Volume", ("Adjusted")
#' @param ... additional parameters
#' @param periodicity periodicity of data to query and return. Must be one of "daily", "weekly", "monthly", ("daily")
get_prices_yahoo <- function(
  yahoo_id,
  column = 'Adjusted',
  from = "2007-01-01",
  to = Sys.Date(),
  ...,
  periodicity = "daily"
) {
  series <- list()
  i <- 0
  
  for(id in yahoo_id) {
    col <- paste(toupper(id), '.', column, sep = '')
    
    ss <- getSymbols.yahoo(id, auto.assign = FALSE,
                           from = from, to = to,
                           periodicity = periodicity, ...)[, col]
    
    i <- i + 1
    series[[i]] <- ss
  }
  
  series <- do.call(cbind, series)
  series <- na.locf(series, na.rm = FALSE, fromLast = TRUE)
  series
}

#' @title get_dlyChg_from_price
#' @description Builds an xts containing daily price changes
#' @param price xts containing prices
get_dlyChg_from_price <- function(
  price
) {
  # calculamos cambios relativos diarios
  dlyChg <- price / lag.xts(price)
  # por conveniencia asignamos 1's al primer renglon de `dlyChg`
  dlyChg[1, ] <- 1
  dlyChg
}

#' @title get_cumRet_from_dlyChg
#' @description Builds an xts containing cumulative returns
#' @param dlyChg xts containing daily changes
get_cumRet_from_dlyChg <- function(
  dlyChg
) {
  cumRet <- xts(order.by=index(dlyChg))
  for (c in 1:ncol(dlyChg)) {
    dlyChg_col <- log(dlyChg[ ,c])
    dlyChg_col[is.na(dlyChg_col)] <- 0
    cumRet <- cbind(cumRet, exp(cumsum(dlyChg_col)))
  }
  cumRet
}

#' @title get_rebalance_t2
#' @description Gets portfolio weights and value at time t2
#' @param dlyChg_t2 daily changes at t2
#' @param portWeightOpen_t2 portfolio weight at the open on t2 
#' @param portValue_t1 portfolio value at t1
get_rebalance_t2 <- function(
  dlyChg_t2,
  portWeightOpen_t2,
  portValue_t1
) {
  portSumTerm_t2 <- dlyChg_t2 * portWeightOpen_t2
  portChg_t2 <- sum(portSumTerm_t2)
  
  portWeightClose_t2 <- portSumTerm_t2 / portChg_t2
  portValue_t2 <- portValue_t1 * portChg_t2
  
  list(portWeight = portWeightClose_t2,
       portValue = portValue_t2)
}

#' @title get_rebalance_
#' @description Computes daily portofolio weights and value 
#' @param dlyChg xts that contains daily changes
#' @param rebWeight xts that contains rebalance weights
get_rebalance_ <- function(
  dlyChg,
  rebWeight
) {
  stopifnot(index(rebWeight)[1] == index(dlyChg)[1])
  
  dates <- as.character(index(dlyChg))
  reb_dates <- as.character(index(rebWeight))
  
  portWeightOpen <- list()
  portWeight <- list()
  portValue <- list()
  
  portWeight[[dates[1]]] <- as.numeric(rebWeight[dates[1]])
  portValue[[dates[1]]] <- 1.0
  
  for(i in 2:length(dates)) {
    if(dates[i-1] %in% reb_dates) {
      portWeightOpen[[dates[i]]] <- as.numeric(rebWeight[dates[i-1]])
    } else {
      portWeightOpen[[dates[i]]] <- portWeight[[dates[i-1]]]
    }
    reb <- get_rebalance_t2(as.numeric(dlyChg[dates[i]]),
                            portWeightOpen[[dates[i]]],
                            portValue[[dates[i-1]]])
    portWeight[[dates[i]]] <- reb$portWeight
    portValue[[dates[i]]] <- reb$portValue
  }
  portWeight <- do.call(rbind, portWeight)
  portValue <- do.call(rbind, portValue)
  portWeight <- xts(portWeight, order.by = index(dlyChg))
  portValue <- xts(portValue, order.by = index(dlyChg))
  names(portWeight) <- names(dlyChg)
  names(portValue) <- 'port'
  
  list(portWeight = portWeight,
       portValue = portValue)
}

#' @title get_rebalance
#' @description Computes daily portofolio weights, value and contributions 
#' @param dlyChg xts that contains daily changes
#' @param rebWeight xts that contains rebalance weights
get_rebalance <- function(
  dlyChg,
  rebWeight
) {
  rebalance <- get_rebalance_(dlyChg, rebWeight)
  
  portContrib <- list()
  for(c in 1:ncol(dlyChg)) {
    dlyChg_c <- dlyChg[, c]
    dlyChg_c$fix <- 1
    rebWeight_c <- rebWeight[, c]
    rebWeight_c$fix <- 1 - rebWeight_c[, 1]
    portContrib[[c]] <- get_rebalance_(dlyChg_c, rebWeight_c)$portValue #validate
  }
  portContrib <- do.call(cbind, portContrib)
  names(portContrib) <- names(dlyChg)
  rebalance[['portContrib']] <- portContrib
  rebalance
}

#' @title get_w_with_geomTruncDecay
#' @description Builds weights with geometric truncated decay 
#' @param T where the geometric distribution is truncated
#' @param halflife number of days that accumulate 0.5 density
#' @param interval as in uniroot
get_w_with_geomTruncDecay <- function(
  T,
  halflife,
  interval = c(0.000001, 0.999999)
) {
  median_p <- function(p) (1-(1-p)^halflife) / (1-(1-p)^T) - 0.5
  s <- uniroot(median_p, interval)
  p <- s$root
  k <- T:1
  w <- p * (1-p)^{k-1} / (1-(1-p)^T)
  w
}

#' @title get_meanRet_from_dlyChg
#' @description Builds weights with geometric truncated decay
#' @param dlyChg xts that contains daily changes
#' @param method either 'arithmetic' or 'geometric'
#' @param halflife number of days that accumulate 0.5 density
#' @param interval as in uniroot
#' @param annualization_factor in days (252 days in a year)
get_meanRet_from_dlyChg <- function(
  dlyChg,
  method = "geometric",
  halflife = NULL,
  interval = c(0.000001, 0.999999),
  annualization_factor = 1
) {
  dlyChg <- dlyChg[-1, ] #pensar como manejar lo del renglon con unos
  dlyChg <- na.trim(dlyChg)
  T <- nrow(dlyChg)
  
  if(is.null(halflife)) {
    w <- 1 / T
  } else {
    w <- get_w_with_geomTruncDecay(T = T, halflife = halflife, interval = interval)
  }
  
  if(method == "arithmetic") meanRet <- annualization_factor * apply(w * dlyChg, 2, sum, na.rm = TRUE)
  if(method == "geometric") meanRet <- exp(annualization_factor * apply(w * log(dlyChg), 2, sum, na.rm = TRUE))
  
  meanRet
}


#' @title get_rollChg_from_dlyChg
#' @description computes rolling returns
#' @param dlyChg xts that contains daily changes
#' @param roll number of days in a rolling window
#' @param method either 'arithmetic' or 'geometric'
#' @param halflife number of days that accumulate 0.5 density
#' @param trim should the first rolling window be eliminated?
#' @param interval as in uniroot
get_rollChg_from_dlyChg <- function(
  dlyChg,
  method = "geometric",
  roll = 30,
  halflife = NULL,
  trim = FALSE,
  interval = c(0.000001, 0.999999)
) {
  dlyChg <- dlyChg[-1, ] #pensar como manejar lo del renglon con unos
  
  if(is.null(halflife)) {
    w <- rep(1/roll, roll)
  } else {
    w <- get_w_with_geomTruncDecay(T = roll, halflife = halflife, interval = interval)
  }
  
  if(method == "arithmetic") rollChg <- rollapply(zoo(dlyChg),
                                                  width = roll,
                                                  align = 'right',
                                                  FUN = weighted.mean,
                                                  w = w,
                                                  na.rm = TRUE,
                                                  fill = NA)
  
  if(method == "geometric") rollChg <- exp(rollapply(zoo(log(dlyChg)),
                                                     width = roll,
                                                     align = 'right',
                                                     FUN = weighted.mean,
                                                     w = w,
                                                     na.rm = TRUE,
                                                     fill = NA))
  
  rollChg <- reclass(rollChg, dlyChg)
  na.trim(rollChg)
}



#' @title weighted.sd
#' @description Compute a weighted standard deviation.
#' @param x a vector or dataframe containing the values
#' whose weighted standard deviation is to be computed.
#' @param w a vector of weights.
#' @param na.rm remove NA values before processing,
#' the default value is FALSE.
#' @details See radiant.data::weighted.sd and
#' stats::weighted.mean.default.
#' @export
weighted.sd <- function (
  x,
  w,
  na.rm = TRUE
){
  if (na.rm) {
    i <- !is.na(x)
    w <- w[i]
    x <- x[i]
  }
  w <- w / sum(w)
  wm <- weighted.mean(x, w)
  n <- length(x)
  sqrt(sum(w * (x - wm)^2) * n / (n - 1))
}

#' @title get_sdRet_from_dlyChg
#' @description computes standard deviation of returns
#' @param dlyChg xts that contains daily changes
#' @param method either 'arithmetic' or 'geometric'
#' @param halflife number of days that accumulate 0.5 density
#' @param interval as in uniroot
#' @param annualization_factor in days (252 days in a year)

get_sdRet_from_dlyChg <- function(
  dlyChg,
  method = "arithmetic",
  halflife = NULL,
  interval = c(0.000001, 0.999999),
  annualization_factor = 1
) {
  dlyChg <- dlyChg[-1, ] #pensar como manejar lo del renglon con unos
  dlyChg <- na.trim(dlyChg)
  T <- nrow(dlyChg)
  
  if(is.null(halflife)) {
    w <- rep(1 / T, T)
  } else {
    w <- get_w_with_geomTruncDecay(T = T, halflife = halflife, interval = interval)
  }
  
  if(method == "arithmetic") sdRet <- apply(sqrt(annualization_factor) * dlyChg,
                                            2,
                                            weighted.sd,
                                            w = w,
                                            na.rm = TRUE)
  if(method == "geometric") sdRet <- exp(apply(sqrt(annualization_factor) * log(dlyChg),
                                               2,
                                               weighted.sd,
                                               w = w,
                                               na.rm = TRUE))
  
  sdRet
}



#' @title plot_riskReward_from_dlyChg
#' @description plots risk-reward
#' @param dlyChg xts that contains daily changes
#' @param roll number of days in a rolling window
#' @param roll_halflife within a rolling window

plot_riskReward_from_dlyChg <- function(
  dlyChg,
  roll,
  roll_halflife = NULL) {
  rollChg <- get_rollChg_from_dlyChg(dlyChg = dlyChg,
                                     roll = roll,
                                     halflife = roll_halflife,
                                     trim = TRUE) ^ roll
  
  reward <- get_meanRet_from_dlyChg(rollChg)
  
  risk <- get_sdRet_from_dlyChg(rollChg) #* sqrt(nrow(dlyChg))
  
  plot_data <- data.frame(
    Instrumento =  names(risk),
    Rendimiento = 100 * (as.numeric(reward) - 1),
    Riesgo = 100 * (as.numeric(risk)),
    check.names = FALSE,
    row.names = NULL
  )
  slopes <- data.frame(
    intercept = 0,
    slope = (reward - 1)/risk
  )
  
  ggplot(data = plot_data, aes(x = Riesgo, y = Rendimiento, color = Instrumento)) +
    geom_point(size = 4, alpha = 0.4) +
    geom_abline(data = slopes, aes(intercept = intercept ,slope = slope), colour = "red", linetype = 2) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(min(0, min(plot_data$Rendimiento)), max(plot_data$Rendimiento))
    ) +
    scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(min(0, min(plot_data$Riesgo)), max(plot_data$Riesgo))
    ) +
    theme_bw() +
    theme(legend.title = element_blank())
}



#' @title get_covRet_from_dlyChg
#' @description builds weighted covariance matrix
#' @param dlyChg xts that contains daily changes
#' @param halflife number of days that accumulate 0.5 density
#' @param ... parameters passed to get_w_with_geomTruncDecay()

get_covRet_from_dlyChg <- function(
  dlyChg,
  halflife = NULL,
  ...
) {
  dlyChg <- dlyChg[-1, ] #pensar como manejar lo del renglon con unos
  dlyChg <- na.omit(dlyChg)
  T <- nrow(dlyChg)
  if(!is.null(halflife)) {
    w <- get_w_with_geomTruncDecay(T, halflife, ...)
    covRet <- cov.wt(dlyChg, w)$cov
  } else {
    covRet <- cov.wt(dlyChg)$cov
  }
  
  covRet
}


#' @title get_naiveParity_weights
#' @description gets risk parity naive weights
#' @param dlyChg xts that contains daily changes
#' @param halflife number of days that accumulate 0.5 density
get_naiveParity_weights <- function(
  dlyChg,
  halflife) {
  
  sdRet <- get_sdRet_from_dlyChg(dlyChg, 
                                 halflife = halflife)
  
  naiveParity_weights <- 1 / sdRet
  naiveParity_weights / sum(naiveParity_weights)
}


#' @title get_contribParity_weights
#' @description gets risk contribution parity weights
#' @param dlyChg xts that contains daily changes
#' @param halflife number of days that accumulate 0.5 density

get_contribParity_weights <- function(
  dlyChg, 
  halflife
) {
  T <- nrow(dlyChg)
  w <- get_w_with_geomTruncDecay(T=T, halflife = halflife)
  Sigma <- cov.wt(dlyChg, w)$cov
  
  N <- ncol(Sigma)
  x0 <- rep(1/N, N)
  
  lower <- rep(0, N)
  upper <- rep(1, N)
  
  heq <- function(x) {
    sum(x) - 1
  }
  
  heqjac <- function(x) {
    rep(1, length(x))
  }
  
  contribParity <- nloptr::slsqp(
    x0 = x0, 
    fn = fn_contribParity, 
    gr = gr_contribParity, 
    lower = lower, 
    upper = upper, 
    heq = heq, 
    heqjac = heqjac, 
    Sigma = Sigma)
  
  contribParity_weights <- contribParity$par
  names(contribParity_weights) <- names(dlyChg)
  contribParity_weights
}

fn_principalParity <- function(x, P, lam, K) {
  y <- t(P) %*% x
  risk_contribs <- lam * y^2
  tot <- 0
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      tot <- tot + (risk_contribs[i] - risk_contribs[j])^2
    }
  }
  1e6 * 2.0 * tot
}

grad_principalParity <- function(x, P, lam, K) {
  y <- t(P) %*% x
  g <- numeric(nrow(P))
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      g <- g + (lam[i] * y[i]^2 - lam[j] * y[j]^2) * (lam[i] * y[i] * P[,i] - lam[j] * y[j] * P[,j])
    }
  }
  1e6 * 8.0 * g
}

#' @title get_principalParity_weights
#' @description gets risk contribution parity weights
#' @param dlyChg xts that contains daily changes
#' @param halflife number of days that accumulate 0.5 density

get_principalParity_weights <- function(
  dlyChg,
  halflife, 
  num_factors = 3
) {
  
  N <- ncol(dlyChg)
  T <- nrow(dlyChg) 
  w <- get_w_with_geomTruncDecay(T=T, halflife = halflife)
  Sigma <- cov.wt(dlyChg, w)$cov 
  eig <- eigen(Sigma, symmetric=TRUE) #eigen vectors - loadings
  P <- eig$vectors
  lam <- eig$values
  
  x0 <- rep(1 / N, N)
  lower <- rep(0.0, N)
  upper <- rep(1.0, N)
  
  principaParity <- nloptr::slsqp(
    x0 = x0, 
    fn = fn_principalParity, 
    gr = grad_principalParity, 
    lower = lower, 
    upper = upper, 
    heq = heq, 
    heqjac = heqjac, 
    P = P, 
    lam = lam, 
    K = num_factors)
  
  principalParity_weights <- principaParity$par
  names(principalParity_weights) <- names(dlyChg)
  principalParity_weights
}


#' @title get_riskParity_dlyWeights
#' @description Builds an `xts` that contains daily risk parity weights
#' @param dlyChg xts that contains daily changes
#' @param fun_riskParity_weights function to compute risk parity weights
#' @param halflife number of days that accumulate 0.5 density
#' @param min_window number of days in the first window
#' @param roll_parity days for the output's the rolling window
#' @param ... fun_riskParity_weights() parameters
get_riskParity_dlyWeights <- function(
  dlyChg, 
  fun_riskParity_weights = get_naiveParity_weights,
  halflife, 
  min_window = ceiling(halflife * 2.3), 
  roll_parity = NULL,
  ...
) {
  dlyChg <- dlyChg[-1, ]
  riskParity_dlyWeights <- data.frame()
  
  for(T in min_window:nrow(dlyChg)){
    riskParity_weights <- fun_riskParity_weights(dlyChg[1:T], 
                                                 halflife = halflife, 
                                                 ...)
    
    riskParity_dlyWeights <- rbind(riskParity_dlyWeights, 
                                   riskParity_weights)
  }
  
  colnames(riskParity_dlyWeights) <- colnames(dlyChg)
  riskParity_dlyWeights <- xts(riskParity_dlyWeights, 
                               order.by = index(dlyChg)[min_window:nrow(dlyChg)])
  
  if(!is.null(roll_parity)) {
    riskParity_dlyWeights <- get_rollChg_from_dlyChg(riskParity_dlyWeights, 
                                                     roll = roll_parity)
    
    scale <- xts(apply(riskParity_dlyWeights, 1, sum), 
                 order.by = index(riskParity_dlyWeights))
    for (c in 1:ncol(riskParity_dlyWeights)) {
      riskParity_dlyWeights[, c] <- riskParity_dlyWeights[, c] / scale
    }
  }
  
  riskParity_dlyWeights
}


#' @title get_dlyImpliedRet_from_dlyWeights
#' @description gets implied returns from portfolio weights
#' @param dlyWeights xts that contains portfolio daily weights
#' @param dlyChg xts that contains daily changes
#' @param halflife number of days that accumulate 0.5 density
#' @param min_window number of days in the first window
#' @param roll_parity days for the output's the rolling window
get_dlyImpliedRet_from_dlyWeights <- function(
  dlyWeights, 
  dlyChg, 
  halflife, 
  min_window = ceiling(halflife * 2.3), 
  roll_parity = NULL
) {
  col_check <- (sum(!names(dlyWeights) == names(dlyChg)) == 0)
  row_check <- (sum(! index(dlyWeights) %in% 
                      index(dlyChg)[min_window:nrow(dlyChg)]) == 0)
  stopifnot(col_check & row_check)
  
  dlyImpliedRet <- data.frame()
  
  for(date in as.character(index(dlyWeights))){
    dates <- index(dlyChg[paste0("/", date)])
    
    w <- get_w_with_geomTruncDecay(T = length(dates),
                                   halflife = halflife)
    Sigma <- cov.wt(dlyChg[dates], w)$cov
    
    h <- as.numeric(dlyWeights[date])
    
    dlyImpliedRet <- rbind(dlyImpliedRet,
                           t(h)  %*% Sigma)
  }
  
  xts(dlyImpliedRet, 
      order.by = index(dlyWeights))
}


#' @title get_optim_dlyWeights_from_dlyImpliedRet
#' @description gets daily weights for optimal portolio from daily implied returns
#' @param dlyImpliedRet xts that contains portfolio daily implied returns
#' @param dlyChg xts that contains daily changes
#' @param halflife number of days that accumulate 0.5 density
#' @param min_window number of days in the first window
#' @param roll_optim days for the output's the rolling window
#' @param delta risk aversion parameter

get_optim_dlyWeights_from_dlyImpliedRet <- function(
  dlyImpliedRet, 
  dlyChg, 
  halflife, 
  min_window = ceiling(halflife * 2.3), 
  roll_optim = NULL, 
  delta = 1
) {
  col_check <- (sum(!names(dlyImpliedRet) == names(dlyChg)) == 0)
  row_check <- (sum(! index(dlyImpliedRet) %in% 
                      index(dlyChg)[min_window:nrow(dlyChg)]) == 0)
  stopifnot(col_check & row_check)
  
  fn_optim <- function(
    x, mu, Sigma, delta
  ) {
    - (t(mu) %*% x) + delta * (t(x) %*% Sigma %*% x)
  }
  
  gr_optim <- function(
    x, mu, Sigma, delta
  ) {
    - mu + 2 * delta * Sigma %*% x
  }
  
  N <- ncol(dlyChg) 
  lower <- rep(0, N)
  upper <- rep(1, N)
  
  heq <- function(x) {
    sum(x) - 1
  }
  
  optim_dlyWeights <- data.frame()
  
  for(date in as.character(index(dlyImpliedRet))){
    dates <- index(dlyChg[paste0("/", date)])
    
    w <- get_w_with_geomTruncDecay(T = length(dates),
                                   halflife = halflife)
    Sigma <- cov.wt(dlyChg[dates], w)$cov
    
    mu <- as.numeric(dlyImpliedRet[date])
    
    x0 <- rep(1/N, N)
    
    optim <- slsqp(x0 = x0, 
                   fn = fn_optim, 
                   gr = gr_optim, 
                   lower = lower, 
                   upper = upper, 
                   heq = heq, 
                   mu = mu, 
                   Sigma = Sigma, 
                   delta = delta)
    
    optim_weights <- optim$par
    
    optim_dlyWeights <- rbind(optim_dlyWeights, optim_weights)
  }
  
  colnames(optim_dlyWeights) <- colnames(dlyChg)
  optim_dlyWeights <- xts(optim_dlyWeights, 
                          order.by = index(dlyImpliedRet))
  
  if(!is.null(roll_optim)) {
    optim_dlyWeights <- get_rollChg_from_dlyChg(optim_dlyWeights, 
                                                roll = roll_optim)
    
    scale <- xts(apply(optim_dlyWeights, 1, sum), 
                 order.by = index(optim_dlyWeights))
    for (c in 1:ncol(optim_dlyWeights)) {
      optim_dlyWeights[, c] <- optim_dlyWeights[, c] / scale
    }
  }
  
  optim_dlyWeights
}