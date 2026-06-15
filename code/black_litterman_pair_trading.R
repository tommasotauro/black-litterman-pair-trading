# =============================================================================
# TESI — Integration of Pairs Trading Signals into the Black-Litterman Model
# Tommaso Tauro — Sapienza Università di Roma — A.Y. 2025/2026
# FILE UNICO OPERATIVO
# =============================================================================

library(dplyr)
library(tidyr)
library(quantmod)
library(tseries)

# =============================================================================
# 1. UNIVERSO DI INVESTIMENTO
# =============================================================================

investment_universe <- data.frame(
  Ticker = c(
    "AAPL","MSFT","IBM","CSCO","INTC",
    "JNJ","PFE","MRK","UNH","AMGN",
    "WMT","PG","KO","PEP","COST",
    "GE","MMM","BA","CAT","HON",
    "NEE","DUK","SO","D","EXC",
    "JPM","BAC","WFC","C","AXP",
    "HD","MCD","DIS","NKE","LOW",
    "XOM","CVX","SLB","OXY","EOG",
    "APD","SHW","ECL","PPG","FCX"
  ),
  Sector = c(
    rep("IT",5), rep("HC",5), rep("CS",5), rep("IND",5), rep("UTL",5),
    rep("FIN",5), rep("CD",5), rep("ENE",5), rep("MAT",5)
  ),
  stringsAsFactors = FALSE
)

# =============================================================================
# 2. DOWNLOAD PREZZI (2004-2018)
# =============================================================================

from_date <- "2004-01-01"
to_date   <- "2018-12-31"
tickers   <- investment_universe$Ticker

price_list <- lapply(tickers, function(t) {
  out <- try(getSymbols(t, src="yahoo", from=from_date, to=to_date,
                        auto.assign=FALSE, warnings=FALSE), silent=TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
})
names(price_list) <- tickers

good_tickers <- tickers[sapply(tickers, function(t) {
  p <- price_list[[t]]
  if (is.null(p)) return(FALSE)
  if (sum(is.na(Ad(p))) > 0) return(FALSE)
  TRUE
})]

wide_prices <- lapply(good_tickers, function(t) {
  data.frame(Date=as.Date(index(price_list[[t]])),
             Adjusted=as.numeric(Ad(price_list[[t]])),
             Ticker=t, stringsAsFactors=FALSE)
}) %>%
  bind_rows() %>%
  pivot_wider(names_from=Ticker, values_from=Adjusted)

wide_prices <- wide_prices[complete.cases(wide_prices), ]

# =============================================================================
# 3. COPPIE INTRA-SETTORE
# =============================================================================

sector_pairs <- investment_universe %>%
  filter(Ticker %in% good_tickers) %>%
  group_by(Sector) %>%
  summarise(pairs=list(combn(Ticker, 2, simplify=FALSE)), .groups="drop") %>%
  unnest(pairs) %>%
  mutate(ticker1=sapply(pairs,`[`,1), ticker2=sapply(pairs,`[`,2)) %>%
  select(Sector, ticker1, ticker2)

# =============================================================================
# 4. PARAMETRI GLOBALI
# =============================================================================

window_size    <- 60
step_size      <- 27
alpha_adf      <- 0.05
max_half_life  <- 30
lambda_view    <- 1.0
delta_t        <- 1/252
gamma_risk     <- 3.0
oos_start_date <- as.Date("2016-01-01")
oos_end_date   <- as.Date("2018-12-31")

# Soglie s-score (Avellaneda & Lee 2010)
s_open  <- 1.3
s_close <- 0.5

price_matrix <- as.matrix(wide_prices[, good_tickers])
all_dates    <- wide_prices$Date

# =============================================================================
# 5. FUNZIONI CORE
# =============================================================================

# --- Test cointegrazione Engle-Granger ---
test_cointegration <- function(p1, p2, alpha=0.05) {
  r1 <- diff(p1)/p1[-length(p1)]
  r2 <- diff(p2)/p2[-length(p2)]
  fit <- lm(r1 ~ r2)
  resid <- residuals(fit)
  X_cumsum <- cumsum(resid)
  adf <- tryCatch(adf.test(X_cumsum, alternative="stationary"), error=function(e) NULL)
  if (is.null(adf)) return(list(cointegrated=FALSE))
  list(
    cointegrated = adf$p.value < alpha,
    adf_pvalue   = adf$p.value,
    a_hat        = as.numeric(coef(fit)[1]),
    b_hat        = as.numeric(coef(fit)[2]),
    residuals    = resid,
    X_cumsum     = X_cumsum
  )
}

# --- Stima parametri OU tramite AR(1) ---
estimate_ou <- function(X, dt=1) {
  n    <- length(X)
  fit  <- lm(X[-1] ~ X[-n])
  c0   <- coef(fit)[1]
  c1   <- coef(fit)[2]
  if (is.na(c1) || c1 <= 0 || c1 >= 1) return(NULL)
  kappa    <- -log(c1)/dt
  m_eq     <- as.numeric(c0/(1-c1))
  sigma    <- as.numeric(sd(residuals(fit))*sqrt(-2*log(c1)/(dt*(1-c1^2))))
  sigma_eq <- sigma/sqrt(2*kappa)
  half_life <- log(2)/kappa
  list(kappa=kappa, m=m_eq, sigma=sigma, sigma_eq=sigma_eq, half_life=half_life)
}

# --- S-score giornaliero OOS ---
compute_sscore <- function(p1, p2, a_hat, b_hat, m, sigma_eq) {
  r1    <- diff(p1)/p1[-length(p1)]
  r2    <- diff(p2)/p2[-length(p2)]
  e_oos <- r1 - (a_hat + b_hat*r2)
  X_oos <- cumsum(e_oos)
  list(X=X_oos, s=(X_oos - m)/sigma_eq)
}

# --- Prior di mercato ---
build_prior <- function(idx_start, idx_end) {
  ret <- apply(price_matrix[idx_start:idx_end, good_tickers], 2,
               function(p) diff(p)/p[-length(p)])
  list(m_hat=colMeans(ret, na.rm=TRUE), Sigma=cov(ret, use="complete.obs"))
}

# --- Costruzione views P, q, W ---
build_bl_views <- function(active_pairs, lambda=1.0, dt=1/252) {
  K <- length(active_pairs)
  N <- length(good_tickers)
  if (K==0) return(NULL)
  P      <- matrix(0, nrow=K, ncol=N, dimnames=list(NULL, good_tickers))
  q      <- numeric(K)
  W_diag <- numeric(K)
  for (k in seq_len(K)) {
    pair       <- active_pairs[[k]]
    P[k, pair$ticker1] <-  1.0
    P[k, pair$ticker2] <- -pair$b_hat
    q[k]      <- lambda * pair$kappa * (pair$m - pair$X_current) * dt
    W_diag[k] <- (pair$sigma_eq^2) * dt
  }
  list(P=P, q=q, W=diag(W_diag))
}

# --- Posterior Black-Litterman ---
compute_bl_posterior <- function(m_hat, Sigma, P, q, W) {
  Si  <- tryCatch(solve(Sigma), error=function(e) NULL)
  Wi  <- tryCatch(solve(W),     error=function(e) NULL)
  if (is.null(Si)||is.null(Wi)) return(NULL)
  PtWi  <- t(P) %*% Wi
  M     <- tryCatch(solve(Si + PtWi %*% P), error=function(e) NULL)
  if (is.null(M)) return(NULL)
  list(m_tilde=as.numeric(M %*% (Si %*% m_hat + PtWi %*% q)), M=M)
}

# --- Mean-Variance optimization ---
mv_optimize <- function(m_tilde, Sigma, gamma=3.0) {
  Si    <- tryCatch(solve(Sigma), error=function(e) NULL)
  if (is.null(Si)) return(rep(1/length(m_tilde), length(m_tilde)))
  w_raw <- (1/gamma) * Si %*% m_tilde
  w_sum <- sum(w_raw)
  if (abs(w_sum) < 1e-10) return(rep(1/length(m_tilde), length(m_tilde)))
  as.numeric(w_raw/w_sum)
}

# --- Segnale trading ---
get_signal <- function(s, pos) {
  if (pos== 0) { if (s < -s_open) return(+1); if (s > s_open) return(-1); return(0) }
  if (pos==+1) { if (s > -s_close) return(0); return(+1) }
  if (pos==-1) { if (s <  s_close) return(0); return(-1) }
  0
}

# =============================================================================
# 6. ROLLING WINDOW COINTEGRATION TEST
# =============================================================================

train_start_idx <- which(all_dates >= as.Date("2006-01-01"))[1]
window_starts   <- seq(train_start_idx, nrow(price_matrix)-window_size, by=step_size)
n_windows       <- length(window_starts)
coint_results   <- vector("list", n_windows)

for (w in seq_along(window_starts)) {
  idx_start  <- window_starts[w]
  idx_end    <- idx_start + window_size - 1
  prices_w   <- price_matrix[idx_start:idx_end, ]
  pairs_found <- list()
  
  for (k in seq_len(nrow(sector_pairs))) {
    t1  <- sector_pairs$ticker1[k]
    t2  <- sector_pairs$ticker2[k]
    sec <- sector_pairs$Sector[k]
    p1  <- prices_w[, t1]
    p2  <- prices_w[, t2]
    if (any(is.na(c(p1,p2)))) next
    res <- test_cointegration(p1, p2, alpha=alpha_adf)
    if (!res$cointegrated) next
    ou  <- estimate_ou(res$X_cumsum, dt=1)
    if (is.null(ou) || ou$half_life > max_half_life) next
    pairs_found[[length(pairs_found)+1]] <- c(
      list(sector=sec, ticker1=t1, ticker2=t2,
           b_hat=res$b_hat, a_hat=res$a_hat,
           adf_pvalue=res$adf_pvalue,
           residuals=list(res$residuals),
           X_cumsum=list(res$X_cumsum)),
      ou
    )
  }
  
  # Best pair per settore: massimo kappa
  best_pairs <- list()
  if (length(pairs_found) > 0) {
    sectors_present <- unique(sapply(pairs_found, `[[`, "sector"))
    for (sec in sectors_present) {
      sec_p  <- Filter(function(p) p$sector==sec, pairs_found)
      kappas <- sapply(sec_p, `[[`, "kappa")
      best_pairs[[sec]] <- sec_p[[which.max(kappas)]]
    }
  }
  
  coint_results[[w]] <- list(
    window     = w,
    date_start = all_dates[idx_start],
    date_end   = all_dates[idx_end],
    date_trade = all_dates[min(idx_end+1, nrow(price_matrix))],
    idx_start  = idx_start,
    idx_end    = idx_end,
    n_pairs    = length(pairs_found),
    best_pairs = best_pairs
  )
}

# =============================================================================
# 7. BACKTESTING OOS
# =============================================================================

spy_raw <- try(getSymbols("SPY", src="yahoo",
                          from=as.character(oos_start_date-90),
                          to=as.character(oos_end_date),
                          auto.assign=FALSE, warnings=FALSE), silent=TRUE)
spy_ret_df <- data.frame(
  Date   = as.Date(index(spy_raw)),
  Return = as.numeric(dailyReturn(Ad(spy_raw)))
)

oos_dates       <- all_dates[all_dates >= oos_start_date & all_dates <= oos_end_date]
all_sectors     <- unique(investment_universe$Sector)
position_state  <- setNames(rep(0L, length(all_sectors)), all_sectors)

oos_window_idx <- which(sapply(coint_results, function(x)
  !is.null(x$date_trade) &&
    x$date_trade >= oos_start_date &&
    x$date_trade <= oos_end_date))

port_ret_vec <- numeric(length(oos_dates))
spy_ret_vec  <- numeric(length(oos_dates))

for (d in seq_along(oos_dates)) {
  today     <- oos_dates[d]
  today_idx <- which(all_dates == today)
  
  spy_today <- spy_ret_df$Return[spy_ret_df$Date == today]
  if (length(spy_today)==0) spy_today <- 0
  spy_ret_vec[d] <- spy_today
  
  # Finestra attiva più recente
  active_w <- NA
  for (w in rev(oos_window_idx)) {
    if (coint_results[[w]]$date_trade <= today) { active_w <- w; break }
  }
  
  if (is.na(active_w) || length(coint_results[[active_w]]$best_pairs)==0 ||
      length(today_idx)==0 || today_idx <= 1) {
    port_ret_vec[d] <- spy_today
    next
  }
  
  w_data     <- coint_results[[active_w]]
  best_pairs <- w_data$best_pairs
  idx_end    <- w_data$idx_end
  
  # Calcola s-score e segnale per ogni settore
  active_pairs <- list()
  for (sec in names(best_pairs)) {
    pair <- best_pairs[[sec]]
    if (today_idx <= idx_end) next
    p1_seq <- price_matrix[idx_end:today_idx, pair$ticker1]
    p2_seq <- price_matrix[idx_end:today_idx, pair$ticker2]
    if (any(is.na(c(p1_seq,p2_seq))) || length(p1_seq)<2) next
    ss      <- compute_sscore(p1_seq, p2_seq, pair$a_hat, pair$b_hat, pair$m, pair$sigma_eq)
    s_today <- tail(ss$s, 1)
    X_today <- tail(ss$X, 1)
    if (is.na(s_today) || is.infinite(s_today)) next
    new_pos <- get_signal(s_today, position_state[sec])
    position_state[sec] <- new_pos
    if (new_pos != 0) {
      p            <- pair
      p$X_current  <- X_today
      p$position   <- new_pos
      active_pairs[[sec]] <- p
    }
  }
  
  n_active <- length(active_pairs)
  
  if (n_active == 0) {
    port_ret_vec[d] <- spy_today
    next
  }
  
  # Prior + Views + Posterior + Ottimizzazione
  prior <- build_prior(w_data$idx_start, w_data$idx_end)
  views <- build_bl_views(active_pairs, lambda=lambda_view, dt=delta_t)
  
  if (is.null(views)) { port_ret_vec[d] <- spy_today; next }
  
  bl <- compute_bl_posterior(prior$m_hat, prior$Sigma, views$P, views$q, views$W)
  
  if (is.null(bl)) { port_ret_vec[d] <- spy_today; next }
  
  w_opt <- mv_optimize(bl$m_tilde, prior$Sigma, gamma=gamma_risk)
  
  # Return realizzato
  p_today     <- price_matrix[today_idx,   good_tickers]
  p_yest      <- price_matrix[today_idx-1, good_tickers]
  asset_ret   <- (p_today - p_yest) / p_yest
  port_active <- sum(w_opt * asset_ret, na.rm=TRUE)
  
  # Blending: settori senza segnale → SPY
  w_active        <- n_active / length(all_sectors)
  w_spy           <- 1 - w_active
  port_ret_vec[d] <- w_active * port_active + w_spy * spy_today
}

# =============================================================================
# 8. TRANSACTION COSTS (sensitivity analysis - Tabella 4.7)
# =============================================================================

apply_transaction_costs <- function(returns, dates, bps_per_trade, short_cost_annual,
                                    position_changes) {
  # bps_per_trade: costo per esecuzione in basis points
  # short_cost_annual: costo annuo posizioni short (es. 0.01 = 1%)
  # position_changes: vettore logico TRUE nei giorni in cui c'è un trade
  
  daily_short_cost <- short_cost_annual / 252
  tc_per_trade     <- bps_per_trade / 10000
  
  net_returns <- returns
  for (d in seq_along(returns)) {
    if (position_changes[d]) {
      net_returns[d] <- net_returns[d] - tc_per_trade
    }
    # Costo short selling proporzionale ai giorni con posizione short
    net_returns[d] <- net_returns[d] - daily_short_cost * 0.5  # approx 50% short
  }
  net_returns
}

# Individua i giorni con cambio posizione (proxy: N_active cambia)
# Qui usiamo una proxy semplificata: ogni giorno OOS con segnale attivo
# nella pratica si traccia ogni apertura/chiusura

scenarios <- list(
  list(label="Gross (0 bps, 0%)",   bps=0,  sc=0.00),
  list(label="Base (5 bps, 1%)",    bps=5,  sc=0.01),
  list(label="Medium (10 bps, 1%)", bps=10, sc=0.01),
  list(label="High (15 bps, 2%)",   bps=15, sc=0.02)
)

# Funzione performance
compute_performance <- function(returns) {
  r        <- na.omit(returns)
  ann      <- 252
  total    <- prod(1+r) - 1
  ann_ret  <- (1+total)^(ann/length(r)) - 1
  vol      <- sd(r) * sqrt(ann)
  sharpe   <- ann_ret / vol
  cum      <- cumprod(1+r)
  max_dd   <- min(cum/cummax(cum) - 1)
  calmar   <- ann_ret / abs(max_dd)
  c(Total=round(total,4), Ann=round(ann_ret,4), Vol=round(vol,4),
    Sharpe=round(sharpe,4), MaxDD=round(max_dd,4), Calmar=round(calmar,4))
}

# Portfolio returns oggetto finale
portfolio_returns <- data.frame(
  Date        = oos_dates,
  Port_Return = port_ret_vec,
  SPY_Return  = spy_ret_vec
)
portfolio_returns$Cum_Strategy <- cumprod(1 + port_ret_vec)
portfolio_returns$Cum_SPY      <- cumprod(1 + spy_ret_vec)

# Tabella performance con transaction costs
tc_results <- lapply(scenarios, function(sc) {
  net_r <- apply_transaction_costs(
    returns          = port_ret_vec,
    dates            = oos_dates,
    bps_per_trade    = sc$bps,
    short_cost_annual= sc$sc,
    position_changes = rep(TRUE, length(port_ret_vec))  # proxy conservativa
  )
  perf <- compute_performance(net_r)
  data.frame(Scenario=sc$label, t(perf), stringsAsFactors=FALSE)
})

performance_table <- do.call(rbind, tc_results)

# SPY benchmark
spy_perf <- compute_performance(spy_ret_vec)
spy_row  <- data.frame(Scenario="SPY (benchmark)", t(spy_perf), stringsAsFactors=FALSE)
performance_table <- rbind(performance_table, spy_row)

# =============================================================================
# OUTPUT FINALE:
#   portfolio_returns  → Date, Port_Return, SPY_Return, Cum_Strategy, Cum_SPY
#   performance_table  → Tabella 4.6 + 4.7 con scenari transaction costs
#   coint_results      → tutti i risultati rolling window
# =============================================================================
# =============================================================================