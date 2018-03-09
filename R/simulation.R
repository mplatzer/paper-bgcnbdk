source('R/load.r')

# Run large-scale simulation study.
# Simulated customers will be stored in `results/sim-cbs.csv.gz`
sim_all <- function() {

  # define parameter scenarios to be simulated
  T.cal   <- 52
  T.stars <- c(4, 16, 52)
  repeats <- 5 # how many time each scenario should be generated
  configs <- data.frame(N=c(4000))
  configs <- merge(configs, data.frame(k=c(1,2,3,4)), all=T)
  configs <- merge(configs, data.frame(r=c(0.25, 0.5, 0.75)), all=T)
  configs <- merge(configs, data.frame(alpha=c(5, 10, 15)), all=T)
  configs$alpha <- configs$alpha / configs$k # re-scale alpha
  configs <- merge(configs, data.frame(a=c(0.5, 0.75, 1.0)), all=T)
  configs <- merge(configs, data.frame(b=c(2.5, 5, 10)), all=T)
  configs$idx_config <- 1:nrow(configs)
  configs <- merge(configs, data.frame(idx_repeat = 1:repeats), all=T)
  configs <- merge(configs, data.frame(dropout_at_zero = c(FALSE, TRUE), all=T))
  configs$idx_total <- 1:nrow(configs)
  setDT(configs)

  set.seed(1)
  cbs_sim <- rbindlist(mclapply(1:nrow(configs), function(i) {
    cat(i, 'of', nrow(configs), '\n')
    config <- as.list(configs[i,])
    params <- c(k=config$k, r=config$r, alpha=config$alpha, a=config$a, b=config$b)
    if (config$dropout_at_zero) {
      cbs <- setDT(mbgcnbd.GenerateData(n=config$N, T.cal=T.cal, T.star=T.stars, params=params)$cbs)
      est.cnbd <- mbgcnbd.EstimateParameters(cbs)
      est.nbd  <- mbgcnbd.EstimateParameters(cbs, k=1)
      for (tstar in T.stars) {
        cbs[[paste0('x.cnbd', tstar)]] <- mbgcnbd.ConditionalExpectedTransactions(est.cnbd, tstar, cbs$x, cbs$t.x, cbs$T.cal)
        cbs[[paste0('x.nbd', tstar)]]  <- mbgcnbd.ConditionalExpectedTransactions(est.nbd,  tstar, cbs$x, cbs$t.x, cbs$T.cal)
      }
    } else {
      cbs <- setDT(bgcnbd.GenerateData(n=config$N, T.cal=T.cal, T.star=T.stars, params=params)$cbs)
      est.cnbd <- bgcnbd.EstimateParameters(cbs)
      est.nbd  <- bgcnbd.EstimateParameters(cbs, k=1)
      for (tstar in T.stars) {
        cbs[[paste0('x.cnbd', tstar)]] <- bgcnbd.ConditionalExpectedTransactions(est.cnbd, tstar, cbs$x, cbs$t.x, cbs$T.cal)
        cbs[[paste0('x.nbd', tstar)]]  <- bgcnbd.ConditionalExpectedTransactions(est.nbd,  tstar, cbs$x, cbs$t.x, cbs$T.cal)
      }
    }
    est.cnbd <- as.list(est.cnbd)
    est.nbd <- as.list(est.nbd)
    names(est.cnbd) <- paste0(names(est.cnbd), '.cnbd')
    names(est.nbd) <- paste0(names(est.nbd), '.nbd')
    cbs <- cbind(cbs, data.frame(config))
    cbs <- cbind(cbs, data.frame(est.cnbd))
    cbs <- cbind(cbs, data.frame(est.nbd))
    return(cbs)
  }, mc.cores = detectCores()))
  write_csv(cbs_sim, 'results/sim-cbs.csv.gz')
  cbs_sim
}

# Calculate summary stats for all simulations.
# Statistics will be stored in `results/sim-stats.csv`
get_stats <- function(cbs) {
  # calculate summary stats for a specific T.star holdout period
  get_stats_for_tstar <- function(cbs, T.star) {
    bias <- function(act, est) (sum(est)-sum(act))/sum(act)
    mae  <- function(act, est) mean(abs(act-est))
    rmse <- function(act, est) sqrt(mean((act-est)^2))
    set(cbs, j = 'x.star', value = cbs[[paste0('x.star', T.star)]])
    set(cbs, j = 'x.cnbd', value = cbs[[paste0('x.cnbd', T.star)]])
    set(cbs, j = 'x.nbd', value = cbs[[paste0('x.nbd', T.star)]])
    stats <- cbs[, .(bias_cnbd = bias(x.star, x.cnbd),
                     bias_nbd  = bias(x.star, x.nbd),
                     mae_cnbd  = mae(x.star, x.cnbd),
                     mae_nbd   = mae(x.star, x.nbd),
                     rmse_cnbd = rmse(x.star, x.cnbd),
                     rmse_nbd  = rmse(x.star, x.nbd),
                     lift_mae  = mae(x.star, x.cnbd) / mae(x.star, x.nbd),
                     lift_rmse = rmse(x.star, x.cnbd) / rmse(x.star, x.nbd),
                     x0=mean(x==0), alive=mean(alive), p=mean(p), lambda=mean(lambda), x_=mean(x[x>0])),
                   by=.(N, k, r, alpha, a, b, dropout_at_zero, idx_total)]
    stats[, T.star := T.star]
    stats[, k := as.ordered(k)]
    stats
  }
  stats <- rbind(get_stats_for_tstar(cbs, 52),
                 get_stats_for_tstar(cbs, 16),
                 get_stats_for_tstar(cbs, 4))
  write_csv(stats, 'results/sim-stats.csv')
  stats
}

# cbs_sim <- sim_all()
# stats <- get_stats(cbs_sim)


report_basic_stats <- function() {
  stats <- read_csv('results/sim-stats.csv')
  stats_paper <- stats[T.star==52 & dropout_at_zero==F]
  round(summary(stats_paper[, (1-x0)*x_]), 2) # number of transactions during calibration
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.47    1.12    1.76    2.02    2.66    6.27
  round(summary(stats_paper[, x0]), 2)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.15    0.31    0.44    0.45    0.59    0.76
}


plot_some <- function() {
  stats <- read_csv('results/sim-stats.csv')

  ss <- stats[T.star==52 & dropout_at_zero==F]
  st <- melt(ss[, .(k, bias_nbd, bias_cnbd, mae_nbd, mae_cnbd)], id.vars='k')
  st[, measure := ifelse(grepl('bias_', variable), 'Bias', 'MAE')]
  st[, model := ordered(ifelse(grepl('_cnbd', variable), 'BG/CNBD-k', 'BG/NBD'), levels=c('BG/NBD', 'BG/CNBD-k'))]
  st[measure=='Bias', median(value), by=.(model, k)]

  p <- qplot(k, value, fill=model, geom='boxplot', data=st[measure=='Bias'],
        main = 'Aggregate Level Error - Bias',
        outlier.size = .5) +
    xlab('Regularity k') + ylab('') +
    scale_fill_manual(name = "Model", values = c("grey90", "grey60")) +
    background_grid(major = "y", minor = "none", colour.major = "grey80") +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-0.22,0.22))
  save_plot("figures/sim-bias.png", p, base_aspect_ratio = 1.5)

  p <- qplot(k, value, fill=model, geom='boxplot', data=st[measure=='MAE'],
        main = 'Individual Level Error - MAE',
        outlier.size = .5) +
    xlab('Regularity k') + ylab('') +
    scale_fill_manual(name = "Model", values = c("grey90", "grey60")) +
    background_grid(major = "y", minor = "none", colour.major = "grey80") +
    coord_cartesian(ylim = c(0.4, 2))
  save_plot("figures/sim-mae.png", p, base_aspect_ratio = 1.5)

  p <- qplot(k, 1-lift_mae, fill=I("grey80"), geom='boxplot', data=ss,
             main = 'Individual Level Error - MAE Lift',
             outlier.size = .5) +
    xlab('Regularity k') + ylab('') +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = c(-0.03, 0.2)) +
    background_grid(major = "y", minor = "none", colour.major = "grey80")
  save_plot("figures/sim-lift.png", p, base_aspect_ratio = 1.3)

  ys <- scale_y_continuous(labels = scales::percent, limit=c(-0.05,0.2))
  bg <- background_grid(major = "xy", minor = "none")
  p1 <- qplot(x0, 1-lift_mae, col=k, data=ss) + xlab('Share of Inactive') + ylab('Lift in MAE') + ys + bg
  p2 <- qplot(alive, 1-lift_mae, col=k, data=ss) + xlab('P(alive)') + ylab('Lift in MAE') + ys + bg
  p3 <- qplot(x_, 1-lift_mae, col=k, data=ss) + xlab('Nr. of Transactions of Actives') + ylab('Lift in MAE') + ys + bg
  p4 <- qplot(p, 1-lift_mae, col=k, data=ss) + xlab('dropout p') + ylab('Lift in MAE') + ys + bg
  p <- plot_grid(p1, p2, p3, p4)
  save_plot("figures/sim-scatterplot.png", p, ncol = 2, nrow = 1, base_aspect_ratio = 1)

  p <- qplot(lambda, p, col=(lift_mae>1), data=ss[k>=2]) + bg
  save_plot("figures/sim-badcases.png", p, ncol = 2, nrow = 1, base_aspect_ratio = 1.3)

  p <- qplot(lambda/as.integer(as.character(k)), 1-lift_mae, data=ss, color=k,
        main = 'Individual Level Error - MAE Lift') +
    xlab(expression(paste('Purchase Frequency ', (lambda/k)))) + ylab('') +
    coord_cartesian(ylim = c(-0.05, 0.2)) +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_grey(start = 0.6, end = 0.2) +
    background_grid(major = "y", minor = "none", colour.major = "grey80")
  save_plot("figures/sim-lift-vs-freq.png", p, base_aspect_ratio = 1.5)

  p <- qplot(p, 1-lift_mae, data=ss, color=k,
             main = 'Individual Level Error - MAE Lift') +
    xlab('Dropout Probability (p)') + ylab('') +
    coord_cartesian(ylim = c(-0.05, 0.2)) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    scale_colour_grey(start = 0.6, end = 0.2) +
    background_grid(major = "y", minor = "none", colour.major = "grey80")
  save_plot("figures/sim-lift-vs-drop.png", p, base_aspect_ratio = 1.5)

  p <- qplot(alive, 1-lift_mae, data=ss, color=k,
             main = 'Individual Level Error - MAE Lift') +
    xlab('Share of Customers still alive at T') + ylab('') +
    coord_cartesian(ylim = c(-0.05, 0.2)) +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(labels = scales::percent) +
    scale_colour_grey(start = 0.6, end = 0.2) +
    background_grid(major = "y", minor = "none", colour.major = "grey80")
  save_plot("figures/sim-lift-vs-palive.png", p, base_aspect_ratio = 1.5)
}


# Basic speed comparison between BG/CNBD-k and Pareto/GGG
# (running on single core)
clock <- function() {
  set.seed(1)
  cbs <- bgcnbd.GenerateData(4000, 52, 52, c(k=2, r=0.25, alpha=2.5, a=0.25, b=2.5))$cbs
  system.time({
    param.est <- bgcnbd.EstimateParameters(cbs)
    cbs$est1 <- bgcnbd.ConditionalExpectedTransactions(param.est, 52, cbs$x, cbs$t.x, cbs$T.cal)
  })
  # 1.6secs
  system.time({
    draws <- pggg.mcmc.DrawParameters(cbs, mc.cores=1, chains=2)
    xstar <- mcmc.DrawFutureTransactions(cbs, draws, T.star=52)
    cbs$est2 <- apply(xstar, 2, mean)
  })
  # 650secs
}


tracking_plot <- function() {
  set.seed(4)
  dt <- bgcnbd.GenerateData(4000, 52, 52, c(k=2, r=0.5, alpha=10/2, a=0.75, b=2.5))
  cbs <- dt$cbs
  est <- bgcnbd.EstimateParameters(cbs)
  est1 <- bgnbd.EstimateParameters(cbs)
  elog <- dt$elog
  elog <- data.table::setDT(elog)
  elog <- elog[, t0 := min(t), by=cust]
  inc.tracking <- elog[t>t0, .N, keyby=ceiling(t)]$N
  T.cal <- max(cbs$T.cal)
  T.tot  <- max(cbs$T.cal+cbs$T.star)
  op <- par(mfrow=c(1,2))
  inc <- bgcnbd.PlotTrackingInc(est, cbs$T.cal, T.tot, inc.tracking)
  nil <- bgcnbd.PlotTrackingInc(c(1, est1), cbs$T.cal, T.tot, inc.tracking, ymax = max(inc) * 1.05)
  par(op)
  dev.print(png, 'figures/sim-trend.png', height=500, width=1000)
  cat('mean p true', 0.75/(0.75+2.5))
  # 23%
  cat('BG/CNBD-k estimate for mean p', est[['a']]/(est[['a']]+est[['b']]))
  # 23.8%
  cat('BG/NBD estimate for mean p', est1[3]/(est1[3]+est1[4]))
  # 14%
}
