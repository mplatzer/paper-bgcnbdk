source('R/load.r')

draw_files <- list.files('pggg-draws', full.names = T)

data_stats <- function() {
  stats <- rbindlist(lapply(draw_files, function(mcmc_file) {
    cat('load', mcmc_file, '\n')
    load(mcmc_file)
    p <- function(x) paste0(r(x)*100, '%')
    r <- function(x) round(mean(x), 2)
    out <- cbs[, .(N=.N,
                   T.cal=max(T.cal), T.star=max(T.star),
                   x0=p(x==0), xx0=p(x.star==0),
                   x=r(x), xx=r(x.star),
                   x_=r(x[x>0]), xx=r(x.star[x.star>0]),
                   x5=p(x>=5), xx5=p(x.star>=5)
                   )]
    out[, file := basename(mcmc_file)]
    out
  }))
  fwrite(stats, 'results/emp-data-stats.csv')
}


param_stats <- function() {
  stats <- rbindlist(lapply(draw_files, function(mcmc_file) {
    cat('load', mcmc_file, '\n')
    load(mcmc_file)
    p1 <- bgcnbd.EstimateParameters(cbs, k=1)
    p2 <- mbgcnbd.EstimateParameters(cbs, k=1)
    p3 <- bgcnbd.EstimateParameters(cbs)
    p4 <- mbgcnbd.EstimateParameters(cbs)
    str <- function(p) paste0(
      paste0(
        c(names(p1), 'itt', 'p', 'LL'),
        ' = ',
        round(c(p, (p[3]/p[2])*p[1], p[4]/(p[4]+p[5])),2)), collapse=' ')
    out <- cbs[, .('BG/NBD' = paste0(str(p1), 'LL=', round(bgcnbd.cbs.LL(p1, cbs))),
                   'MBG/NBD' = paste0(str(p2), 'LL=', round(mbgcnbd.cbs.LL(p2, cbs))),
                   'BG/CNBD-k' = paste0(str(p3), 'LL=', round(bgcnbd.cbs.LL(p3, cbs))),
                   'MBG/CNBD-k' = paste0(str(p4), 'LL=', round(mbgcnbd.cbs.LL(p4, cbs))))]
    out[, file := basename(mcmc_file)]
    out
  }))
  fwrite(stats, 'results/emp-param-stats.csv')
}


bias <- function(act, est) (sum(est)-sum(act))/sum(act)
mae  <- function(act, est) mean(abs(act-est))
rmse <- function(act, est) sqrt(mean((act-est)^2))
lift <- function(act, old, new) 1-mae(act, new)/mae(act, old)

bench_mle <- function(cbs) {
  # (M)BG/CNBD-k
  params.bgnbd <- bgcnbd.EstimateParameters(cbs, k=1)
  params.bgcnbd <- bgcnbd.EstimateParameters(cbs)
  params.mbgnbd <- mbgcnbd.EstimateParameters(cbs, k=1)
  params.mbgcnbd <- mbgcnbd.EstimateParameters(cbs)
  add_params <- function(dt, params, prefix='') for (p in names(params)) dt[, (paste0(prefix, p)) := params[p]]
  add_params(cbs, params.bgnbd, 'bgnbd.')
  add_params(cbs, params.bgcnbd, 'bgcnbd.')
  add_params(cbs, params.mbgnbd, 'mbgnbd.')
  add_params(cbs, params.mbgcnbd, 'mbgcnbd.')
  cbs[, ll.bgnbd := bgcnbd.cbs.LL(params.bgnbd, cbs)]
  cbs[, ll.bgcnbd := bgcnbd.cbs.LL(params.bgcnbd, cbs)]
  cbs[, ll.mbgnbd := mbgcnbd.cbs.LL(params.mbgnbd, cbs)]
  cbs[, ll.mbgcnbd := mbgcnbd.cbs.LL(params.mbgcnbd, cbs)]
  cbs[, x.est.bgnbd := bgcnbd.ConditionalExpectedTransactions(params.bgnbd, T.star, x, t.x, T.cal)]
  cbs[, x.est.bgcnbd := bgcnbd.ConditionalExpectedTransactions(params.bgcnbd, T.star, x, t.x, T.cal)]
  cbs[, x.est.mbgnbd := mbgcnbd.ConditionalExpectedTransactions(params.mbgnbd, T.star, x, t.x, T.cal)]
  cbs[, x.est.mbgcnbd := mbgcnbd.ConditionalExpectedTransactions(params.mbgcnbd, T.star, x, t.x, T.cal)]
  cbs
}

bench_mcmc <- function(cbs, draws_1k, draws) {
  # Pareto/GGG
  xstar_1k <- mcmc.DrawFutureTransactions(cbs, draws_1k)
  xstar <- mcmc.DrawFutureTransactions(cbs, draws)
  cbs[, x.est.pnbd  := apply(xstar_1k, 2, mean)]
  cbs[, x.est.pggg := apply(xstar, 2, mean)]
  cbs
}

sum_stats <- function(cbs) {
  out <- rbind(
    data.table(measure = 'MAE',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(mae(cbs$x.star, cbs$x.est.bgnbd), mae(cbs$x.star, cbs$x.est.bgcnbd), mae(cbs$x.star, cbs$x.est.mbgnbd), mae(cbs$x.star, cbs$x.est.mbgcnbd), mae(cbs$x.star, cbs$x.est.pnbd), mae(cbs$x.star, cbs$x.est.pggg)))
    ,
    data.table(measure = 'RMSE',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(rmse(cbs$x.star, cbs$x.est.bgnbd), rmse(cbs$x.star, cbs$x.est.bgcnbd), rmse(cbs$x.star, cbs$x.est.mbgnbd), rmse(cbs$x.star, cbs$x.est.mbgcnbd), rmse(cbs$x.star, cbs$x.est.pnbd), rmse(cbs$x.star, cbs$x.est.pggg)))
    ,
    data.table(measure = 'BIAS',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(bias(cbs$x.star, cbs$x.est.bgnbd), bias(cbs$x.star, cbs$x.est.bgcnbd), bias(cbs$x.star, cbs$x.est.mbgnbd), bias(cbs$x.star, cbs$x.est.mbgcnbd), bias(cbs$x.star, cbs$x.est.pnbd), bias(cbs$x.star, cbs$x.est.pggg)))
    ,
    data.table(measure = 'LL',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(cbs$ll.bgnbd[1], cbs$ll.bgcnbd[1], cbs$ll.mbgnbd[1], cbs$ll.mbgcnbd[1], NA, NA)
               )
    ,
    data.table(measure = 'r',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(cbs$bgnbd.r[1], cbs$bgcnbd.r[1], cbs$mbgnbd.r[1], cbs$mbgcnbd.r[1], NA, NA)
               )
    ,
    data.table(measure = 'alpha',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(cbs$bgnbd.alpha[1], cbs$bgcnbd.alpha[1], cbs$mbgnbd.alpha[1], cbs$mbgcnbd.alpha[1], NA, NA)
               )
    ,
    data.table(measure = 'a',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(cbs$bgnbd.a[1], cbs$bgcnbd.a[1], cbs$mbgnbd.a[1], cbs$mbgcnbd.a[1], NA, NA)
               )
    ,
    data.table(measure = 'b',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(cbs$bgnbd.b[1], cbs$bgcnbd.b[1], cbs$mbgnbd.b[1], cbs$mbgcnbd.b[1], NA, NA)
               )
    ,
    data.table(measure = 'k',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = c(cbs$bgnbd.k[1], cbs$bgcnbd.k[1], cbs$mbgnbd.k[1], cbs$mbgcnbd.k[1], NA, NA)
               )
    ,
    data.table(measure = 'N',
               model = c('BG/NBD', 'BG/CNBD-k', 'MBG/NBD', 'MBG/CNBD-k', 'P/NBD', 'P/GGG'),
               value = rep(nrow(cbs), 6)
               )
  )
}

gen_stats <- function() {
  stats <- rbindlist(lapply(draw_files, function(mcmc_file) {
    cat('load', mcmc_file, '\n')
    load(mcmc_file)
    cbs <- bench_mle(cbs)
    cbs <- bench_mcmc(cbs, draws_1k, draws)
    out <- sum_stats(cbs)
    out[, file := basename(mcmc_file)]
    out
  }))
  fwrite(stats, 'results/emp-stats.csv')
}

calc_lift <- function(stats, m='MAE', models=c('MBG/NBD', 'BG/CNBD-k', 'MBG/CNBD-k', 'P/GGG'), base_model='BG/NBD') {
  lift <- merge(stats[measure==m & model %in% models, .(file, model, value)],
                    stats[measure==m & model == base_model, .(file, base=value)], by='file')
  lift <- lift[order(file, model)][, lift:=round(1-(value/base), 3)]
  lift[, data := gsub('-.*', '', file)][, model:=ordered(model, levels=models)]
}

plot_lift <- function(dt, title) {
  qplot(x=model, y=lift, data=dt, size=I(5), main=title) +
    scale_y_continuous(labels = scales::percent, limit=c(-0.05,0.2)) +
    background_grid(major = "y", minor = "none") +
    theme(axis.text.x = element_text(angle=20, vjust=0.5, size=12), axis.title = element_blank(), axis.line.x=element_blank()) +
    geom_abline(slope=0, intercept = 0)
}

plot_lifts <- function(stats) {
  lift <- calc_lift(stats, 'MAE')
  p1 <- plot_lift(lift[data=='cdnow'], 'CDs (k=1)')
  p2 <- plot_lift(lift[data=='m18'], 'Apparel & Accessories (k=1)')
  p3 <- plot_lift(lift[data=='donations'], 'Donations (k=2)')
  p4 <- plot_lift(lift[data=='grocery'], 'Groceries (k=2)')
  p5 <- plot_lift(lift[data=='dietary'], 'Dietary Supplements (k=2)')
  p6 <- plot_lift(lift[data=='office'], 'Office Supply (k=2)')
  pg <- plot_grid(p1, p2, p3, p4, p5, p6)
  save_plot('figures/emp-mae-lift.png', pg, ncol=3, nrow=2, base_aspect_ratio=0.8)
  pg
}

show_stats <- function() {
  (bias <- dcast(stats[measure=='BIAS'], file~model, fun=function(x) round(mean(x), 3)))
  #                     file BG/CNBD-k BG/NBD MBG/CNBD-k MBG/NBD  P/GGG  P/NBD
  # 1:     cdnow-draws.rdata    -0.121 -0.121     -0.162  -0.162 -0.123 -0.068 - BG/NBD
  # 2:   dietary-draws.rdata     0.310  0.582      0.288   0.575  0.494  0.508 - MBG/CNBD-2
  # 3: donations-draws.rdata     0.067  0.196      0.030   0.190 -0.109  0.147 - MBG/CNBD-2
  # 4:   grocery-draws.rdata     0.133  0.182      0.174   0.191  0.125  0.188 - BG/CNBD-2
  # 5:       m18-draws.rdata     0.084  0.084      0.057   0.057  0.140  0.116 - MBG/CNBD-2
  # 6:    office-draws.rdata     0.008 -0.008     -0.034  -0.006 -0.070 -0.010 - MBG/NBD
  (ll <- dcast(stats[measure=='LL'], file~model, fun=function(x) round(mean(x), 3)))
  #                     file   BG/CNBD-k      BG/NBD  MBG/CNBD-k     MBG/NBD P/GGG P/NBD
  # 1:     cdnow-draws.rdata   -9582.429   -9582.429   -9582.136   -9582.136    NA    NA - MBG/NBD
  # 2:   dietary-draws.rdata   -1910.976   -1911.859   -1910.065   -1911.710    NA    NA - MBG/CNBD-2
  # 3: donations-draws.rdata -144713.168 -146900.696 -144540.942 -146899.381    NA    NA - MBG/CNBD-2
  # 4:   grocery-draws.rdata  -15054.173  -15836.566  -14978.211  -15781.939    NA    NA - MBG/CNBD-2
  # 5:       m18-draws.rdata   -5284.613   -5284.613   -5285.782   -5285.782    NA    NA - BG/NBD
  # 6:    office-draws.rdata   -2898.236   -2954.237   -2898.805   -2954.356    NA    NA - BG/CNBD-2
  (rmse <- dcast(stats[measure=='RMSE'], file~model, fun=function(x) round(mean(x), 3)))
  #                     file BG/CNBD-k BG/NBD MBG/CNBD-k MBG/NBD P/GGG P/NBD
  # 1:     cdnow-draws.rdata     1.608  1.608      1.607   1.607 1.610 1.606 - MBG/NBD    & P/NBD
  # 2:   dietary-draws.rdata     0.281  0.284      0.280   0.284 0.291 0.288 - MBG/CNBD-2
  # 3: donations-draws.rdata     0.629  0.646      0.632   0.647 0.629 0.650 - BG/CNBD-2  & P/GGG
  # 4:   grocery-draws.rdata     3.101  3.153      3.192   3.212 3.142 3.206 - BG/CNBD-2
  # 5:       m18-draws.rdata     0.863  0.863      0.861   0.861 0.854 0.850 - MBG/NBD    & P/NBD
  # 6:    office-draws.rdata     0.622  0.635      0.621   0.635 0.635 0.596 - MBG/CNBD-2 & P/NBD
}


gen_stats_grocery <- function() {
  load('../paper-data/grocery-elog.rdata') # elog
  stopifnot(exists('elog'))
  cats <- names(elog)[-c(1:4)]
  stats_cats <- (lapply(cats, function(cat) {
    cat(cat, '\n')
    cohort_start <- as.Date("2006-01-01")
    cohort_end   <- as.Date("2006-06-30")
    T.cal <- as.Date("2006-12-31")
    T.tot <- as.Date("2007-12-31")
    try({
      elog.cat <- elog[elog[[cat]]==1, list(cust, date)]
      # customers are assumed 'new' if they haven't been active in first two years i.e. 2004/2005
      # but had a transaction in first quarter of 2006
      elog.cat[, first:=min(date), by="cust"]
      elog.cat <- unique(elog.cat)
      elog.cat <- elog.cat[first >= cohort_start & first <= cohort_end]
      elog.cat <- elog.cat[date <= T.tot]
      cbs <- elog2cbs(elog.cat, per="week", T.cal=T.cal, T.tot=T.tot)
      cbs <- bench_mle(cbs)
      sum_stats(cbs)[, category := cat]
    })
  }))
  stats_cats <- rbindlist(Filter(function(x) is.data.table(x), stats_cats))
  stats_cats <- stats_cats[model %in% c('BG/NBD', 'MBG/NBD', 'BG/CNBD-k', 'MBG/CNBD-k')]
  k_cats <- sapply(cats, function(cat) try(estimateRegularity(elog[elog[[cat]]==1 & date <= T.tot, list(cust, date)])))
  stats_cats <- merge(stats_cats, data.table(category=names(k_cats), k_wheat=k_cats), by='category', all.x=T)
  fwrite(stats_cats, 'results/emp-stats-cats.csv')

  stats_cats <- fread('results/emp-stats-cats.csv')
  cats500 <- stats_cats[measure=='N' & value>=500, unique(category)]
  stats_cats <- stats_cats[category %in% cats500]
  lift <- merge(stats_cats[measure=='MAE' & model %in% c('MBG/NBD', 'BG/CNBD-k', 'MBG/CNBD-k'), .(category, model, value, k_wheat=as.numeric(k_wheat))],
                stats_cats[measure=='MAE' & model == 'BG/NBD', .(category, base=value)], by='category')
  lift <- merge(lift, stats_cats[measure=='k' & model %in% c('BG/CNBD-k', 'MBG/CNBD-k'), .(category, model, k=value)], by=c('category', 'model'), all.x=T)
  lift <- merge(lift, stats_cats[measure=='N', .(category, model, N=value)], by=c('category', 'model'), all.x=T)
  lift <- lift[order(category, model)][, lift:=round(1-(value/base), 3)]
  lift[, model:=ordered(model, levels=c('MBG/NBD', 'BG/CNBD-k', 'MBG/CNBD-k'))]

  plot(lift ~ k_wheat, data = lift[model=='MBG/CNBD-k' & N>=500, .(lift, k_wheat)], xlim=c(1,3))
  abline(lm(lift ~ k_wheat, data = lift[model=='MBG/CNBD-k' & N>=500, .(lift, k_wheat)]))
  beanplot(lift ~ k, data = lift[model=='BG/CNBD-k' & N>=500, .(lift, k)], col = "lightgray", border = "grey", overallline = "median", what=c(0,1,1,1), ylim=c(0,0.25), axes=T,
           main='BG/CNBD-k Lift', las=1, ylab='Lift', xlab='k')
  beanplot(lift ~ k, data = lift[model=='MBG/CNBD-k' & N>=500, .(lift, k)], col = "lightgray", border = "grey", overallline = "median", what=c(0,1,1,1), ylim=c(0,0.25), axes=T,
           main='MBG/CNBD-k Lift', las=1, ylab='Lift', xlab='k')
  (p <- ggplot(lift[model=='BG/CNBD-k' & N>=500, .(lift, k)], aes(x=factor(k), y=lift)) + geom_boxplot(fill='grey') +
    ggtitle('BG/CNBD-k Lift') +
    scale_y_continuous(labels = scales::percent, limit=c(-0.05,0.2)) +
    background_grid(major = "y", minor = "none") +
    theme(axis.text.x = element_text(angle=20, vjust=0.5, size=12), axis.title = element_blank(), axis.line.x=element_blank()) +
    geom_abline(slope=0, intercept = 0))
  (p <- ggplot(lift[model=='MBG/CNBD-k' & N>=500, .(lift, k)], aes(x=factor(k), y=lift)) + geom_boxplot(fill='grey') +
    #geom_jitter(width=0.25) +
    ggtitle('MBG/CNBD-k Lift') +
    scale_y_continuous(labels = scales::percent, limit=c(-0.05,0.3)) +
    background_grid(major = "y", minor = "none") +
    theme(axis.text.x = element_text(angle=0, vjust=0.5, size=12), axis.title = element_blank(), axis.line.x=element_blank()) +
    geom_abline(slope=0, intercept = 0))
  save_plot("figures/emp-category-lift.png", p, base_aspect_ratio = 1.2)
}

plotTracDonations <- function() {
  elog <- data.table(read.csv('../paper-pareto-ggg/data/donations-elog.csv', stringsAsFactors = F))
  elog[, date := as.Date(date, format='%m/%d/%Y')]
  elog[, first := min(date), by='cust']
  cbs <- elog2cbs(elog, per='week', T.cal = as.Date('2005/07/01'), T.tot = as.Date('2006/07/01'))
  params.mbgcnbd <- mbgcnbd.EstimateParameters(cbs)
  params.mbgnbd  <- mbgcnbd.EstimateParameters(cbs, k=1)
  elog[, t:=as.numeric(date)]
  elog[, t0:=as.numeric(first)]
  inc.tracking <- elog[t>t0, .N, keyby=ceiling(t/14)]$N
  T.tot  <- max(ceiling(cbs$T.cal)+cbs$T.star)
  png('figures/emp-donations-track.png', width = 800, height=450)
  op <- par(mfrow=c(1,2))
  inc <- mbgcnbd.PlotTrackingInc(params.mbgcnbd, ceiling(cbs$T.cal), T.tot, inc.tracking)
  abline(v=max(cbs$T.cal)/2)
  nil <- mbgcnbd.PlotTrackingInc(params.mbgnbd, ceiling(cbs$T.cal), T.tot, inc.tracking, ymax = max(inc) * 1.05)
  abline(v=max(cbs$T.cal)/2)
  dev.off()
}

timing <- function() {
  times <- rbindlist(lapply(draw_files, function(mcmc_file) {
    cat('load', mcmc_file, '\n')
    load(mcmc_file)
    # BG/NBD
    time.para1 <- system.time({params <- bgcnbd.EstimateParameters(cbs, k=1)})[["elapsed"]]
    time.cond1 <- system.time({est <- bgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    time.aliv1 <- system.time({est <- bgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    # MBG/NBD
    time.para2 <- system.time({params <- mbgcnbd.EstimateParameters(cbs, k=1)})[["elapsed"]]
    time.cond2 <- system.time({est <- mbgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    time.aliv2 <- system.time({est <- mbgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    # BG/CNBD-k
    time.para3 <- system.time({params <- bgcnbd.EstimateParameters(cbs)})[["elapsed"]]
    time.cond3 <- system.time({est <- bgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    time.aliv3 <- system.time({est <- bgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    # MBG/CNBD-k
    time.para4 <- system.time({params <- mbgcnbd.EstimateParameters(cbs)})[["elapsed"]]
    time.cond4 <- system.time({est <- mbgcnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    time.aliv4 <- system.time({est <- mbgcnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    # Pareto/NBD (MLE)
    time.para0 <- 0#system.time({params <- pnbd.EstimateParameters(cbs)})[["elapsed"]]
    time.cond0 <- 0#system.time({est <- pnbd.ConditionalExpectedTransactions(params, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    time.aliv0 <- 0#system.time({est <- pnbd.PAlive(params, cbs$x, cbs$t.x, cbs$T.cal)})[["elapsed"]]
    # Pareto/NBD (HB)
    # draw params: 1s for 1000 users x 1000 steps
    # draw future: 4s for 1000 users x 1000 steps
    # Pareto/GGG
    # draw params: 30s for 1000 users x 1000 steps
    # draw future: 6s for 1000 users x 1000 steps
    # out
    data.table(file = basename(mcmc_file), n = nrow(cbs),
               model = c("Pareto/NBD", "BG/NBD", "MBG/NBD", "BG/CNBD-k", "MBG/CNBD-k", "Pareto/NBD (HB)", "Pareto/GGG"),
               para = c(time.para0, time.para1, time.para2, time.para3, time.para4, nrow(cbs) * 8 * 0.001, nrow(cbs) * 8 * 0.030),
               cond = c(time.cond0, time.cond1, time.cond2, time.cond3, time.cond4, nrow(cbs) * 8 * 0.004, nrow(cbs) * 8 * 0.006),
               aliv = c(time.aliv0, time.aliv1, time.aliv2, time.aliv3, time.aliv4, 0, 0))
  }))
  # fix donations: we used 32k MCMC steps instead of 8k
  times[file=="donations-draws.rdata" & model=="Pareto/NBD (HB)", para := para * 4]
  times[file=="donations-draws.rdata" & model=="Pareto/GGG", para := para * 4]
  wide <- dcast(times, file ~ model, value.var = c("para", "cond"))
  fwrite(wide, 'results/timing.csv')
  times
}

get_fit <- function() {

  censor <- 5

  fit <- rbindlist(lapply(draw_files, function(mcmc_file) {
    load(mcmc_file)
    data <- gsub("-(.*)", "", basename(mcmc_file))

    par_b1 <- bgcnbd.EstimateParameters(cbs, k=1)
    est_b1 <- bgcnbd.ConditionalExpectedTransactions(par_b1, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
    cal_b1 <- apply(mbgcnbd.pmf(par_b1, t = cbs$T.cal, x = 0:(censor-1)), 1, mean)

    par_m1 <- mbgcnbd.EstimateParameters(cbs, k=1)
    est_m1 <- mbgcnbd.ConditionalExpectedTransactions(par_m1, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
    cal_m1 <- apply(mbgcnbd.pmf(par_m1, t = cbs$T.cal, x = 0:(censor-1)), 1, mean)

    par_bk <- bgcnbd.EstimateParameters(cbs)
    est_bk <- bgcnbd.ConditionalExpectedTransactions(par_bk, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
    cal_bk <- apply(bgcnbd.pmf(par_bk, t = cbs$T.cal, x = 0:(censor-1)), 1, mean)

    par_mk <- mbgcnbd.EstimateParameters(cbs)
    est_mk <- mbgcnbd.ConditionalExpectedTransactions(par_mk, cbs$T.star, cbs$x, cbs$t.x, cbs$T.cal)
    cal_mk <- apply(mbgcnbd.pmf(par_mk, t = cbs$T.cal, x = 0:(censor-1)), 1, mean)

    est_pg <- apply(xstar, 2, mean)
    cal_pg <- apply(mcmc.pmf(draws, t = cbs$T.cal, x = 0:(censor-1)), 1, mean)

    est <- data.table(x = cbs$x, "Actuals" = cbs$x.star, "Pareto/GGG" = est_pg,
                      "BG/NBD" = est_b1, "MBG/NBD" = est_m1, "BG/CNBD-k" = est_bk, "MBG/CNBD-k" = est_mk)
    est[, bin := pmin(x, censor)]
    val <- est[, lapply(.SD, mean), keyby=bin, .SDcols=c("Actuals", "Pareto/GGG", "BG/NBD", "MBG/NBD", "BG/CNBD-k", "MBG/CNBD-k")]
    val <- melt(as.data.table(val), id.var = "bin")

    N <- cbs[, .N, keyby=.(bin=pmin(x, censor))]$N
    cal <- cbind("Pareto/GGG" = cal_pg, "BG/NBD" = cal_b1, "MBG/NBD" = cal_m1, "BG/CNBD-k" = cal_bk, "MBG/CNBD-k" = cal_mk)
    cal <- rbind(cal, 1 - apply(cal, 2, sum))
    cal <- cbind(bin = 0:censor, "N" = N, "Actuals" = N/sum(N), cal)
    cal <- melt(as.data.table(cal), id.var = "bin")

    qp <- function(tall, fn, ylab) {
      ymax <- max(tall$value)
      ymax <- ifelse(ymax < 1, ceiling(max(tall$value) * 10) / 10, ceiling(ymax))
      p <- ggplot(tall, aes(x=bin, y=value, group = variable, linetype = variable)) +
        geom_line() +
        geom_point(shape = 1) +
        scale_x_continuous(labels = c(0:(censor-1), paste0(censor, "+")), breaks = 0:censor, name = "Calibration period transactions") +
        background_grid(major = "y", minor = "none", colour.major = "grey80") +
        guides(linetype = guide_legend(title = NULL)) +
        theme(legend.position = c(0.3, 0.88), legend.background = element_rect(fill = "white"))
      if (ymax > 1) p <- p + scale_y_continuous(breaks = 0:ymax, name = ylab, limits = c(0, ymax))
      else p <- p + scale_y_continuous(name = ylab, limits = c(0, ymax), labels = scales::percent)
      save_plot(fn, p, base_aspect_ratio = 1.5)
    }
    qp(val[variable %in% c("Actuals", "MBG/NBD", "MBG/CNBD-k")],
      paste0("figures/emp-fit-val-", data, "-1.png"),
      "Holdout period transactions")
    qp(val[variable %in% c("Actuals", "Pareto/GGG", "MBG/CNBD-k")],
       paste0("figures/emp-fit-val-", data, "-2.png"),
       "Holdout period transactions")
    qp(cal[variable %in% c("Actuals", "MBG/NBD", "MBG/CNBD-k")],
       paste0("figures/emp-fit-cal-", data, "-1.png"),
       "Calibration period share")
    qp(cal[variable %in% c("Actuals", "Pareto/GGG", "MBG/CNBD-k")],
       paste0("figures/emp-fit-cal-", data, "-2.png"),
       "Calibration period share")

    cal[, type := "calibration"]
    val[, type := "validation"]
    out <- rbind(cal, val, fill=TRUE, use.names=TRUE)
    out[, data := data]
    out
  }))
  fwrite(dcast(fit, data + bin ~ type + variable), 'results/emp-fit.csv')
  fit
}
