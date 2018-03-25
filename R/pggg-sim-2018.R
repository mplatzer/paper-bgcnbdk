source('R/load.R')

# simulation scenarios
configs <- data.frame(T.cal=52, T.star=52)
configs <- merge(configs, data.frame(t     = c(1.6, 5, 6, 8, 17),
                                     gamma = c(0.4, 2.5, 4, 8, 20)), all=T)
configs <- merge(configs, data.frame(r=c(0.25, 0.75)), all=T)
configs <- merge(configs, data.frame(alpha=c(5, 15)), all=T)
configs <- merge(configs, data.frame(s=c(0.25, 0.75)), all=T)
configs <- merge(configs, data.frame(beta=c(5, 15)), all=T)
configs <- merge(configs, data.frame(N=c(1000, 4000)), all=T)
fwrite(configs, 'results/sim-2018-configs.csv')

nil <- mclapply(1:nrow(configs), function(idx) {
  fn <- paste0('pggg-draws/sim-2018/simulation-', idx, '.rdata')
  if (file.exists(fn)) return()
  cat(idx, '\n', file = 'pggg-draws.log', append = T)
  # generate data
  params <- configs[idx, ]
  set.seed(ifelse(params$N==4000, 1, 2)) # use different seeds for N=1000/4000
  cbs <- pggg.GenerateData(params$N, T.cal=params$T.cal, T.star=params$T.star, params=params)$cbs
  # estimate P/GGG
  set.seed(1)
  cat(idx, '-', format(Sys.time(), '%H:%M:%S'), '\n', file = 'pggg-draws.log', append = T)
  param_init <- list(t=1, gamma=1, r=1, alpha=1, s=1, beta=1)
  mcmc <- 6000
  burnin <- 2000
  thin <- 200
  draws_1k <- pnbd.mcmc.DrawParameters(cbs, mcmc=mcmc, burnin=burnin, thin=thin, chains=4, mc.cores=1, param_init=param_init) # 1s for 1000 users x 1000 steps
  xstar_1k <- mcmc.DrawFutureTransactions(cbs, draws_1k)
  cat(idx, '-', format(Sys.time(), '%H:%M:%S'), '\n', file = 'pggg-draws.log', append = T)
  draws <- pggg.mcmc.DrawParameters(cbs, mcmc=mcmc, burnin=burnin, thin=thin, chains=4, mc.cores=1, param_init=param_init) # 30s for 1000 users x 1000 steps
  xstar <- mcmc.DrawFutureTransactions(cbs, draws)
  cat(idx, '-', format(Sys.time(), '%H:%M:%S'), '\n', file = 'pggg-draws.log', append = T)
  cbs$xstar_pggg <- apply(xstar, 2, mean)
  cbs$xstar_pnbd <- apply(xstar_1k, 2, mean)
  # estimate BG/CNBD-k
  mle.params <- list()
  mle.params$bgcnbd.params <- bgcnbd.EstimateParameters(cbs)
  cbs$xstar_bgcnbd <- bgcnbd.ConditionalExpectedTransactions(
                        params = mle.params$bgcnbd.params,
                        T.star = cbs$T.star,
                        x      = cbs$x,
                        t.x    = cbs$t.x,
                        T.cal  = cbs$T.cal)
  # save to disk
  true.params <- params
  save(cbs, true.params, mle.params, draws, draws_1k, xstar, xstar_1k, file=fn)
  return()
}, mc.cores = detectCores())

# extract cbs
fns <- list.files('pggg-draws/sim-2018/', '*.rdata')
cbs_all <- rbindlist(lapply(fns, function(fn) {
  load(fn)
  cbs$fn <- fn
  return(cbs)
}))
fst::write_fst(cbs_all, 'results/sim-2018-cbs.fst')
