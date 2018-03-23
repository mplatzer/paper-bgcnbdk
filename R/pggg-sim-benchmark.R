library(data.table)
library(BTYDplus)
fns <- list.files('pggg-draws/simulation/', pattern = '*.rdata$', full.names = TRUE)

# extract simulation parameters
configs <- rbindlist(lapply(fns, function(fn) {
  load(fn)
  true.params
}))

# estimate BG/CNBD-k for each simulation
nil <- lapply(fns, function(fn) {
  cat(fn, '\n')
  load(fn)
  if ('xstar_bgcnbd' %in% names(cbs)) next()
  mle.params$bgcnbd.params <- bgcnbd.EstimateParameters(cbs)
  cbs$xstar_bgcnbd <- bgcnbd.ConditionalExpectedTransactions(
                        params = mle.params$bgcnbd.params,
                        T.star = cbs$T.star,
                        x      = cbs$x,
                        t.x    = cbs$t.x,
                        T.cal  = cbs$T.cal)
  save(cbs, draws, draws_1k, mle.params, true.params, xstar, xstar_1k,
       file = fn)
  return()
})

# add P/GGG mean estimates
nil <- lapply(fns, function(fn) {
  cat(fn, '\n')
  load(fn)
  if ('xstar_pggg' %in% names(cbs)) next()
  cbs$xstar_pggg <- apply(xstar, 2, mean)
  save(cbs, draws, draws_1k, mle.params, true.params, xstar, xstar_1k,
       file = fn)
  return()
})

# calculate performance across simulations
perf <- rbindlist(lapply(fns, function(fn) {
  load(fn)
  cbind(as.data.table(true.params),
        bgcnbd_mape = mean(abs(cbs$xstar_bgcnbd - cbs$x.star)),
        bgcnbd_bias = mean(cbs$xstar_bgcnbd - cbs$x.star),
        pggg_mape = mean(abs(cbs$xstar_pgg - cbs$x.star)),
        pggg_bias = mean(cbs$xstar_pgg - cbs$x.star))
}))
fwrite(perf, 'results/pggg_sim-benchmark.csv')

plot(density(perf$bgcnbd_mape / perf$pggg_mape, adjust = 0.5), main = 'MAPE - BG/CNBD-k vs P/GGG')
