source('R/load.R')

# calculate performance across simulations
cbs_all <- as.data.table(fst::read_fst('results/sim-2018-cbs.fst'))
perf <- cbs_all[,
          .(bgcnbd_mae = mean(abs(xstar_bgcnbd - x.star)),
            bgcnbd_bias = mean(xstar_bgcnbd - x.star),
            pggg_mae = mean(abs(xstar_pggg - x.star)),
            pggg_bias = mean(xstar_pggg - x.star)),
          by = fn]

# add configs & timings
configs <- fread('results/sim-2018-configs.csv')
configs[, fn := paste0('simulation-', .I, '.rdata')]
lines <- readLines('results/pggg-draws.log')
timings <- sapply(1:160, function(idx) {
  times <- grep(paste0('^', idx, ' - '), lines, value = T)
  start <- gsub(paste0('^', idx, ' - '), '', head(times, 1))
  end   <- gsub(paste0('^', idx, ' - '), '', tail(times, 1))
  dur <- as.numeric(hms(end)) - as.numeric(hms(start))
  if (dur<0) dur <- dur + 24 * 3600
  if (dur>12*3600) dur <- 24 * 3600 - dur
  return(dur)
})
configs$compute_time <- timings

# store to disk
perf <- merge(perf, configs, by = 'fn')
perf[, idx := as.integer(gsub('simulation-([0-9]+)\\.rdata', '\\1', fn))]
setkey(perf, idx)
fwrite(perf, 'results/sim-2018-perf.csv')
perf[, (mean(compute_time)) / 3600, by = N]

# plot ratio
plot(density(perf$bgcnbd_mae / perf$pggg_mae, adjust = 0.5), main = 'MAE - BG/CNBD-k vs P/GGG')
