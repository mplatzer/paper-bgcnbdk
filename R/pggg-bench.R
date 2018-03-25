source('R/load.R')

# calculate performance across simulations
cbs_all <- as.data.table(fst::read_fst('results/sim-2018-cbs.fst'))
perf <- cbs_all[,
          .(bgcnbd_mae = mean(abs(xstar_bgcnbd - x.star)),
            bgcnbd_bias = mean(xstar_bgcnbd - x.star),
            pggg_mae = mean(abs(xstar_pggg - x.star)),
            pggg_bias = mean(xstar_pggg - x.star)),
          by = fn]

configs <- fread('results/sim-2018-configs.csv')
configs[, fn := paste0('simulation-', .I, '.rdata')]
perf <- merge(perf, configs, by = 'fn')
fwrite(perf, 'results/sim-2018-perf.csv')

# plot ratio
plot(density(perf$bgcnbd_mae / perf$pggg_mae, adjust = 0.5), main = 'MAE - BG/CNBD-k vs P/GGG')
