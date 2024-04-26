# Generation of simulated datasets
simulate_hashtag_data <- function(ncells,
                                  nhashtags,
                                  multiplet_pct,
                                  quality){
  
  sim_neg_distrib <- function(ncells){
    return(rnbinom(n = ncells, mu = 20, size = 10))
  }
  
  sim_pos_distrib <- function(ncells, quality){
    if(quality == "high"){
      return(rnbinom(n = ncells, mu = 100, size = 10))
    } else if (quality == "medium"){
      return(rnbinom(n = ncells, mu = 80, size = 10))
    } else if (quality == "low"){
      return(rnbinom(n = ncells, mu = 60, size = 10))
    } else {
      stop("Unknown \"quality\" input! Should be either \"high\", \"medium\", or \"low\"")
    }
  }
  
  cells_per_hash <- rep(ncells / nhashtags, nhashtags)
  names(cells_per_hash) <- 1:nhashtags
  n_multiplets <- ceiling(ncells * multiplet_pct / 100)
  simulate_multiplets <- sapply(X = 1:n_multiplets, function(x){sample(1:nhashtags, 2, replace = T)})
  multiplets_table <- table(simulate_multiplets)
  final_cells_per_hash <- cells_per_hash - floor(as.vector(multiplets_table) / 2)
  
  counts <- lapply(X = 1:nhashtags, function(x){
    ct <- matrix(nrow = final_cells_per_hash[x], ncol = nhashtags)
    rownames(ct) <- paste0(x, "-", 1:final_cells_per_hash[x])
    colnames(ct) <- paste0("HTO-", 1:nhashtags)
    for(j in 1:ncol(ct)){
      if(x == j){
        ct[, j] <- sim_pos_distrib(final_cells_per_hash[x], quality)
      } else {
        ct[, j] <- sim_neg_distrib(final_cells_per_hash[x])
      }
    }
    return(ct)
  })
  counts <- do.call(rbind, counts)
  
  multiplet_counts <- t(sapply(X = 1:ncol(simulate_multiplets), function(i){
    cell1 <- sim_neg_distrib(nhashtags)
    cell2 <- sim_neg_distrib(nhashtags)
    
    cell1[simulate_multiplets[1, i]] <- sim_pos_distrib(1, quality)
    cell2[simulate_multiplets[2, i]] <- sim_pos_distrib(1, quality)
    return(cell1 + cell2)
  }))
  multiplet_names <- sapply(X = 1:ncol(simulate_multiplets), function(x){
    paste0("mp", x, ".", simulate_multiplets[1, x], "-", simulate_multiplets[2, x])
  })
  rownames(multiplet_counts) <- multiplet_names
  colnames(multiplet_counts) <- paste0("HTO-", 1:nhashtags)
  
  final_counts <- rbind(counts, multiplet_counts)
  return(final_counts)
}

low_4k <- simulate_hashtag_data(ncells = 4000, nhashtags = 8, multiplet_pct = 0.43, quality = "low")
low_8k <- simulate_hashtag_data(ncells = 8000, nhashtags = 8, multiplet_pct = 0.91, quality = "low")
low_12k <- simulate_hashtag_data(ncells = 12000, nhashtags = 8, multiplet_pct = 1.47, quality = "low")

med_4k <- simulate_hashtag_data(ncells = 4000, nhashtags = 8, multiplet_pct = 0.43, quality = "medium")
med_8k <- simulate_hashtag_data(ncells = 8000, nhashtags = 8, multiplet_pct = 0.91, quality = "medium")
med_12k <- simulate_hashtag_data(ncells = 12000, nhashtags = 8, multiplet_pct = 1.47, quality = "medium")

high_4k <- simulate_hashtag_data(ncells = 4000, nhashtags = 8, multiplet_pct = 0.43, quality = "high")
high_8k <- simulate_hashtag_data(ncells = 8000, nhashtags = 8, multiplet_pct = 0.91, quality = "high")
high_12k <- simulate_hashtag_data(ncells = 12000, nhashtags = 8, multiplet_pct = 1.47, quality = "high")

# Implementation of BFF
library(cellhashR)

run_demux <- function(data, 
                      metricsFile,
                      methods,
                      doTSNE = F,
                      bff_cluster.doublet_thresh = 0.05,
                      bff_cluster.neg_thresh = 0.05,
                      bff_cluster.dist_frac = 0.1){
  cellbarcodeWhitelist <- colnames(data)
  output <- GenerateCellHashingCalls(barcodeMatrix = data, methods = methods, 
                                     cellbarcodeWhitelist = cellbarcodeWhitelist, 
                                     metricsFile = metricsFile, doTSNE = doTSNE,
                                     bff_cluster.doublet_thresh=bff_cluster.doublet_thresh, 
                                     bff_cluster.neg_thresh=bff_cluster.neg_thresh,
                                     bff_cluster.dist_frac=bff_cluster.dist_frac,
                                     htodemux.positive.quantile = 0.99)
  return(output)
}
dir.create("bff_raw")
bff_raw_high_12k <- run_demux(t(high_12k), "bff_raw/high_12k_metrics.txt", methods = "bff_raw")
bff_raw_high_8k <- run_demux(t(high_8k), "bff_raw/high_8k_metrics.txt", methods = "bff_raw")
bff_raw_high_4k <- run_demux(t(high_4k), "bff_raw/high_4k_metrics.txt", methods = "bff_raw")
bff_raw_med_12k <- run_demux(t(med_12k), "bff_raw/med_12k_metrics.txt", methods = "bff_raw")
bff_raw_med_8k <- run_demux(t(med_8k), "bff_raw/med_8k_metrics.txt", methods = "bff_raw")
bff_raw_med_4k <- run_demux(t(med_4k), "bff_raw/med_4k_metrics.txt", methods = "bff_raw")
bff_raw_low_12k <- run_demux(t(low_12k), "bff_raw/low_12k_metrics.txt", methods = "bff_raw")
bff_raw_low_8k <- run_demux(t(low_8k), "bff_raw/low_8k_metrics.txt", methods = "bff_raw")
bff_raw_low_4k <- run_demux(t(low_4k), "bff_raw/low_4k_metrics.txt", methods = "bff_raw")

dir.create("bff_cluster")
bff_cluster_high_12k <- run_demux(t(high_12k), "bff_cluster/high_12k_metrics.txt", methods = "bff_cluster")
bff_cluster_high_8k <- run_demux(t(high_8k), "bff_cluster/high_8k_metrics.txt", methods = "bff_cluster")
bff_cluster_high_4k <- run_demux(t(high_4k), "bff_cluster/high_4k_metrics.txt", methods = "bff_cluster")
bff_cluster_med_12k <- run_demux(t(med_12k), "bff_cluster/med_12k_metrics.txt", methods = "bff_cluster")
bff_cluster_med_8k <- run_demux(t(med_8k), "bff_cluster/med_8k_metrics.txt", methods = "bff_cluster")
bff_cluster_med_4k <- run_demux(t(med_4k), "bff_cluster/med_4k_metrics.txt", methods = "bff_cluster")
bff_cluster_low_12k <- run_demux(t(low_12k), "bff_cluster/low_12k_metrics.txt", methods = "bff_cluster")
bff_cluster_low_8k <- run_demux(t(low_8k), "bff_cluster/low_8k_metrics.txt", methods = "bff_cluster")
bff_cluster_low_4k <- run_demux(t(low_4k), "bff_cluster/low_4k_metrics.txt", methods = "bff_cluster")

# Implementation of HashedDrops
library(DropletUtils)
hashed_drops_high_12k <- hashedDrops(t(high_12k), confident.min = 0.5, doublet.min = 0.5)
hashed_drops_high_8k <- hashedDrops(t(high_8k), confident.min = 0.5, doublet.min = 0.5)
hashed_drops_high_4k <- hashedDrops(t(high_4k), confident.min = 0.5, doublet.min = 0.5)
hashed_drops_med_12k <- hashedDrops(t(med_12k), confident.min = 0.5, doublet.min = 0.5)
hashed_drops_med_8k <- hashedDrops(t(med_8k), confident.min = 0.5, doublet.min = 0.5)
hashed_drops_med_4k <- hashedDrops(t(med_4k), confident.min = 0.5, doublet.min = 0.5)
hashed_drops_low_12k <- hashedDrops(t(low_12k), confident.min = 0.5, doublet.min = 0.5)
hashed_drops_low_8k <- hashedDrops(t(low_8k), confident.min = 0.5, doublet.min = 0.5)
hashed_drops_low_4k <- hashedDrops(t(low_4k), confident.min = 0.5, doublet.min = 0.5)

# Implementation of GMM-Demux
dir.create("gmm_demux")
gmm_demux_high_12k <- run_demux(t(high_12k), "gmm_demux/high_12k_metrics.txt", methods = "gmm_demux")
gmm_demux_high_8k <- run_demux(t(high_8k), "gmm_demux/high_8k_metrics.txt", methods = "gmm_demux")
gmm_demux_high_4k <- run_demux(t(high_4k), "gmm_demux/high_4k_metrics.txt", methods = "gmm_demux")
gmm_demux_med_12k <- run_demux(t(med_12k), "gmm_demux/med_12k_metrics.txt", methods = "gmm_demux")
gmm_demux_med_8k <- run_demux(t(med_8k), "gmm_demux/med_8k_metrics.txt", methods = "gmm_demux")
gmm_demux_med_4k <- run_demux(t(med_4k), "gmm_demux/med_4k_metrics.txt", methods = "gmm_demux")
gmm_demux_low_12k <- run_demux(t(low_12k), "gmm_demux/low_12k_metrics.txt", methods = "gmm_demux")
gmm_demux_low_8k <- run_demux(t(low_8k), "gmm_demux/low_8k_metrics.txt", methods = "gmm_demux")
gmm_demux_low_4k <- run_demux(t(low_4k), "gmm_demux/low_4k_metrics.txt", methods = "gmm_demux")

# Implementation of HTO demux
dir.create("htodemux")
htodemux_high_12k <- run_demux(t(high_12k), "htodemux/high_12k_metrics.txt", methods = "htodemux")
htodemux_high_8k <- run_demux(t(high_8k), "htodemux/high_8k_metrics.txt", methods = "htodemux")
htodemux_high_4k <- run_demux(t(high_4k), "htodemux/high_4k_metrics.txt", methods = "htodemux")
htodemux_med_12k <- run_demux(t(med_12k), "htodemux/med_12k_metrics.txt", methods = "htodemux")
htodemux_med_8k <- run_demux(t(med_8k), "htodemux/med_8k_metrics.txt", methods = "htodemux")
htodemux_med_4k <- run_demux(t(med_4k), "htodemux/med_4k_metrics.txt", methods = "htodemux")
htodemux_low_12k <- run_demux(t(low_12k), "htodemux/low_12k_metrics.txt", methods = "htodemux")
htodemux_low_8k <- run_demux(t(low_8k), "htodemux/low_8k_metrics.txt", methods = "htodemux")
htodemux_low_4k <- run_demux(t(low_4k), "htodemux/low_4k_metrics.txt", methods = "htodemux")

# Implementation of demuxmix
library(demuxmix)
demuxmix_high_12k <- demuxmix(hto = t(high_12k), model = "naive")
demuxmix_high_8k <- demuxmix(hto = t(high_8k), model = "naive")
demuxmix_high_4k <- demuxmix(hto = t(high_4k), model = "naive")
demuxmix_med_12k <- demuxmix(hto = t(med_12k), model = "naive")
demuxmix_med_8k <- demuxmix(hto = t(med_8k), model = "naive")
demuxmix_med_4k <- demuxmix(hto = t(med_4k), model = "naive")
demuxmix_low_12k <- demuxmix(hto = t(low_12k), model = "naive")
demuxmix_low_8k <- demuxmix(hto = t(low_8k), model = "naive")
demuxmix_low_4k <- demuxmix(hto = t(low_4k), model = "naive")

# Classification plots
library(ggplot2)
library(viridisLite)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
bff_raw_assignments <- rbind(prop.table(table(bff_raw_high_12k$consensuscall.global)), 
                             prop.table(table(bff_raw_high_8k$consensuscall.global)),
                             prop.table(table(bff_raw_high_4k$consensuscall.global)),
                             prop.table(table(bff_raw_med_12k$consensuscall.global)), 
                             prop.table(table(bff_raw_med_8k$consensuscall.global)),
                             prop.table(table(bff_raw_med_4k$consensuscall.global)),
                             c(0, 1, 0), 
                             c(0, 1, 0), 
                             c(0, 1, 0))
bff_raw_assignments <- bff_raw_assignments[, c(3, 1, 2)]
rownames(bff_raw_assignments) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                   "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(bff_raw_assignments), las = 2, col = rev(viridis(3)), main = "BFF Raw")
legend(x = "topright", inset=c(-0.3,0), legend = c("Singlet", "Doublet", "Negative"), 
       col = rev(viridis(3)), pch = 15)

bff_cluster_assignments <- rbind(prop.table(table(bff_cluster_high_12k$consensuscall.global)), 
                                 prop.table(table(bff_cluster_high_8k$consensuscall.global)),
                                 prop.table(table(bff_cluster_high_4k$consensuscall.global)),
                                 prop.table(table(bff_cluster_med_12k$consensuscall.global)), 
                                 prop.table(table(bff_cluster_med_8k$consensuscall.global)),
                                 prop.table(table(bff_cluster_med_4k$consensuscall.global)),
                                 c(0, 1, 0), 
                                 c(0, 1, 0), 
                                 c(0, 1, 0))
bff_cluster_assignments <- bff_cluster_assignments[, c(3, 1, 2)]
rownames(bff_cluster_assignments) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                       "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(bff_cluster_assignments), las = 2, col = rev(viridis(3)), main = "BFF Cluster")
legend(x = "topright", inset=c(-0.3,0), legend = c("Singlet", "Doublet", "Negative"), 
       col = rev(viridis(3)), pch = 15)


gmm_demux_assignments <- rbind(prop.table(table(gmm_demux_high_12k$consensuscall.global)), 
                               prop.table(table(gmm_demux_high_8k$consensuscall.global)),
                               prop.table(table(gmm_demux_high_4k$consensuscall.global)),
                               prop.table(table(gmm_demux_med_12k$consensuscall.global)), 
                               prop.table(table(gmm_demux_med_8k$consensuscall.global)),
                               prop.table(table(gmm_demux_med_4k$consensuscall.global)),
                               prop.table(table(gmm_demux_low_12k$consensuscall.global)), 
                               prop.table(table(gmm_demux_low_8k$consensuscall.global)),
                               prop.table(table(gmm_demux_low_4k$consensuscall.global)))
gmm_demux_assignments <- gmm_demux_assignments[, c(3, 1, 2)]
rownames(gmm_demux_assignments) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                     "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(gmm_demux_assignments), las = 2, col = rev(viridis(3)), main = "GMM Demux")
legend(x = "topright", inset=c(-0.3,0), legend = c("Singlet", "Doublet", "Negative"), 
       col = rev(viridis(3)), pch = 15)

htodemux_assignments <- rbind(prop.table(table(htodemux_high_12k$consensuscall.global)), 
                              prop.table(table(htodemux_high_8k$consensuscall.global)),
                              prop.table(table(htodemux_high_4k$consensuscall.global)),
                              prop.table(table(htodemux_med_12k$consensuscall.global)), 
                              prop.table(table(htodemux_med_8k$consensuscall.global)),
                              prop.table(table(htodemux_med_4k$consensuscall.global)),
                              prop.table(table(htodemux_low_12k$consensuscall.global)), 
                              prop.table(table(htodemux_low_8k$consensuscall.global)),
                              prop.table(table(htodemux_low_4k$consensuscall.global)))
htodemux_assignments <- htodemux_assignments[, c(3, 1, 2)]
rownames(htodemux_assignments) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                    "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(htodemux_assignments), las = 2, col = rev(viridis(3)), main = "HTODemux")
legend(x = "topright", inset=c(-0.3,0), legend = c("Singlet", "Doublet", "Negative"), 
       col = rev(viridis(3)), pch = 15)

hashed_drops_assignment <- rbind(prop.table(c(sum(hashed_drops_high_12k$Confident), sum(hashed_drops_high_12k$Doublet), 
                                              sum(!hashed_drops_high_12k$Confident & !hashed_drops_high_12k$Confident))),
                                 prop.table(c(sum(hashed_drops_high_8k$Confident), sum(hashed_drops_high_8k$Doublet), 
                                              sum(!hashed_drops_high_8k$Confident & !hashed_drops_high_8k$Confident))),
                                 prop.table(c(sum(hashed_drops_high_4k$Confident), sum(hashed_drops_high_4k$Doublet), 
                                              sum(!hashed_drops_high_4k$Confident & !hashed_drops_high_4k$Confident))),
                                 prop.table(c(sum(hashed_drops_med_12k$Confident), sum(hashed_drops_med_12k$Doublet), 
                                              sum(!hashed_drops_med_12k$Confident & !hashed_drops_med_12k$Confident))),
                                 prop.table(c(sum(hashed_drops_med_8k$Confident), sum(hashed_drops_med_8k$Doublet), 
                                              sum(!hashed_drops_med_8k$Confident & !hashed_drops_med_8k$Confident))),
                                 prop.table(c(sum(hashed_drops_med_4k$Confident), sum(hashed_drops_med_4k$Doublet), 
                                              sum(!hashed_drops_med_4k$Confident & !hashed_drops_med_4k$Confident))),
                                 prop.table(c(sum(hashed_drops_low_12k$Confident), sum(hashed_drops_low_12k$Doublet), 
                                              sum(!hashed_drops_low_12k$Confident & !hashed_drops_low_12k$Confident))),
                                 prop.table(c(sum(hashed_drops_low_8k$Confident), sum(hashed_drops_low_8k$Doublet), 
                                              sum(!hashed_drops_low_8k$Confident & !hashed_drops_low_8k$Confident))),
                                 prop.table(c(sum(hashed_drops_low_4k$Confident), sum(hashed_drops_low_4k$Doublet), 
                                              sum(!hashed_drops_low_4k$Confident & !hashed_drops_low_4k$Confident))))
colnames(hashed_drops_assignment) <- c("Singlet", "Doublet", "Negative")
rownames(hashed_drops_assignment) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                       "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(hashed_drops_assignment), las = 2, col = rev(viridis(3)), main = "HashedDrops")
legend(x = "topright", inset=c(-0.3,0), legend = c("Singlet", "Doublet", "Negative"), 
       col = rev(viridis(3)), pch = 15)


demuxmix_assignments <- rbind(prop.table(table(dmmClassify(demuxmix_high_12k)$Type)),
                              prop.table(table(dmmClassify(demuxmix_high_8k)$Type)),
                              prop.table(table(dmmClassify(demuxmix_high_4k)$Type)),
                              prop.table(table(dmmClassify(demuxmix_med_12k)$Type)),
                              prop.table(table(dmmClassify(demuxmix_med_8k)$Type)),
                              prop.table(table(dmmClassify(demuxmix_med_4k)$Type)),
                              prop.table(table(dmmClassify(demuxmix_low_12k)$Type)),
                              prop.table(table(dmmClassify(demuxmix_low_8k)$Type)),
                              prop.table(table(dmmClassify(demuxmix_low_4k)$Type)))
rownames(demuxmix_assignments) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                    "Med_4k", "Low_12k", "Low_8k", "Low_4k")
demuxmix_assignments <- demuxmix_assignments[, c(3, 1, 2, 4)]
barplot(t(demuxmix_assignments), las = 2, col = c(rev(viridis(3)), "grey"), main = "demuxmix")
legend(x = "topright", inset=c(-0.3,0), legend = c("Singlet", "Doublet", "Negative", "Uncertain"), 
       col = c(rev(viridis(3)), "grey"), pch = 15)

actual_assignments <- rbind(prop.table(c(grep("mp", rownames(high_12k))[1] - 1, length(grep("mp", rownames(high_12k))))),
                            prop.table(c(grep("mp", rownames(high_8k))[1] - 1, length(grep("mp", rownames(high_8k))))),
                            prop.table(c(grep("mp", rownames(high_4k))[1] - 1, length(grep("mp", rownames(high_4k))))),
                            prop.table(c(grep("mp", rownames(med_12k))[1] - 1, length(grep("mp", rownames(med_12k))))),
                            prop.table(c(grep("mp", rownames(med_8k))[1] - 1, length(grep("mp", rownames(med_8k))))),
                            prop.table(c(grep("mp", rownames(med_4k))[1] - 1, length(grep("mp", rownames(med_4k))))),
                            prop.table(c(grep("mp", rownames(low_12k))[1] - 1, length(grep("mp", rownames(low_12k))))),
                            prop.table(c(grep("mp", rownames(low_8k))[1] - 1, length(grep("mp", rownames(low_8k))))),
                            prop.table(c(grep("mp", rownames(low_4k))[1] - 1, length(grep("mp", rownames(low_4k))))))
rownames(actual_assignments) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                  "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(actual_assignments), las = 2, col = rev(viridis(3))[1:2], main = "Actual")
legend(x = "topright", inset=c(-0.3,0), legend = c("Singlet", "Doublet", "Negative", "Uncertain"), 
       col = c(rev(viridis(3)), "grey"), pch = 15)


#Singlet Assignment
high_12k_singlets <- high_12k[-grep("mp", rownames(high_12k)), ]
high_12k_multiplets <- high_12k[grep("mp", rownames(high_12k)), ]
high_8k_singlets <- high_8k[-grep("mp", rownames(high_8k)), ]
high_8k_multiplets <- high_8k[grep("mp", rownames(high_8k)), ]
high_4k_singlets <- high_4k[-grep("mp", rownames(high_4k)), ]
high_4k_multiplets <- high_4k[grep("mp", rownames(high_4k)), ]
med_12k_singlets <- med_12k[-grep("mp", rownames(med_12k)), ]
med_12k_multiplets <- med_12k[grep("mp", rownames(med_12k)), ]
med_8k_singlets <- med_8k[-grep("mp", rownames(med_8k)), ]
med_8k_multiplets <- med_8k[grep("mp", rownames(med_8k)), ]
med_4k_singlets <- med_4k[-grep("mp", rownames(med_4k)), ]
med_4k_multiplets <- med_4k[grep("mp", rownames(med_4k)), ]
low_12k_singlets <- low_12k[-grep("mp", rownames(low_12k)), ]
low_12k_multiplets <- low_12k[grep("mp", rownames(low_12k)), ]
low_8k_singlets <- low_8k[-grep("mp", rownames(low_8k)), ]
low_8k_multiplets <- low_8k[grep("mp", rownames(low_8k)), ]
low_4k_singlets <- low_4k[-grep("mp", rownames(low_4k)), ]
low_4k_multiplets <- low_4k[grep("mp", rownames(low_4k)), ]

error_rates_singlet <- function(x,
                                x_calls,
                                method) {
  true_assignment <- paste0("HTO-", sapply(strsplit(rownames(x), "-"), '[', 1))
  
  if(method == "bff_raw" | method == "bff_cluster" | method == "gmm_demux" | method == "htodemux"){
    n_wrong <- sum(true_assignment != x_calls$consensuscall[1:nrow(x)])
    n_doublet <- sum(x_calls$consensuscall[1:nrow(x)] == "Doublet")
    n_nocall <- sum(x_calls$consensuscall[1:nrow(x)] == "Negative")
    n_singlet <- n_wrong - n_doublet - n_nocall
  } else if(method == "demuxmix"){
    n_wrong <- sum(true_assignment != x_calls$HTO[1:nrow(x)])
    n_doublet <- length(grep(",", x_calls$HTO[1:nrow(x)]))
    n_nocall <- sum(x_calls$HTO[1:nrow(x)] == "uncertain" | x_calls$HTO[1:nrow(x)] == "negative")
    n_singlet <- n_wrong - n_doublet - n_nocall
  } else if(method == "hashed_drops"){
    n_doublet <- sum(x_calls$Doublet[1:nrow(x)])
    n_nocall <- sum(!x_calls$Confident[1:nrow(x)])
    n_singlet <- sum((true_assignment != paste0("HTO-", x_calls$Best[1:nrow(x)]))[!x_calls$Doublet[1:nrow(x)] & x_calls$Confident[1:nrow(x)]])
  }
  
  all_wrong <- c(n_singlet, n_doublet, n_nocall)
  names(all_wrong) <- c("Wrong Barcode", "Doublet", "No call")
  
  return(all_wrong/nrow(x) * 100)
}

bff_raw_error <- rbind(error_rates_singlet(x = high_12k_singlets, x_calls = bff_raw_high_12k, method = "bff_raw"),
                       error_rates_singlet(x = high_8k_singlets, x_calls = bff_raw_high_8k, method = "bff_raw"),
                       error_rates_singlet(x = high_4k_singlets, x_calls = bff_raw_high_4k, method = "bff_raw"),
                       error_rates_singlet(x = med_12k_singlets, x_calls = bff_raw_med_12k, method = "bff_raw"),
                       error_rates_singlet(x = med_8k_singlets, x_calls = bff_raw_med_8k, method = "bff_raw"),
                       error_rates_singlet(x = med_4k_singlets, x_calls = bff_raw_med_4k, method = "bff_raw"))
rownames(bff_raw_error) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                             "Med_4k")
barplot(t(bff_raw_error), col = c("red", "grey", "black"), las = 2, ylim = c(0, 40), main = "BFF Raw")


bff_cluster_error <- rbind(error_rates_singlet(x = high_12k_singlets, x_calls = bff_cluster_high_12k, method = "bff_cluster"),
                           error_rates_singlet(x = high_8k_singlets, x_calls = bff_cluster_high_8k, method = "bff_cluster"),
                           error_rates_singlet(x = high_4k_singlets, x_calls = bff_cluster_high_4k, method = "bff_cluster"),
                           error_rates_singlet(x = med_12k_singlets, x_calls = bff_cluster_med_12k, method = "bff_cluster"),
                           error_rates_singlet(x = med_8k_singlets, x_calls = bff_cluster_med_8k, method = "bff_cluster"),
                           error_rates_singlet(x = med_4k_singlets, x_calls = bff_cluster_med_4k, method = "bff_cluster"))
rownames(bff_cluster_error) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                 "Med_4k")
barplot(t(bff_cluster_error), col = c("red", "grey", "black"), las = 2, ylim = c(0, 40), main = "BFF Cluster")

gmm_demux_error <- rbind(error_rates_singlet(x = high_12k_singlets, x_calls = gmm_demux_high_12k, method = "gmm_demux"),
                         error_rates_singlet(x = high_8k_singlets, x_calls = gmm_demux_high_8k, method = "gmm_demux"),
                         error_rates_singlet(x = high_4k_singlets, x_calls = gmm_demux_high_4k, method = "gmm_demux"),
                         error_rates_singlet(x = med_12k_singlets, x_calls = gmm_demux_med_12k, method = "gmm_demux"),
                         error_rates_singlet(x = med_8k_singlets, x_calls = gmm_demux_med_8k, method = "gmm_demux"),
                         error_rates_singlet(x = med_4k_singlets, x_calls = gmm_demux_med_4k, method = "gmm_demux"),
                         error_rates_singlet(x = low_12k_singlets, x_calls = gmm_demux_low_12k, method = "gmm_demux"),
                         error_rates_singlet(x = low_8k_singlets, x_calls = gmm_demux_low_8k, method = "gmm_demux"),
                         error_rates_singlet(x = low_4k_singlets, x_calls = gmm_demux_low_4k, method = "gmm_demux"))
rownames(gmm_demux_error) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                               "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(gmm_demux_error), col = c("red", "grey", "black"), las = 2, ylim = c(0, 40), main = "GMM Demux")

htodemux_error <- rbind(error_rates_singlet(x = high_12k_singlets, x_calls = htodemux_high_12k, method = "htodemux"),
                        error_rates_singlet(x = high_8k_singlets, x_calls = htodemux_high_8k, method = "htodemux"),
                        error_rates_singlet(x = high_4k_singlets, x_calls = htodemux_high_4k, method = "htodemux"),
                        error_rates_singlet(x = med_12k_singlets, x_calls = htodemux_med_12k, method = "htodemux"),
                        error_rates_singlet(x = med_8k_singlets, x_calls = htodemux_med_8k, method = "htodemux"),
                        error_rates_singlet(x = med_4k_singlets, x_calls = htodemux_med_4k, method = "htodemux"),
                        error_rates_singlet(x = low_12k_singlets, x_calls = htodemux_low_12k, method = "htodemux"),
                        error_rates_singlet(x = low_8k_singlets, x_calls = htodemux_low_8k, method = "htodemux"),
                        error_rates_singlet(x = low_4k_singlets, x_calls = htodemux_low_4k, method = "htodemux"))
rownames(htodemux_error) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                              "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(htodemux_error), col = c("red", "grey", "black"), las = 2, ylim = c(0, 40), main = "HTO Demux")

demuxmix_error <- rbind(error_rates_singlet(x = high_12k_singlets, x_calls = dmmClassify(demuxmix_high_12k), method = "demuxmix"),
                        error_rates_singlet(x = high_8k_singlets, x_calls = dmmClassify(demuxmix_high_8k), method = "demuxmix"),
                        error_rates_singlet(x = high_4k_singlets, x_calls = dmmClassify(demuxmix_high_4k), method = "demuxmix"),
                        error_rates_singlet(x = med_12k_singlets, x_calls = dmmClassify(demuxmix_med_12k), method = "demuxmix"),
                        error_rates_singlet(x = med_8k_singlets, x_calls = dmmClassify(demuxmix_med_8k), method = "demuxmix"),
                        error_rates_singlet(x = med_4k_singlets, x_calls = dmmClassify(demuxmix_med_4k), method = "demuxmix"),
                        error_rates_singlet(x = low_12k_singlets, x_calls = dmmClassify(demuxmix_low_12k), method = "demuxmix"),
                        error_rates_singlet(x = low_8k_singlets, x_calls = dmmClassify(demuxmix_low_8k), method = "demuxmix"),
                        error_rates_singlet(x = low_4k_singlets, x_calls = dmmClassify(demuxmix_low_4k), method = "demuxmix"))
rownames(demuxmix_error) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                              "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(demuxmix_error), col = c("red", "grey", "black"), las = 2, ylim = c(0, 40), main = "Demuxmix")


hashed_drops_error <- rbind(error_rates_singlet(x = high_12k_singlets, x_calls = hashed_drops_high_12k, method = "hashed_drops"),
                            error_rates_singlet(x = high_8k_singlets, x_calls = hashed_drops_high_8k, method = "hashed_drops"),
                            error_rates_singlet(x = high_4k_singlets, x_calls = hashed_drops_high_4k, method = "hashed_drops"),
                            error_rates_singlet(x = med_12k_singlets, x_calls = hashed_drops_med_12k, method = "hashed_drops"),
                            error_rates_singlet(x = med_8k_singlets, x_calls = hashed_drops_med_8k, method = "hashed_drops"),
                            error_rates_singlet(x = med_4k_singlets, x_calls = hashed_drops_med_4k, method = "hashed_drops"),
                            error_rates_singlet(x = low_12k_singlets, x_calls = hashed_drops_low_12k, method = "hashed_drops"),
                            error_rates_singlet(x = low_8k_singlets, x_calls = hashed_drops_low_8k, method = "hashed_drops"),
                            error_rates_singlet(x = low_4k_singlets, x_calls = hashed_drops_low_4k, method = "hashed_drops"))
rownames(hashed_drops_error) <- c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k",
                                  "Med_4k", "Low_12k", "Low_8k", "Low_4k")
barplot(t(hashed_drops_error), col = c("red", "grey", "black"), las = 2, ylim = c(0, 40), main = "Hashed Drops")
legend(x = "topright", inset=c(-0.3,0), legend = c("Wrong Singlet Call", "Doublet", "No Call"),
       col = c("red", "grey", "black"), pch = 15)


# Multiplet Assignment Accuracy
multiplet_assignment <- function(x_calls,
                                 method){
  
  assignment_same <- sapply(X = 1:length(x_calls), function(x){
    if(method == "bff_raw" | method == "bff_cluster" | method == "gmm_demux" | method == "htodemux"){
      x_calls_mp <- x_calls[[x]][grep("mp", x_calls[[x]]$cellbarcode), ]
      matching <- sapply(strsplit(sapply(strsplit(x_calls_mp$cellbarcode, "\\."), '[', 2), "-"), '[', 1) ==
        sapply(strsplit(sapply(strsplit(x_calls_mp$cellbarcode, "\\."), '[', 2), "-"), '[', 2)
      
      assignment_same <- prop.table(table(x_calls_mp$consensuscall.global[matching]))
    } else if(method == "demuxmix"){
      x_calls_mp <- dmmClassify(x_calls[[x]])[grep("mp", rownames(dmmClassify(x_calls[[x]]))), ]
      matching <- sapply(strsplit(sapply(strsplit(rownames(x_calls_mp), "\\."), '[', 2), "-"), '[', 1) ==
        sapply(strsplit(sapply(strsplit(rownames(x_calls_mp), "\\."), '[', 2), "-"), '[', 2)
      
      x_calls_mp$Type[x_calls_mp$Type == "uncertain"] <- "negative"
      calls_factor <- factor(x = x_calls_mp$Type, levels = c("multiplet", "negative", "singlet"))
      
      assignment_diff <- prop.table(table(calls_factor[matching]))
    } else if (method == "hashed_drops"){
      x_calls_mp <- x_calls[[x]][grep("mp", rownames(x_calls[[x]])), ]
      matching <- sapply(strsplit(sapply(strsplit(rownames(x_calls_mp), "\\."), '[', 2), "-"), '[', 1) ==
        sapply(strsplit(sapply(strsplit(rownames(x_calls_mp), "\\."), '[', 2), "-"), '[', 2)
      
      calls <- rep("negative", nrow(x_calls_mp))
      calls[x_calls_mp$Doublet] <- "doublet"
      calls[x_calls_mp$Confident] <- "singlet"
      
      calls <- factor(x = calls, levels = c("doublet", "singlet", "negative"))
      assignment_same <- prop.table(table(calls[matching]))
    }
  })
  
  assignment_diff <- sapply(X = 1:length(x_calls), function(x){
    if(method == "bff_raw" | method == "bff_cluster" | method == "gmm_demux" | method == "htodemux"){
      x_calls_mp <- x_calls[[x]][grep("mp", x_calls[[x]]$cellbarcode), ]
      matching <- sapply(strsplit(sapply(strsplit(x_calls_mp$cellbarcode, "\\."), '[', 2), "-"), '[', 1) ==
        sapply(strsplit(sapply(strsplit(x_calls_mp$cellbarcode, "\\."), '[', 2), "-"), '[', 2)
      
      assignment_diff <- prop.table(table(x_calls_mp$consensuscall.global[!matching]))
    } else if(method == "demuxmix"){
      x_calls_mp <- dmmClassify(x_calls[[x]])[grep("mp", rownames(dmmClassify(x_calls[[x]]))), ]
      matching <- sapply(strsplit(sapply(strsplit(rownames(x_calls_mp), "\\."), '[', 2), "-"), '[', 1) ==
        sapply(strsplit(sapply(strsplit(rownames(x_calls_mp), "\\."), '[', 2), "-"), '[', 2)
      
      x_calls_mp$Type[x_calls_mp$Type == "uncertain"] <- "negative"
      calls_factor <- factor(x = x_calls_mp$Type, levels = c("multiplet", "negative", "singlet"))
      
      assignment_diff <- prop.table(table(calls_factor[!matching]))
    } else if (method == "hashed_drops"){
      x_calls_mp <- x_calls[[x]][grep("mp", rownames(x_calls[[x]])), ]
      matching <- sapply(strsplit(sapply(strsplit(rownames(x_calls_mp), "\\."), '[', 2), "-"), '[', 1) ==
        sapply(strsplit(sapply(strsplit(rownames(x_calls_mp), "\\."), '[', 2), "-"), '[', 2)
      
      calls <- rep("negative", nrow(x_calls_mp))
      calls[x_calls_mp$Doublet] <- "doublet"
      calls[x_calls_mp$Confident] <- "singlet"
      
      calls <- factor(x = calls, levels = c("doublet", "singlet", "negative"))
      assignment_diff <- prop.table(table(calls[!matching]))
    }
  })
  print(assignment_diff)
  
  assignment_same <- assignment_same[, ncol(assignment_same):1]
  assignment_diff <- assignment_diff[, ncol(assignment_diff):1]
  
  tables <- list(assignment_same, assignment_diff)
  names(tables) <- c("Same", "Diff")
  return(tables)
}

bff_raw_mp_assign <- multiplet_assignment(list(bff_raw_high_12k, bff_raw_high_8k, bff_raw_high_4k,
                                               bff_raw_med_12k, bff_raw_med_8k, bff_raw_med_4k), method = "bff_raw")
bff_raw_mp_assign[["Same"]] <- cbind(c(0, 1, 0), c(0, 1, 0), c(0, 1, 0), bff_raw_mp_assign[["Same"]])
bff_raw_mp_assign[["Diff"]] <- cbind(c(0, 0, 0), c(0, 1, 0), c(0, 1, 0), c(0, 1, 0), bff_raw_mp_assign[["Diff"]])

colnames(bff_raw_mp_assign[["Same"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                               "Med_4k", "Low_12k", "Low_8k", "Low_4k"))
colnames(bff_raw_mp_assign[["Diff"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                               "Med_4k", "Low_12k", "Low_8k", "Low_4k", ""))

par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
barplot(cbind(bff_raw_mp_assign[["Same"]], bff_raw_mp_assign[["Diff"]]), las = 2,  
        col = c("red", "grey", "black"), horiz = T, main = "BFF_Raw")
legend(x = "topright", inset=c(-0.3,0), legend = c("Doublet", "Negative", "Singlet"), 
       col = c("red", "grey", "black"), pch = 15)


bff_cluster_mp_assign <- multiplet_assignment(list(bff_cluster_high_12k, bff_cluster_high_8k, bff_cluster_high_4k,
                                                   bff_cluster_med_12k, bff_cluster_med_8k, bff_cluster_med_4k), method = "bff_cluster")
bff_cluster_mp_assign[["Same"]] <- cbind(c(0, 1, 0), c(0, 1, 0), c(0, 1, 0), bff_cluster_mp_assign[["Same"]])
bff_cluster_mp_assign[["Diff"]] <- cbind(c(0, 0, 0), c(0, 1, 0), c(0, 1, 0), c(0, 1, 0), bff_cluster_mp_assign[["Diff"]])

colnames(bff_cluster_mp_assign[["Same"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                   "Med_4k", "Low_12k", "Low_8k", "Low_4k"))
colnames(bff_cluster_mp_assign[["Diff"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                   "Med_4k", "Low_12k", "Low_8k", "Low_4k", ""))

par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
barplot(cbind(bff_cluster_mp_assign[["Same"]], bff_cluster_mp_assign[["Diff"]]), las = 2,  
        col = c("red", "grey", "black"), horiz = T, main = "BFF_Cluster")
legend(x = "topright", inset=c(-0.3,0), legend = c("Doublet", "Negative", "Singlet"), 
       col = c("red", "grey", "black"), pch = 15)


gmm_demux_mp_assign <- multiplet_assignment(list(gmm_demux_high_12k, gmm_demux_high_8k, gmm_demux_high_4k,
                                                 gmm_demux_med_12k, gmm_demux_med_8k, gmm_demux_med_4k, 
                                                 gmm_demux_low_12k, gmm_demux_low_8k, gmm_demux_low_4k), method = "gmm_demux")
gmm_demux_mp_assign[["Diff"]] <- cbind(c(0, 0, 0), gmm_demux_mp_assign[["Diff"]])

colnames(gmm_demux_mp_assign[["Same"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                 "Med_4k", "Low_12k", "Low_8k", "Low_4k"))
colnames(gmm_demux_mp_assign[["Diff"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                 "Med_4k", "Low_12k", "Low_8k", "Low_4k", ""))

par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
barplot(cbind(gmm_demux_mp_assign[["Same"]], gmm_demux_mp_assign[["Diff"]]), las = 2,  
        col = c("red", "grey", "black"), horiz = T, main = "gmm_demux")
legend(x = "topright", inset=c(-0.3,0), legend = c("Doublet", "Negative", "Singlet"), 
       col = c("red", "grey", "black"), pch = 15)


htodemux_mp_assign <- multiplet_assignment(list(htodemux_high_12k, htodemux_high_8k, htodemux_high_4k,
                                                htodemux_med_12k, htodemux_med_8k, htodemux_med_4k, 
                                                htodemux_low_12k, htodemux_low_8k, htodemux_low_4k), method = "htodemux")
htodemux_mp_assign[["Diff"]] <- cbind(c(0, 0, 0), htodemux_mp_assign[["Diff"]])

colnames(htodemux_mp_assign[["Same"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                "Med_4k", "Low_12k", "Low_8k", "Low_4k"))
colnames(htodemux_mp_assign[["Diff"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                "Med_4k", "Low_12k", "Low_8k", "Low_4k", ""))

par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
barplot(cbind(htodemux_mp_assign[["Same"]], htodemux_mp_assign[["Diff"]]), las = 2,  
        col = c("red", "grey", "black"), horiz = T, main = "htodemux")
legend(x = "topright", inset=c(-0.3,0), legend = c("Doublet", "Negative", "Singlet"), 
       col = c("red", "grey", "black"), pch = 15)


demuxmix_mp_assign <- multiplet_assignment(list(demuxmix_high_12k, demuxmix_high_8k, demuxmix_high_4k,
                                                demuxmix_med_12k, demuxmix_med_8k, demuxmix_med_4k, 
                                                demuxmix_low_12k, demuxmix_low_8k, demuxmix_low_4k), method = "demuxmix")
demuxmix_mp_assign[["Diff"]] <- cbind(c(0, 0, 0), demuxmix_mp_assign[["Diff"]])

colnames(demuxmix_mp_assign[["Same"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                "Med_4k", "Low_12k", "Low_8k", "Low_4k"))
colnames(demuxmix_mp_assign[["Diff"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                "Med_4k", "Low_12k", "Low_8k", "Low_4k", ""))

par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
barplot(cbind(demuxmix_mp_assign[["Same"]], demuxmix_mp_assign[["Diff"]]), las = 2,  
        col = c("red", "grey", "black"), horiz = T, main = "demuxmix")
legend(x = "topright", inset=c(-0.3,0), legend = c("Doublet", "Negative", "Singlet"), 
       col = c("red", "grey", "black"), pch = 15)


hashed_drops_mp_assign <- multiplet_assignment(list(hashed_drops_high_12k, hashed_drops_high_8k, hashed_drops_high_4k,
                                                    hashed_drops_med_12k, hashed_drops_med_8k, hashed_drops_med_4k, 
                                                    hashed_drops_low_12k, hashed_drops_low_8k, hashed_drops_low_4k), method = "hashed_drops")
hashed_drops_mp_assign[["Diff"]] <- cbind(c(0, 0, 0), hashed_drops_mp_assign[["Diff"]])

colnames(hashed_drops_mp_assign[["Same"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                    "Med_4k", "Low_12k", "Low_8k", "Low_4k"))
colnames(hashed_drops_mp_assign[["Diff"]]) <- rev(c("High_12k", "High_8k", "High_4k", "Med_12k", "Med_8k", 
                                                    "Med_4k", "Low_12k", "Low_8k", "Low_4k", ""))

par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
barplot(cbind(hashed_drops_mp_assign[["Same"]], hashed_drops_mp_assign[["Diff"]]), las = 2,  
        col = c("red", "grey", "black"), horiz = T, main = "hashed_drops")
legend(x = "topright", inset=c(-0.3,0), legend = c("Doublet", "Negative", "Singlet"), 
       col = c("red", "grey", "black"), pch = 15)


# Calculate F-score
f_score <- function(x, x_calls, method){
  if(method == "bff_raw" | method == "bff_cluster" | method == "gmm_demux" | method == "htodemux"){
    true_assignment <- c(paste0("HTO-", sapply(strsplit(rownames(x)[!grepl("mp", rownames(x))], "-"), '[', 1)),
                         rep("Doublet", length(grep("mp", rownames(x)))))
    fn <- sum(x_calls$consensuscall == "Negative")
    fp <- sum(x_calls$consensuscall != true_assignment) - fn
  } else if(method == "demuxmix"){
    true_assignment <- c(paste0("HTO-", sapply(strsplit(rownames(x)[!grepl("mp", rownames(x))], "-"), '[', 1)),
                         rep("multiplet", length(grep("mp", rownames(x)))))
    fn <- sum(dmmClassify(x_calls)$HTO == "negative" | dmmClassify(x_calls)$HTO == "uncertain")
    singlet_index <- !grepl("mp", rownames(dmmClassify(x_calls)))
    fp <- sum(dmmClassify(x_calls)$HTO[singlet_index] != true_assignment[singlet_index], 
              dmmClassify(x_calls)$Type[!singlet_index] != true_assignment[!singlet_index]) - fn
  } else if(method == "hashed_drops"){
    true_assignment <- c(paste0("HTO-", sapply(strsplit(rownames(x)[!grepl("mp", rownames(x))], "-"), '[', 1)),
                         rep("Doublet", length(grep("mp", rownames(x)))))
    
    hashed_assignment <- rep("negative", nrow(x_calls))
    hashed_assignment[x_calls$Doublet] <- "Doublet"
    hashed_assignment[x_calls$Confident] <- paste0("HTO-", x_calls$Best[x_calls$Confident])
    
    fn <- sum(hashed_assignment == "negative")
    fp <- sum(hashed_assignment != true_assignment) - fn
  }
  tp <- length(true_assignment) - fn - fp
  f <- tp / (tp + 0.5 * (fp + fn))
  
  return(f)
}

par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)

high_12k_fscore <- c(f_score(high_12k, bff_raw_high_12k, method = "bff_raw"),
                     f_score(high_12k, bff_cluster_high_12k, method = "bff_cluster"),
                     f_score(high_12k, gmm_demux_high_12k, method = "gmm_demux"),
                     f_score(high_12k, htodemux_high_12k, method = "htodemux"),
                     f_score(high_12k, hashed_drops_high_12k, method = "hashed_drops"),
                     f_score(high_12k, demuxmix_high_12k, method = "demuxmix"))
names(high_12k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(high_12k_fscore))
cols[which.max(high_12k_fscore)] <- "red"

barplot(high_12k_fscore, ylim = c(0,1), las = 2, col = cols, main = "High 12k F-scores")

high_8k_fscore <- c(f_score(high_8k, bff_raw_high_8k, method = "bff_raw"),
                     f_score(high_8k, bff_cluster_high_8k, method = "bff_cluster"),
                     f_score(high_8k, gmm_demux_high_8k, method = "gmm_demux"),
                     f_score(high_8k, htodemux_high_8k, method = "htodemux"),
                     f_score(high_8k, hashed_drops_high_8k, method = "hashed_drops"),
                     f_score(high_8k, demuxmix_high_8k, method = "demuxmix"))
names(high_8k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(high_8k_fscore))
cols[which.max(high_8k_fscore)] <- "red"

barplot(high_8k_fscore, ylim = c(0,1), las = 2, col = cols, main = "High 8k F-scores")

high_4k_fscore <- c(f_score(high_4k, bff_raw_high_4k, method = "bff_raw"),
                    f_score(high_4k, bff_cluster_high_4k, method = "bff_cluster"),
                    f_score(high_4k, gmm_demux_high_4k, method = "gmm_demux"),
                    f_score(high_4k, htodemux_high_4k, method = "htodemux"),
                    f_score(high_4k, hashed_drops_high_4k, method = "hashed_drops"),
                    f_score(high_4k, demuxmix_high_4k, method = "demuxmix"))
names(high_4k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(high_4k_fscore))
cols[which.max(high_4k_fscore)] <- "red"

barplot(high_4k_fscore, ylim = c(0,1), las = 2, col = cols, main = "High 4k F-scores")


med_12k_fscore <- c(f_score(med_12k, bff_raw_med_12k, method = "bff_raw"),
                    f_score(med_12k, bff_cluster_med_12k, method = "bff_cluster"),
                    f_score(med_12k, gmm_demux_med_12k, method = "gmm_demux"),
                    f_score(med_12k, htodemux_med_12k, method = "htodemux"),
                    f_score(med_12k, hashed_drops_med_12k, method = "hashed_drops"),
                    f_score(med_12k, demuxmix_med_12k, method = "demuxmix"))
names(med_12k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(med_12k_fscore))
cols[which.max(med_12k_fscore)] <- "red"

barplot(med_12k_fscore, ylim = c(0,1), las = 2, col = cols, main = "Med 12k F-scores")

med_8k_fscore <- c(f_score(med_8k, bff_raw_med_8k, method = "bff_raw"),
                    f_score(med_8k, bff_cluster_med_8k, method = "bff_cluster"),
                    f_score(med_8k, gmm_demux_med_8k, method = "gmm_demux"),
                    f_score(med_8k, htodemux_med_8k, method = "htodemux"),
                    f_score(med_8k, hashed_drops_med_8k, method = "hashed_drops"),
                    f_score(med_8k, demuxmix_med_8k, method = "demuxmix"))
names(med_8k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(med_8k_fscore))
cols[which.max(med_8k_fscore)] <- "red"

barplot(med_8k_fscore, ylim = c(0,1), las = 2, col = cols, main = "Med 8k F-scores")


med_4k_fscore <- c(f_score(med_4k, bff_raw_med_4k, method = "bff_raw"),
                   f_score(med_4k, bff_cluster_med_4k, method = "bff_cluster"),
                   f_score(med_4k, gmm_demux_med_4k, method = "gmm_demux"),
                   f_score(med_4k, htodemux_med_4k, method = "htodemux"),
                   f_score(med_4k, hashed_drops_med_4k, method = "hashed_drops"),
                   f_score(med_4k, demuxmix_med_4k, method = "demuxmix"))
names(med_4k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(med_4k_fscore))
cols[which.max(med_4k_fscore)] <- "red"

barplot(med_4k_fscore, ylim = c(0,1), las = 2, col = cols, main = "Med 4k F-scores")


low_12k_fscore <- c(f_score(low_12k, bff_raw_low_12k, method = "bff_raw"),
                   0,
                   f_score(low_12k, gmm_demux_low_12k, method = "gmm_demux"),
                   f_score(low_12k, htodemux_low_12k, method = "htodemux"),
                   f_score(low_12k, hashed_drops_low_12k, method = "hashed_drops"),
                   f_score(low_12k, demuxmix_low_12k, method = "demuxmix"))
names(low_12k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(low_12k_fscore))
cols[which.max(low_12k_fscore)] <- "red"

barplot(low_12k_fscore, ylim = c(0,1), las = 2, col = cols, main = "Low 12k F-scores")


low_8k_fscore <- c(f_score(low_8k, bff_raw_low_8k, method = "bff_raw"),
                    0,
                    f_score(low_8k, gmm_demux_low_8k, method = "gmm_demux"),
                    f_score(low_8k, htodemux_low_8k, method = "htodemux"),
                    f_score(low_8k, hashed_drops_low_8k, method = "hashed_drops"),
                    f_score(low_8k, demuxmix_low_8k, method = "demuxmix"))
names(low_8k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(low_8k_fscore))
cols[which.max(low_8k_fscore)] <- "red"

barplot(low_8k_fscore, ylim = c(0,1), las = 2, col = cols, main = "Low 8k F-scores")


low_4k_fscore <- c(f_score(low_4k, bff_raw_low_4k, method = "bff_raw"),
                   0,
                   f_score(low_4k, gmm_demux_low_4k, method = "gmm_demux"),
                   f_score(low_4k, htodemux_low_4k, method = "htodemux"),
                   f_score(low_4k, hashed_drops_low_4k, method = "hashed_drops"),
                   f_score(low_4k, demuxmix_low_4k, method = "demuxmix"))
names(low_4k_fscore) <- c("bff_raw", "bff_cluster", "gmm_demux", "htodemux", "hashed_drops", "demuxmix")
cols <- rep("black", length(low_4k_fscore))
cols[which.max(low_4k_fscore)] <- "red"

barplot(low_4k_fscore, ylim = c(0,1), las = 2, col = cols, main = "Low 4k F-scores")


fscores_combined <- rbind(high_12k_fscore, high_8k_fscore, high_4k_fscore, med_12k_fscore, med_8k_fscore, med_4k_fscore,
                          low_12k_fscore, low_8k_fscore, low_4k_fscore)
fscores_means <- colMeans(fscores_combined)

cols <- rep("black", length(fscores_means))
cols[which.max(fscores_means)] <- "red"
barplot(fscores_means, ylim = c(0,1), las = 2, col = cols, main = "Mean F scores")

# PPT Plots

layout(t(rbind(c(1, 2, 3))))
p1 <- hist(high_12k[, 1], breaks = 100, xlim = c(0, 200), xlab = "HTO-1", main = "High Quality", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
p2 <- hist(med_12k[, 1], breaks = 100, xlim = c(0, 200), xlab = "HTO-1", main = "Medium Quality", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
p3 <- hist(low_12k[, 1], breaks = 100, xlim = c(0, 200), xlab = "HTO-1", main = "Low Quality", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)

plot(high_12k[, 1], high_12k[, 2], pch = 16, cex = 0.5, main = "High quality - 12,000 cells", xlab = "HTO-1", ylab = "HTO-2")
plot(log(high_12k[, 1]), log(high_12k[, 2]), pch = 16, cex = 0.5, main = "High quality - 12,000 cells, log-scaled", xlab = "HTO-1", ylab = "HTO-2")
tail(high_12k, n = 20)

cols_h12k <- rep("black", nrow(high_12k))
cols_h12k[c(grep("mp.*1-2", rownames(high_12k)), grep("mp.*1-1", rownames(high_12k)), grep("mp.*2-2", rownames(high_12k)))] <- "red"

plot(high_12k[, 1], high_12k[, 2], pch = 16, cex = 0.5, main = "High quality - 12,000 cells", 
     xlab = "HTO-1", ylab = "HTO-2", col = cols_h12k)

layout(matrix(1:9, nrow = 3, ncol = 3, byrow = T))
par(mar = c(1, 1, 1, 1))
plot(log(high_12k[, 1]), log(high_12k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
plot(log(high_8k[, 1]), log(high_8k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
plot(log(high_4k[, 1]), log(high_4k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
plot(log(med_12k[, 1]), log(med_12k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
plot(log(med_8k[, 1]), log(med_8k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
plot(log(med_4k[, 1]), log(med_4k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
plot(log(low_12k[, 1]), log(low_12k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
plot(log(low_8k[, 1]), log(low_8k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
plot(log(low_4k[, 1]), log(low_4k[, 2]), pch = 16, cex = 0.2, 
     xlab = "HTO-1", ylab = "HTO-2", xlim = c(0, 6), ylim = c(0, 6), labels = F, yaxt = 'n', xaxt = 'n')
dev.off()

layout(rbind(1, 2))
hist(rnbinom(n = 12000, mu = 100, size = 10), col = c(rep("red", 12000)), xlim = c(0, 250), breaks = 50, main = "Positive distribution")
hist(rnbinom(n = 12000, mu = 20, size = 10), col = c(rep("blue", 12000)), xlim = c(0, 250), breaks = 20, main = "Negative distribution")
dev.off()

cell1 <- c(rnbinom(n = 1, mu = 100, size = 10), rnbinom(n = 1, mu = 20, size = 10), 
           rnbinom(n = 1, mu = 20, size = 10), rnbinom(n = 1, mu = 20, size = 10), 
           rnbinom(n = 1, mu = 20, size = 10), rnbinom(n = 1, mu = 20, size = 10), 
           rnbinom(n = 1, mu = 20, size = 10), rnbinom(n = 1, mu = 20, size = 10))

cell2 <- c(rnbinom(n = 1, mu = 20, size = 10), rnbinom(n = 1, mu = 20, size = 10), 
           rnbinom(n = 1, mu = 20, size = 10), rnbinom(n = 1, mu = 20, size = 10), 
           rnbinom(n = 1, mu = 20, size = 10), rnbinom(n = 1, mu = 100, size = 10), 
           rnbinom(n = 1, mu = 20, size = 10), rnbinom(n = 1, mu = 20, size = 10))

names(cell1) <- paste0("HTO-", 1:8)
names(cell2) <- paste0("HTO-", 1:8)

layout(rbind(1, 2))
barplot(cell1, las = 2, main = "Cell positive for HTO-1")
barplot(cell2, las = 2, main = "Cell positive for HTO-6")
dev.off()

barplot(cell1 + cell2, las = 2, main = "Doublet positive for both HTO-1 and HTO-6")




