library(dplyr)
c('ggplot2', 'HGC', 'reshape2', 'Seurat', 'zellkonverter') %>%
  pacman::p_load(char = .)

adata.wmc <- readH5AD('220516-ABM.velo.wmc.h5ad')
wmc.clustering.trees <- sapply(
  adata.wmc@assays %>% names %>% .[grepl('wmc_', .) | . == 'X'],
  function(wmc) {
    as.Seurat(adata.wmc, counts = wmc, data = wmc) %>%
    FindVariableFeatures(., nfeatures = 2000) %>%
    ScaleData(., features = rownames(.)) %>%
    RunPCA(., features = VariableFeatures(object = .)) %>%
    # Determine dimensionality
    JackStraw(num.replicate = 100) %>%
    ScoreJackStraw(dims = 1:20) %>%
    # Construct the graph and cluster the cells with HGC
    FindNeighbors(dims = 1:10) %>%
    FindClusteringTree(graph.type = "SNN")
  }
)
saveRDS(wmc.clustering.trees, 'wmc_clustering_trees.rds')
wmc.clustering.trees <- readRDS('wmc_clustering_trees.rds')

# Save dendrograms and ARI plot
sapply(names(wmc.clustering.trees), function(wmc) {
  wmc.clustering.trees[[wmc]] %>% {
    .@graphs$ClusteringTree$height <- log(.@graphs$ClusteringTree$height + 1)
    tree <- .@graphs$ClusteringTree
    meta.data <- .@meta.data[c('phase', 'ann0608')] %>%
      rename(all_of(c('Cell Phase' = 'phase', 'Cell Type' = 'ann0608')))
    plot_fname <- ifelse(wmc == 'X', 'NoWMC', sub('^[^_]*_', '', wmc))
    res <- 300
    length <- 480 / 72 * res
    png(
      file.path('DendrogramPlots', sprintf('%s.png', plot_fname)),
      res = res, width = length, height = length
    )
    par(mar = c(3, 3, 0, 0))
    HGC.PlotDendrogram(tree, plot.label = TRUE, labels = meta.data)
    dev.off()
    HGC.PlotARIs(tree, labels = meta.data, return.ARI = TRUE) %>%
      data.frame %>% mutate(k = 1:nrow(.) + 1,
        WMC = ifelse(wmc == 'X', 'No WMC',
                     wmc %>% sub('wmc_', '', .) %>% sub('band', '-Band', .) %>%
                       sub('_low', ' Low', .) %>% sub('_high', ' High ', .)))
  }
}, simplify = FALSE) %>% bind_rows %>%
  melt(id.vars = c('k', 'WMC'), variable.name = 'Cell', value.name = 'ARI') %>%
  mutate(Cell = sub('Cell.', '', Cell)) %>% {
    ggplot(., aes(x = k, y = ARI, color = WMC, shape = Cell)) +
      geom_line() + geom_point()
  } %>% ggsave('ARI.pdf', .)