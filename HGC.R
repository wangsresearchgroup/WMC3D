library(dplyr)
c('data.table', 'ggplot2', 'grid', 'HGC', 'ggpubr', 'reshape2', 'Seurat') %>%
  pacman::p_load(char = .)
source('WMC3d.R')

rename.files <- function(dir) {
  sapply(list.files(dir), function(file) {
    file %>%
      gsub('.*barcodes', 'barcodes', .) %>%
      gsub('.*genes', 'genes', .) %>%
      gsub('.*.mtx', 'matrix.mtx', .) %>%
      gsub('.*.csv', 'metadata.csv', .) %>%
      {file.rename(paste(dir, file, sep = '/'), paste(dir, ., sep = '/'))}
  })
}

read.seurat <- function(dir) {
  Read10X(
    data.dir = dir, strip.suffix = TRUE,
    gene.column = read.csv(file.path(dir, 'genes.tsv'), sep = '\t') %>% ncol) %>%
    CreateSeuratObject(
      project = dir, min.cells = 3, min.features = 200,
      meta.data = read.csv(file.path(dir, 'metadata.csv')))
}

quality.control <- function(seurat) {
  seurat %>%
    subset(nFeature_RNA > 200 & nFeature_RNA < 2500) %>% {
      if ('percent.mito' %in% names(.@meta.data)) subset(., percent.mito < 5)
      else .
    }
}

preprocess <- function(seurat) {
  seurat %>%
    NormalizeData(normalization.method = "LogNormalize",
                  scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
}

wmc <- function(seurat, M.range = 2:4) {
  matrix <- seurat@assays$RNA@layers$data %>% as.matrix
  wmc_list <- list('No WMC' = seurat)
  for (M in M.range) {
    wmc_list <- mDWT(matrix, M) %>% {
      names(.) <- paste0(M, '-Band ', names(.))
      sapply(., function(x) {
        wmc_seurat <- seurat
        wmc_seurat@assays$RNA@layers$data <- x
        wmc_seurat@assays$RNA@layers$counts <- x
        wmc_seurat
      })
    } %>% append(wmc_list, .)
  }
  return(wmc_list)
}

dim.reduce <- function(seurat) {
  seurat %>%
    ScaleData(features = rownames(.)) %>%
    RunPCA(features = VariableFeatures(object = .)) %>%
    JackStraw(num.replicate = 100) %>%
    ScoreJackStraw(dims = 1:20)
}

clustering <- function(seurat) {
  seurat %>%
    FindNeighbors(dims = 1:10) %>%
    RunUMAP(n.components = 3, dims = 1:10) %>%
    FindClusteringTree(graph.type = "SNN")
}

dirs <- c('CID3921', 'CID4463', 'CID4495', 'CID4523', 'RT1')
sapply(dirs, function(dir) {
  rename.files(dir)
  dir %>% read.seurat %>% quality.control %>% preprocess %>% wmc %>%
    sapply(function(seurat) seurat %>% dim.reduce %>% clustering) %>%
    saveRDS(file.path(dir, 'wmc_clustering_trees.rds'))
})

# UMAP
dirs.umap <- dirs %>% subset(grepl('CID', .))
sapply(dirs, function(dir) {
  wmc.clustering.trees <- readRDS(file.path(dir, 'wmc_clustering_trees.rds'))
  sapply(names(wmc.clustering.trees), function(wmc) {
    dt <- wmc.clustering.trees[[wmc]]@reductions[["umap"]]@cell.embeddings %>%
      data.table
    setnames(dt, names(dt), names(dt) %>% toupper)
    dt[, WMC := wmc]
    wmc.clustering.trees[[wmc]]@meta.data %>% {
      dt[, `:=` (subtype = .[['subtype']], celltype = .[['celltype_major']])]
    }
    return(dt)
  }, simplify = FALSE) %>% rbindlist %>% fwrite(file.path(dir, 'UMAP.csv'))
})

plot.umap <- function(dir) {
  fread(file.path(dir, 'UMAP.csv')) %>% {
    setnames(., c('celltype', 'subtype'), c('Cell Type', 'Subtype'))
    # Set subplot order
    gaps <- paste('Gap', c(2, 5))
    dummy <- .[rep(1, length(gaps))]
    dummy[, WMC := gaps]
    . <- rbind(., dummy)
    plot.order <- function(x) {
      if (x == 'No WMC') 1
      else if (grepl('Gap', x)) as.numeric(gsub('Gap ', '', x))
      else 3 * (as.numeric(substr(x, 0, 1)) - 1) + ifelse(
        grepl('Low', x), 0, as.numeric(gsub('.*High ', '', x))
      )
    }
    .[, WMC := WMC %>% {
      factor(., levels = unique(.)[order(sapply(unique(.), plot.order))])
    }]
    p <- lapply(
      .[['WMC']] %>% {unique(.)[order(sapply(unique(.), plot.order))]},
      function(wmc) {
        if (grepl('Gap', wmc)) {
          ggplot() + theme_void()
        } else {
          .[WMC == wmc] %>% {
            ggplot(., aes(x = UMAP_1, y = UMAP_2, color = `Cell Type`)) +
              geom_point() + theme_bw() +
              ggtitle(wmc) + theme(
                axis.title.x = element_blank(), axis.title.y = element_blank(),
                plot.title = element_text(hjust = 0.5))
          }
        }
      }) %>% ggarrange(plotlist = ., ncol = 4, nrow = 3, common.legend = TRUE,
                       legend = 'right')
    annotate_figure(p, bottom = text_grob('UMAP1'),
                    left = text_grob('UMAP2', rot = 90))
  } %>%
    ggsave(sprintf('UMAP %s.pdf', dir), ., width = 12, height = 9)
}
sapply(dirs.umap, plot.umap)

# Save ARI data frames
dirs.hgc <- dirs %>% .[-c(4, length(.))]
sapply(dirs.hgc, function(dir) {
  wmc.clustering.trees <- readRDS(file.path(dir, 'wmc_clustering_trees.rds'))
  sapply(names(wmc.clustering.trees), function(wmc) {
    wmc.clustering.trees[[wmc]] %>% {
      tree <- .@graphs$ClusteringTree
      tree$height <- log(tree$height + 1)
      subtype <- .@meta.data[1, 'subtype']
      meta.data <- .@meta.data %>% select(celltype_major) %>%
        rename('Cell Type' = 'celltype_major')
      dir.create(file.path(dir, 'DendrogramPlots'), showWarnings = FALSE)
      res <- 300
      length <- 480 / 72 * res
      png(
        file.path(dir, 'DendrogramPlots', paste0(wmc, '.png')),
        res = res, width = length, height = length
      )
      par(mar = c(3, 3, 0, 0))
      HGC.PlotDendrogram(tree, plot.label = TRUE, labels = meta.data)
      dev.off()
      HGC.PlotARIs(tree, labels = meta.data, return.ARI = TRUE) %>%
        data.frame %>% mutate(k = 1:nrow(.) + 1, WMC = wmc, subtype = subtype)
    }
  }, simplify = FALSE) %>% bind_rows %>% rename('ARI' = 'Cell Type') %>%
    write.csv(file.path(dir, 'ARI.csv'), row.names = FALSE)
})

# Combine ARIs and ggplot
sapply(dirs.hgc, function(dir) {
  read.csv(file.path(dir, 'ARI.csv'))
}, simplify = FALSE) %>% bind_rows %>%
  {
  ggplot(., aes(x = k, y = ARI, color = WMC)) + geom_line() + geom_point() +
    ylab('ARI for Cell Type') + facet_wrap(~subtype)
} %>% ggsave(dirs.for.hgc %>% paste(collapse = ' ') %>%
               sprintf('ARI Subtype %s.pdf', .), .)