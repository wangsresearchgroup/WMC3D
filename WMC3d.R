library(dplyr)
library(Seurat)
library(patchwork)
library(wavelets)
library(Matrix)
library(scCATCH)
library(hdf5r)
library(readr)
library(UpSetR)
library(scatterplot3d)
library(plotly)
library(dittoSeq)

# Step 1: Data pre-processing

# input 
# pbmc.data <- Read10X(...) 

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

lb_feature=200
ub_feature=4200
max.percent.mt=5

pbmc <- subset(pbmc, subset = nFeature_RNA > lb_feature & nFeature_RNA < ub_feature & percent.mt < max.percent.mt)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_clean <- pbmc@assays$RNA@data
working_matrix=as.matrix(pbmc_clean)


## Step2: Discrete Wavelet transform
two_band <- function(X){
  n=dim(X)[1]
  nn=2*ceiling(n/2)
  m=dim(X)[2]
  VV=matrix(NA,nn,m)
  WW=matrix(NA,nn,m)
  XX=matrix(0,nn,m)
  XX[1:n,]=X
  
  for (i in 1:m){
    p<-XX[,i]
    B<-dwt(p,'d4',n.level=1) #d8 d12 ... for alternative 
    B1<-B
    B2<-B
    B1@W$W1<-matrix(0,nn/2,1)
    B2@V$V1<-matrix(0,nn/2,1)
    VV[,i]<-idwt(B1)
    WW[,i]<-idwt(B2)
  }
  VV=VV[1:n,]
  WW=WW[1:n,]
  result= list("low" = VV, "high1" = WW)
  return (result)
}

three_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0=c(0.33838609728386,0.53083618701374,0.72328627674361,0.23896417190576,0.04651408217589,-0.14593600755399)
  h1=c(-0.11737701613483,0.54433105395181,-0.01870574735313,-0.69911956479289,-0.13608276348796,0.42695403781698)
  h2=c(0.40363686892892,-0.62853936105471,0.46060475252131,-0.40363686892892,-0.07856742013185,0.24650202866523)
  #n=number of genes
  #m=number of cells
  nn=3*ceiling(n/3)
  
  vv1=matrix(0,nn,m)
  ww1=matrix(0,nn,m)
  ww2=matrix(0,nn,m)
  for (k in 1:m){
    ss=X[,k]
    if(n%%3){ss=c(ss,rep(0,3-n%%3))}
    ss=c(ss,ss[1:3])
    v1=rep(0,nn/3)
    w1=rep(0,nn/3)
    w2=rep(0,nn/3)
    j=1
    for (i in 1:(nn/3)){
      v1[i]=sum(h0*ss[j:(j+5)])
      w1[i]=sum(h1*ss[j:(j+5)])
      w2[i]=sum(h2*ss[j:(j+5)])
      j=j+3
    }
    v1=c(v1[nn/3],v1)
    w1=c(w1[nn/3],w1)
    w2=c(w2[nn/3],w2)
    
    for (i in 1:(nn/3)){
      vv1[3*i-2,k]=h0[4]*v1[i]+h0[1]*v1[i+1]
      vv1[3*i-1,k]=h0[5]*v1[i]+h0[2]*v1[i+1]
      vv1[3*i,k]=h0[6]*v1[i]+h0[3]*v1[i+1]
      ww1[3*i-2,k]=h1[4]*w1[i]+h1[1]*w1[i+1]
      ww1[3*i-1,k]=h1[5]*w1[i]+h1[2]*w1[i+1]
      ww1[3*i,k]=h1[6]*w1[i]+h1[3]*w1[i+1]
      ww2[3*i-2,k]=h2[4]*w2[i]+h2[1]*w2[i+1]
      ww2[3*i-1,k]=h2[5]*w2[i]+h2[2]*w2[i+1]
      ww2[3*i,k]=h2[6]*w2[i]+h2[3]*w2[i+1]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  ww2=ww2[1:n,]
  result= list("low" = vv1, "high1" = ww1,'high2'=ww2)
  return (result)
}


four_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0<-(c(-0.067371764,0.094195111,0.40580489,0.567371764,0.567371764,0.40580489,0.094195111,-0.067371764))
  h1<-(c(-0.094195111,0.067371764, 0.567371764 ,0.40580489,-0.40580489,-0.567371764,-0.067371764,0.094195111))
  h2<-(c(-0.094195111,-0.067371764,0.567371764,-0.40580489,-0.40580489,0.567371764,-0.067371764,-0.094195111))
  h3<-(c(-0.067371764,-0.094195111,0.40580489,-0.567371764,0.567371764,-0.40580489,0.094195111,0.067371764))
  
  nn=4*ceiling(n/4)
  
  vv1=matrix(0,nn,m)
  ww1=matrix(0,nn,m)
  ww2=matrix(0,nn,m)
  ww3=matrix(0,nn,m)
  for (k in 1:m){
    ss=X[,k]
    if(n%%4){ss=c(ss,rep(0,4-n%%4))}
    ss=c(ss,ss[1:4])
    v1=rep(0,nn/4)
    w1=rep(0,nn/4)
    w2=rep(0,nn/4)
    w3=rep(0,nn/4)
    j=1
    for (i in 1:(nn/4)){
      v1[i]=sum(h0*ss[j:(j+7)])
      w1[i]=sum(h1*ss[j:(j+7)])
      w2[i]=sum(h2*ss[j:(j+7)])
      w3[i]=sum(h3*ss[j:(j+7)])
      j=j+4
    }
    v1=c(v1[nn/4],v1)
    w1=c(w1[nn/4],w1)
    w2=c(w2[nn/4],w2)
    w3=c(w3[nn/4],w3)
    for (i in 1:(nn/4)){
      vv1[4*i-3,k]=h0[5]*v1[i]+h0[1]*v1[i+1]
      vv1[4*i-2,k]=h0[6]*v1[i]+h0[2]*v1[i+1]
      vv1[4*i-1,k]=h0[7]*v1[i]+h0[3]*v1[i+1]
      vv1[4*i,k]=h0[8]*v1[i]+h0[4]*v1[i+1]
      ww1[4*i-3,k]=h1[5]*w1[i]+h1[1]*w1[i+1]
      ww1[4*i-2,k]=h1[6]*w1[i]+h1[2]*w1[i+1]
      ww1[4*i-1,k]=h1[7]*w1[i]+h1[3]*w1[i+1]
      ww1[4*i,k]=h1[8]*w1[i]+h1[4]*w1[i+1]
      ww2[4*i-3,k]=h2[5]*w2[i]+h2[1]*w2[i+1]
      ww2[4*i-2,k]=h2[6]*w2[i]+h2[2]*w2[i+1]
      ww2[4*i-1,k]=h2[7]*w2[i]+h2[3]*w2[i+1]
      ww2[4*i,k]=h2[8]*w2[i]+h2[4]*w2[i+1]
      ww3[4*i-3,k]=h3[5]*w3[i]+h3[1]*w3[i+1]
      ww3[4*i-2,k]=h3[6]*w3[i]+h3[2]*w3[i+1]
      ww3[4*i-1,k]=h3[7]*w3[i]+h3[3]*w3[i+1]
      ww3[4*i,k]=h3[8]*w3[i]+h3[4]*w3[i+1]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  ww2=ww2[1:n,]
  ww3=ww3[1:n,]
  
  result= list("low" = vv1, "high1" = ww1,'high2'=ww2,'high3'=ww3)
  return (result)
}

mDWT<-function(X,M){
  if (M==2){
    return(two_band(X))
  }
  if (M==3){
    return(three_band(X))
  }
  if (M==4){
    return(four_band(X))
  }
  warning('invalid band number')
}

bandnumber=4
wavelet_result=mDWT(working_matrix,bandnumber)

#Step 3: dimesion reducing and clustering

dimreduce<-function(data_,ob,X){
  pbmc_vv<-ob
  colnames(data_)<-colnames(X)
  rownames(data_)<-rownames(X)
  vv_sprase<-Matrix(data_, sparse = TRUE)
  pbmc_vv@assays$RNA@data<-vv_sprase
  pbmc_vv <- FindVariableFeatures(pbmc_vv, selection.method = "vst", nfeatures = 5000)
  top10_vv <- head(VariableFeatures(pbmc_vv), 10)
  
  # plot1 <- VariableFeaturePlot(pbmc_vv)
  # plot2 <- LabelPoints(plot = plot1, points = top10_vv, repel = TRUE)
  # plot1 + plot2
  
  all.genes <- rownames(pbmc_vv)
  pbmc_vv <- ScaleData(pbmc_vv, features = all.genes)
  
  pbmc_vv <- RunPCA(pbmc_vv, features = VariableFeatures(object = pbmc_vv))
  
  VizDimLoadings(pbmc_vv, dims = 1:2, reduction = "pca")
  DimPlot(pbmc_vv, reduction = "pca")
  pbmc_vv <- JackStraw(pbmc_vv, num.replicate = 100)
  pbmc_vv <- ScoreJackStraw(pbmc_vv, dims = 1:20)
  
  JackStrawPlot(pbmc_vv, dims = 1:15)
  ElbowPlot(pbmc_vv)
  
  pbmc_vv <- FindNeighbors(pbmc_vv, dims = 1:10)
  pbmc_vv <- FindClusters(pbmc_vv, resolution = 0.5)
  head(Idents(pbmc_vv), 10)
  pbmc_vv <- RunUMAP(pbmc_vv, dims = 1:10)
  
  pbmc_vv.markers <- FindAllMarkers(pbmc_vv, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  result= list("data" = pbmc_vv, "markers" = pbmc_vv.markers)
  return(result)
}

reduce_low<-dimreduce(wavelet_result$low,pbmc,working_matrix)
reduce_h1<-dimreduce(wavelet_result$high1,pbmc,working_matrix)
reduce_h2<-dimreduce(wavelet_result$high2,pbmc,working_matrix)
reduce_h3<-dimreduce(wavelet_result$high3,pbmc,working_matrix)

#Step 4: Multi-view plotting and annotation

withoutdwt<-dimreduce(working_matrix,pbmc,working_matrix)
DimPlot(withoutdwt$data, reduction = "umap")

data_new=reduce_low$data@reductions$umap@cell.embeddings
z1=as.numeric(reduce_h1$data$seurat_clusters)
z2=as.numeric(reduce_h2$data$seurat_clusters)
z3=as.numeric(reduce_h3$data$seurat_clusters)
clu=reduce_low$data$seurat_clusters
data_new=cbind(data_new,clu,z1,z2,z3)
data_new=as.data.frame(data_new)
colnames(data_new)[3:6]=c('cluster','high1','high2','high3')
plot_ly(data=data_new,x=~UMAP_1,y=~UMAP_2,z=~high3,color=~cluster,colors=dittoColors(),
        type='scatter3d',mode='markers',marker=list(size=2,width=1))


#cell identification
cell_identity<-function(data,tissue_name,match_CellMatch_,cancer){
  clu_markers <- findmarkergenes(object = data,species = 'Human',cluster = 'All',match_CellMatch = match_CellMatch_,
                                 cancer = cancer,
                                 tissue = tissue_name,
                                 cell_min_pct = 0.1,
                                 logfc = 0.1,
                                 pvalue = 0.01)
  
  clu_ann <- scCATCH(object = clu_markers$clu_markers,
                     species = 'Human',
                     cancer = cancer,
                     tissue = tissue_name)
  
  result= list("clu_markers " = clu_markers , "indenity" =clu_ann )
  return (result)
}

cell_tppe_name<-function(d1,d2){
  name=d1$indenity$cell_type
  d3=d2$data@meta.data$seurat_clusters
  len_=length(name)
  for(i in 1:len_){
    if (is.na(name[i])){
      name[i]='Unknown'
    }
  }
  a=as.numeric(unique(d3))-1
  new_name=1:max(a)
  b=as.numeric(d1$indenity$cluster)
  c=setdiff(a,b)+1
  
  if (length(c)!=0){
  for (j in c){
    new_name[j]='Unknown'
  }
  
  k=1
  for (j in (b+1)){
    new_name[j]=name[k]
    k=k+1
  }
  return (new_name)
  }
  else{
    return (name)
  }
}

annotation<-function(reduced,tissu_name,mtch_CellMatch_name,cancer_name){
  cellid<-cell_identity(reduced$data,tissu_name,mtch_CellMatch_name,cancer_name)
  return(cell_tppe_name(cellid,reduced))
}

tissu_name<-'Breast'
mtch_CellMatch_name<-FALSE
cancer_name<-NULL

annotation(reduce_low,tissu_name,mtch_CellMatch_name,cancer_name)
