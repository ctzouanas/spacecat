getwd()

out_dir = paste0(getwd(),'/output/')
if (!file.exists(out_dir)){dir.create(out_dir)}

#source some functions Genshaft likes
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/Interacting Cells/T_MDDC_RNAseq/HIV Analysis/LoveHIVfunctions2.R')
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/Interacting Cells/T_MDDC_RNAseq/PMA Study/Pairs_PMA_Analysis/LovePMApipeFunctions.R')
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/Interacting Cells/T_MDDC_RNAseq/All_Pairs_Data/All_Pairs_Reannotate_Data.R')
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/Interacting Cells/T_MDDC_RNAseq/HIV Analysis/MAST_wrapper_function.R')
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab/Projects/Interacting Cells (Pairs)/Rita & Britt Interacting APC.Tcell/Functions_ASG.R')
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/General Code/ChooseClusterResolutionDownsample_Seurat_v4.R')


load(file = '~/Dropbox (MIT)/shalek lab top of tree/Shalek Lab/Projects/SpaceCat/Experiments/191021_Tumor_1/Analysis/Deep_both/20200222_LungKPTumor.RData')
ls()

#the original seurat object from Carly
data.all.seurat

#convert the seurat object to the most up-to-date seurat object class w accessible functions, etc
data.all.seurat.update = UpdateSeuratObject(data.all.seurat)
data.all.seurat.update = SetIdent(data.all.seurat.update, value = data.all.seurat.update$Clusters.Names)
table(data.all.seurat.update$Sort,data.all.seurat.update$Tumor)
table(data.all.seurat.update$Clusters.Names,data.all.seurat.update$Sort)
# get some wes anderson colors going.
require(wesanderson)
wes.colors = c(wes_palette("Zissou1", n=25, type="continuous"), wes_palette("Darjeeling2", n=5, type="discrete"), wes_palette("GrandBudapest1", n=4, type="discrete"), wes_palette("GrandBudapest2", n=4, type="discrete"))

# remove Core, remove doublet/low quality cluster
data.new.all.seurat = subset(data.all.seurat.update, 
                         idents = names(table(data.all.seurat.update@active.ident))[c(1:10,12:21)], 
                         subset = `Sort` %in% c('Outer','Rest'))
data.new.all.seurat
table(data.new.all.seurat@active.ident)
# do some house keeping to rrename things
data.new.all.seurat@meta.data$Clusters.Names = as.character(data.new.all.seurat@meta.data$Clusters.Names)
data.new.all.seurat@meta.data$Clusters.Names[which(data.new.all.seurat@meta.data$Clusters.Names == "DCs Ccr7 high Ccl22 high")] = "Ccr7+ Ccl22+ DCs"
data.new.all.seurat@meta.data$Clusters.Names[which(data.new.all.seurat@meta.data$Clusters.Names == "Areg+ DCs")] = "Areg+ Macrophages"
data.new.all.seurat@meta.data$Clusters.Names[which(data.new.all.seurat@meta.data$Clusters.Names == "Macrophages")] = "Monocytes or Macrophages"
data.new.all.seurat@meta.data$Clusters.Names = as.factor(data.new.all.seurat@meta.data$Clusters.Names)

data.new.all.seurat@meta.data$Tumor = as.character(data.new.all.seurat@meta.data$Tumor)
data.new.all.seurat@meta.data$Tumor[which(data.new.all.seurat@meta.data$Tumor == "Tumor1")] = "Tumor 2"
data.new.all.seurat@meta.data$Tumor[which(data.new.all.seurat@meta.data$Tumor == "Tumor2")] = "Tumor 1"
data.new.all.seurat@meta.data$Tumor = as.factor(data.new.all.seurat@meta.data$Tumor)

data.new.all.seurat@meta.data$Sort = as.character(data.new.all.seurat@meta.data$Sort)
data.new.all.seurat@meta.data$Sort[which(data.new.all.seurat@meta.data$Sort == "Outer")] = "Healthy Tumor Border"
data.new.all.seurat@meta.data$Sort[which(data.new.all.seurat@meta.data$Sort == "Rest")] = "Whole Tumor"
data.new.all.seurat@meta.data$Sort = as.factor(data.new.all.seurat@meta.data$Sort)

table(data.new.all.seurat@active.ident)
table(data.new.all.seurat$Sort,data.new.all.seurat$Tumor)
table(data.new.all.seurat$Tumor)
table(data.new.all.seurat$Clusters.Names)

data.new.all.seurat = SetIdent(data.new.all.seurat, value = data.new.all.seurat$Clusters.Names)

#rerun generic pipeline on this subset of data
sub_dir = paste0(out_dir,'All_Data_no_core_no_low_quality/')
if (!file.exists(sub_dir)){dir.create(sub_dir)}

data.new.all.seurat = FindVariableFeatures(data.new.all.seurat, nfeatures = 2000)
data.new.all.seurat = ScaleData(object = data.new.all.seurat, do.scale=TRUE, do.center=TRUE)
range(data.new.all.seurat@assays$RNA@scale.data)

g = VariableFeaturePlot(data.new.all.seurat, assay = 'RNA', log = 2)
g
pdf(paste0(sub_dir,"MeanVarPlot_nonames.pdf"), useDingbats = FALSE)
g
dev.off()

top <- head(x = VariableFeatures(object = data.new.all.seurat), 
            n = 40) # 40 appeared to be the right number, but feel free to toggle

# Plot variable features with labels
g = LabelPoints(plot = g, 
                points = top, 
                repel = TRUE)
g

pdf(paste0(sub_dir,"MeanVarPlot_names.pdf"), useDingbats = FALSE)
g
dev.off()

#run pcs over variable genes
#pcs.range = 1:100
data.new.all.seurat = RunPCA(data.new.all.seurat, features = NULL, npcs = 100) #run on variable genes, 100 PCs


ElbowPlot(data.new.all.seurat, ndims=100) + geom_vline(xintercept = 25, col = 'red') +geom_vline(xintercept = 25, col = 'red')
pdf(paste0(sub_dir,"PCjackstraw_elbow.pdf"), useDingbats = FALSE)
ElbowPlot(data.new.all.seurat, ndims=100) + geom_vline(xintercept = 25, col = 'red') +geom_vline(xintercept = 25, col = 'red')
dev.off()

top = c(1:3)
bottom = c(100)
sets.of = 9

t = list()
t$first = min(top,bottom)
t$last = max(top,bottom)
t$start.of.sets = seq(from = t$first, to = t$last, by = sets.of)
t$start.of.sets

g = lapply(t$start.of.sets , function(x){
  if(x+sets.of>dim(data.new.all.seurat@reductions$pca)[2]){
    temp = dim(data.new.all.seurat@reductions$pca)[2]
  } else{
    temp = x + sets.of - 1 # adjustments to final heatmap to plot, should do 3x3 grid
  }
  g = DimHeatmap(data.new.all.seurat, dims = c(x:temp), cells = 500, balanced = TRUE)
  return(g)
})

print(g[[1]])

pdf(paste0(sub_dir,"PCHeatmap_",x,"-",temp,".pdf"), useDingbats = FALSE)
g
dev.off()




sig.pcs = c(1:25)

#UMAP etc to get new figure 4B
set.seed(123)
data.new.all.seurat <- FindNeighbors(data.new.all.seurat, reduction = "pca", dims = sig.pcs, nn.eps = 0.5)
data.new.all.seurat <- RunUMAP(data.new.all.seurat, dims = sig.pcs, min.dist = 0.5, n.neighbors = 30, seed.use = 123)

install.packages("wesanderson")
require(wesanderson)

set.seed(123)
pal = sample(c(wes_palette('Darjeeling1', type = 'continuous', n = 10),
               wes_palette("IsleofDogs1", type = 'continuous', n = 10)),
             size = 20, replace = FALSE)


g = DimPlot(data.new.all.seurat, reduction = 'umap', pt.size=.6, cols = pal, label=T, group.by="ident", repel = T) + 
  theme(legend.position = "none", 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

g  

pdf(paste0(sub_dir,"UMAP_All_Cells_Figure_4B.pdf"), useDingbats = FALSE)
g
dev.off()

pdf(paste0(sub_dir,"UMAP_All_Cells_Figure_4B_legend.pdf"), useDingbats = FALSE)
g+theme(legend.position = 'right')
dev.off()

data.new.all.seurat = SetIdent(data.new.all.seurat, value = data.new.all.seurat$Clusters.Names)
data.new.all.seurat@active.ident = factor(data.new.all.seurat@active.ident, levels = levels(data.new.all.seurat@active.ident)[c(5,7,1,20,4,2,11,6,13,18,15,8,19,3,14,10,17,16,9,12)])
data.subset@ident = temp.3

#rescale to make sure you get all the genes you want
fig.4c.genes = c("H2-Ab1", "Cd74", "Ccl17", "Flt3", "Ccr7", "Ccl22","Areg", "Pdpn", "Mmp12", "Dcstamp","Vdr", "Tnfrsf9", "C1qa", "C1qb", "Apoe", "C4b", "Marco", "Arg1","Chi3l3", "Lpl", "Il18",  "Vcan", "Plac8", "Ccr2","S100a9", "Csf3r", "G0s2",  "Rorc", "Foxp3", "Il23r", "Ifngr1", "Lef1", "S1pr1", "Nkg7", "Gzma", "Prf1", "Il5", "Il1rl1", "Gata3", "Pou2af1", "Igkc", "Cd19", "Gm10722", "Gm10717","Gm21738", "Krt18", "Clu", "Krt19", "Wnt7b", "Itga2","Sox4",  "Sftpb", "Ager", "Sftpc", "Col3a1", "Dcn", "Mgp","Pecam1", "Podxl", "Lyve1")
data.new.all.seurat = ScaleData(object = data.new.all.seurat, features = fig.4c.genes, do.scale=TRUE, do.center=TRUE)

saveRDS(data.new.all.seurat, file = '20210324_Revision_Tumor_Seurat_obj.rds')

g = DotPlot(data.new.all.seurat, features = rev(fig.4c.genes), 
            cols =c("lightgrey", "black"), col.min=0, dot.min=.01) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 2, hjust = 0),
        legend.position="right") + 
  scale_x_discrete(position = 'top')
  

pdf(paste0(sub_dir,"DotPlot Final Subcluster Both Tumors 3 genes per cluster Figure 4C.pdf"), width=16, height=5, useDingbats=FALSE)
g
dev.off()

Tumor.1.Colors = wes_palette("Cavalcanti1")[c(1,2,4,5)]
g = DimPlot(data.new.all.seurat, 
            cells = which(data.new.all.seurat$Tumor == 'Tumor 1'),
            reduction = 'umap', pt.size=.6, cols = Tumor.1.Colors, group.by="Array_Sort") + 
  theme(legend.position = "none", 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

g  

pdf(paste0(sub_dir,"UMAP_All_Cells_Figure_4D_tumor_1.pdf"), useDingbats = FALSE)
g
dev.off()

pdf(paste0(sub_dir,"UMAP_All_Cells_Figure_4D_tumor_1_legend.pdf"), useDingbats = FALSE)
g+theme(legend.position = 'right')
dev.off()

data.new.all.seurat = data.all.seurat
wes_palette("Darjeeling2")[5]
Tumor.2.Colors = c(wes_palette("Zissou1")[c(1,4,5)],wes_palette("Darjeeling2")[5])
g = DimPlot(data.new.all.seurat, 
            cells = which(data.new.all.seurat$Tumor == 'Tumor 2'),
            reduction = 'umap', pt.size=.6, cols = Tumor.2.Colors, group.by="Array_Sort") + 
  theme(legend.position = "none", 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

g  

pdf(paste0(sub_dir,"UMAP_All_Cells_Figure_4D_tumor_2.pdf"), useDingbats = FALSE)
g
dev.off()

pdf(paste0(sub_dir,"UMAP_All_Cells_Figure_4D_tumor_2_legend.pdf"), useDingbats = FALSE)
g+theme(legend.position = 'right')
dev.off()

######### done with the easy TSNE plotting, time for the stacked barplots
head(data.new.all.seurat@meta.data)
table(data.new.all.seurat$Clusters.Names,data.new.all.seurat$Array_Sort,data.new.all.seurat$Tumor)

#calculate order for the clusters by looking at overall, whole tumor frequency of each cluster
# this is cosmetic so it doesnt have to be 'perfect' but we want the smallest values to the 
# lowest part of hte graph because we will log the y axis
overall.sort.df = table(data.new.all.seurat@active.ident,data.new.all.seurat$Sort)
overall.sort.df = t(t(overall.sort.df)/rowSums(t(overall.sort.df)))

#cluster.order = names(sort(overall.sort.df[,'Whole Tumor'], decreasing = FALSE))

cluster.order = levels(data.new.all.seurat@active.ident)

overall.array.df = table(data.new.all.seurat@active.ident,data.new.all.seurat$Array)
overall.array.df = t(t(overall.array.df)/rowSums(t(overall.array.df)))[,-8]

head(melt(tumor.1.df.freq))

data.new.all.seurat = AddMetaData(data.new.all.seurat,
                                  metadata = paste(data.new.all.seurat$Sort,data.new.all.seurat$Tumor,sep='_'),
                                  col.name = 'Sort_Tumor')

sample.org.df = table(data.new.all.seurat$Array,data.new.all.seurat$Sort_Tumor)
sample.org.df[which(sample.org.df>0)]=1


cluster.freq.array.list=list()
for (i in colnames(sample.org.df)){
  cluster.freq.array.list[[i]] = list()
  temp = row.names(sample.org.df)[which(sample.org.df[,i]==1)]
  for(j in cluster.order){
    cluster.freq.array.list[[i]][[j]] = overall.array.df[j,temp]
    names(cluster.freq.array.list[[i]][[j]]) = temp
  }
}


cluster.freq.array.list.summary = lapply(cluster.freq.array.list , function(x){
  return(lapply(x,function(y){
    temp = c(mean(y))
    names(temp)=c('mean')
    return(t(as.matrix(temp)))
  }))
})
#fna.freq.list.summary

require(reshape2)
cluster.freq.df.summary = reshape2::melt(cluster.freq.array.list.summary)
colnames(cluster.freq.df.summary)[3] = 'mean' 

cluster.freq.list.summary = lapply(cluster.freq.array.list , function(x){
  return(lapply(x,function(y){
    temp = c(sd(y))
    names(temp)=c('sd')
    return(t(as.matrix(temp)))
  }))
})
#fna.freq.df.summary
#reshape2::melt(fna.freq.list.summary,)
cluster.freq.df.summary = cbind(cluster.freq.df.summary, sd = reshape2::melt(cluster.freq.list.summary)$value)

cluster.freq.df.summary
names(cluster.freq.array.list)
cluster.freq.df.summary$adjusted_mean = 0
#this for loop only works because they were already in the correct order.....
for (i in names(cluster.freq.array.list)){
  cluster.freq.df.summary[which(cluster.freq.df.summary$L1 == i),'adjusted_mean'] = cumsum(cluster.freq.df.summary[which(cluster.freq.df.summary$L1 == i),'mean'])
}
cluster.freq.df.summary$L2 = factor(cluster.freq.df.summary$L2, levels = rev(cluster.order))
#you may want to set the order of these groups in a certain way... right now, lets do tumor 1 first then tumor 2
cluster.freq.df.summary$L1 = factor(cluster.freq.df.summary$L1, levels = c('Whole Tumor_Tumor 1','Healthy Tumor Border_Tumor 1',
                                                                           'Whole Tumor_Tumor 2','Healthy Tumor Border_Tumor 2'))
Clusters.Colors.2 = pal
names(Clusters.Colors.2) = sort(levels(data.new.all.seurat@active.ident))

cluster.freq.df.summary$colors.chosen = factor(Clusters.Colors.2[cluster.freq.df.summary$L2])

cluster.freq.df.summary$L2


g = ggplot(cluster.freq.df.summary, aes(x = L1, y = mean, fill = L2)) + 
  scale_fill_manual(values = Clusters.Colors.2[cluster.order]) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymax = adjusted_mean + sd/3, ymin=adjusted_mean - sd/3), 
                position = "identity", width = .4, alpha = 1, size = .1) + 
  scale_color_manual(values = Clusters.Colors.2) + 
  theme(axis.text.x = element_text(angle = -30, hjust = 0))

g = make_ggplot_pretty(g, ylab.custom = 'Frequency of Cells', xlab.custom = 'Sample')

g

pdf(paste0(sub_dir,"UMAP_All_Cells_Figure_4E_both.pdf"), useDingbats = FALSE)
g
dev.off()

pdf(paste0(sub_dir,"UMAP_All_Cells_Figure_4E_both_no_legend.pdf"), useDingbats = FALSE)
g + theme(legend.position = 'none')
dev.off()

cluster.freq.df.summary$tumor = unlist(lapply(strsplit(as.character(cluster.freq.df.summary$L1),
                                                       split = '_', fixed = TRUE),
                                              function(x){
                                                return(x[2])
                                              }))

# lets do this as a side by side bar plot instead
library(scales)
g = ggplot(cluster.freq.df.summary, aes(y = (mean+0.0001)*10000, x = L2, fill = L1)) + 
  scale_fill_manual(values = wes_palette("BottleRocket2")[c(1,4,2,3)]) +
  geom_bar(stat = "identity", position = 'dodge') + 
  geom_errorbar(aes(ymax = (mean+sd/3+0.0001)*10000, ymin= (mean-sd/3+0.0001)*10000), 
                position = "dodge",  alpha = 1, size = .1) + 
  scale_color_manual(values = Clusters.Colors.2) + 
  #coord_trans(y = 'log10', ylim = c(0.0001,1)) +
  scale_y_log10(lim = c(1,10000), 
                breaks = c(1,10,100,1000,10000), 
                labels = comma(c(0.0001,0.001,0.01,0.1,1),accuracy = 0.0001))+
  facet_wrap(~tumor) + coord_flip() + labs(y = 'Cell Type Frequency', x = 'Cell Type') +
  theme_bw() +
  theme(panel.spacing.x = unit(1,'lines'), legend.title = element_blank()) 
  #ylim() 

pdf(paste0(sub_dir,"cell_frequency_bar_plot.pdf"), useDingbats = FALSE, width = 10, height = 5)
g
dev.off()

#g+theme_bw()

min(cluster.freq.df.summary$mean[cluster.freq.df.summary$mean!=0])

##### do fisher exact test on these compositions to determine significance

table(data.new.all.seurat$Clusters.Names,data.new.all.seurat$Sort_Tumor)
table(data.new.all.seurat$Sort_Tumor, data.new.all.seurat$Array)

table(data.new.all.seurat$Clusters.Names,data.new.all.seurat$Sort_Tumor)
### loop through to do both tumors
tumor.vec = c('Tumor 1', 'Tumor 2')
cell.freq.confusion.mat = list()
t.mat = matrix(data = NA, nrow = 2, ncol = 2)
for(i in 1:length(tumor.vec)){
  t.seurat = subset(data.new.all.seurat, subset = Tumor == tumor.vec[i])
  t.table.overall = table(t.seurat$Sort)
  t.table.cluster = table(t.seurat$Clusters.Names,t.seurat$Sort)
  cell.freq.confusion.mat[[i]] = list()
  for (j in 1:dim(t.table.cluster)[1]){
    cell.freq.confusion.mat[[i]][[j]] = t.mat
    cell.freq.confusion.mat[[i]][[j]][1,] = c(t.table.cluster[j,'Healthy Tumor Border'], t.table.overall['Healthy Tumor Border']-t.table.cluster[j,'Healthy Tumor Border'])
    cell.freq.confusion.mat[[i]][[j]][2,] = c(t.table.cluster[j,'Whole Tumor'], t.table.overall['Whole Tumor']-t.table.cluster[j,'Whole Tumor'])
  }
  names(cell.freq.confusion.mat[[i]]) = row.names(t.table.cluster)
}
names(cell.freq.confusion.mat) = tumor.vec

cell.freq.fisher.list = lapply(cell.freq.confusion.mat, function(x){
  lapply(x, function(y){
    return(fisher.test(y)$p.value)})})

cell.freq.fisher.melt = melt(cell.freq.fisher.list)
cell.freq.fisher.melt$p_val_adj = p.adjust(cell.freq.fisher.melt$value,method = 'BH')

cell.freq.fisher.melt.sig = subset(cell.freq.fisher.melt, subset = p_val_adj < 0.001)

table(cell.freq.fisher.melt.sig$L1,cell.freq.fisher.melt.sig$L2)

table(data.new.all.seurat$Tumor)

#################################################
########                                #########
########  Lets do Mono/Mac Analysis     #########
########       adapted from CGKZ        #########
########                                #########
#################################################

### note - when CGKZ did analysis, tumor 2 was tumor 1 and vice versa. 
###   thus their analysis did not align w the figure/paper annotation
###   this analysis does match the figure/paper annotation
### thus here we will focus on tumor 1, where CGKZ focused on tumor 2. they are the same tumor - samples 5E-9I
###    one other change is that the tumor core sample 8H has been removed from this analysis
###    that will change the heatmap (fewer cells)
###     it may also effect the genes revealed....... as certain cells will be removed from the correlation analyses 
###                                                     (both null model & signal)

##### start w tumor 2
data.subset = subset(data.new.all.seurat, cells = intersect(which(data.new.all.seurat$Tumor == 'Tumor 2'),
                                                            which(data.new.all.seurat$Clusters.Names == 'Monocytes or Macrophages')))
table(data.subset$Clusters.Names_Tumor)
table(data.subset$Array)

data.subset = SetIdent(data.subset, value = as.vector(data.subset@meta.data$Sort))
table(data.subset@active.ident)

out.de = FindMarkers(data.subset, ident.1 = "Whole Tumor", ident.2 = "Healthy Tumor Border", logfc.threshold=.01, only.pos=FALSE, min.pct=.01, test.use="bimod")

#### plot volcano plot
range(out.de$avg_log2FC)
range(-log10(out.de$p_val_adj))
out.de$color.set = 'grey'
out.de$color.set[intersect(which(out.de$p_val_adj < 0.05),
                           which(out.de$avg_log2FC < 0))] = 'black'
out.de$color.set[intersect(which(out.de$p_val_adj < 0.05),
                           which(out.de$avg_log2FC > 0))] = 'dodgerblue2'
out.de$gene = row.names(out.de)
require(ggrepel)
g = ggplot(out.de, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color.set)) +
  geom_hline(aes(yintercept=-log10(.05)), linetype="dashed", color="grey") +
  geom_vline(aes(xintercept =0), linetype="dashed", color="grey") + 
  xlim(-2.2,2.2) + ylim(0,12) + 
  geom_point(size = .1) +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.05,as.character(gene),'')),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50', size = 4, force = 1, max.overlaps = 10) +
  scale_color_manual(values=c("black", "firebrick2", "grey")) + 
  theme_bw() + 
  theme(legend.position = 'none') 

pdf(paste0(sub_dir,"mono_mac_de_whole_vs_border_Tumor_2.pdf"), width=7, height=6, useDingbats=FALSE)
g
dev.off()

#test for il1a
g = ggplot(out.de, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color.set)) +
  geom_hline(aes(yintercept=-log10(.05)), linetype="dashed", color="grey") +
  geom_vline(aes(xintercept =0), linetype="dashed", color="grey") + 
  xlim(-2.2,2.2) + ylim(0,12) + 
  geom_point(size = .1) +
  geom_text_repel(data=subset(out.de, gene == 'Il1a'),
                  aes(label = gene),
                  box.padding   = 0.35, 
                  point.padding = 0.5,
                  nudge_y = 10+log10(subset(out.de, gene == 'Il1a')$p_val_adj),
                  segment.color = 'grey50', size = 4, force = 1, max.overlaps = 10) +
  scale_color_manual(values=c("black", "firebrick2", "grey")) + 
  theme_bw() + 
  theme(legend.position = 'none') 

##### now tumor 1

data.subset = subset(data.new.all.seurat, cells = intersect(which(data.new.all.seurat$Tumor == 'Tumor 1'),
                                                            which(data.new.all.seurat$Clusters.Names == 'Monocytes or Macrophages')))
table(data.subset$Clusters.Names_Tumor)
table(data.subset$Array)

data.subset = SetIdent(data.subset, value = as.vector(data.subset@meta.data$Sort))
table(data.subset@active.ident)

out.de = FindMarkers(data.subset, ident.1 = "Whole Tumor", ident.2 = "Healthy Tumor Border", logfc.threshold=.01, only.pos=FALSE, min.pct=.01, test.use="bimod")
sig.outer = out.de[which(out.de$avg_logFC<0),]
sig.outer = sig.outer[which(sig.outer $p_val_adj<0.05),]
dim(sig.outer)

#### plot volcano plot
out.de$color.set = 'grey'
out.de$color.set[intersect(which(out.de$p_val_adj < 0.05),
                           which(out.de$avg_log2FC < 0))] = 'black'
out.de$color.set[intersect(which(out.de$p_val_adj < 0.05),
                           which(out.de$avg_log2FC > 0))] = 'dodgerblue2'
out.de$gene = row.names(out.de)
require(ggrepel)
g = ggplot(out.de, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color.set)) +
  geom_hline(aes(yintercept=-log10(.05)), linetype="dashed", color="grey") +
  geom_vline(aes(xintercept =0), linetype="dashed", color="grey") + 
  xlim(-3,3) + ylim(0,110) + 
  geom_point(size = .1) +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.05,as.character(gene),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', size = 4) +
  scale_color_manual(values=c("black", "dodgerblue2", "grey")) + 
  theme_bw() + 
  theme(legend.position = 'none') 

pdf(paste0(sub_dir,"mono_mac_de_whole_vs_border_Tumor_1.pdf"), width=7, height=6, useDingbats=FALSE)
g
dev.off()

# im pretty sure we got the same results based on the number of genes and how this should be an exact replication... but 

#read in previously performed DE... just to be sure.
sig.outer.read = read.table(file=paste0(getwd(),'/../Deep_both/Tumor2_Macrophages/Sig DE genes by bimod Outer up Tumor2_Macrophages .txt'),
                            sep="\t")

length(which(row.names(sig.outer) == row.names(sig.outer.read)))
# yup full match, good.

# proceed to the heatmap

sig.rest = out.de[which(out.de$avg_logFC>0),]
sig.rest = sig.rest[which(sig.rest$p_val_adj<0.05),]
dim(sig.rest)

# take the tumor 2 macrophages, make a "score" of the top DE genes, and rank cells by their expression of each 
dim(sig.outer)
dim(sig.rest)

data.subset = ScaleData(data.subset,features = c(rownames(sig.outer),rownames(sig.rest)))

score.outer = colSums(data.subset@assays$RNA@scale.data[rownames(sig.outer),])
score.rest = colSums(data.subset@assays$RNA@scale.data[rownames(sig.rest),])


score.outer

plot(score.outer, score.rest, pch = 20, cex = 0.3)
model = lm(score.outer ~ score.rest)
abline(model, col = 'red')


plot(rank(score.outer), rank(score.rest), pch = 20, cex = 0.3)
model = lm(rank(score.outer) ~ rank(score.rest))
abline(model, col = 'red')

plot(y= 1:30,x = 1:30, pch = 1:30)


# now want to order cells as most score outer high vs most score rest high 


# start just with the ranked score.outer

plot(rank(score.outer), score.outer)



data.subset = AddMetaData(data.subset, metadata = score.outer, col.name = 'Outer_Score')
data.subset = AddMetaData(data.subset, metadata = rank(score.outer), col.name = 'Outer_Rank')
                          
head(data.subset@meta.data)


VlnPlot(data.subset, c("Outer_Score"))


# do a heatmap over the score.outer genes wiht the columns ranked by the rankscore.outer

temp = data.subset@assays$RNA@scale.data[rownames(sig.outer),order(score.outer)]

head(temp)

dim(temp)

temp = as.matrix(temp)
require(gplots)
heatmap.2(temp, Rowv=NA, Colv=NA, 
          scale="none", density="none", trace="none", 
          col=colorRampPalette(c("dodgerblue3","dodgerblue3", "white", "firebrick3", "firebrick3"))(100))


# add the color map on the colside colors for the location
table(data.subset@meta.data$Sort)
colormap.temp = c("black", "dodgerblue2")[as.integer(factor(data.subset$Sort))]

colormap.temp = colormap.temp[order(score.outer)]




pdf(paste("Heatmap Score on Significant Genes Outer Plot Sig Genes Outer vs Inner colbySort", ParentCluster, ".pdf"), width=8, height=6, useDingbats=FALSE)
heatmap.2(temp, Rowv=NA, Colv=NA, scale="row", density="none", trace="none", col=colorRampPalette(c("dodgerblue3","dodgerblue3","dodgerblue3", "white", "firebrick3", "firebrick3", "firebrick3"))(100), ColSideColors = colormap.temp, labCol=NA)
dev.off()




plot(density(score.outer[which(colormap.temp=="dodgerblue2")]), col="dodgerblue2", xlim=c(-10,25))
lines(density(score.outer[which(colormap.temp=="black")]), col="black")



# then do correlation across this
head(data.subset@meta.data)
dim(as.matrix(data.subset@data))
length(as.vector(data.subset$Outer_Score))


cor.out = cor(t(as.matrix(data.subset@assays$RNA@data)), as.vector(data.subset$Outer_Score), method = 'pearson')
dim(cor.out)



cor.out = cbind(cor.out, cor.out)
head(cor.out)

cor.out = cor.out[which(cor.out[,1]!="NA"),]


cor.out = cor.out[order(cor.out[,1]),]
head(cor.out)
tail(cor.out)




# get the significant background, cut on this

null.mat = sapply(1:10000,function(x){
  return(as.vector(data.subset$Outer_Score)[sample(1:dim(data.subset)[2], dim(data.subset)[2])])
})
dim(null.mat)

cor.out.null.mat = cor(t(as.matrix(data.subset$RNA@data)), null.mat)
dim(cor.out.null.mat)

null.mean = mean(cor.out.null.mat, na.rm = TRUE)
null.sd = sd(cor.out.null.mat, na.rm = TRUE)
null.mean
null.sd

cor.out.outer = cor.out[which(cor.out[,1]<null.sd*-5),]
cor.out.rest = cor.out[which(cor.out[,1]>null.sd*5),]

head(cor.out.outer)
head(cor.out.rest)

data.subset = ScaleData(object = data.subset, 
                        features = c(rownames(cor.out.outer), rownames(cor.out.rest)), 
                        do.scale=TRUE, do.center=TRUE)

temp = data.subset$RNA@scale.data[rev(c(rownames(cor.out.outer), rownames(cor.out.rest))),order(score.outer)]


head(temp)

dim(temp)

temp = as.matrix(temp)


# add the color map on the colside colors for the location
table(data.subset@meta.data$Sort)
colormap.temp = c("black", "dodgerblue2")[as.integer(factor(data.subset@meta.data$Sort))]

colormap.temp = colormap.temp[order(score.outer)]

pdf(paste0(sub_dir,"Heatmap Score on Significant Genes Outer vs Inner Plot Sig Genes Correlation up and down colbySort tumor 1 monocyte and macrophage.pdf"), width=10, height=7, useDingbats=FALSE)
heatmap.2(temp, Rowv=NA, Colv=NA, scale="row", density="none", trace="none", col=colorRampPalette(c("dodgerblue4","dodgerblue4","dodgerblue3","dodgerblue3", "white", "firebrick3", "firebrick3", "firebrick4", "firebrick4"))(100), ColSideColors = colormap.temp, labCol=NA)
dev.off()


pdf(paste0(sub_dir,"Heatmap Score on Significant Genes Outer vs Inner  Plot Sig Genes Correlation up and down colbySort lowerer thresh tumor 1 monocyte and macrophage.pdf"), width=10, height=7, useDingbats=FALSE)
heatmap.2(temp, Rowv=NA, Colv=NA, scale="row", density="none", trace="none", col=colorRampPalette(c("dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue3", "white", "firebrick3", "firebrick3","firebrick3","firebrick3", "firebrick4", "firebrick4", "firebrick4", "firebrick4"))(100), ColSideColors = colormap.temp, labCol=NA)
dev.off()

pdf(paste0(sub_dir,"Heatmap Score on Significant Genes Outer vs Inner  Plot Sig Genes Correlation up and down colbySort lowererer thresh tumor 1 monocyte and macrophage.pdf"), width=10, height=7, useDingbats=FALSE)
heatmap.2(temp, Rowv=NA, Colv=NA, scale="row", density="none", trace="none", col=colorRampPalette(c("dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue3","dodgerblue2", "dodgerblue2", "white", "firebrick2", "firebrick2", "firebrick3", "firebrick3","firebrick3","firebrick3", "firebrick4", "firebrick4", "firebrick4", "firebrick4"))(100), ColSideColors = colormap.temp, labCol=NA)
dev.off()

head(row.names(temp))

#install.packages(ComplexHeatmap)
BiocManager::install("ComplexHeatmap")
require(ComplexHeatmap)

genes.of.interest = c('Ccl3','Cxcl2','Cd274','Tnf','Ccl6','Ccl4','Ccr1','F10',
                      'Clec4e','Cxcl3','F3','Tnfrsf1b','Havcr2','Il23a',
                      'Hspa1a','Hspd1','HspH1','Idh2','Siglec5','Idh3g','Il1r2',
                      'Mdh2','Mdh1','Psmd2','Psmd3','Mtss1','Lcp1','Actn1','Xbp1')
row.names.temp = row.names(temp)
row.names(temp) = NA
colnames(temp) = NA
ha  = rowAnnotation(foo = anno_mark(at = which(row.names(temp) %in% genes.of.interest), 
                                    labels = row.names(temp)[which(row.names(temp) %in% genes.of.interest)]))
h = Heatmap(temp,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        right_annotation = ha,
        use_raster = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE)
pdf(paste0(sub_dir,"Heatmap Score on Significant Genes Outer vs Inner  Plot Sig Genes Correlation up and down colbySort lowerer thresh tumor 1 monocyte and macrophage_for_annotations.pdf"), width=10, height=7, useDingbats=FALSE)
h
dev.off()


write.table(cor.out.outer, file=paste0(sub_dir,"Genes Negatively Correlated with Healthy Tumor Border Score outer vs inner.txt"), quote=FALSE, sep="\t")
write.table(row.names(cor.out.outer), file=paste0(sub_dir,"Genes Negatively Correlated with Healthy Tumor Border Score outer vs inner_just_genes.txt"), quote=FALSE, sep="\t",row.names = FALSE, col.names = FALSE)

write.table(cor.out.rest, file=paste0(sub_dir,"Genes Positively Correlated with Healthy Tumor Border Score outer vs inner.txt"), quote=FALSE, sep="\t")
write.table(row.names(cor.out.rest), file=paste0(sub_dir,"Genes Positively Correlated with Healthy Tumor Border Score outer vs inner_just_genes.txt"), quote=FALSE, sep="\t",row.names = FALSE, col.names = FALSE)

write.table(row.names(data.subset), 
            file=paste0(sub_dir,"Genes all_just_genes.txt"), quote=FALSE, sep="\t",row.names = FALSE, col.names = FALSE)


head(cor.out.outer)


##### create GO dataframe lifting information from the saved GO results in :
###     Shalek Lab/Projects/SpaceCat/Experiments/191021_Tumor_1/Analysis/2020_Revision_Analysis/output/All_Data_no_core_no_low_quality/GO_enrichment_tumor_border.xlsx

df = data.frame(GO.Num = NA, GO.Term = NA, Count = NA, PValue = NA, Up.In.Border = NA)
df[1,] = c('GO:0002376','Inflammatory Response',30,as.numeric(6.70E-11), 1)
df[2,] = c('GO:0005615','Extracellular Space',52,as.numeric(1.13E-06), 1)
df[3,] = c('mmu04620','TLR Signaling Pathway',12,as.numeric(3.58E-05), 1)
df[4,] = c('mmu04060','Cytokine-Cytokine Receptor Interaction',17,as.numeric(3.06E-04), 1)
df[5,] = c('IPR000975','Interleukin-1',4,as.numeric(3.10E-04), 1)
df[6,] = c('GO:0007159','Leukocyte Cell-Cell Adhesion',5,as.numeric(8.10E-04), 1)
df[7,] = c('GO:0030670','Phagocytic Vesicle Membrane',11,as.numeric(4.62E-09), 1)

df[8,] = c('GO:0000166','Nucleotide Binding',66,as.numeric(1.15E-08), -1)
df[9,] = c('GO:0005739','Mitochondrion',58,as.numeric(2.22E-08), -1)
df[10,] = c('GO:0098641','Cadherin Binding Involved in Cell-Cell Adhesion',18,as.numeric(3.80E-06), -1)
df[11,] = c('GO:0006412','Translation',22,as.numeric(1.88E-06), -1)
df[12,] = c('GO:0006457','Protein Folding',11,as.numeric(3.79E-05), -1)
df[13,] = c('GO:0005856','Cytoskeleton',34,as.numeric(2.34E-04), -1)
df[14,] = c('GO:0005753','Mitochondrial Proton-Transporting ATP Synthase Complex',7,as.numeric(1.99E-07), -1)
df[15,] = c('mmu00020','TCA Cycle',6,as.numeric(5.39E-04), -1)

df$PValue = as.numeric(df$PValue)
df$Up.In.Border = as.numeric(df$Up.In.Border)
df$Transformed.PValue = -log10(df$PValue)*df$Up.In.Border
df = df[order(df$Transformed.PValue, decreasing = FALSE),]
df$GO.Term = factor(df$GO.Term, levels = df$GO.Term)
df$Up.In.Border = as.factor(df$Up.In.Border)

g = ggplot(df, aes(x = GO.Term, y = Transformed.PValue, fill = Up.In.Border)) + 
  geom_hline(yintercept = c(-log10(0.05),log10(0.05)), linetype = 'dashed', color = 'gray70') + 
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid = element_blank()) + 
  coord_flip() +
  #ylab(expression(Signed~log["10"](p-value))) + 
  scale_fill_manual(values=c("dodgerblue2", "firebrick2")) + 
  labs(title = 'Gene Ontologies Correlated with Healthy/Tumor Border Score', 
       y = expression(Signed~log["10"](p-value)),
       x = NULL) + 
  scale_x_discrete(position = "top") 
g

#install.packages('tidyverse')
require(tidyverse)
pdf(paste0(sub_dir,"GO_Term_enrichment_tumor_border_vs_rest.pdf"), width=6, height=7, useDingbats=FALSE)
g
dev.off()



table(data.subset@active.ident)
table(data.subset@meta.data$Tumor)


#################################################
########                                #########
########  Lets compare correlation to   #########
########       Monocle Pseudotime       #########
########                                #########
#################################################

list.files(out_dir)
sub_dir = paste0(out_dir,'monocyte_macrophage_subset/')
if(!file.exists(sub_dir)){dir.create(sub_dir)}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

remotes::install_github('satijalab/seurat-wrappers')
library(monocle3)
library(SeuratWrappers)

table(data.subset$Sort)
table(data.subset$Array)

data.subset.cds <- as.cell_data_set(data.subset)
data.subset.cds <- cluster_cells(cds = data.subset.cds, reduction_method = "UMAP")
data.subset.cds <- learn_graph(data.subset.cds, use_partition = TRUE)

# load the pre-selected HSCs
#hsc <- readLines("../vignette_data/hsc_cells.txt") 

data.subset.cds <- order_cells(data.subset.cds, reduction_method = "UMAP")

DimPlot(data.subset)

# plot trajectories colored by pseudotime
plot_cells(
  cds = data.subset.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

data.subset <- AddMetaData(
  object = data.subset,
  metadata = data.subset.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime.Original"
)

data.subset = FindVariableFeatures(data.subset, nfeatures = 2000)
data.subset = ScaleData(object = data.subset, do.scale=TRUE, do.center=TRUE)
range(data.new.all.seurat@assays$RNA@scale.data)

g = VariableFeaturePlot(data.subset, assay = 'RNA', log = 2)
g
pdf(paste0(sub_dir,"MeanVarPlot_nonames.pdf"), useDingbats = FALSE)
g
dev.off()

top <- head(x = VariableFeatures(object = data.subset), 
            n = 40) # 40 appeared to be the right number, but feel free to toggle

# Plot variable features with labels
g = LabelPoints(plot = g, 
                points = top, 
                repel = TRUE)
g

pdf(paste0(sub_dir,"MeanVarPlot_names.pdf"), useDingbats = FALSE)
g
dev.off()

#run pcs over variable genes
#pcs.range = 1:100
data.subset = RunPCA(data.subset, features = NULL, npcs = 100) #run on variable genes, 100 PCs


ElbowPlot(data.subset, ndims=100) + geom_vline(xintercept = 15, col = 'red')
pdf(paste0(sub_dir,"PCjackstraw_elbow.pdf"), useDingbats = FALSE)
ElbowPlot(data.subset, ndims=100) + geom_vline(xintercept = 15, col = 'red')
dev.off()

top = c(1:3)
bottom = c(100)
sets.of = 9

t = list()
t$first = min(top,bottom)
t$last = max(top,bottom)
t$start.of.sets = seq(from = t$first, to = t$last, by = sets.of)
t$start.of.sets

g = lapply(t$start.of.sets , function(x){
  if(x+sets.of>dim(data.subset@reductions$pca)[2]){
    temp = dim(data.subset@reductions$pca)[2]
  } else{
    temp = x + sets.of - 1 # adjustments to final heatmap to plot, should do 3x3 grid
  }
  g = DimHeatmap(data.subset, dims = c(x:temp), cells = 500, balanced = TRUE)
  return(g)
})

print(g[[1]])

pdf(paste0(sub_dir,"PCHeatmap_",x,"-",temp,".pdf"), useDingbats = FALSE)
g
dev.off()




sig.pcs = c(1:9)

#UMAP etc to get new figure 4B
set.seed(123)
data.subset <- FindNeighbors(data.subset, reduction = "pca", dims = sig.pcs, nn.eps = 0.5)
data.subset <- RunUMAP(data.subset, dims = sig.pcs, min.dist = 0.5, n.neighbors = 30, seed.use = 123)

# cols = Clusters.Colors.2,
g = DimPlot(data.subset, reduction = 'umap', pt.size=.6, label=T, group.by="Sort", repel = T) + 
  theme(legend.position = "none", 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

pdf(paste0(sub_dir,"umap_monocyte_macrophage_by_sort.pdf"), useDingbats = FALSE)
g  
dev.off()

g = DimPlot(data.subset, reduction = 'umap', pt.size=.6, label=T, group.by="ident", repel = T) + 
  theme(legend.position = "none", 
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

pdf(paste0(sub_dir,"umap_monocyte_macrophage_by_cluster.pdf"), useDingbats = FALSE)
g  
dev.off()

data.subset.cds.2 <- as.cell_data_set(data.subset)
data.subset.cds.2 <- cluster_cells(cds = data.subset.cds.2, reduction_method = "UMAP")
data.subset.cds.2 <- learn_graph(data.subset.cds.2, use_partition = TRUE)

# load the pre-selected HSCs
#hsc <- readLines("../vignette_data/hsc_cells.txt") 
# this is just a list of the cells to start with - we can either supply nothing or use something liek the core or border?

data.subset.cds.2 <- order_cells(data.subset.cds.2, reduction_method = "UMAP")

g = plot_cells(
  cds = data.subset.cds.2,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

data.subset <- AddMetaData(
  object = data.subset,
  metadata = data.subset.cds.2@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime.Second.New.UMAP"
)


data.subset 
data.subset = ChooseClusterResolutionDownsample(data.subset, n.pcs = sig.pcs, 
                                                sample.name = "LC_Methods_Top", bias="under", 
                                                res.low = .1, res.high = 2,
                                                res.n=20, 
                                                proportion.downsample = 1,
                                                dir = sub_dir)

data.subset <- FindClusters(data.subset, resolution = .117)

table(data.subset@meta.data$RNA_snn_res.0.117)
DimPlot(data.subset,group.by = 'RNA_snn_res.0.117')

data.subset <- FindClusters(data.subset, resolution = .2)
DimPlot(data.subset,group.by = 'RNA_snn_res.0.2')

data.subset <- FindClusters(data.subset, resolution = .5)
DimPlot(data.subset,group.by = 'RNA_snn_res.0.5')

data.subset <- FindClusters(data.subset, resolution = .3)
DimPlot(data.subset,group.by = 'RNA_snn_res.0.3')

table(data.subset@active.ident)


data.subset.cds.3 <- as.cell_data_set(data.subset)
data.subset.cds.3 <- cluster_cells(cds = data.subset.cds.3, reduction_method = "UMAP")
data.subset.cds.3 <- learn_graph(data.subset.cds.3, use_partition = TRUE)

data.subset <- AddMetaData(
  object = data.subset,
  metadata = data.subset.cds.3@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime.Third.New.UMAP.SeuratCluster"
)

which(order(data.subset$Pseudotime.Second.New.UMAP) == order(data.subset$Pseudotime.Third.New.UMAP.SeuratCluster))

which(order(data.subset$Pseudotime.Second.New.UMAP) == order(data.subset$Pseudotime.Third.New.UMAP.SeuratCluster))

# basically the third & second look very similar with some inversions
plot(data.subset$Pseudotime.Second.New.UMAP, data.subset$Pseudotime.Third.New.UMAP.SeuratCluster)
#the first & second look sorta similar and are certainly corrrelated, but do not show the same structure
plot(data.subset$Pseudotime.Second.New.UMAP, data.subset$Pseudotime.Original)
abline(lm(data.subset$Pseudotime.Original~ data.subset$Pseudotime.Second.New.UMAP), col = 'red')


# load the pre-selected HSCs
#hsc <- readLines("../vignette_data/hsc_cells.txt") 
# this is just a list of the cells to start with - we can either supply nothing or use something liek the core or border?

data.subset.cds.3 <- order_cells(data.subset.cds.3, reduction_method = "UMAP")

g = plot_cells(
  cds = data.subset.cds.3,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE) +   
  #geom_point(size = 1, aes(color = pseudotime)) +
  scale_color_gradient(low="dodgerblue4", high ="firebrick2")


pdf(paste0(sub_dir,"umap_monocyte_macrophage_by_pseudotime.pdf"), useDingbats = FALSE)
g  
dev.off()


plot(x = rank(data.subset$Pseudotime.Third.New.UMAP.SeuratCluster), 
     y = rank(data.subset$Outer_Score),
     pch = 20, xlab = 'Monocle Pseudotime Rank' , ylab = 'Healthy/Tumor Border Score Rank')

pdf(paste0(sub_dir,"healthy_tumor_border_rank_vs_pseudotime_rank.pdf"), useDingbats = FALSE)
plot(x = rank(data.subset$Pseudotime.Third.New.UMAP.SeuratCluster), 
     y = rank(data.subset$Outer_Score),
     pch = 20, xlab = 'Monocle Pseudotime Rank' , ylab = 'Healthy/Tumor Border Score Rank')
dev.off()

x = 30
plot(1:x,1:x,pch=1:x)


##############################################
##############################################
######  DE between mono/mac tumor 1 vs mono/mac tumor 2 ############
##############################################
##############################################

mono.mac.tumor.1.2.de = FindMarkers(data.new.all.seurat,
                                    test.use = 'MAST',
                                    ident.1 = WhichCells(data.new.all.seurat,
                                                         ident = 'Monocytes or Macrophages',
                                                         cells = which(data.new.all.seurat$Tumor == 'Tumor 1')),
                                    ident.2 = WhichCells(data.new.all.seurat,
                                                         ident = 'Monocytes or Macrophages',
                                                         cells = which(data.new.all.seurat$Tumor == 'Tumor 2')))

temp = WhichCells(data.new.all.seurat,
                  ident = 'Monocytes or Macrophages',
                  cells = which(data.new.all.seurat$Tumor == 'Tumor 1'))
temp3 = WhichCells(data.new.all.seurat,
                   ident = 'Monocytes or Macrophages',
                   cells = which(data.new.all.seurat$Tumor == 'Tumor 2'))
df.temp = data.frame(cell.id = c(temp,temp3), group.id = factor(c(rep(1,length(temp)),
                                                                  rep(2,length(temp3)))))

df.temp2 = data.new.all.seurat@assays$RNA[row.names(mono.mac.tumor.1.2.de),
                                      df.temp$cell.id]

temp4 = apply(df.temp2,1,function(x){
  return(cohen.d(x ~ df.temp$group.id)$estimate)})

mono.mac.tumor.1.2.de.cohen.d = cbind(mono.mac.tumor.1.2.de,
                                      cohens.d = temp4)



mono.mac.tumor.1.2.de.cohen.d$color = 'black'
mono.mac.tumor.1.2.de.cohen.d$color[intersect(which(mono.mac.tumor.1.2.de.cohen.d$p_val_adj < 0.001),
          which(abs(mono.mac.tumor.1.2.de.cohen.d$cohens.d) > 0.5))] = 'red'

mono.mac.tumor.1.2.de.cohen.d$name = NA

mono.mac.tumor.1.2.de.cohen.d$name[which(mono.mac.tumor.1.2.de.cohen.d$color == 'red')] = row.names(mono.mac.tumor.1.2.de.cohen.d)[which(mono.mac.tumor.1.2.de.cohen.d$color == 'red')]

length(intersect(which(mono.mac.tumor.1.2.de.cohen.d$p_val_adj < 0.001),
                 which(abs(mono.mac.tumor.1.2.de.cohen.d$cohens.d) > 0.5)))

require(ggrepel)
g = ggplot(data = mono.mac.tumor.1.2.de.cohen.d, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
  geom_point() + 
  scale_color_manual(values = c('black', 'red')) + 
  geom_label_repel(aes(label=name),label.padding = 0.1, ) + 
  theme_bw() +
  theme(legend.position = 'none') + 
  labs(x = 'Average Log2 Fold Change', y = '-Log10 (Adjusted P-value)', title = 'Differential Expression between Monocytes or Macrophages from Tumor 1 vs Tumor 2')

pdf(paste0(sub_dir,"monocyte_macrophage_tumor_DE_volcano_plot.pdf"), width=12, height=7, useDingbats=FALSE)
g
dev.off()

pdf(paste0(sub_dir,"Xist_expression_by_tumor.pdf"), width=4, height=4, useDingbats=FALSE)
VlnPlot(data.all.seurat, features = 'Xist', group.by = 'Tumor')
dev.off()

pdf(paste0(sub_dir,"Xist_expression_by_cluster_tumor.pdf"), width=12, height=4, useDingbats=FALSE)
VlnPlot(data.all.seurat, features = 'Xist', group.by = 'Clusters.Names', split.by = 'Tumor')
dev.off()


