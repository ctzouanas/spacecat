#####
#Create a Joint eSet of PMA/CEFT/TT/DMSO/etc + HIV Data so we can do all analysis in cleaner  way
####
collect = NULL
out_dir = "/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab/Projects/SpaceCat/Experiments/180731_DMZ/SS2/DMZ_SS2_Analysis/"
lib_dir = "/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab/Collaborations/Nir and Michael/Pipeline 2.0/RCODE"
housekeeping_list = "/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab/Collaborations/Nir and Michael/house_keeping_human_names.txt"

#### source necessary functions
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/Interacting Cells/T_MDDC_RNAseq/HIV Analysis/LoveHIVfunctions2.R')
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/Interacting Cells/T_MDDC_RNAseq/PMA Study/Pairs_PMA_Analysis/LovePMApipeFunctions.R')
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/Interacting Cells/T_MDDC_RNAseq/All_Pairs_Data/All_Pairs_Reannotate_Data.R')
source('/Users/Alex/Dropbox (MIT)/Shalek Lab Top of Tree/Shalek Lab Genshaft Files/Interacting Cells/T_MDDC_RNAseq/HIV Analysis/MAST_wrapper_function.R')
# load all required packages
loadLibraries(lib_dir)

require(openxlsx)
meta.data.list = list()
meta.data.list[['HEK_PA']] = read.xlsx(file = '../RawData/DMZ_Master_Data.xlsx',sheetName = 'Plate1_HEK')
meta.data.list[['3T3_PA']] = read.xlsx(file = '../RawData/DMZ_Master_Data.xlsx',sheetName = 'Plate3_3T3')
meta.data.list[['Neg_PA']] = read.xlsx(file = '../RawData/DMZ_Master_Data.xlsx',sheetName = 'Plate5_unsorted')
head(meta.data.list[['HEK_PA']] )
head(meta.data.list[['3T3_PA']] )
head(meta.data.list[['Neg_PA']] )

meta.data.list.classified = lapply(meta.data.list, function(x){
  temp = rep(NA,dim(x)[1])
  temp2 = apply(x[,2:3],1,max)
  temp[which(temp2 > 100,000)] = 0
  temp[which(temp == 0)] = x[which(temp==0),2] - x[which(temp==0),3]
  temp[which(temp<0)] = 'HUMAN'
  temp[which(temp>0)] = 'MOUSE'
  return(cbind(x,temp))
})

head(meta.data.list.classified[[2]])

### source some other functions
# subsetting eSet
subset_eSet_list = function(eSet_list,
                            cell_subset = colnames(exprs(eSet_list[[1]])),
                            gene_subset = row.names(exprs(eSet_list[[1]]))){
  return(lapply(eSet_list,function(eSet){
    return(eSet[gene_subset,cell_subset])
  }))
}

# QC stats by some basic groupings
box_plot_qc_multiplot = function(eSet = NULL, 
                                 qc.variable.vector = NULL, 
                                 condition.variable = NULL,
                                 save.file = NULL){
  g = list()
  g[[1]]  = box_plot_by_condition(variable.vec = pData(eSet)[,qc.variable.vector[1]],
                                  condition.vec = pData(eSet)[,condition.variable],
                                  ylab.custom = 'Percent Transcrpitome',
                                  title.custom = '% Transcriptome by Library Type',
                                  xlab.custom = FALSE,
                                  points = FALSE)
  g[[2]]  = box_plot_by_condition(variable.vec = pData(eSet)[,qc.variable.vector[2]],
                                  condition.vec = pData(eSet)[,condition.variable],
                                  ylab.custom = 'Reads',
                                  title.custom = 'Reads by Library Type',
                                  xlab.custom = FALSE,
                                  points = FALSE)
  g[[3]]  = box_plot_by_condition(variable.vec = pData(eSet)[,qc.variable.vector[3]],
                                  condition.vec = pData(eSet)[,condition.variable],
                                  ylab.custom = 'Gene Count',
                                  title.custom = 'Genes by Library Type',
                                  xlab.custom = FALSE,
                                  points = FALSE)
  g[[4]]  = box_plot_by_condition(variable.vec = pData(eSet)[,qc.variable.vector[4]],
                                  condition.vec = pData(eSet)[,condition.variable],
                                  ylab.custom = 'Percent RNA',
                                  title.custom = '% RNA by Library Type',
                                  xlab.custom = FALSE,
                                  points = FALSE)
  g[[5]]  = box_plot_by_condition(variable.vec = pData(eSet)[,qc.variable.vector[5]],
                                  condition.vec = pData(eSet)[,condition.variable],
                                  ylab.custom = 'Percent Genome',
                                  title.custom = '% Genome by Library Type',
                                  xlab.custom = FALSE,
                                  points = FALSE)
  #### make this fraction of cells that make it through / have fastqs rather than dont....
  #g[[6]]
  
  figure <- ggarrange(plotlist = g,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol = 3, nrow = 2, 
                      common.legend = TRUE, legend = 'bottom')
  
  if(!is.null(save.file)){
    pdf(file = save.file, useDingbats = FALSE,
        width = 12, height = 7)
    print(figure)
    dev.off()
  }
  
  return(figure)
}

###### lets look at gene characteristics (mean / var)
mean_var_multiplot = function(eSet_list = NULL, 
                              save.file = NULL){
  
  #if(length(eSet_list) > 2){
  #  stop('not expecting eSet list this long, modify code to make it work.')
  #}
  
  g = list()
  g = lapply(eSet_list, function(x){
    df = data.frame(mean = rowMeans(log2(1+exprs(x))), 
                    var = apply(log2(1+exprs(x)),1,var),
                    genes = get_simple_genes(row.names(exprs(x))))
    panels = list()
    panels[[1]] = ggplot(df, aes(x = mean)) + 
      geom_histogram(bins = floor(sqrt(length(df$mean)))) +
      labs(title = 'Histogram of Mean Log2(1+Expression)')
    panels[[2]] = ggplot(df, aes(x = var)) + 
      geom_histogram(bins = floor(sqrt(length(df$var)))) +
      labs(title = 'Histogram of Var Log2(1+Expression)')
    panels[[3]] =  ggplot(df,aes(y=var,x=mean,label=as.character(genes))) +
      geom_point() + 
      theme(axis.ticks = element_line(colour = 'black', size = 0.2),
            axis.text = element_text(colour = 'black', size = 8)) +
      labs(title = 'Variance vs Mean of Log2(1+Expression)')
    panels[[4]] =  ggplot(df,aes(y=var,x=mean,label=as.character(genes))) +
      geom_text() + 
      theme(axis.ticks = element_line(colour = 'black', size = 0.2),
            axis.text = element_text(colour = 'black', size = 8)) +
      labs(title = 'Variance vs Mean of Log2(1+Expression) [names]')
    
    panels = lapply(panels,make_ggplot_pretty)
    
    figure = ggarrange(ggarrange(plotlist = panels[1:2],
                                 labels = c("A", "B"),
                                 nrow = 2),
                       panels[[3]],
                       panels[[4]],
                       ncol = 3,
                       labels = c('','C','D'))
    
    #make them ggplots pretty & get on with it
    return(figure)})
  
  figure <- ggarrange(plotlist = g, 
                      nrow = 2)
  
  if(!is.null(save.file)){
    pdf(file = save.file, useDingbats = FALSE,
        width = 12, height = 7)
    print(figure)
    dev.off()
  }
  
  return(figure)
}

