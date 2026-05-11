setwd("C:/Users/mcooke5/Desktop/JDW_NGS/__SCRIPTS/CRISPRI_SEQ/")
name = 'OEpykF_MinGlu'
counts_initial = read.table('1.0_count_table_OEpykF_MinGlu.txt')

#DESEQ2 1.36.1 ; ggplot2 3.4.2 ; stringr 1.5.0
required_packages <- c("DESeq2", "ggplot2", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)}
  library(pkg, character.only = TRUE)}

##Block 1 - Load Counts and Filter Low Coverage
counts=counts_initial[rowSums(counts_initial)>20,]
samples = data.frame(label=colnames(counts))
samples$replicate <- as.factor(sub(".*_", "", samples$label))
samples$generation <- factor(ifelse(grepl("_T0_", samples$label), 0, ifelse(grepl("_TF_", samples$label), 1, NA)))
samples <- samples[order(samples$replicate, samples$generation), ]
rownames(samples) = samples$label
counts=counts[,samples$label]
remove(counts_initial,pkg,required_packages)
#------------------------------------------------------
##Block 2 - Deseq2 Computations
dds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design = ~ generation + replicate)
controlguides <- readLines("./_controlguides.txt")
targetControl <- as.numeric(colMedians(counts(dds)[controlguides,1:ncol(counts(dds))]))
sizeFactors(dds) <- targetControl/mean(targetControl)
#dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = 'local')
dds <- nbinomWaldTest(dds, betaPrior=FALSE)
res <- results(dds,contrast=c('generation',1,0))
res$target=rownames(res)
writeLines(capture.output(summary(res)), paste0("2.0_",name,"_DESEQstats.txt"))

##Block 2 Metrics
pdf(paste0("2.1_",name,"_DispEst.pdf"),width=6,height=4,colormodel='rgb',paper = 'A4')
plotDispEsts(dds,main = "Dispersion Estimates")
dev.off()

pdf(paste0("2.2_",name,"_MA.pdf"),width=8,height=6,colormodel='rgb',paper = 'A4')
plotMA(dds,alpha=0.05,main="LFC",colNonSig="gray60",colLine="gray40",ylim=c(-4,4))
dev.off()

forpca <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")
pdf(paste0("2.3_",name,"_PCA.pdf"),width=6,height=6,colormodel='rgb',paper = 'A4')
plotPCA(forpca, intgroup = "label",ntop = 500)
dev.off()

remove(counts,dds,forpca,samples,controlguides,targetControl)
#--------------------------------------------------------------------
##Block 3 - Gene Annotation of Guides
guide=read.table('./_EcoWG1_library.txt',header=T)
data=merge(guide[,c("target","pos","gene","ori","coding","gene_right","gene_left"),],
           as.data.frame(res[,c('target','baseMean','log2FoldChange','padj')]),by ='target')
data <- data[order(data$log2FoldChange) & data$baseMean > 10, ]
operons=read.table('./_OPERONS_regulonDB.txt')

score_gene<-function(data){
  temp=data[!is.na(data$gene_left),]
  medians=unique(temp[!is.na(temp[temp$coding!='promoter','gene']),c('gene','gene_left','gene_right','ori','coding')])
  medians[c("median", "mad", "nb", "pos")] <- c(1, 1, 0, 0)
  for (i in 1:dim(medians)[1]){
    if (length(operons[grepl(as.character(medians$gene[i]),operons$genes),]$name)==1){
      medians$operon[i]=as.character(operons[grepl(medians$gene[i],operons$genes),]$name)
      medians$nb[i]=operons[grepl(medians$gene[i],operons$genes),]$number
      g=operons[grepl(medians$gene[i],operons$genes),]$genes
      medians$pos[i]=str_count(substr(g,start = 1,stop=grepRaw(medians$gene[i],g)),pattern = ',')+1
    }else{
      medians$operon[i]=as.character(medians$gene[i])
      medians$nb[i]=1}}
  for (i in 1:dim(medians)[1]){
    medians$median[i]=median(data[!is.na(data$log2FoldChange)& data$gene==medians$gene[i],]$log2FoldChange)
    medians$mad[i]=mad(data[!is.na(data$log2FoldChange)& data$gene==medians$gene[i],]$log2FoldChange,constant=1)}
  medians=medians[order(medians$median),]
  return(medians)}
scores=score_gene(data)
scores <- transform(scores, rown = rownames(scores))
data$rown=rownames(data)
scores <- cbind(scores, data[match(scores$rown, data$rown), c("target", "padj")])

#Final Output
guides=subset(data,select = c('gene','log2FoldChange','padj','target','pos','gene_left','gene_right','coding'))
genes=subset(scores, select = c('gene','median','mad','padj','target','operon','pos','ori','gene_left','gene_right','coding'))
genes_signif=subset(scores, padj<=0.1 & median<=-2 | padj<=0.1 & median>=2, select = c('gene','median','mad','padj','target','operon','pos','ori','gene_left','gene_right','coding'))
write.csv(guides,paste0('2.4_',name,'_Guide.csv'))
write.csv(genes,paste0('2.5_',name,'_Gene.csv'))
write.csv(genes_signif,paste0('2.6_',name,'_Gene_Signif.csv'))

remove(data,genes,genes_signif,guide,guides,operons,res,score_gene)
#-----------------------------------------------------------------------------------

scores=subset(scores, mad!=0)
scores$rank=1:dim(scores)[1]
p=ggplot(scores,aes(rank,median,colour=padj<=0.1 & median<=-2 | padj<=0.1 & median>=2))+
  geom_errorbar(data=scores,aes(ymin=median-mad,ymax=median+mad),color='#999999',width=0.001,linewidth=0.1)+
  geom_point(alpha=0.6,size=2,stroke=0)+
  geom_hline(yintercept = 0)+
  scale_color_manual(values=c('#666666','#CC0000'))+
  geom_hline(yintercept = -2,linetype='dotted')+
  annotate("text", x = c(3900,3900), y = c(-10,-9.2), label = c("Insignificant Change", "Signficant Change"),size=3,hjust = 1)+
  annotate("point", x = c(4000,4000), y = c(-10,-9.2),colour=c('#666666','#CC0000'),alpha=0.6)+
  ylab('Median log2FC')+
  xlab('Rank-ordered genes')+
  scale_y_continuous(breaks=seq(-10, 6, 2))+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),panel.background = element_blank())
pdf(paste0("2.7_",name,"_Rankorder.pdf"),width=8,height=6,colormodel='rgb',paper = 'A4')
print(p)
dev.off()
remove(p,scores,name)
