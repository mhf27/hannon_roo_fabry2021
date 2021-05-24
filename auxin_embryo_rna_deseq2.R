### Analysis of eMF73 auxin treated embryos

library("ggplot2", lib.loc="~/Library/R/3.6/library")
library("DESeq2")
library("pheatmap")
library(ggrepel)
library(stringr)
library(svglite)
library("ggpubr", lib.loc="~/Library/R/3.6/library")

directory <- "~/Dropbox (hannonlab)/Sequencing Data/Roo/rnaseq/emf91/individual/count/comb_new/"


sampleFiles <- grep(".txt",list.files(directory),value=TRUE)
sampleCondition <- as.character(replicate(3,c("ctrl", "auxin")))
sampleBatch <- as.character(c("a","a","b","b", "c", "c"))
sampleName <- gsub(".comb.txt*", "", sampleFiles)
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          batch = sampleBatch)



ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ batch + condition)

ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = c("ctrl","auxin"))

dds <- DESeq(ddsHTSeq)
res <- results(dds)

resLFC <- lfcShrink(dds, coef="condition_auxin_vs_ctrl", type="apeglm")
#write.table(resLFC, "~/Dropbox (hannonlab)/Sequencing Data/Roo/rnaseq/emf91/emf91_3_rep_resLFC.txt", quote = F, sep = "\t")
#write.table(res, "~/Dropbox (hannonlab)/Sequencing Data/Roo/rnaseq/emf91/emf91_3_rep_res.txt", quote = F, sep = "\t")
resOrdered <- res[order(res$pvalue),]

plotMA(res, ylim=c(-4,4))
plotMA(resLFC, ylim=c(-2,2))

results <- as.data.frame(resLFC)
results$te <- grepl("FBgn", rownames(results))
results$gene <- str_split_fixed(rownames(results), "_", 2)[,2]
results[,1] <- log10(results[,1])

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

plotPCA(vsd, intgroup=c("condition"))



ggplot(results, aes(x=baseMean, y=log2FoldChange)) +
  geom_point(data=results[which(results[,6]=='FALSE'& results[,5] > 0.05),], aes(x=baseMean, y=log2FoldChange), colour="grey",stat="identity", size=2, shape=16) +
  geom_point(data=results[which(results[,6]=='FALSE'& is.na(results[,5])),], aes(x=baseMean, y=log2FoldChange), colour="grey",stat="identity", size=2, shape=16) +
  geom_point(data=results[which(results[,6]=='FALSE'& results[,5] < 0.05),], aes(x=baseMean, y=log2FoldChange), colour="purple",stat="identity", size=2, shape=16) +
  geom_point(data=results[which(results[,6]=='FALSE'& results[,5] < 0.05 & results[,2] < -1),], aes(x=baseMean, y=log2FoldChange), colour="purple",stat="identity", size=2, shape=16) +
  geom_point(data=results[which(results[,6]=='FALSE'& results[,5] < 0.05 & results[,2] > 1),], aes(x=baseMean, y=log2FoldChange), colour="purple",stat="identity", size=2, shape=16) +
  geom_point(data=results[which(results[,6]=='TRUE'& results[,5] > 0.05),], aes(x=baseMean, y=log2FoldChange), colour="darkgreen",stat="identity", size=2, shape=16) +
  geom_point(data=results[which(results[,6]=='TRUE'& is.na(results[,5])),], aes(x=baseMean, y=log2FoldChange), colour="darkgreen",stat="identity", size=2, shape=16) +
  geom_point(data=results[which(results[,6]=='TRUE'& results[,5] < 0.05 & results[,2] < 1),], aes(x=baseMean, y=log2FoldChange), colour="darkgreen",stat="identity", size=2, shape=16) +
  geom_point(data=results[which(results[,6]=='TRUE'& results[,5] < 0.05 & results[,2] > 1),], aes(x=baseMean, y=log2FoldChange), colour="red",stat="identity", size=2, shape=16) +
  geom_text_repel(data=results[which(results[,6]=='TRUE'& results[,5] < 0.05 & results[,2] < 1),], aes(label=gene,hjust=0, vjust=0)) +
  geom_text_repel(data=results[which(results[,6]=='TRUE'& results[,5] < 0.05 & results[,2] > 1),], aes(label=gene,hjust=0, vjust=0)) +
  geom_abline(intercept = 1, slope = 0, linetype="dotted", colour="red") +
  geom_abline(intercept = 0, slope = 0, linetype="dotted", colour="black") +
  geom_abline(intercept = -1, slope = 0, linetype="dotted", colour="blue") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-1.6,1.6),xlim = c(1,5)) +
  ylab("log2 Fold Change auxin vs ctrl") +
  xlab("log10 mean expression") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"))

ggsave(filename="~/Dropbox (hannonlab)/Fabry_2020/updated_figures_feb21/figure4/update_fig4_c.pdf", plot = last_plot(), width=5, height=5)


ggplot(results, aes(x=baseMean, y=log2FoldChange)) +
  geom_point(data=results[which(results[,6]=='FALSE'& results[,5] > 0.05),], aes(x=baseMean, y=log2FoldChange), colour="grey", alpha=0.2) +
  geom_point(data=results[which(results[,6]=='FALSE'& is.na(results[,5])),], aes(x=baseMean, y=log2FoldChange), colour="grey", alpha=0.2) +
  geom_point(data=results[which(results[,6]=='FALSE'& results[,5] < 0.05),], aes(x=baseMean, y=log2FoldChange), colour="grey", alpha=0.25) +
  geom_point(data=results[which(results[,6]=='FALSE'& results[,5] < 0.05 & results[,2] < -1),], aes(x=baseMean, y=log2FoldChange), colour="grey", alpha=0.25) +
  geom_point(data=results[which(results[,6]=='FALSE'& results[,5] < 0.05 & results[,2] > 1),], aes(x=baseMean, y=log2FoldChange), colour="grey", alpha=0.25) +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-1.6,1.6),xlim = c(1,5)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), line = element_blank())

#ggsave(filename="~/Dropbox (hannonlab)/Sequencing Data/Roo/rnaseq/emf91/emf91_3_rep_dm6_te_deseq_genes_only.png", plot = last_plot(), width=5, height=5)


ggplot(results, aes(x=baseMean, y=log2FoldChange)) +
  geom_point(data=results[which(results[,6]=='TRUE'& results[,5] > 0.05),], aes(x=baseMean, y=log2FoldChange), colour="darkgreen", size=2, alpha=0.5) +
  geom_point(data=results[which(results[,6]=='TRUE'& is.na(results[,5])),], aes(x=baseMean, y=log2FoldChange), colour="darkgreen", size=2, alpha=0.5) +
  geom_point(data=results[which(results[,6]=='TRUE'& results[,5] < 0.05 & results[,2] < 1),], aes(x=baseMean, y=log2FoldChange), colour="darkgreen", size=2, alpha=0.5) +
  geom_point(data=results[which(results[,6]=='TRUE'& results[,5] < 0.05 & results[,2] > 1),], aes(x=baseMean, y=log2FoldChange), colour="red", size=2) +
  geom_text_repel(data=results[which(results[,6]=='TRUE'& results[,5] < 0.05 & results[,2] < 1),], aes(label=gene,hjust=0, vjust=0)) +
  geom_text_repel(data=results[which(results[,6]=='TRUE'& results[,5] < 0.05 & results[,2] > 1),], aes(label=gene,hjust=0, vjust=0)) +
  geom_abline(intercept = 1, slope = 0, linetype="dotted", colour="red") +
  geom_abline(intercept = 0, slope = 0, linetype="dotted", colour="black") +
  geom_abline(intercept = -1, slope = 0, linetype="dotted", colour="blue") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-1.6,1.6),xlim = c(1,5)) +
  ylab("log2 Fold Change auxin vs ctrl") +
  xlab("log10 mean expression") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"))

#ggsave(filename="~/Dropbox (hannonlab)/Sequencing Data/Roo/rnaseq/emf91/emf91_3_rep_dm6_te_deseq_tes_only.svg", plot = last_plot(), width=5, height=5)
