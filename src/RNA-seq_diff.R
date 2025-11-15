library(dplyr)
library(DESeq2)
library(data.table)
library(ggrepel)
library(GEOquery)

# ----------------------------------------- Liu ------------------------------------------

count_dat <- read.table("Data/Liu/rnaseq_rawcounts.txt",row.names = 1,header = T)
gene <- rownames(count_dat)
count_dat <- apply(count_dat,c(1,2),round)

sam_dat <- data.table::fread("Data/Liu/sample_info.txt",header = T,fill = T)
sam <- dplyr::filter(sam_dat,BR %in% c("CR","PR","PD"),
                     `biopsyContext (1=Pre-Ipi; 2=On-Ipi; 3=Pre-PD1; 4=On-PD1)`==3)

sam$BR <- ifelse(sam$BR %in% c("CR","PR"), "Responders", "NonResponders")
sam <- sam[sam$V1 %in% colnames(count_dat),]

count <- count_dat[,match(sam$V1,colnames(count_dat))]

colData= data.frame(group=as.factor(sam$BR),
                    purity = sam$purity
                    )

dds <- DESeqDataSetFromMatrix(countData = count, 
                              colData = colData, design = ~ purity + group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", c("NonResponders","Responders")))
liu_res_dat <- as.data.frame(res)

plot_dat <- res[complete.cases(res),]
plot_dat$fill <- "Not-siginaficant"
plot_dat$fill[plot_dat$padj < 0.2 & plot_dat$log2FoldChange < 0 ] <- "Down-regulated" 
plot_dat$fill[plot_dat$padj < 0.2 & plot_dat$log2FoldChange > 0 ] <- "Up-regulated"


plot_txt <- plot_dat[rownames(plot_dat) == "CLU",]
plot_txt$gene <- "CLU"

res_plot <- ggplot(data=plot_dat,aes(x=log2FoldChange,y=-log10(padj),fill = fill))+
  geom_point(alpha=1, size=2,shape = 21,show.legend = T) +
  scale_fill_manual(values = c(`Up-regulated` = c("#FF4040"), `Not-siginaficant` = c("#ABABAB") , `Down-regulated` = c("#1874CD")))+
  geom_hline(yintercept = -log10(0.2),lty=4,color="#DC0000FF",linewidth=0.8) +
  #scale_x_continuous(breaks = c(-5:5))+
  coord_cartesian(xlim = c(-10, 10),ylim = c(0,10))+
  labs(x=expression(log[2]*(Fold~change)), y=expression(log[10]*(P.adjust)))+
  theme_bw() +
  theme(axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=12),
        panel.background = element_blank(),
        legend.title=element_blank()
        #panel.grid =element_line(colour = "gray"),
        #legend.position = "none"
        )
res_plot <- res_plot + 
  geom_text_repel(data = plot_txt,
                  aes(x=log2FoldChange,y=-log10(padj),label = gene),
                  color="#2CA02CFF",
                  size = 3,
                  point.padding = unit(0.8, "lines"),
                  xlim = c(-2,1),ylim = c(4,10),
                  segment.color = "#2CA02CFF",
                  show.legend = FALSE
                  )


png(filename = "Results/Liu_dotplot.png",width = 12,height = 8,units = "cm",res = 300)
res_plot
dev.off()


# ----------------------------------------- Riaz ------------------------------------------
rm(list = ls())

mat <- read.delim("Data/Riaz/CountData.BMS038.txt")
gene <- mat$HUGO

sam<- read.csv("Data/Riaz/SampleTableCorrected.9.19.16.csv",                         row.names=1)
sam <- dplyr::filter(sam,Response %in% c("PRCR","PD"),
                     PreOn == "Pre"
                     )
sam$Response <- ifelse(sam$Response %in% c("PRCR"), "NonResponders","Responders" )

sam <- sam[sam$Sample %in% names(mat),]
mat <- mat[,match(sam$Sample,names(mat))]



colData= data.frame(case=as.factor(sam$PatientID),
                    group = as.factor(sam$Response))

dds <- DESeqDataSetFromMatrix(countData = mat, 
                              colData = colData, design = ~ group)

dds <- DESeq(dds)
res <- results(dds, contrast = c("group", c("NonResponders","Responders")))
res$gene <- gene
Riaz_res_dat <- as.data.frame(res)

plot_dat <- res[complete.cases(res),]
plot_dat$fill <- "Not-siginaficant"
plot_dat$fill[plot_dat$padj < 0.2 & plot_dat$log2FoldChange < 0 ] <- "Down-regulated" 
plot_dat$fill[plot_dat$padj < 0.2 & plot_dat$log2FoldChange > 0 ] <- "Up-regulated"

plot_txt <- plot_dat[plot_dat$gene == "CLU",]

res_plot <- ggplot(data=plot_dat,aes(x=log2FoldChange,y=-log10(padj),fill = fill))+
  geom_point(alpha=1, size=2,shape = 21,show.legend = T) +
  scale_fill_manual(values = c(`Up-regulated` = c("#FF4040"), `Not-siginaficant` = c("#ABABAB") , `Down-regulated` = c("#1874CD")))+
  geom_hline(yintercept = -log10(0.2),lty=4,color="#DC0000FF",linewidth=0.8) +
  #scale_x_continuous(breaks = c(-5:5))+
  coord_cartesian(xlim = c(-10, 10),ylim = c(0,5))+
  labs(x=expression(log[2]*(Fold~change)), y=expression(log[10]*(P.adjust)))+
  theme_bw() +
  theme(axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=12),
        panel.background = element_blank(),
        legend.title=element_blank()
        #panel.grid =element_line(colour = "gray"),
        #legend.position = "none"
        )
res_plot <- res_plot + 
  geom_text_repel(data = plot_txt,
                  aes(x=log2FoldChange,y=-log10(padj),label = gene),
                  color="#2CA02CFF",
                  size = 3,
                  point.padding = unit(0.8, "lines"),
                  xlim = c(5,10),ylim = c(2,3),
                  segment.color = "#2CA02CFF",
                  show.legend = FALSE
                  )


png(filename = "Results/Riaz_dotplot.png",width = 12,height = 8,units = "cm",res = 300)
res_plot
dev.off()


