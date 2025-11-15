library(tidyverse)
library(readr)
library(ggrepel)
library(ggstatsplot)
library(ggpubr)
library(pROC)
library(GSVA)

dat <- read_tsv(file = "Data/Annotation_combine.txt")
inten_dat <- dplyr::select(dat,starts_with("LFQ intensity"))
RP_dat <- dat[,149:150]
info_dat <- read.csv(file = "Data/clinical_information.csv")
info_dat$ID <- paste0("LFQ intensity GC", info_dat$ID)
info_dat <- info_dat[match(colnames(inten_dat),info_dat$ID),]
info_dat$responder <- ifelse(info_dat$`Clinical.outcomes`=="Progressors","Non-responders","Responders")

invalid_dat <- inten_dat[,info_dat$`Clinical.outcomes` == "Progressors"]
valid_dat <- inten_dat[,info_dat$`Clinical.outcomes` == "Responders"]

diff_fun <- function(d1,d2){
  ratio <- c()
  pvalue <- c()
  for(i in 1:nrow(d1)){
    x1 <- as.numeric(d1[i,])
    x2 <- as.numeric(d2[i,])
    x1 <- x1[x1!=0]
    x2 <- x2[x2!=0]
    if(length(x1)<2 | length(x2)<2){
      r <- NA
      p <- NA
    }else{
      r <- mean(x1)/mean(x2)
      p <- t.test(log2(x1),log2(x2),var.equal = T)$p.value
    }
    ratio <- c(ratio,r)
    pvalue <- c(pvalue,p)
  }
  return(data.frame(ratio,pvalue))
}
diff_res <- diff_fun(invalid_dat,valid_dat)
res <- cbind(dat[,1:3],diff_res)


## volcano plot

plot_dat <- res[complete.cases(res),]
plot_dat$full <- "NS"
plot_dat$full[plot_dat$pvalue < 0.05 & plot_dat$ratio < (1/1.5)] <- "Down"
plot_dat$full[plot_dat$pvalue < 0.05 & plot_dat$ratio > 1.5] <- "Up"

plot_txt <- plot_dat[plot_dat$`Gene name` == "CLU",]

res_plot <- ggplot(data=plot_dat,aes(x=log2(ratio),y=-log10(pvalue),fill = full))+
  geom_point(alpha=1, size=3,shape = 21) +
  scale_fill_manual(values = c(Up = "#C21520", NS = "#80796BFF" , Down = "#154B8C"))+
  #geom_hline(yintercept = -log10(0.05),lty=4,color="#DC0000FF",size=0.8) +
  #scale_x_continuous(breaks = c(-5:5))+
  xlim(c(-5,5))+
  labs(x="log2(Fold change)",
       y="-log10 (P value)")+
  theme_bw()+
  theme(axis.text = element_text(size=12,color = "black"),
        axis.title = element_text(size=12),
        panel.background = element_blank(),
        legend.title=element_blank()
        #panel.grid =element_line(colour = "gray"),
        #legend.position = "none"
        )
res_plot <- res_plot + geom_text_repel(data = plot_txt,
                  aes(x=log2(ratio),y=-log10(pvalue),label = `Gene name`),
                  colour="red",
                  point.padding = unit(0.8, "lines"),
                  #segment.color = "black",
                  show.legend = F)
png(filename = "../result/dotplot.png",width = 15,height = 10,units = "cm",res = 300)
res_plot
dev.off()

## KEGG

library(clusterProfiler)
library(org.Hs.eg.db)
sig_diff <- dplyr::filter(res,ratio< (1/1.5) | ratio > 1.5,pvalue<0.05)

ekegg=enrichKEGG(gene=sig_diff$`Protein accession`,
                 organism="hsa",
                 keyType="uniprot",
                 pAdjustMethod = "BH")
                 #pvalueCutoff = 0.05,
                # qvalueCutoff = 0.1)

kegg_dat <- ekegg@result %>% arrange(pvalue) %>% .[1:10,]
kegg_dat$Description <- factor(kegg_dat$Description,levels = rev(kegg_dat$Description))

kegg_plot <- ggplot(data=kegg_dat,aes(y = Description, x = -log10(pvalue))) +
  geom_col(aes(fill=Count))+
  scale_fill_gradient(low = "#154B8C",high = "#C51321")+
  labs(x = "-log10(P value)",y = "") +
  theme(
    axis.title = element_text(size=12),
    axis.text.y = element_text(size=12,color = "black"),
    axis.text.x = element_text(size=12,color = "black"),
    axis.line = element_line(size=0.8,color="black"),
    axis.ticks= element_line(size=0.8,colour = "black"),
    panel.grid =element_blank(),
    panel.background = element_blank()
    #plot.margin = unit(c(0,0,0,7),units = "cm")
    )

png(filename = "Results/kegg.png",width = 20,height = 10,units = "cm",res = 300)
kegg_plot
dev.off()


# ------------------------ ROC AUC ---------------------------------

roc_dat <- data.frame(ex=t(inten_dat[dat$`Gene name`=="CLU",]),
                      group = info_dat$`Clinical.outcomes` == "Progressors")
roc_dat <- roc_dat[roc_dat$ex > 0,]
roc_res <- roc(group~ex,roc_dat,plot=T)


png(filename = "Results/CLU_roc.png",
    width =12 ,height =12 ,units = "cm",res = 300)
plot(roc_res, main = "",
     col = "#9A32CD",#"#3CB371", 
     lwd = 3, 
     print.auc=TRUE,
     print.auc.y=0.4
     )
dev.off()

