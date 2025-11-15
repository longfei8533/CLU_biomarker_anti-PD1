library(data.table)
library(tidyverse)

anno_dat <- fread(file = "Data/Annotation_combine.txt") %>% as.data.frame()
sampleDat <- read.csv("Data/sample_data.csv")
dat <- fread("Data/NAguideR_seqknn.csv") %>% as.data.frame()
genes <- anno_dat$`Gene name`[match(dat[,1],anno_dat$`Protein accession`)]
dat <- dat[,-1]
dat <- log2(dat)
dat <- data.frame(genes = genes,dat)

# ------------------------------------- ssGSEA --------------------------------
library(GSVA)

cell_dat <- readLines("Data/CellReports.txt")

geneSet <- list()
for(i in cell_dat){
  itms <- strsplit(i,"\t") %>% .[[1]]
  geneSet[[itms[1]]] <- itms[3:length(itms)]
  print(itms[1])
}

ssdat <- dat[!duplicated(dat$genes),]
ssmat <- ssdat[,-1] %>% as.matrix()
rownames(ssmat) <- ssdat$genes

ss_res <- gsva(ssmat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=T)


# -------------------------------- plot ----------------------------------------------
# t-test

pvalue_less <- c()
pvalue_greater <- c()

for(i in 1:nrow(ss_res)){
  tRes_less <- t.test(ss_res[i,][sampleDat$Label == 0],
                      ss_res[i,][sampleDat$Label == 1],alternative = "less")
  tRes_greater <- t.test(ss_res[i,][sampleDat$Label == 0],
                 ss_res[i,][sampleDat$Label == 1],alternative = "greater")
  pvalue_less <- c(pvalue_less,tRes_less$p.value)
  pvalue_greater <- c(pvalue_greater,tRes_greater$p.value)
}

label_ <- c()

for(i in 1:length(pvalue_greater)){
  l = ""
  if(pvalue_less[i] < 0.05 | pvalue_greater[i] < 0.05){l = "*"}
  if(pvalue_less[i] < 0.01 | pvalue_greater[i] < 0.01){l = "**"}
  if(pvalue_less[i] < 0.001 | pvalue_greater[i] < 0.001){l = "***"}
  label_ <- c(label_,l)
}

# label_dat

label_dat <- data.frame(x = rownames(ss_res), y = 2.8 , label = label_)

ssplot_dat <- ss_res %>% t() %>% scale() %>% as.data.frame()
ssplot_dat$Group <- c(rep("Non-responders",11),rep("Responders",17))
ssplot_dat <- tidyr::pivot_longer(ssplot_dat, -Group,
                    names_to = "Cell type", values_to = "Score")


sel_label_dat <- label_dat[label_dat$label != "",]
sel_ssplot_dat <- ssplot_dat[ssplot_dat$`Cell type` %in% sel_label_dat$x,]


sel_p_box <- ggplot(sel_ssplot_dat,aes(x = `Cell type`, y = Score))+
  geom_boxplot(aes(fill = Group),position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) +
  scale_fill_manual(values = c("#fded95", "#49c2d9"))+
  coord_cartesian(ylim = c(-2.5,3) )+
  labs(x = "")+
  theme_bw()+
  theme( axis.text.x = element_text(angle = 45, hjust = 1 ),
         legend.position = "top"
         )

sel_p_box <- sel_p_box + 
    geom_text(data = sel_label_dat, aes(x = x, y = y, label = label),
            color = "#B33D27",size = 4,fontface = "bold"
            )
png(file = "Results/tme_boxplot.png",width = 10,height = 10, units = "cm",res = 300)
sel_p_box
dev.off()
