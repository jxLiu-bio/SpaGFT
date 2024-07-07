# preload function 
jaccard_similarity <- function(A, B) {
  intersection = length(intersect(A, B))
  union = length(A) + length(B) - intersection
  return (intersection/union)
}
library(ggplot2)
#
setwd("e:/spaGFT/")
# set up markers
my.marker <- openxlsx::read.xlsx("marker_cortical layer marker.xlsx")
table(my.marker$species)
human_marker <- unique(my.marker$`gene_symbol(H)`)
mouse_marker <- unique(my.marker$`gene_symbol(M)`)
# set up sample names
#  check data
my.sample.list <- c("MOB-Rep-10_ST_MouseBrain_expression.h5seurat",
                    "GSE152506-N02-D1_ST_MouseBrain_expression.h5seurat",
                    "HE_Visum_MouseBrain_expression.h5seurat",
                    "Puck-200115-08_Slideseq-v2_MouseBrain_expression.h5seurat",
                    "2-5_Visum_HumanBrain_expression.h5seurat",
                    "MOB-Rep-1_ST_MouseBrain_expression.h5seurat",
                    "MOB-Rep-11_ST_MouseBrain_expression.h5seurat")

my.sample.name <- gsub("_expression.h5seurat","",my.sample.list)
all.sample.list <- list.files(path = "e:/spaGFT/combination/",pattern = ".csv")
# setup software
my.software <- c("MERINGUE","SpaGCN","SPARK","SPARK-X","SpatialDE","SpaGFT")
# Jaccard, Szymkiewicz–Simpson coefficient, Tversky index, Sørensen coefficient, F1

my.final.matrix <- list()

for (i in 1:length(my.software)){
  tmp.software <- my.software[i]
  tmp.res.list <- grep(tmp.software,all.sample.list,value = T)
  print(i)
  my.score.matrix <- c()
  for (j in 1:length(tmp.res.list)){
    tmp.sample.name <- tmp.res.list[j]
    file.name <- file.path("e:/spaGFT/combination/",tmp.res.list[j])
    # file.existence <- file.exists(file.name)
    # if(!file.existence){
    #   my.score_Jaccard[tmp.sample.name,tmp.software] <- "NA"
    #   my.score_OverlapCoef[tmp.sample.name,tmp.software]<- "NA"
    #   my.score_TverskyIndex[tmp.sample.name,tmp.software] <- "NA"
    #   next()
    # }
    # read results
    tmp_results <- read.csv(file.name,row.names = 1)
    # filter 
    if (tmp.software == "MERINGUE"){
      out_results <- tmp_results[tmp_results$p.adj < 0.05,]
    } 
    if (tmp.software == "SpaGCN"){
      out_results <- tmp_results[tmp_results$pvals_adj <0.05,]
      out_results <- out_results[!duplicated(out_results$genes),]
      rownames(out_results) <- out_results$genes
    } 
    if (tmp.software == "SPARK"){
      out_results <- tmp_results[tmp_results$adjusted_pvalue <0.05,]
    } 
    if (tmp.software == "SPARK-X"){
      out_results <- tmp_results[tmp_results$adjustedPval <0.05,]
    } 
    if (tmp.software == "SpatialDE"){
      out_results <- tmp_results[tmp_results$BIC > 0 & tmp_results$qval <0.05,]
      rownames(out_results) <- out_results$g
    }
    if (tmp.software == "SpaGFT"){
      tmp_results$qvalue <- p.adjust(tmp_results$pvalue)
      out_results <- tmp_results[tmp_results$cutoff_gft_score == "True" & tmp_results$qvalue < 0.05,]
      rownames(out_results) <- out_results$features
    }
    predicted.gene <- rownames(out_results)
    # select human and mouse marker
    if(grepl("human",file.name,ignore.case = T)){gt_marker = human_marker}
    if(grepl("mouse",file.name,ignore.case = T)){gt_marker = mouse_marker}
    score_Jaccard <- jaccard_similarity(predicted.gene,gt_marker)
    score_Overlap_coef <- ribiosUtils::overlapCoefficient(predicted.gene,gt_marker)
    score_TverskyIndex <- tcR::tversky.index(predicted.gene,gt_marker)
    score_FisherStatistic <- GeneOverlap::newGeneOverlap(predicted.gene,gt_marker,genome.size = 30000)
    score_FisherStatistic <- GeneOverlap::testGeneOverlap(score_FisherStatistic)
    score_FisherStatistic <- score_FisherStatistic@odds.ratio
    tmp.sampleID <- sapply(strsplit(tmp.sample.name,"_"), "[",1)
    tmp.tech <-  sapply(strsplit(tmp.sample.name,"_"), "[",2)
    tmp.parameter <-  sapply(strsplit(tmp.sample.name,"_"), "[",5)
    
    tmp.score.matrix <- c(tmp.sampleID ,
                          tmp.software,
                          tmp.tech ,
                          tmp.parameter,
                          score_Jaccard,score_Overlap_coef,score_TverskyIndex,score_FisherStatistic)
    my.score.matrix <- rbind(my.score.matrix,tmp.score.matrix)
  }
  colnames(my.score.matrix) <- c("sampleID",	
                                 "software",
                                 "tech",
                                 "parameter",
                                 "score_Jaccard",
                                 "score_Overlap_coef",
                                 "score_TverskyIndex",
                                 "score_FisherStatistic")
  my.final.matrix[[i]] <- my.score.matrix
  names(my.final.matrix)[i] <- tmp.software
}
my.software

xlsx::write.xlsx(my.final.matrix[[1]],file = "e:/spaGFT/Benchmark_combinationres.xlsx",sheetName = "MERINGUE")
xlsx::write.xlsx(my.final.matrix[[2]],file = "e:/spaGFT/Benchmark_combinationres.xlsx",sheetName = "SpaGCN" ,append = T)
xlsx::write.xlsx(my.final.matrix[[3]],file = "e:/spaGFT/Benchmark_combinationres.xlsx",sheetName = "SPARK",append = T)
xlsx::write.xlsx(my.final.matrix[[4]],file = "e:/spaGFT/Benchmark_combinationres.xlsx",sheetName = "SPARK-X",append = T)
xlsx::write.xlsx(my.final.matrix[[5]],file = "e:/spaGFT/Benchmark_combinationres.xlsx",sheetName = "SpatialDE",append = T)
xlsx::write.xlsx(my.final.matrix[[6]],file = "e:/spaGFT/Benchmark_combinationres.xlsx",sheetName = "SpaGFT",append = T)

all.ressult <- rbind(my.final.matrix[[1]],
                     my.final.matrix[[2]],
                     my.final.matrix[[3]],
                     my.final.matrix[[4]],
                     my.final.matrix[[5]],
                     my.final.matrix[[6]])
all.ressult <- as.data.frame(all.ressult)
colnames(all.ressult) <- c("sampleID",	
                           "software",
                           "tech",
                           "parameter",
                           "score_Jaccard",
                           "score_Overlap_coef",
                           "score_TverskyIndex",
                           "score_FisherStatistic")
my.sample.list
HE_Visium.df <- all.ressult[all.ressult$sampleID == "HE",]
GSE152506_ST.df <- all.ressult[all.ressult$sampleID == "GSE152506-N02-D1",]
MOB_ST.df <- all.ressult[all.ressult$sampleID == "MOB-Rep-10",]
Puck_slideseq2.df <- all.ressult[all.ressult$sampleID == "Puck-200115-08",]
S2_5_Visium.df <- all.ressult[all.ressult$sampleID == "2-5",]
MOB_ST_1.df <- all.ressult[all.ressult$sampleID == "MOB-Rep-1",]
MOB_ST_11.df <- all.ressult[all.ressult$sampleID == "MOB-Rep-11",]

my.df <- HE_Visium.df
plt.ID <- paste0(unique(my.df$sampleID),"_",unique(my.df$tech))
p_HE <- ggplot(my.df,aes(x = software,y =as.numeric(score_Jaccard),fill = software ))+
  geom_boxplot()+ggtitle(plt.ID)+theme_classic()+ylim(0,0.18) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = as.character(jcolors::jcolors(palette = "pal8"))[1:6],breaks = c("MERINGUE","SpaGCN","SPARK","SPARK-X","SpatialDE","SpaGFT"))


my.df <- GSE152506_ST.df
plt.ID <- paste0(unique(my.df$sampleID),"_",unique(my.df$tech))
p_GSE152506_ST <- ggplot(my.df,aes(x = software,y =as.numeric(score_Jaccard),fill = software ))+
  geom_boxplot()+ggtitle(plt.ID)+theme_classic()+ylim(0,0.18)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = as.character(jcolors::jcolors(palette = "pal8"))[1:6],breaks = c("MERINGUE","SpaGCN","SPARK","SPARK-X","SpatialDE","SpaGFT"))

my.df <- MOB_ST.df
plt.ID <- paste0(unique(my.df$sampleID),"_",unique(my.df$tech))
p.MOB_ST<- ggplot(my.df,aes(x = software,y =as.numeric(score_Jaccard),fill = software ))+
  geom_boxplot()+ggtitle(plt.ID)+theme_classic()+ylim(0,0.18)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = as.character(jcolors::jcolors(palette = "pal8"))[1:6],breaks = c("MERINGUE","SpaGCN","SPARK","SPARK-X","SpatialDE","SpaGFT"))

my.df <- MOB_ST_1.df
plt.ID <- paste0(unique(my.df$sampleID),"_",unique(my.df$tech))
p.MOB_ST.1<- ggplot(my.df,aes(x = software,y =as.numeric(score_Jaccard),fill = software ))+
  geom_boxplot()+ggtitle(plt.ID)+theme_classic()+ylim(0,0.18)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = as.character(jcolors::jcolors(palette = "pal8"))[1:6],breaks = c("MERINGUE","SpaGCN","SPARK","SPARK-X","SpatialDE","SpaGFT"))

my.df <- MOB_ST_11.df
plt.ID <- paste0(unique(my.df$sampleID),"_",unique(my.df$tech))
p.MOB_ST.11<- ggplot(my.df,aes(x = software,y =as.numeric(score_Jaccard),fill = software ))+
  geom_boxplot()+ggtitle(plt.ID)+theme_classic()+ylim(0,0.18)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = as.character(jcolors::jcolors(palette = "pal8"))[1:6],breaks = c("MERINGUE","SpaGCN","SPARK","SPARK-X","SpatialDE","SpaGFT"))


my.df <- Puck_slideseq2.df
plt.ID <- paste0(unique(my.df$sampleID),"_",unique(my.df$tech))
p.puck <- ggplot(my.df,aes(x = software,y =as.numeric(score_Jaccard),fill = software ))+
  geom_boxplot()+ggtitle(plt.ID)+theme_classic()+ylim(0,0.18)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = as.character(jcolors::jcolors(palette = "pal8"))[1:6],breaks = c("MERINGUE","SpaGCN","SPARK","SPARK-X","SpatialDE","SpaGFT"))


my.df <- S2_5_Visium.df
plt.ID <- paste0(unique(my.df$sampleID),"_",unique(my.df$tech))
p.S2_5_Visium <- ggplot(my.df,aes(x = software,y =as.numeric(score_Jaccard),fill = software ))+
  geom_boxplot()+ggtitle(plt.ID)+theme_classic() +ylim(0,0.18)+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_fill_manual(values = as.character(jcolors::jcolors(palette = "pal8"))[1:6],breaks = c("MERINGUE","SpaGCN","SPARK","SPARK-X","SpatialDE","SpaGFT"))

cowplot::plot_grid(p_HE,p_GSE152506_ST,p.MOB_ST,p.puck,p.S2_5_Visium,p.MOB_ST.1,p.MOB_ST.11)


