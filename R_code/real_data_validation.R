setwd("e:/Dropbox/my published paper/SVG_jixin/results/03022022-read-data/SpecDomain/")
library(Seurat)
library(scales)
sample_name <- "151673"
spec_range <- c(2:50)
file_name_spectrum <- paste0("e:/Dropbox/my published paper/SVG_jixin/results/03022022-read-data/SpecDomain/",sample_name, "_spec_low400_high100.csv")
file_name_gene_siml <- paste0("e:/Data and Code/Source data/original/",sample_name)
expression_df <- Load10X_Spatial(data.dir = file_name_gene_siml, filename = "filtered_feature_bc_matrix.h5")
spectrum_df <- read.csv(file_name_spectrum,header = T,row.names = 1)



# GPS1

gene_name <- "MOBP"
#
x_y <- expression_df@images$slice1@coordinates$row


plot_dim_df <- data.frame(row.names = colnames(expression_df),
                          x = expression_df@images$slice1@coordinates$row,
                          y = expression_df@images$slice1@coordinates$col,
                          Gene = expression_df@assays$Spatial@data[gene_name,]
)

#
library(ggplot2)
library(cowplot)
p_SVG <- ggplot(plot_dim_df, aes(x = x, y = y, color = log1p(Gene))) + 
  geom_point(size  = 2) + 
  scale_color_gradientn(colours = viridis::viridis(10)) + theme_void()+theme(legend.position = "none",axis.title=element_blank(),
                                                              axis.text=element_blank(),
                                                              axis.ticks=element_blank())

SVG_fourier <- spectrum_df[,gene_name] 

max_spec <- max(spectrum_df)
SVG_fourier <- round(c(SVG_fourier/max_spec),2)
min_spec <- min(SVG_fourier)
#SVG_fourier <- (SVG_fourier - min(SVG_fourier))/ (max(SVG_fourier) - min(SVG_fourier))
# softmax normalization
# SVG_fourier <- exp(SVG_fourier)/sum(exp(SVG_fourier))
# softmax <- function(par){
#   n.par <- length(par)
#   par1 <- sort(par, decreasing = TRUE)
#   Lk <- par1[1]
#   for (k in 1:(n.par-1)) {
#     Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
#   }
#   val <- exp(par - Lk)
#   return(val)
# }
# SVG_fourier <- softmax(SVG_fourier)


gg_df <- as.data.frame(cbind(ID = 1:length(as.numeric(SVG_fourier)[spec_range]), value = as.numeric(SVG_fourier)[spec_range]))
gg_df$colors <-  c(rep("A",15),
                   rep("B",nrow(gg_df)-15))

barplot <- ggplot(gg_df, aes( x = ID, y = value, fill = colors)) +
  scale_fill_manual(values =c("#c81e1e","#385994") ,breaks = c("A","B"))+ 
  geom_col() + theme_classic()+ theme(axis.text.y = element_text(size = 18,face = "bold",vjust = 0),
                                      axis.title.y=element_blank(),
                                      axis.title.x=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.ticks.x=element_blank())+
  scale_x_continuous(limits = c(1,max(gg_df$ID)), expand = c(0, 0))+
  scale_y_continuous(limits = c(0,max(gg_df$value)), expand = expansion(mult = c(0, .1)))+
  coord_cartesian(ylim = c(0, 0.3))+
  guides(fill="none") 


# gg_df$sig_frag <- c(rep("topology",floor(0.05* nrow(gg_df))),
#                       rep("signal",ceiling(0.95* nrow(gg_df))))
# gg_df$sig_frag_col <- gg_df$sig_frag 
# gg_df$sig_frag_col[gg_df$sig_frag_col == "topology"] <- colorRampPalette(c("#FFFFFF","#008AFC"))(n= length(which(gg_df$sig_frag_col == "topology")))
# gg_df$sig_frag_col[gg_df$sig_frag_col == "signal"] <- colorRampPalette(c("#FFFFFF","#FF3A3A"))(n= length(which(gg_df$sig_frag_col == "signal")))

# bar_code <- ggplot(gg_df, aes( x = ID, y = 1, fill = sig_frag)) + 
#   geom_tile()+ 
#   theme_void()+ scale_fill_manual(values = gg_df$sig_frag_col) + theme(legend.position = "none")
bar_code <- ggplot(gg_df, aes( x = ID, y = 1, fill = value)) + 
  geom_tile()+ 
  theme_void()+ scale_fill_gradient(low = "white",high = "black") + theme(legend.position = "none")


library(cowplot)
plot_grid(p_SVG,barplot, bar_code,ncol = 1,rel_heights = c(2,1,1))

