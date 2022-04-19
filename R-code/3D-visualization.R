library(plotly)
Plot_df <- plot_dim_df
Plot_df$x <- round(Plot_df$x,digits = 0)
Plot_df$y <- round(Plot_df$y,digits = 0)
max(Plot_df$x)
max(Plot_df$y)
min(Plot_df$x)
min(Plot_df$y)

Plot_df$x <- Plot_df$x+1
Plot_df$y <- Plot_df$y+1

my.matrix <- matrix(rep(0,30*30),nrow = 30)
for (i in 1:30){
  for (j in 1:30){
    tmp_value <- max(Plot_df$SVG_expression[Plot_df$x==i & Plot_df$y==j])
    my.matrix[i,j] <- tmp_value
  }
}

my.matrix <- ifelse(is.infinite(my.matrix),0,my.matrix)

library(rayimage)
library(Seurat)
x <- Load10X_Spatial(data.dir = "e:/Human_brain/151673/",filename = "filtered_feature_bc_matrix.h5")
x <- NormalizeData(x)
range(x@images$slice1@coordinates$imagecol)
range(x@images$slice1@coordinates$imagerow)
x_y <- data.frame(x = x@images$slice1@coordinates$imagerow+1, 
                  y = x@images$slice1@coordinates$imagecol+1, 
                  gene = x@assays$Spatial@counts["GPS1",])

ggplot(x_y,aes(x = x, y=y, color = gene))+
  geom_point()+ 
  scale_color_gradientn(colors  = viridis::magma(10))
  

r_spot <- floor(x@images$slice1@scale.factors$fiducial/3)

gene_count_total <- sum(x_y$gene)
gene_count_total <- gene_count_total + length(which(x_y$gene == 0))

x_y_simulation <- matrix(rep(NA, gene_count_total*2),nrow = gene_count_total)
colnames(x_y_simulation) <- c("x","y")

i = 1
j = 1
while (i <= nrow(x_y_simulation)) {
  tmp_counts <- x_y$gene[i]
  if(tmp_counts == 0){
    tmp_x <- x_y$x[i]
    tmp_y <- x_y$y[i]
    new_tmp_x <- tmp_x 
    new_tmp_y <- tmp_y 
    
    x_y_simulation[c(j:(j+1)),1] <- new_tmp_x
    x_y_simulation[c(j:(j+1)),2] <- new_tmp_y
    
    i=i+1
    j = j+1
  }else{
    tmp_x <- x_y$x[i]
    tmp_y <- x_y$y[i]
    
    theta_factor <- sample(seq(from = 0,to = 2*pi,length.out = 10*tmp_counts),
                           tmp_counts)
    new_tmp_x <- tmp_x + floor(1.2*r_spot*cos(theta_factor))
    new_tmp_y <- tmp_y + floor(1.2*r_spot*sin(theta_factor))
    
    x_y_simulation[c(j:(j+tmp_counts-1)),1] <- new_tmp_x
    x_y_simulation[c(j:(j+tmp_counts-1)),2] <- new_tmp_y
    i=i+1
    j = j+tmp_counts
    }
}
library(rayshader)
library(ggplot2)
library(tidyverse)
x_y_simulation <- as.data.frame(x_y_simulation)

gg = ggplot(x_y_simulation) +
  stat_density_2d(aes(x = x, y = y, fill = stat(nlevel)),
                  geom = "polygon",
                  n = 50,bins = 30,contour = T) +
  scale_fill_viridis_c(option = "A") + theme_void()+ theme(legend.position="none")

gg  

# gg = ggplot(x_y_simulation,aes(x = x,y=y)) +
#   geom_hex(bins = 70) +
#   scale_fill_continuous(type = "viridis") +
#   theme_bw()
# gg
plot_gg(gg, width = 5, height = 5, multicore = TRUE, scale = 250, 
        zoom = 0.7, theta = 10, phi = 30, windowsize = c(800, 800))

# plot network










