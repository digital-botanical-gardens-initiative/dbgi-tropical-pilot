library(vegan)
library(ape)
library(ggtree)
library(wesanderson)
library(ggalt)
library(randomForest)
library(rfPermute)
library(igraph)
library(rcdk)
library(dplyr)
library(plotly)

######################################################################
####################### load data

metabo_matrix <- read.csv("G:/My Drive/taf/git_repository/gnps_job/dbgi_pilot/quantification_table/quantification_table-00000.csv")
metadata <- read.csv("G:/My Drive/taf/git_repository/gnps_job/dbgi_pilot/metadata_table/metadata_table-00000.tsv",sep="\t")
annotation <- read.csv("G:/My Drive/taf/git_repository/dbgi-tropical-pilot/docs/results/met_annot_enhancer/dbgi-pilot-tropical/dbgi-pilot-tropical_spectral_match_results_repond.tsv", sep ="\t")
metadata_site <- read.csv("./docs/MAPP_batch_00018/metadata_site.csv", sep =",")


filename <- metadata$filename
sample_name <- metadata$ATTRIBUTE_query_otol_species
discriminating_factor_1 <- metadata$ATTRIBUTE_sample_prep_mode
attribute_1 <- metadata$ATTRIBUTE_query_otol_family
attribute_2 <- metadata$ATTRIBUTE_sample_prep_mode

metadata_format <- data.frame(filename,sample_name,attribute_1,attribute_2,discriminating_factor_1)


metabo_matrix_format <- t(metabo_matrix[,grep("DBGI",colnames(metabo_matrix))])
sample_id <- rownames(metabo_matrix_format)
sample_id <-gsub("X","", gsub(".Peak.area","",sample_id))
feature_id <- metabo_matrix$row.ID 
colnames(metabo_matrix_format) <- feature_id

merge_file_sample_name<- data.frame(gsub("X","", gsub(".Peak.area","",row.names(metabo_matrix_format))))
colnames(merge_file_sample_name)[1] <- "filename"


name_format <-  merge(merge_file_sample_name,metadata_format, by="filename" )
rownames(metabo_matrix_format) <- name_format$sample_name


metabo_matrix_format <- metabo_matrix_format[row.names(metabo_matrix_format) %in% metadata_format$sample_name,]
metadata_format <- metadata_format[ metadata_format$filename %in% sample_id,]


###################
###3 remove_qc

dist_metabo<- vegdist(metabo_matrix_format,method="bray") # method="man" # is a bit better
hclust_metabo<- hclust(dist_metabo, method = "complete")

dendro_metabo <- as.phylo(hclust_metabo)
#### remove_QC 


plot(dendro_metabo)


############ circular plot 

library(scales)
g2 <- split(as.factor(metadata_format$sample_name), as.factor(metadata_format$attribute_1))
tree_plot2 <- groupOTU(dendro_metabo, g2,group_name = "grouper")

cols <- sample(c(wes_palette("Cavalcanti1"),wes_palette("GrandBudapest1")))
numcol  <-length(unique(metadata_format$attribute_1))
colfunc <- colorRampPalette(c(cols))
cols<- colfunc(numcol+1)


circ <- ggtree(tree_plot2, aes(color=grouper,fill=grouper),layout='circular',size=1)+ #
  geom_tiplab(size=4, offset=0.08)  +
  scale_colour_manual(values=cols ) + theme(legend.position="none")

#df <- data.frame(metadata_noqc$attribute_1)

#rownames(df) <- paste(c(1:length(tree_plot2$tip.label)),tree_plot2$tip.label)


#ggheat_plot<-gheatmap(circ, df[, "metadata_noqc.attribute_1", drop=F], offset=0, width=0.1,colnames = FALSE,
#                     colnames_angle=90, colnames_offset_y = 0) + scale_fill_manual(values=cols,breaks=0:2)

#ggheat_plot

####3 meta mds


nmds <- metaMDS(metabo_matrix_format,distance="gower")
scores <- vegan::scores(nmds)
sample_name_nmds <-  row.names(scores)
scores<- data.frame(scores)
scores$sample_name_nmds  <- sample_name_nmds

nmds_plot<- scores %>%
  left_join(metadata_format, by = c("sample_name_nmds" = "sample_name"))%>%
  data.frame()%>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_encircle(aes(group = attribute_2, color = attribute_2, fill = attribute_2), alpha = 0.4,s_shape=0.8) +
  geom_point(aes(color = attribute_2)) +
  #annotate("text", x = -0.5, y = 0.3, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw() + scale_fill_manual(values=c("darkgreen","blue","darkred")) + scale_color_manual(values=c("darkgreen","blue","darkred"))

nmds_plot

### data supervizer

metabo_matrix_format <- data.frame(metabo_matrix_format)
metabo_matrix_format$attribute <- metadata_format$
  
  
  ##################### plsda 
  
  liney <- MASS::lda(attribute_2~., data = metabo_matrix_format)
result <- predict(liney, metabo_matrix_format)
lhat <- 1 - sum(result$class == metabo_matrix_format$attribute_2)/length(metabo_matrix_format$attribute_2)

data <- data.frame(x1=result$x[,1], y=metabo_matrix_format$attribute_2)
data$y <- factor(data$y)
LDA <- ggplot(data, aes(x=x1, fill=y)) +
  geom_density(adjust=5, alpha=0.6) +
  xlab("x1") +
  ylab("Density") +
  ggtitle(sprintf("PLS-LDA, L = %.2f", lhat)) + scale_fill_manual(values=cols2)+
  theme_bw()

