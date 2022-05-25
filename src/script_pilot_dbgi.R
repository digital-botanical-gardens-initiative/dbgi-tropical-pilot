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
metabo_matrix_format$attribute <- as.factor(metadata_format$attribute_2)
  
  

model_1 = randomForest(x = metabo_matrix_format[, colnames(metabo_matrix_format) != "attribute"],
  y = metabo_matrix_format$attribute, importance = TRUE,ntree=100)
  
ozone.rp <- rfPermute(x = metabo_matrix_format[, colnames(metabo_matrix_format) != "attribute"],
                      y = metabo_matrix_format$attribute, na.action = na.omit, ntree = 100, num.rep = 5)

# plotImportance(ozone.rp, scale = TRUE,size = 3)


importance =  randomForest::importance(model_1)
varImportance = data.frame(Variables = row.names(importance),
                           Importance =round(importance[, "MeanDecreaseAccuracy"],2))

rankImportance= varImportance %>% mutate(Rank=paste("#",dplyr::dense_rank(dplyr::desc(Importance))))
rankImportance2 <- rankImportance[rankImportance$Importance > 1.1 ,]#| rankImportance$Importance < -0.5

var_imp<- ggplot(rankImportance2,aes(x=reorder(Variables,Importance),y=Importance,fill=Importance))+ 
  geom_bar(stat="identity") + 
  geom_text(aes(x = Variables, y = 0.5, label = Rank),hjust=0, vjust=0.55, size = 4, colour = "white") +
  labs(x = "Variables") +
  coord_flip() + 
  theme_classic()


setwd("G:/My Drive/taf/git_repository/dbgi-tropical-pilot/docs/results")
sink("table.txt")
print(summary(ozone.rp))
sink() 

Rr_perm <- readLines("G:/My Drive/taf/git_repository/dbgi-tropical-pilot/docs/results/table.txt")
Rr_perm <- Rr_perm[-length(Rr_perm)]



pdf("G:/My Drive/taf/git_repository/dbgi-tropical-pilot/docs/results/result_pilot_dbgi.pdf",width=15, height=15)
title <- "result_pilot_dbgi"
grid::grid.text(title,x = (0.5), y = (0.6))
circ
nmds_plot
var_imp
par(mar = c(0.1,1,1,0.1))
plot(x=c(1,21),y=c(1,29.7),type="n", axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
liner <- rev(seq(1,29,0.4)) ## interline
for (i in c(1:length(Rr_perm))) {
  text(x=2, y=liner[i], labels=Rr_perm[i],cex = 0.5,pos = 4) 
}
dev.off()

