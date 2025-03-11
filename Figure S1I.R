library(Seurat)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(anndata)
library(dplyr)
library(tidyr)
library(janitor)
library(stringr)
library(alluvial)
library(ggalluvial)

# Converting .rds to make compatible with map my results
LS.integrated<-readRDS(file = "LS_integrated.rds")
DefaultAssay(LS.integrated) <- "RNA"

LS.integrated[['integrated']] <- NULL
LS.integrated[['SCT']] <- NULL

SaveH5Seurat(LS.integrated, filename = "LSintegrated2.h5Seurat")
Convert("LSintegrated2.h5Seurat", dest = "h5ad")

mapping <- read.csv("correlation_mapping_results.csv",comment.char="#")
head(data.frame(mapping))
data.frame(Cell_counts=head(sort(table(mapping$class_name),decreasing=T),8))

dataQC_h5ad <- read_h5ad('LSintegrated.h5ad')
dataQC <- t(as.matrix(dataQC_h5ad$X))
rownames(dataQC) <- rownames(dataQC_h5ad$var)
colnames(dataQC) <- rownames(dataQC_h5ad$obs)

# Loading and formatting data
LS.integrated <- NormalizeData(LS.integrated)
table(Idents(LS.integrated))
new.ident <- c("Gaba1","Gaba2","Gaba3","Gaba4","Gaba5","Gaba6","Gaba7","Glu1","Gaba8","Gaba9","Gaba10","Glu2","Gaba11","Gaba12")
names(x = new.ident) <- levels(x =LS.integrated)
LS.integrated<- RenameIdents(object =LS.integrated, new.ident)
table(Idents(LS.integrated))

celltype <- data.frame(Idents(LS.integrated))
LS.integrated@meta.data$celltype <- celltype

# read in mapping results
mapping <- read.csv("Hierarchical_mapping_results.csv",comment.char="#")
head(data.frame(mapping))

sub_mapped <- mapping[c('cell_id','subclass_name','class_name')]
celltype <- as.character(Idents(LS.integrated))
tot <- cbind(sub_mapped, celltype)

allu <- tot %>% 
  group_by(class_name, subclass_name, celltype,) %>%   # grouping
  summarise(Freq = n())

allu <- allu[order(allu$Freq, decreasing=TRUE),]
allu$subclass_name <- factor(allu$subclass_name, levels=unique(allu$subclass_name))

allu <- allu[order(allu$Freq, decreasing=TRUE),]
allu$class_name <- factor(allu$class_name, levels=unique(allu$class_name))


allu <- allu[order(allu$Freq, decreasing=TRUE),]
allu$subclass_name <- factor(allu$subclass_name, levels=unique(allu$subclass_name))

allu <- allu[order(allu$Freq, decreasing=TRUE),]
allu$celltype <- factor(allu$celltype, levels=new.ident)

ggplot(as.data.frame(allu),
       aes(y = Freq, axis1 = celltype, axis2 = class_name, axis3=subclass_name)) +
  geom_alluvium(aes(fill = Freq), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Cell type", "Cell class"), expand = c(.05, .05)) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Correlation")
