# Clusters-and-Dimensional-Reduction-of-Seurat-Objects

# Results


# Load necessary libraries for data manipulation and visualization
library(Seurat)
library(scCustFx)
library(pheatmap)   # Heatmap visualization
library(ggplot2)    # General plotting



knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
knitr::opts_knit$set(root.dir = getwd())  # Set working directory for knitting


## Load Seurat Object(s)


Seurat_Object_1 = readRDS("File_Path_to_Seurat_Object")

head(Seurat_Object_1@meta.data)


DimPlot(Seurat_Object_1)
ImageDimPlot(Seurat_Object_1)


VariableFeatures(Seurat_Object_1)[1:50]



# Seurat_Object_1_new <- ConvertFOVToSpatialSeurat(Seurat_Object_1)



# Idents(Seurat_Object_1_new) = "ClusterNames_0.4"
# Idents(Seurat_Object_1) = "ClusterNames_0.4"
# SpatialDimPlot(Seurat_Object_1_new) | DimPlot(Seurat_Object_1)



p1 <- FeaturePlot(Seurat_Object_1, features = "Biomarker")

# The cell types will be the second list of cell types in the metadata. The other lists with cell types do not work. 

p2 <- DimPlot(Seurat_Object_1, group.by = "Cell type in metadata")

p1+p2

FeaturePlot(Seurat_Object_1,
            features="PC_2",
            raster = F,
            raster.dpi = c(2000, 2000),
            pt.size = .1,
            max.cutoff = 'q95',
            min.cutoff = 'q01',
            order = T) + 
  # coord_flip()  + scale_y_reverse() +
  theme_classic(base_size = 14) + 
  # NoLegend()+
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank() #,plot.title = element_blank()
  )   &
  ggplot2::scale_colour_gradientn(colours = c("navy", "darkblue", "dodgerblue", "gold", "red", "maroon"))




FeatureScatter(Seurat_Object_1, feature1 = "nef", "env") |
  FeatureScatter(Seurat_Object_1, feature1 = "gag", "pol")


sort( abs(Seurat_Object_1@reductions$pca@feature.loadings[,1]), decreasing = T)[1:20]
sort( abs(Seurat_Object_1@reductions$pca@feature.loadings[,2]), decreasing = T)[1:20]
