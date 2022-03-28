# Load in required packages
library("affy")
library("affyPLM")
library("sva")
library("AnnotationDbi")
library("hgu133plus2.db")
library("ggplot2")

# Set working directory to CEL file locations
setwd("/projectnb/bf528/users/hedgehog_2022/project_1/sample_symlink")

# Save list of files as character vector
files <- list.files()

# Read CEL files
data <- ReadAffy(filenames = files)

# Create Expression Set
Eset <- rma(data)

# Create Pset
Pset <- fitPLM(data, normalize=TRUE, background=TRUE)

# Compute RLE
RLE_data <- data.frame(t(RLE(Pset,type="stats")))

# Check for RLE outliers (abs(RLE) > 0.1)
RLE_outliers <- RLE_data[which(abs(RLE_data$median)>0.1),]

# Parameter that allows two plots to be diplayed on the same image
par(mfrow=c(1,2))

# Create Histogram from the median RLE values
hist(RLE_data$median,
     main = "Median RLE by Sample",
     xlab = "Median RLE",
     breaks = 10,
     col = 'orangered')
mtext("(A)", font = 2, side = 3, adj = -0.35)

# Compute NUSE
NUSE_data <- data.frame(t(NUSE(Pset,type="stats")))

# Check for NUSE outliers
NUSE_outliers <- NUSE_data[which(NUSE_data$median>1.05),]

# Create Histogram from the median NUSE values
hist(NUSE_data$median,
     main = "Median NUSE by Sample",
     xlab = "Median NUSE",
     breaks = 10,
     col = 'skyblue3')
mtext("(B)", font = 2, side = 3, adj = -0.35)

# Set variables for ComBat
edata <- exprs(Eset)
# Read in metadata table
metadata <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
batch <- metadata$normalizationcombatbatch
mod <-  as.factor(metadata$normalizationcombatmod)
# Run ComBat with parameters
final <- ComBat(edata, batch = batch, mod = mod)

# Write expression data to
write.csv(final, "/projectnb/bf528/users/hedgehog_2022/project_1/outputs/hedgehog_intensity.csv",
          row.names = TRUE)

# Scale data, transposing pre scale and post scale
scaled_final <- t(scale(t(final))) %>%
  # Perform PCA on prescaled data
  PCA <- prcomp(center = FALSE, scale=FALSE, rank = 200)
# Save summary of PCA data
PCA_summary <- summary(PCA)
# Save explained variance values as percentages for PC1 and PC2
PC1_EV <- paste(round(test$importance[2,1], digits = 2) * 100, "%", sep = "")
PC2_EV <- paste(round(test$importance[2,2], digits = 2) * 100, "%", sep = "")

# Create subset data frame from the rotation data
PCA_rotation <- as.data.frame(PCA$rotation)
# Plot PCA
ggplot(PCA_rotation, aes(x = PC1, y = PC2, color = metadata$cit.coloncancermolecularsubtype)) +
  geom_point() +
  theme_light() +
  labs(title = "PCA of Cancer Subtypes in Expression Data",
       x = paste("PC1:", PC1_EV, "Explained Variance"),
       y = paste("PC2:", PC2_EV, "Explained Variance"),
       color = "Cancer Subtype") +
  theme(legend.position = "bottom")