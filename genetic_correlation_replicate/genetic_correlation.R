library(data.table)

setwd("/projectnb/bs859/students/mindyt5/GWAS-Analysis-Endometriosis/genetic_correlation_replicate")

# Read the file
dat <- fread("zcat GCST90478722_dysmenorrhea.tsv.gz", sep = "\t", header = TRUE, na.strings = c("NA", "#NA", ""))

head(dat)

# Force-convert p_value to numeric (non-numeric becomes NA)
dat[, p_value := as.numeric(p_value)]

# Check problematic rows (optional)
print(dat[is.na(p_value), .N])  # Count NA p-values



sapply(dat, class)  # Ensure p_value is "numeric"
fwrite(dat, "GCST90478722_dysmenorrhea_clean.tsv", sep = "\t")
system("gzip GCST90478722_dysmenorrhea_clean.tsv")  # Compress



library(data.table)
dat <- fread("zcat GCST90478722_dysmenorrhea_clean.tsv.gz | sed 's/\r//g'", sep = "\t", header = TRUE)
fwrite(dat, "GCST90478722_fixed.tsv", sep = "\t", quote = FALSE, na = "NA")
system("gzip -f GCST90478722_fixed.tsv")



### BARPLOT
# Load required library
library(ggplot2)
# Create a data frame with p-values included
df <- data.frame(
  Trait = c("Asthma", "Uterine Fibroids", "Migraine", "Osteoarthritis", "Dysmenorrhea", "Menorrhagia"),
  rg = c(0.0047, 0.5194, 0.1898, 0.0294, 0.4105, 0.6446),
  se = c(0.0592, 0.0879, 0.0602, 0.0557, 0.2834, 0.0547),
  p = c(9.3649e-01, 3.4764e-09, 1.6201e-03, 5.9717e-01, 1.4754e-01, 4.2567e-32)
)

# Format p-values for display
df$p_label <- formatC(df$p, format = "e", digits = 2)

# Plot with p-values as labels
ggplot(df, aes(x = reorder(Trait, rg), y = rg)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = rg - se, ymax = rg + se), width = 0.2) +
  geom_text(aes(label = paste0("p = ", p_label)), 
            vjust = -0.8, size = 5, color = "black") +
  coord_flip() +
  theme_minimal(base_size = 16) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(face = "bold", color = "black", size = 18),
    plot.title = element_text(face = "bold", color = "black", size = 20, hjust = 0.5)
  ) +
  labs(
    title = "Genetic Correlation with Endometriosis",
    x = "Trait",
    y = "Genetic Correlation (rg)"
  )


# Load required library
library(ggplot2)

# Create a data frame with p-values included
df <- data.frame(
  Trait = c("Asthma", "Uterine Fibroids", "Migraine", "Osteoarthritis", "Dysmenorrhea", "Menorrhagia"),
  rg = c(0.0047, 0.5194, 0.1898, 0.0294, 0.4105, 0.6446),
  se = c(0.0592, 0.0879, 0.0602, 0.0557, 0.2834, 0.0547),
  p = c(9.3649e-01, 3.4764e-09, 1.6201e-03, 5.9717e-01, 1.4754e-01, 4.2567e-32)
)

# Format p-values for display
df$p_label <- formatC(df$p, format = "e", digits = 2)

# Plot with p-values as labels
ggplot(df, aes(x = reorder(Trait, rg), y = rg)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = rg - se, ymax = rg + se), width = 0.2) +
  geom_text(aes(label = paste0("p = ", p_label)), 
            vjust = -3.5, size = 5, color = "red") +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(face = "bold", color = "black", size = 14),
    plot.title = element_text(face = "bold", color = "black", size = 16, hjust = 0.5)
  ) +
  labs(
    title = "Genetic Correlation with Endometriosis",
    x = "Trait",
    y = "Genetic Correlation (rg)"
  )





