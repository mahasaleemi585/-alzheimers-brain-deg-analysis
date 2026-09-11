# ===== Phase 2: Load and label the data =====

# Load the expression matrix
expr_data <- read.delim("GSE104704_RNA-Seq_Table.txt.gz", header = TRUE, row.names = 1)
dim(expr_data)
head(expr_data[, 1:5])

# Separate the 30 real samples from the 3 extra comparison columns
sample_cols <- expr_data[, 1:30]
sample_names <- colnames(sample_cols)

# Match each column to its correct group using a lookup table
# (built from GEO sample metadata, matched by sample ID rather than column order)
group_lookup <- c(
  "2-12A"="Young", "3-17T"="Young", "4-13A"="Young", "5-18T"="Young",
  "6-14A"="Young", "7-19T"="Young", "8-15A"="Young", "9-16A"="Young",
  "10-8A"="Old", "11-10T"="Old", "12-6A"="Old", "13-11T"="Old",
  "14-7A"="Old", "15-13T"="Old", "16-14T"="Old", "17-9A"="Old",
  "18-10A"="Old", "19-11A"="Old",
  "20-1T"="AD", "21-1A"="AD", "22-2T"="AD", "23-2A"="AD",
  "24-3T"="AD", "25-5T"="AD", "26-3A"="AD", "27-5A"="AD",
  "28-8T"="AD", "29-6T"="AD", "30-9T"="AD", "31-7T"="AD"
)

# Recover the original sample IDs (R renamed them slightly) and match to groups
sample_ids <- sub("^X", "", sample_names)
sample_ids <- sub("^([0-9]+\\.[0-9]+[A-Z])\\..*$", "\\1", sample_ids)
sample_ids <- sub("\\.", "-", sample_ids)
group <- group_lookup[sample_ids]

# Sanity check - verify each sample name matches its assigned group
data.frame(sample_names, sample_ids, group)
