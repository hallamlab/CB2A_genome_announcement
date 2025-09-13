# Install required packages
install.packages("BiocManager")
BiocManager::install(c("genbankr", "Biostrings"))
install.packages("dplyr")


setwd("//wsl.localhost/Ubuntu/home/beth1313")

library(genbankr)
library(Biostrings)
library(dplyr)

# File paths
genome1_file <- "CB2.gb"
genome2_file <- "CB15.gb"

# Parse GenBank files and extract CDS features
genome1 <- readGenBank(genome1_file)
genome2 <- readGenBank(genome2_file)

# Function to extract CDS features and sequences
extract_cds <- function(genome) {
  # Extract CDS features
  cds_features <- cds(genome)
  
  # Extract sequences for CDS regions
  genome_seq <- unlist(genome@sequence)  # Extract full genome sequence
  cds_sequences <- lapply(seq_along(cds_features), function(i) {
    region <- cds_features[i]
    extractAt(genome_seq, ranges(region))
  })
  
  # Create a data frame with extracted information
  data.frame(
    start = start(cds_features),
    end = end(cds_features),
    strand = strand(cds_features),
    gene = ifelse(is.na(mcols(cds_features)$gene), "Unknown", mcols(cds_features)$gene),
    product = ifelse(is.na(mcols(cds_features)$product), "Unknown", mcols(cds_features)$product),
    sequence = sapply(cds_sequences, toString)
  )
}

genome1_cds <- extract_cds(genome1)
genome2_cds <- extract_cds(genome2)


# Create a new table with FASTA headers and sequences
genome1_cds_2col <- genome1_cds %>%
  mutate(gene = paste0("gene", row_number())) %>%  # Number the genes as gene1, gene2, etc.
  select(gene, sequence)  # Keep only the FASTA header and sequence columns

genome2_cds_2col <- genome2_cds %>%
  mutate(gene = paste0("gene", row_number())) %>%  # Number the genes as gene1, gene2, etc.
  select(gene, sequence)  # Keep only the FASTA header and sequence columns

# Write FASTA file
write_fasta <- function(df, output_file) {
  sink(output_file)
  for (i in seq_len(nrow(df))) {
    cat(paste0(">", df$gene[i], "\n"))  # FASTA header
    cat(paste0(df$sequence[i], "\n"))  # Sequence
  }
  sink()
}

# Save to a FASTA file
write_fasta(genome1_cds_2col, "genome1_cds.fa")
write_fasta(genome2_cds_2col, "genome2_cds.fa")

# Create BLAST databases for both genomes
system("makeblastdb -in genome1_cds.fa -dbtype nucl -out genome1_db")
system("makeblastdb -in genome2_cds.fa -dbtype nucl -out genome2_db")
# Query genome 1 CDS against genome 2
system("blastn -query genome1_cds.fa -db genome2_db -out genome1_vs_genome2.out -outfmt 6")
# Query genome 2 CDS against genome 1
system("blastn -query genome2_cds.fa -db genome1_db -out genome2_vs_genome1.out -outfmt 6")

library(data.table)

# Load BLAST results
genome1_vs_genome2 <- fread("genome1_vs_genome2.out", header = FALSE, sep = "\t")
genome2_vs_genome1 <- fread("genome2_vs_genome1.out", header = FALSE, sep = "\t")

# BLAST column names
colnames(genome1_vs_genome2) <- colnames(genome2_vs_genome1) <- c(
  "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
  "qstart", "qend", "sstart", "send", "evalue", "bitscore"
)


# Extract gene IDs from FASTA headers
gene_ids_genome1 <- genome1_cds_2col$gene
gene_ids_genome2 <- genome2_cds_2col$gene

# Extract matched query and subject gene IDs
matched_query_ids <- unique(genome1_vs_genome2$qseqid)
matched_subject_ids <- unique(genome1_vs_genome2$sseqid)

# Find genes in genome1 not matched in genome2
unique_genome1 <- setdiff(gene_ids_genome1, matched_query_ids)

# Find genes in genome2 not matched in genome1
unique_genome2 <- setdiff(gene_ids_genome2, matched_subject_ids)

# View unique genes
print(unique_genome1)
print(unique_genome2)

# Filter unique sequences from the original data frames
unique_sequences_genome1 <- genome1_cds_2col %>% filter(gene %in% unique_genome1)
unique_sequences_genome2 <- genome2_cds_2col %>% filter(gene %in% unique_genome2)

# Save unique sequences to new FASTA files
write_fasta(unique_sequences_genome1, "unique_genome1.fa")
write_fasta(unique_sequences_genome2, "unique_genome2.fa")


