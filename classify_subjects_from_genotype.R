# setwd(my_dir)

#genotype_table is the name of input txt file.
genotype_table = "brca_rad51_germline_genotype_table.txt"
library(data.table)

# data table is melted: RecordID, Sample, Variable, Value
melted <- fread(genotype_table)

# Each record has CHROM, POS, REF, ALT, GQ, GT fields for each sample, etc., in different rows
dt = data.table(RecordID = unique(melted$RecordID))

# Get the CHROM for each record in the melted table, then merge into summary dt
dt = dt[melted[Variable == "CHROM", .(chr = Value), by = "RecordID"], on = "RecordID"]

# Repeat with other values
dt = dt[melted[Variable == "POS", .(pos = Value), by = "RecordID"], on = "RecordID"]
dt = dt[melted[Variable == "REF", .(ref = Value), by = "RecordID"], on = "RecordID"]
dt = dt[melted[Variable == "ALT", .(alt = Value), by = "RecordID"], on = "RecordID"]

# To-do: What is the "Pooled" sample? (not making a special column for it for now)
# Get several GQ statistics for each record
dt = dt[melted[Variable == "GQ" & ! is.na(Value), .(prop_GQ_equal_99 = mean(Value == 99),
                                                    prop_GQ_over_40 = mean(Value > 40),
                                                    prop_GQ_over_20 = mean(Value > 20),
                                                    GQ_90_quantile = quantile(as.numeric(Value), probs = .90),
                                                    GQ_10_quantile = quantile(as.numeric(Value), probs = .10),
                                                    GQ_mean = mean(as.numeric(Value)),
                                                    GQ_median = median(as.numeric(Value))), 
                                                    by = "RecordID"], on = "RecordID"]

# Saving table (looks like a threshold of 40 will keep most data; 20 is also a common threshold)
#fwrite(dt, file = "brca_rad51_site_gq_summary.txt", sep = "\t")
tcga_MAF=fread("TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf")
tcga_MAF$tcgaID=substr(tcga_MAF$Matched_Norm_Sample_Barcode,1,16)
## To-do: Make the input argument rsID rather than record number (possibly easiest way is to add rsIDs to GQ summary table)
## This function takes the "melted" genotype table has already been read in with fread, and the RecordID
classify_samples = function(melted, record_num, gq_cutoff = 40) {
  # get ref and alt alleles
#  all_chr = melted[Variable == "CHROM", .(chr = Value), by = RecordID]
#  all_pos = melted[Variable == "POS", .(pos = Value), by = RecordID]
#  all_loci = all_chr[all_pos, on = "RecordID"]
#  record_num = all_loci[chr == requested_chr & pos == requested_pos, RecordID]
  ref_allele = melted[RecordID == record_num & Variable == "REF", Value]
  alt_allele = melted[RecordID == record_num & Variable == "ALT", Value]
  
  # get current record's GT and GQ values for each sample
  gt_data = melted[Variable == "GT" & RecordID == record_num][, .(Sample, GT = Value)]
  gq_data = melted[Variable == "GQ" & RecordID == record_num][, .(Sample, GQ = Value)]
  
  # merge together into one genotype table that has one row per sample
  genotypes = gt_data[gq_data, on = "Sample"]
  
  # convert table's Sample field to TCGA sample ID
  genotypes[, tcga_name := sapply(Sample, function(x) paste0("TCGA-", paste(strsplit(x, "-")[[1]][c(2:4)], collapse = "-")))]
  # six TCGA samples had 2 BAM files instead of one, so each has an extra row here
  # if the genotypes disagree, mark as ambiguous; otherwise, take the record with higher GQ
  ambiguous_samples = genotypes[, .(distinct_gt = length(unique(GT))), by = "tcga_name"][distinct_gt > 1, tcga_name]
  if (length(ambiguous_samples) > 0) {
    genotypes = genotypes[! tcga_name %in% ambiguous_samples]
  }
  # sort table by GQ descending and take just the first associated with each TCGA ID 
  tcga_samples = unique(genotypes$tcga_name)
  genotypes = genotypes[order(GQ, decreasing = T)]
  setkey(genotypes, "tcga_name")
  genotypes = genotypes[tcga_samples, mult = "first"]
  
  # remove samples with low genotype quality
  low_quality = genotypes[GQ < gq_cutoff, tcga_name]
  genotypes = genotypes[! tcga_name %in% low_quality]

  # split GT into individual alleles
  # assume all biallelic SNPs (a couple aren't, but we won't use these)
  # interesting: a very small number of records use | instead of / to separate GT; these are apparently
  # phased calls (that is, the chromosome of each allele has been identified)
  genotypes[, c("allele1", "allele2") := tstrsplit(GT, split = '[/|]')]
  
  homozygous_ref = genotypes[allele1 == ref_allele & allele2 == ref_allele, tcga_name]
  homozygous_alt = genotypes[allele1 == alt_allele & allele2 == alt_allele, tcga_name]
  heterozygous = genotypes[(allele1 == ref_allele & allele2 == alt_allele) | (allele1 == alt_allele & allele2 == ref_allele), tcga_name]
  
  called_samples = c(homozygous_ref, homozygous_alt, heterozygous)
  bad_gt = setdiff(genotypes$tcga_name, called_samples)
  if(length(bad_gt) > 0) {
    warning(sprintf("%d samples had non-biallelic or unparseable GT values for record %d: %s.", length(bad_gt), record_num, paste(bad_gt, collapse = ", ")))
  }
  return(list(homozygous_ref = homozygous_ref, homozygous_alt = homozygous_alt, heterozygous = heterozygous, 
              low_quality = low_quality, ambiguous = ambiguous_samples, bad_genotype = bad_gt))
}

# Example of amibiguous site, corresponding to TCGA-A7-A26E-10A
# melted[RecordID == 21 & grepl('A7-A26E-10A', melted$Sample) == T]


genotype_summary=list()
for(i in 1:28){
  genotype_summary[[i]] = classify_samples(melted, record_num=i, gq_cutoff = 40)
}
saveRDS(genotype_summary,"genotype_summary.rds")
