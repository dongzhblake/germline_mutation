library(cancereffectsizeR)
library(data.table)

setwd("/Users/dongzihan/Desktop/germline-somatic mutation")
a=fread("TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf")
a$sampleID=substr(a$Matched_Norm_Sample_Barcode,1,16)
analysis = CESAnalysis(genome = "hg19")
analysis = load_maf(analysis, maf = "TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf", chain_file = "hg38ToHg19.over.chain")
sig_to_remove = suggest_cosmic_v3_signatures_to_remove("BRCA", treatment_naive = T, quiet = T)
sig_to_remove =c(sig_to_remove,"SBS84","SBS85")
analysis = trinuc_mutation_rates(cesa=analysis,signatures_to_remove = sig_to_remove)
saveRDS(analysis,"tri_nuc.rds")
# save ces result for trinuc_rates, so that we don't need to re-compute it for every site
analysis = set_trinuc_rates(analysis,trinuc_rates)
analysis = gene_mutation_rates(analysis, covariate_file = "breast_pca")
analysis = annotate_variants(analysis)
analysis = ces_snv(analysis)


mut_sig=data.frame(readRDS("tri_nuc.rds")$mutational_signatures)
mut_sig=mut_sig[mut_sig$group_avg_blended==F,]
rownames(mut_sig)=mut_sig$Unique_Patient_Identifier
mut_sig=mut_sig[,c(-1,-2,-3)]
mut_sig=mut_sig[,which((colSums(mut_sig)!=0))]
rownames(mut_sig)=a[match(rownames(mut_sig),a$Tumor_Sample_Barcode),]$sampleID
mut_sig$group="undefined"

# function to calculate the total absolute effect on mutational signatures of alt VS ref on each site.
cal_total_effect<-function(site){
  samples=list(sub_samples[[site]][1]$homozygous_ref,sub_samples[[site]][2]$homozygous_alt,sub_samples[[site]][3]$heterozygous,c(sub_samples[[site]][2]$homozygous_alt,sub_samples[[site]][3]$heterozygous))
  if(min(length(samples[[1]]),length(samples[[4]]))<=10){
    return(c(0,NA,NA))
  }
  else{
  mut_sig[rownames(mut_sig) %in% samples[[1]],]$group="homo_ref"
  mut_sig[rownames(mut_sig) %in% samples[[4]],]$group="any_alt"
  total_effect=0
  for(i in 1:(ncol(mut_sig)-1)){
    group_ref_mean=mean(mut_sig[mut_sig$group=="homo_ref",i])
    group_alt_mean=mean(mut_sig[mut_sig$group=="any_alt",i])
    total_effect=total_effect+abs(group_ref_mean-group_alt_mean)
  }
  return(c(total_effect,nrow(mut_sig[rownames(mut_sig) %in% samples[[1]],]),nrow(mut_sig[rownames(mut_sig) %in% samples[[4]],])))
}}
q=cal_total_effect(1)
for(i in 2:28){
  q=rbind(q,cal_total_effect(i))
}
colnames(q)=c("total_abs_effect_size","samplesize_ref","samplesize_alt")
rownames(q)= paste0("site",c(1:28))
# site 3: 0.14229560  328  188
# site 7: 0.14586869  279  174
# site 12: 0.13440323  287  200
# site 8, 10, 14 have the most effect on signature mutation processes. 





