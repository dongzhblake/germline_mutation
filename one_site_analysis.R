  library(crayon)
  library(cancereffectsizeR)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  setwd("/Users/dongzihan/Desktop/germline-somatic mutation")
  a=fread("TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf")
  a$sampleID=substr(a$Matched_Norm_Sample_Barcode,1,16)
  sub_samples=readRDS("genotype_summary.rds")
# select which site to analyze
i=1
samples=list(sub_samples[[i]][1]$homozygous_ref,sub_samples[[i]][2]$homozygous_alt,sub_samples[[i]][3]$heterozygous,c(sub_samples[[i]][2]$homozygous_alt,sub_samples[[i]][3]$heterozygous))
a$genotype="undefined"
a[a$sampleID %in% samples[[1]],]$genotype="homo_ref"
a[a$sampleID %in% samples[[4]],]$genotype="any_alt"
a=a[a$genotype != "undefined",]
ces_anal = CESAnalysis(genome = "hg19", 
                       progression_order = c("homo_ref","any_alt"))
analysis <- load_maf(ces_anal, maf = a, progression_col = "genotype",chain_file = "hg38ToHg19.over.chain")
sig_to_remove = suggest_cosmic_v3_signatures_to_remove("BRCA", treatment_naive = T, quiet = T)
sig_to_remove =c(sig_to_remove,"SBS84","SBS85")
analysis = trinuc_mutation_rates(cesa=analysis,signatures_to_remove = sig_to_remove)
analysis = gene_mutation_rates(analysis, covariate_file = "breast_pca")
analysis = annotate_variants(analysis)
analysis = ces_snv(analysis)
saveRDS(analysis,paste0("site_",i,"_SI.rds"))



#tri-nuc mutatation rate heatmap 
kirc_02=analysis=readRDS(paste0("site_",i,"_SI.rds"))
a=analysis@samples
stage1=a[a$progression_name=="homo_ref",]$Unique_Patient_Identifier
stage234=a[a$progression_name=="any_alt",]$Unique_Patient_Identifier
KIRC_trinuc_stage <- 
  kirc_02@trinucleotide_mutation_weights$trinuc_proportion_matrix
tumor_stage1 <- KIRC_trinuc_stage[stage1, ]
tumor_stage234 <- KIRC_trinuc_stage[stage234, ]

tumor_stage1_avg <- 
  apply(tumor_stage1, 2, mean)
tumor_stage234_avg <- 
  apply(tumor_stage234, 2, mean)

tumor_stage1_avg_ordered <- 
  data.frame(average = tumor_stage1_avg,
             mutation = names(tumor_stage1_avg)) %>%
  mutate(upstream = substr(mutation, 1, 1),
         downstream = substr(mutation, 7, 7),
         mut_from = substr(mutation, 3, 3),
         mut_to = substr(mutation, 5, 5),
         mutation_name = paste0(mut_from, "\u2192", mut_to)) %>%
  arrange(downstream, upstream)

tumor_stage234_avg_ordered <- 
  data.frame(average = tumor_stage234_avg,
             mutation = names(tumor_stage1_avg)) %>%
  mutate(upstream = substr(mutation, 1, 1),
         downstream = substr(mutation, 7, 7),
         mut_from = substr(mutation, 3, 3),
         mut_to = substr(mutation, 5, 5),
         mutation_name = paste0(mut_from, "\u2192", mut_to)) %>%
  arrange(downstream, upstream)
trinuc.mutation_data <- tumor_stage1_avg_ordered[, -1]

KIRC_trinuc_stage_heatmap_data <- 
  data.frame(deconstrucSig=trinuc.mutation_data$mutation,
             Upstream=trinuc.mutation_data$upstream,
             Downstream=trinuc.mutation_data$downstream,
             mutated_from=trinuc.mutation_data$mut_from,
             mutated_to=trinuc.mutation_data$mut_to,
             trinuc_Size_3_and_less=tumor_stage1_avg_ordered$average,
             trinuc_Size_3.1_and_up=tumor_stage234_avg_ordered$average) %>%
  mutate(mutation = paste0(mutated_from, "\u2192", mutated_to))

KIRC_trinuc_stage1_heatmap <- 
  ggplot(data=KIRC_trinuc_stage_heatmap_data, aes(x=Downstream, Upstream)) +
  geom_tile(aes(fill = trinuc_Size_3_and_less*100), color="white") + 
  scale_fill_gradient(low="white", high="dark green", name="Percent") +
  facet_grid(.~mutation)+
  geom_text(aes(label = round(trinuc_Size_3_and_less, 4)*100), size=2) +
  labs(title=paste0("homo_ref_site_",i), x="Downstream", y="Upstream") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_text(size=15),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12))


KIRC_trinuc_stage234_heatmap <- 
  ggplot(data=KIRC_trinuc_stage_heatmap_data, aes(x=Downstream, Upstream)) +
  geom_tile(aes(fill = trinuc_Size_3.1_and_up*100), color="white") + 
  scale_fill_gradient(low="white", high="dark green", name="Percent") +
  facet_grid(.~mutation)+
  geom_text(aes(label = round(trinuc_Size_3.1_and_up, 4)*100), size=2) +
  labs(title=paste0("any_alt_site_",i), x="Downstream", y="Upstream") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.ticks = element_blank(),
                     strip.text = element_text(size=15),
                     axis.title.x = element_text(size=15),
                     axis.title.y = element_text(size=15),
                     axis.text.x = element_text(size=12),
                     axis.text.y = element_text(size=12))

KIRC_trinuc_stage1_heatmap
ggsave(paste0("germline_site_",i,"_trinuc_heatmap_homoref.png"), width=8, height=2.25)

KIRC_trinuc_stage234_heatmap
ggsave(paste0("germline_site_",i,"_trinuc_heatmap_anyalt.png"), width=8, height=2.25)

combined_trinuc_heatmap_stage <- 
  plot_grid(KIRC_trinuc_stage1_heatmap, KIRC_trinuc_stage234_heatmap, 
            labels = c("A", "B"), align="h", ncol=1)
combined_trinuc_heatmap_stage
ggsave(paste0("germline_site_",i,"_trinuc_heatmap_total.png"), width=8, height=4.5)


#### gene summary table
kirc_02@selection_results$gene=apply(data.frame(kirc_02@selection_results[,1]),1,FUN = function(x){strsplit(x,"_")[[1]][1]})
summary_multivar <- function(data, group_name) {
  data_clean <- data %>% 
    filter(group == group_name) %>%
    arrange(desc(selection_intensity)) %>%
    filter(selection_intensity > 1)
  
  # Summarise information of gene with multiple variants
  info1 <- data_clean %>% group_by(gene) %>%
    summarise(cum_si = sum(selection_intensity), # change sum to mean and sd
              mean_si = mean(selection_intensity),
              sd = sd(selection_intensity),
              max_si = max(selection_intensity),
              n_variant = n_distinct(variant)) %>%
    filter(n_variant > 1) 
  
  top_variant <- data_clean %>%
    group_by(gene) %>% filter(row_number() == 1)
  
  gene_variants <- data_clean %>%
    group_by(gene) %>%
    summarise(all_variants=paste(unique(variant),collapse = ","))
  
  merge_info <- merge(info1, top_variant[, -3], by.x = "gene") %>%
    arrange(desc(cum_si), desc(n_variant))
  
  merge_info <- merge(merge_info, gene_variants, by.x = "gene") %>%
    arrange(desc(cum_si), desc(n_variant))
  
  return(merge_info)
}

stage_data <- 
  data.frame(variant = kirc_02@selection_results$variant,
             gene = kirc_02@selection_results$gene,
             selection_intensity = kirc_02@selection_results$selection_intensity,
             group = kirc_02@selection_results$progression,stringsAsFactors = F)

stage1_info <- summary_multivar(stage_data, group_name = "homo_ref")
stage234_info <- summary_multivar(stage_data, group_name = "any_alt")
allsnvs=data.frame(analysis$mutations$snv)
#gene="TP53"
#variants_ref_gene=stage_data[stage_data$group=="homo_ref" & stage_data$variant%in%strsplit(stage1_info[stage1_info$gene==gene,]$all_variants,",")[[1]],]
#allsnvs=allsnvs[allsnvs$assoc_aa_mut %in% variants_ref_gene$variant,]
#names(allsnvs)[6] = "variant"
#variants_ref_gene=merge(variants_ref_gene,allsnvs,by="variant")

fill_na <- function(x, fill = 0) {
  x = ifelse(is.na(x), fill, x)
  return(x)
}


sd_info <- stage_data %>% 
  filter(selection_intensity > 1) %>%
  group_by(gene) %>% summarise(sd = sd(selection_intensity))

stage_merge <- merge(stage1_info, stage234_info, by = "gene", all = T, 
                     suffixes = c(".1", ".234")) %>%
  mutate_at(c("cum_si.1", "cum_si.234", 
              "mean_si.1", "mean_si.234", 
              "sd.1", "sd.234",
              "n_variant.1", "n_variant.234"), fill_na) %>%
  mutate(n_variant_total = n_variant.1 + n_variant.234, 
         mean_si_total = (cum_si.1 + cum_si.234) / n_variant_total) %>%
  arrange(desc(n_variant_total))

stage_merge <- stage_merge %>%
  mutate(n_variant = paste(as.character(n_variant.1), 
                           as.character(n_variant.234), sep = "|"),
         topvar.1 = ifelse(is.na(variant.1), "NA", 
                           paste0(variant.1, "(", 
                                  round(max_si.1, 1),"x10^0)")),
         topvar.234 = ifelse(is.na(variant.234), "NA", 
                             paste0(variant.234, "(", 
                                    round(max_si.234, 1),"x10^0)")),
         topvar = paste(topvar.1, topvar.234, sep = "|"),
         mean_si = paste(round(mean_si.1, 2), round(mean_si.234, 2), sep = "|"),
         mean_si_total = round(mean_si_total, 2),
         sd = paste(round(sd.1, 2), round(sd.234, 2), sep = "|")
  )
stage_merge_sub <- stage_merge %>%
  select("gene", "n_variant", "mean_si", "sd", "topvar") %>%
  # arrange(desc(cum_si_total)) %>%
  mutate(format = "homo_ref | any_alt")

names(stage_merge_sub) <- 
  c("gene", "#variant", "mean SI", "sd", "TopVariant", "format")

write.csv(stage_merge_sub, file = paste0("germline_site_",i,"_gene_summary.csv"), quote = F, row.names = F)

library(dplyr)
library(ggplot2)
library(RColorBrewer)

# boxplot of SI

stage_data <- 
  data.frame(variant = kirc_02@selection_results$variant,
             gene = kirc_02@selection_results$gene,
             selection_intensity = kirc_02@selection_results$selection_intensity,
             group = kirc_02@selection_results$progression) %>%
  filter(selection_intensity > 1)

SR=data.frame(kirc_02@selection_results)
stage1_info <- summary_multivar(stage_data, group_name = "homo_ref")
stage234_info <- summary_multivar(stage_data, group_name = "any_alt")

genes_Pri=as.character(stage1_info[order(stage1_info$mean_si,decreasing = T),]$gene[1:10])
genes_Meta=as.character(stage234_info[order(stage234_info$mean_si,decreasing = T),]$gene[1:10])

top_gene_stage=genes_Pri
top_gene_stage=genes_Meta

stage_plot_data <- stage_data %>% 
  filter(gene %in% top_gene_stage) %>%
  mutate(gene = factor(gene, levels = top_gene_stage))

si_boxplot <- function(data, group_names, genes, colormap, yticks,main) {
  palette <- brewer.pal(6, colormap)[c(2, 5)]
  myplt <-
    boxplot(selection_intensity ~ group*gene, data = data, boxwex=0.4,
            col = palette, xlab = "", ylab = "", 
            xaxt="n", yaxt="n")
  title(ylab = expression(paste("Selection Intensity /", "10"^"4")), 
        mgp = c(2, 0, 0))
  title(xlab = "Gene", mgp = c(2, 0, 0))
  title(main=main)
  
  axis(1, mgp = c(0, 0.2, 0), 
       at = seq(1.5 , 20 , 2), 
       labels = genes, 
       tick=FALSE , cex.axis=0.5)
  
  axis(2, at = yticks * 1e4, las=2,
       labels = yticks, cex.axis=0.6)
  
  # Add the grey vertical lines
  for(i in seq(0.5 , 21 , 2)){ 
    abline(v=i, lty=1, col="grey")
  }
  
  # Add a legend
  legend("topright", legend = group_names, 
         col=palette, 
         pch = 15, bty = "n", pt.cex = 2, cex = 1,  horiz = F)
}

pdf(paste0("germline_site",i,"_topref_genes.pdf"), width = 8, height = 6)
si_boxplot(stage_plot_data, c("any_alt","homo_ref"), 
           top_gene_stage, "Greens", yticks = seq(0, 40, 2.5),main="homo_ref_genes")
dev.off()



pdf(paste0("germline_site",i,"_topalt_genes.pdf"), width = 8, height = 6)
si_boxplot(stage_plot_data, c( "any_alt","homo_ref"), 
           top_gene_stage, "Greens", yticks = seq(0, 40, 2.5),main="any_alt_genes")
dev.off()











