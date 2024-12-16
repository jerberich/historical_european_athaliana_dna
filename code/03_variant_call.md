# Historical European *A. thaliana* DNA processing
## Variant calling

### Filter samples based on sequencing statistics
```R
library(tidyverse)


# Load QC statistics
stats <- read_tsv('/path/to/grenepipe/output/stats/erberich_summary_stats.tsv')

# Load sample table
samples <- read_tsv('/path/to/grenepipe/output/sample_table.tsv')


# Merge
df <- stats %>% 
  select(sample, qc) %>% 
  mutate(sample = str_replace(sample, ' ', '_')) %>%
  right_join(samples) %>% 
  mutate(qc = replace_na(qc, 'keep'))

df_clean <- df %>% 
  filter(qc == 'keep') %>% 
  select(sample, unit, platform, fq1, fq2)

df_clean %>% write_tsv('/path/to/grenepipe/output/sample_table.tsv')
```

### Run Grenepipe
```bash
snakemake \
    --conda-frontend mamba \
    --conda-prefix /path/to/grenepipe/output/conda-envs \
    --executor slurm \
    --profile /path/to/grenepipe/workflow/profiles/slurm/ \
    --directory /gpath/to/grenepipe/output/
```

### Variant calling statistics
Calculate with bcftools. Here is a helpful guide to what each of these statistics means: [https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/](https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/)
```bash
bcftools query genotyped/all.vcf.gz -f '%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n' > qc/all.vcf_FS.SOR.MQRS.RPRS.QD.MQ.DP.txt
```

Plot with R and ggplot2
```R

# Load qc stats
df <- read_tsv('all.vcf_FS.SOR.MQRS.RPRS.QD.MQ.DP.txt',
               col_names = c('FS','SOR','MQRS','RPRS','QD','MQ','DP'), 
               na = c('','NA','.'))

colMeans(is.na(df))
df <- df %>% 
  mutate(type = 'historical')

# Plot FS (GATK recommends removing >60)

cutoff = 0.3
lbl <- df %>% 
  summarise(removed = sum(FS > 0.09) *100 / n()) %>% 
  pull(removed) %>% 
  round(digits = 1)

df %>% 
  ggplot(aes(sqrt(FS), fill=type)) +
  geom_histogram(position = "identity",binwidth = .2,alpha = 0.5) +
  # scale_y_log10() +
  # theme_bw() +
  geom_vline(xintercept = cutoff, color ='blue', linetype="dotted")+ 
  # annotate("label", x=Inf, y=Inf, label = paste0('Cutoff: ',as.character(cutoff)), hjust = 2, vjust = 2) +
  # annotate("label", x=Inf, y=Inf, label = paste0('Removes ',as.character(lbl)[1],'% GSVD\nRemoves ',as.character(lbl)[2],'% nonGSVD'), hjust = 1.5, vjust = 3) +
  # geom_text("label", x=Inf, y=Inf, label = paste0('Cutoff: ',as.character(cutoff),'\n','Filters:\n',as.character(lbl)[1],'% GSVD\n',as.character(lbl)[2],'% nonGSVD')) +
  scale_fill_discrete(labels=c(paste0(as.character(lbl)[1],'% historical')), type = c('black')) +
  scale_y_continuous(label = unit_format(unit = "M", scale = 1e-6, sep = "")) +
  labs(y = 'SNPs', fill = paste0('Cutoff > ',as.character(cutoff),'\n','filters:'), title = 'Fisher Stand') +
  theme_classic()+
  theme(aspect.ratio=1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        legend.position = c(0.80, 0.75),
        legend.box.background = element_rect(color = "grey20", linewidth = 1),
        legend.key = element_blank()) +
  guides(fill = guide_legend(override.aes = aes(label = "", fill = NA, size = NA)))

plot_file_path <- '/path/to/plots/erberich_data_vcf_filter_'
ggsave(filename = paste0(plot_file_path, 'FS.png'),
       dpi = 300,
       width = 200,
       height = 200,
       scale = 7,
       units = "px",
)
# # Plot SOR (GATK recommends removing > 3)
cutoff = 1.7
lbl <- df %>%
  summarise(removed = sum(SOR > cutoff ) *100 / n()) %>%
  pull(removed) %>%
  round(digits = 1)

df %>%
  ggplot(aes(SOR)) +
  geom_histogram(binwidth = 0.12) ++
  theme_bw()+
  geom_vline(xintercept = cutoff, color ='blue', linetype="dotted")+
  annotate("label", x=Inf, y=Inf, label = paste0('Cutoff: ',as.character(cutoff)), hjust = 2, vjust = 2) +
  annotate("label", x=Inf, y=Inf, label = paste0('Removes ',as.character(lbl),'%'), hjust = 1.5, vjust = 3)

# Plot MQ (GATK recommends remobing < 40)
cutoff = 40
lbl <- df %>% 
  summarise(removed = sum(MQ < cutoff ) *100 / n()) %>% 
  pull(removed) %>% 
  round(digits = 1)

cutoff = 40
lbl <- df %>% 
  group_by(type) %>% 
  summarise(removed = sum(MQ < cutoff ) *100 / n()) %>% 
  pull(removed) %>% 
  round(digits = 1)

df %>% 
  ggplot(aes(MQ, fill=type)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  theme_bw() +
  geom_vline(xintercept = cutoff, color ='blue', linetype="dotted")+ 
  scale_fill_discrete(labels=c(paste0(as.character(lbl)[1],'% historical')), type = c('black')) +
  scale_y_continuous(label = unit_format(unit = "M", scale = 1e-6, sep = "")) +
  labs(y = 'SNPs', fill = paste0('Cutoff > ',as.character(cutoff),'\n','filters:'), title = 'Mapping Qualtiy') +
  theme_classic()+
  theme(aspect.ratio=1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        legend.position = c(0.25, 0.75),
        legend.box.background = element_rect(color = "grey20", linewidth = 1),
        legend.key = element_blank()) +
  guides(fill = guide_legend(override.aes = aes(label = "", fill = NA, size = NA)))

ggsave(filename = paste0(plot_file_path, 'MQ.png'),
       dpi = 300,
       width = 200,
       height = 200,
       scale = 7,
       units = "px",
)

# Plot MQRS (Note mostly missing data, GATK recommends removing < -12)
cutoff = -2
cutoff2 = 2
lbl <- df %>% 
  summarise(removed = sum(MQRS < cutoff |
                          MQRS >cutoff2, na.rm = TRUE) *100 / n()) %>% 
  pull(removed) %>% 
  round(digits = 1)

df %>% 
  ggplot(aes(MQRS, fill=type)) +
  geom_histogram(binwidth = 0.5, position = "identity", alpha = 0.5) +
  theme_bw() +
  geom_vline(xintercept = cutoff, color ='blue', linetype="dotted")+ 
  geom_vline(xintercept = cutoff2, color ='blue', linetype="dotted")+ 
  scale_fill_discrete(labels=c(paste0(as.character(lbl)[1],'% historical')), type = c('black')) +
  scale_y_continuous(label = unit_format(unit = "M", scale = 1e-6, sep = "")) +
  labs(y = 'SNPs', fill = paste0(as.character(cutoff2),' > Cutoff > ',as.character(cutoff),'\n','filters:'), title = 'MQ Rank Sum') +
  theme_classic()+
  theme(aspect.ratio=1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        legend.position = c(0.25, 0.75),
        legend.box.background = element_rect(color = "grey20", linewidth = 1),
        legend.key = element_blank()) +
  guides(fill = guide_legend(override.aes = aes(label = "", fill = NA, size = NA)))

ggsave(filename = paste0(plot_file_path, 'MQRS.png'),
       dpi = 300,
       width = 200,
       height = 200,
       scale = 7,
       units = "px",
)

# Plot QD (GATK recommends removing < 2)
cutoff = 22
lbl <- df %>% 
  summarise(removed = sum(QD < cutoff, na.rm = TRUE) *100 / n()) %>% 
  pull(removed) %>% 
  round(digits = 1)

df %>% 
  ggplot(aes(QD, fill=type)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  theme_bw()+
  geom_vline(xintercept = cutoff, color ='blue', linetype="dotted")+ 
  scale_fill_discrete(labels=c(paste0(as.character(lbl)[1],'% historical')), type = c('black')) +
  scale_y_continuous(label = unit_format(unit = "M", scale = 1e-6, sep = "")) +
  labs(y = 'SNPs', fill = paste0('Cutoff < ',as.character(cutoff),'\n','filters:'), title = 'Quality by Depth') +
  theme_classic()+
  theme(aspect.ratio=1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        legend.position = c(0.25, 0.75),
        legend.box.background = element_rect(color = "grey20", linewidth = 1),
        legend.key = element_blank()) +
  guides(fill = guide_legend(override.aes = aes(label = "", fill = NA, size = NA)))

# plot_file_path <- '/global/home/users/jerberic/plots/erberich_data_qc_'
ggsave(filename = paste0(plot_file_path, 'QD.png'),
       dpi = 300,
       width = 200,
       height = 200,
       scale = 7,
       units = "px",
)

# Plot ReadPosRankSum (Note mostly missing data, GATK recommends removing > -8.0)
cutoff = -2
cutoff2 = 2
lbl <- df %>% 
  summarise(removed = sum(RPRS < cutoff |
                            RPRS >cutoff2, na.rm = TRUE) *100 / n()) %>% 
  pull(removed) %>% 
  round(digits = 1)
  
df %>% 
  ggplot(aes(RPRS, fill=type)) +
  geom_histogram(binwidth = 0.5, position = "identity", alpha = 0.5) +
  theme_bw()+
  geom_vline(xintercept = cutoff, color ='blue', linetype="dotted")+ 
  geom_vline(xintercept = cutoff2, color ='blue', linetype="dotted")+ 
  scale_fill_discrete(labels=c(paste0(as.character(lbl)[1],'% historical')), type = c('black')) +
  scale_y_continuous(label = unit_format(unit = "M", scale = 1e-6, sep = "")) +
  labs(y = 'SNPs', fill = paste0(as.character(cutoff2),' > Cutoff > ',as.character(cutoff),'\n','filters:'), title = 'Read Position Rank Sum') +
  theme_classic()+
  theme(aspect.ratio=1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        legend.position = c(0.80, 0.75),
        legend.box.background = element_rect(color = "grey20", linewidth = 1),
        # legend.background = element_rect(fill = 'white'),
        # panel.background = element_rect(fill = "grey90"),
        # legend.key = element_rect(fill="white"),
        legend.key = element_blank()) +
  guides(fill = guide_legend(override.aes = aes(label = "", fill = NA, size = NA)))

ggsave(filename = paste0(plot_file_path, 'RPRS.png'),
       dpi = 300,
       width = 200,
       height = 200,
       scale = 7,
       units = "px",
)

# Plot cumulative read depth (guide recommends to removed those that have extremely high DP)
print(paste(mean(df$DP)," " ,median(df$DP)))

cutoff = 22
cutoff2 = 140000
lbl <- df %>% 
  summarise(removed = sum(DP < cutoff |
                          DP >cutoff2, na.rm = TRUE) *100 / n()) %>% 
  pull(removed)%>% 
  round(digits = 1)


# DP
df %>% 
  filter(DP < 70000) %>%
  ggplot(aes(DP)) +
  geom_histogram(binwidth = 100) +
  theme_bw() +
  scale_y_log10() #+
  # scale_x_log10()

# Total snps filtered
df2 <- df %>% 
  mutate(qc = case_when(FS > 0.09  | MQ < 40 | MQRS < -2 | MQRS > 2 | QD < 22 | RPRS < -2 | RPRS > 2 | DP < 22 | DP > 140000 ~ "rm", TRUE ~ "keep"))

df2 %>% 
  summarise(removed = sum(qc == 'rm')/ n())
```

### Filter for pop-gen analysis
We removed SNPs that don't have reliable calling statistics as visualized above and SNPs that are unique to a single sample with the minor allele count. Since *A. thaliana* is a selfing plant, most SNPs are homozygous. We then filter samples based on genotype missingness with PLINK, this also forces sites to be biallelic.
```bash
# Apply hard filter
bcftools filter -e 'FS>0.04 || MQ<40 || MQRankSum<-2.0 || MQRankSum>2.0 || QD<22.0 || ReadPosRankSum<-2.0 || ReadPosRankSum>2.0 || INFO/DP<22 || INFO/DP>140000' -O z -o genotyped-all.hf.vcf.gz genotyped-all.vcf.gz

# Filter to just SNPs
vcftools --gzvcf genotyped-all.hf.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip -c > genotyped-all.hf.SNPonly.vcf.gz

# Filter sites on missingness
vcftools --gzvcf genotyped-all.hf.SNPonly.vcf.gz --max-missing 0.15 --mac 3 --recode --recode-INFO-all --stdout | gzip -c > genotyped-all.hf.g15mac3.SNPonly.vcf.gz

# Filter samples on missingness and force biallelic variants
plink2 --vcf genotyped-all.hf.g15mac3.SNPonly.vcf --mind .2 --allow-extra-chr --const-fid --keep-allele-order --recode vcf --out genotyped-all.hf.g15mac3.SNPonly_bi
```

### Principal component analysis (PCA) to visualize population structure
Calculate using PLINK
```bash
plink2 --vcf genotyped-all.hf.g15mac3.SNPonly_bi.vcf --mind .2 --pca --allow-extra-chr --const-fid --keep-allele-order --out genotyped-all.hf.g15mac3.SNPonly_bi
plink2 --vcf genotyped-all.hf.SNPonly.vcf --mind .2 --pca --allow-extra-chr --const-fid --keep-allele-order --out genotyped-all.hf.SNPonly_bi
```

Visualize in R
```R
df_clean_eigen_vec <- read_tsv('genotyped-all.hf.g15mac3.SNPonly_bi.eigenvec') %>% 
  rename(sample_plant = `#IID`)
df_clean_eigen_val <- read_tsv('genotyped-all.hf.g15mac3.SNPonly_bi.eigenval', col_names = c('perc'))  %>% 
  mutate(pve = perc/sum(perc), PC = paste0('PC',row_number()))

df_eigen_vec <- read_tsv('genotyped-all.hf.SNPonly_bi.eigenvec') %>% 
  rename(sample_plant = `#IID`)
df_eigen_val <- read_tsv('genotyped-all.hf.SNPonly_bi.eigenval', col_names = c('perc')) %>% 
  mutate(pve = perc/sum(perc), PC = paste0('PC',row_number()))

# Plot pve
df_clean_eigen_val %>% 
  ggplot(aes(reorder(PC, pve), pve)) + geom_bar(stat = "identity") +
  theme_classic()+
  labs(title = 'Filtered VCF')

df_clean_eigen_vec %>% 
  left_join(samples_df) %>% 
  ggplot(aes(PC1, PC2, color = year)) +
  scale_color_viridis(direction = -1) +
  geom_point() +
  labs(title = 'PCA after filtering') +
  theme_classic() +
  theme(aspect.ratio=1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        legend.position = 'none')
plot_file_path <- '/path/to/plots/erberich_data_vcf_filter_'

ggsave(filename = paste0(plot_file_path, 'pca_post_filtering.png'),
       dpi = 300,
       width = 200,
       height = 200,
       scale = 7,
       units = "px",
)
 
# Plot pve
df_eigen_val %>% 
  ggplot(aes(reorder(PC, pve), pve)) + geom_bar(stat = "identity") +
  theme_classic()+
  labs(title = 'Hard_filtered VCF')

df_eigen_vec %>% left_join(samples_df) %>% 
  ggplot(aes(PC1, PC2, color = year)) +
  scale_color_viridis(direction = -1) +
  geom_point() +
  labs(title = 'PCA before filtering') +
  theme_classic() +
  theme(aspect.ratio=1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=14),
        legend.position = 'none')

ggsave(filename = paste0(plot_file_path, 'pca_pre_filtering.png'),
       dpi = 300,
       width = 200,
       height = 200,
       scale = 7,
       units = "px",
)                     
```

#### [Pseudoage long-read sequencing for variant call QC](code/04_pseudoage_long_read_athaliana.md)