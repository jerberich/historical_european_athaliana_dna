
## Setup and run Grenepipe

| Name | Source |
| ----------- | ----------- |
| Grenepipe (v0.13.0) | [https://github.com/moiexpositoalonsolab/grenepipe/wiki](https://github.com/moiexpositoalonsolab/grenepipe/wiki) |
| config.yaml | [code/config_historical.yaml](code/config_historical.yaml) |
| Paragraph | Text |

## Install Grenepipe
```bash
git clone https://github.com/moiexpositoalonsolab/grenepipe.git
cd grenepipe
mamba env create -f envs/grenepipe.yaml
```

Use Grenepipe tool to create sample table for each sequencing run
```bash
./tools/generate-table.py  path/to/sequencing/run1/  path/to/sequencing/run1/run1_sample_table.tsv
./tools/generate-table.py  path/to/sequencing/run2/  path/to/sequencing/run2/run2_sample_table.tsv
# ...
```

Merge sample tables across sequencing runs
```R
library(tidyverse)

df_1 <- read_tsv('path/to/sequencing/run1/run1_sample_table.tsv')
df_2 <- read_tsv('path/to/sequencing/run2/run2_sample_table.tsv')
# ...

# Set one of the datasets as unit '2' (add more as needed)
df_2 <- df_2 %>% 
  mutate(unit = 2)

# Merge sample tables and export
df <- rbind(df_1, df_2) %>% 
  arrange(sample, unit)

df %>% write_tsv('path/to/grenepipe/output/sample_table.tsv')
```

Move config file and run grenepipe
```bash
cp config_historical.yaml /path/to/grenepipe/output/config.yaml

snakemake \
    --conda-frontend mamba \
    --conda-prefix /path/to/grenepipe/output/conda-envs \
    --executor slurm \
    --profile /path/to/grenepipe/workflow/profiles/slurm/ \
    --directory /gpath/to/grenepipe/output/
```