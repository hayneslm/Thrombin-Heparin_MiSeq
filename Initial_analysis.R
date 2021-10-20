install.packages("tidyverse")
install.packages("dplyr")
library(tidyverse)
library(readr)
library(scales)
library(dplyr)
library(ggplot2)


#import datasets
Input_Thrombin = read_tsv('Total_Input_Throm.txt') %>% 
  select(!"...1") %>% 
  filter(!grepl('-', Variant), !grepl('S0', Variant), !grepl('B', Variant)) %>%
  na.omit() %>% rename(type = Phenotype) %>% 
  mutate(score = if_else(padj<0.05, "pass", "fail"), 
         log2FoldChange = -log2FoldChange)

Input_Heparin = read_tsv('Total_Input_Heparin.txt') %>% 
  select(!"...1") %>% 
  filter(!grepl('-', Variant), !grepl('S0', Variant), !grepl('B', Variant)) %>%
  na.omit() %>% rename(type = Phenotype) %>% 
  mutate(score = if_else(padj<0.05, "pass", "fail"), 
         log2FoldChange = -log2FoldChange)

Input_uPA = read_tsv('Total_Input_uPA.txt') %>% 
  select(!"...1") %>% 
  filter(!grepl('-', Variant), !grepl('S0', Variant), !grepl('B', Variant)) %>%
  na.omit() %>% rename(type = Phenotype) %>% 
  mutate(score = if_else(padj<0.05, "pass", "fail"), 
         log2FoldChange = -log2FoldChange)

#Make MA Plots
ggplot()+
  scale_x_log10()+
  geom_hline(yintercept = 0, color = "black") +
  ggtitle('Thrombin MA Plot')+
  theme_classic()+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_Thrombin, type == 'Missense'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_Thrombin, type == 'Nonsense'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_Thrombin, type == 'wt'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_vline(xintercept = 10, color = "black", linetype = "dashed")

ggplot()+
  scale_x_log10()+
  geom_hline(yintercept = 0, color = "black") +
  ggtitle('Thrombin+Heparin MA Plot')+
  theme_classic()+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_Heparin, type == 'Missense'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_Heparin, type == 'Nonsense'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_Heparin, type == 'wt'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_vline(xintercept = 10, color = "black", linetype = "dashed")

ggplot()+
  scale_x_log10()+
  geom_hline(yintercept = 0, color = "black") +
  ggtitle('uPA MA Plot')+
  theme_classic()+
  scale_shape_manual(values = c(1,19))+
  geom_jitter(data = subset(Input_uPA, type == 'Missense'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_uPA, type == 'Nonsense'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_jitter(data = subset(Input_uPA, type == 'wt'),
              aes(x = baseMean, y = log2FoldChange, 
                  color = type, shape = score),
              stroke = 0.5, size = 1.2)+
  geom_vline(xintercept = 10, color = "black", linetype = "dashed")

#Filter data, basemean = 10, padj < = 0.05
Input_Thrombin_filtered= Input_Thrombin %>% filter(baseMean > 100, padj <= 0.05)
Input_uPA_filtered = Input_uPA %>% filter(baseMean > 100, padj <= 0.05)
Input_Heparin_filtered = Input_Heparin %>% filter(baseMean > 100, padj <= 0.05)

write_tsv(Input_Thrombin_filtered, "Input_thormbin_filtered.xls")
write_tsv(Input_uPA_filtered, "Input_uPA_filtered.xls")
write_tsv(Input_Heparin_filtered, "Input_heparin_filtered.xls")