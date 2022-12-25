library(ggpubr)
library(rbiom)
library(tidyverse)
library(phytools)
library(ggtree)
library(scatterpie)
library(ggnewscale)
library(ggthemes)

### Sample metadata
srr_dat <- read.delim("./data/sample_metadata/srr_with_host_taxonomy.txt",
                      stringsAsFactors = F)
srr_dat %>% # fix order for conifers
  mutate(Order = case_when(Host_genus %in% c("Juniperus", "Cupressus", "Thuja") ~ "Cupressales",
                           Host_genus %in% c("Pinus", "Tsuga", "Abies", "Picea") ~ "Pinales",
                           TRUE ~ Order)) -> srr_dat

import_reads_tax<- function(target, region, layout, split=NA, np=FALSE, cp=FALSE, rank=NA){
  if (is.na(split)){
    biom_path <- paste0("./data/otu_tables/feature-table_", target, "_", region, "_", layout, ".biom")
    tax_path <- paste0("./data/constax/constax_taxonomy_", target, "_", region, "_", layout, ".txt")
    fname_path <- paste0("./data/sample_metadata/fnames_", target, "_", region, "_", layout, ".txt")
  }
  else {
    biom_path <- paste0("./data/otu_tables/feature-table_", target, "_", region, "_", layout, "_", split, ".biom")
    tax_path <- paste0("./data/constax/constax_taxonomy_", target, "_", region, "_", layout, "_", split, ".txt")
    fname_path <- paste0("./data/sample_metadata/fnames_", target, "_", region, "_", layout, "_", split, ".txt")
    
  }
  print(biom_path)
  print(tax_path)
  print(fname_path)
  
  otu_dat <- read.biom(biom_path)$counts %>%
    as.matrix(.) %>%
    as.data.frame(.) %>%
    rownames_to_column(var = "OTU_ID") %>%
    mutate(OTU_name = paste("ASV", row_number(), sep="_"))
  if (np){
    otu_dat %>%
      dplyr::select(starts_with("JL_")) %>%
      colSums(.) -> no_passing_reads
    no_passing_reads[no_passing_reads == 0] %>%
      names(.) -> np
    return(np)
  }
  else {
    tax_dat <- read.delim(tax_path, stringsAsFactors = F)
    comb_dat <- left_join(otu_dat, tax_dat, by = "OTU_ID")
    fnames_dat <- read.delim(fname_path, skip = 1, header = F, stringsAsFactors = F) %>%
      mutate(sample_name = str_extract(V1, "JL_\\d*"),
             SRR_acc = str_extract(V1, ".RR\\d*")) %>%
      left_join(srr_dat, by = "SRR_acc") %>%
      rename(Host_order=Order)

    long_dat <- comb_dat %>%
      filter(! High_level_taxonomy %in% c("Mitochondria", "Chloroplast"),
             ! str_detect(High_level_taxonomy,
                          paste(c("Alveolata", "Archaeplastida", "Archaea"),
                                collapse = "|"))) %>%
      pivot_longer(cols = starts_with("JL_"),
                   names_to = "sample_name", values_to = "Reads") %>%
      left_join(fnames_dat %>%
                  dplyr::select(sample_name, SRR_acc, Compartment, Host_genus, Host_order),
                by = "sample_name")
    if (! cp && is.na(rank)){
      return(long_dat)
    }
    else if(cp){
      long_dat %>%
        group_by(Compartment, Host_order) %>%
        summarise(total_reads_com_host = sum(Reads)) -> out_dat
      return(out_dat)
    }
    else {
      long_dat %>%
        group_by_at(vars(c("Compartment", "Host_order", rank))) %>%
        summarise(total_reads = sum(Reads)) -> out_dat
      return(out_dat)
    }
  }
}



### No reads
# np_its1.its2 <- import_reads_tax("Fungi", "ITS1-ITS2", "pe", np = T)
# np_its1 <- import_reads_tax("Fungi", "ITS1", "pe", np = T)
# np_its2 <- import_reads_tax("Fungi", "ITS2", "pe", np = T)
# 
# no_reads_df <- data.frame(Target = "Fungi",
#                           Region = c(rep("ITS1-ITS2", length(np_its1.its2)),
#                                      rep("ITS1", length(np_its1)),
#                                      rep("ITS2", length(np_its2))),
#                           Sample = c(np_its1.its2, np_its1, np_its2))
# no_reads_df
# write_delim(no_reads_df, "./data/sample_metadata/no_reads_fungi.txt")# Single end reads for those which had no passing reads from PE



### Combine data

fungi_its1.its2_long_pe <- import_reads_tax("Fungi", "ITS1-ITS2", "pe")
fungi_its1_long_pe <- import_reads_tax("Fungi", "ITS1", "pe")
fungi_its2_long_pe <- import_reads_tax("Fungi", "ITS2", "pe")
fungi_its1.its2_long_se <- import_reads_tax("Fungi", "ITS1-ITS2", "se")
fungi_its1_long_se <- import_reads_tax("Fungi", "ITS1", "se")
fungi_its2_long_se <- import_reads_tax("Fungi", "ITS2", "se")

fungi_its1.its2_long_se %>%
  filter(is.na(Host_order)) %>%
  dplyr::select(SRR_acc) -> a
a
unique(a$SRR_acc)

all_fungi_long <- rbind(fungi_its1.its2_long_pe,
                        fungi_its1_long_pe,
                        fungi_its2_long_pe,
                        fungi_its1.its2_long_se,
                        fungi_its1_long_se,
                        fungi_its2_long_se) %>%
  filter(Phylum != "")

rm(fungi_its1.its2_long_pe,
   fungi_its1_long_pe,
   fungi_its2_long_pe,
   fungi_its1.its2_long_se,
   fungi_its1_long_se,
   fungi_its2_long_se)

unique(all_fungi_long$SRR_acc) -> srr_acc_fungi
srr_dat %>%
  filter(SRR_acc %in% srr_acc_fungi) -> srr_dat_fungi
length(unique(srr_dat_fungi$BioProject))

all_fungi_long %>%
  .$High_level_taxonomy -> HL_others
unique(HL_others)

all_fungi_long %>%
  group_by(Compartment, Host_order) %>%
  summarise(total_reads_com_host = sum(Reads)) -> all_fungi_com_host_totals
all_fungi_com_host_totals
write.table(all_fungi_com_host_totals, "./data/intermediate/all_fung_com_host_totals.txt", sep = "\t")
all_fungi_com_host_totals <- read.table("./data/intermediate/all_fungi_com_host_totals.txt", sep = "\t")

all_fungi_long %>%
  group_by(Compartment, Host_order, Class) %>%
  summarise(total_reads = sum(Reads)) -> all_fungi_class_totals
all_fungi_class_totals
write.table(all_fungi_class_totals, "./data/intermediate/all_fungi_class_totals.txt", sep = "\t")
all_fungi_class_totals <- read.table("./data/intermediate/all_fun_class_totals.txt", sep = "\t")


left_join(all_fungi_class_totals,
          all_fungi_com_host_totals,
          by = c("Compartment", "Host_order")) %>%
  mutate(prop = total_reads/total_reads_com_host,
         Class = str_remove(Class, "_.*")) -> all_fungi_read_props
all_fungi_read_props
write.table(all_fungi_read_props, "./data/intermediate/all_fungi_read_props.txt", sep = "\t")
all_fungi_read_props <- read.table("./data/intermediate/all_fungi_read_props.txt", sep = "\t")


all_fungi_read_props %>%
  group_by(Class) %>%
  summarize(prop_reads_total = sum(prop, na.rm = T)) %>%
  arrange(desc(prop_reads_total)) -> all_fungi_class_props_ordered

if (sum(all_fungi_class_props_ordered[1:12,"Class"] == "") > 0){
  top_tax_fungi <- all_fungi_class_props_ordered[1:13,] %>%
    filter(Class != "") %>%
    pull(Class)
  others_fungi <- all_fungi_class_props_ordered %>%
    filter(! Class %in% top_tax_fungi) %>%
    pull(Class)
} else {
  top_tax_fungi <- all_fungi_class_props_ordered[1:12, "Class"]
  others_fungi <- all_fungi_class_props_ordered %>%
    filter(! Class %in% top_tax_fungi) %>%
    pull(Class)
}
others_fungi

# all_fungi_read_props %>%
#   mutate(Class = case_when(str_detect(Class, paste(top_tax_fungi, collapse = "|")) ~ Class,
#                            TRUE ~ "Other")) %>%
#   # filter(Compartment == "Endophytic") %>%
#   filter(! is.na(Host_order),
#          ! is.na(Compartment)) %>%
#   mutate(ypos = cumsum(prop) - 0.5*prop) %>%
#   ggplot(aes(x = "", y = prop, fill = Class))+
#   geom_bar(stat="identity") + 
#   theme_void() +
#   coord_polar("y", start=0) +
#   # scale_fill_discrete(labels = c(top_tax[order(top_tax)], "Other")) +
#   facet_wrap(vars(Host_order, Compartment)) -> g
# g  
# ggsave("./figures/fungi_piecharts.png", g, height = 6, width = 4)
# ggsave("./figures/fungi_piecharts.pdf", g, height = 6, width = 10)


all_fungi_long %>%
  filter(! is.na(Host_order),
         ! is.na(Compartment)) %>%
  arrange(desc(Reads)) %>%
  dplyr::select(Phylum, Class, Reads, SRR_acc) -> a
a
unique(a$SRR_acc)
### Bacteria
bact_v1.v3_long_se <- import_reads_tax("Bacteria", "v1-v3", "se")
bact_v3.v4_long_se <- import_reads_tax("Bacteria", "v3-v4", "se")
bact_v4.1_long_se <- import_reads_tax("Bacteria", "v4", "se", split = 1)
bact_v4.2_long_se <- import_reads_tax("Bacteria", "v4", "se", split = 2)
bact_v4.3_long_se <- import_reads_tax("Bacteria", "v4", "se", split = 3)
bact_v4.v5_long_se <- import_reads_tax("Bacteria", "v4-v5", "se")
bact_v5.v6_1_long_se <- import_reads_tax("Bacteria", "v5-v6", "se", split = 1)
bact_v5.v6_2_long_se <- import_reads_tax("Bacteria", "v5-v6", "se", split = 2)
bact_v5.v6_3_long_se <- import_reads_tax("Bacteria", "v5-v6", "se", split = 3)
bact_v5.v6_4_long_se <- import_reads_tax("Bacteria", "v5-v6", "se", split = 4)
bact_v5.v6_5_long_se <- import_reads_tax("Bacteria", "v5-v6", "se", split = 5)
bact_v6.v7_long_se <- import_reads_tax("Bacteria", "v6-v7", "se")

unique(bact_v4.v5_long_se$Host_genus)

all_bact_long <- rbind(bact_v1.v3_long_se,
                       bact_v3.v4_long_se,
                       bact_v4.1_long_se,
                       bact_v4.2_long_se,
                       bact_v4.3_long_se,
                       bact_v4.v5_long_se,
                       bact_v5.v6_1_long_se,
                       bact_v5.v6_2_long_se,
                       bact_v5.v6_3_long_se,
                       bact_v5.v6_4_long_se,
                       bact_v5.v6_5_long_se,
                       bact_v6.v7_long_se) %>%
  filter(Rank_2 != "")

rm(bact_v1.v3_long_se,
   bact_v3.v4_long_se,
   bact_v4.1_long_se,
   bact_v4.2_long_se,
   bact_v4.3_long_se,
   bact_v4.v5_long_se,
   bact_v5.v6_1_long_se,
   bact_v5.v6_2_long_se,
   bact_v5.v6_3_long_se,
   bact_v5.v6_4_long_se,
   bact_v5.v6_5_long_se,
   bact_v6.v7_long_se)

write.table(all_bact_long, "./data/intermediate/all_bacteria_long.txt", sep = "\t")
all_bact_long <- read.table("./data/intermediate/all_bacteria_long.txt", sep = "\t")

all_bact_long %>%
  .$High_level_taxonomy -> HL_others
unique(HL_others)

all_bact_long %>%
  group_by(Compartment, Host_order) %>%
  summarise(total_reads_com_host = sum(Reads)) -> all_bact_com_host_totals
all_bact_com_host_totals
write.table(all_bact_com_host_totals, "./data/intermediate/all_bact_com_host_totals.txt", sep = "\t")
all_bact_com_host_totals <- read.table("./data/intermediate/all_bact_com_host_totals.txt", sep = "\t")

all_bact_long %>%
  group_by(Compartment, Host_order, Rank_2) %>%
  summarise(total_reads = sum(Reads)) -> all_bact_r2_totals
all_bact_r2_totals %>%
  mutate(Rank_2 = str_to_sentence(Rank_2)) %>%
  group_by(Compartment, Host_order, Rank_2) %>%
  summarise(total_reads = sum(total_reads))  -> all_bact_r2_totals

write.table(all_bact_r2_totals, "./data/intermediate/all_bact_r2_totals.txt", sep = "\t")
all_bact_r2_totals <- read.table("./data/intermediate/all_bact_r2_totals.txt", sep = "\t")


left_join(all_bact_r2_totals,
          all_bact_com_host_totals,
          by = c("Compartment", "Host_order")) %>%
  mutate(prop = total_reads/total_reads_com_host,
         Rank_2 = str_remove(Rank_2, "_.*")) -> all_bact_read_props
view(all_bact_read_props)

write.table(all_bact_read_props, "./data/intermediate/all_bact_read_props.txt", sep = "\t")
all_bact_read_props <- read.table("./data/intermediate/all_bact_read_props.txt", sep = "\t")

all_bact_read_props %>%
  group_by(Rank_2) %>%
  summarize(prop_reads_total = sum(prop, na.rm = T)) %>%
  arrange(desc(prop_reads_total)) -> all_bact_class_props_ordered
view(all_bact_class_props_ordered)
n <- 20
if (sum(all_bact_class_props_ordered[1:n,"Rank_2"] == "") > 0){
  top_tax_bact <- all_bact_class_props_ordered[1:n+1,] %>%
    filter(Rank_2 != "") %>%
    pull(Rank_2)
  others_bact <- all_bact_class_props_ordered %>%
    filter(! Rank_2 %in% top_tax_bact) %>%
    pull(Rank_2)
} else {
  top_tax_bact <- all_bact_class_props_ordered[1:n,] %>%
    pull(Rank_2)
  others_bact <- all_bact_class_props_ordered %>%
    filter(! Rank_2 %in% top_tax_bact) %>%
    pull(Rank_2)
}
others_bact

# bact_v1.v3_long_se_cp <- import_reads_tax("Bacteria", "v1-v3", "se", cp=T)
# bact_v3.v4_long_se_cp <- import_reads_tax("Bacteria", "v3-v4", "se", cp=T)
# bact_v4.1_long_se_cp <- import_reads_tax("Bacteria", "v4", "se", split = 1, cp=T)
# bact_v4.2_long_se_cp <- import_reads_tax("Bacteria", "v4", "se", split = 2, cp=T)
# bact_v4.3_long_se_cp <- import_reads_tax("Bacteria", "v4", "se", split = 3, cp=T)
# bact_v4.v5_long_se_cp <- import_reads_tax("Bacteria", "v4-v5", "se", cp=T)
# bact_v5.v6_1_long_se_cp <- import_reads_tax("Bacteria", "v5-v6", "se", split = 1, cp=T)
# bact_v5.v6_2_long_se_cp <- import_reads_tax("Bacteria", "v5-v6", "se", split = 2, cp=T)
# bact_v5.v6_3_long_se_cp <- import_reads_tax("Bacteria", "v5-v6", "se", split = 3, cp=T)
# bact_v5.v6_4_long_se_cp <- import_reads_tax("Bacteria", "v5-v6", "se", split = 4, cp=T)
# bact_v5.v6_5_long_se_cp <- import_reads_tax("Bacteria", "v5-v6", "se", split = 5, cp=T)
# bact_v6.v7_long_se_cp <- import_reads_tax("Bacteria", "v6-v7", "se", cp=T)


# all_bact_com_host_totals <- rbind(bact_v1.v3_long_se_cp,
#                                   bact_v3.v4_long_se_cp,
#                                   bact_v4.1_long_se_cp,
#                                   bact_v4.2_long_se_cp,
#                                   bact_v4.3_long_se_cp,
#                                   bact_v4.v5_long_se_cp,
#                                   bact_v5.v6_1_long_se_cp,
#                                   bact_v5.v6_2_long_se_cp,
#                                   bact_v5.v6_3_long_se_cp,
#                                   bact_v5.v6_4_long_se_cp,
#                                   bact_v5.v6_5_long_se_cp,
#                                   bact_v6.v7_long_se_cp)

### Import plant phylogenetic tree from https://doi.org/10.1038/s41586-019-1693-2
small_tree <- ape::read.tree("./data/phylo/pruned_tree.tre")
plot(small_tree)
mrca_node <- findMRCA(small_tree, c("SZYG", "HMHL"))
rooted_small_tree <- reroot(small_tree, mrca_node)
rooted_small_tree
rooted_small_tree <- drop.tip(rooted_small_tree, "EWXK")
plot.new()
plot(rooted_small_tree)
nodelabels(frame = NULL, cex = 0.4)
tiplabels(frame = NULL, cex = 0.4)
rooted_small_tree

moss_node <- findMRCA(rooted_small_tree, c("DHWX", "GOWD"))
fern_node <- 33
liv_node <- 35
lyco_node <- 34
gymno_node <- findMRCA(rooted_small_tree, c("XMGP", "DZQM"))
mono_node <- findMRCA(rooted_small_tree, c("BYQM", "BPKH"))
mag.chlo_node <- findMRCA(rooted_small_tree, c("OSHQ", "CSSK"))
eudi_node <- findMRCA(rooted_small_tree, c("AUGV", "ACFP"))
# library(ggthemes)


ggtree(rooted_small_tree) +
  geom_text(aes(label=node), hjust=-.3)
# write.tree(rooted_small_tree, "./data/phylo/pruned_tree_w_annot_rooted.tre")
ggplot(rooted_small_tree) + 
  geom_tree() +
  theme_tree() +
  geom_tiplab(align=TRUE, linesize=.5) -> g_tip_labs
g_tip_labs
order_tip_names <- read.delim("./data/phylo/order_annot_tip_label.txt", sep = "\t")
colnames(order_tip_names) <- c("Host_order", "Tip_label")

# # g <- g_tip_labs
# tip_df <- data.frame(x = g$data$x[g$data$isTip], 
#                             y = g$data$y[g$data$isTip],
#                             Tip_label = g$data$label[g$data$isTip]) %>%
#   left_join(order_tip_names, by="Tip_label")
# tip_df

all_fungi_read_props %>%
  mutate(Class = case_when(str_detect(Class, paste(top_tax_fungi, collapse = "|")) ~ Class,
                           TRUE ~ "Other")) %>%
  filter(! is.na(Host_order),
         ! is.na(Compartment),
         total_reads != 0) %>%
  group_by(Compartment, Host_order, Class) %>%
  summarise(prop = sum(prop)) -> fungi_prop_df
fungi_prop_df
# view(all_bact_read_props)
all_bact_read_props %>%
  mutate(Rank_2 = case_when(str_detect(Rank_2, paste(top_tax_bact, collapse = "|")) ~ Rank_2,
                           TRUE ~ "Other")) %>%
  filter(! is.na(Host_order),
         ! is.na(Compartment),
         total_reads != 0) %>%
  group_by(Compartment, Host_order, Rank_2) %>%
  summarise(prop = sum(prop)) -> bact_prop_df
bact_prop_df

ggplot(rooted_small_tree) +
  coord_equal() +
  geom_tree() +
  theme_tree() -> g_tree_alone
# view(order_tip_names)
g_tree_alone
g_tree_alone$data$y <- g_tree_alone$data$y*2.2
g_tree_alone$data$x <- g_tree_alone$data$x*0.3

highlight_df <- data.frame(node = c(liv_node,moss_node,lyco_node,fern_node,gymno_node,
                                    mono_node,mag.chlo_node,eudi_node),
                           group = c("Liverworts",
                                     "Mosses",
                                     "Lycophytes",
                                     "Ferns",
                                     "Gymnosperms",
                                     "Monocots",
                                     "Magnoliids and\nChloranthales",
                                     "Eudicots"),
                           color = c("#c6eac3",
                                     "#f2afdb",
                                     "#bfe9e5",
                                     "#e5c887",
                                     "#b4d584",
                                     "#d6bde3",
                                     "#bdd6f1",
                                     "#fbcbc7"))

g_tree_alone <- g_tree_alone + 
  geom_highlight(data = highlight_df,
                 aes(node = node, fill = group), alpha = 1, to.bottom = T) +
  scale_fill_manual(guide = "none",
                    breaks = highlight_df$group,
                    values = highlight_df$color)
g_tree_alone$layers[[1]]$data$xmin <- max(g_tree_alone$layers[[1]]$data$xmin)
g_tree_alone$layers[[1]]$data$xmax <- max(g_tree_alone$layers[[1]]$data$xmax)
g_tree_alone

g_tree_alone$layers[[1]]$data
highlight_df_annot <- left_join(highlight_df, g_tree_alone$data, by="node")
highlight_df_annot$x <- c(rep(-1, 5), 0.5, 2, 2)
highlight_df_annot

tip_df <- data.frame(x = g_tree_alone$data$x, 
                     y = g_tree_alone$data$y,
                     Tip_label = g_tree_alone$data$label)[g_tree_alone$data$isTip,] %>%
  left_join(order_tip_names, by="Tip_label")
tip_df

fungi_prop_df %>%
  pivot_wider(id_cols = 1:2, names_from = Class, values_from = prop) %>%
  mutate(across(c(top_tax_fungi, "Other"), ~ replace_na(.x, replace=0))) %>%
  left_join(tip_df, by = "Host_order") %>%
  mutate(Target = "Fungi",
         x_offset = case_when(Compartment == "Total_phyllosphere" ~ 0,
                              Compartment == "Endophytic" ~ 1,
                              Compartment == "Epiphytic" ~ 2)) %>%
  mutate(x_end = max(g_tree_alone$data$x) + 2 + x_offset*2,
         radius = 0.9) %>%
  filter(! is.na(x)) -> fungi_pie_df
fungi_pie_df
bact_prop_df %>%
  pivot_wider(id_cols = 1:2, names_from = Rank_2, values_from = prop) %>%
  mutate(across(c(top_tax_bact, "Other"), ~ replace_na(.x, replace=0))) %>%
  left_join(tip_df, by = "Host_order") %>%
  mutate(Target = "Bacteria",
         x_offset = case_when(Compartment == "Total_phyllosphere" ~ 4,
                              Compartment == "Endophytic" ~ 5,
                              Compartment == "Epiphytic" ~ 6)) %>%
  mutate(x_end = max(g_tree_alone$data$x) + 2 + x_offset*2,
         radius = 0.9) %>%
  filter(! is.na(x)) -> bact_pie_df
bact_pie_df

annot_data <- data.frame(x_start = g_tree_alone$data$x,
                       y_start = g_tree_alone$data$y,
                       x_end = max(g_tree_alone$data$x),
                       y_end = g_tree_alone$data$y,
                       Tip_label = g_tree_alone$data$label)[g_tree_alone$data$isTip,] %>%
  left_join(tip_df, by = "Tip_label")
annot_data
label_data <- data.frame(x = c(0,1,2,4,5,6)*2 + (max(g_tree_alone$data$x) + 2),
                         y = 0,
                         label = c("Total Phyllosphere", "Endophytic", "Epiphytic"))
label_data <- rbind(label_data, data.frame(x = c(1,5)*2 + (max(g_tree_alone$data$x) + 2),
                             y = max(g_tree_alone$data$y) + 2,
                             label = c("Fungi", "Bacteria")))

label_data

g_tree_alone + 
  new_scale_fill() + 
  geom_scatterpie(aes(x = x_end, y = y, r = radius), data = fungi_pie_df,
                  cols = c(top_tax_fungi, "Other"), color = NA) +
  labs(fill = "Fungal Class") +
  scale_fill_discrete(guide = guide_legend(order = 1)) +
  new_scale_fill() + 
  geom_scatterpie(aes(x = x_end, y = y, r = radius), data = bact_pie_df,
                  cols = c(top_tax_bact), color = NA) +
  scale_fill_discrete(guide = guide_legend(order = 2)) +
  geom_segment(data = annot_data,
               aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               lty = 3) +
  geom_segment(aes(x = max(fungi_pie_df$x_end) + 2,
                   xend = max(fungi_pie_df$x_end) + 2,
                   y = min(g_tree_alone$data$y),
                   yend = max(g_tree_alone$data$y)+2)) +
  geom_text(data = annot_data,
            aes(x = x_end+16, y = y_end, label = Host_order),
            hjust = 0, vjust = 0.5, angle = 0, size = 3) +
  geom_text(data = label_data[1:6,],
            aes(x = x, y = y, label = label),
            hjust = 0, vjust = 0.5, angle = -60, size = 3) +
  geom_text(data = label_data[7:8,],
            aes(x = x, y = y, label = label),
            hjust = 0.5, vjust = 0, angle = 0, size = 3) +
  geom_text(data = highlight_df_annot,
            aes(x = x, y = y, label = group, color = group),
            hjust = 1, vjust = 0.5, size = 3)+
  scale_color_manual(guide = "none",
                    breaks = highlight_df$group,
                    values = highlight_df$color) +
  theme(legend.position = c(1.17, 0.5),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(9, "pt"),
        plot.margin = margin(2, 2, 2, 2, "pt")) +
  labs(fill = "Bacterial Phylum") +
  guides(fill=guide_legend(ncol=1)) +
  scale_x_continuous(limits = c(-15, max(annot_data$x)+30)) + 
  scale_y_continuous(limits = c(-10, max(annot_data$y)+2)) -> g
g
ggsave("./figures/tree_with_pies.png", g, height = 8, width = 8)
ggsave("./figures/tree_with_pies.svg", g, height = 8, width = 8)
