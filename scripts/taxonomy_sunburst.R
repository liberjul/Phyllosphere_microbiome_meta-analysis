library(ggpubr)
library(tidyverse)
library(igraph)
library(ggraph)
library(svglite)

host_taxonomy_wide_matched_leaf <- read.csv("./data/biosample_data/host_table_matched_annot_leaf.csv", stringsAsFactors = F)
host_taxonomy_wide_matched_phyllo <- read.csv("./data/biosample_data/host_table_matched_annot_phyllo.csv", stringsAsFactors = F)

host_taxonomy_wide_matched <- rbind(host_taxonomy_wide_matched_leaf,
                                    host_taxonomy_wide_matched_phyllo)

host_taxonomy_wide_matched %>%
  dplyr::select(count, Name, rank_1:Genus) %>%
  filter(rank_9 != "<NA>") ->
  host_tax_filt
host_tax_filt %>%
  group_by(rank_7, rank_8) %>%
  rename(source = rank_7, target = rank_8) %>%
  summarize(sample_count = sum(count)) %>%
  mutate(rank="rank_8") -> rank_8_sum
host_tax_filt %>%
  group_by(rank_8, rank_9) %>%
  rename(source = rank_8, target = rank_9) %>%
  summarize(sample_count = sum(count)) %>%
  mutate(rank="rank_9") -> rank_9_sum
rank_9_sum
host_tax_filt %>%
  group_by(rank_9, rank_10) %>%
  rename(source = rank_9, target = rank_10) %>%
  summarize(sample_count = sum(count)) %>%
  mutate(rank="rank_10")  -> rank_10_sum
rank_10_sum
host_tax_filt %>%
  filter(rank_10 != "Lycopodiaceae")%>%
  group_by(rank_10, Order) %>%
  rename(source = rank_10, target = Order) %>%
  summarize(sample_count = sum(count)) %>%
  mutate(rank="Order")  -> order_sum
order_sum

host_tax_filt %>%
  group_by(Order, Family) %>%
  rename(source = Order, target = Family) %>%
  summarize(sample_count = sum(count)) %>%
  mutate(rank="Family")  -> ord_to_fam_count
ord_to_fam_count
host_tax_filt %>%
  group_by(Family, Genus) %>%
  rename(source = Family, target = Genus) %>%
  summarize(sample_count = sum(count)) %>%
  mutate(rank="Genus") -> fam_to_gen_count
fam_to_gen_count
comb_source_target <- rbind(rank_8_sum,
                            rank_9_sum,
                            rank_10_sum,
                            order_sum,
                            ord_to_fam_count,
                            fam_to_gen_count)

# comb_source_target %>% filter(!rank %in% c("Genus", "Family")) -> comb_source_target

comb_source_target %>%
  replace_na(list(source="Not identified", target="Not identified")) %>%
  filter(!(source == "Not identified" & target == "Not identified")) -> comb_source_target
comb_source_target %>%
  group_by(source) %>%
  summarize(source_sum = sum(sample_count)) %>%
  mutate(source_count_str = paste0(source,
                                   " (", format(source_sum, trim = T, big.mark=","), ")")) -> source_names
comb_source_target %>%
  filter(!target %in%
           comb_source_target$source) %>%
  group_by(target) %>%
  summarize(target_sum = sum(sample_count)) %>%
  mutate(target_count_str = paste0(target,
                                   " (", format(target_sum, trim = T, big.mark=","), ")")) -> target_names

source_names$source_count_str[match(comb_source_target$source,
                                    source_names$source)] ->
  comb_source_target$source_count_str
target_names$target_count_str[match(comb_source_target$target,
                                    target_names$target)] ->
  comb_source_target$target_count_str
comb_source_target$target_count_str[is.na(comb_source_target$target_count_str)] <-
  source_names$source_count_str[match(comb_source_target$target[is.na(comb_source_target$target_count_str)],
                                      source_names$source)]
comb_source_target

nodes <- data.frame(
  name=c(as.character(comb_source_target$source_count_str),
         as.character(comb_source_target$target_count_str)) %>% unique()
)
nodes$name

# With networkD3, connection must be provided using id, not using real name like in the comb_source_target dataframe.. So we need to reformat it.
comb_source_target$IDsource <- match(comb_source_target$source_count_str, nodes$name)-1
comb_source_target$IDtarget <- match(comb_source_target$target_count_str, nodes$name)-1


comb_source_target %>%
  filter(!target %in% c("", "Not identified")) %>%
  mutate(from = source,
         to = target) -> tax_edges
comb_source_target %>%
  group_by(source, source_count_str) %>%
  summarize(size = sum(sample_count)) %>%
  mutate(name = source,
         rank = "Not genus",
         annot_name = source_count_str) %>%
  dplyr::select(name, size, rank, annot_name) -> tax_vert_top

comb_source_target %>%
  filter(! target %in% c(comb_source_target$source, "Not identified", "")) %>%
  group_by(rank, target, target_count_str) %>%
  summarize(size = sum(sample_count)) %>%
  mutate(name = target,
         annot_name = target_count_str) %>%
  dplyr::select(name, size, rank, annot_name) -> tax_vert_bottom
tax_vert <- full_join(tax_vert_top, tax_vert_bottom)
tax_vert[!duplicated(tax_vert$name),] -> tax_vert
tax_vert[,2:ncol(tax_vert)] -> tax_vert

tax_vert_top[duplicated(tax_vert_top$name),]
tax_vert[duplicated(tax_vert$name),]
tax_edges[! tax_edges$to %in% tax_vert$name,]
tax_vert$Genus <- tax_vert$name
tax_vert$Genus[tax_vert$rank != "Genus"] <- ""
tax_vert$Genus[tax_vert$size < 100] <- ""
tax_vert$higherName <- tax_vert$name
tax_vert$higherName[tax_vert$rank != "Not genus"] <- ""
tax_vert$higherName[tax_vert$size < 500] <- ""
tax_vert$plot_segment <- NA
tax_vert$plot_segment[tax_vert$Genus != ""] <- 1

tax_vert[,1]

mygraph <- graph_from_data_frame(tax_edges, vertices = tax_vert)
set.seed(1234)

out_layer_lab_pos <- function(coord, r, scaling, plot = 1){
  return(coord + sign(coord)*r*scaling)
}

ggraph(mygraph, 'partition', circular = TRUE, weight = size) +
  geom_node_arc_bar(aes(fill = as.factor(depth))) +
  theme_void() +
  scale_fill_viridis(discrete = T, end = 0.8) +
  # scale_fill_manual(values = c("0" = magma(4)[2], "1" = magma(4)[3], "2" = "#ffffff")) +
  theme(legend.position="none") -> tax_sunburst

tax_sunburst

annotation_df <- tax_sunburst$data %>%
  dplyr::select(x, y, r0, r, plot_segment, Genus, higherName, size) %>%
  mutate(xend = out_layer_lab_pos(x, r, 0.1),
         yend = out_layer_lab_pos(y, r, 0.1),
         hjust_o = (sign(x)*-1+1)*0.5,
         vjust_o = (sign(y)*-1+1)*0.5,
         hjust_i = (sign(x)+1)*0.5,
         vjust_i = (sign(y)+1)*0.5,
         Genus_exp = str_c("italic('",
                           str_extract(Genus, "^[:alpha:]*"),
                           "')~ (", size, ")"))
annotation_df$Genus_exp[annotation_df$Genus == ""] <- ""
annotation_df$yend[annotation_df$Genus == "Pinus"] <- 2.8
annotation_df$yend[annotation_df$Genus == "Solanum"] <- 2.5
annotation_df$yend[annotation_df$Genus == "Sorghum"] <- 0.2
annotation_df %>%
  filter(higherName != "") -> annot_df_hn

annot_df_hn %>%
  mutate(xend = c(0, 0, 0, 0.5, 0, 1.8, 0, 0.7, -2.2, -1.2, 1.5, 0.6, -2.4, 0, -1.3),
         yend = c(0, -1.52, -1.8, 1.6, -2.1, -0.7, -2.3, 1.95, -0.2, 1.9, -1.1, 2.25, -0.5, -2.52, 2.1),
         hjust_i = c(0.5, 0.5, 0.5, 1, 0.5, 1, 0.5, 1, 0, 0, 1, 1, 0, 0.5, 0),
         vjust_i = c(0.5, 0.5, 0.5, 1, 0.5, 1, 0.5, 1, 1, 1, 1, 1, 1, 0.5, 1)) -> annot_df_hn

# set.seed(98467)
fs <- 1.5
tax_sunburst +
  geom_text(data = annotation_df%>% filter(Genus != ""),
            aes(x = xend,
                y = yend,
                # label=Genus_exp,
                label=str_c("italic('", Genus, "') (", format(size, trim = T, big.mark=","), ")"),
                hjust=hjust_o, vjust = vjust_o),
            color = "black", size = fs,
            parse = T) +
  geom_segment(data = annotation_df %>%filter(plot_segment == 1),
               aes(x = x, y = y,
                   xend = xend,
                   yend = yend,
                   alpha = plot_segment), color = "black") +
  geom_text(data = annot_df_hn,
            aes(x = xend,
                y = yend,
                hjust=hjust_i, vjust = vjust_i,
                label=higherName),
            color = "white", size=fs) +
  geom_segment(data = annot_df_hn %>% filter(y > -1.5),
               aes(x = x, y = y,
                   xend = xend,
                   yend = yend), color = "white") +
  geom_text(data = annotation_df %>%filter(higherName == "Tracheophyta"),
            aes(x = x,
                y = y,
                label=str_c(higherName, " \n(n = ", format(size, trim = T, big.mark=","), ")")),
            hjust = 0.5,
            color = "black",
            size=fs) +
  scale_x_continuous(expand = c(.3, .3)) -> tx

tx
ggsave("./figures/taxonomy_sunburst.png", tx, units = "cm",
       width = 8, height = 5)
ggsave("./figures/taxonomy_sunburst.pdf", tx, units = "in",
       width = 12, height = 7.5)
ggsave("./figures/taxonomy_sunburst.svg", tx, units = "cm",
       width = 8, height = 5)
