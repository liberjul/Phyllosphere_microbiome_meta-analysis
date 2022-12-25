library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(classInt)
library(viridis)
library(tidyverse)
library(svglite)
library(measurements)

country_from_geo_str <- function(in_str){
  return(unlist(strsplit(in_str, ":"))[1])
}
loc_from_geo_str <- function(in_str){
  return(trimws(unlist(strsplit(in_str, ":"))[2]))
}
to_from_replace_str <- function(x, to, from){
  map = setNames(to, from)
  if (x %in% from){
    return(map[x])
  }
  else {
    return(x)
  }
}
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  # from https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}


country_region_table <- read.csv("./data/biosample_data/country_region_table.csv")
sra_meta_leaf_metag_long <- read.delim("./data/biosample_data/leaf_metagenome_ORGN_sample_metadata_2021_11_09_long.txt", stringsAsFactors = F)
sra_meta_phyllo_metag_long <- read.delim("./data/biosample_data/phyllosphere_metagenome_ORGN_sample_metadata_2021_11_08_long.txt", stringsAsFactors = F)
comb_meta_long <- rbind(sra_meta_leaf_metag_long, sra_meta_phyllo_metag_long)

sra_meta_leaf_metag_long %>%
  distinct() %>%
  pivot_wider(names_from = Attribute_name,
              values_from = Attribute_value,
              # values_fn = length,
              id_cols = Sample_ID) %>%
  dplyr::select(-starts_with(c("jgi", "samn", "phyllo"))) %>%
  mutate(geo_loc_name = replace(geo_loc_name,geo_loc_name=="Russia: Crimea, Yalta", "Ukraine: Crimea, Yalta")) %>%
  mutate(country = sapply(.$geo_loc_name, country_from_geo_str),
         locality = sapply(.$geo_loc_name, loc_from_geo_str)) -> sra_meta_leaf_metag_wide


sra_meta_phyllo_metag_long %>%
  distinct() %>%
  pivot_wider(names_from = Attribute_name,
              values_from = Attribute_value,
              # values_fn = length,
              id_cols = Sample_ID) %>%
  dplyr::select(-starts_with(c("jgi", "samn", "phyllo"))) %>%
  #Cleaning
  mutate(geo_loc_name = replace(geo_loc_name,geo_loc_name=="Salinas Valley", "USA: California, Salinas Valley")) %>%
  mutate(country = sapply(.$geo_loc_name, country_from_geo_str),
         locality = sapply(.$geo_loc_name, loc_from_geo_str)) -> sra_meta_phyllo_metag_wide


colnames(sra_meta_phyllo_metag_wide)

comb_wide <- bind_rows(sra_meta_leaf_metag_wide,
          sra_meta_phyllo_metag_wide)
comb_wide %>%
  group_by(country) %>%
  summarise(count = n()) %>%
  arrange(-count) %>%
  mutate(country = replace(country, country == "USA", "United States"))-> country_table
country_table


convert_lat_long <- function(ll_string, index){
  dir_list <- list(N = 1, S = -1, E = 1, W = -1)
  dec_coord_vec <- c()
  coord_vec <- c()
  dir_vec <- c()
  if (str_detect(ll_string, ".*?\u0081.*")){ # degree sign
    ll_string <- str_replace_all(ll_string, "?\u0081", "")
    ll_string <- str_replace_all(ll_string, "'", "")
    str_lat <- unlist(str_split(ll_string, ", "))[1]
    str_long <- unlist(str_split(ll_string, ", "))[2]
    lat <- conv_unit(str_lat, from = "deg_min_sec", to = "dec_deg")
    long <- conv_unit(str_long, from = "deg_min_sec", to = "dec_deg")
    dec_coord_vec <- as.numeric(c(lat, long))
  }
  else if (str_detect(ll_string, ".*?\u00C1.*")){ # such as "38Á 31.899', -121Á 45.204'"
    ll_string <- str_replace_all(ll_string, "\u00C1", "")
    ll_string <- str_replace_all(ll_string, "'", "")
    str_lat <- unlist(str_split(ll_string, ", "))[1]
    str_long <- unlist(str_split(ll_string, ", "))[2]
    if (str_detect(ll_string, "\\.")){ # degree minute
      lat <- conv_unit(str_lat, from = "deg_dec_min", to = "dec_deg")
      long <- conv_unit(str_long, from = "deg_dec_min", to = "dec_deg")
    }
    else { # degree minute second
      lat <- conv_unit(str_lat, from = "deg_min_sec", to = "dec_deg")
      long <- conv_unit(str_long, from = "deg_min_sec", to = "dec_deg") 
    }
    dec_coord_vec <- as.numeric(c(lat, long))
  }
  else if (! str_detect(ll_string, "\\.")){ 
    if (str_detect(ll_string, "[NS]\\d")){
      lat_dir <- str_extract(ll_string, "[NS]")
      lon_dir <- str_extract(ll_string, "[EW]")
      ll_split <- str_extract_all(ll_string, "[NSEW]\\d[\\d\\s]*\\d", simplify = T)
      for (i in ll_split){
        coord_str <- str_extract(i, "\\d[\\d\\s]*\\d")
        coord <- conv_unit(coord_str, from = "deg_min_sec", to = "dec_deg")
        dir <- str_extract(i, "[NSEW]")
        if (! is.na(coord)){
          coord_vec <- c(coord_vec, coord)
        }
        if (! is.na(dir)){
          dir_vec <- c(dir_vec, dir)
        }
      }
    }
    else { # string without decimals, like 44 N 83 W
      ll_split <- unlist(str_split(ll_string, "\\s"))
      for (i in ll_split){
        coord <- str_extract(i, "\\d{1,3}")
        dir <- str_extract(i, "[NSEW]")
        if (! is.na(coord)){
          coord_vec <- c(coord_vec, coord)
        }
        if (! is.na(dir)){
          dir_vec <- c(dir_vec, dir)
        }
      }
    }
    for (i in 1:2){
      dec_coord <- dir_list[[dir_vec[i]]]*as.numeric(coord_vec[i])
      dec_coord_vec <- c(dec_coord_vec, dec_coord)
    }
    return(dec_coord_vec[index])
  }
  else { # Decimal string like 44.134 N 83.221 W or with directions first
    if (str_detect(ll_string, "[NS]\\d")){
      lat_dir <- str_extract(ll_string, "[NS]")
      lon_dir <- str_extract(ll_string, "[EW]")
      ll_split <- str_extract_all(ll_string, "[NSEW]\\d[\\d\\s\\.]*\\d", simplify = T)
      for (i in ll_split){
        coord_str <- str_extract(i, "\\d[\\d\\s\\.]*\\d")
        coord <- conv_unit(coord_str, from = "deg_min_sec", to = "dec_deg")
        dir <- str_extract(i, "[NSEW]")
        if (! is.na(coord)){
          coord_vec <- c(coord_vec, coord)
        }
        if (! is.na(dir)){
          dir_vec <- c(dir_vec, dir)
        }
      }
      for (i in 1:2){
        dec_coord <- dir_list[[dir_vec[i]]]*as.numeric(coord_vec[i])
        dec_coord_vec <- c(dec_coord_vec, dec_coord)
      }
    }
    else {
      ll_split <- unlist(str_split(ll_string, "\\s"))
      for (i in ll_split){
        coord <- str_extract(i, "\\d{1,3}\\.\\d*")
        dir <- str_extract(i, "[NSEW]")
        if (! is.na(coord)){
          coord_vec <- c(coord_vec, coord)
        }
        if (! is.na(dir)){
          dir_vec <- c(dir_vec, dir)
        }
      }
      for (i in 1:2){
        dec_coord <- dir_list[[dir_vec[i]]]*as.numeric(coord_vec[i])
        dec_coord_vec <- c(dec_coord_vec, dec_coord)
      }
    }
  }
  return(dec_coord_vec[index])
}

ll_string <- "N60 7 44.0 E14 32 1.0"
convert_lat_long(ll_string, 1)
ll_string <- "38Á 31.899', -121Á 45.204'"
convert_lat_long(a, 1)

comb_meta_long %>%
  filter(Attribute_name == "lat_lon",
         ! is.na(Attribute_value),
         ! Attribute_value %in% c("unknown", "missing", "not collected", "Not collected",
                                  "not applicable", "na", "Not applicable", "N/A"),
         str_detect(Attribute_value, ".*\\d.*")) -> lat_lon_data

lat_lon_data$Lat <- sapply(lat_lon_data$Attribute_value, convert_lat_long, index = 1)
lat_lon_data$Long <- sapply(lat_lon_data$Attribute_value, convert_lat_long, index = 2)
lat_lon_data %>%
  filter(is.na(Long))

plot(lat_lon_data$Long, lat_lon_data$Lat)
# unique(lat_lon_data$Attribute_value)

world <- ne_countries(scale = "medium", returnclass = "sf")
world$sample_num <- rep(0,length(world$name))

samples_counts <- country_table %>% filter(!is.na(country))
# colnames(samples_counts) <- c("Country", "Count")
# samples_counts
for (i in samples_counts$country){
  world$sample_num[world$name == i] <- samples_counts$count[samples_counts$country == i]
}
# world$name

world$sample_num[world$sample_num == 0] <- NA
# Code adapted from https://github.com/milos-agathon/neet_2019/blob/main/R/neet_europe.r
ni = classIntervals(world$sample_num, 
                    n = 6, 
                    style = 'quantile')$brks
ni
labels <- c()
for(i in 1:length(ni)){
  labels <- c(labels, paste0(round(ni[i], 0), 
                             "-", 
                             round(ni[i + 1], 0)))
}
labels <- labels[1:length(labels)-1]
labels

world$cat <- cut(world$sample_num, 
                 breaks = ni, 
                 labels = labels, 
                 include.lowest = T)
levels(world$cat) # let's check how many levels it has (6)
# let's check how many levels it has (6)

# label NAs, too
lvl <- levels(world$cat)
lvl[length(lvl) + 1] <- "No data"
world$cat <- factor(world$cat, levels = lvl)
world$cat[is.na(world$cat)] <- "No data"
levels(world$cat)

# cbar_breaks <- c(1, 10, 20, 100, 200, 400, 800, 1600, 3200)
crs_use <- st_crs(world)$input
sample_coords <- st_as_sf(lat_lon_data,
                          coords = c("Long", "Lat"),
                          crs = crs_use) %>%
  st_transform(crs = crs_use)
# check <- st_intersects(x = sample_coords, y = world,
#                        sparse = F)
# 
# lat_lon_data[check[,world$sovereignt == "Russia"], ]

plot.new()
world %>%
  ggplot() + geom_sf(aes(fill=cat), color = NA,
                     alpha = 0.8, na.rm = T) +
  geom_sf(data = sample_coords, 
             color = "black", alpha = 0.2, size = 0.6) +
  scale_fill_manual(name = "Number of collected samples",
                    values = c(viridis(n=6, end = 0.8, direction = -1), "grey80"),
                    labels = c("1-5", "5-14", "14-47", "47-81", "81-184", "184-4697", "No data"),
                    na.value="gray90", drop=F) + 
  guides(fill=guide_legend(
    direction = "horizontal",
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )
  ) +
  theme_minimal() + labs(fill ="Number of samples\ncollected") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(),
        legend.position = c(.55, .04),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "white", size = 0.2),
        plot.title = element_text(size=20, color="#6c2d83", hjust=0.5, vjust=-10),
        plot.subtitle = element_text(size=14, color="#bd5288", hjust=0.5, vjust=-15, face="bold"),
        plot.caption = element_text(size=9, color="grey60", hjust=0.5, vjust=9),
        axis.title.x = element_text(size=7, color="grey60", hjust=0.5, vjust=5),
        legend.text = element_text(size=10, color="grey20"),
        legend.title = element_text(size=11, color="grey20"),
        strip.text = element_text(size=12),
        plot.margin = unit(c(t=0, r=0, b=0, l=0),"lines"), #added these narrower margins to enlarge map
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) -> country_map
country_map
ggsave("./figures/Countries_samples_chloropleth.png", country_map, width=7, height=8.5, dpi = 600, device='png')
ggsave("./figures/Countries_samples_chloropleth.pdf", country_map, width=7, height=8.5, dpi = 600, device='pdf')
ggsave("./figures/Countries_samples_chloropleth.svg", country_map, width=7, height=8.5, dpi = 600, device='svg')
