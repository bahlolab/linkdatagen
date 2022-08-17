library(tidyverse)
library(glue)

setwd("~/Documents/linkdatagen/")

hapmap_version <- "2" # or "3"

if(hapmap_version == "2") {
  bed <- read_tsv("hglft_genome_d0ad_be980.bed", 
                  col_names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand"))
}
if(hapmap_version == "3") {
  bed <- read_tsv("hglft_genome_33f8_655670.bed", 
                col_names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand"))
  annotation <- read_tsv("data/annotHapMap3U.txt", comment = "#")
}

bed

annotation <- read_tsv(glue("data/annotHapMap{hapmap_version}U.txt"), comment = "#")


hmu <- inner_join(annotation, bed, by = c("rs_name" = "name"))

# Check when things are different
hmu %>% filter(Chrom != chrom)

hmu %>% filter(Strand != strand)

# filter and mutate to what we can use
flip <- function(x) {
  # the allele on the opposite strand
  case_when(
    x == "A" ~ "T",
    x == "C" ~ "G",
    x == "G" ~ "C",
    x == "T" ~ "A",
  )
}

hmu_1 <- hmu %>% 
  filter(Chrom == chrom & Strand == strand) %>% 
  rename(physical_position_build38 = chromEnd) %>%
  #TODO: filter out non-standard chromosomes
  #mutate(
  #  TOPBOT = if_else(Strand == strand, TOPBOT, if_else(TOPBOT == "T", "B", "T")),
  #  A =      if_else(Strand == strand, A, flip(A)),
  #  B =      if_else(Strand == strand, B, flip(B))
  #) %>% 
  #filter(Chrom == "chr1" & (physical_position_build38 >= 120381544 & physical_position_build38 <= 121605956 ))
  #filter(physical_position_build38 < lead(physical_position_build38) | Chrom != lead(Chrom)) %>% 
  select(Chrom, physical_position_build38, everything())

ggplot(hmu_1, aes(x = physical_position_build37 - lag(physical_position_build37), 
                  y = physical_position_build38 - lag(physical_position_build38))) +
  geom_point()

ggplot(hmu_1, aes(x = physical_position_build37 - lag(physical_position_build37), 
                  y = physical_position_build38 - lag(physical_position_build38))) +
  geom_hex()


ggplot(hmu_1, aes(x = if_else(Chrom == lag(Chrom), physical_position_build37 - lag(physical_position_build37), NA_real_), 
                  y = if_else(Chrom == lag(Chrom), physical_position_build38 - lag(physical_position_build38), NA_real_))) +
  geom_point() +
  xlab("hg19 adjacent SNP physical difference") +
  ylab("hg38 adjacent SNP physical difference")

ggplot(hmu_1, aes(x = if_else(Chrom == lag(Chrom), physical_position_build37 - lag(physical_position_build37), NA_real_), 
                  y = if_else(Chrom == lag(Chrom), physical_position_build38 - lag(physical_position_build38), NA_real_))) +
  geom_point() +
  xlim(0, 1e6) + ylim(-1e6, 1e6) +
  xlab("hg19 adjacent SNP physical difference") +
  ylab("hg38 adjacent SNP physical difference")

ggplot(hmu_1, aes(x = if_else(Chrom == lag(Chrom), physical_position_build37 - lag(physical_position_build37), NA_real_), 
                  y = if_else(Chrom == lag(Chrom), physical_position_build38 - lag(physical_position_build38), NA_real_))) +
  geom_hex() +
  xlim(-10, 1e6) + ylim(-1e6, 1e6)

ggplot(hmu_1, aes(x = if_else(Chrom == lag(Chrom), physical_position_build37 - lag(physical_position_build37), NA_real_), 
                  y = if_else(Chrom == lag(Chrom), physical_position_build38 - lag(physical_position_build38), NA_real_))) +
  geom_hex() +
  xlim(-10, 1e6) + ylim(-1e6, 1e6) +
  scale_fill_gradient(trans = "log") +
  xlab("hg19 adjacent SNP physical difference") +
  ylab("hg38 adjacent SNP physical difference")

ggplot(hmu_1, aes(x = if_else(Chrom == lag(Chrom), physical_position_build37 - lag(physical_position_build37), NA_real_), 
                  y = if_else(Chrom == lag(Chrom), physical_position_build38 - lag(physical_position_build38), NA_real_))) +
  geom_hex() +
  xlim(-10, 1e5) + ylim(-1e5, 1e5) +
  scale_fill_gradient(trans = "log") +
  xlab("hg19 adjacent SNP physical difference") +
  ylab("hg38 adjacent SNP physical difference")


ggplot(hmu_1, aes(x = if_else(Chrom == lag(Chrom), abs(physical_position_build37 - lag(physical_position_build37) - 
                                (physical_position_build38 - lag(physical_position_build38))) + 1, NA_real_))
       ) +
  xlab("hg19 to hg38 absolulte difference + 1") + 
  geom_density() +
  scale_y_log10() +
  scale_x_log10()
  

# Filter out the bad ordering of SNPs. This is probably nowhere near the fastest
# way to do this, but won't have to run this many times (hopefully).
#hmu_2 <- hmu_1
# last_nrow <- nrow(hmu_2) + 1
# while(last_nrow != nrow(hmu_2)) {
#   last_nrow <- nrow(hmu_2)
#   hmu_2 <- hmu_2 %>% 
#     filter(physical_position_build38 < lead(physical_position_build38) | Chrom != lead(Chrom))
#   message("nrow: ", nrow(hmu_2), "\n")
# }

keep <- rep(TRUE, nrow(hmu_1))
max_chrom_position <- 1
for(i in seq(2, nrow(hmu_1), 1)) {
  if(i %% 10000 == 0) {
    message(i)
  }
  if(hmu_1[i-1,]$Chrom != hmu_1[i,]$Chrom) {
    max_chrom_position <- 1
    next
  }
  new_chrom_position <- hmu_1[i,]$physical_position_build38
  if(new_chrom_position < max_chrom_position) {
    keep[i] <- FALSE
  } else {
    max_chrom_position <- new_chrom_position
  }
}
hmu_2 <- hmu_1[keep,]
nrow(hmu_2)

hmu_2 %>% filter((physical_position_build38 > lead(physical_position_build38) | lag(physical_position_build38) > physical_position_build38) & Chrom == lead(Chrom) & lag(Chrom) == Chrom) %>% select(Chrom, physical_position_build38, Strand, A, B, strand) %>% View()
hmu_2 %>% filter((deCODE_genetic_map_position > lead(deCODE_genetic_map_position) | lag(deCODE_genetic_map_position) > deCODE_genetic_map_position) & Chrom == lead(Chrom) & lag(Chrom) == Chrom) %>% select(Chrom, physical_position_build38, Strand, A, B, strand) %>% View()

hmu_3 <- hmu_2 %>%
  select(-physical_position_build37, -chrom, -chromStart, -score) %>% 
  select(-strand)

write_tsv(hmu_3, glue("annotHapMap{hapmap_version}U_hg38.txt"))
