library(tidyverse)
library(GGally)

New_Cor_By_Month_Polity <- read_delim("sars-cov2_genomes_clade_region_v3.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(date = lubridate::mdy(date), Month = lubridate::month(date)) %>% 
  group_by(Month, country, location, protein) %>% 
  tally() %>% 
  pivot_wider(names_from = protein, values_from = n, values_fill = 0, names_prefix = "prot_") %>% 
  ungroup() %>% 
  dplyr::filter(Month != 12) %>% 
  dplyr::select(starts_with("prot"))%>% 
  rowwise() %>% 
  mutate(TotalCases = sum(c_across(prot_1A:prot_2C))) %>%
  as.data.frame() %>% 
  mutate_at(vars(prot_1A:prot_2C),list(~round(./TotalCases,4)))  %>% 
  mutate(Type = "Data", Alpha = 1)


Prots <- readRDS("~/Documents/Response_Curves/Simulation_response_probs.rds")%>% 
  mutate(Type = "Simulation", Alpha = 0.3) %>% bind_rows(New_Cor_By_Month_Polity) %>% 
  mutate(Type = fct_relevel(Type, "Simulation"))



ggpairs(Prots,
        lower = list(continuous = GGally::wrap("smooth", 
                                               size=0.1)
                     ),
        columns = c(1:5),
        mapping=ggplot2::aes(colour = Type, alpha = Alpha)
        ) + theme_bw()

ggpairs(Prots,
        columns = c(1:5),
        mapping=ggplot2::aes(colour = Type, alpha = Alpha, size =0.1),
        size = 0.1
) + theme_bw()


