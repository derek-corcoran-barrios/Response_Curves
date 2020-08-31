library(tidyverse)
library(corrr)

New_Cor_By_Month_Polity <- read_delim("sars-cov2_genomes_clade_region_v3.txt", 
                                                             "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(date = lubridate::mdy(date), Month = lubridate::month(date)) %>% 
  group_by(Month, country, location, protein) %>% 
  tally() %>% 
  pivot_wider(names_from = protein, values_from = n, values_fill = 0, names_prefix = "prot_") %>% 
  ungroup() %>% 
  dplyr::filter(Month != 12) %>% 
  dplyr::select(starts_with("prot")) %>% 
  correlate(diagonal = 1) %>% 
  dplyr::select(-rowname) %>% 
  as.matrix()

rownames(New_Cor_By_Month_Polity) <- colnames(New_Cor_By_Month_Polity)


Cor3 <-read_delim("sars-cov2_genomes_clade_region_v3.txt", 
                  "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(date = lubridate::mdy(date), Month = lubridate::month(date)) %>% 
  group_by(Month, country, location, protein) %>% 
  tally() %>% 
  rename(adm0InCovid =  country) %>% 
  pivot_wider(names_from = protein, values_from = n, values_fill = 0) %>% 
  ungroup() %>%  
  pivot_longer(`1A`:`2C`, names_to = "protein", values_to = "n") %>% 
  dplyr::select(n, protein) %>% 
  table() %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = n, values_from = Freq) %>% 
  rowwise() %>% 
  mutate(TotalCases = sum(c_across(`0`:`20`))) %>% 
  mutate_at(vars(`0`:`20`),list(~./TotalCases)) %>% 
  pivot_longer(`0`:`20`, names_to = "Count", values_to = "Probability") %>% 
  group_by(protein) %>% 
  mutate(Probability = cumsum(Probability)) %>% 
  dplyr::select(-TotalCases) %>% 
  pivot_wider(names_from = protein, values_from = Probability)


library(GenOrd)


prob_1A <- Cor3$`1A`[Cor3$`1A` < 1]
prob_1B <- Cor3$`1B`[Cor3$`1B` < 1]
prob_2A <- Cor3$`2A`[Cor3$`2A` < 1]
prob_2B <- Cor3$`2B`[Cor3$`2B` < 1]
prob_2C <- Cor3$`2C`[Cor3$`2C` < 1]

marginal <- list(prob_1A, prob_1B, prob_2A, prob_2B, prob_2C)

Sup <-  list(0:(length(prob_1A)),0:(length(prob_1B)),0:(length(prob_2A)),0:(length(prob_2B)), 0:(length(prob_2C)))

corrcheck(marginal, support = Sup)


ordsample2 <- function (n, marginal, Sigma, support = list(), Spearman = FALSE, 
          cormat = "discrete") 
{
  if (!all(unlist(lapply(marginal, function(x) (sort(x) == 
                                                x & min(x) > 0 & max(x) < 1))))) 
    stop("Error in assigning marginal distributions!")
  if (!isSymmetric(Sigma) | min(eigen(Sigma)$values) < 0 | 
      !all(diag(Sigma) == 1)) 
    stop("Correlation matrix not valid!")
  k <- length(marginal)
  kj <- numeric(k)
  len <- length(support)
  for (i in 1:k) {
    kj[i] <- length(marginal[[i]]) + 1
    if (len == 0) {
      support[[i]] <- 1:kj[i]
    }
  }
  if (cormat == "discrete") {
    Sigmac <- ordcont(marginal = marginal, Sigma = Sigma, 
                      support = support, Spearman = Spearman)[[1]]
    Sigma <- Sigmac
  }
  valori <- mvrnorm(n, rep(0, k), Sigma)
  if (n == 1) 
    valori <- matrix(valori, nrow = 1)
  for (i in 1:k) {
    valori[, i] <- as.integer(cut(valori[, i], breaks = c(min(valori[, 
                                                                     i]) - 1, unique(qnorm(marginal[[i]])), max(valori[, i]) + 
                                                            1)))
    valori[, i] <- support[[i]][valori[, i]]
  }
  return(valori)
}

m <- ordsample2(n = 10000, marginal, Sigma = New_Cor_By_Month_Polity, support = Sup)

m2 <- m %>% as.data.frame()  %>% rowwise() %>% 
  mutate(TotalCases = sum(c_across(prot_1A:prot_2C)))%>% 
  dplyr::filter(TotalCases != 0) %>% 
  mutate_at(vars(prot_1A:prot_2C),list(~./TotalCases)) %>% 
  ungroup() 

m2[,-6] %>% GGally::ggpairs(., 
                            lower = list(continuous = GGally::wrap("smooth", alpha = 0.1, size=0.1)))

saveRDS(m2, "Simulation_response_probs.rds")
