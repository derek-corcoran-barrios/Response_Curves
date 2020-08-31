library(tidyverse)
library(caret)
library(gbm)
library(MuMIn)
## Read in and generate climatic means

df <- read_rds("df2_corrected_prot_up_low_bound_small.rds") %>% 
  dplyr::select("TempLag14mean", "rhLag14mean", "uvLag14max", "propOld", "lambda7") %>% 
  summarise_all(mean)

## Read in and join with Simulation of prots

Prots <- read_rds("Simulation_response_probs.rds") %>% 
  rename(Est_1A = prot_1A, Est_1B = prot_1B, Est_2A = prot_2A, Est_2B = prot_2B, Est_2C = prot_2C)

df <- bind_cols(df, Prots)

## PRedict GLM

df2 <- df

GLM_Mod <- read_rds("avgmod.95p2.rds")

df2$Pred <- predict(GLM_Mod, newdata = df)

df2 <- df2 %>% 
  pivot_longer(starts_with("Est"), names_to = "Protein", values_to = "Proportion") %>% 
  mutate(Protein = str_remove_all(Protein, "Est_")) %>% 
  dplyr::select(Protein, Proportion, Pred) %>% 
  group_by(Protein, Proportion) %>% 
  summarise_at("Pred", .funs = list(SD = sd, Lambda = mean)) %>% 
  dplyr::filter(!is.na(SD))

## Graph for glm
## Noisy dynamics

ggplot(df2, aes(x = Proportion, y = Lambda)) +
  geom_ribbon(aes(ymax = Lambda + SD, ymin = Lambda - SD), fill = "red", alpha = 0.5) + 
  geom_path() +
  facet_wrap(~Protein) +
  theme_bw()
##########################################################
####GBM##################################################
##########################################################

