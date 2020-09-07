##########################################
#############Version 2####################
##########################################

library(tidyverse)
library(caret)
library(gbm)
library(MuMIn)
library(furrr)
library(mgcv)


## Read in and generate climatic means

df <- read_rds("df2_corrected_prot_up_low_bound_small.rds")


Prots <- read_rds("Simulation_response_probs.rds") %>% 
  rename(Est_1A = prot_1A, Est_1B = prot_1B, Est_2A = prot_2A, Est_2B = prot_2B, Est_2C = prot_2C)
Weight <-  read_rds("DF_Prop_gbm_Est.rds") %>% 
  mutate(ID = 1:n(),Weight = 1/(R2 - R2_Post)^3, Total = sum(Weight), Std_Weight = Weight/Total) %>% 
  pull(Std_Weight)


Models <- read_rds("models_Prop_gbm_Est.rds") 

Prot = "2C"

Smoothed_Heat <- function(Data, n, Variable, Weight, Models, ncores = 3, Prot, Simulation){
  Sequence <- seq(from = min(df[Variable]), max(df[Variable]), length.out = n)
  
  Vars <- c("rhLag14mean", "TempLag14mean", "uvLag14max", "propOld", "lambda7")
  
  Vars <- Vars[Vars != Variable]
  
  dfMean <- Data %>% 
    dplyr::select(all_of(Vars)) %>% 
    summarise_all(mean)
  
  
  dfVar <- data.frame(V1 = Sequence)
  colnames(dfVar) <- Variable
  
  df <- bind_cols(dfVar, dfMean) %>% 
    group_split(.data[[Variable]]) %>% 
    purrr::map(~bind_cols(.x, Simulation)) %>% 
    reduce(bind_rows)
  
  message("Starting predictions, this might take a while")
  
  plan("multiprocess", workers = ncores)
  
  df$Lambda<- Models %>% 
    furrr::future_map(~predict(.x, df)) %>% 
    furrr::future_map2(.y =Weight,~.x*.y) %>% 
    reduce(`+`)
  
  gc()
  
  message("Finished predictions!")
  
  df <-  df %>% 
    pivot_longer(starts_with("Est"), names_to = "Protein", values_to = "Proportion") %>% 
    mutate(Protein = str_remove_all(Protein, "Est_"))  %>% 
    dplyr::select(Protein, Proportion, Lambda, all_of(Variable)) %>% 
    dplyr::filter(Protein == Prot)
  
  Formula <- paste0("Lambda ~ s(Proportion) + s(", Variable, ")")
  Loess3d <- gam(as.formula(Formula), data = df)
  
  
  newDF <- expand.grid(Proportion = seq(from = min(df$Proportion), to = max(df$Proportion), length.out = 200),
                       V1 = seq(from = min(df[Variable]), to = max(df[Variable]), length.out = 200))
  
  colnames(newDF)[2] <- Variable
  
  newDF$Pred <- predict(Loess3d, newDF)
  
  G <- ggplot(newDF, aes_string(x = "Proportion", y = Variable)) + geom_raster(aes(fill = Pred)) + theme_bw() + scale_fill_viridis_c()
  print(G)
  return(list(DF = newDF, Graph = G, GAM = Loess3d))
}


Test <- Smoothed_Heat(Data = df, 
                      n = 20,
                      Variable = "rhLag14mean",
                      Weight = Weight,
                      Models = Models,
                      ncores = 3,
                      Prot = "2C",
                      Simulation = Prots)

beepr::beep(8)

Test2 <- Smoothed_Heat(Data = df, 
                      n = 20,
                      Variable = "propOld",
                      Weight = Weight,
                      Models = Models,
                      ncores = 3,
                      Prot = "1A",
                      Simulation = Prots)

beepr::beep(8)
