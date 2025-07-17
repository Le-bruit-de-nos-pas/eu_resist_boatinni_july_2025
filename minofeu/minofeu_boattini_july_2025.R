library(readxl)
library(tidyverse)
library(data.table)
options(scipen = 999)


# Look-up species-antibacterials combinations resistance cutoffs ------------


# MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)
# 
# names(MinoFeu_Complete_dataset)
# 
# Look_up <- MinoFeu_Complete_dataset %>% select(`Species identification`, `Ampicillin MIC`:`Doxycycline inhib zone diam`) 
#   
# Look_up <- Look_up %>% gather(Abx, Rx, `Ampicillin MIC`:`Doxycycline inhib zone diam`) %>% 
#   drop_na() %>% select(-Rx) %>% distinct()
# 
# length(unique(Look_up$`Species identification`))
# 
# length(unique(Look_up$Abx))
# 
# Look_up <- Look_up %>% mutate(Abx=str_replace(Abx, " MIC", "")) %>%
#   mutate(Abx=str_replace(Abx, " inhib zone diam", "")) %>% distinct()
# 
# fwrite(Look_up, "Look_up_MinoFeu_Complete_dataset.csv")

# ----------------------      

# Resistant/susceptible flags using MIC --------------------------------

Look_up <- fread("../data/Look_up_MinoFeu_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MinoFeu_Complete_dataset$`Species identification`))

sort(unique(MinoFeu_Complete_dataset$`Species identification`))

names(MinoFeu_Complete_dataset)


MIC_data <- MinoFeu_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MinoFeu_Complete_dataset %>% select(contains(" MIC"))
  )

names(MIC_data)

MIC_data <- MIC_data %>% select(-c(`AST by broth microdilution=1`)) %>%
  gather(Abx, MIC, `Ampicillin MIC`:`Doxycycline MIC`) 


MIC_data <- MIC_data %>% mutate(Abx=str_replace(Abx, " MIC", "")) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>%
  mutate(MIC = as.numeric(MIC))

names(Look_up)

MIC_data <- MIC_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_disc_lower, CLSI_disc_lower_m100, CLSI_disc_lower_m45, CASFM_disc_lower)))

length(unique(MIC_data$`Code event`))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_m100_Resist=ifelse(MIC>=CLSI_mic_bigeq_m100 ,1,0)) %>%
  mutate(CLIST_mic_m45_Resist=ifelse(MIC>=CLSI_mic_bigeq_m45 ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


summary_MIC_concent <- MIC_data %>% filter(!is.na(MIC)) %>%
  group_by(`Species identification`, Abx) %>% 
  summarise(mean=mean(MIC, na.rm=T), 
            sd=sd(MIC), 
            mean=mean(MIC), 
            median=median(MIC), 
            Q90=quantile(MIC, 0.90),
            n=n()) 

fwrite(summary_MIC_concent, "../out/summary_MIC_concent_Jul_16_MIC90.csv")



EUCAST_resist_counts <- MIC_data %>% filter(!is.na(EUCAST_mic_Resist )) %>%
  group_by(`Species identification`, Abx, EUCAST_mic_Resist) %>% count() %>%
  spread(key=EUCAST_mic_Resist, value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_eucast=`1`+`0`, perc_r_eucast=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

CLSI_m100_resist_counts <- MIC_data %>% filter(!is.na(CLIST_mic_m100_Resist )) %>%
  group_by(`Species identification`, Abx, CLIST_mic_m100_Resist  ) %>% count() %>%
  spread(key=CLIST_mic_m100_Resist  , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_clsi=`1`+`0`, perc_r_clsi=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

CLSI_m45_resist_counts <- MIC_data %>% filter(!is.na(CLIST_mic_m45_Resist )) %>%
  group_by(`Species identification`, Abx, CLIST_mic_m45_Resist  ) %>% count() %>%
  spread(key=CLIST_mic_m45_Resist  , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_clsi=`1`+`0`, perc_r_clsi=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

casfm_resist_counts <- MIC_data %>% filter(!is.na(CASFM_mic_Resist  )) %>%
  group_by(`Species identification`, Abx, CASFM_mic_Resist  ) %>% count() %>%
  spread(key=CASFM_mic_Resist  , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_casfm=`1`+`0`, perc_r_casfm=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))


fwrite(MIC_data, "../out/MIC_data_July_16.csv")
fwrite(MIC_data, "MIC_data_Jan29.csv")

fwrite(EUCAST_resist_counts, "../out/EUCAST_resist_counts_Jul_16.csv")
fwrite(CLSI_m100_resist_counts, "../out/CLSI_m100_resist_counts_Jul_16.csv")
fwrite(CLSI_m45_resist_counts, "../out/CLSI_m45_resist_counts_Jul_16.csv")
fwrite(casfm_resist_counts, "../out/casfm_resist_counts_Jul_16.csv")

# ---------------------------
# Resistant/susceptible flags using inhibitory zone diameter --------------------------------

Look_up <- fread("../data/Look_up_MinoFeu_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MinoFeu_Complete_dataset$`Species identification`))

sort(unique(MinoFeu_Complete_dataset$`Species identification`))

names(MinoFeu_Complete_dataset)


ZONE_data <- MinoFeu_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MinoFeu_Complete_dataset %>% select(contains(" inhib zone"))
  )

names(ZONE_data)

ZONE_data <- ZONE_data %>%
  gather(Abx, ZONE, `Ampicillin inhib zone diam`:`Doxycycline inhib zone diam`) 

ZONE_data <- ZONE_data %>% mutate(Abx=str_replace(Abx, " inhib zone diam", "")) 

ZONE_data$ZONE <- as.numeric(ZONE_data$ZONE)

unique(ZONE_data$ZONE)

ZONE_data <- ZONE_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_mic_big , CLSI_mic_bigeq_m100, CLSI_mic_bigeq_m45 , CASFM_mic_big )))


ZONE_data <- ZONE_data %>% 
  mutate(EUCAST_Diam_Resist=ifelse(ZONE<EUCAST_disc_lower,1,0)) %>%
  mutate(CLIST_Diam_m100_Resist=ifelse(ZONE<CLSI_disc_lower_m100  ,1,0)) %>%
  mutate(CLIST_Diam_m45_Resist=ifelse(ZONE<CLSI_disc_lower_m45  ,1,0)) %>%
  mutate(CASFM_Diam_Resist=ifelse(ZONE<CASFM_disc_lower ,1,0)) 


summary_ZONE_concent <- ZONE_data %>% filter(!is.na(ZONE)) %>%
  group_by(`Species identification`, Abx) %>% 
  summarise(mean=mean(ZONE, na.rm=T), 
            sd=sd(ZONE), 
            mean=mean(ZONE), 
            median=median(ZONE), 
            Q1=quantile(ZONE, 0.25),
            Q3=quantile(ZONE, 0.75),
            n=n()) 


EUCAST_resist_counts_zone <- ZONE_data %>% filter(!is.na(EUCAST_Diam_Resist)) %>%
  group_by(`Species identification`, Abx, EUCAST_Diam_Resist) %>% count() %>%
  spread(key=EUCAST_Diam_Resist, value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_eucast=`1`+`0`, perc_r_eucast=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

clsi_m100_resist_counts_zone <- ZONE_data %>% filter(!is.na(CLIST_Diam_m100_Resist   )) %>%
  group_by(`Species identification`, Abx, CLIST_Diam_m100_Resist  ) %>% count() %>%
  spread(key=CLIST_Diam_m100_Resist  , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_clsi=`1`+`0`, perc_r_clsi=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

clsi_m45_resist_counts_zone <- ZONE_data %>% filter(!is.na(CLIST_Diam_m45_Resist   )) %>%
  group_by(`Species identification`, Abx, CLIST_Diam_m45_Resist  ) %>% count() %>%
  spread(key=CLIST_Diam_m45_Resist  , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_clsi=`1`+`0`, perc_r_clsi=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))

casfm_resist_counts_zone <- ZONE_data %>% filter(!is.na(CASFM_Diam_Resist  )) %>%
  group_by(`Species identification`, Abx, CASFM_Diam_Resist ) %>% count() %>%
  spread(key=CASFM_Diam_Resist , value=n) %>%
  mutate(`1`=ifelse(is.na(`1`),0,`1`)) %>%
  mutate(`0`=ifelse(is.na(`0`),0,`0`)) %>%
  mutate(n_casfm=`1`+`0`, perc_r_casfm=`1`/(`1`+`0`)) %>% select(-c(`1`,`0`))


fwrite(summary_ZONE_concent, "../out/summary_ZONE_concent_July_16.csv")

fwrite(ZONE_data, "../out/ZONE_data_July_16.csv")

fwrite(EUCAST_resist_counts_zone, "../out/EUCAST_resist_counts_zone_July_16.csv")
fwrite(clsi_m100_resist_counts_zone, "../out/clsi_m100_resist_counts_zone_July_16.csv")
fwrite(clsi_m45_resist_counts_zone, "../out/clsi_m45_resist_counts_zone_July_16.csv")
fwrite(casfm_resist_counts_zone, "../out/casfm_resist_counts_zone_July_16.csv")





# ---------------------


# Overall figures ---------------------------

MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

Look_up <- fread("../data/Look_up_MinoFeu_Complete_dataset.csv")

# isolates per year 

MinoFeu_Complete_dataset %>% select(`Code event`, `Species identification`, `2020=1`:`2024=1`) %>%
  filter(!is.na(`Species identification`)) %>%
  gather(Group, exp, `2020=1`:`2024=1`) %>%
  filter(exp==1) %>% group_by(`Code event`) %>% filter(Group==max(Group)) %>% ungroup() %>%
  group_by(Group) %>% count() %>% mutate(n=n/4402)

#   Group      n
# 1 2020=1 0.159
# 2 2021=1 0.216
# 3 2022=1 0.214
# 4 2023=1 0.211
# 5 2024=1 0.200

# Isolates per country

data.frame(MinoFeu_Complete_dataset %>% select(`Code event`, `Species identification`, `Austria`:`Turkey`) %>%
             filter(!is.na(`Species identification`)) %>%
             gather(Group, exp, `Austria`:`Turkey`) %>%
             filter(exp==1) %>%
             group_by(Group) %>% count() %>% 
             arrange(n) %>% mutate(n=n/4420) ) 

# 1        Slovakia 0.0004524887
# 2         Iceland 0.0065610860
# 3        Bulgaria 0.0138009050
# 4         Austria 0.0156108597
# 5         Romania 0.0169683258
# 6        Portugal 0.0201357466
# 7          Norway 0.0210407240
# 8     Switzerland 0.0212669683
# 9         Hungary 0.0233031674
# 10         Sweden 0.0253393665
# 11         Turkey 0.0253393665
# 12        Denmark 0.0260180995
# 13         Poland 0.0289592760
# 14        Belgium 0.0309954751
# 15       Slovenia 0.0314479638
# 16        Ireland 0.0330316742
# 17    Netherlands 0.0355203620
# 18 Czech Republic 0.0391402715
# 19          Spain 0.0416289593
# 20        Croatia 0.0615384615
# 21         Greece 0.0710407240
# 22        Germany 0.0762443439
# 23          Italy 0.0907239819
# 24         France 0.2438914027

# Medical ward

data.frame(MinoFeu_Complete_dataset %>% select(`Code event`, `Species identification`,  `Emergency=1`:`Neonatal ICU=1`) %>%
             filter(!is.na(`Species identification`)) %>%
             gather(Group, exp, `Emergency=1`:`Neonatal ICU=1`) %>%
             filter(exp==1) %>% distinct() %>% 
             group_by(Group) %>%
             count() %>%
             arrange(n) %>% mutate(n=n/4313))


#                    Group           n
# 1 Surgical paediatrics=1 0.003941572
# 2         Neonatal ICU=1 0.004869001
# 3       Paediatric ICU=1 0.022722003
# 4  Medical Paediatrics=1 0.056573151
# 5        Surgical ward=1 0.099698586
# 6            Emergency=1 0.112682588
# 7                  ICU=1 0.207744030
# 8         Medical ward=1 0.491769070



# Most Common species

Look_up

TOP50 <- data.frame(MinoFeu_Complete_dataset %>%  filter(!is.na(`Species identification`)) %>% 
                      group_by(`Species identification`) %>% count() %>%
                      arrange(-n) %>% ungroup() %>% mutate(cum=cumsum(n)) %>% mutate(n=n/4420)) 
  

top_20_species <- TOP50 %>% slice(1:20) %>% select(Species.identification) 

top_20_species <- list(top_20_species$Species.identification)[[1]]

# Most Common species by country TOP10 per country

temp <- data.frame(MinoFeu_Complete_dataset %>%  filter(!is.na(`Species identification`)) %>% 
                     gather(Group, exp, `Austria`:`Turkey`) %>%
                     filter(exp==1) %>%
                     mutate(`Species identification`=ifelse(`Species identification`%in%top_20_species, `Species identification`, "Other")) %>%
                     group_by(Group, `Species identification`) %>% count() %>%
                     rename("sub"="n") %>% ungroup() %>%
                     group_by(Group) %>% mutate(Tot=sum(sub)) %>%
                     mutate(perc=sub/Tot) %>%
                     arrange(Group, -perc ) %>%
                     group_by(Group) %>%  select(-c(sub, Tot)) %>%
                     spread(key=Group, value=perc))


fwrite(temp, "../out/temp.csv")

temp <- data.frame(MinoFeu_Complete_dataset %>%  filter(!is.na(`Species identification`)) %>%
              gather(Group, exp, `Emergency=1`:`Neonatal ICU=1`) %>%
             filter(exp==1) %>% distinct() %>% 
             mutate(`Species identification`=ifelse(`Species identification`%in%top_20_species, `Species identification`, "Other")) %>%
             group_by(Group,`Species identification`) %>%
             count() %>%
             rename("sub"="n") %>% ungroup() %>%
             group_by(Group) %>% mutate(Tot=sum(sub)) %>%
             mutate(perc=sub/Tot) %>%
             arrange(Group, -perc ) %>%
             group_by(Group) %>% select(-c(sub, Tot)) %>%
             spread(key=Group, value=perc))

fwrite(temp, "../out/temp.csv")



# --------------------

# Summary Table MIC and Diameters ----------------

MIC_workbook_jul_16 <- read_excel(path = "../out/MIC_workbook_jul_16.xlsx",  sheet="Summary MIC Values", skip = 1)
unique(MIC_workbook_jul_16$`Species identification`)
TOP <-  MIC_workbook_jul_16 

MIC_workbook_jul_16 <- read_excel(path = "../out/MIC_workbook_jul_16.xlsx",  sheet="EUCAST Resist Thresh", skip = 0)
MIC_workbook_jul_16 <- MIC_workbook_jul_16 %>% rename("N_EUCASAT_thre"="# Nr Isolates")

TOP <- TOP %>% rename("Abx"="Abx MIC tested") %>% left_join(MIC_workbook_jul_16)

MIC_workbook_jul_16 <- read_excel(path = "../out/MIC_workbook_jul_16.xlsx",  sheet="CLSI 100 Resist Thresh", skip = 0)
MIC_workbook_jul_16 <- MIC_workbook_jul_16 %>% rename("N_CLSI_thre_m100"="# Nr Isolates") 

TOP <- TOP %>% left_join(MIC_workbook_jul_16)

MIC_workbook_jul_16 <- read_excel(path = "../out/MIC_workbook_jul_16.xlsx",  sheet="CLSI 45 Resist Thresh", skip = 0)
MIC_workbook_jul_16 <- MIC_workbook_jul_16 %>% rename("N_CLSI_thre_m45"="# Nr Isolates") 

TOP <- TOP %>% left_join(MIC_workbook_jul_16)

MIC_workbook_jul_16 <- read_excel(path = "../out/MIC_workbook_jul_16.xlsx",  sheet="CASFM Resist Thresh", skip = 0)
MIC_workbook_jul_16 <- MIC_workbook_jul_16 %>% rename("N_CASFM_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(MIC_workbook_jul_16 )

fwrite(TOP, "../out/MIC_Summary_All_Jul_16.csv")




diam_inhib_workbook_jul_16 <- read_excel(path = "../out/diam_inhib_workbook_jul_16.xlsx",  sheet="Summary Zone Diam Values", skip = 1)

TOP <- diam_inhib_workbook_jul_16

diam_inhib_workbook_jul_16  <- read_excel(path = "../out/diam_inhib_workbook_jul_16.xlsx",  sheet="EUCAST Resist Thresh", skip = 0)
diam_inhib_workbook_jul_16 <- diam_inhib_workbook_jul_16 %>% rename("N_EUCASAT_thre"="# Nr Isolates")

TOP <- TOP %>% rename("Abx"="Abx Zone tested") %>% left_join(diam_inhib_workbook_jul_16 )

diam_inhib_workbook_jul_16  <- read_excel(path = "../out/diam_inhib_workbook_jul_16.xlsx",  sheet="CLSI m100 Resist Thresh", skip = 0)
diam_inhib_workbook_jul_16 <- diam_inhib_workbook_jul_16 %>% rename("N_CLSI_thre_m100"="# Nr Isolates")

TOP <- TOP %>%left_join(diam_inhib_workbook_jul_16 )


diam_inhib_workbook_jul_16  <- read_excel(path = "../out/diam_inhib_workbook_jul_16.xlsx",  sheet="CLSI m45 Resist Thresh", skip = 0)
diam_inhib_workbook_jul_16 <- diam_inhib_workbook_jul_16 %>% rename("N_CLSI_thre_m45"="# Nr Isolates")

TOP <- TOP %>%left_join(diam_inhib_workbook_jul_16 )

diam_inhib_workbook_jul_16 <- read_excel(path = "../out/diam_inhib_workbook_jul_16.xlsx",  sheet="CASFM Resist Thresh", skip = 0)
diam_inhib_workbook_jul_16 <- diam_inhib_workbook_jul_16 %>% rename("N_CASFM_thre"="# Nr Isolates")

TOP <- TOP %>% left_join(diam_inhib_workbook_jul_16)

data.frame(TOP)

fwrite(TOP, "../out/Diams_Summary_All_Jul_16.csv")



# -------------
# EUCAST Plot species vs resistance rate MIC ---------------------


Look_up <- fread("../data/Look_up_MinoFeu_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MinoFeu_Complete_dataset$`Species identification`))

sort(unique(MinoFeu_Complete_dataset$`Species identification`))

names(MinoFeu_Complete_dataset)


MIC_data <- MinoFeu_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MinoFeu_Complete_dataset %>% select(contains(" MIC"))
  )

names(MIC_data)

MIC_data <- MIC_data %>% select(-c(`AST by broth microdilution=1`)) %>%
  gather(Abx, MIC, `Ampicillin MIC`:`Doxycycline MIC`) 


MIC_data <- MIC_data %>% mutate(Abx=str_replace(Abx, " MIC", "")) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>%
  mutate(MIC = as.numeric(MIC))

names(Look_up)

MIC_data <- MIC_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_disc_lower, CLSI_disc_lower_m100, CLSI_disc_lower_m45, CASFM_disc_lower)))

length(unique(MIC_data$`Code event`))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_m100_Resist=ifelse(MIC>=CLSI_mic_bigeq_m100 ,1,0)) %>%
  mutate(CLIST_mic_m45_Resist=ifelse(MIC>=CLSI_mic_bigeq_m45 ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


temp <- MIC_data %>% select(`Code event`, `Species identification`, Abx, EUCAST_mic_Resist)

library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_mic_Resist) * 100
  ) %>% filter(n_samples>10) 


heatmap_data <- resistance_summary %>% ungroup() %>% 
  select(`Species identification`, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$`Species identification`

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "white",  # Color for missing values
                 number_color = "black",  # Set label numbers to black
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "EUCAST MIC \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >10 Samples] \n")
                 
ggsave(file="../out/dendo_10plus_mic_eucast.svg", plot=plot, width=7, height=7)

resistance_summary <- temp %>%
  select(`Code event`, `Species identification`, Abx, EUCAST_mic_Resist) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_mic_Resist) * 100
  ) %>% filter(n_samples>10) 


plot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "EUCAST MIC \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >10 Samples]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="../out/bubble_10plus_mic_eucast.svg", plot=plot, width=8, height=8)




# ----------------
# CLSI m100 Plot species vs resistance rate MIC ---------------------


Look_up <- fread("../data/Look_up_MinoFeu_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MinoFeu_Complete_dataset$`Species identification`))

sort(unique(MinoFeu_Complete_dataset$`Species identification`))

names(MinoFeu_Complete_dataset)


MIC_data <- MinoFeu_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MinoFeu_Complete_dataset %>% select(contains(" MIC"))
  )

names(MIC_data)

MIC_data <- MIC_data %>% select(-c(`AST by broth microdilution=1`)) %>%
  gather(Abx, MIC, `Ampicillin MIC`:`Doxycycline MIC`) 


MIC_data <- MIC_data %>% mutate(Abx=str_replace(Abx, " MIC", "")) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>%
  mutate(MIC = as.numeric(MIC))

names(Look_up)

MIC_data <- MIC_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_disc_lower, CLSI_disc_lower_m100, CLSI_disc_lower_m45, CASFM_disc_lower)))

length(unique(MIC_data$`Code event`))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_m100_Resist=ifelse(MIC>=CLSI_mic_bigeq_m100 ,1,0)) %>%
  mutate(CLIST_mic_m45_Resist=ifelse(MIC>=CLSI_mic_bigeq_m45 ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


temp <- MIC_data %>% select(`Code event`, `Species identification`, Abx, CLIST_mic_m100_Resist)

library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CLIST_mic_m100_Resist) * 100
  ) %>% filter(n_samples>10) 


heatmap_data <- resistance_summary %>% ungroup() %>% 
  select(`Species identification`, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$`Species identification`

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "white",  # Color for missing values
                 number_color = "black",  # Set label numbers to black
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "CLSI m100 MIC \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >10 Samples] \n")
                 
ggsave(file="../out/dendo_10plus_mic_clsi_m100.svg", plot=plot, width=7, height=7)

resistance_summary <- temp %>%
  select(`Code event`, `Species identification`, Abx, CLIST_mic_m100_Resist) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CLIST_mic_m100_Resist) * 100
  ) %>% filter(n_samples>10) 


plot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "CLSI m100 MIC \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >10 Samples]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="../out/bubble_10plus_mic_clsi_m100.svg", plot=plot, width=8, height=8)




# ----------------
# CASFM  Plot species vs resistance rate MIC ---------------------


Look_up <- fread("../data/Look_up_MinoFeu_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MinoFeu_Complete_dataset$`Species identification`))

sort(unique(MinoFeu_Complete_dataset$`Species identification`))

names(MinoFeu_Complete_dataset)


MIC_data <- MinoFeu_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MinoFeu_Complete_dataset %>% select(contains(" MIC"))
  )

names(MIC_data)

MIC_data <- MIC_data %>% select(-c(`AST by broth microdilution=1`)) %>%
  gather(Abx, MIC, `Ampicillin MIC`:`Doxycycline MIC`) 


MIC_data <- MIC_data %>% mutate(Abx=str_replace(Abx, " MIC", "")) 

unique(MIC_data$MIC)

MIC_data <- MIC_data %>%
  mutate(MIC = as.numeric(MIC))

names(Look_up)

MIC_data <- MIC_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_disc_lower, CLSI_disc_lower_m100, CLSI_disc_lower_m45, CASFM_disc_lower)))

length(unique(MIC_data$`Code event`))

MIC_data <- MIC_data %>% 
  mutate(EUCAST_mic_Resist=ifelse(MIC>EUCAST_mic_big,1,0)) %>%
  mutate(CLIST_mic_m100_Resist=ifelse(MIC>=CLSI_mic_bigeq_m100 ,1,0)) %>%
  mutate(CLIST_mic_m45_Resist=ifelse(MIC>=CLSI_mic_bigeq_m45 ,1,0)) %>%
  mutate(CASFM_mic_Resist=ifelse(MIC>CASFM_mic_big,1,0)) 


temp <- MIC_data %>% select(`Code event`, `Species identification`, Abx, CASFM_mic_Resist)

library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CASFM_mic_Resist) * 100
  ) %>% filter(n_samples>10) 


heatmap_data <- resistance_summary %>% ungroup() %>% 
  select(`Species identification`, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$`Species identification`

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "white",  # Color for missing values
                 number_color = "black",  # Set label numbers to black
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "CASFM MIC \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >10 Samples] \n")
                 
ggsave(file="../out/dendo_10plus_mic_casfm.svg", plot=plot, width=7, height=7)

resistance_summary <- temp %>%
  select(`Code event`, `Species identification`, Abx, CASFM_mic_Resist) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CASFM_mic_Resist) * 100
  ) %>% filter(n_samples>10) 


plot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "CASFM MIC \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >10 Samples]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="../out/bubble_10plus_mic_casfm.svg", plot=plot, width=8, height=8)




# ----------------
# EUCAST Plot species vs resistance rate zone inhibition ---------------------


Look_up <- fread("../data/Look_up_MinoFeu_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MinoFeu_Complete_dataset$`Species identification`))

sort(unique(MinoFeu_Complete_dataset$`Species identification`))

names(MinoFeu_Complete_dataset)


ZONE_data <- MinoFeu_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MinoFeu_Complete_dataset %>% select(contains(" inhib zone"))
  )

names(ZONE_data)

ZONE_data <- ZONE_data %>%
  gather(Abx, ZONE, `Ampicillin inhib zone diam`:`Doxycycline inhib zone diam`) 

ZONE_data <- ZONE_data %>% mutate(Abx=str_replace(Abx, " inhib zone diam", "")) 

ZONE_data$ZONE <- as.numeric(ZONE_data$ZONE)



ZONE_data <- ZONE_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_mic_big , CLSI_mic_bigeq_m100, CLSI_mic_bigeq_m45 , CASFM_mic_big )))


ZONE_data <- ZONE_data %>% 
  mutate(EUCAST_Diam_Resist=ifelse(ZONE<EUCAST_disc_lower,1,0)) %>%
  mutate(CLIST_Diam_m100_Resist=ifelse(ZONE<CLSI_disc_lower_m100  ,1,0)) %>%
  mutate(CLIST_Diam_m45_Resist=ifelse(ZONE<CLSI_disc_lower_m45  ,1,0)) %>%
  mutate(CASFM_Diam_Resist=ifelse(ZONE<CASFM_disc_lower ,1,0)) 


temp <- ZONE_data %>% select(`Code event`, `Species identification`, Abx, EUCAST_Diam_Resist )

library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_Diam_Resist) * 100
  ) %>% filter(n_samples>10) 


heatmap_data <- resistance_summary %>% ungroup() %>% 
  select(`Species identification`, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$`Species identification`

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "white",  # Color for missing values
                 number_color = "black",  # Set label numbers to black
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "EUCAST Zone Inhibition \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >10 Samples] \n")
                 
ggsave(file="../out/dendo_10plus_zone_inhib_eucast.svg", plot=plot, width=7, height=7)


resistance_summary <- temp %>%
  select(`Code event`, `Species identification`, Abx, EUCAST_Diam_Resist) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(EUCAST_Diam_Resist) * 100
  ) %>% filter(n_samples>10) 


plot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "EUCAST Zone Inhibition \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >10 Samples]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="../out/bubble_10plus_zone_inhib_eucast.svg", plot=plot, width=8, height=8)




# ----------------
# CASFM  Plot species vs resistance rate zone inhibition ---------------------


Look_up <- fread("../data/Look_up_MinoFeu_Complete_dataset.csv")

sort(unique(Look_up$`Species identification`))

length(unique(Look_up$`Species identification`))

MinoFeu_Complete_dataset <- read_xlsx(path="../data/2.MINOFEu_dataset_matteo_only.xlsx",sheet = "3. Complete dataset", skip=0, col_types = "text", trim_ws = TRUE)

length(unique(MinoFeu_Complete_dataset$`Species identification`))

sort(unique(MinoFeu_Complete_dataset$`Species identification`))

names(MinoFeu_Complete_dataset)


ZONE_data <- MinoFeu_Complete_dataset %>% select(`Code event`,  `Species identification`) %>%
  bind_cols(
    MinoFeu_Complete_dataset %>% select(contains(" inhib zone"))
  )

names(ZONE_data)

ZONE_data <- ZONE_data %>%
  gather(Abx, ZONE, `Ampicillin inhib zone diam`:`Doxycycline inhib zone diam`) 

ZONE_data <- ZONE_data %>% mutate(Abx=str_replace(Abx, " inhib zone diam", "")) 

ZONE_data$ZONE <- as.numeric(ZONE_data$ZONE)



ZONE_data <- ZONE_data %>% 
  left_join(Look_up %>% select(-c( EUCAST_mic_big , CLSI_mic_bigeq_m100, CLSI_mic_bigeq_m45 , CASFM_mic_big )))


ZONE_data <- ZONE_data %>% 
  mutate(EUCAST_Diam_Resist=ifelse(ZONE<EUCAST_disc_lower,1,0)) %>%
  mutate(CLIST_Diam_m100_Resist=ifelse(ZONE<CLSI_disc_lower_m100  ,1,0)) %>%
  mutate(CLIST_Diam_m45_Resist=ifelse(ZONE<CLSI_disc_lower_m45  ,1,0)) %>%
  mutate(CASFM_Diam_Resist=ifelse(ZONE<CASFM_disc_lower ,1,0)) 


temp <- ZONE_data %>% select(`Code event`, `Species identification`, Abx, CASFM_Diam_Resist)

library(pheatmap)
library(dplyr)
library(tidyr)

resistance_summary <- temp %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CASFM_Diam_Resist) * 100
  ) %>% filter(n_samples>10) 


heatmap_data <- resistance_summary %>% ungroup() %>% 
  select(`Species identification`, Abx, resistance_rate) %>%
  pivot_wider(names_from = Abx, values_from = resistance_rate, values_fill = NA)

heatmap_matrix <- as.matrix(heatmap_data[, -1])

rownames(heatmap_matrix) <- heatmap_data$`Species identification`

heatmap_matrix[is.na(heatmap_matrix)] <- -99

display_matrix <- round(heatmap_matrix, 0)

plot <- pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("white", "lightcyan1", "royalblue4"))(50),  
                 cluster_rows = TRUE,  # Cluster species
                 cluster_cols = TRUE,  # Cluster antibiotics
                 na_col = "white",  # Color for missing values
                 number_color = "black",  # Set label numbers to black
                 fontsize_row = 10,  # Font size for species
                 fontsize_col = 10,  # Font size for antibiotics
                 display_numbers = display_matrix,  # Show exact resistance rates
                 main = "CASFM Zone Inhibition \n % Antibiotic Resistance Clustering \n [Species-Abx Combinations With >10 Samples] \n")
                 
ggsave(file="../out/dendo_10plus_zone_inhibition_casfm.svg", plot=plot, width=7, height=7)

resistance_summary <- temp %>%
  select(`Code event`, `Species identification`, Abx, CASFM_Diam_Resist) %>% drop_na() %>%
  group_by(`Species identification`, Abx) %>%
  summarise(
    n_samples = n(),
    resistance_rate = mean(CASFM_Diam_Resist) * 100
  ) %>% filter(n_samples>10) 


plot <- ggplot(resistance_summary, aes(x = Abx, y = `Species identification`)) +
  geom_point(aes(size = n_samples , color = resistance_rate), alpha = 0.9) +
  scale_size(range = c(1, 20)) +  # Adjust the bubble size range
  scale_color_gradient(low = "lightgray", high = "midnightblue") +  # Color scale from susceptible (green) to resistant (red)
  labs(title = "CASFM Zone Inhibition \n % Antibiotic Resistance by Species and Antibiotic \n [Species-Abx Combinations With >10 Samples]",
       x = "Antibiotic",
       y = "Species",
       size = "Sample Size",
       color = "Proportion Resistant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  


ggsave(file="../out/bubble_10plus_zone inhib_casfm.svg", plot=plot, width=8, height=8)




# ----------------