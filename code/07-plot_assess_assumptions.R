rm(list=ls())
graphics.off()

# Specify paths ----------------------------------------------------------------

source("code/specify_paths.R")

# Load results -----------------------------------------------------------------

results <- data.table::fread("output/assess_assumptions.csv", data.table = FALSE)

# Format exposure labels -------------------------------------------------------

gwas <- data.table::fread("raw/gwas.csv",
                          stringsAsFactors = FALSE,
                          select = c("trait","trait_long"),
                          data.table = FALSE)

colnames(gwas) <- c("exposure","exposure_name")

results <- merge(results, gwas, by = c("exposure"))

# Format outcome labels -------------------------------------------------------

finn <- data.frame(rbind(c(1,"finn-a-I_INFECT_PARASIT","I","Certain infectious and parasitic diseases"),
                         c(2,"finn-a-II_NEOPLASM","II","Neoplasms"),
                         c(3,"finn-a-III_BLOOD_IMMUN","III","Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism"),
                         c(4,"finn-a-IV_ENDOCRIN_NUTRIT","IV","Endocrine, nutritional and metabolic diseases"),
                         c(5,"finn-a-V_MENTAL_BEHAV","V","Mental and behavioural disorders"),
                         c(6,"finn-a-VI_NERVOUS","VI","Diseases of the nervous system"),
                         c(7,"finn-a-VII_EYE_ADNEXA","VII","Diseases of the eye and adnexa"),
                         c(8,"finn-a-VIII_EAR_MASTOID","VIII","Diseases of the ear and mastoid process"),
                         c(9,"finn-a-IX_CIRCULATORY","IX","Diseases of the circulatory system"),
                         c(10,"finn-a-X_RESPIRATORY","X","Diseases of the respiratory system"),
                         c(11,"finn-a-XI_DIGESTIVE","XI","Diseases of the digestive system"),
                         c(12,"finn-a-XII_SKIN_SUBCUTAN","XII","Diseases of the skin and subcutaneous tissue"),
                         c(13,"finn-a-XIII_MUSCULOSKELET","XIII","Diseases of the musculoskeletal system and connective tissue"),
                         c(14,"finn-a-XIV_GENITOURINARY","XIV","Diseases of the genitourinary system"),
                         c(15,"finn-a-XV_PREGNANCY_BIRTH","XV","Pregnancy, childbirth and the puerperium"),
                         c(16,"finn-a-XVI_PERINATAL","XVI","Certain conditions originating in the perinatal period"),
                         c(17,"finn-a-XVII_MALFORMAT_ABNORMAL","XVII","Congenital malformations, deformations and chromosomal abnormalities"),
                         c(18,"finn-a-XVIII_MISCFINDINGS","XVIII","Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified"),
                         c(19,"finn-a-XIX_INJURY_POISON","XIX","Injury, poisoning and certain other consequences of external causes"),
                         c(20,"finn-a-XX_EXTERNAL_MORB_MORT","XX","External causes of morbidity and mortality"),
                         c(21,"finn-a-XXI_HEALTHFACTORS","XXI","Factors influencing health status and contact with health services")),
                   stringsAsFactors = FALSE)

colnames(finn) <- c("icd10_order","id.outcome","icd10_chapter","outcome_name")
finn$outcome_name <- paste0(finn$icd10_chapter,". ",finn$outcome_name)

results <- merge(results, finn, by = c("id.outcome"))

# Plot -------------------------------------------------------------------------

ggplot2::ggplot(data = results, 
                mapping = ggplot2::aes(y = forcats::fct_rev(exposure_name), 
                                       x = reorder(icd10_chapter,as.numeric(icd10_order)),
                                       color = pval,
                                       size = abs(b))) +
  ggplot2::labs(size = "Absolute value\nof estimate", 
                color = "P-value\n",
                x = "ICD Chapter",
                caption = stringr::str_wrap(paste0("ICD Chapters: ",paste(finn$outcome_name,collapse = "; ")),width = 150)) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::scale_color_gradient(low = "#084594",
                                high = "#9ecae1",
                                na.value = "#eff3ff",
                                guide = "colourbar", 
                                lim = c(0,0.05)) +
  ggplot2::scale_size(breaks = c(0.5,1,2),
                      labels = c(0.5,1,2)) +
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                 panel.grid.minor = ggplot2::element_blank(),
                 axis.text = ggplot2::element_text(size=8),
                 axis.title.y = ggplot2::element_blank(),
                 text = ggplot2::element_text(size=8),
                 legend.position = "bottom",
                 legend.text = ggplot2::element_text(size=8),
                 legend.title = ggplot2::element_text(size=8),
                 plot.caption = ggplot2::element_text(hjust = 0)) 

ggplot2::ggsave(filename = "output/assess_assumptions.jpeg",
                dpi = 300, width = 210, height = 297, unit = "mm", scale = 1)
