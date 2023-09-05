VisualizingDataBrowser
================
Cera Fisher
2023-08-15

This is the Github MD version of the RMarkdown Notebook in the `src` directory of this repo. The below figures can be fully reproduced by cloning this repo and running the Rmd. 

# Making visualizations of hosted datasets

The biggest difficulty with summarizing the data in the CGC’s connected
datasets is the difficulty that always exists with making summaries over
disparate datasets: the datasets are not harmonized and do not use the
same metadata.

Here, I’m trying to get a representation of data volume that will both
be useful to scientists and will also be reproducible as we ingest more
metadata. It seems to me that the things most useful to researchers who
are trying to use our hosted data are \[Fig1\] Disease Type (what kind
of cancer is this data about?) and \[Fig2\] Primary Sites (what part of
the body did the sample come from?).

<figure>
<img src="out/Top25DiseaseTypesInConnectedData.png"
alt="Figure 1 - Summarized Datasets by Disease Type" />
<figcaption aria-hidden="true">Figure 1 - Summarized Datasets by Disease
Type</figcaption>
</figure>

And we didn't use this on the new website, but I also made the primary site figure. 

<figure>
<img src="out/Top25PrimarySitesInConnectedData.png"
alt="Figure 2 - Summarized Datasets by Primary Site" />
<figcaption aria-hidden="true">Figure 2 - Summarized Datasets by Primary Site</figcaption>
</figure>



### Part 1 - Import and Aggregate Data

``` r
library(tidyverse)
library(ggplot2)
```

#### Where did the data tables come from?

For seven of our connected data sources, the metadata were accessed by
using Data Browser on the CGC. These sources are: CPTAC, CPTAC3, FM,
ICGC, TARGET, TCIA, and TCGA. To reproduce these searches and re-pull
the metadata, you must affirmatively select all the options in the given
metadata key you want to pull, and use the ‘Export’ option at the bottom
right of the screen to export to TSV. For example for TCGA, use Data
Browser and affirmatively select every available ‘Primary Site’ and
‘Disease Type’ (Figure 3 below)

<figure>
<img src="img/TCGADataBrowser.png"
alt="Figure 3 - TCGA Data Browser" />
<figcaption aria-hidden="true">Figure 3 - TCGA Data Browser</figcaption>
</figure>


#### What are the data sources represented here?

*CPTAC* and *CPTAC-3* - the Clinical Proteomic Tumor Analysis
Consortium. Data consists of genomics and proteomics.

*FM* - Foundation Medicine Adult Cancer Clinical Dataset.

*ICDC* - Integrated Canine Data Commons

*ICGC* - International Cancer Genome Consortium

*PDC* - Proteomic Data Commons

*TARGET* - Therapeutically Applicable Research to Generate Effective
Treatments

*TCGA* - The Cancer Genome Atlas

Not currently included, but to be incorporated:

*Cancer Data Services* datasets – not included because of difficulties
in collapsing the disease types categories in those datasets with the
disease types present in these larger sources. Generally issues of much
higher specificity in the CDS datasets; for example, one disease type
from MCI is “Pilocytic Astrocytoma Kiaa1549::Braf Fusion-Positive (Who
Grade 1)”.

*Genomic Data Commons* – not included because the GENIE dataset of over
40,000 cases so badly swamps everything else; the graph ends up being a
GENIE graph rather than showing all of the other datasets.

*TCIA* - The Cancer Imaging Archive – not included because there is no
field in the TCIA data analogous to “disease type”.

#### Import metadata files

These metadata files were created by using Data Browser.

``` r
filesdir = "data/"


# This parse tables function doesn't do much besides make a dataframe and add
# the dataset name

parseTables <- function(tsv, name) {
   df <- read.table(tsv, header=TRUE, sep = "\t")
   df <- df %>% mutate(dataset = name)
   return(df)
}

# Read in the datasets' metadata 

cptac.file <- paste(filesdir, "cptac_20230815_0.tsv", sep="")
cptac3.file <- paste(filesdir, "cptac_3_20230815_0.tsv", sep="")
tcga.file <- paste(filesdir, "tcga_grch38_20230815_0.tsv", sep="")
tcia.file <- paste(filesdir, "tcia_20230815_0.tsv", sep="")
target.file <- paste(filesdir, "target_20230815_0.tsv", sep="")
fm.file <- paste(filesdir, "fm_20230815_0.tsv", sep="")
icgc.file <- paste(filesdir, "icgc_20230815_0.tsv", sep="")


cptac.df <- parseTables(cptac.file, "cptac")
cptac3.df <- parseTables(cptac3.file, "cptac3")
tcga.df <- parseTables(tcga.file, "tgca")
fm.df <- parseTables(fm.file, "fm")
tcia.df <- parseTables(tcia.file, "tcia")
target.df <- parseTables(target.file, "target")
icgc.df <- parseTables(icgc.file, "icgc")
```

Now we’ve got these datasets, we’ve got to harmonize their column names.

``` r
colnames(cptac.df)
```

    ## [1] "case"               "case_diseasetype_1" "case_primarysite_1"
    ## [4] "investigation"      "dataset"

``` r
colnames(tcia.df)
```

    ## [1] "case"               "case_diseasetype_1" "investigation"     
    ## [4] "dataset"

``` r
colnames(icgc.df)
```

    ## [1] "donor"                 "donor_primarysite_1"   "project"              
    ## [4] "project_projectname_1" "dataset"

``` r
# let's just keep the caseid/donor, primary site, and disease type. And rename them as we go. 

cptac.df <- cptac.df %>% select(case, disease_type = case_diseasetype_1, primary_site = case_primarysite_1, dataset)
cptac3.df <- cptac3.df %>% select(case, disease_type = case_diseasetype_1, primary_site = case_primarysite_1, dataset)
fm.df <- fm.df %>% select(case, disease_type=case_diseasetype_1, primary_site=case_primarysite, dataset)
icgc.df <- icgc.df %>% mutate(disease_type=NA) %>% select(case=donor, disease_type, primary_site=donor_primarysite_1, dataset)
target.df <- target.df %>% select(case, disease_type = case_diseasetype_1, primary_site=case_primarysite_1, dataset)
tcga.df <- tcga.df %>% select(case, disease_type=case_diseasetype_1, primary_site=case_primarysite_1, dataset)
tcia.df <- tcia.df %>% mutate(primary_site=NA) %>% select(case, disease_type=case_diseasetype_1, primary_site, dataset)
```

Now they all have the same column names in the same order. We can row
bind them, which we’ll do after importing PDC and ICDC.

#### Add Proteomic Data Commons

Metadata from PDC and ICDC were downloaded from their respective data
repositories. Both of them permit download of metadata as a manifest; in
the case of PDC, the “biospecimens” table was downloaded from
<https://proteomic.datacommons.cancer.gov/pdc/browse>.

``` r
pdc.file <- paste(filesdir, "PDC_metadata_082223.txt", sep="")
pdc.df <- parseTables(pdc.file, "pdc")
colnames(pdc.df)
```

    ##  [1] "Aliquot.ID"                         "Aliquot.Submitter.ID"              
    ##  [3] "Sample.ID"                          "Sample.Submitter.ID"               
    ##  [5] "Case.ID"                            "Case.Submitter.ID"                 
    ##  [7] "Project.Name"                       "Sample.Type"                       
    ##  [9] "Primary.Site"                       "Disease.Type"                      
    ## [11] "Aliquot.Is.Ref"                     "Aliquot.Status"                    
    ## [13] "Aliquot.Quantity"                   "Aliquot.Volume"                    
    ## [15] "Amount"                             "Analyte.Type"                      
    ## [17] "Concentration"                      "Case.Status"                       
    ## [19] "Sample.Status"                      "Sample.Is.Ref"                     
    ## [21] "Biospecimen.Anatomic.Site"          "Composition"                       
    ## [23] "Current.Weight"                     "Days.To.Collection"                
    ## [25] "Days.To.Sample.Procurement"         "Diagnosis.Pathologically.Confirmed"
    ## [27] "Freezing.Method"                    "Initial.Weight"                    
    ## [29] "Intermediate.Dimension"             "Longest.Dimension"                 
    ## [31] "Method.Of.Sample.Procurement"       "Pathology.Report.UUID"             
    ## [33] "Preservation.Method"                "Sample.Type.id"                    
    ## [35] "Sample.Ordinal"                     "Shortest.Dimension"                
    ## [37] "Time.Between.Clamping.And.Freezing" "Time.Between.Excision.and.Freezing"
    ## [39] "Tissue.Collection.Type"             "Tissue.Type"                       
    ## [41] "Tumor.Code"                         "Tumor.Code.ID"                     
    ## [43] "Tumor.Descriptor"                   "Program.Name"                      
    ## [45] "dataset"

``` r
pdc.df <- pdc.df %>%
  select(case=Case.ID, disease_type=Disease.Type, primary_site=Primary.Site, dataset)
```

#### Canine data commons

In the case of the ICDC, the “Cases” table was downloaded from
<https://caninecommons.cancer.gov/#/explore> - Select the “Cases” tab,
look for the icon of a cloud with an arrow in it right above the right
hand side of the table, and click it to “Download table contents as
csv”.

Extensive manual data sanitizing was done for ICDC, due to
capitalization issues and comma insertions.

``` r
icdc.file <- paste(filesdir, "ICDC_Cases_082223_metadata_sanitized.txt", sep="")
icdc.df <- parseTables(icdc.file, "icdc")

icdc.df <- icdc.df %>% select(case=Case.ID, disease_type=Diagnosis, primary_site=Disease.Site, dataset) %>%
  unique()

## there are 678 unique cases, which matches what is on the website. 
```

#### Bind the small dataframes into a large one

``` r
datasets.df <- bind_rows(cptac.df, cptac3.df, fm.df, icgc.df, target.df, tcga.df, tcia.df, pdc.df, icdc.df)
dim(datasets.df)
```

    ## [1] 19948     4

#### Count the disease types and determine the top twenty

We want two things – we want to count the disease types in each dataset,
and we also want the top disease types over all.

``` r
disease_types <- datasets.df %>% group_by(dataset) %>% count(disease_type) %>% arrange(desc(n)) %>% ungroup()
# Get rid of OTHER and NOS
disease_types <- disease_types %>% filter(!(disease_type == "Other"))
disease_types$disease_type <- gsub(", NOS", "", disease_types$disease_type, ignore.case = TRUE)
disease_types$disease_type <- gsub("Pediatric/AYA Brain Tumors", "Pediatric Brain Tumors", disease_types$disease_type)
write.table(disease_types, "../out/datasets_disease_types.txt", sep="\t", row.names = FALSE, quote = FALSE, eol="\n")

datasets_disease <- disease_types
top_diseases <- datasets_disease %>% group_by(disease_type) %>% summarize(n=sum(n)) %>% arrange(desc(n))
top_25 <- top_diseases %>% filter(!(is.na(disease_type))) %>% slice_head(n=25)

knitr::kable(top_25, caption="The top 25 disease types in our hosted data regardless of source", format = 'markdown')
```

| disease_type                                |    n |
|:--------------------------------------------|-----:|
| Adenomas and Adenocarcinomas                | 2509 |
| Ductal and Lobular Neoplasms                | 1526 |
| Myeloid Leukemias                           | 1369 |
| Lymphoid Leukemias                          | 1185 |
| Ovarian Serous Cystadenocarcinoma           |  946 |
| Squamous Cell Neoplasms                     |  789 |
| Lung Adenocarcinoma                         |  592 |
| Epithelial Neoplasms                        |  540 |
| Clear Cell Renal Cell Carcinoma             |  483 |
| Transitional Cell Papillomas and Carcinomas |  480 |
| Breast Invasive Carcinoma                   |  479 |
| Colon Adenocarcinoma                        |  446 |
| Pancreatic Ductal Adenocarcinoma            |  421 |
| Head and Neck Squamous Cell Carcinoma       |  373 |
| Uterine Corpus Endometrial Carcinoma        |  333 |
| Hepatocellular Carcinoma                    |  330 |
| Neuroepitheliomatous Neoplasms              |  306 |
| Osteosarcoma                                |  281 |
| Kidney Renal Clear Cell Carcinoma           |  267 |
| Lung Squamous Cell Carcinoma                |  265 |
| Pediatric Brain Tumors                      |  264 |
| Glioblastoma Multiforme                     |  262 |
| Gliomas                                     |  258 |
| Cystic, Mucinous and Serous Neoplasms       |  222 |
| Low Grade Glioma                            |  199 |

The top 25 disease types in our hosted data regardless of source

This is a nice table, but we want it in a plot.

### Part 2 - Create the dataframe for plotting and make the chart

``` r
# Let's collect these categories of disease_types from the dataset-grouped list, 
# but add the Level information to it from top_25$disease_type so we 
# can use that for making Factors. 
disease_plot_data <- datasets_disease %>% filter(disease_type %in% top_25$disease_type) %>%
   mutate(level = match(disease_type, top_25$disease_type))

 disease_plot_data <- disease_plot_data %>% arrange(by=level)
 
 disease_plot_data$disease_type <- factor(disease_plot_data$disease_type, unique(disease_plot_data$disease_type))
 levels(disease_plot_data$disease_type)
```

    ##  [1] "Adenomas and Adenocarcinomas"               
    ##  [2] "Ductal and Lobular Neoplasms"               
    ##  [3] "Myeloid Leukemias"                          
    ##  [4] "Lymphoid Leukemias"                         
    ##  [5] "Ovarian Serous Cystadenocarcinoma"          
    ##  [6] "Squamous Cell Neoplasms"                    
    ##  [7] "Lung Adenocarcinoma"                        
    ##  [8] "Epithelial Neoplasms"                       
    ##  [9] "Clear Cell Renal Cell Carcinoma"            
    ## [10] "Transitional Cell Papillomas and Carcinomas"
    ## [11] "Breast Invasive Carcinoma"                  
    ## [12] "Colon Adenocarcinoma"                       
    ## [13] "Pancreatic Ductal Adenocarcinoma"           
    ## [14] "Head and Neck Squamous Cell Carcinoma"      
    ## [15] "Uterine Corpus Endometrial Carcinoma"       
    ## [16] "Hepatocellular Carcinoma "                  
    ## [17] "Neuroepitheliomatous Neoplasms"             
    ## [18] "Osteosarcoma"                               
    ## [19] "Kidney Renal Clear Cell Carcinoma"          
    ## [20] "Lung Squamous Cell Carcinoma"               
    ## [21] "Pediatric Brain Tumors"                     
    ## [22] "Glioblastoma Multiforme"                    
    ## [23] "Gliomas"                                    
    ## [24] "Cystic, Mucinous and Serous Neoplasms"      
    ## [25] "Low Grade Glioma"

``` r
disease_plot <- ggplot(disease_plot_data, aes(x=disease_type, y=n, group=dataset)) +
  theme_minimal() +
  geom_col(position="stack", aes(fill=dataset)) + 
  coord_flip() + 
  scale_x_discrete(limits=rev(levels(disease_plot_data$disease_type))) +
  theme(panel.grid = element_blank(),
        axis.text.y=element_text(size=8)) +
  labs(x="", y="Number of cases", title="Top 25 Disease Types in Connected Datasets")

palette_cgc8 <- c("#1da1f2", "#1b94df", "#1887cb", "#167ab8", "#146da5", "#116191", "#0f547e", "#0d476a")

disease_plot +
  scale_color_manual(values=rev(palette_cgc8), aesthetics = c("color", "fill"))
```

![](DatasetsVisualizations_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
# png(filename = "../out/Top25DiseaseTypesInConnectedData.png", width = 720)
# disease_plot +
#   scale_color_manual(values=rev(palette_cgc8), aesthetics = c("color", "fill"))
# dev.off()
# 
# pdf(file="../out/Top25DiseaseTypesInConnectedData.pdf")
# disease_plot +
#   scale_color_manual(values=rev(palette_cgc8), aesthetics = c("color", "fill"))
# dev.off()
# 
# disease_plot +
#   scale_color_manual(values=rev(palette_cgc9), aesthetics = c("color", "fill"))
# dev.off()
```

The below block was used to make higher resolution version with bigger
fonts and capitalize the datasets – this is the website version. I don’t
evaluate this code in this notebook because it will look terrible in a
markdown doc.

``` r
disease_plot_data <- disease_plot_data %>% mutate(dataset = toupper(dataset))

disease_plot <- ggplot(disease_plot_data, aes(x=disease_type, y=n, group=dataset)) +
  theme_minimal() +
  geom_col(position="stack", aes(fill=dataset)) + 
  coord_flip() + 
  scale_x_discrete(limits=rev(levels(disease_plot_data$disease_type))) +
  theme(panel.grid = element_blank(),
        axis.text.y=element_text(size=28),
        plot.title=element_text(size=34, hjust = 0),
        legend.key.size=unit(1.5, 'cm'),
        legend.text=element_text(size=22),
        legend.title=element_blank(),
        axis.text.x=element_text(size=24),
        axis.title.x=element_text(size=24)) +
  labs(x="", y="Number of cases", title="Top disease types in the connected CRDC datasets")

disease_plot +
  scale_color_manual(values=rev(palette_cgc9), aesthetics = c("color", "fill"))
# 
# png(filename = "out/Top25DiseaseTypesInConnectedData_ForWeb.png", width = 2160, height = 1440 )
# disease_plot +
#   scale_color_manual(values=rev(palette_cgc8), aesthetics = c("color", "fill"))
# dev.off()

# svg(filename="out/Top25DiseaseTypes_forWeb.svg", width=720)
# disease_plot_svg +
#   scale_color_manual(values=rev(palette_cgc8), aesthetics = c("color", "fill"))
# dev.off()
```

# Primary Site

We start from our same datasets.df, but to make things a little more
manageable, I am doing some synonymizing over these anatomical sites.

``` r
datasets.df$primary_site <- gsub("Liver and intrahepatic bile ducts", "Liver", datasets.df$primary_site)
# and just let's remove the "NOS" after Uterus
datasets.df$primary_site <- gsub(", NOS", "", datasets.df$primary_site)

datasets.df$primary_site <- gsub("Hematopoietic and reticuloendothelial systems", "Bone marrow", datasets.df$primary_site)
datasets.df$primary_site <- gsub("Spinal cord, cranial nerves, and other parts of central nervous system", "Central nervous system", datasets.df$primary_site)
datasets.df$primary_site <- gsub("Bronchus and lung", "Lung", datasets.df$primary_site)
datasets.df$primary_site <- gsub("Cervix uteri", "Cervix and uterus", datasets.df$primary_site)
datasets.df$primary_site <- gsub("Uterus", "Cervix and uterus", datasets.df$primary_site)
datasets.df$primary_site <- gsub("Head and Neck", "Head and neck", datasets.df$primary_site)
datasets.df$primary_site <- gsub("Rectum", "Colorectal", datasets.df$primary_site)
datasets.df$primary_site <- gsub("Peripheral nerves and autonomic nervous system", "Autonomic nervous system", datasets.df$primary_site)
datasets.df$primary_site <- gsub("Heart, mediastinum, and pleura", "Heart", datasets.df$primary_site)




top_primary_sites <- datasets.df %>%
  count(primary_site) %>%
  arrange(desc(n)) %>%
  filter(!is.na(primary_site))



top_25_ps <- datasets.df %>%
  count(primary_site) %>%
  arrange(desc(n)) %>%
  filter(!is.na(primary_site)) %>%
  filter(!(primary_site == "Other and ill-defined sites"),
         !(primary_site == "Unknown"), 
         !(primary_site == "Not Reported")) %>%
  slice_head(n=25)
  


knitr::kable(top_25_ps, caption="Primary sites ranked by number of cases, regardless of source")
```

| primary_site                   |    n |
|:-------------------------------|-----:|
| Bone marrow                    | 2800 |
| Breast                         | 2049 |
| Lung                           | 1727 |
| Kidney                         | 1230 |
| Pancreas                       | 1227 |
| Colon                          | 1107 |
| Ovary                          | 1049 |
| Liver                          |  783 |
| Cervix and uterus              |  761 |
| Esophagus                      |  665 |
| Bladder                        |  570 |
| Brain                          |  565 |
| Head and neck                  |  353 |
| Bone                           |  332 |
| Stomach                        |  226 |
| Adrenal gland                  |  203 |
| Blood                          |  203 |
| Colorectal                     |  178 |
| Central nervous system         |  175 |
| Prostate                       |  149 |
| Skin                           |  128 |
| Lymph Node                     |  116 |
| Heart                          |   58 |
| Retroperitoneum and peritoneum |   52 |
| Autonomic nervous system       |   39 |

Primary sites ranked by number of cases, regardless of source

Now we’ll break them out by dataset:

``` r
primary_sites_datasets <- datasets.df %>% 
  filter(!(is.na(primary_site)),
         !(primary_site == "Other and ill-defined sites"),
         !(primary_site == "Unknown"), 
         !(primary_site == "Not reported")) %>%
  group_by(dataset) %>% 
  count(primary_site) %>% 
  arrange(desc(n)) 
```

Alright. Now we just have to repeat what we did to get the top 25 and
factor them correctly

``` r
prsite_plot_data <- primary_sites_datasets %>% filter(primary_site %in% top_25_ps$primary_site) %>%
  mutate(level=match(primary_site, top_25_ps$primary_site))


prsite_plot_data <- prsite_plot_data %>% arrange(by=level)
prsite_plot_data$primary_site <- factor(prsite_plot_data$primary_site, unique(prsite_plot_data$primary_site))

levels(prsite_plot_data$primary_site)
```

    ##  [1] "Bone marrow"                    "Breast"                        
    ##  [3] "Lung"                           "Kidney"                        
    ##  [5] "Pancreas"                       "Colon"                         
    ##  [7] "Ovary"                          "Liver"                         
    ##  [9] "Cervix and uterus"              "Esophagus"                     
    ## [11] "Bladder"                        "Brain"                         
    ## [13] "Head and neck"                  "Bone"                          
    ## [15] "Stomach"                        "Adrenal gland"                 
    ## [17] "Blood"                          "Colorectal"                    
    ## [19] "Central nervous system"         "Prostate"                      
    ## [21] "Skin"                           "Lymph Node"                    
    ## [23] "Heart"                          "Retroperitoneum and peritoneum"
    ## [25] "Autonomic nervous system"

And make the plot with the colors:

``` r
site_plot <- ggplot(prsite_plot_data, aes(x=primary_site, y=n, group=dataset)) +
  theme_minimal() +
  geom_col(position="stack", aes(fill=dataset)) + 
  coord_flip() + 
  scale_x_discrete(limits=rev(levels(prsite_plot_data$primary_site))) +
  scale_y_continuous(breaks=c(0, 500, 1000, 1500, 2000, 3000)) +
  # scale_x_discrete(limits=rev(levels(disease_plot_data$disease_type)), labels = function(x) str_wrap(x, width=40)) +
  theme(panel.grid = element_blank(),
        axis.text.y=element_text(size=8)) +
  labs(x="", y="Number of cases", title="Top 25 Sampled Body Sites in Connected Datasets")
  


site_plot + 
  scale_color_manual(values=rev(palette_cgc8), aesthetics = c("color", "fill"))
```

![](DatasetsVisualizations_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->
\##### Save Rdata and print session info

``` r
save.image(file="SummarizeConnectedDatasets.RData")
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 22621)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.2 forcats_1.0.0   stringr_1.5.0   dplyr_1.1.0    
    ##  [5] purrr_1.0.1     readr_2.1.4     tidyr_1.3.0     tibble_3.2.0   
    ##  [9] ggplot2_3.4.2   tidyverse_2.0.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] highr_0.10       pillar_1.9.0     compiler_4.2.2   tools_4.2.2     
    ##  [5] digest_0.6.31    timechange_0.2.0 evaluate_0.21    lifecycle_1.0.3 
    ##  [9] gtable_0.3.3     pkgconfig_2.0.3  rlang_1.1.0      cli_3.6.0       
    ## [13] rstudioapi_0.14  yaml_2.3.7       xfun_0.37        fastmap_1.1.1   
    ## [17] withr_2.5.0      knitr_1.42       generics_0.1.3   vctrs_0.6.0     
    ## [21] hms_1.1.3        grid_4.2.2       tidyselect_1.2.0 glue_1.6.2      
    ## [25] R6_2.5.1         fansi_1.0.4      rmarkdown_2.21   farver_2.1.1    
    ## [29] tzdb_0.3.0       magrittr_2.0.3   scales_1.2.1     htmltools_0.5.4 
    ## [33] colorspace_2.1-0 labeling_0.4.2   utf8_1.2.3       stringi_1.7.12  
    ## [37] munsell_0.5.0
