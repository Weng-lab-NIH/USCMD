#install.packages("shinydashboard")

library(shiny)
library(shinydashboard)
library(tidyverse)
library(argparse)
library(biomaRt)
library(MASS)
library(lme4)

ensembl_id_conversion <- function(ensembl_id){

  mart <- readRDS("./mart.Rds")
#   print("********")
#   print(ensembl_id)
  converted <- getBM(values=ensembl_id,
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id","external_gene_name", "description"),
    mart= mart)
  if(length(converted == 0)){
    return(c("ensembl_gene_id" = NA,"external_gene_name" = NA, 
      "description" = NA))
  }
  # print(converted)
  # print(class(converted))
  # print("********")
  return(converted)
}

parser <- ArgumentParser()
parser$add_argument("--donor_list_csv")
parser$add_argument("--port")
args <- parser$parse_args()
options(shiny.port = as.numeric(args$port))

donor_df <- read.csv(args$donor_list_csv)

constructed_df <- data.frame(
  sample_name = character(), 
  pre_UMI_num = numeric(),
  post_UMI_num = numeric(),
  false_neg_num = numeric()
)

all_umi_pass <- list()

data_per_cell <- data.frame(
  sample = character(),
  log2_mut = numeric(),
  UMI = numeric(),
  COV = numeric(),
  COV_EXOME = numeric()
)

for (i in 1:nrow(donor_df)){
  pipeline_dir <- donor_df[i,4]
  sample_name <- donor_df[i, 1]
  print(sample_name)

  step2_dir <- file.path(pipeline_dir, "step2_out", "*", "*.bai")
  num_cells <- Sys.glob( step2_dir) 
  num_cells <- length(num_cells)

  step11_path <- file.path(pipeline_dir, "step11_out", "mutations.csv")
  step11_csv <- read.csv(step11_path, header=F)
  colnames(step11_csv) <- c(
    "Sample",
    "READ_EXOME",
    "READ_SAM",
    "UMI",
    "COV_EXOME",
    "COV")
  pattern_to_remove <- paste0("^.*/step2_out/+.*/")
  print(pattern_to_remove)
  step11_csv <- mutate_at(step11_csv, "Sample", str_replace, pattern_to_remove, "") %>%
    mutate_at("Sample", str_replace, ".bam", "") %>%
    mutate_at("Sample", str_sub, -18)
  print(head(step11_csv))
  umi_counts <- step11_csv$V4
  umi_count <- sum(umi_counts)

  step9_path <- file.path(pipeline_dir, "step9_out", "filtered_ScoredMutations.csv")
  step9_csv <- read.csv(step9_path, header=T)
  # print("step9_csv")
  # print(head(step9_csv))
  mutect_pass <- step9_csv %>%
    filter(mutect_filter == 'pass')
  pre_UMI_num <- nrow(mutect_pass)

  umi_pass <- step9_csv %>%
    filter(
      (reads_in_umi == 'pass' & 
      umi_fraction_filter == 'pass' & 
      num_variant_filter == 'pass') | 
    recovered_double==T) %>%
    dplyr::select(ENSEMBL_GENE_ID, Chr, POS, REF, ALT, AA_CHANGE, bc) 
  post_UMI_num <- nrow(umi_pass)

  # print("umi_pass")
  # print(nrow(umi_pass))
  if (nrow(umi_pass)==0){
    next
  }

  all_umi_pass[[i]] <- umi_pass

  false_neg <- step9_csv %>%
    filter(recovered_double==T)  
  false_neg_num <- nrow(false_neg)

  # false negatives are just all the ones that has.two.variant
  #normalized counts
  num_mut_per_cell <- umi_pass %>% 
    group_by(sample, bc) %>% tally()
  # print("num_mut_per_cell")
  # print(num_mut_per_cell)
  num_mut_per_cell$log2_mut <- log2(num_mut_per_cell$n + 1)
  new_data_per_cell <- full_join(num_mut_per_cell, step11_csv, 
    by = c('bc', 'sample')) %>%
    replace_na(list("n"=1, "log2_mut"=0))
  new_data_per_cell$sample <- sample_name
  
  data_per_cell <- rbind(data_per_cell, new_data_per_cell)


  constructed_df <- rbind(constructed_df, 
    data.frame (sample_name = sample_name,
      pre_UMI_num = pre_UMI_num,
      post_UMI_num = post_UMI_num,
      false_neg_num = false_neg_num,
      num_cells = num_cells,
      umi_count = umi_count)
    )
}
# print(sapply(donor_df, class))
# print(sapply(constructed_df, class))
combined_df <- inner_join(donor_df, constructed_df, by="sample_name")
all_umi_pass <- bind_rows(all_umi_pass)
num_umi_pass <- dim(all_umi_pass)[1]

top_genes <- all_umi_pass %>% 
  count(ENSEMBL_GENE_ID, Chr, POS, REF, ALT,AA_CHANGE, sort = TRUE) %>%
  filter(!grepl( "-", ENSEMBL_GENE_ID)) %>%
  #filter(!is.na(AA_CHANGE)) %>%  
  top_n(min(10, num_umi_pass))

top_gene_summary <- bind_rows(lapply(top_genes$ENSEMBL_GENE_ID, ensembl_id_conversion))
top_gene_summary <- top_gene_summary %>%
  full_join(top_genes, by=c("ensembl_gene_id" = "ENSEMBL_GENE_ID"))

print("data_per_cell")
print(data_per_cell)

nbGLM <- glm.nb(log2_mut ~ UMI + COV + COV_EXOME, 
  data=data_per_cell)
data_per_cell$norm_mut = (nbGLM$residuals + nbGLM[["coefficients"]][["(Intercept)"]])

adj_quantile <- quantile(data_per_cell$norm_mut, probs = c(0.001,0.999))
data_per_cell[data_per_cell$norm_mut < min(adj_quantile),]$norm_mut <- min(adj_quantile)
data_per_cell[data_per_cell$norm_mut > max(adj_quantile),]$norm_mut <- max(adj_quantile)
value <- data_per_cell$norm_mut
data_per_cell$scale_mut <- (value - min(value)) / (max(value) - min(value)) * max(data_per_cell$n) 

scaled_mut_num <- group_by(data_per_cell, sample) %>%
  tally() %>% rename(sample_name = sample)
combined_df$scaled_mut_num <- scaled_mut_num$n

## app.R ##
ui <- dashboardPage(
  dashboardHeader(title = "USCMD Output Summary"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Widgets", tabName = "widgets", icon = icon("th"))
    )
  ),
  ## Body content
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              fluidRow(
                # # A static valueBox
                # valueBox(33, "Number of Samples", icon = icon("align-right")),
                # # A static valueBox
                # valueBox(3450, "Number of Cells / Sample", icon = icon("align-right")),
                # # A static valueBox
                # valueBox(100, "Quality", icon = icon("align-right"))
              ),
              
              fluidRow(
                box(plotOutput("plot1", height = 250), 
                    background = 'black'),
                box(plotOutput("plot2", height = 250), 
                    background = 'black')
                
                # box(
                #   title = "Controls",
                #   background = 'black',
                #   sliderInput("slider", "Max Age of Donor:", 1, 100, 50)
                # )
              ),
              
              fluidRow(
                box(tableOutput('table2'))
                ),
              fluidRow(
                box(tableOutput('table1')),
                
                # box(
                #   title = "Controls",
                #   background = 'black',
                #   sliderInput("slider", "Max Age of Donor:", 1, 100, 50)
                # )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "widgets",
              h2("Widgets tab content")
      )
    )
  )
)


server <- function(input, output) {
  set.seed(122)
  


  output$plot1 <- renderPlot({
    data <- combined_df %>%
      pivot_longer(
        cols = ends_with("_num"),
        names_to = "correction_state",
        values_to = "num_mutation"
        ) %>%
      mutate(correction_state = factor(correction_state, 
        levels=c("pre_UMI_num", "post_UMI_num", "scaled_mut_num"))) %>%
      drop_na()
    #print("**********************, data")
    #print(data)
    ggplot(data, aes(x=sample_name, y=num_mutation, group = correction_state)) +
      theme_minimal() + theme_bw() +
      ylab('Mutation Number (no UMI correction)') + xlab('Sample') +
      geom_bar(stat = 'identity', position = 'dodge2', alpha = 0.5, 
        aes(color = correction_state, fill = correction_state))
  })
  output$plot2 <- renderPlot({
    data <- combined_df
    ggplot(data, aes(x=sample_name, y=false_neg_num)) +
      theme_minimal() + theme_bw() +
      ylab('Number of False Negatives') + xlab('Sample') +
      geom_bar(stat = 'identity', position = 'dodge2', alpha = 0.5) +
      scale_fill_manual(values = c('grey30', 'blue', 'skyblue'))
  })
  output$table1 <- renderTable({
    top_gene_summary
    })
  output$table2 <- renderTable({
    data <- combined_df
    })
}
shinyApp(ui, server)
