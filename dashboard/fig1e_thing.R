#install.packages("shinydashboard")

library(shiny)
library(shinydashboard)
library(tidyverse)
library(argparse)
library(biomaRt)

ensembl_id_conversion <- function(ensembl_id){

  mart <- readRDS("./mart.Rds")
#   print("********")
#   print(ensembl_id)
  converted <- getBM(values=ensembl_id,
    filters= "ensembl_gene_id", 
    attributes= c("external_gene_name", "description"),
    mart= mart)
  # print(converted)
  # print(class(converted))
  # print("********")
  return(converted)
}

parser <- ArgumentParser()
parser$add_argument("--donor_list_csv")
parser$add_argument("--port")
args <- parser$parse_args()
print("args parsed.")
print(paste("port:", args$port))
options(shiny.port = as.numeric(args$port))

donor_df <- read.csv(args$donor_list_csv)
print("donor_df")
print(donor_df)

constructed_df <- data.frame(
  sample_name = character(), 
  pre_UMI_num = numeric(),
  post_UMI_num = numeric(),
  false_neg_num = numeric()
)

all_umi_pass <- list()

print('right before loop')
for (i in 1:nrow(donor_df)) {
  pipeline_dir <- donor_df[i,4]
  sample_name <- donor_df[i, 1]
  # print(pipeline_dir)
  print(sample_name)

  step2_dir <- file.path(pipeline_dir, "step2_out", "*", "*.bai")
  num_cells <- Sys.glob( step2_dir) 
  num_cells <- length(num_cells)

  # step9_path <- file.path(pipeline_dir, "step9", "filtered_ScoredMutations.csv")
  step9_path <- file.path(pipeline_dir, "step9", "ScoredMutations.csv")
  step9_csv <- read.csv(step9_path, header=T)
  print("step9_csv")
  print(dim(step9_csv))
  print(table(step9_csv$mutect_filter))
  mutect_pass <- step9_csv %>%
    filter(mutect_filter == 'pass')
  pre_UMI_num <- nrow(mutect_pass)
  print("pre_UMI_num")
  print(pre_UMI_num)

  # umi_pass <- step9_csv %>%
  #   filter(( reads_in_umi_filter == 'pass') | recovered_double==T) %>%
  #   dplyr::select(ENSEMBL_GENE_ID, Chr, POS, REF, ALT, AA_CHANGE, bc) 
  umi_pass <- step9_csv %>%
    filter(( reads_in_umi_filter == 'pass') ) %>%
    dplyr::select(ENSEMBL_GENE_ID, Chr, POS, REF, ALT, AA_CHANGE, bc) 
  # umi_pass <- step9_csv %>%
  #   filter((reads_in_umi_filter == 'pass' & umi_fraction_filter == 'pass') | recovered_double==T) %>%
  #   dplyr::select(ENSEMBL_GENE_ID, Chr, POS, REF, ALT, AA_CHANGE, bc) 
  post_UMI_num <- nrow(umi_pass)

  print("umi_pass")
  print(nrow(umi_pass))
  print(unique(step9_csv$reads_in_umi_filter))
  print(unique(step9_csv$umi_fraction_filter))
  print(table(step9_csv$reads_in_umi_filter, step9_csv$umi_fraction_filter))
  if (nrow(umi_pass)==0){
    next
  }

  all_umi_pass[[i]] <- umi_pass

  false_neg <- step9_csv %>%
    filter(recovered_double==T)  
  false_neg_num <- nrow(false_neg)

  constructed_df <- rbind(constructed_df, 
    data.frame (sample_name = sample_name,
      pre_UMI_num = pre_UMI_num,
      post_UMI_num = post_UMI_num,
      false_neg_num = false_neg_num,
      num_cells = num_cells)
    )
}
# print(donor_df)
# print(constructed_df)
combined_df <- inner_join(donor_df, constructed_df, by="sample_name")
#print("here")
#print(all_umi_pass)
all_umi_pass <- bind_rows(all_umi_pass)
print("all_umi_pass")
print(dim(all_umi_pass))
# print("all_umi_pass")
# print(all_umi_pass)

top_genes <- all_umi_pass %>% 
  count(ENSEMBL_GENE_ID, Chr, POS, REF, ALT,AA_CHANGE, sort = TRUE) %>%
  filter(!grepl( "-", ENSEMBL_GENE_ID)) %>%
  filter(!is.na(AA_CHANGE)) %>%  
  head(n=10L)
print(top_genes)

top_gene_summary <- bind_rows(lapply(top_genes$ENSEMBL_GENE_ID, ensembl_id_conversion))
top_gene_summary <- top_gene_summary %>%
  mutate(num_cells = top_genes$n,
    Chromosome = top_genes$Chr,
    Position = top_genes$POS,
    NT_mutated_from = top_genes$REF,
    NT_mutated_to = top_genes$ALT,
    AA_CHANGE = top_genes$AA_CHANGE) %>%
  rename(gene_description = description,
    gene_id = external_gene_name)

fig1e_table <- combined_df %>%
  dplyr::select(sample_name, pre_UMI_num, post_UMI_num)
write.csv(fig1e_table, "fig1e_table.csv", row.names=F)

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
        cols = ends_with("UMI_num"),
        names_to = "correction_state",
        values_to = "num_mutation"
        ) %>%
      mutate(correction_state = factor(correction_state, 
        levels=c("pre_UMI_num", "post_UMI_num")))
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
