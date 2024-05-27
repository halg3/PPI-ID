library(shiny)
library(DT)
library(readr)
library(dplyr)

ppidm = read.csv("/Users/hvygoodwin/Downloads/DifB_Project/app/PPIDM_FullSortedDataset_84K_GSB.csv")
compiled_interactions = read.csv("/Users/hvygoodwin/Downloads/DifB_Project/app/ddi_interactions.csv")
compiled_interactions = na.omit(compiled_interactions)
# write.csv(compiled_interactions, "/Users/hvygoodwin/Downloads/DifB_Project/app/ddi_interactions.csv", row.names = FALSE)

interaction_id <- read_tsv("/Users/hvygoodwin/Downloads/DifB_Project/app/domine-tables-2.0/elm_interaction_domains.tsv")
colnames(interaction_id) <- c("elm_identifier", "interaction_domain_id", "interaction_domain_description", "interaction_domain_name")

dif_interpro <- read_tsv("/Users/hvygoodwin/Downloads/dif_interpro_sample.tsv")
dif_elm <- read_tsv("/Users/hvygoodwin/Downloads/dif_elm_sample - Sheet1.tsv")
cact_interpro <- read_tsv("/Users/hvygoodwin/Downloads/cact_interpro_sample.tsv")
cact_elm <- read_tsv("/Users/hvygoodwin/Downloads/cact_elm_sample - Sheet2.tsv")
table_titles <- read_tsv("/Users/hvygoodwin/Downloads/table_titles - Sheet2.tsv")

slim_pfam <- function(protein){
  pfam_id = c()
  
  # Only keeps SLiM names that are also found in interaction_id. interaction_id is from the ELM database, and it 
  # provides identifying information for SLiMs.
  keep = intersect(protein$slim, interaction_id$elm_identifier)
  protein_filtered = protein[protein$slim %in% keep, ]
  
  slim = protein_filtered$slim
  range = protein_filtered$range
  prob = protein_filtered$probability
  
  for(motif in protein_filtered$slim){
    for(i in 1:nrow(interaction_id)){
      if(motif == interaction_id$elm_identifier[i]){
        pfam_id = c(pfam_id, interaction_id$interaction_domain_id[i])
        break
      }
    }
  }
  return(data.frame(slim, ID = pfam_id, range, prob))
}
 
elm_algo <- function(protein1, protein2){
  
  table_names = c("slim", "sequence", "range", "misc", "description", "location", "pattern", "phi", "structure", "probability")
  default_table_names = c("Elm Name", "Instances (Matched Sequence)", "Positions", "View in Jmol", "Elm Description", "Cell Compartment", "Pattern", "PHI-Blast Instance Mapping", "Structural Filter Info", "Probability")
  
  if(ncol(protein1) != 10 | ncol(protein2) != 10){
    return("Please ensure that ELM results for both proteins are submitted")
  } else {
    if(sum(colnames(protein1) == default_table_names) == 10){
      colnames(protein1) = table_names
    } else if(colnames(protein1)[1] %in% interaction_id$elm_identifier){
      firstrow1 = colnames(protein1)
      protein1 = rbind(firstrow1, protein1)
      colnames(protein1) <- table_names
    } else {
      return("Please ensure that ELM data for Protein 1 is formatted correctly")
    }
    if(sum(colnames(protein2) == default_table_names) == 10){
      colnames(protein2) = table_names
    } else if(colnames(protein2)[1] %in% interaction_id$elm_identifier){
      firstrow2 = colnames(protein2)
      protein2 = rbind(firstrow2, protein2)
      colnames(protein2) <- table_names
    } else {
      return("Please ensure that ELM data for Protein 2 is formatted correctly")
    }
    
    protein1_range = gsub(" \\[A\\]", "", protein1$range)
    protein1$range = NULL
    protein1$range = protein1_range
    protein2_range = gsub(" \\[A\\]", "", protein2$range)
    protein2$range = NULL
    protein2$range = protein2_range
    
    protein1_cleaned = slim_pfam(protein1)
    protein2_cleaned = slim_pfam(protein2)
    
    domain_name_1 = c()
    domain_name_2 = c()
    domain_range_1 = c()
    domain_range_2 = c()
    combined_prob = c()
    
    for(domain1 in 1:nrow(protein1_cleaned)){
      for(domain2 in 1:nrow(protein2_cleaned)){
        if(protein1_cleaned$ID[domain1] %in% compiled_interactions$domain_1){
          poi1 = compiled_interactions[compiled_interactions$domain_1 == protein1_cleaned$ID[domain1], ]
          if(protein2_cleaned$ID[domain2] %in% poi1$domain_2){
            prob_1 = c()
            prob_2 = c()
            
            domain_name_1 = c(domain_name_1, protein1_cleaned$slim[domain1])
            domain_range_1 = c(domain_range_1, protein1_cleaned$range[domain1])
            prob_1 = c(prob_1, as.numeric(protein1_cleaned$prob[domain1]))
            
            domain_name_2 = c(domain_name_2, protein2_cleaned$slim[domain2])
            domain_range_2 = c(domain_range_2, protein2_cleaned$range[domain2])
            prob_2 = c(prob_2, as.numeric(protein2_cleaned$prob[domain2]))
            
            combined_prob = c(combined_prob, prob_1 * prob_2)
          }
        }
        if(protein1_cleaned$ID[domain1] %in% compiled_interactions$domain_2){
          poi1 = compiled_interactions[compiled_interactions$domain_2 == protein1_cleaned$ID[domain1], ]
          if(protein2_cleaned$ID[domain2] %in% poi1$domain_1){
            prob_1 = c()
            prob_2 = c()
            
            domain_name_1 = c(domain_name_1, protein1_cleaned$slim[domain1])
            domain_range_1 = c(domain_range_1, protein1_cleaned$range[domain1])
            prob_1 = c(prob_1, as.numeric(protein1_cleaned$prob[domain1]))
            
            domain_name_2 = c(domain_name_2, protein2_cleaned$slim[domain2])
            domain_range_2 = c(domain_range_2, protein2_cleaned$range[domain2])
            prob_2 = c(prob_2, as.numeric(protein2_cleaned$prob[domain2]))
            
            combined_prob = c(combined_prob, prob_1 * prob_2)
          }
        }
      }
    }
  }
  unsorted_table = unique(data.frame(domain_name_1, domain_name_2, domain_range_1, domain_range_2, combined_prob))
  
  indices = order(unsorted_table$combined_prob, decreasing = FALSE)
  final_table = unsorted_table[indices,]
  colnames(final_table) <- c("SLiM 1", "SLiM 2", "SLiM 1 Range(s)", "SLiM 2 Range(s)", "Combined Probability")
  row.names(final_table) = NULL
  return(final_table)
}

interpro_algo <- function(protein1, protein2){
  domain_name_1 = c()
  domain_name_2 = c()
  domain_start_1 = c()
  domain_end_1 = c()
  domain_start_2 = c()
  domain_end_2 = c()
  
  if(ncol(protein1) != 15 | ncol(protein2) != 15){
    return("Please ensure that Interpro results for both proteins are submitted.")
  } else {
    firstrow1 = colnames(protein1)
    protein1 = rbind(firstrow1, protein1)
    colnames(protein1) = c("Sequence", "Name", "Length", "Database", "ID", "Description", "Start", "End", "Probability", "T/F", "Date", "Interpro_ID", "DomainName", "GO", "Desc1")
    protein1 = protein1[protein1$Database == "Pfam", ]
  
    firstrow2 = colnames(protein2)
    protein2 = rbind(firstrow2, protein2)
    colnames(protein2) = c("Sequence", "Name", "Length", "Database", "ID", "Description", "Start", "End", "Probability", "T/F", "Date", "Interpro_ID", "DomainName", "GO", "Desc1")
    protein2 = protein2[protein2$Database == "Pfam", ]
  
    for(domain1 in 1:nrow(protein1)){
      for(domain2 in 1:nrow(protein2)){
        if(protein1$ID[domain1] %in% compiled_interactions$domain_1){
          poi1 = compiled_interactions[compiled_interactions$domain_1 == protein1$ID[domain1], ]
          if(protein2$ID[domain2] %in% poi1$domain_2){
            domain_name_1 = c(domain_name_1, protein1$DomainName[domain1])
            domain_start_1 = c(domain_start_1, protein1$Start[domain1])
            domain_end_1 = c(domain_end_1, protein1$End[domain1])
            
            domain_name_2 = c(domain_name_2, protein2$DomainName[domain2])
            domain_start_2 = c(domain_start_2, protein2$Start[domain2])
            domain_end_2 = c(domain_end_2, protein2$End[domain2])
          }
        }
        if (protein1$ID[domain1] %in% compiled_interactions$domain_2){
          poi1 = compiled_interactions[compiled_interactions$domain_2 == protein1$ID[domain1], ]
          if(protein2$ID[domain2] %in% poi1$domain_1){
            domain_name_1 = c(domain_name_1, protein1$DomainName[domain1])
            domain_start_1 = c(domain_start_1, protein1$Start[domain1])
            domain_end_1 = c(domain_end_1, protein1$End[domain1])
            
            domain_name_2 = c(domain_name_2, protein2$DomainName[domain2])
            domain_start_2 = c(domain_start_2, protein2$Start[domain2])
            domain_end_2 = c(domain_end_2, protein2$End[domain2])
          }
        }
      }
    }
  }
  range1 = paste(domain_start_1, domain_end_1, sep = "-")
  range2 = paste(domain_start_2, domain_end_2, sep = "-")
  final_table = data.frame(domain_name_1, domain_name_2, range1, range2)
  colnames(final_table) <- c("Domain 1", "Domain 2", "Domain 1 Range", "Domain 2 Range")
  return(final_table)
}


ui <- fluidPage(
  titlePanel("Domain-Domain Interaction Prediction"),
  tabsetPanel(
    tabPanel("Predict",
             sidebarLayout(
               sidebarPanel(
                 fileInput("file1", "Choose TSV File for Protein 1:", accept = ".tsv"),
                 fileInput("file2", "Choose TSV File for Protein 2:", accept = ".tsv"),
                 radioButtons("database", "Select Database",
                              choices = list("InterPro" = "interpro", "ELM" = "elm")),
                 downloadButton("downloadData", "Download Results")
               ),
               mainPanel(
                 DTOutput("dataTable")
               )
             )
    ),
    tabPanel("About",
             h3("About"),
             h4("This work?"),
             h3("Help"),
             h3("References"),
             p("Test"),
             p("Test")
    )
  )
)

server <- function(input, output, session){
  results <- reactive({
    req(input$file1)
    req(input$file2)
    
    protein1 <- read_tsv(input$file1$datapath)
    protein2 <- read_tsv(input$file2$datapath)
    
    if (input$database == "interpro") {
      interpro_algo(protein1, protein2)
    } else if (input$database == "elm") {
      elm_algo(protein1, protein2)
    }
  })
  
  output$dataTable <- renderDT({
    results()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("ddip-prediction", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      write_tsv(results(), file)
    }
  )
}

shinyApp(ui = ui, server = server)
