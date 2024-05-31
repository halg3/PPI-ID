library(shiny)
library(DT)
library(readr)
library(dplyr)
library(r3dmol)
library(bio3d)

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

# Function that identifies which tool (InterPro or ELM) was used to generate protein domain/slim information.
# InterPro is assigned 1, ELM is assigned 0, and data that fits neither InterPro nor ELM is assigned -5.
db_id = function(protein){
  if(ncol(protein) == 15){
    db = 1
  } else if(ncol(protein) == 10){
    db = 0
  } else {
    db = -5
  }
  
  return(db)
}

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
    
    protein1_range = gsub(" \\[A\\]", ",", protein1$range)
    for(i in 1:length(protein1_range)){
      if(substr(protein1_range[i], nchar(protein1_range[i]), nchar(protein1_range[i])) == ","){
        protein1_range[i] = substr(protein1_range[i], 1, nchar(protein1_range[i])-1)
      }
    }
    protein1$range = NULL
    protein1$range = protein1_range
    protein2_range = gsub(" \\[A\\]", ",", protein2$range)
    for(i in 1:length(protein2_range)){
      if(substr(protein2_range[i], nchar(protein2_range[i]), nchar(protein2_range[i])) == ","){
        protein2_range[i] = substr(protein2_range[i], 1, nchar(protein2_range[i])-1)
      }
    }
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
  colnames(final_table) <- c("SLiM 1", "SLiM 2", "Range 1", "Range 2", "Combined Probability")
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
  combined_prob = c()
  
  if(ncol(protein1) != 15 | ncol(protein2) != 15){
    return("Please ensure that Interpro results for both proteins are submitted.")
  } else {
    firstrow1 = colnames(protein1)
    protein1 = rbind(firstrow1, protein1)
    colnames(protein1) = c("Sequence", "Name", "Length", "Database", "ID", "Description", "Start", "End", "prob", "T/F", "Date", "Interpro_ID", "DomainName", "GO", "Desc1")
    protein1 = protein1[protein1$Database == "Pfam", ]
    
    firstrow2 = colnames(protein2)
    protein2 = rbind(firstrow2, protein2)
    colnames(protein2) = c("Sequence", "Name", "Length", "Database", "ID", "Description", "Start", "End", "prob", "T/F", "Date", "Interpro_ID", "DomainName", "GO", "Desc1")
    protein2 = protein2[protein2$Database == "Pfam", ]
    
    for(domain1 in 1:nrow(protein1)){
      for(domain2 in 1:nrow(protein2)){
        if(protein1$ID[domain1] %in% compiled_interactions$domain_1){
          poi1 = compiled_interactions[compiled_interactions$domain_1 == protein1$ID[domain1], ]
          if(protein2$ID[domain2] %in% poi1$domain_2){
            prob_1 = c()
            prob_2 = c()
            
            domain_name_1 = c(domain_name_1, protein1$DomainName[domain1])
            domain_start_1 = c(domain_start_1, protein1$Start[domain1])
            domain_end_1 = c(domain_end_1, protein1$End[domain1])
            prob_1 = c(prob_1, as.numeric(protein1$prob[domain1]))
            
            domain_name_2 = c(domain_name_2, protein2$DomainName[domain2])
            domain_start_2 = c(domain_start_2, protein2$Start[domain2])
            domain_end_2 = c(domain_end_2, protein2$End[domain2])
            prob_2 = c(prob_2, as.numeric(protein2$prob[domain2]))
            
            combined_prob = c(combined_prob, prob_1 * prob_2)
          }
        }
        if (protein1$ID[domain1] %in% compiled_interactions$domain_2){
          poi1 = compiled_interactions[compiled_interactions$domain_2 == protein1$ID[domain1], ]
          if(protein2$ID[domain2] %in% poi1$domain_1){
            prob_1 = c()
            prob_2 = c()
            
            domain_name_1 = c(domain_name_1, protein1$DomainName[domain1])
            domain_start_1 = c(domain_start_1, protein1$Start[domain1])
            domain_end_1 = c(domain_end_1, protein1$End[domain1])
            prob_1 = c(prob_1, as.numeric(protein1$prob[domain1]))
            
            domain_name_2 = c(domain_name_2, protein2$DomainName[domain2])
            domain_start_2 = c(domain_start_2, protein2$Start[domain2])
            domain_end_2 = c(domain_end_2, protein2$End[domain2])
            prob_2 = c(prob_2, as.numeric(protein2$prob[domain2]))
            
            combined_prob = c(combined_prob, prob_1 * prob_2)
          }
        }
      }
    }
  }
  range1 = paste(domain_start_1, domain_end_1, sep = "-")
  range2 = paste(domain_start_2, domain_end_2, sep = "-")
  
  unsorted_table = unique(data.frame(domain_name_1, domain_name_2, range1, range2, combined_prob))
  
  indices = order(unsorted_table$combined_prob, decreasing = FALSE)
  final_table = unsorted_table[indices,]
  colnames(final_table) <- c("Domain 1", "Domain 2", "Range 1", "Range 2", "Combined Probability")
  row.names(final_table) = NULL
  return(final_table)
}

mixed_algo = function(protein1, protein2, prot_db_1){
  if(prot_db_1 == 1){
    firstrow1 = colnames(protein1)
    protein1 = rbind(firstrow1, protein1)
    colnames(protein1) = c("Sequence", "Name", "Length", "Database", "ID", "Description", "Start", "End", "prob", "T/F", "Date", "Interpro_ID", "DomainName", "GO", "Desc1")
    protein1 = protein1[protein1$Database == "Pfam", ]
    
    elm_table_names = c("slim", "sequence", "range", "misc", "description", "location", "pattern", "phi", "structure", "probability")
    elm_default_table_names = c("Elm Name", "Instances (Matched Sequence)", "Positions", "View in Jmol", "Elm Description", "Cell Compartment", "Pattern", "PHI-Blast Instance Mapping", "Structural Filter Info", "Probability")
    
    if(sum(colnames(protein2) == elm_default_table_names) == 10){
      colnames(protein2) = elm_table_names
    } else if(colnames(protein2)[1] %in% interaction_id$elm_identifier){
      firstrow2 = colnames(protein2)
      protein2 = rbind(firstrow2, protein2)
      colnames(protein2) <- elm_table_names
    }
    
    protein2_range = gsub(" \\[A\\]", ",", protein2$range)
    for(i in 1:length(protein2_range)){
      if(substr(protein2_range[i], nchar(protein2_range[i]), nchar(protein2_range[i])) == ","){
        protein2_range[i] = substr(protein2_range[i], 1, nchar(protein2_range[i])-1)
      }
    }
    protein2$range = NULL
    protein2$range = protein2_range
    
    protein2_cleaned = slim_pfam(protein2)
    
    domain_name_1 = c()
    domain_name_2 = c()
    domain_start_1 = c()
    domain_end_1 = c()
    domain_range_2 = c()
    combined_prob = c()
    
    for(domain1 in 1:nrow(protein1)){
      for(domain2 in 1:nrow(protein2_cleaned)){
        if(protein1$ID[domain1] %in% compiled_interactions$domain_1){
          poi1 = compiled_interactions[compiled_interactions$domain_1 == protein1$ID[domain1], ]
          if(protein2_cleaned$ID[domain2] %in% poi1$domain_2){
            prob_1 = c()
            prob_2 = c()
            
            domain_name_1 = c(domain_name_1, protein1$DomainName[domain1])
            domain_start_1 = c(domain_start_1, protein1$Start[domain1])
            domain_end_1 = c(domain_end_1, protein1$End[domain1])
            prob_1 = c(prob_1, as.numeric(protein1$prob[domain1]))
            
            domain_name_2 = c(domain_name_2, protein2_cleaned$slim[domain2])
            domain_range_2 = c(domain_range_2, protein2_cleaned$range[domain2])
            prob_2 = c(prob_2, as.numeric(protein2_cleaned$prob[domain2]))
            
            combined_prob = c(combined_prob, prob_1 * prob_2)
          }
        }
        if(protein1$ID[domain1] %in% compiled_interactions$domain_2){
          poi1 = compiled_interactions[compiled_interactions$domain_2 == protein1$ID[domain1], ]
          if(protein2_cleaned$ID[domain2] %in% poi1$domain_1){
            prob_1 = c()
            prob_2 = c()
            
            domain_name_1 = c(domain_name_1, protein1$DomainName[domain1])
            domain_start_1 = c(domain_start_1, protein1$Start[domain1])
            domain_end_1 = c(domain_end_1, protein1$End[domain1])
            prob_1 = c(prob_1, as.numeric(protein1$prob[domain1]))
            
            domain_name_2 = c(domain_name_2, protein2_cleaned$slim[domain2])
            domain_range_2 = c(domain_range_2, protein2_cleaned$range[domain2])
            prob_2 = c(prob_2, as.numeric(protein2_cleaned$prob[domain2]))
            
            combined_prob = c(combined_prob, prob_1 * prob_2)
          }
        }
      }
    }
    
    domain_range_1 = paste(domain_start_1, domain_end_1, sep = "-")
    unsorted_table = unique(data.frame(domain_name_1, domain_name_2, domain_range_1, domain_range_2, combined_prob))
    indices = order(unsorted_table$combined_prob, decreasing = FALSE)
    final_table = unsorted_table[indices,]
    colnames(final_table) <- c("Domain 1", "SLiM 2", "Range 1", "Range 2", "Combined Probability")
    row.names(final_table) = NULL
    return(final_table)
    
  } else {
    firstrow2 = colnames(protein2)
    protein2 = rbind(firstrow2, protein2)
    colnames(protein2) = c("Sequence", "Name", "Length", "Database", "ID", "Description", "Start", "End", "prob", "T/F", "Date", "Interpro_ID", "DomainName", "GO", "Desc1")
    protein2 = protein2[protein2$Database == "Pfam", ]
    
    elm_table_names = c("slim", "sequence", "range", "misc", "description", "location", "pattern", "phi", "structure", "probability")
    elm_default_table_names = c("Elm Name", "Instances (Matched Sequence)", "Positions", "View in Jmol", "Elm Description", "Cell Compartment", "Pattern", "PHI-Blast Instance Mapping", "Structural Filter Info", "Probability")
    
    if(sum(colnames(protein1) == elm_default_table_names) == 10){
      colnames(protein1) = elm_table_names
    } else if(colnames(protein1)[1] %in% interaction_id$elm_identifier){
      firstrow1 = colnames(protein1)
      protein1 = rbind(firstrow1, protein1)
      colnames(protein1) <- elm_table_names
    }
    
    protein1_range = gsub(" \\[A\\]", ",", protein1$range)
    for(i in 1:length(protein1_range)){
      if(substr(protein1_range[i], nchar(protein1_range[i]), nchar(protein1_range[i])) == ","){
        protein1_range[i] = substr(protein1_range[i], 1, nchar(protein1_range[i])-1)
      }
    }
    protein1$range = NULL
    protein1$range = protein1_range
    
    protein1 = slim_pfam(protein1)
    
    domain_name_1 = c()
    domain_name_2 = c()
    domain_start_2 = c()
    domain_end_2 = c()
    domain_range_1 = c()
    combined_prob = c()
    
    for(domain1 in 1:nrow(protein1)){
      for(domain2 in 1:nrow(protein2)){
        if(protein1$ID[domain1] %in% compiled_interactions$domain_1){
          poi1 = compiled_interactions[compiled_interactions$domain_1 == protein1$ID[domain1], ]
          if(protein2$ID[domain2] %in% poi1$domain_2){
            prob_1 = c()
            prob_2 = c()
            
            domain_name_1 = c(domain_name_1, protein1$slim[domain1])
            domain_range_1 = c(domain_range_1, protein1$range[domain1])
            prob_1 = c(prob_1, as.numeric(protein1$prob[domain1]))
            
            domain_name_2 = c(domain_name_2, protein2$DomainName[domain2])
            domain_start_2 = c(domain_start_2, protein2$Start[domain2])
            domain_end_2 = c(domain_end_2, protein2$End[domain2])
            prob_2 = c(prob_2, as.numeric(protein2$prob[domain2]))
            
            combined_prob = c(combined_prob, prob_1 * prob_2)
          }
        }
        if(protein1$ID[domain1] %in% compiled_interactions$domain_2){
          poi1 = compiled_interactions[compiled_interactions$domain_2 == protein1$ID[domain1], ]
          if(protein2$ID[domain2] %in% poi1$domain_1){
            prob_1 = c()
            prob_2 = c()
            
            domain_name_1 = c(domain_name_1, protein1$slim[domain1])
            domain_range_1 = c(domain_range_1, protein1$range[domain1])
            prob_1 = c(prob_1, as.numeric(protein1$prob[domain1]))
            
            domain_name_2 = c(domain_name_2, protein2$DomainName[domain2])
            domain_start_2 = c(domain_start_2, protein2$Start[domain2])
            domain_end_2 = c(domain_end_2, protein2$End[domain2])
            prob_2 = c(prob_2, as.numeric(protein2$prob[domain2]))
            
            combined_prob = c(combined_prob, prob_1 * prob_2)
          }
        }
      }
    }
    
    domain_range_2 = paste(domain_start_2, domain_end_2, sep = "-")
    unsorted_table = unique(data.frame(domain_name_1, domain_name_2, domain_range_1, domain_range_2, combined_prob))
    indices = order(unsorted_table$combined_prob, decreasing = FALSE)
    final_table = unsorted_table[indices,]
    colnames(final_table) <- c("Domain 1", "SLiM 2", "Range 1", "Range 2", "Combined Probability")
    row.names(final_table) = NULL
    return(final_table)
  }
}

range_conversion = function(range){
  parts = unlist(strsplit(range, "-"))
  numeric_parts = as.numeric(parts)
  range_string = paste(parts[1], parts[2], sep = ":")
  range_correct = eval(parse(text = range_string))
  
  return(range_correct)
}

multi_range_conversion = function(range){
  ranges_correct = c()
  parts = unlist(strsplit(range, ", "))
  
  for(i in 1:length(parts)){
    range_part = range_conversion(parts[i])
    ranges_correct = c(ranges_correct, range_part)
  }
  
  return(ranges_correct)
}

ui <- fluidPage(
  titlePanel("Domain-Domain Interaction Prediction"),
  tabsetPanel(
    tabPanel("Predict",
             sidebarLayout(
               sidebarPanel(
                 fileInput("file1", "Choose TSV File for Protein 1:", accept = ".tsv"),
                 fileInput("file2", "Choose TSV File for Protein 2:", accept = ".tsv"),
                 downloadButton("downloadData", "Download Results"),
                 tags$hr(style = "border-top: 1px solid #444444;"),
                 numericInput("numericFilter", "Only Show Probabilities Less Than...", min = 0, max = 1, value = NULL, step = 1e-10),
                 tags$hr(style = "border-top: 1px solid #444444;"),
                 fileInput("pdb", "Choose PDB file for Multimer", accept = ".pdb"),
                 actionButton("refreshPdb", "Revert Model to Original"),
                 tags$hr(style = "border-top: 1px solid #444444;"),
                 tags$h5("3D Molecular Model Color Key:"),
                 tags$h6(tags$span(style = "color: hotpink;", "Protein 1"),
                         tags$span(style = "color: #00cc96;", "Protein 2"))
               ),
               mainPanel(
                 DTOutput("dataTable"),
                 uiOutput("message"),
                 r3dmolOutput("structure3d")
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
  error_message = reactiveVal(NULL)
  
  results <- reactive({
    req(input$file1)
    req(input$file2)
    
    protein1 <- read_tsv(input$file1$datapath)
    protein2 <- read_tsv(input$file2$datapath)
    
    prot_db_1 <- db_id(protein1)
    prot_db_2 <- db_id(protein2)
    
    if(prot_db_1 + prot_db_2 == 2){
      error_message(NULL)
      interpro_algo(protein1, protein2)
    } else if(prot_db_1 + prot_db_2 == 0){
      error_message(NULL)
      elm_algo(protein1, protein2)
    } else if(prot_db_1 + prot_db_2 == 1){
      error_message(NULL)
      mixed_algo(protein1, protein2, prot_db_1)
    } else {
      error_message("Please ensure domain/SLiM information for both proteins are from either InterPro or the ELM Database")
      return(NULL)
    }
  })
  
  output$dataTable <- renderDT({
    datatable(results(), selection = 'single')
  })
  
  observeEvent(input$numericFilter, {
    if (!is.null(input$numericFilter) && !is.na(input$numericFilter) && input$numericFilter != "") {
      filteredData <- results()
      if (!is.null(filteredData)) {
        filteredData <- filteredData[filteredData[["Combined Probability"]] < input$numericFilter, ]
        output$dataTable <- renderDT({
          datatable(filteredData, selection = 'single')
        })
      }
    } else {
      output$dataTable <- renderDT({
        resultsData <- results()
        if (!is.null(resultsData)) {
          datatable(resultsData, selection = 'single')
        }
      })
    }
  })
  
  output$message <- renderUI({
    msg <- error_message()
    if (!is.null(msg)) {
      div(style = "color: red;", msg)
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("ddip-prediction", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      write_tsv(results(), file)
    }
  )
  
  pdb_data = reactiveVal(NULL)
  chains = reactiveVal(NULL)
  
  observeEvent(input$pdb, {
    req(input$pdb)
    pdb_file <- input$pdb$datapath
    pdb_data(pdb_file)
    
    pdb_object = read.pdb(pdb_file)
    chains(unique(pdb_object$atom$chain))
    
    output$structure3d <- renderR3dmol({
      r3dmol(  
        viewer_spec = m_viewer_spec(
          cartoonQuality = 10,
          lowerZoomLimit = 50,
          upperZoomLimit = 350
        )
      ) %>%
        m_add_model(
          data = pdb_data(), 
          format = "pdb"
        ) %>%
        m_set_style(
          sel = m_sel(chain = chains()[1]),
          style = m_style_cartoon(
            color = "hotpink"
          )
        ) %>%
        m_set_style(
          sel = m_sel(chain = chains()[2]),
          style = m_style_cartoon(
            color = "#00cc96"
          )
        ) %>%
        m_zoom_to()
    })
  })
  
  observeEvent(input$dataTable_rows_selected, {
    selected_row <- input$dataTable_rows_selected
    if(length(selected_row) > 0){
      selected_data <- results() [selected_row, ]
      
      range1 = selected_data[["Range 1"]]
      range2 = selected_data[["Range 2"]]
      
      if(length(unlist(strsplit(range1, ", "))) > 1){
        range1 = multi_range_conversion(range1)
      } else {
        range1 = range_conversion(range1)
      }
      
      if(length(unlist(strsplit(range2, ", "))) > 1){
        range2 = multi_range_conversion(range2)
      } else {
        range2 = range_conversion(range2)
      }
      
      output$structure3d <- renderR3dmol({
        r3dmol(  
          viewer_spec = m_viewer_spec(
            cartoonQuality = 10,
            lowerZoomLimit = 50,
            upperZoomLimit = 350
          )
        ) %>%
          m_add_model(
            data = pdb_data(), 
            format = "pdb"
          ) %>%
          m_set_style(
            sel = m_sel(chain = chains()[1]),
            style = m_style_cartoon(
              color = "hotpink"
            )
          ) %>%
          m_set_style(
            sel = m_sel(chain = chains()[2]),
            style = m_style_cartoon(
              color = "#00cc96"
            )
          ) %>%
          m_set_style(
            sel = m_sel(chain = chains()[1], resi = range1),
            style = m_style_cartoon(
              color = "blue"
            )
          ) %>%
          m_set_style(
            sel = m_sel(chain = chains()[2], resi = range2),
            style = m_style_cartoon(
              color = "blue"
            )
          ) %>%
          m_zoom_to()
      })
    }
  })
  
  observeEvent(input$refreshPdb, {
    req(input$pdb)
    
    output$structure3d <- renderR3dmol({
      r3dmol(  
        viewer_spec = m_viewer_spec(
          cartoonQuality = 10,
          lowerZoomLimit = 50,
          upperZoomLimit = 350
        )
      ) %>%
        m_add_model(
          data = pdb_data(), 
          format = "pdb"
        ) %>%
        m_set_style(
          sel = m_sel(chain = chains()[1]),
          style = m_style_cartoon(
            color = "hotpink"
          )
        ) %>%
        m_set_style(
          sel = m_sel(chain = chains()[2]),
          style = m_style_cartoon(
            color = "#00cc96"
          )
        ) %>%
        m_zoom_to()
    })
  })
}

shinyApp(ui = ui, server = server)
