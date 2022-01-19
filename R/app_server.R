#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd

app_server <- function( input, output, session ) {
#shinyServer(function(input, output, session) {
  #################################################################################################
  r <- reactiveValues(dataset = reactiveValues())
  p <- reactiveValues(dataset = reactiveValues())
  
  data <- reactive({singlet})
  feature <- reactive({rownames(singlet)})
  patient <- reactive({input$patient})

  observe({
    r$dataset <- data()
    r$feature <- feature()
  })
  observe({p$patient <- patient()})
  observeEvent(input$actBtnLoadPatient,{
    shinybusy::show_modal_spinner(spin = "semipolar",color = "deepskyblue",text = "Please wait...")
    
    rm(list = ls());gc();gc();gc()
    
    r$dataset <<- get(load(file = paste0("datasets/Patient/",input$Patient_Condition,"_",input$Patient_Type,"_",input$patient,".RData")))
    #r$dataset <<- get(load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/",input$patient,"/",input$Patient_Condition,"_",input$Patient_Type,"_",input$patient,".RData")))
    
    Seurat::Idents(r$dataset)<-'Phénotype.fine'
    r$dataset <- subset(r$dataset, idents = names(table(r$dataset@meta.data$Phénotype.fine)[which(table(r$dataset@meta.data$Phénotype.fine) > 10)]))
    Seurat::Idents(r$dataset)<-'Sample'
    
    shinyWidgets::sendSweetAlert(session = session,title = "Done !",type = "success")
    shinyjs::runjs("window.scrollTo(0, 50)")
    shinybusy::remove_modal_spinner()
  })
  observeEvent(input$actBtnLoadDataset,{
    shinybusy::show_modal_spinner(spin = "semipolar",color = "deepskyblue",text = "Please wait...")
    
    rm(list = ls());gc();gc();gc()
    r$dataset <<- get(load(file = paste0("datasets/All/All_",input$Dataset_Condition,"_",input$Dataset_Type,".RData")))
    #r$dataset <<- get(load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_meta/All_",input$Dataset_Condition,"_",input$Dataset_Type,".RData")))
    
    Seurat::Idents(r$dataset)<-'Sample'
    
    shinyWidgets::sendSweetAlert(session = session,title = "Done !",type = "success")
    shinyjs::runjs("window.scrollTo(0, 50)")
    shinybusy::remove_modal_spinner()
  })
  
  
  #################################################################################################
  callModule(mod_Presentation_server, "Presentation_ui_1")
  output$nb_clonotype_cell <- renderText({if (is.null(r$dataset@tools$vloupe)) {return("Sélectionnez un patient !")} else {return(paste0("Clonotype majoritaire :","\n Proportion : ", round(r$dataset@tools$vloupe$proportion[1],3)*100,"% \n Type : ", r$dataset@tools$vloupe$type[1],"\n Isotype : ", r$dataset@tools$vloupe$igh_c_genes[1],"\n Heavy : ", r$dataset@tools$vloupe$V_lourde[1], "/", r$dataset@tools$vloupe$D_lourde[1], "/", r$dataset@tools$vloupe$J_lourde[1],"\n Light : ", r$dataset@tools$vloupe$V_legere[1], "/", r$dataset@tools$vloupe$J_legere[1]))}})
  
  callModule(mod_Patient_Metadata_server, "Patient_Metadata_ui_1", r = r)
  callModule(mod_Patient_Reduction_Dimension_server, "Patient_Reduction_Dimension_ui_1", r = r, p = p)
  callModule(mod_Patient_VDJ_server, "Patient_VDJ_ui_1", r = r)
  callModule(mod_Patient_Parti_server, "Patient_Parti_ui_1", r = r, p = p)
  callModule(mod_Patient_QC_server, "Patient_QC_ui_1", r = r)
  callModule(mod_Patient_GOT_server,"Patient_GOT_ui_1", r = r)
  #callModule(mod_Patient_Expression_server, "Patient_Expression_ui_1", r = r)
  #callModule(mod_Patient_Subset_server, "Patient_Subset_ui_1")#, r = r)
  
  #################################################################################################
  output$information <- renderText({if (length(table(r$dataset@meta.data$orig.ident)) > 1) {return("Information : All")} else {return(paste("Information :", input$patient))}})
  output$nb_tot_cell <- renderText({return(paste("Number of cells : \t", as.character(as.numeric(sum(table(r$dataset@meta.data$Phénotype))))))})
  output$nb_b_cell <- renderText({return(paste("B cells : \t", as.character(as.numeric(table(r$dataset@meta.data$Phénotype)["B cells"]))))})
  output$nb_other_cell <- renderText({return(paste("\nOther cells : \t", as.character(as.numeric(sum(table(r$dataset@meta.data$Phénotype))) - as.numeric(table(r$dataset@meta.data$Phénotype)["B cells"]))))})
  output$nb_Controle_cell <- renderText({return(paste("\nControle : \t", as.character(as.numeric(table(r$dataset@meta.data$Condition)["Pré-greffe"]))))})
  output$nb_Excipient_cell <- renderText({return(paste("\nExcipient : \t", as.character(as.numeric(table(r$dataset@meta.data$Condition)["Excipient"]))))})
  output$nb_RCHOP_cell <- renderText({return(paste("\nRCHOP : \t", as.character(as.numeric(table(r$dataset@meta.data$Condition)["RCHOP"]))))})
  
  callModule(mod_All_Metadata_server, "All_Metadata_ui_1")
  callModule(mod_All_Reduction_Dimension_server, "All_Reduction_Dimension_ui_1", r = r)
  callModule(mod_All_DE_server, "All_DE_ui_1", r = r)
  #callModule(mod_All_Velocity_server, "All_Velocity_ui_1")
  #callModule(mod_All_scRepertoir_server, "All_scRepertoir_ui_1")
  #callModule(mod_All_Enrichissement_server, "All_Enrichissement_ui_1")
  
  
  #################################################################################################
  ## -- Rapport Patient -- ##
  output$PDF_report <- downloadHandler(
    filename = function() {paste0("rapport_", input$patient,".pdf")},
    content = function(file) {file.copy(paste0("inst/app/www/",input$patient,"/rapport_", input$patient,".pdf"), file)}
  )
  output$HTML_report <- downloadHandler(
    filename = function() {paste0("rapport_", input$patient,".html")},
    content = function(file) {file.copy(paste0("inst/app/www/",input$patient,"/rapport_", input$patient,".html"), file)}
  )
  output$Word_report <- downloadHandler(
    filename = function() {paste0("rapport_", input$patient,".docx")},
    content = function(file) {file.copy(paste0("inst/app/www/",input$patient,"/rapport_", input$patient,".docx"), file)}
  )
}