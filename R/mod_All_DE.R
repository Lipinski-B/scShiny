mod_All_DE_ui <- function(id) {
  ns <- NS(id)
  
  # metadata DE
  tabItem(tabName = "All_DE",
          navbarPage("Effet : ",
                      tabPanel("VlnPlot",
                        fluidRow(
                          column(8,plotOutput(ns("Vln_plot"), width = "100%",  height = "600px")),
                          column(4,box(width = 12,
                            column(12,h4("Parameters :")),
                            column(12,pickerInput(inputId = ns("Vln_metadata"),label = "Métadata" , choices = c("Sample", "Heavy", "Light", "Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST","CDR3_AA_Heavy_TRUST",
                                                                                                                                                              "V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST","CDR3_AA_Light_TRUST",
                                                       "State","Condition","Phénotype","Phénotype.fine","Phase","Reponse","Greffe","Clonotype","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light","Chaine_Heavy", 
                                                                                                                                                              "V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy",
                                                                                                                                          "BCL2_L23L","BCL2_K22K","CD79B_Y196H","EZH2_A682G","EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S"), multiple = F, options = list(`actions-box` = TRUE))),
                            column(12,pickerInput(inputId = ns("Vln_metadata_split"),label = "Métadata" , choices = c("Sample", "Heavy", "Light", "Phénotype_TRUST","V_Heavy_TRUST","D_Heavy_TRUST","J_Heavy_TRUST","Heavy_TRUST","CDR3_DNA_Heavy_TRUST","CDR3_AA_Heavy_TRUST",
                                                                                                                "V_Light_TRUST","D_Light_TRUST","J_Light_TRUST","Light_TRUST","CDR3_DNA_Light_TRUST","CDR3_AA_Light_TRUST",
                                                                                                                "State","Condition","Phénotype","Phénotype.fine","Phase","Reponse","Greffe","Clonotype","Chaine_Light","V_Light","D_Light","J_Light","C_Light","CDR3_Light","CDR3_nt_Light","Chaine_Heavy", 
                                                                                                                "V_Heavy","D_Heavy","J_Heavy","C_Heavy","CDR3_Heavy","CDR3_nt_Heavy",
                                                                                                                "BCL2_L23L","BCL2_K22K","CD79B_Y196H","EZH2_A682G","EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S"), multiple = T, options = list(`actions-box` = TRUE))),
                            
                            column(12,selectizeInput(inputId = ns('Vln_feature'), label= 'Feature : ', choices = NULL, multiple=TRUE)),
                            column(12,selectizeInput(inputId = ns('Vln_ident'), label= 'Ident : ', choices = NULL, multiple=TRUE)),
                            column(12,awesomeRadio(inputId = ns('Vln_assay'), label= 'Assay : ', choices = c('RNA','SCT','integrated'),inline = T)),
                            column(12,awesomeRadio(inputId = ns('Vln_slot'), label= 'Slot : ', choices = c('data','count','scale.data'),inline = T)),
                          ))
                        )
                      ),
                      tabPanel("Réponse",h3("RC vs RP :"),
                        fluidRow(
                          column(6,
                            tagList(DT::DTOutput(ns("result_Reponse"))),br(),
                            tagList(plotOutput(ns("DE_Heatmap_Reponse"), width = "95%",  height = "600px")),
                          ),
                          column(6,
                            tagList(plotOutput(ns("DE_RidgePlot_Reponse"), width = "95%",  height = "600px")),br(),
                            tagList(plotOutput(ns("DE_VlnPlot_Reponse"), width = "95%",  height = "600px"))
                          )
                        )
                      ),
                      tabPanel("RCHOP",h3("RCHOP vs Excipient :"),            
                        fluidRow(
                          column(6,
                            tagList(DT::DTOutput(ns("result_RE"))),br(),
                            tagList(plotOutput(ns("DE_Heatmap_RE"), width = "95%",  height = "600px")),
                          ),
                          column(6,
                            tagList(plotOutput(ns("DE_RidgePlot_RE"), width = "95%",  height = "600px")),br(),
                            tagList(plotOutput(ns("DE_VlnPlot_RE"), width = "95%",  height = "600px"))
                          )
                        )
                      ),
                              
          
                      tabPanel("Pré-greffe",h3("Pré-greffe vs Excipient :"),
                        fluidRow(
                          column(6,
                            tagList(DT::DTOutput(ns("result_PE"))),br(),
                            tagList(plotOutput(ns("DE_Heatmap_PE"), width = "95%",  height = "600px")),
                          ),
                          column(6,
                            tagList(plotOutput(ns("DE_RidgePlot_PE"), width = "95%",  height = "600px")),br(),
                            tagList(plotOutput(ns("DE_VlnPlot_PE"), width = "95%",  height = "600px"))
                          )
                        )
                      )
            

                      
                      # tabPanel("Resultats",
                      #        splitLayout(cellWidths = c("40%", "60%"),
                                       #tagList(DT::DTOutput(ns("result"))),
                                       #tagList(plotOutput(ns("intersect"), width = "90%",  height = "650px")))),
                     
                      # tabPanel("Intersection", 
                      #         splitLayout(cellWidths = c("50%", "50%"),
                      #                  plotOutput(ns("heatmap1"), width = "90%",  height = "850px"),
                      #                  plotOutput(ns("heatmap2"), width = "90%",  height = "850px")))
    
  ))
}

mod_All_DE_server <- function(input, output, session, r = r) {
  ns <- session$ns
  

  
  output$result_RE <-  DT::renderDataTable({DT::datatable(as.data.frame(format(r$dataset@tools$DE_RE,scientific =T)))})
  output$DE_Heatmap_RE <- renderPlot({Seurat::Idents(r$dataset)<-"Condition"; Seurat::DoHeatmap(r$dataset, cells = rownames(r$dataset@meta.data)[which(r$dataset@meta.data$Condition==c("RCHOP","Excipient"))], features = rownames(r$dataset@tools$DE_RE)[1:50], size = 3, assay = 'integrated', slot = "data")})
  output$DE_RidgePlot_RE <- renderPlot({Seurat::Idents(r$dataset)<-"Condition"; Seurat::RidgePlot(r$dataset, idents = c("RCHOP","Excipient"), features = rownames(r$dataset@tools$DE_RE)[1:12], ncol = 3, assay = 'integrated', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")})
  output$DE_VlnPlot_RE <- renderPlot({Seurat::Idents(r$dataset)<-"Condition"; Seurat::VlnPlot(r$dataset, idents = c("RCHOP","Excipient"), features = rownames(r$dataset@tools$DE_RE)[1:12], sort = T, ncol = 3, assay = 'integrated', slot = "data") & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")})
  
  output$result_PE <-  DT::renderDataTable({DT::datatable(as.data.frame(format(r$dataset@tools$DE_PE,scientific =T)))})
  output$DE_Heatmap_PE <- renderPlot({Seurat::Idents(r$dataset)<-"Condition"; Seurat::DoHeatmap(r$dataset, cells = rownames(r$dataset@meta.data)[which(r$dataset@meta.data$Condition==c("Pré-greffe","Excipient"))], features = rownames(r$dataset@tools$DE_PE)[1:50], size = 3, assay = 'integrated', slot = "data")})
  output$DE_RidgePlot_PE <- renderPlot({Seurat::Idents(r$dataset)<-"Condition"; Seurat::RidgePlot(r$dataset, idents = c("Pré-greffe","Excipient"), features = rownames(r$dataset@tools$DE_PE)[1:12], ncol = 3, assay = 'integrated', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")})
  output$DE_VlnPlot_PE <- renderPlot({Seurat::Idents(r$dataset)<-"Condition"; Seurat::VlnPlot(r$dataset, idents = c("Pré-greffe","Excipient"), features = rownames(r$dataset@tools$DE_PE)[1:12], sort = T, ncol = 3, assay = 'integrated', slot = "data") & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")})
  
  
  output$result_Reponse <-  DT::renderDataTable({DT::datatable(as.data.frame(format(r$dataset@tools$DE_Reponse,scientific =T)))})
  output$DE_Heatmap_Reponse <- renderPlot({Seurat::Idents(r$dataset)<-"Reponse" ; Seurat::DoHeatmap(r$dataset, cells = rownames(r$dataset@meta.data)[which(r$dataset@meta.data$Reponse==c("RC","RP"))], features = rownames(r$dataset@tools$DE_Reponse)[1:50], size = 3, assay = 'integrated', slot = "data")})
  output$DE_RidgePlot_Reponse <- renderPlot({Seurat::Idents(r$dataset)<-"Reponse" ; Seurat::RidgePlot(r$dataset, idents = c("RC","RP"), features = rownames(r$dataset@tools$DE_Reponse)[1:12], ncol = 3, assay = 'integrated', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")})
  output$DE_VlnPlot_Reponse <- renderPlot({Seurat::Idents(r$dataset)<-"Reponse" ; Seurat::VlnPlot(r$dataset, idents = c("RC","RP"), features = rownames(r$dataset@tools$DE_Reponse)[1:12], sort = T, ncol = 3, assay = 'integrated', slot = "data") & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")})
  
  
  reactive({input$Vln_metadata})
  reactive({input$Vln_assay})
  reactive({input$Vln_slot})
  reactive({input$Vln_feature})
  reactive({input$Vln_ident})
  reactive({input$Vln_metadata_split})
  
  observe({updateSelectizeInput(session, 'Vln_feature', choices = c("RCHOP_SCT_Score1", "RCHOP_RNA_Score1",  "RCHOP_i_Score1", "RCHOP_AUC_Score","percent.mt","percent.ig","percent.rb","percent.rb.meta","percent.ig.meta","percent.mt.meta","nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT","nCount_HTO","nFeature_HTO", r$feature), selected = "RCHOP_RNA_Score1", server = TRUE)})
  observe({updateSelectizeInput(session, 'Vln_ident', choices = names(table(r$dataset[[as.character(input$Vln_metadata)]])), server = TRUE)})

  output$Vln_plot = renderPlot({
    Seurat::Idents(r$dataset) <- input$Vln_metadata ; Seurat::VlnPlot(r$dataset, idents = input$Vln_ident, split.by= input$Vln_metadata_split, features = input$Vln_feature, assay = input$Vln_assay, slot = input$Vln_slot) & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
  })
  
  #load(file="inst/app/www/DE.RData")
  #count <- cluster_counts[which(rownames(cluster_counts) %in% rownames(res)),]
  # output$intersect <- renderPlot({
  #   sign <- list(Bulk = signature_RE_bulk, Single_cell = signature_RE_sc)
  #   ggVennDiagram::ggVennDiagram(sign, label_alpha = 1, label = "count") + ggplot2::scale_fill_gradient(low="white",high = "white") + ggplot2::theme(legend.position = "none") + ggplot2::ggtitle("Expression différentielle :\nRCHOP vs Excipient") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # })
  # output$heatmap1 <- renderPlot({pheatmap::pheatmap(count[ , 1:length(colnames(count))], color = RColorBrewer::brewer.pal(6, "YlOrRd") , cluster_rows = T, show_rownames = F, annotation = cluster_metadata[, c("patient_id", "group_id", "cluster_id")], cluster_cols = F, border_color = NA, fontsize = 10, scale = "row", fontsize_row = 10, height = 20)})  
  # output$heatmap2 <- renderPlot({
  #   count <- count[, c(1,3,5,7,9,11,13,2,4,6,8,10,12,14)]
  #   pheatmap::pheatmap(count[ , 1:length(colnames(count))], color = RColorBrewer::brewer.pal(6, "YlOrRd") , cluster_rows = T, show_rownames = F, annotation = cluster_metadata[, c("patient_id", "group_id", "cluster_id")], cluster_cols = F, border_color = NA, fontsize = 10, scale = "row", fontsize_row = 10, height = 20)
  # })

}

## To be copied in the UI
# mod_All_DE_ui("All_DE_ui_1")

## To be copied in the server
# callModule(mod_All_DE_server, "All_DE_ui_1")