mod_All_DE_ui <- function(id) {
  ns <- NS(id)
  
  # metadata DE
  tabItem(tabName = "All_DE",
          navbarPage("Analyses : ",

                      tabPanel("Heatmap",h3("RCHOP vs Excipient :"),
                               splitLayout(cellWidths = c("40%", "60%"),
                                           tagList(DT::DTOutput(ns("result"))),
                                           tagList(plotOutput(ns("DE_Heatmap"), width = "95%",  height = "600px")))
                               ),
                              
          
                      tabPanel("Plot",h3("RCHOP vs Excipient :"),
                               splitLayout(cellWidths=c("50%","50%"),
                                           tagList(plotOutput(ns("DE_RidgePlot"), width = "95%",  height = "1200px")),
                                           tagList(plotOutput(ns("DE_VlnPlot"), width = "95%",  height = "1200px")))
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
    
  )
}

mod_All_DE_server <- function(input, output, session, r = r) {
  ns <- session$ns
  
  observe(Seurat::Idents(r$dataset)<-"Condition")
  
  output$result <-  DT::renderDataTable({DT::datatable(as.data.frame(format(r$dataset@tools$DE_RE,scientific =T)))})
  output$DE_Heatmap <- renderPlot({Seurat::DoHeatmap(r$dataset, cells = rownames(r$dataset@meta.data)[which(r$dataset@meta.data$Condition==c("Pré-greffe","Excipient"))], features = rownames(r$dataset@tools$DE_RE)[1:50], size = 3, assay = 'integrated', slot = "data")})
  output$DE_RidgePlot <- renderPlot({Seurat::RidgePlot(r$dataset, idents = c("Pré-greffe","Excipient"), features = rownames(r$dataset@tools$DE_RE)[1:12], ncol = 3, assay = 'integrated', slot = "data") & theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")})
  output$DE_VlnPlot <- renderPlot({Seurat::VlnPlot(r$dataset, idents = c("Pré-greffe","Excipient"), features = rownames(r$dataset@tools$DE_RE)[1:12], sort = T, ncol = 3, assay = 'integrated', slot = "data") & theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom")})
  
  
  
  
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