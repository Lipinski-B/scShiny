mod_Patient_Reduction_Dimension_ui <- function(id) {
  ns <- NS(id)

  tabItem(tabName = "Patient_Reduction_Dimension",          
          #a(href = "https://www.dublinbus.ie/RTPI/Sources-of-Real-Time-Information/", "Bus stop numbers can be found here.", target = "_blank"),
          #br(),h3("Custom URL"),
          #p("A custom URL can be used to pre select choices when loading the app. Use the button below to create a URL for the choices currently selected."),
          splitLayout(cellWidths = c("80%", "20%"),
            tagList(
              conditionalPanel(condition = "input.Analyse_type == 'RNA'", ns = ns,
                plotOutput(ns("PCA"), width = "100%",  height = "650px"), uiOutput(ns("feature_pca")),conditionalPanel(condition = "input.genes23 == true",plotOutput(ns("FeaturePlot_23_PCA"), width = "100%",  height = "650px")),),
              conditionalPanel(condition = "input.Analyse_type == 'Velocity'", ns = ns,
                plotOutput(ns("Velocity"), width = "100%",  height = "1000px"))),
            
            tagList(
              column(10,shinyBS::bsCollapse(id = "collapse", open = "Analyse : ",
                shinyBS::bsCollapsePanel("Analyse : ",br(),
                column(12,awesomeRadio(inputId = ns("Analyse_type"), label = NULL, choices = c("RNA","Velocity"), selected = "RNA", inline = F, checkbox = TRUE))
              ))),
              
              column(10,shinyBS::bsCollapse(id = ns("collapse"), open = "Paramètres : ",
                shinyBS::bsCollapsePanel("Paramètres : ",
                
                conditionalPanel(condition = sprintf("input['%s'] == 'RNA'", "Analyse_type"), ns = ns,
                  column(12,h5("Reduction :")),
                  column(12,pickerInput(inputId = ns("Reduction"),label = NULL , choices = c("umap","pca","tsne"), multiple = F, options = list(`actions-box` = TRUE))),fluidRow(),
                  column(12,h5("Group by : ")),
                  column(12,pickerInput(ns("Groupes"), inline = FALSE, choices = c("Phénotype", "Phénotype.fine","Phase","Condition", "Clonotype","Chaine","V","D","J","C","CDR3","BCL2_L23L", "BCL2_K22K", "CD79B_Y196H", "EZH2_A682G", "EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S"), multiple = TRUE, options = list(`actions-box` = TRUE))),fluidRow(),
                  column(12,h5("Split by : ")),
                  column(12,pickerInput(ns("Splites"), inline = FALSE, choices = c("seurat_clusters", "Phénotype", "Phénotype.fine","Phase","Condition", "Clonotype","Chaine","V","D","J","C","CDR3","BCL2_L23L", "BCL2_K22K", "CD79B_Y196H", "EZH2_A682G", "EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S"), multiple = TRUE, options = list(`actions-box` = TRUE))),fluidRow(),
                  
                  column(12,sliderTextInput(ns("Point_Size"),"Taille des points :",choices = seq(from = 0,to = 2,by = 0.2),grid = TRUE,selected = 0.8)),fluidRow(),
                  column(12,awesomeCheckbox(ns("Label_element"),"Show Label",TRUE)),fluidRow(),
                  column(12,selectizeInput(inputId = ns('variables'), label= 'Expression : ', choices = NULL, multiple=TRUE)),
                  column(12,materialSwitch(inputId = "genes23", label = "23 genes", status = "primary", value = FALSE)),br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br()
                ),
                
                conditionalPanel(condition = sprintf("input['%s'] == 'Velocity'", "Analyse_type"), ns = ns,
                  column(12,h5("Analyse : ")),
                  column(12,pickerInput(inputId = ns("Velo_type"),label = NULL , choices = c("Reads", "Spliced"), multiple = F, options = list(`actions-box` = TRUE))),fluidRow(),
                  column(12,h5("Reduction :")),
                  column(12,pickerInput(inputId = ns("Reduction_Velo"),label = NULL , choices = c("umap","tsne"), multiple = F, options = list(`actions-box` = TRUE))),fluidRow(),
                  column(12,h5("Condition :")),
                  column(12,pickerInput(inputId = ns("Velo_condition"),label = NULL , choices = c("RCHOP","Excipient","Pré-greffe"), multiple = F, options = list(`actions-box` = TRUE))),fluidRow(),
                  column(12,h5("Group by : ")),
                  column(12,pickerInput(inputId = ns("Groupes_Velo"),label = NULL , choices = c("seurat_clusters", "Phénotype", "Phénotype.fine","Phase","Condition", "Clonotype","Chaine","V","D","J","C","CDR3"), multiple = F, options = list(`actions-box` = TRUE))),fluidRow(),
                  
                  column(12,sliderTextInput(ns("nFeature_spliced_Velo"),"Subset :",choices = seq(from = 10,to = 200,by = 10),grid = TRUE,selected = 10)),fluidRow(),
                  
                  column(12,sliderTextInput(ns("Arrow_scale_Velo"),"Taille des flèches :",choices = seq(from = 0,to = 1,by = 0.1),grid = TRUE,selected = 0.8)),fluidRow(),
                  column(12,sliderTextInput(ns("nb_grid_Velo"),"Nombre de flèche :",choices = seq(from = 10,to = 100,by = 10),grid = TRUE,selected = 40)),fluidRow(),
                  column(12,sliderTextInput(ns("arrow_Velo"),"Poids des flèches :",choices = seq(from = 0,to = 2,by = 0.2),grid = TRUE,selected = 1)),fluidRow(),fluidRow(),
                  column(12,sliderTextInput(ns("border_cell_Velo"),"Bordures des points :",choices = seq(from = 0,to = 1,by = 0.1),grid = TRUE,selected = 0.1)),fluidRow(),br(),
                  column(12,actionButton(inputId = ns("actBtnVelocity"), label = "Submit",icon = icon("play"), width ='87%')),fluidRow(), br()
                )))))))
}

mod_Patient_Reduction_Dimension_server <- function(input, output, session, r, p) {
  ns <- session$ns

  reactive({input$Groupes})
  reactive({input$Splites})
  reactive({input$variables})
  reactive({input$Reduction})
  reactive({input$Point_Size})
  
  observe({updateSelectizeInput(session, 'variables', choices = c("percent.mt","percent.ig","percent.rb","nCount_RNA","nFeature_RNA","nCount_SCT","nFeature_SCT","nCount_HTO","nFeature_HTO", r$feature), server = TRUE)})
  
  Groupes <- eventReactive(input$actBtnVelocity,{input$Groupes_Velo})
  Reduction <- eventReactive(input$actBtnVelocity,{input$Reduction_Velo})
  nFeature_spliced <- eventReactive(input$actBtnVelocity,{input$nFeature_spliced_Velo})
  
  Point_Size <- eventReactive(input$actBtnVelocity,{input$Point_Size_Velo})
  Arrow_scale <- eventReactive(input$actBtnVelocity,{input$Arrow_scale_Velo})
  nb_grid <- eventReactive(input$actBtnVelocity,{input$nb_grid_Velo})
  arrow <- eventReactive(input$actBtnVelocity,{input$arrow_Velo})
  border_cell <- eventReactive(input$actBtnVelocity,{input$border_cell_Velo})
  Velo_condition <- eventReactive(input$actBtnVelocity,{input$Velo_condition})
  
  ## -- Réduction de dimensions -- ##
  output$PCA <- renderPlot({Seurat::DimPlot(object = r$dataset, group.by = input$Groupes, split.by = input$Splites, label.size = 5, pt.size = input$Point_Size, reduction = input$Reduction, label = input$Label_element)   & Seurat::NoLegend() & xlab(label = paste0(input$Reduction, " / PCA 1 : ", round(Seurat::Stdev(r$dataset[["pca"]])[1],2), " %")) & ylab(label = paste0(input$Reduction, " / PCA 2 : ", round(Seurat::Stdev(r$dataset[["pca"]])[2],2), " %"))})
  
  ## -- FeaturePlot -- ##
  output$feature_pca = renderUI({if(length(input$variables)>0){plotOutput(ns("FeaturePlot_PCA"), width = "100%",  height = "500px")}})
  output$FeaturePlot_PCA <- renderPlot({Seurat::FeaturePlot(r$dataset, features = input$variables, reduction = input$Reduction, split.by = input$Splites, pt.size = input$Point_Size, combine = T , ncol = 3) & Seurat::NoAxes()})
  output$FeaturePlot_23_PCA <- renderPlot({Seurat::FeaturePlot(r$dataset, features = "genes231", reduction = input$Reduction, split.by = input$Splites, pt.size = 1, combine = T, label = TRUE, repel = TRUE) & Seurat::NoAxes() & scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))})
  
  ## -- Vélocity -- ##
  observeEvent(input$actBtnVelocity, {
    if(input$Analyse_type == "Velocity"){
      
      if(input$Velo_type == "Spliced"){
        rm(list = ls());gc();gc();gc()
        if(input$Velo_condition == "RCHOP"){
          object <<- get(load(file = paste0("datasets/Velocity/",p$patient,"_RCHOP_velocyto.RData")))
        } else { 
          object <<- get(load(file = paste0("datasets/Velocity/",p$patient,"_Excipient_velocyto.RData")))
        }
        
        velo <- Seurat::Tool(object = object, slot = "RunVelocity")
        velo$current <- velo$current[,which(!is.na(Matrix::colSums(velo$current)))] 
        velo$projected <- velo$projected[,which(!is.na(Matrix::colSums(velo$projected)))]
        
        Seurat::Idents(object) <- Groupes()
        ident.colors <- (scales::hue_pal())(n = length(x = levels(x = object))) 
        names(x = ident.colors) <- levels(x = object)
        cell.colors <- ident.colors[Seurat::Idents(object = object)] 
        names(x = cell.colors) <- colnames(x = object)
        emb <<- Seurat::Embeddings(object = subset(object, nFeature_spliced > nFeature_spliced()), reduction = Reduction())
      }
      

      
      if(input$Velo_type == "Reads"){
        velo <- Seurat::Tool(object = r$dataset, slot = "RunVelocity")
        velo$current <- velo$current[,which(!is.na(Matrix::colSums(velo$current)))] 
        velo$projected <- velo$projected[,which(!is.na(Matrix::colSums(velo$projected)))]
        
        Seurat::Idents(r$dataset) <- Groupes()
        ident.colors <- (scales::hue_pal())(n = length(x = levels(x = r$dataset))) 
        names(x = ident.colors) <- levels(x = r$dataset)
        cell.colors <- ident.colors[Seurat::Idents(object = r$dataset)] 
        names(x = cell.colors) <- colnames(x = r$dataset)
        emb <<- Seurat::Embeddings(object = subset(r$dataset, nFeature_spliced > nFeature_spliced()), reduction = Reduction())
      } 
      
      
      
      output$Velocity <-renderPlot({
        velocyto.R::show.velocity.on.embedding.cor(emb = emb, 
                                                   vel = velo, n = 200, scale = "sqrt", 
                                                   cell.colors = velocyto.R::ac(x = cell.colors, alpha = 0.5), 
                                                   show.grid.flow = TRUE, do.par = T, min.grid.cell.mass = 0.5, xlab = Velo_condition(),
                                                   cex = Point_Size(), #1, 
                                                   arrow.scale = Arrow_scale(), #0.8,
                                                   grid.n = nb_grid(), #50, 
                                                   arrow.lwd = arrow(), #1,  
                                                   cell.border.alpha = border_cell()) #0.1)
        })
      }
  })
}

## To be copied in the UI
# mod_Patient_Reduction_Dimension_ui("Reduction_Dimension_ui_1")

## To be copied in the server
# callModule(mod_Patient_Reduction_Dimension_server, "Reduction_Dimension_ui_1")