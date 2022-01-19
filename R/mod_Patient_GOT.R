#' Patient_GOT UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_Patient_GOT_ui <- function(id){
  ns <- NS(id)

  tabItem(tabName = "Patient_GOT",
    fluidRow(
      plotOutput(ns("histo1"), width = "100%",  height = "600px"),
      plotOutput(ns("histo2"), width = "100%",  height = "600px"),
    )
  )
}
    
#' Patient_GOT Server Functions
#'
#' @noRd 
#' 
mod_Patient_GOT_server <-function(input, output, session, r) {
    ns <- session$ns
 
    
    output$histo1 <- renderPlot({
      histo1 <- data.frame(
        hotspot=c(rep("BCL2_L23L",length(table(useNA = "always", r$dataset[["BCL2_L23L"]]))),  rep("BCL2_K22K",length(table(useNA = "always", r$dataset[["BCL2_K22K"]]))),  rep("CD79B_Y196H",length(table(useNA = "always", r$dataset[["CD79B_Y196H"]]))),  rep("EZH2_A682G",length(table(useNA = "always", r$dataset[["EZH2_A682G"]]))),  rep("EZH2_A692V",length(table(useNA = "always", r$dataset[["EZH2_A692V"]])))),  
        Génotype=c(names(table(useNA = "always", r$dataset[["BCL2_L23L"]])),    names(table(useNA = "always", r$dataset[["BCL2_K22K"]])),    names(table(useNA = "always", r$dataset[["CD79B_Y196H"]])),    names(table(useNA = "always", r$dataset[["EZH2_A682G"]])),   names(table(useNA = "always", r$dataset[["EZH2_A692V"]]))),     
        frequence=c(round(table(useNA = "always", r$dataset[["BCL2_L23L"]])/length(colnames(r$dataset))*100,2),    round(table(useNA = "always", r$dataset[["BCL2_K22K"]])/length(colnames(r$dataset))*100,2),    round(table(useNA = "always", r$dataset[["CD79B_Y196H"]])/length(colnames(r$dataset))*100,2),    round(table(useNA = "always", r$dataset[["EZH2_A682G"]])/length(colnames(r$dataset))*100,2),   round(table(useNA = "always", r$dataset[["EZH2_A692V"]])/length(colnames(r$dataset))*100,2))
      )
      
      ggplot(data=histo1, aes(x=hotspot, y=frequence, fill=Génotype)) +   geom_bar(stat="identity", color="black", position=position_dodge()) + labs(title="Génotypage du transcriptome : Hotspots BCL2, CD79B & EZH2 A") +
        xlab("Mutation") + ylab("Fréquence des génotypes") + theme_minimal() + scale_fill_manual(values=c('#250000','#E69F00', '#999999')) + scale_fill_brewer(palette="Blues") +   
        theme(axis.text.x = element_text(angle = 45, size = 12), axis.title.y = element_text(size = 18), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    })
    
    output$histo2 <- renderPlot({
      histo2 <- data.frame(
        hotspot=c(rep("EZH2_Y646C",length(table(useNA = "always", r$dataset[["EZH2_Y646C"]]))),  rep("EZH2_Y646F",length(table(useNA = "always", r$dataset[["EZH2_Y646F"]]))),  rep("EZH2_Y646H",length(table(useNA = "always", r$dataset[["EZH2_Y646H"]]))),  rep("EZH2_Y646N",length(table(useNA = "always", r$dataset[["EZH2_Y646N"]]))),  rep("EZH2_Y646S",length(table(useNA = "always", r$dataset[["EZH2_Y646S"]])))),
        Génotype=c(names(table(useNA = "always", r$dataset[["EZH2_Y646C"]])),    names(table(useNA = "always", r$dataset[["EZH2_Y646F"]])),    names(table(useNA = "always", r$dataset[["EZH2_Y646H"]])),    names(table(useNA = "always", r$dataset[["EZH2_Y646N"]])),    names(table(useNA = "always", r$dataset[["EZH2_Y646S"]]))  ),
        frequence=c(round(table(useNA = "always", r$dataset[["EZH2_Y646C"]])/length(colnames(r$dataset))*100,2),    round(table(useNA = "always", r$dataset[["EZH2_Y646F"]])/length(colnames(r$dataset))*100,2),    round(table(useNA = "always", r$dataset[["EZH2_Y646H"]])/length(colnames(r$dataset))*100,2),    round(table(useNA = "always", r$dataset[["EZH2_Y646N"]])/length(colnames(r$dataset))*100,2),    round(table(useNA = "always", r$dataset[["EZH2_Y646S"]])/length(colnames(r$dataset))*100,2)  )
      )
      
      ggplot(data=histo2, aes(x=hotspot, y=frequence, fill=Génotype)) +   geom_bar(stat="identity", color="black", position=position_dodge()) + labs(title="Génotypage du transcriptome : Hotspots EZH2 Y") +
        xlab("Mutation") + ylab("Fréquence des génotypes") + theme_minimal() + scale_fill_manual(values=c('#250000','#E69F00', '#999999')) + scale_fill_brewer(palette="Blues") +   
        theme(axis.text.x = element_text(angle = 45, size = 12), axis.title.y = element_text(size = 18), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    })
    
}
    
## To be copied in the UI
# mod_Patient_GOT_ui("Patient_GOT_1")
    
## To be copied in the server
# mod_Patient_GOT_server("Patient_GOT_1")
