library(stringr)
library(plyr)
library(ggplot2)
library(plotly)
library(dplyr)

patient="FL12C1888"
load(file=paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/VDJ.RData"))

# Create test data.
data0 <- data.frame(
  clonotype=vloupe$clonotype_id,
  frequency=vloupe$frequency
)

data$fraction = data$frequency / sum(data$frequency) # Compute percentages
data$ymax = cumsum(data$fraction) # Compute the cumulative percentages (top of each rectangle)
data$ymin = c(0, head(data$ymax, n=-1)) # Compute the bottom of each rectangle
data$labelPosition <- (data$ymax + data$ymin) / 2 # Compute label position

# Make the plot
ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=3, xmin=1, fill=clonotype)) +
  geom_rect() +
  geom_text(x=4, aes(y=labelPosition, label=frequency, color=clonotype), size=3) + # x here controls label position (inner / outer)
  coord_polar(theta="y") + xlim(c(-1, 4)) +
  theme_void() + theme(legend.position = "none")


data %>% plot_ly(labels = ~clonotype, values = ~frequency) %>% 
         add_pie(hole = 0.6) %>% 
         layout(title = "Donut charts using Plotly",  showlegend = F,
                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


siege <- c("FL140304","FL12C1888","FL09C1164","FL08G0293")
siege <- c("FL02G095","FL05G0330")

for (patient in siege) {
  load(file=paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/VDJ.RData"))
  e <- plot_ly(hole = 0.85 ,type = "pie",labels = vloupe$clonotype_id, values = vloupe$frequency, showlegend = FALSE,
          textposition = 'side',textinfo = 'percent',
          hovertemplate = paste0('Clonotype : ', vloupe$clonotype_id, "\nProportion : ", round(vloupe$proportion,3)*100, "% \nType : ", vloupe$type, "\nIsotype : ", vloupe$igh_c_genes,"\nHeavy : ", vloupe$V_lourde, " / ", vloupe$D_lourde, " / ", vloupe$J_lourde, "\nLight : ", vloupe$V_legere, " / ", vloupe$J_legere),
          domain = list(x = c(0, 1), y = c(0, 1)) ) %>% layout(title = patient, autosize = T) %>%
    subplot(nrows = 3, shareX = F, shareY = F, margin = 0.04, heights = c(0.2,0.2,0.2), widths = c(0.2,0.2),
        plot_ly(x = V[[1]], y = V[[2]], name = "V", type = "bar", showlegend = T) %>% layout(xaxis= list(showticklabels = FALSE)),
        plot_ly(x = Heavy[[1]], y = Heavy[[2]], name = "Heavy", type = "bar", showlegend = FALSE,
                hovertemplate = paste0("Locus : %{x}\nProportion : ", round(Heavy[[2]]/sum(Heavy[[2]]),3)*100,"%")) %>% layout(xaxis= list(showticklabels = T, tickangle = 45)),
        plot_ly(x = Light[[1]], y = Light[[2]], name = "Light", type = "bar", showlegend = FALSE,
                hovertemplate = paste0("Locus : %{x}\nProportion : ", round(Light[[2]]/sum(Light[[2]]),3)*100,"%")) %>% layout(xaxis= list(showticklabels = T, tickangle = 45)),
        plot_ly(x = D[[1]], y = D[[2]], name = "D", type = "bar", showlegend = T) %>% layout(xaxis= list(showticklabels = FALSE)),
        plot_ly(x = J[[1]], y = J[[2]], name = "J", type = "bar", showlegend = T) %>% layout(xaxis= list(showticklabels = FALSE))
      )
  
  
  load(file = paste0("/home/boris/Documents/analyse/singlet_",patient,".RData"))
  f <-plot_ly(all@tools$sunburst, ids = ~ids ,labels = ~labels, parents = ~parents, values = ~values, marker = list(colors = c( "#BEBADA", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462")),
          type = 'sunburst', branchvalues = 'total', hoverinfo = "text", hovertext = paste(all@tools$sunburst$labels, ":", round((all@tools$sunburst$values/all@tools$sunburst$values[1])*100,2),"%", "\nTotal : " ,all@tools$sunburst$values)) 
  
  
  save(e, f, file = paste0("/home/boris/Bureau/scShiny/shiny/www/",patient,"_VDJ.RData"))
}

