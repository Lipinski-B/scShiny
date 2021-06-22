library(plotly)
library(stringr)
count <- function(column){
  names<-c();number<-c()
  for (k in 1:length(column)) {
    for (i in 1:length(names(table(vloupe[column[k]])))) {
      if(names(table(vloupe[column[k]]))[i]!="NA"){
        names <- c(names, names(table(vloupe[column[k]]))[i])
        cell <- subset(vloupe, vloupe[column[k]] == names(table(vloupe[column[k]]))[i])
        number <- c(number, sum(cell$frequency))
      }
    }
  }
  return(c(list(names), list(number)))
}

siege <- c("FL140304","FL12C1888","FL09C1164","FL08G0293")
patient <- siege[2]

vloupe <- as.data.frame(read.table( paste0("/home/boris/Documents/lipinskib/boris/VDJ/", patient,"_vloupe-clonotypes.csv"), sep = ",", header = T))
vloupe$clonotype_id <- as.numeric(str_replace(vloupe$clonotype_id,"clonotype",""))

vloupe$lourde <- ifelse(vloupe$igh_c_genes!="", "IGH", "NA")
vloupe$legere[vloupe$igl_c_genes=="" & vloupe$igk_c_genes==""] <- "NA"
vloupe$legere[vloupe$igl_c_genes=="" & vloupe$igk_c_genes!=""] <- "IGK"
vloupe$legere[vloupe$igl_c_genes!="" & vloupe$igk_c_genes==""] <- "IGL"
vloupe$legere[vloupe$igl_c_genes!="" & vloupe$igk_c_genes!=""] <- "IGK-IGL"
vloupe$type <- paste0(vloupe$lourde,"-",vloupe$legere)

vloupe$igh_c_genes <- ifelse(vloupe$igh_c_genes!="", vloupe$igh_c_genes, "NA")
vloupe$V_lourde <- ifelse(vloupe$igh_v_genes!="", vloupe$igh_v_genes, "NA")
vloupe$D_lourde <- ifelse(vloupe$igh_d_genes!="", vloupe$igh_d_genes, "NA")
vloupe$J_lourde <- ifelse(vloupe$igh_j_genes!="", vloupe$igh_j_genes, "NA")

vloupe$V_legere <- ifelse(vloupe$igk_v_genes!="", vloupe$igk_v_genes, ifelse(vloupe$igl_v_genes!="", vloupe$igl_v_genes, "NA"))
vloupe$J_legere <- ifelse(vloupe$igk_j_genes!="", vloupe$igk_j_genes, ifelse(vloupe$igl_j_genes!="", vloupe$igl_j_genes, "NA"))
vloupe$heavy <-str_sub(vloupe$igh_c_genes, 1, 4)

### sauvegarde des objets
Isotype <- count("igh_c_genes")
V <- count(c("V_lourde","V_legere"))
D <- count("D_lourde")
J <- count(c("J_lourde","J_legere"))
Type <- count("type")
Light <- count("legere")
Heavy <- count("heavy")
Clonotype <- c(list(vloupe$clonotype_id[1:5]), list(vloupe$frequency[1:5]))

save(Clonotype, Type, Isotype, V, D, J, vloupe, Heavy, file=paste0("/home/boris/Documents/lipinskib/flinovo/result/",patient,"/R/VDJ.RData"))





### élaboration des figures
fig <- plot_ly(x = names, y = number, name = "Clonotype", type = "bar") %>% layout(xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))
Type <- plot_ly(x = names(table(vloupe$type)), y = table(vloupe$type),name = "Clonotype", type = "bar")
Clonotype <- plot_ly(x = vloupe$clonotype_id[1:5], y = vloupe$frequency[1:5], name = "Clonotype", type = "bar", 
                     text = paste0("Proportion : ", round(vloupe$proportion[1:5],3)*100, "% \nType : ", vloupe$type[1:5], "\nIsotype : ", vloupe$igh_c_genes[1:5], 
                                   "\nHeavy : ", vloupe$V_lourde[1:5], " / ", vloupe$D_lourde[1:5], " / ", vloupe$J_lourde[1:5], "\nLight : ", vloupe$V_legere[1:5], " / ", vloupe$J_legere[1:5]))







for (patient in c("FL140304","FL12C1888","FL09C1164","FL08G0293")) {
  load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))
  
  vloupe <- as.data.frame(read.table( paste0("/home/boris/Documents/lipinskib/boris/VDJ/", patient,"_vloupe-clonotypes.csv"), sep = ",", header = T))
  vloupe$clonotype_id <- as.numeric(str_replace(vloupe$clonotype_id,"clonotype",""))
  
  vloupe$lourde <- ifelse(vloupe$igh_c_genes!="", "IGH", "NA")
  vloupe$legere[vloupe$igl_c_genes=="" & vloupe$igk_c_genes==""] <- "NA"
  vloupe$legere[vloupe$igl_c_genes=="" & vloupe$igk_c_genes!=""] <- "IGK"
  vloupe$legere[vloupe$igl_c_genes!="" & vloupe$igk_c_genes==""] <- "IGL"
  vloupe$legere[vloupe$igl_c_genes!="" & vloupe$igk_c_genes!=""] <- "IGK-IGL"
  vloupe$type <- paste0(vloupe$lourde,"-",vloupe$legere)
  
  vloupe$igh_c_genes <- ifelse(vloupe$igh_c_genes!="", vloupe$igh_c_genes, "NA")
  vloupe$V_lourde <- ifelse(vloupe$igh_v_genes!="", vloupe$igh_v_genes, "NA")
  vloupe$D_lourde <- ifelse(vloupe$igh_d_genes!="", vloupe$igh_d_genes, "NA")
  vloupe$J_lourde <- ifelse(vloupe$igh_j_genes!="", vloupe$igh_j_genes, "NA")
  
  vloupe$V_legere <- ifelse(vloupe$igk_v_genes!="", vloupe$igk_v_genes, ifelse(vloupe$igl_v_genes!="", vloupe$igl_v_genes, "NA"))
  vloupe$J_legere <- ifelse(vloupe$igk_j_genes!="", vloupe$igk_j_genes, ifelse(vloupe$igl_j_genes!="", vloupe$igl_j_genes, "NA"))
  vloupe$heavy <-str_sub(vloupe$igh_c_genes, 1, 4)

  all@tools$Heavy <- count("heavy")
  save(all, file = paste0("/home/boris/Documents/analyse/singlet_",patient,".RData"))
}
