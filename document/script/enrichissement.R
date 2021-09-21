## -- TEST : DEenrichRPlot -- ##
ER <- list()

for (patient in siege) {
  load(file = paste0("/home/boris/Bureau/scShiny/www/", patient,"/", patient,".RData"))
  ER[[patient]]$KEGG$pos <- singlet@tools$KEGG$pos[1]
  ER[[patient]]$GO_Biological$pos <- singlet@tools$GO_Biological$pos[1]
  ER[[patient]]$GO_Cellular$pos <- singlet@tools$GO_Cellular$pos[1]
  ER[[patient]]$GO_Molecular$pos <- singlet@tools$GO_Molecular$pos[1]
    
  ER[[patient]]$KEGG$neg <- singlet@tools$KEGG$neg[1]
  ER[[patient]]$GO_Biological$neg <- singlet@tools$GO_Biological$neg[1]
  ER[[patient]]$GO_Cellular$neg <- singlet@tools$GO_Cellular$neg[1]
  ER[[patient]]$GO_Molecular$neg <- singlet@tools$GO_Molecular$neg[1]
  
  rm(singlet) ; gc() ; gc() ; gc()
}

df <- tibble(Groupe = rep(c("Complete","Partiel"), each = 8),Level = rep(c("UP","DOWN"),8), Fonction = "1", Enrichissement = rep(c(rep("GO : Biological Process",2),rep("GO : Molecular Function",2),rep("GO : Cellular Component",2),rep("KEGG",2)),2))

df$Fonction[1] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$pos[1], ER[["FL12C1888"]]$GO_Biological$pos[1], ER[["FL08G0293"]]$GO_Biological$pos[1])) ))>0,paste(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$pos[1], ER[["FL12C1888"]]$GO_Biological$pos[1], ER[["FL08G0293"]]$GO_Biological$pos[1]))), collapse = "\n"), "NA")
df$Fonction[2] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$neg[1], ER[["FL12C1888"]]$GO_Biological$neg[1], ER[["FL08G0293"]]$GO_Biological$neg[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$neg[1], ER[["FL12C1888"]]$GO_Biological$neg[1], ER[["FL08G0293"]]$GO_Biological$neg[1]))), collapse = "\n"), "NA")
df$Fonction[3] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$pos[1], ER[["FL12C1888"]]$GO_Molecular$pos[1], ER[["FL08G0293"]]$GO_Molecular$pos[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$pos[1], ER[["FL12C1888"]]$GO_Molecular$pos[1], ER[["FL08G0293"]]$GO_Molecular$pos[1]))), collapse = "\n"), "NA")
df$Fonction[4] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$neg[1], ER[["FL12C1888"]]$GO_Molecular$neg[1], ER[["FL08G0293"]]$GO_Molecular$neg[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$neg[1], ER[["FL12C1888"]]$GO_Molecular$neg[1], ER[["FL08G0293"]]$GO_Molecular$neg[1]))), collapse = "\n"), "NA")
df$Fonction[5] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$pos[1], ER[["FL12C1888"]]$GO_Cellular$pos[1], ER[["FL08G0293"]]$GO_Cellular$pos[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$pos[1], ER[["FL12C1888"]]$GO_Cellular$pos[1], ER[["FL08G0293"]]$GO_Cellular$pos[1]))), collapse = "\n"), "NA")
df$Fonction[6] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$neg[1], ER[["FL12C1888"]]$GO_Cellular$neg[1], ER[["FL08G0293"]]$GO_Cellular$neg[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$neg[1], ER[["FL12C1888"]]$GO_Cellular$neg[1], ER[["FL08G0293"]]$GO_Cellular$neg[1]))), collapse = "\n"), "NA")
df$Fonction[7] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$KEGG$pos[1], ER[["FL12C1888"]]$KEGG$pos[1], ER[["FL08G0293"]]$KEGG$pos[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL140304"]]$KEGG$pos[1], ER[["FL12C1888"]]$KEGG$pos[1], ER[["FL08G0293"]]$KEGG$pos[1]))), collapse = "\n"), "NA")
df$Fonction[8] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL140304"]]$KEGG$neg[1], ER[["FL12C1888"]]$KEGG$neg[1], ER[["FL08G0293"]]$KEGG$neg[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL140304"]]$KEGG$neg[1], ER[["FL12C1888"]]$KEGG$neg[1], ER[["FL08G0293"]]$KEGG$neg[1]))), collapse = "\n"), "NA")

df$Fonction[9] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$pos[1], ER[["FL02G095"]]$GO_Biological$pos[1], ER[["FL05G0330"]]$GO_Biological$pos[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$pos[1], ER[["FL02G095"]]$GO_Biological$pos[1], ER[["FL05G0330"]]$GO_Biological$pos[1]))), collapse = "\n"), "NA")
df$Fonction[10] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$neg[1], ER[["FL02G095"]]$GO_Biological$neg[1], ER[["FL05G0330"]]$GO_Biological$neg[1])) ))>0,paste(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$neg[1], ER[["FL02G095"]]$GO_Biological$neg[1], ER[["FL05G0330"]]$GO_Biological$neg[1]))), collapse = "\n"), "NA")
df$Fonction[11] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$pos[1], ER[["FL02G095"]]$GO_Molecular$pos[1], ER[["FL05G0330"]]$GO_Molecular$pos[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$pos[1], ER[["FL02G095"]]$GO_Molecular$pos[1], ER[["FL05G0330"]]$GO_Molecular$pos[1]))), collapse = "\n"), "NA")
df$Fonction[12] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$neg[1], ER[["FL02G095"]]$GO_Molecular$neg[1], ER[["FL05G0330"]]$GO_Molecular$neg[1])) ))>0,paste(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$neg[1], ER[["FL02G095"]]$GO_Molecular$neg[1], ER[["FL05G0330"]]$GO_Molecular$neg[1]))), collapse = "\n"), "NA")
df$Fonction[13] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$pos[1], ER[["FL02G095"]]$GO_Cellular$pos[1], ER[["FL05G0330"]]$GO_Cellular$pos[1])) ))>0,paste(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$pos[1], ER[["FL02G095"]]$GO_Cellular$pos[1], ER[["FL05G0330"]]$GO_Cellular$pos[1]))), collapse = "\n"), "NA")
df$Fonction[14] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$neg[1], ER[["FL02G095"]]$GO_Cellular$neg[1], ER[["FL05G0330"]]$GO_Cellular$neg[1])) ))>0,paste(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$neg[1], ER[["FL02G095"]]$GO_Cellular$neg[1], ER[["FL05G0330"]]$GO_Cellular$neg[1]))), collapse = "\n"), "NA")
df$Fonction[15] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$pos[1], ER[["FL02G095"]]$KEGG$pos[1], ER[["FL05G0330"]]$KEGG$pos[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$pos[1], ER[["FL02G095"]]$KEGG$pos[1], ER[["FL05G0330"]]$KEGG$pos[1]))), collapse = "\n"), "NA")
df$Fonction[16] <- ifelse(length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$neg[1], ER[["FL02G095"]]$KEGG$neg[1], ER[["FL05G0330"]]$KEGG$neg[1])) ))>0, paste(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$neg[1], ER[["FL02G095"]]$KEGG$neg[1], ER[["FL05G0330"]]$KEGG$neg[1]))), collapse = "\n"), "NA")


df$Nombre[1] <- length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$pos[1], ER[["FL12C1888"]]$GO_Biological$pos[1], ER[["FL08G0293"]]$GO_Biological$pos[1])) ))
df$Nombre[2] <--length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Biological$neg[1], ER[["FL12C1888"]]$GO_Biological$neg[1], ER[["FL08G0293"]]$GO_Biological$neg[1])) ))
df$Nombre[3] <- length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$pos[1], ER[["FL12C1888"]]$GO_Molecular$pos[1], ER[["FL08G0293"]]$GO_Molecular$pos[1])) ))
df$Nombre[4] <--length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Molecular$neg[1], ER[["FL12C1888"]]$GO_Molecular$neg[1], ER[["FL08G0293"]]$GO_Molecular$neg[1])) ))
df$Nombre[5] <- length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$pos[1], ER[["FL12C1888"]]$GO_Cellular$pos[1], ER[["FL08G0293"]]$GO_Cellular$pos[1])) ))
df$Nombre[6] <--length(unlist(Reduce(intersect, list(ER[["FL140304"]]$GO_Cellular$neg[1], ER[["FL12C1888"]]$GO_Cellular$neg[1], ER[["FL08G0293"]]$GO_Cellular$neg[1])) ))
df$Nombre[7] <- length(unlist(Reduce(intersect, list(ER[["FL140304"]]$KEGG$pos[1], ER[["FL12C1888"]]$KEGG$pos[1], ER[["FL08G0293"]]$KEGG$pos[1])) ))
df$Nombre[8] <--length(unlist(Reduce(intersect, list(ER[["FL140304"]]$KEGG$neg[1], ER[["FL12C1888"]]$KEGG$neg[1], ER[["FL08G0293"]]$KEGG$neg[1])) ))

df$Nombre[9] <- length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$pos[1], ER[["FL02G095"]]$GO_Biological$pos[1], ER[["FL05G0330"]]$GO_Biological$pos[1])) ))
df$Nombre[10] <--length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Biological$neg[1], ER[["FL02G095"]]$GO_Biological$neg[1], ER[["FL05G0330"]]$GO_Biological$neg[1])) ))
df$Nombre[11] <- length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$pos[1], ER[["FL02G095"]]$GO_Molecular$pos[1], ER[["FL05G0330"]]$GO_Molecular$pos[1])) ))
df$Nombre[12] <--length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Molecular$neg[1], ER[["FL02G095"]]$GO_Molecular$neg[1], ER[["FL05G0330"]]$GO_Molecular$neg[1])) ))
df$Nombre[13] <- length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$pos[1], ER[["FL02G095"]]$GO_Cellular$pos[1], ER[["FL05G0330"]]$GO_Cellular$pos[1])) ))
df$Nombre[14] <--length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$GO_Cellular$neg[1], ER[["FL02G095"]]$GO_Cellular$neg[1], ER[["FL05G0330"]]$GO_Cellular$neg[1])) ))
df$Nombre[15] <- length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$pos[1], ER[["FL02G095"]]$KEGG$pos[1], ER[["FL05G0330"]]$KEGG$pos[1])) ))
df$Nombre[16] <--length(unlist(Reduce(intersect, list(ER[["FL09C1164"]]$KEGG$neg[1], ER[["FL02G095"]]$KEGG$neg[1], ER[["FL05G0330"]]$KEGG$neg[1])) ))


plot_ly(df, x = df$Enrichissement[1:8], y = df$Nombre[1:8], name = "Complete", type = "bar", hovertemplate = paste("Groupe :", df$Groupe[1:8],"\nLevel :", df$Level[1:8],"\nEnrichissement :", df$Enrichissement[1:8],"\nNombre :", df$Nombre[1:8],"\nFonction :\n", df$Fonction[1:8] )) %>% 
  add_trace(df, x = df$Enrichissement[9:16], y = df$Nombre[9:16], name = 'Partiel',hovertemplate = paste("Groupe :", df$Groupe[9:16],"\nLevel :", df$Level[9:16],"\nEnrichissement :", df$Enrichissement[9:16],"\nNombre :", df$Nombre[9:16],"\nFonction :\n", df$Fonction[9:16] ))

df2 <- tibble(Groupe = c("Complete", "Complete", "Complete", "Complete", "Partiel", "Partiel", "Partiel", "Partiel"), Enrichissement = c("GO : Biological Process", "GO : Molecular Function", "GO : Cellular Component", "KEGG", "GO : Biological Process", "GO : Molecular Function", "GO : Cellular Component", "KEGG"),
             UP = c(df$Fonction[1], df$Fonction[3],df$Fonction[5],df$Fonction[7],df$Fonction[9],df$Fonction[11],df$Fonction[13],df$Fonction[15]), DOWN = c(df$Fonction[2], df$Fonction[4],df$Fonction[6],df$Fonction[8],df$Fonction[10],df$Fonction[12],df$Fonction[14],df$Fonction[16]))

df2 %>% datatable(rownames = FALSE)

