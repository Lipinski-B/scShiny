library(stringr)
library(plyr)
library(ggplot2)

patient <- "FL09C1164"
setwd(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/"))

got <- function(hotspot, colname, patient){
  #BC <- read.table(paste0("/home/boris/Documents/lipinskib/BTG1/data/GOT/barcodes/patients/BC_AC.txt"))
  BC <- read.table(paste0("/home/boris/Documents/lipinskib/flinovo/data/",patient,"/GOT/barecodes/barcodes.txt"))
  result <- cbind(hotspot[,17],hotspot[,18],hotspot[,19], hotspot[,19])
  rownames(result) <- unique(unlist(str_split(hotspot[,1],";")))
  colnames(result) <- c("WT","MUT", colname,"SUM")
  
  for (i in 1:length(rownames(result))){
    result[i,4] <- sum(as.integer(result[i,1]),as.integer(result[i,2]),as.integer(result[i,3]))
    if(result[i,1]=="0" & result[i,2]!="0"){result[i,3]  <- colnames(result)[2]}
    else if(result[i,1]!="0" & result[i,2]=="0"){ result[i,3]  <- colnames(result)[1]} 
    else if(result[i,1]=="0" & result[i,2]=="0"){ result[i,3]  <-  "AMB"} 
    else {
      WT <- as.double(result[i,1]) / (as.double(result[i,1]) +  as.double(result[i,2]))
      MUT <- as.double(result[i,2]) / (as.double(result[i,1]) +  as.double(result[i,2]))
      genotype <- max(WT,MUT)
    
      if(genotype  < 0.8){result[i,3]  <- "AMB"} 
      else {
        if(genotype == WT){result[i,3]  <- colnames(result)[1]
        } else {result[i,3]  <- colnames(result)[2]}
      }
    }
  }
  final <- list()
  final[["tomerge"]] <- data.frame(tomerge=result[,3])
  final[["sum"]] <- data.frame(sum=result[,4], hotspot=rep(colname, length(result[,4])))
  rownames(final[["sum"]])<-NULL
  a <- data.frame(tomerge=setdiff(BC$V1, rownames(final[["tomerge"]])))
  rownames(a) <- a$tomerge
  a[,1] <- NA
  
  final[["tomerge"]] <- rbind(final[["tomerge"]],a)

  return(final)
}
vlnp <- function(hotspot,mutation){
  PAT <- c("FL140304","FL12C1888","FL09C1164","FL08G0293")
  e <- NULL
  for (i in 1:length(PAT)) {
    patient <- PAT[i]
    gt <- got(read.table(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/",mutation,".out/",mutation,".summTable.txt"), header = T), mutation, patient)
    
    figure1 <- read.table(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient,"/R/figure_GoT/", hotspot,"/", hotspot,"_nb_umi_figure1.txt"), header = T)
    figure2 <- read.table(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient,"/R/figure_GoT/", hotspot,"_hotspot/", hotspot,"_hotspot_nb_umi_figure1.txt"), header = T)
    figure3 <- data.frame(name=rep(hotspot,length(gt[["sum"]][,1])), sum=as.numeric(gt[["sum"]][,1]))
    
    m <- data.frame(
      name=c(rep("10X",length(figure1[,1])),rep("10X Targeted",length(figure2[,1])),rep("GoT",length(figure3[,1]))), 
      sum=c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum)),
      patien=rep(paste0(patient,": ",round((length(figure3$name)/length(gt$tomerge[,1]))*100,2)," % genotyping"),length(c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum)))),
      #patien=rep(patient,length(c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum)))),
      #hotspot=rep(mutation, length(c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum)))),
      hotspot=rep(paste0(mutation,": ",round((length(figure3$name)/length(gt$tomerge[,1]))*100,2)," % genotyping"), length(c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum))))
    )
    
    e <- rbind(e, m)
  }
  return(e)
}

##################################################################################################################################################################
######### Flinovo #########
###########################
BCL2_L23L <- got(read.table("BCL2_L23L.out/BCL2_L23L.summTable.txt", header = T), "BCL2_L23L", patient)
BCL2_K22K <- got(read.table("BCL2_K22K.out/BCL2_K22K.summTable.txt", header = T), "BCL2_K22K", patient)
CD79B <- got(read.table("CD79B.out/CD79B.summTable.txt", header = T), "CD79B_Y696H", patient)
#EZH2_A682G_1 <- got(read.table("EZH2_A682G_1.out/EZH2_A682G_1.summTable.txt", header = T), "EZH2_A682G_1", patient)
EZH2_A682G_2 <- got(read.table("EZH2_A682G_2.out/EZH2_A682G_2.summTable.txt", header = T), "EZH2_A682G_2", patient)
EZH2_A692V_1 <- got(read.table("EZH2_A692V_1.out/EZH2_A692V_1.summTable.txt", header = T), "EZH2_A692V_1", patient)
#EZH2_A692V_2 <- got(read.table("EZH2_A692V_2.out/EZH2_A692V_2.summTable.txt", header = T), "EZH2_A692V_2", patient)
EZH2_Y646C <- got(read.table("EZH2_Y646C.out/EZH2_Y646C.summTable.txt", header = T), "EZH2_Y646C", patient)
EZH2_Y646F <- got(read.table("EZH2_Y646F.out/EZH2_Y646F.summTable.txt", header = T), "EZH2_Y646F", patient)
EZH2_Y646H <- got(read.table("EZH2_Y646H.out/EZH2_Y646H.summTable.txt", header = T), "EZH2_Y646H", patient)
EZH2_Y646N <- got(read.table("EZH2_Y646N.out/EZH2_Y646N.summTable.txt", header = T), "EZH2_Y646N", patient)
EZH2_Y646S <- got(read.table("EZH2_Y646S.out/EZH2_Y646S.summTable.txt", header = T), "EZH2_Y646S", patient)

GOT<- cbind(BCL2_L23L[["tomerge"]], BCL2_K22K[["tomerge"]], CD79B[["tomerge"]], EZH2_A682G_2[["tomerge"]], EZH2_A692V_1[["tomerge"]], EZH2_Y646C[["tomerge"]],EZH2_Y646F[["tomerge"]], EZH2_Y646H[["tomerge"]], EZH2_Y646N[["tomerge"]], EZH2_Y646S[["tomerge"]])
rownames(GOT) <-paste0(rownames(GOT),"-1")
colnames(GOT) <- c( "BCL2_L23L","BCL2_K22K", "CD79B","EZH2_A682G","EZH2_A692V","EZH2_Y646C","EZH2_Y646F","EZH2_Y646H","EZH2_Y646N","EZH2_Y646S")
save(GOT, file=paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/", patient, "_GoT.Rdata"))

##################################################################################################################################################################
######### Violin plot #########
###############################
cd79B <- vlnp("CD79B","CD79B")
bcl2 <- vlnp("BCL2", "BCL2_K22K")
ezh2 <- vlnp("EZH2","EZH2_A682G_1")

a <- rbind(cd79B, bcl2, ezh2)

ggplot(bcl2, aes(x=name,y=sum)) +
  geom_violin(trim = FALSE) +
  annotate("rect", xmin=1.5, xmax=Inf, ymin=0, ymax=Inf, alpha=0.2, fill="red") + annotate("rect", xmin=0.41, xmax=1.5, ymin=0, ymax=Inf, alpha=0.2, fill="blue")+
  geom_violin(trim = FALSE) + geom_boxplot(width=.03, outlier.size=0, fill="grey60") + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2) +
  xlab("") + ylab("UMI counts per cell") + scale_y_continuous(breaks = seq(0,250,20), trans = "log2") + labs(title="HOTSPOT GoT: ") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text = element_text(colour = "black"))+
  theme(panel.background = element_rect(colour = "grey70", size = 2, linetype = "solid"),panel.grid.major = element_line(size = 0.9, linetype = 'solid', colour = "white"),panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "white"))+
  facet_grid(~hotspot)


##################################################################################################################################################################
######### Histograme #########
##############################
histo <- data.frame(hotspot=c(rep("BCL2_L23L",length(table(useNA = "always", BCL2_L23L[[1]]))),rep("BCL2_K22K",length(table(useNA = "always", BCL2_K22K[[1]]))),
                              rep("CD79B_Y196H",length(table(useNA = "always", CD79B[[1]]))),rep("EZH2_A682G",length(table(useNA = "always", EZH2_A682G_2[[1]]))),
                              rep("EZH2_A692V",length(table(useNA = "always", EZH2_A692V_1[[1]]))),
                              rep("EZH2_Y646C",length(table(useNA = "always", EZH2_Y646C[[1]]))),
                              rep("EZH2_Y646F",length(table(useNA = "always", EZH2_Y646F[[1]]))),rep("EZH2_Y646H",length(table(useNA = "always", EZH2_Y646H[[1]]))),
                              rep("EZH2_Y646N",length(table(useNA = "always", EZH2_Y646N[[1]]))),rep("EZH2_Y646S",length(table(useNA = "always", EZH2_Y646S[[1]])))),
                    Génotype=c(names(table(useNA = "always", BCL2_L23L[[1]])),names(table(useNA = "always", BCL2_K22K[[1]])),names(table(useNA = "always", CD79B[[1]])),names(table(useNA = "always", EZH2_A682G_2[[1]])),
                                names(table(useNA = "always", EZH2_A692V_1[[1]])), names(table(useNA = "always", EZH2_Y646C[[1]])),
                                names(table(useNA = "always", EZH2_Y646F[[1]])),names(table(useNA = "always", EZH2_Y646H[[1]])),names(table(useNA = "always", EZH2_Y646N[[1]])),names(table(useNA = "always", EZH2_Y646S[[1]]))), 
                    frequence=c(round((table(useNA = "always", GOT[,1])/length(GOT[,1]))*100,2),round((table(useNA = "always", GOT[,2])/length(GOT[,1]))*100,2),round((table(useNA = "always", GOT[,3])/length(GOT[,1]))*100,2),round((table(useNA = "always", GOT[,4])/length(GOT[,1]))*100,2),round((table(useNA = "always", GOT[,5])/length(GOT[,1]))*100,2),
                                round((table(useNA = "always", GOT[,6])/length(GOT[,1]))*100,2),round((table(useNA = "always", GOT[,7])/length(GOT[,1]))*100,2),round((table(useNA = "always", GOT[,8])/length(GOT[,1]))*100,2),round((table(useNA = "always", GOT[,9])/length(GOT[,1]))*100,2),round((table(useNA = "always", GOT[,10])/length(GOT[,1]))*100,2)))

ggplot(data=histo, aes(x=hotspot, y=frequence, fill=Génotype)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + labs(title="Génotypage du transcriptome") +
  xlab("Mutation") + ylab("Fréquence des génotypes") +
  theme_minimal() + scale_fill_manual(values=c('#250000','#E69F00', '#999999'))# + scale_fill_brewer(palette="Blues")
