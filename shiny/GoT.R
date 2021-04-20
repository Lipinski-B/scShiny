library(stringr)
library(plyr)
library(ggplot2)

patient <- "FL12C1888"
setwd(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/"))

got <- function(hotspot, colname, patient){
  BC <- read.table(paste0("/home/boris/Documents/lipinskib/flinovo/data/", patient, "/GOT/barcodes/barcodes.txt"))
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
BCL2_L23L <- got(read.table("BCL2_L23L.out/BCL2_L23L.summTable.txt", header = T), "BCL2_L23L", patient)
BCL2_K22K <- got(read.table("BCL2_K22K.out/BCL2_K22K.summTable.txt", header = T), "BCL2_K22K", patient)
CD79B <- got(read.table("CD79B.out/CD79B.summTable.txt", header = T), "CD79B_Y696H", patient)
EZH2_A682G_1 <- got(read.table("EZH2_A682G_1.out/EZH2_A682G_1.summTable.txt", header = T), "EZH2_A682G_1", patient)
EZH2_A682G_2 <- got(read.table("EZH2_A682G_2.out/EZH2_A682G_2.summTable.txt", header = T), "EZH2_A682G_2", patient)
EZH2_A692V_1 <- got(read.table("EZH2_A692V_1.out/EZH2_A692V_1.summTable.txt", header = T), "EZH2_A692V_1", patient)
EZH2_A692V_2 <- got(read.table("EZH2_A692V_2.out/EZH2_A692V_2.summTable.txt", header = T), "EZH2_A692V_2", patient)
EZH2_Y646C <- got(read.table("EZH2_Y646C.out/EZH2_Y646C.summTable.txt", header = T), "EZH2_Y646C", patient)
EZH2_Y646F <- got(read.table("EZH2_Y646F.out/EZH2_Y646F.summTable.txt", header = T), "EZH2_Y646F", patient)
EZH2_Y646H <- got(read.table("EZH2_Y646H.out/EZH2_Y646H.summTable.txt", header = T), "EZH2_Y646H", patient)
EZH2_Y646N <- got(read.table("EZH2_Y646N.out/EZH2_Y646N.summTable.txt", header = T), "EZH2_Y646N", patient)
EZH2_Y646S <- got(read.table("EZH2_Y646S.out/EZH2_Y646S.summTable.txt", header = T), "EZH2_Y646S", patient)

GOT<- cbind(BCL2_L23L[["tomerge"]], BCL2_K22K[["tomerge"]], CD79B[["tomerge"]], EZH2_A682G_1[["tomerge"]], EZH2_A682G_2[["tomerge"]], EZH2_A692V_1[["tomerge"]], EZH2_A692V_2[["tomerge"]], EZH2_Y646C[["tomerge"]],EZH2_Y646F[["tomerge"]], EZH2_Y646H[["tomerge"]], EZH2_Y646N[["tomerge"]], EZH2_Y646S[["tomerge"]])
rownames(GOT) <-paste0(rownames(GOT),"-1")
#save(GOT, file=paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/", patient, "_GoT.Rdata"))

##################################################################################################################################################################
######### Violin plot #########
###############################
vlnp <- function(hotspot,mutation){
  PAT <- c("FL12C1888", "FL140304")
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
      hotspot=rep(mutation, length(c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum))))
    )

    e <- rbind(e, m)
  }
  return(e)
}

#patien=rep(patient,length(c(as.numeric(figure1[,2]),as.numeric(figure2[,2]),as.numeric(figure3$sum)))),
cd79B <- vlnp("CD79B","CD79B")
bcl2 <- vlnp("BCL2", "BCL2_K22K")
ezh2 <- vlnp("EZH2","EZH2_A682G_1")

a <- rbind(cd79B, bcl2, ezh2)

ggplot(ezh2, aes(x=name,y=sum)) +
  geom_violin(trim = FALSE) +
  annotate("rect", xmin=1.5, xmax=Inf, ymin=0, ymax=Inf, alpha=0.2, fill="red")+
  annotate("rect", xmin=0.41, xmax=1.5, ymin=0, ymax=Inf, alpha=0.2, fill="blue")+
  geom_violin(trim = FALSE) + geom_boxplot(width=.03, outlier.size=0, fill="grey60") +
  stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2) +
  xlab("") + ylab("UMI counts per cell") + scale_y_continuous(breaks = seq(0,64,4), trans = "log2") +
  labs(title="EZH2 GoT: ") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text = element_text(colour = "black"))+
  theme(
    panel.background = element_rect(colour = "grey70", size = 2, linetype = "solid"),
    panel.grid.major = element_line(size = 0.9, linetype = 'solid', colour = "white"), 
    panel.grid.minor = element_line(size = 0.5, linetype = 'solid', colour = "white"))+
  facet_grid(hotspot ~ patien)


##################################################################################################################################################################
######### Histograme 1 #########
################################
VLNP<- rbind(BCL2_L23L[["sum"]], BCL2_K22K[["sum"]], CD79B_Y696H[["sum"]], EZH2_A682G_1[["sum"]], 
             EZH2_A682G_2[["sum"]], EZH2_A692V_1[["sum"]], EZH2_A692V_2[["sum"]], EZH2_Y646C[["sum"]], 
             EZH2_Y646F[["sum"]], EZH2_Y646H[["sum"]], EZH2_Y646N[["sum"]], EZH2_Y646S[["sum"]])

ggplot(VLNP, aes(x=hotspot,y=sum)) + 
  geom_violin(trim = FALSE) +
  geom_boxplot(width=.1, outlier.size=0, fill="grey50") +
  stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=4) +
  xlab("Hotspot") +
  ylab("Number of UMI/Cell") +
  theme(axis.text = element_text(colour = "black"))# Violin plot basicp



df2 <- data.frame(hotspot=rep(c("EZH2_A682G_1", "EZH2_A682G_2", "EZH2_A692V_1",  "EZH2_A692V_2", "EZH2_Y646C", "EZH2_Y646F", "EZH2_Y646N", "EZH2_Y646S"), each=4), 
                  condition=rep(c("AMB", "MUT", "None", "WT"),8),
                  frequence=c(round((table(GOT[,4])/length(GOT[,1]))*100,2),
                              round((table(GOT[,5])/length(GOT[,1]))*100,2), 
                              round((table(GOT[,6])/length(GOT[,1]))*100,2), 
                              round((table(GOT[,7])/length(GOT[,1]))*100,2), 
                              round((table(GOT[,8])/length(GOT[,1]))*100,2), 
                              round((table(GOT[,9])/length(GOT[,1]))*100,2),
                              round((table(GOT[,11])/length(GOT[,1]))*100,2),
                              round((table(GOT[,12])/length(GOT[,1]))*100,2)))

p <- ggplot(data=df2, aes(x=hotspot, y=frequence, fill=condition)) + geom_bar(stat="identity", color="black", position=position_dodge()) + theme_minimal()
p + scale_fill_manual(values=c('#999999','#E69F00'))
p + scale_fill_brewer(palette="Blues")

#Histograme 2
df3 <- data.frame(hotspot=rep(c("BCL2_L23L", "BCL2_K22K", "CD79B_Y696H",  "EZH2_Y646H"), each=3),
                  condition=rep(c("AMB", "None", "WT"),4),
                  frequence=c(
                    round((table(GOT[,1])/length(GOT[,1]))*100,2),
                    round((table(GOT[,2])/length(GOT[,1]))*100,2),
                    round((table(GOT[,3])/length(GOT[,1]))*100,2),
                    round((table(GOT[,10])/length(GOT[,1]))*100,2)))

p <- ggplot(data=df3, aes(x=hotspot, y=frequence, fill=condition)) + geom_bar(stat="identity", color="black", position=position_dodge())+ theme_minimal()
p + scale_fill_manual(values=c('#999999','#E69F00'))
p + scale_fill_brewer(palette="Blues")
