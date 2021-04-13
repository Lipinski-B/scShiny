library(stringr)
library(plyr)
library(ggplot2)

patient <- "FL140304"
setwd(paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/"))
BC <- read.table(paste0("/home/boris/Documents/lipinskib/flinovo/data/", patient, "/GOT/barcodes/barcodes.txt"))

got <- function(hotspot, colname){
  result <- cbind(hotspot[,17],hotspot[,18],hotspot[,18])
  rownames(result) <- unique(unlist(str_split(hotspot[,1],";")))
  colnames(result) <- c("WT","MUT","BCL2_L23L")
  
  for (i in 1:length(rownames(result))){
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
  
  final <- as.data.frame(result[,3])
  count(final)
  
  a <- as.data.frame(setdiff(BC$V1, rownames(final)))
  rownames(a) <- a$`setdiff(BC$V1, rownames(final))`
  colnames(a) <- c("tomerge")
  colnames(final) <- c("tomerge")
  a[,1] <- 'None'
  
  final <- rbind(final,a)
  colnames(final) <- colname
  
  return(final)
}
BCL2_L23L <- got(read.table("BCL2_L23L.out/BCL2_L23L.summTable.txt", header = T), "BCL2_L23L")
BCL2_K22K <- got(read.table("BCL2_K22K.out/BCL2_K22K.summTable.txt", header = T), "BCL2_K22K")
CD79B_Y696H <- got(read.table("CD79B_Y696H.out/CD79B_Y696H.summTable.txt", header = T), "CD79B_Y696H")
EZH2_A682G_1 <- got(read.table("EZH2_A682G_1.out/EZH2_A682G_1.summTable.txt", header = T), "EZH2_A682G_1")
EZH2_A682G_2 <- got(read.table("EZH2_A682G_2.out/EZH2_A682G_2.summTable.txt", header = T), "EZH2_A682G_2")
EZH2_A692V_1 <- got(read.table("EZH2_A692V_1.out/EZH2_A692V_1.summTable.txt", header = T), "EZH2_A692V_1")
EZH2_A692V_2 <- got(read.table("EZH2_A692V_2.out/EZH2_A692V_2.summTable.txt", header = T), "EZH2_A692V_2")
EZH2_Y646C <- got(read.table("EZH2_Y646C.out/EZH2_Y646C.summTable.txt", header = T), "EZH2_Y646C")
EZH2_Y646F <- got(read.table("EZH2_Y646F.out/EZH2_Y646F.summTable.txt", header = T), "EZH2_Y646F")
EZH2_Y646H <- got(read.table("EZH2_Y646H.out/EZH2_Y646H.summTable.txt", header = T), "EZH2_Y646H")
EZH2_Y646N <- got(read.table("EZH2_Y646N.out/EZH2_Y646N.summTable.txt", header = T), "EZH2_Y646N")
EZH2_Y646S <- got(read.table("EZH2_Y646S.out/EZH2_Y646S.summTable.txt", header = T), "EZH2_Y646S")

GOT <- cbind(BCL2_L23L, BCL2_K22K, CD79B_Y696H, EZH2_A682G_1, EZH2_A682G_2, EZH2_A692V_1, EZH2_A692V_2, EZH2_Y646C, EZH2_Y646F, EZH2_Y646H, EZH2_Y646N, EZH2_Y646S)

rownames(GOT) <-paste0(rownames(GOT),"-1")
save(GOT, file=paste0("/home/boris/Documents/lipinskib/flinovo/result/", patient, "/GOT/result/", patient, "_GoT.Rdata"))
dim(GOT)

#Histograme 1
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
