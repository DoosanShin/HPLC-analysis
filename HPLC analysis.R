#Packages installation
install.packages("readxl")
installed.packages("stringr")
install.packages('dplyr')
install.packages('readr')
install.packages('ggplot2')
install.packages('patchwork')
install.packages('gridExtra')
install.packages("ggforce")
library(ggforce)
library(readxl)
library(stringr)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(gridExtra)

setwd("C:\\Users\\creac\\Desktop\\KFB_UGT_KO_leaf") # folder with HPLC excel files
filenames <- dir()
allData <- list()

#Retrive all Exel file from the directory folder 
for(i in 1:length(filenames)){
  fullSheet <- read_excel(filenames[i], sheet="Integration")
  startRow <- which(fullSheet[ , 1] == "Integration Results") + 4
  endRow <- nrow(fullSheet) - 1
  intResults <- fullSheet[startRow:endRow, ]
  colnames(intResults) <- paste(fullSheet[startRow - 3, ])
  allData[[i]] <- intResults
}

names(allData) <- filenames

# need to find sample with most rows - account for all possible metabolite peaks
peakCount <- vector(mode="numeric", length=length(allData))
for(i in 1:length(allData)){
  peakCount[i] <- nrow(allData[[i]])
}
allPeakNames <- allData[[which.max(peakCount)]]$'Peak Name'
allPeakNames[is.na(allPeakNames)] <- "Unknown"

# create new data frames to store height and area data
areaData <- data.frame(matrix(ncol=length(allPeakNames), nrow=length(allData)))

colnames(areaData) <- allPeakNames

# pull out area data from all samples 
for(i in 1:length(allData)){
  sampleName <- names(allData[i])
  peakNames <- allData[[i]]$`Peak Name`
  peakAreas <- allData[[i]]$Area
  for(peak in peakNames){
    areaData[[peak]][i] <- peakAreas[which(peakNames == peak)]
  }
  rownames(areaData)[i] <- sampleName
}

areaData[areaData=="n.a."] <- NA



#Data analysis

areaData <- tibble::rownames_to_column(areaData, "samples")
areaData2 <- select(areaData, -contains("Unknown"))

#Calculate average of each genotypes
Component_num <- as.numeric(sum(str_count(allPeakNames, "Component")))  # The number of detected component
Average <- data.frame(matrix(ncol=Component_num, nrow=6)) #nrow=length(genotype number)
colnames(Average) <- grep("Component ", allPeakNames,value=TRUE)
row.names(Average) <- c("VC", "CR2-2","CR2-5", "CR4-1","CR4-2", "UC8-2" )

names(Average)<-str_replace_all(names(Average), c(" " = "_"))
names(areaData2)<-str_replace_all(names(areaData2), c(" " = "_"))

for(i in 1:Component_num) {
  Component_avr_VC <- areaData2 %>%
    filter(grepl('VC', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  avr <- mean(Component_avr_VC[,1], na.rm = TRUE)
  Average[1, i] <- avr
  
  Component_avr_CR2_2 <- areaData2 %>%
    filter(grepl('CR2-2', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  avr <- mean(Component_avr_CR2_2[,1], na.rm = TRUE)
  Average[2, i] <- avr
  
  Component_avr_CR2_5 <- areaData2 %>%
    filter(grepl('CR2-5', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  avr <- mean(Component_avr_CR2_5[,1], na.rm = TRUE)
  Average[3, i] <- avr
  
  Component_avr_CR4_1 <- areaData2 %>%
    filter(grepl('CR4-1', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  avr <- mean(Component_avr_CR4_1[,1], na.rm = TRUE)
  Average[4, i] <- avr
  
  Component_avr_CR4_2 <- areaData2 %>%
    filter(grepl('CR4-2', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  avr <- mean(Component_avr_CR4_2[,1], na.rm = TRUE)
  Average[5, i] <- avr
  
  Component_avr_UC8_2 <- areaData2 %>%
    filter(grepl('UC8-2', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  avr <- mean(Component_avr_UC8_2[,1], na.rm = TRUE)
  Average[6, i] <- avr
}

#Calculate standard error of each genotypes
Component_num <- as.numeric(sum(str_count(allPeakNames, "Component")))  # The number of detected component
STE <- data.frame(matrix(ncol=Component_num, nrow=6)) #nrow=length(genotype number)
colnames(STE) <- grep("Component ", allPeakNames,value=TRUE)
row.names(STE) <- c("VC", "CR2-2","CR2-5", "CR4-1","CR4-2", "UC8-2" )

names(STE)<-str_replace_all(names(STE), c(" " = "_"))
names(areaData2)<-str_replace_all(names(areaData2), c(" " = "_"))

stderr <- function(x, na.rm=TRUE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

for(i in 1:Component_num) {
  Component_avr_VC <- areaData2 %>%
    filter(grepl('VC', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  SE <- stderr(Component_avr_VC[,1])
  STE[1, i] <- SE
  
  Component_avr_CR2_2 <- areaData2 %>%
    filter(grepl('CR2-2', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  SE <- stderr(Component_avr_CR2_2[,1])
  STE[2, i] <- SE
  
  Component_avr_CR2_5 <- areaData2 %>%
    filter(grepl('CR2-5', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  SE <- stderr(Component_avr_CR2_5[,1])
  STE[3, i] <- SE
  
  Component_avr_CR4_1 <- areaData2 %>%
    filter(grepl('CR4-1', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  SE <- stderr(Component_avr_CR4_1[,1])
  STE[4, i] <- SE
  
  Component_avr_CR4_2 <- areaData2 %>%
    filter(grepl('CR4-2', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  SE <- stderr(Component_avr_CR4_2[,1])
  STE[5, i] <- SE
  
  Component_avr_UC8_2 <- areaData2 %>%
    filter(grepl('UC8-2', areaData2$samples)) %>%
    select(1+i) %>%
    mutate_all(as.numeric) 
  SE <- stderr(Component_avr_UC8_2[,1])
  STE[6, i] <- SE

}



#Bar graph drawing
Average <- tibble::rownames_to_column(Average, "samples")
#Turn your 'treatment' column into a character vector
Average$samples <- as.character(Average$samples)
#Then turn it back into a factor with the levels in the correct order <- this is for x-axis order
Average$samples <- factor(Average$samples, levels=unique(Average$samples))

# To remove specific samples
Average <- Average[-4, ]
STE <- STE[-4, ]

#Plotting test
p1 <- ggplot(Average) +
  geom_bar( aes(x=Average$samples, y=Average[,2]), stat="identity", fill="grey2", alpha=0.9, width=0.6) +
  geom_errorbar( aes(x=Average$samples, ymin=Average[,2]-STE[,1], ymax=Average[,2]+STE[,1]), width=0.2, 
                 colour="black", alpha=0.9, size=0.9) + 
  ggtitle(colnames(Average[2])) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Peak area(mAU*min)")

STE[,1]
Average[,2]
p2 <- ggplot(Average) +
  geom_bar( aes(x=Average$samples, y=Average[,3]), stat="identity", fill="grey2", alpha=0.9, width=0.6) +
  geom_errorbar( aes(x=Average$samples, ymin=Average[,3]-STE[,2], ymax=Average[,3]+STE[,2]), width=0.2, 
                 colour="black", alpha=0.9, size=0.9) + 
  ggtitle(colnames(Average[3])) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Peak area(mAU*min)")

p1 + p2


#Generation mutiple plot
plot_list = list()
for(i in 1:Component_num){
  graphy <- print(ggplot(Average, aes_string(x=Average$samples, y=Average[ ,i+1])) +
                    geom_bar(stat="identity", fill="grey2", alpha=0.9, width=0.6) +
                    geom_errorbar( aes_string(x=Average$samples, ymin=Average[ ,i+1]-STE[ ,i], ymax=Average[ ,i+1]+STE[ ,i]), width=0.2, 
                                   colour="black", alpha=0.9, size=0.9) + 
                    ggtitle(colnames(Average[i+1])) +
                    theme(axis.title.x=element_blank(),
                          axis.ticks.x=element_blank()) +
                    ylab("Peak area(mAU*min)"))
  plot_list[[i]] <- graphy
}

#Plot paring
paired_plot <- list()
plot_mat=matrix(1:6, nrow = 3, ncol = 2)
paired_plot <- grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],  plot_list[[5]], plot_list[[6]],
             layout_matrix=plot_mat)

ggsave(path = "C:\\Users\\creac\\Desktop\\KFB_UGT_KO_leaf plot", file="whatever.pdf", paired_plot)

#Plot paring looping and saving
paired_plot <- list()
plot_mat=matrix(1:6, nrow = 3, ncol = 2)

for (i in seq(from=1, to=(length(plot_list)), by=6)) {
  paired <- grid.arrange(plot_list[[i]], plot_list[[i+1]], plot_list[[i+2]], plot_list[[i+3]],  plot_list[[i+4]], plot_list[[i+5]],
                              layout_matrix=plot_mat)
  paired_plot[[i]] <- paired
  
  ggsave(path = "C:\\Users\\creac\\Desktop\\KFB_UGT_KO_leaf plot\\Original", 
         file=paste('SlKFB1 and SlUGT78 double KO leaf',  i, ".pdf", sep=""), 
         paired_plot[[i]])
  
}

  
  
  