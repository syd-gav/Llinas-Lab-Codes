#Imhoi Koo and Gabriel W. Rangel

#Note from Sydney M Gavula: make sure that Java JDK is installed and set as an environment variable on your computer and set as new environment variable in settings

library(pacman)
library(rJava)
library(rcdk)
pacman::p_load(readxl, dplyr, ggplot2, xlsx)


#Load the MS-Dial exported Area, reformatted
posFileName <- "raw-pos-alignment.xlsx"
negFileName <- "raw-neg-alignment.xlsx"


#Load in the Amino Acid Reference Chart .xlsx file
fileName <- "./AAreferenceChart.xlsx"

#Load in Peptide List
peptides <- read.table("./Peptides_01.txt", stringsAsFactors = FALSE, col.names = "shortName")

maxLengthPep <- 13

#Set Output File Name
outputFileName <- paste("hem-peptides",maxLengthPep,".xlsx")



calFt <- function(eltHb) {
  setHb <- NULL
  for (i in 1:length(eltHb)-1){
    for (j in 1:maxLengthPep){
      if (i+j > length(eltHb)){
        break
      } else {
        setHb <- c(setHb, paste0(eltHb[i:(i+j)], collapse = ""))
      }
    }
  }
  
  setHb <- data.frame(AAs = unique(setHb), stringsAsFactors = FALSE)
  
  mass <- apply(AAdata, 1, function(x) {
    y <- rcdk::get.formula(x["Molecular Formula"][[1]], charge=0)@mass
    return(cbind(x["1-Letter Symbol"], x["Molecular Formula"], y))
  })
  H2Omass <- rcdk::get.formula("H2O", charge = 0)@mass
  AAMass <- data.frame(AA=mass[1,],
                       formula = (mass[2,]),
                       mass = as.numeric(mass[3,]),
                       stringsAsFactors = FALSE)
  
  for (i in 1:length(setHb$AAs)){
    temp1 <- strsplit(setHb$AAs[i], split="")[[1]]
    for (j in 1:length(temp1)){
      indx <- which(temp1[j]==AAdata$`1-Letter Symbol`)
    }
    temp2 <- sapply(temp1,
                    function(x) {
                      indx <- which(x==AAMass$AA)
                      AAMass[indx,]
                    })
    temp3 <- paste0(temp2["formula",], collapse="")
    temp4 <- sum(unlist(temp2["mass",])) -(length(temp1)-1)*H2Omass
    # temp4 <- sum(unlist(temp2["mass",])) -(stringr::str_length(setHb$AAs[i])-1)*H2Omass
    setHb[i,2] <- temp3
    setHb[i,3] <- temp4
    rm(list=c("temp1", 'temp2', 'temp3', 'temp4'))
  }
  colnames(setHb) <- c("shortName","formula",'mass')
  Hmass <- 1.007276
  
  setHb <- setHb %>% mutate(mzPos = mass+Hmass, mzNeg = mass-Hmass)
  return(setHb)
}

calFtList <- function(list=peptides$shortName) {
  
  setHb <- data.frame(AAs = unique(list), stringsAsFactors = FALSE)
  
  mass <- apply(AAdata, 1, function(x) {
    y <- rcdk::get.formula(x["Molecular Formula"][[1]], charge=0)@mass
    return(cbind(x["1-Letter Symbol"], x["Molecular Formula"], y))
  })
  H2Omass <- rcdk::get.formula("H2O", charge = 0)@mass
  AAMass <- data.frame(AA=mass[1,],
                       formula = (mass[2,]),
                       mass = as.numeric(mass[3,]),
                       stringsAsFactors = FALSE)
  
  for (i in 1:length(setHb$AAs)){
    temp1 <- strsplit(setHb$AAs[i], split="")[[1]]
    for (j in 1:length(temp1)){
      indx <- which(temp1[j]==AAdata$`1-Letter Symbol`)
    }
    temp2 <- sapply(temp1,
                    function(x) {
                      indx <- which(x==AAMass$AA)
                      AAMass[indx,]
                    })
    temp3 <- paste0(temp2["formula",], collapse="")
    temp4 <- sum(unlist(temp2["mass",])) -(length(temp1)-1)*H2Omass
    # temp4 <- sum(unlist(temp2["mass",])) -(stringr::str_length(setHb$AAs[i])-1)*H2Omass
    setHb[i,2] <- temp3
    setHb[i,3] <- temp4
    rm(list=c("temp1", 'temp2', 'temp3', 'temp4'))
  }
  colnames(setHb) <- c("shortName","formula",'mass')
  Hmass <- 1.007276
  
  setHb <- setHb %>% mutate(mzPos = mass+Hmass, mzNeg = mass-Hmass)
  return(setHb)
}

selFt <- function(fileName, setH=setHb, mz='mzNeg') {
  sheetsName <- excel_sheets(fileName)
  
  data <- lapply(sheetsName, function(x) as.data.frame(read_excel(fileName, sheet = x)))
  names(data) <- sheetsName
  rownames(data$sample_info) <- data$sample_info$fileName
  
  rownames(data$feature_info) <- data$feature_info$ID
  rownames(data$feature_matrix) <- data$feature_info$ID
  # data$feature_matrix$ID <- NULL
  # rownames(data$TIC_normalize) <- data$feature_info$ID
  # data$TIC_normalize$ID <- NULL
  # rownames(data$IS_normal) <- data$feature_info$ID
  # data$IS_normal$ID <- NULL
  
  uniqMZ <- unique(setH[,mz])
  ppmThrhd = 15
  selFeature <- NULL
  for (i in 1:length(uniqMZ)){
    temp1 <- paste0(setH$shortName[which(setH[,mz]==uniqMZ[i])], collapse=";")
    temp2 <- which(abs(data$feature_info$`Average Mz`-uniqMZ[i])<ppmThrhd*uniqMZ[i]*10^-6)
    if (sum(temp2)!=0){
      temp3 <- data.frame(unsetldName=temp1,
                          calmz = uniqMZ[i],
                          data$feature_info[temp2,])
      selFeature <- rbind(selFeature, temp3)
    }
  }
  return(selFeature)
}

Hbeta  <- 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
Halpha <- 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR'

eltHb <- strsplit(Hbeta, split="")[[1]]
eltHa <- strsplit(Halpha, split="")[[1]]

AAdata   <- as.data.frame(read_excel(fileName))

setHb <- calFt(eltHb)
setHa <- calFt(eltHa)



selFeatureHbPos <- selFt(posFileName, setHb, 'mzPos')
selFeatureHaPos <- selFt(posFileName, setHa, 'mzPos')


write.xlsx2(selFeatureHbPos, file=outputFileName,
            sheetName="selFeatureHbPos",
            col.names=TRUE, row.names=TRUE, append=FALSE)
write.xlsx2(selFeatureHaPos, file=outputFileName,
            sheetName="selFeatureHaPos",
            col.names=TRUE, row.names=TRUE, append=TRUE)

selFeatureHbNeg <- selFt(negFileName, setHb, 'mzNeg')
selFeatureHaNeg <- selFt(negFileName, setHa, 'mzNeg')

merdata <- rbind(selFeatureHaNeg, selFeatureHaNeg)


write.xlsx2(selFeatureHbNeg, file=outputFileName,
            sheetName="selFeatureHbNeg",
            col.names=TRUE, row.names=TRUE, append=TRUE)
write.xlsx2(selFeatureHaNeg, file=outputFileName,
            sheetName="selFeatureHaNeg",
            col.names=TRUE, row.names=TRUE, append=TRUE)
