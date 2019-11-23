# Developers: [
# Cherlynn Silva <cdnsprado@gmail.com https://github.com/CDsilvaA>, 
# Thiago Roberto do Prado <trprado@outlook.com https://github.com/trprado>
# ]

# Example of data
# ----------- Crossover output -----------
#   I = Individual
#   S C{} = Sire Choromosome{Positions}
#   D C{} = Dam  Choromosome{Positions}
#   
#   I 115702
#   S C1{ 48.019651}
#   I 115704
#   D C1{ 38.691369 54.26158 56.290077 69.561402}
#   I 115706
#   D C1{ 85.880854}
#   I 115707
#   S C1{ 151.579322}


setwd("your_directory") 
dados <- read.delim("name_of_your_file_crosso_005.txt", header = T)

Crossovers <- function(set, chr) {
  # Returns a list containing matrices ordered with information about the animal, sire and dam with the positions 
  # on the chromosome where crossovers occurred in the genome. Adapted for the output of the QMSim crossover 
  # file.
  lista <- list()
  j = 0
  
  for (j in 1:chr) {
    ID <- NULL
    i = 0
    for (i in 1:length(set)) {
      sire <- NULL
      dam <- NULL
      vari <- strsplit(as.character(set[i]), "[ \t=]+")[[1]]
      
      if (vari[1] == "I" && vari[2] != "Individual")
      {
        ID <- as.numeric(vari[2])
        
        if (isTRUE(strsplit(as.character(set[i + 1]), "[ \tC{}]+")[[1]][2] == j) ||
            isTRUE(strsplit(as.character(set[i + 2]), "[ \tC{}]+")[[1]][2] == j)) {
          if (isTRUE(strsplit(as.character(set[i + 1]), "[ \tC{}]+")[[1]][1] == "S")
              &&
              isTRUE(strsplit(as.character(set[i + 2]), "[ \tC{}]+")[[1]][1] == "D"))
          {
            temp <- strsplit(as.character(set[i + 1]), "[ \tC{}]+")[[1]]
            temp2 <-
              strsplit(as.character(set[i + 2]), "[ \tC{}]+")[[1]]
            sire <- as.numeric(temp[3:length(temp)])
            dam <- as.numeric(temp2[3:length(temp2)])
          }
          
          else if (isTRUE(strsplit(as.character(set[i + 1]), "[ \tC{}]+")[[1]][1] == "S")
                   ||
                   isTRUE(strsplit(as.character(set[i + 1]), "[ \tC{}]+")[[1]][1] == "D"))
          {
            temp <- strsplit(as.character(set[i + 1]), "[ \tC{}]+")[[1]]
            
            if (temp[1] == "S")
              sire <- as.numeric(temp[3:length(temp)])
            
            else
              dam <- as.numeric(temp[3:length(temp)])
          }
          
          if (length(sire) > length(dam))
          {
            mm <- matrix(0, nrow = length(sire), ncol = 4)
            colnames(mm) <- c("ID", "SIRE", "DAM", "CHR")
            mm[1, 1] <- ID
            mm[1, 4] <- j
            mm[1:length(sire), 2] <- sire
            if (!is.null(dam))
              mm[1:length(dam), 3] <- dam
          }
          
          else
          {
            mm <- matrix(0, nrow = length(dam), ncol = 4)
            colnames(mm) <- c("ID", "SIRE", "DAM", "CHR")
            mm[1, 1] <- ID
            mm[1, 4] <- j
            mm[1:length(dam), 3] <- dam
            if (!is.null(sire))
              mm[1:length(sire), 2] <- sire
          }
          lista[[length(lista) + 1]] <- mm
        }
      }
    }
  }
  return(lista)
}
