library(tidyr)
library(dplyr)
library(ggplot2)
library(plyr)
library(zoo)

##Summary Statistics for Individual Protein##
aa_vector <- table(final_df$`Amino Acid`) #creates a table of all the amino acids from the protein
unlist(aa_vector) #creates a count of every amino acid
all_amino_acids <- aa_vector[c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")] #assigns the count of every amino acid to a table all_amino_acids
hydrophobic_amino_acids <- aa_vector[c("Ala", "Cys", "Ile", "Leu", "Met", "Phe", "Val")] #assigns the hydrophobic amino acids to a table hydrophobic_amino_acids


#Average hydrophobic value of entire protein
Mean_HydrophobicValue <- mean(final_df$`Hydrophobic Value`, na.rm = TRUE) #calculates the average hydrophobic value of the protein

#Percentage of all amino acids in protein
Percentage_AminoAcids <- (all_amino_acids/nrow(final_df)) *100 #calculates the percentage of every amino acid in the protein

#Percentage of all hydrophobic amino acids in protein
Percentage_HydrophobicAminoAcids <- ((sum(hydrophobic_amino_acids))/nrow(final_df)) *100 #calculates the sum percentage of the all amino acids with a positive hydrophobic value

#Number of times Index is greater or equal to 1.6
Total_Index <- sum(final_df$Index >= 1.6, na.rm = TRUE) #calculates the total amount of times the index value is greater than 1.6 which likely signifies a transmembrane region 

#Clean up data frame
ProteinSummary <- data.frame(Mean_HydrophobicValue, Percentage_AminoAcids, Percentage_HydrophobicAminoAcids, Total_Index) #places the above statistics into a data frame
ProteinSummary <- spread(ProteinSummary, Var1, Freq) #places all features into separate columns
ProteinSummary["Type"] <- "Non-membrane" #creates a column for the protein type; this is researched beforehand 
ProteinSummary["Protein"] <- "AQP5" #creates a column for the protein name
ProteinSummary <- ProteinSummary[colnames(ProteinSummary)[c(25, 24, 1:23)]] #reorganizes the columns
colnames(ProteinSummary)[3] <- "Average_Hydrophobic" #change column name
colnames(ProteinSummary)[4] <- "Percentage_Hydrophobic" #change column name
colnames(ProteinSummary)[5] <- "Index_Count" #change column name


#Add data frame to final data frame and save
ProteinComparison_df <- full_join(ProteinSummary, ProteinComparison_df) #adds recently completed protein data frame to final, aggregated protein data frame
write.csv(ProteinComparison_df, file = "ProteinComparison_df") #saves final data frame 