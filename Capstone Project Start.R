library(tidyr)
library(dplyr)

#Convert FASTA file into a dataframe
df <- data.frame(matrix(unlist(AQP5), nrow = 1, byrow = T), stringsAsFactors = FALSE) 

#Unite dataframe into single column for analysis
df1 <- unite(df, X1:X27) 
df2 <- gsub("_", "", df1)

#Split nucleotides into codons
df3 <- strsplit(df2, '(?<=.{3})', perl = TRUE)
df4 <- data.frame(matrix(unlist(df3), nrow = 1, byrow = T), stringsAsFactors = FALSE) 
final_df <- as.data.frame(t(df4))

#Create custom function to convert codons into amino acid counterpart
aa_conversion <- function(final_df) {
  #Alanine
  final_df = gsub("GCT", "Ala", final_df)
  final_df = gsub("GCC", "Ala", final_df)
  final_df = gsub("GCA", "Ala", final_df)
  final_df = gsub("GCG", "Ala", final_df)
  #Arginine
  final_df = gsub("CGG", "Arg", final_df)
  final_df = gsub("CGC", "Arg", final_df)
  final_df = gsub("CGT", "Arg", final_df)
  final_df = gsub("CGA", "Arg", final_df)
  final_df = gsub("AGA", "Arg", final_df)
  final_df = gsub("AGG", "Arg", final_df)
  #Asparagine
  final_df = gsub("AAT", "Asn", final_df)
  final_df = gsub("AAC", "Asn", final_df)
  #Aspartic Acid
  final_df = gsub("GAT", "Asp", final_df)
  final_df = gsub("GAC", "Asp", final_df)
  #Cysteine
  final_df = gsub("TGT", "Cys", final_df)
  final_df = gsub("TGC", "Cys", final_df)
  #Glutamic Acid
  final_df = gsub("GAA", "Glu", final_df)
  final_df = gsub("GAG", "Glu", final_df)
  #Glutamine
  final_df = gsub("CAA", "Gln", final_df)
  final_df = gsub("CAG", "Gln", final_df)
  #Glycine
  final_df = gsub("GGT", "Gly", final_df)
  final_df = gsub("GGC", "Gly", final_df)
  final_df = gsub("GGA", "Gly", final_df)
  final_df = gsub("GGG", "Gly", final_df)
  #Histidine
  final_df = gsub("CAT", "His", final_df)
  final_df = gsub("CAC", "His", final_df)
  #Isoleucine
  final_df = gsub("ATT", "Ile", final_df)
  final_df = gsub("ATC", "Ile", final_df)
  final_df = gsub("ATA", "Ile", final_df)
  #Leucine
  final_df = gsub("CTT", "Leu", final_df)
  final_df = gsub("CTC", "Leu", final_df)
  final_df = gsub("CTA", "Leu", final_df)
  final_df = gsub("CTG", "Leu", final_df)
  final_df = gsub("TTa", "Leu", final_df)
  final_df = gsub("TTG", "Leu", final_df)
  #Lysine
  final_df = gsub("AAA", "Lys", final_df)
  final_df = gsub("AAG", "Lys", final_df)
  #Methionine
  final_df = gsub("ATG", "Met", final_df)
  #Phenylalanine
  final_df = gsub("TTT", "Phe", final_df)
  final_df = gsub("TTC", "Phe", final_df)
  #Proline
  final_df = gsub("CCT", "Pro", final_df)
  final_df = gsub("CCC", "Pro", final_df)
  final_df = gsub("CCA", "Pro", final_df)
  final_df = gsub("CCG", "Pro", final_df)
  #Serine
  final_df = gsub("TCT", "Ser", final_df)
  final_df = gsub("TCC", "Ser", final_df)
  final_df = gsub("TCA", "Ser", final_df)
  final_df = gsub("TCG", "Ser", final_df)
  final_df = gsub("AGT", "Ser", final_df)
  final_df = gsub("AGC", "Ser", final_df)
  #Threonine
  final_df = gsub("ACC", "Thr", final_df)
  final_df = gsub("ACT", "Thr", final_df)
  final_df = gsub("ACA", "Thr", final_df)
  final_df = gsub("ACG", "Thr", final_df)
  #Tyrosine
  final_df = gsub("TAT", "Tyr", final_df)
  final_df = gsub("TAC", "Tyr", final_df)
  #Tryptophan
  final_df = gsub("TGG", "Trp", final_df)
  #Valine
  final_df = gsub("GTT", "Val", final_df)
  final_df = gsub("GTC", "Val", final_df)
  final_df = gsub("GTA", "Val", final_df)
  final_df = gsub("GTG", "Val", final_df)
  #Stop Codons
  final_df = gsub("TAA", "Stop", final_df)
  final_df = gsub("TAG", "Stop", final_df)
  final_df = gsub("TGA", "Stop", final_df)
  return(final_df)
}

#Convert codons into amino acids 
final_df$Amino_Acid <- sapply(final_df, aa_conversion)                 

#Create custom function to apply hydrophobic values to amino acids
hyd_values <- function(final_df){
  final_df = gsub("Ile", "4.5", final_df)
  final_df = gsub("Val", "4.2", final_df)
  final_df = gsub("Leu", "3.8", final_df)
  final_df = gsub("Phe", "2.8", final_df)
  final_df = gsub("Cys", "2.5", final_df)
  final_df = gsub("Met", "1.9", final_df)
  final_df = gsub("Ala", "1.8", final_df)
  final_df = gsub("Gly", "-0.4", final_df)
  final_df = gsub("Thr", "-0.7", final_df)
  final_df = gsub("Ser", "-0.8", final_df)
  final_df = gsub("Trp", "-0.9", final_df)
  final_df = gsub("Tyr", "-1.3", final_df)
  final_df = gsub("Pro", "-1.6", final_df)
  final_df = gsub("His", "-3.2", final_df)
  final_df = gsub("Glu", "-3.5", final_df)
  final_df = gsub("Gln", "-3.5", final_df)
  final_df = gsub("Asp", "-3.5", final_df)
  final_df = gsub("Asn", "-3.5", final_df)
  final_df = gsub("Lys", "-3.9", final_df)
  final_df = gsub("Arg", "-4.5", final_df)
  final_df = gsub("Stop", "NA", final_df)
}

#Assign amino acids hydrophobic values
final_df$Hyd_Value <- sapply(final_df$Amino_Acid, hyd_values)
  
#Change column names
colnames(final_df)[1] <- "Codon"
colnames(final_df)[2] <- "Amino Acid"
colnames(final_df)[3] <- "Hydrophobic Value"



