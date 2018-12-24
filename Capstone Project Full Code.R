library(tidyr)
library(dplyr)
library(ggplot2)
library(plyr)
library(zoo)
library(ROCR)
library(C50)

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
  final_df = gsub("TTA", "Leu", final_df)
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
  final_df = gsub("Stop", "0", final_df)
}

#Convert FASTA file into a dataframe
df <- data.frame(matrix(unlist(AQP5), nrow = 1, byrow = T), stringsAsFactors = FALSE) #converts base code into data frame with one row and multiple columns
df <- df[ , colSums(is.na(df)) == 0] #removes all NAs from the data frame

#Unite dataframe into single column for analysis
df1 <- unite(df, X1) #unites all columns into one single column 
df2 <- gsub("_", "", df1) #creates character string vector and removes and spaces

#Split nucleotides into codons
df3 <- strsplit(df2, '(?<=.{3})', perl = TRUE) #splits the character string vector into a list of triplets i.e. codons
df4 <- data.frame(matrix(unlist(df3), nrow = 1, byrow = T), stringsAsFactors = FALSE) #creates a data frame from the df3 list with 1 row and each column containing a codon 
final_df <- as.data.frame(t(df4)) #moves each column from df4 into a separate row; result is one column with a row containing every codon

#Convert codons into amino acids 
final_df$Amino_Acid <- sapply(final_df, aa_conversion) #applies the aa_conversion function to every row of the final_df resulting in a new column with every amino acid of the protein                

#Assign amino acids hydrophobic values
final_df$Hyd_Value <- sapply(final_df$Amino_Acid, hyd_values) #applies the  hyd_values function to every row of the Amino_Acid column resulting in a new column with a hydrophobic value specific to each amino acid
  
#Change column names
colnames(final_df)[1] <- "Codon"
colnames(final_df)[2] <- "Amino Acid"
colnames(final_df)[3] <- "Hydrophobic Value"

#Constructing Hydropathy Plot
final_df$AA_position <- 1:nrow(final_df) #creates a new column that orders amino acids by position
final_df$`Hydrophobic Value` <- as.numeric(final_df$`Hydrophobic Value`) #makes hydrophobic values numeric
final_df$`Amino Acid` <- factor(final_df$`Amino Acid`) #makes amino acids factors

final_df$Index <- rollapplyr(final_df$`Hydrophobic Value`, 20, mean, partial=TRUE) #creates a hydrophobic index value for each amino acid position by calculating a rolling average of 20 hydrophobic values at a time

ggplot(final_df, aes(x = AA_position, y = Index)) + #creates hydropathy plot, areas above the blue line are representative of likely transmembrane regions
  geom_line() +
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 1.6, color = "blue") +
  ylim(-2, 2) +
  coord_fixed(ratio = 50)

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


##Exploratory Data Analysis##
##Summary Statistics for all Proteins##

#Calculate mean, variance and standard deviation for average hydrophobic and percentage hydrophic values
A.H.mean <- ddply(ProteinComparison_df, "Type", summarise, Average_Hydrophobic.mean = mean(Average_Hydrophobic)) #average hydrophobic mean for both protein types
P.H.mean <- ddply(ProteinComparison_df, "Type", summarise, Percentage_Hydrophobic.mean = mean(Percentage_Hydrophobic)) #percentage hydrophobic mean for both protein types
A.H.var <- ddply(ProteinComparison_df, "Type", summarise, Average_Hydrophobic.var = var(Average_Hydrophobic)) #average hydrophobic variance for both protein types
P.H.var <- ddply(ProteinComparison_df, "Type", summarise, Percentage_Hydrophobic.var = var(Percentage_Hydrophobic)) #percentage hydrophobci variance for both protein types
A.H.sd <- ddply(ProteinComparison_df, "Type", summarise, Average_Hydrophobic.sd = sd(Average_Hydrophobic)) #average hydrophobic standard deviation for both protein types
P.H.sd <- ddply(ProteinComparison_df, "Type", summarise, Percentage_Hydrophobic.sd = sd(Percentage_Hydrophobic)) #percentage hydrophobic standard deviation for both protein types

#Calculate mean index value for both types  
ProteinComparison_df$Index_Count <- as.numeric(ProteinComparison_df$Index_Count) #makes Index_Count numeric
Index.mean <- ddply(ProteinComparison_df, "Type", summarise, Index_Count.mean = mean(Index_Count)) #calculates average index count for both protein types

#Create data frame with all summary statistics and save
Summary <- full_join(Summary, P.H.sd) #combines all summary statistics into one data frame
write.csv(Summary, file = "SummaryStats") #saves summary statistics data frame 

#Average Hydrophobic Density Plot
ggplot(ProteinComparison_df, aes(x = Average_Hydrophobic, fill = Type)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = A.H.mean, aes(xintercept = Average_Hydrophobic.mean, colour = Type),
             linetype = "dashed", size = 1)

#Percentage Hydrophobic Density Plot
ggplot(ProteinComparison_df, aes(x = Percentage_Hydrophobic, fill = Type)) +
  geom_density(alpha = 0.3) +
  geom_vline(data = P.H.mean, aes(xintercept = Percentage_Hydrophobic.mean, colour = Type),
             linetype = "dashed", size = 1)

#Calculating mean for each amino acid based on type
  Ala.mean <- ddply(ProteinComparison_df, "Type", summarise, Ala = mean(Ala))
  Arg.mean <- ddply(ProteinComparison_df, "Type", summarise, Arg = mean(Arg))
  Asn.mean <- ddply(ProteinComparison_df, "Type", summarise, Asn = mean(Asn))
  Asp.mean <- ddply(ProteinComparison_df, "Type", summarise, Asp = mean(Asp))
  Cys.mean <- ddply(ProteinComparison_df, "Type", summarise, Cys = mean(Cys))
  Gln.mean <- ddply(ProteinComparison_df, "Type", summarise, Gln = mean(Gln))  
  Glu.mean <- ddply(ProteinComparison_df, "Type", summarise, Glu = mean(Glu))
  Gly.mean <- ddply(ProteinComparison_df, "Type", summarise, Gly = mean(Gly))
  His.mean <- ddply(ProteinComparison_df, "Type", summarise, His = mean(His))
  Ile.mean <- ddply(ProteinComparison_df, "Type", summarise, Ile = mean(Ile))
  Leu.mean <- ddply(ProteinComparison_df, "Type", summarise, Leu = mean(Leu))
  Met.mean <- ddply(ProteinComparison_df, "Type", summarise, Met = mean(Met))
  Lys.mean <- ddply(ProteinComparison_df, "Type", summarise, Lys = mean(Lys))  
  Phe.mean <- ddply(ProteinComparison_df, "Type", summarise, Phe = mean(Phe))
  Pro.mean <- ddply(ProteinComparison_df, "Type", summarise, Pro = mean(Pro))
  Ser.mean <- ddply(ProteinComparison_df, "Type", summarise, Ser = mean(Ser))
  Thr.mean <- ddply(ProteinComparison_df, "Type", summarise, Thr = mean(Thr))
  Trp.mean <- ddply(ProteinComparison_df, "Type", summarise, Trp = mean(Trp))
  Tyr.mean <- ddply(ProteinComparison_df, "Type", summarise, Tyr = mean(Tyr))
  Val.mean <- ddply(ProteinComparison_df, "Type", summarise, Val = mean(Val))
  
#Create data frame for all amino acid means
AA_mean <- full_join(AA_mean, Val.mean) #join all amino acid means into one data frame

#Bar chart of amino acid means based on type   
AA_mean2 <- gather(AA_mean, Amino_Acid, Mean, 2:21) #rearrange data frame for bar chart

ggplot(AA_mean2, aes(fill = Type, x = Amino_Acid, y = Mean)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust =1))

#T-tests
t.test(ProteinComparison_df$Leu[ProteinComparison_df$Type=="Transmembrane"], ProteinComparison_df$Leu[ProteinComparison_df$Type=="Non-membrane"])
t.test(ProteinComparison_df$Lys[ProteinComparison_df$Type=="Transmembrane"], ProteinComparison_df$Lys[ProteinComparison_df$Type=="Non-membrane"])
t.test(ProteinComparison_df$Cys[ProteinComparison_df$Type=="Transmembrane"], ProteinComparison_df$Cys[ProteinComparison_df$Type=="Non-membrane"])
t.test(ProteinComparison_df$Glu[ProteinComparison_df$Type=="Transmembrane"], ProteinComparison_df$Glu[ProteinComparison_df$Type=="Non-membrane"])
t.test(ProteinComparison_df$Trp[ProteinComparison_df$Type=="Transmembrane"], ProteinComparison_df$Trp[ProteinComparison_df$Type=="Non-membrane"])
t.test(ProteinComparison_df$Percentage_Hydrophobic[ProteinComparison_df$Type=="Transmembrane"], ProteinComparison_df$Percentage_Hydrophobic[ProteinComparison_df$Type=="Non-membrane"])
t.test(ProteinComparison_df$Average_Hydrophobic[ProteinComparison_df$Type=="Transmembrane"], ProteinComparison_df$Average_Hydrophobic[ProteinComparison_df$Type=="Non-membrane"])


##Model Creation and Prediction##
#Logistic Regression#

#Create new column for prediction
ProteinComparison_df$Type2 <- ProteinComparison_df$Type #create new column Type2
ProteinComparison_df$Type2 <- gsub("Transmembrane", "0", ProteinComparison_df$Type2) #assign 0 to transmembrane proteins
ProteinComparison_df$Type2 <- gsub("Non-membrane", "1", ProteinComparison_df$Type2) #assign 1 to non-membrane proteins
ProteinComparison_df$Type2 <- factor(ProteinComparison_df$Type2) #make Type2 column factors

#Split final data frame into training and testing sets
index <- sample(1:nrow(ProteinComparison_df),36) #split final data frame randomly 
train_df <- ProteinComparison_df[index, ] #place ~70% of data frame into training set
test_df <- ProteinComparison_df[-index, ] #place ~30% of data frame into testing set

#Create logistic regression model
model1 <- glm(Type2 ~ Average_Hydrophobic + Lys + Met + Trp, family = "binomial", data = train_df)
summary(model1)
anova(model1, test = "Chisq") #measure of significance for each feature

#Model prediction with ROC curve and auc
prob <- predict(model1, newdata = test_df, type = "response")
pred <- prediction(prob, test_df$Type2)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf) #plots ROC curve 
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc #value represents percentage of accurately predicted proteins

#Create Decision Tree model
train_tree <- C5.0(Type2 ~ Average_Hydrophobic + Met + Trp, data = train_df)
plot(train_tree)

#Decision Tree model predicition 
results <- predict(object = train_tree, newdata = test_df, type = "class")
table(results, test_df$Type2)

