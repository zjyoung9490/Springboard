library(tidyr)
library(dplyr)
library(ggplot2)
library(plyr)
library(zoo)
library(ROCR)
library(C50)

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

