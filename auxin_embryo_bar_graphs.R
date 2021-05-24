library(plyr)
library(readr)
library(ggplot2)
library("ggpubr", lib.loc="~/Library/R/3.6/library")
library(dplyr)
library(reshape2)

emf97 <- read.table("~/Dropbox (hannonlab)/Sequencing Data/roo/chipseq/emf97/individual/count_het_chr4_roo/emf97_ind_count_roo_het_chr4.txt", header=T, fill = T, as.is = T)
emf97_roo <- as.data.frame(emf97[,2])
emf97_roo[,2] <-emf97[,6]*1000000/emf97[,3]

colnames(emf97_roo) <- c("timepoint", "value")
emf97_roo$timepoint <- factor(emf97_roo$timepoint,levels = unique(emf97_roo$timepoint))

emf97_het <- as.data.frame(emf97[,2])
emf97_het[,2] <-emf97[,4]*1000000/emf97[,3]

colnames(emf97_het) <- c("timepoint", "value")
emf97_het$timepoint <- factor(emf97_het$timepoint,levels = unique(emf97_het$timepoint))

t.test(emf97_het[grep("ctrl",emf97_het[,1]),2], emf97_het[grep("auxin",emf97_het[,1]),2],paired = FALSE, alternative = "two.sided")


emf97_chr4 <- as.data.frame(emf97[,2])
emf97_chr4[,2] <-emf97[,5]*1000000/emf97[,3]

colnames(emf97_chr4) <- c("timepoint", "value")
emf97_chr4$timepoint <- factor(emf97_chr4$timepoint,levels = unique(emf97_chr4$timepoint))

t.test(emf97_chr4[grep("ctrl",emf97_chr4[,1]),2], emf97_chr4[grep("auxin",emf97_chr4[,1]),2],paired = FALSE, alternative = "two.sided")


# Bar graphs --------------------------------------------------------------


# chr4

emf97_chr4_reps <- emf97_chr4 %>% # the names of the new data frame and the data frame to be summarised
  group_by(timepoint) %>%   # the grouping variable
  summarise(mean_PL = mean(value),  # calculates the mean of each group
            sd_PL = sd(value), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(value)/sqrt(n())) # calculates the standard error of each group

emf97_chr4$timepoint <- factor(emf97_chr4$timepoint, levels = unique(emf97_chr4$timepoint))


ggplot(emf97_chr4_reps, aes( y=mean_PL, x=timepoint)) + 
  geom_bar(position="dodge", stat="identity", fill = NA, colour = "red") +
  geom_errorbar(aes(ymin = mean_PL - sd_PL, ymax = mean_PL + sd_PL), width=0.2) +
  geom_point(data = emf97_chr4, aes(y=value, x=timepoint),
             stat="identity",
             size=2,
             shape=16,
             position = position_dodge(width = .9)) +
  geom_hline(yintercept=0) +
  scale_y_continuous(expand = c(0,0),limits = c(0,20000)) +
  scale_fill_manual(values = c("brown1", "skyblue2")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))  

t.test(emf97_chr4[grep("ctrl", emf97_chr4[,1]),2],emf97_chr4[grep("auxin", emf97_chr4[,1]),2],paired = TRUE, alternative = "two.sided")

ggsave(filename="~/Desktop/figure4J_chr4_bar_k9.pdf", plot = last_plot(), width=2, height=4)

# het

emf97_het_reps <- emf97_het %>% # the names of the new data frame and the data frame to be summarised
  group_by(timepoint) %>%   # the grouping variable
  summarise(mean_PL = mean(value),  # calculates the mean of each group
            sd_PL = sd(value), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(value)/sqrt(n())) # calculates the standard error of each group

emf97_het$timepoint <- factor(emf97_het$timepoint, levels = unique(emf97_het$timepoint))


ggplot(emf97_het_reps, aes( y=mean_PL, x=timepoint)) + 
  geom_bar(position="dodge", stat="identity", fill = NA, colour = "red") +
  geom_errorbar(aes(ymin = mean_PL - sd_PL, ymax = mean_PL + sd_PL), width=0.2) +
  geom_point(data = emf97_het, aes(y=value, x=timepoint),
             stat="identity",
             size=2,
             shape=16,
             position = position_dodge(width = .9)) +
  geom_hline(yintercept=0) +
  scale_y_continuous(expand = c(0,0),limits = c(0,500000)) +
  scale_fill_manual(values = c("brown1", "skyblue2")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))  

t.test(emf97_het[grep("ctrl", emf97_het[,1]),2],emf97_het[grep("auxin", emf97_het[,1]),2],paired = TRUE, alternative = "two.sided")

ggsave(filename="~/Desktop/figure4I_het_bar_k9.pdf", plot = last_plot(), width=2, height=4)

# het

emf97_roo_reps <- emf97_roo %>% # the names of the new data frame and the data frame to be summarised
  group_by(timepoint) %>%   # the grouping variable
  summarise(mean_PL = mean(value),  # calculates the mean of each group
            sd_PL = sd(value), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(value)/sqrt(n())) # calculates the standard error of each group

emf97_roo$timepoint <- factor(emf97_roo$timepoint, levels = unique(emf97_roo$timepoint))


ggplot(emf97_roo_reps, aes( y=mean_PL, x=timepoint)) + 
  geom_bar(position="dodge", stat="identity", fill = NA, colour = "red") +
  geom_errorbar(aes(ymin = mean_PL - sd_PL, ymax = mean_PL + sd_PL), width=0.2) +
  geom_point(data = emf97_roo, aes(y=value, x=timepoint),
             stat="identity",
             size=2,
             shape=16,
             position = position_dodge(width = .9)) +
  geom_hline(yintercept=0) +
  scale_y_continuous(expand = c(0,0),limits = c(0,25000)) +
  scale_fill_manual(values = c("brown1", "skyblue2")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))  

t.test(emf97_het[grep("ctrl", emf97_roo[,1]),2],emf97_roo[grep("auxin", emf97_roo[,1]),2],paired = TRUE, alternative = "two.sided")

ggsave(filename="~/Desktop/figure4G_roo_bar_k9.pdf", plot = last_plot(), width=2, height=4)
