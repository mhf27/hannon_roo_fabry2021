library(plyr)
library(readr)
library(ggplot2)
library("ggpubr", lib.loc="~/Library/R/3.6/library")
library(dplyr)
library(reshape2)

emf55 <- read.table("~/Dropbox (hannonlab)/Sequencing Data/roo/chipseq/emf55/new_sep20/individual/count/emf55_ind_count_roo_het_chr4.txt", header=T, fill = T, as.is = T)
emf55_roo <- as.data.frame(emf55[,2])
emf55_roo[,2] <-emf55[,4]*1000000/emf55[,3]

colnames(emf55_roo) <- c("timepoint", "value")
emf55_roo$timepoint <- factor(emf55_roo$timepoint,levels = unique(emf55_roo$timepoint))

ggboxplot(emf55_roo, x = "timepoint", y = "value", fill = "timepoint") +
  ylab("Signal Intensity in RPM") +
  xlab("time points") +
  scale_y_continuous(breaks=c(c(10000,20000, 30000, 40000, 50000)), limits = c(0,50000)) +
  scale_fill_brewer(palette="Set3") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"))

ggsave(filename="~/Dropbox (hannonlab)/Sequencing Data/roo/chipseq/emf55/new_sep20/individual/count/emf55_boxplot_roo.pdf", plot = last_plot(), width=4, height=4)


emf55_het <- as.data.frame(emf55[,2])
emf55_het[,2] <-emf55[,5]*1000000/emf55[,3]

colnames(emf55_het) <- c("timepoint", "value")
emf55_het$timepoint <- factor(emf55_het$timepoint,levels = unique(emf55_het$timepoint))

ggboxplot(emf55_het, x = "timepoint", y = "value", fill = "timepoint") +
  ylab("Signal Intensity in RPM") +
  xlab("time points") +
  scale_y_continuous(breaks=c(c(200000,400000,600000)), limits = c(0,600000)) +
  scale_fill_brewer(palette="Set3") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"))

ggsave(filename="~/Dropbox (hannonlab)/Sequencing Data/roo/chipseq/emf55/new_sep20/individual/count/emf55_boxplot_het.pdf", plot = last_plot(), width=4, height=4)
emf55_chr4 <- as.data.frame(emf55[,2])
emf55_chr4[,2] <-emf55[,6]*1000000/emf55[,3]

colnames(emf55_chr4) <- c("timepoint", "value")
emf55_chr4$timepoint <- factor(emf55_chr4$timepoint,levels = unique(emf55_chr4$timepoint))

ggboxplot(emf55_chr4, x = "timepoint", y = "value", fill = "timepoint") +
  ylab("Signal Intensity in RPM") +
  xlab("time points") +
  scale_y_continuous(breaks=c(c(10000,20000,30000,40000,50000)), limits = c(0,50000)) +
  scale_fill_brewer(palette="Set3") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.x = element_text(size = 15, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 15, angle = 90, hjust = .5, vjust = .5, face = "plain"))

ggsave(filename="~/Dropbox (hannonlab)/Sequencing Data/roo/chipseq/emf55/new_sep20/individual/count/emf55_boxplot_chr4.pdf", plot = last_plot(), width=4, height=4)




### Bar graphs

# roo

emf55_roo_reps <- emf55_roo %>% # the names of the new data frame and the data frame to be summarised
  group_by(timepoint) %>%   # the grouping variable
  summarise(mean_PL = mean(value),  # calculates the mean of each group
            sd_PL = sd(value), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(value)/sqrt(n())) # calculates the standard error of each group

emf55_roo_reps$timepoint <- factor(emf55_roo_reps$timepoint, levels = unique(emf55_roo$timepoint))

emf55_roo_reps <- emf55_roo_reps[1:6,]
emf55_roo <- emf55_roo[1:12,]

ggplot(emf55_roo_reps, aes( y=mean_PL, x=timepoint)) + 
  geom_bar(position="dodge", stat="identity", fill = NA, colour = "red") +
  geom_errorbar(aes(ymin = mean_PL - sd_PL, ymax = mean_PL + sd_PL), width=0.2) +
  geom_point(data = emf55_roo, aes(y=value, x=timepoint),
             stat="identity",
             size=2,
             shape=16,
             position = position_dodge(width = .9)) +
  geom_hline(yintercept=0) +
  scale_y_continuous(expand = c(0,0),limits = c(0,45000)) +
  scale_fill_manual(values = c("brown1", "skyblue2")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, colour="black"), axis.text.y = element_text( colour="black"), axis.ticks = element_line(colour="black"))

ggsave(filename="~/Desktop/figure3F_roo_bar_k9.pdf", plot = last_plot(), width=4, height=4)

#ggsave(filename="~/Dropbox (hannonlab)/Sequencing Data/roo/chipseq/emf55/new_sep20/individual/count/emf55_barplot_roo.pdf", plot = last_plot(), width=4, height=4)

# chr4

emf55_chr4_reps <- emf55_chr4 %>% # the names of the new data frame and the data frame to be summarised
  group_by(timepoint) %>%   # the grouping variable
  summarise(mean_PL = mean(value),  # calculates the mean of each group
            sd_PL = sd(value), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(value)/sqrt(n())) # calculates the standard error of each group

emf55_chr4_reps$timepoint <- factor(emf55_chr4_reps$timepoint, levels = unique(emf55_roo$timepoint))

emf55_chr4_reps <- emf55_chr4_reps[1:6,]
emf55_chr4 <- emf55_chr4[1:12,]

ggplot(emf55_chr4_reps, aes( y=mean_PL, x=timepoint)) + 
  geom_bar(position="dodge", stat="identity", fill = NA, colour = "red") +
  geom_errorbar(aes(ymin = mean_PL - sd_PL, ymax = mean_PL + sd_PL), width=0.2) +
  geom_point(data = emf55_chr4, aes(y=value, x=timepoint),
             stat="identity",
             size=2,
             shape=16,
             position = position_dodge(width = .9)) +
  geom_hline(yintercept=0) +
  scale_y_continuous(expand = c(0,0),limits = c(0,50000)) +
  scale_fill_manual(values = c("brown1", "skyblue2")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, colour="black"), axis.text.y = element_text( colour="black"), axis.ticks = element_line(colour="black"))

ggsave(filename="~/Desktop/figure3E_chr4_bar_k9.pdf", plot = last_plot(), width=4, height=4)

#ggsave(filename="~/Dropbox (hannonlab)/Sequencing Data/roo/chipseq/emf55/new_sep20/individual/count/emf55_barplot_chr4.pdf", plot = last_plot(), width=4, height=4)


# het

emf55_het_reps <- emf55_het %>% # the names of the new data frame and the data frame to be summarised
  group_by(timepoint) %>%   # the grouping variable
  summarise(mean_PL = mean(value),  # calculates the mean of each group
            sd_PL = sd(value), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(value)/sqrt(n())) # calculates the standard error of each group

emf55_het_reps$timepoint <- factor(emf55_het_reps$timepoint, levels = unique(emf55_roo$timepoint))

emf55_het_reps <- emf55_het_reps[1:6,]
emf55_het <- emf55_het[1:12,]

ggplot(emf55_het_reps, aes( y=mean_PL, x=timepoint)) + 
  geom_bar(position="dodge", stat="identity", fill = NA, colour = "red") +
  geom_errorbar(aes(ymin = mean_PL - sd_PL, ymax = mean_PL + sd_PL), width=0.2) +
  geom_point(data = emf55_het, aes(y=value, x=timepoint),
             stat="identity",
             size=2,
             shape=16,
             position = position_dodge(width = .9)) +
  geom_hline(yintercept=0) +
  scale_y_continuous(expand = c(0,0),limits = c(0,400000)) +
  scale_fill_manual(values = c("brown1", "skyblue2")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1))  

ggsave(filename="~/Desktop/figure3D_het_bar_k9.pdf", plot = last_plot(), width=4, height=4)

#ggsave(filename="~/Dropbox (hannonlab)/Sequencing Data/roo/chipseq/emf55/new_sep20/individual/count/emf55_barplot_het.pdf", plot = last_plot(), width=4, height=4)
