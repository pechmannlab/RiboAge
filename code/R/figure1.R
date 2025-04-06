# RiboAge | Figure 1
# S Pechmann (2025)

library(ggplot2)
library(reshape2)
library(cowplot)
library(EnvStats)

setwd("~/M2/RiboAge")



# --------------------- Figure 1 A2/B2 ---------------------
scer_rd <- read.table("data/processed/scer_ribodensities.txt", header=T)
scer_allRD <- as.matrix(read.table("data/processed/scer_allRD.txt"))
scer_allRD <- as.vector(scer_allRD[is.na(scer_allRD)==F])
scer_allRD <- scer_allRD[scer_allRD<50]       # remove most extreme outliers

scer_young <- demp(scer_rd$RD_young, scer_allRD)
scer_middle <- demp(scer_rd$RD_middle, scer_allRD)
scer_old <- demp(scer_rd$RD_old, scer_allRD) 

scer_df = data.frame(entropy=c(-sum(scer_young*log(scer_young)), -sum(scer_middle*log(scer_middle)), -sum(scer_old*log(scer_old))),
                     class=c("young", "middle", "old"), age=c(0, 2, 4))
scer_df$rel_entropy = scer_df$entropy / scer_df$entropy[1]




entropy_noise <- function(RDdf, allRD, err){
  addnoise <- function(rd, err, N){
    res <- matrix(NA, nrow=length(rd), ncol=N)
    for (i in 1:length(rd)){
      x <- rd[i]
      res[i,] <- runif(N, (1-(2*err))*x, (1+(2*err))*x)
    }
    res
  }
  N <- 1000
  rd_young <- RDdf$RD_young
  rd_middle <- RDdf$RD_middle
  rd_old <- RDdf$RD_old
  noise_young <- addnoise(rd_young, err, N)
  noise_middle <- addnoise(rd_middle, err, N)
  noise_old <- addnoise(rd_old, err, N)
  P_young <- array( demp(as.vector(noise_young), allRD), dim(noise_young) )
  P_middle <- array( demp(as.vector(noise_middle), allRD), dim(noise_middle) )
  P_old <- array( demp(as.vector(noise_old), allRD), dim(noise_old) )
  res <- matrix(NA, nrow=N, ncol=3)
  for (i in 1:N){
    S_young <- -sum(P_young[,i]*log(P_young[,i]))
    S_middle <- -sum(P_middle[,i]*log(P_middle[,i]))
    S_old <- -sum(P_old[,i]*log(P_old[,i]))
    S <- c(S_young, S_middle, S_old)
    res[i,] <- S / S[1]
  }
  res
}



scer_r_5 <- entropy_noise(scer_rd, scer_allRD, 0.05)
scer_r_10 <- entropy_noise(scer_rd, scer_allRD, 0.1)
scer_r_20 <- entropy_noise(scer_rd, scer_allRD, 0.2)
scer_err <- data.frame(age=c(0, 2, 4), err5=colMeans(scer_r_5), err10=colMeans(scer_r_10), err20=colMeans(scer_r_20))
scer_err.m <- melt(scer_err, id="age")


svg("figures/figure1/scer_entropy.svg", height=5, width=5)

ggplot(scer_df, aes(x=age, y=rel_entropy)) + 
  geom_line(data=scer_err.m, inherit.aes = F, aes(x=age, y=value, color=variable), linetype=5, size=1.1) + 
  geom_line(color="#228b22", size=2) + 
  geom_point(size=2.5, shape=15) + 
  scale_color_manual(values=c("#777777", "#AAAAAA", "#DDDDDD")) + 
  labs(x="Age [d]", y="Rel. information content") + 
  scale_x_continuous(breaks=scer_df$age, labels=scer_df$age) + 
  theme_classic() +
  theme(
    text = element_text(size=20),
    axis.title = element_text(size=28),
    axis.text = element_text(size=24),
    legend.title = element_blank(),
    legend.text = element_text(size=18),
    legend.position = c(0.2, 0.2)
  )

dev.off()



#~~~~~~~~~~~

cele_rd <- read.table("data/processed/cele_ribodensities.txt", header=T)
cele_allRD <- as.matrix(read.table("data/processed/cele_allRD.txt"))
cele_allRD <- as.vector(cele_allRD[is.na(cele_allRD)==F])
cele_allRD <- cele_allRD[cele_allRD<50]

cele_young <- demp(cele_rd$RD_young, cele_allRD)
cele_middle <- demp(cele_rd$RD_middle, cele_allRD)
cele_old <- demp(cele_rd$RD_old, cele_allRD) 

cele_df = data.frame(entropy=c(-sum(cele_young*log(cele_young)), -sum(cele_middle*log(cele_middle)), -sum(cele_old*log(cele_old))),
                     class=c("young", "middle", "old"), age=c(1, 6, 12))
cele_df$rel_entropy = cele_df$entropy / cele_df$entropy[1]


cele_r_5 <- entropy_noise(cele_rd, cele_allRD, 0.05)
cele_r_10 <- entropy_noise(cele_rd, cele_allRD, 0.1)
cele_r_20 <- entropy_noise(cele_rd, cele_allRD, 0.2)

cele_err <- data.frame(age=c(1, 6, 12), err5=colMeans(cele_r_5), err10=colMeans(cele_r_10), err20=colMeans(cele_r_20))
cele_err.m <- melt(cele_err, id="age")


svg("figures/figure1/cele_entropy.svg", height=5, width=5)

ggplot(cele_df, aes(x=age, y=rel_entropy)) + 
  geom_line(data=cele_err.m, inherit.aes = F, aes(x=age, y=value, color=variable), linetype=5, size=1.1) + 
  geom_line(color="#9826bb", size=2) + 
  geom_point(size=2.5, shape=15) + 
  scale_color_manual(values=c("#777777", "#AAAAAA", "#DDDDDD")) + 
  labs(x="Age [d]", y="Rel. information content") + 
  scale_x_continuous(breaks=cele_df$age, labels=cele_df$age) + 
  theme_classic() +
  theme(
    text = element_text(size=20),
    axis.title = element_text(size=28),
    axis.text = element_text(size=24),
    legend.title = element_blank(),
    legend.text = element_text(size=18),
    legend.position = c(0.2, 0.2)
  )

dev.off()



# --------------------- Figure 1 A1/B1 ---------------------
scer_rd <- read.table("data/processed/scer_ribodensities.txt", header=T)

svg("figures/figure1/scer_rd.svg", height=5, width=5)

ggplot(scer_rd, aes(x=RD_young)) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.3, color="#777777", linetype="twodash" ) +
  geom_vline(aes(xintercept = summary(RD_young)[2]), size=0.3, color="#222222", alpha=0.5 ) +
  geom_vline(aes(xintercept = summary(RD_young)[5]), size=0.3, color="#222222", alpha=0.5 ) +
  #geom_hline(aes(yintercept = summary(RD_middle)[2]), size=0.3, col="#2cb42c", alpha=0.5) +
  #geom_hline(aes(yintercept = summary(RD_middle)[5]), size=0.3, col="#2cb42c", alpha=0.5) +
  geom_hline(aes(yintercept = summary(RD_old)[2]), size=0.3, color="#186218", alpha=0.5) +
  geom_hline(aes(yintercept = summary(RD_old)[5]), size=0.3, color="#186218", alpha=0.5) +
  geom_point(aes(y=RD_middle), col="#2cb42c", size=2.5, shape=15, alpha=0.7) + 
  geom_point(aes(y=RD_old), col="#186218", size=2.5) + 
  labs(x="Codon RD (0d)", y="Codon RD (aged)") + 
  scale_x_continuous(breaks=c(0.5,1,1.5,2), labels=c(0.5,1,1.5,2), limits=c(0.5, 2.38)) + 
  scale_y_continuous(breaks=c(0.5,1,1.5,2), labels=c(0.5,1,1.5,2), limits=c(0.5, 2.38)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.title = element_text(size=28),
    axis.text = element_text(size=24)
  )

dev.off()



cele_rd <- read.table("data/processed/cele_ribodensities.txt", header=T)

svg("figures/figure1/cele_rd.svg", height=5, width=5)

ggplot(cele_rd, aes(x=RD_young)) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.3, color="#777777", linetype="twodash" ) +
  geom_vline(aes(xintercept = summary(RD_young)[2]), size=0.3, color="#222222", alpha=0.5 ) +
  geom_vline(aes(xintercept = summary(RD_young)[5]), size=0.3, color="#222222", alpha=0.5 ) +
  #geom_hline(aes(yintercept = summary(RD_middle)[2]), size=0.3, col="#b33dd7", alpha=0.5) +
  #geom_hline(aes(yintercept = summary(RD_middle)[5]), size=0.3, col="#b33dd7", alpha=0.5) +
  geom_hline(aes(yintercept = summary(RD_old)[2]), size=0.3, color="#761d91", alpha=0.5) +
  geom_hline(aes(yintercept = summary(RD_old)[5]), size=0.3, color="#761d91", alpha=0.5) +
  geom_point(aes(y=RD_middle), col="#b33dd7", size=2.5, shape=15, alpha=0.7) + 
  geom_point(aes(y=RD_old), col="#761d91", size=2.5) + 
  labs(x="Codon RD (1d)", y="Codon RD (aged)") + 
  scale_x_continuous(breaks=c(1,1.5, 2, 2.5), labels=c(1,1.5, 2, 2.5), limits=c(0.5, 2.5)) + 
  scale_y_continuous(breaks=c(1,1.5, 2, 2.5), labels=c(1,1.5, 2, 2.5), limits=c(0.5, 2.5)) + 
  theme_classic() + 
  theme(
    text = element_text(size=20),
    axis.title = element_text(size=28),
    axis.text = element_text(size=24)
  )

dev.off()



scer_corr <- read.table("data/processed/scer_corr.txt", header=T)

s1 <- scer_corr$young - scer_corr$nyoung
s2 <- scer_corr$middle - scer_corr$nmiddle
s3 <- scer_corr$old - scer_corr$nold
boxplot(s1,s2,s3)
abline(h=0)




cele_corr <- read.table("data/processed/cele_corr.txt", header=T)

s1 <- cele_corr$young - cele_corr$nyoung
s2 <- cele_corr$middle - cele_corr$nmiddle
s3 <- cele_corr$old - cele_corr$nold
boxplot(s1,s2,s3)
abline(h=0)







# --------------------- Figure 1C ---------------------
scer_rd <- read.table("data/processed/scer_ribodensities.txt", header=T)
cele_rd <- read.table("data/processed/cele_ribodensities.txt", header=T)

par(mfrow=c(2,1))
barplot(t(as.matrix(scer_rd[,c(3:5)])), beside=T, col=c(2, 3, 4), names=scer_rd$aa)
barplot(t(as.matrix(cele_rd[,c(3:5)])), beside=T, col=c(2, 3, 4), names=cele_rd$aa)


scer_rd2 <- rbind( scer_rd[scer_rd$aa=="T",], scer_rd[scer_rd$aa=="A",], scer_rd[scer_rd$aa=="G",] )
scer_rd2m <- melt(scer_rd2, index=codon)

plot.scer <- ggplot(scer_rd2m, aes(x=codon, y=value)) + 
  geom_col(aes(fill=variable), position=position_dodge2()) +
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) +
  labs(x="", y="Ribosome density") + 
  theme_classic() + 
  theme(
    text = element_text(size=22),
    axis.text.x = element_text(size=18, angle=90, hjust=1, vjust=0.5), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.title = element_blank()
  )



cele_rd2 <- rbind( cele_rd[cele_rd$aa=="T",], cele_rd[cele_rd$aa=="A",], cele_rd[cele_rd$aa=="G",] )
cele_rd2m <- melt(cele_rd2, index=codon)

plot.cele <- ggplot(cele_rd2m, aes(x=codon, y=value)) + 
  geom_col(aes(fill=variable), position=position_dodge2()) +
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) +
  labs(x="", y="Ribosome density") + 
  theme_classic() + 
  theme(
    text = element_text(size=22),
    axis.text.x = element_text(size=18, angle=90, hjust=1, vjust=0.5), 
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.title = element_blank()
  )


svg("figures/figure1/highlight_rd.svg", height=5, width=7)

plot_grid(plot.scer, plot.cele, labels ="", ncol = 1, align = 'v')

dev.off()




# --------------------- Figure 1D ---------------------
inhibitory = c('CGACCG', 'CGAGCG', 'CGAATA', 'CTCCCG', 'CGACGA', 'CGACGG', 'CGACTG', 'GTACCG', 'GTGCGA', 'GTACGA', 'CTGCCG', 'CTGCGA', 'ATACGA', 'ATACGG', 'AGGCGG', 'CTGATA', 'AGGCGA')
good_inhibitory = c("CGAATA", "GTACCG", "GTGCGA", "GTACGA", "CTGCCG", "CTGCGA", "ATACGA", "ATACGG", "AGGCGG", "CTGATA", "AGGCGA")

scer_pair <- read.table("data/processed/scer_pairdensities.txt", header=T)

idx_scer <- c()
for (i in 1:length(inhibitory)){
  idx_scer[i] <- which(scer_pair$codon == inhibitory[i])
}


colnames(scer_pair) <- c("codon", " all(0d)", " all(2d)", " all(4d)")
pair_all_scer <- melt(scer_pair)
inhib_scer <- scer_pair[idx_scer,]
colnames(inhib_scer) <- c("codon", " inh(0d)", " inh(2d)", " inh(4d)")
pair_inhib_scer <- melt(inhib_scer)
pair_merged_scer <- rbind(pair_all_scer, pair_inhib_scer)

svg("figures/figure1/scer_pair.svg", height=5, width=3)

ggplot(pair_merged_scer) + 
  geom_boxplot(aes(x=variable, y=value, fill=variable) ) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218", "#888888", "#555555", "#222222")) + 
  scale_y_continuous(limits=c(0, 8.5)) + 
  theme_classic() + 
  labs(x="", y="Codon pair RD") + 
  theme(
    text = element_text(size=28), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    legend.position = 'none'
  )

dev.off()






cele_pair <- read.table("data/processed/cele_pairdensities.txt", header=T)

idx_cele <- c()
for (i in 1:length(inhibitory)){
  idx_cele[i] <- which(cele_pair$codon == inhibitory[i])
}

colnames(cele_pair) <- c("codon", "all(1d)", "all(6d)", "all(12d)")
pair_all_cele <- melt(cele_pair)
inhib_cele <- cele_pair[idx_cele,]
colnames(inhib_cele) <- c("codon", "inh(1d)", "inh(6d)", "inh(12d)")
pair_inhib_cele <- melt(inhib_cele)
pair_merged_cele <- rbind(pair_all_cele, pair_inhib_cele)

svg("figures/figure1/cele_pair.svg", height=5, width=3)

ggplot(pair_merged_cele) + 
  geom_boxplot(aes(x=variable, y=value, fill=variable) ) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91", "#888888", "#555555", "#222222")) + 
  theme_classic() + 
  labs(x="", y="Codon pair RD") + 
  theme(
    text = element_text(size=28), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
    legend.position = 'none'
  )

dev.off()




# --------------------- Figure 1D2 ---------------------
# reload as changed colnames above (!!)
scer_pair <- read.table("data/processed/scer_pairdensities.txt", header=T)
scer_deltaDF <- data.frame(dYO=scer_pair$RD_young-scer_pair$RD_old, dYM=scer_pair$RD_young-scer_pair$RD_middle)

svg("figures/figure1/scer_pair_delta.svg", height=5, width=8)
ggplot(scer_deltaDF) + 
  geom_vline(aes(xintercept = 2*sd(dYO, na.rm=T)), size=1.2, color="#222222",  linetype='dashed') +
  geom_vline(aes(xintercept =-2*sd(dYO, na.rm=T)), size=1.2, color="#222222", linetype='dashed') +
  geom_histogram(aes(x=dYO), binwidth=0.02, fill="#228b22") + 
  scale_x_continuous(limits=c(-1, 1)) + 
  labs(x=paste("\u0394", "RD = RD(0d) - RD(4d)", sep=""), y="Count") + 
  theme_classic() + 
  theme(
    text = element_text(size=48)
  )
dev.off()


scer_extrDF <- data.frame(  dYO=c(sum(scer_deltaDF$dYO < -2*sd(scer_deltaDF$dYO, na.rm=T), na.rm=T), sum(scer_deltaDF$dYO > 2*sd(scer_deltaDF$dYO, na.rm=T), na.rm=T))/3723*100,
                            dYM=c(sum(scer_deltaDF$dYM < -2*sd(scer_deltaDF$dYM, na.rm=T), na.rm=T), sum(scer_deltaDF$dYM > 2*sd(scer_deltaDF$dYM, na.rm=T), na.rm=T))/3723*100, 
                            cat=c("up", "down") )

scer_extr.m <- melt(scer_extrDF, id='cat')
scer_extr.m$cat <- factor(scer_extr.m$cat, levels=c("up", "down"))
scer_extr.m$variable <- factor(scer_extr.m$variable, levels=c("dYM", "dYO"))
scer_extr.m$label <- c("4d", "4d", "2d", "2d")

svg("figures/figure1/scer_pair_updown.svg", height=5, width=5)

ggplot(scer_extr.m) + 
  geom_col(aes(x=label, y=value, fill=cat), position=position_dodge2()) + 
  scale_fill_manual(values=c("#739bd0", "#0f4c81")) + 
  labs(x="", y="% codon pairs", title="Scer") + 
  theme_classic() + 
  theme(
    text = element_text(size=48), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.28, 0.85)
  )
dev.off()


### cele ###

cele_pair <- read.table("data/processed/cele_pairdensities.txt", header=T)
cele_deltaDF <- data.frame(dYO=cele_pair$RD_young-cele_pair$RD_old, dYM=cele_pair$RD_young-cele_pair$RD_middle)


svg("figures/figure1/cele_pair_delta.svg", height=5, width=8)
ggplot(cele_deltaDF) + 
  geom_vline(aes(xintercept = 2*sd(dYO, na.rm=T)), size=1.2, color="#222222", linetype='dashed') +
  geom_vline(aes(xintercept =-2*sd(dYO, na.rm=T)), size=1.2, color="#222222", linetype='dashed') +
  geom_histogram(aes(x=dYO), binwidth=0.02, fill="#9826bb") + 
  scale_x_continuous(limits=c(-1, 1)) + 
  labs(x=paste("\u0394", "RD = RD(1d) - RD(12d)", sep=""), y="Count") + 
  theme_classic() + 
  theme(
    text = element_text(size=48)
  )
dev.off()


cele_extrDF <- data.frame(  dYO=c(sum(cele_deltaDF$dYO < -2*sd(cele_deltaDF$dYO, na.rm=T), na.rm=T), sum(cele_deltaDF$dYO > 2*sd(cele_deltaDF$dYO, na.rm=T), na.rm=T))/3723*100,
                            dYM=c(sum(cele_deltaDF$dYM < -2*sd(cele_deltaDF$dYM, na.rm=T), na.rm=T), sum(cele_deltaDF$dYM > 2*sd(cele_deltaDF$dYM, na.rm=T), na.rm=T))/3723*100, 
                            cat=c("up", "down") )

cele_extr.m <- melt(cele_extrDF, id='cat')
cele_extr.m$cat <- factor(cele_extr.m$cat, levels=c("up", "down"))
cele_extr.m$variable <- factor(cele_extr.m$variable, levels=c("dYM", "dYO"))
cele_extr.m$label <- factor(c("12d", "12d", "6d", "6d"), levels=c("6d", "12d") )


svg("figures/figure1/cele_pair_updown.svg", height=5, width=5)

ggplot(cele_extr.m) + 
  geom_col(aes(x=label, y=value, fill=cat), position=position_dodge2()) + 
  scale_fill_manual(values=c("#739bd0", "#0f4c81")) + 
  labs(x="", y="% codon pairs", title="Cele") + 
  theme_classic() + 
  theme(
    text = element_text(size=48), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.85)
  )
dev.off()


# --------------------- Figure 1E ---------------------
scer_kmer <- read.table("data/processed/kmer_scer_counts.txt", header=T, sep='\t')
scer_kmer <- scer_kmer[scer_kmer$count >= 20,]

#sel <- ((scer_kmer$class_young == "DT") + (scer_kmer$class_middle == "DT") + (scer_kmer$class_old == "DT")) > 0
sel <- ((scer_kmer$class_young == "DT") + (scer_kmer$class_old == "DT")) > 0
scer_kmer <- scer_kmer[sel,]




res <- matrix(NA, nrow=nrow(scer_kmer), ncol=2)
for (i in 1:nrow(scer_kmer)){

  a <- scer_kmer$hi_young[i]
  b <- scer_kmer$lo_young[i]
  c <- scer_kmer$hi_old[i]
  d <- scer_kmer$lo_old[i]
  
  
  arr <- array(c(a,b,c,d), dim=c(2,2))
  f <- fisher.test(arr)
  res[i,1] <- f$estimate
  res[i,2] <- f$p.val
}

scer_DT <- as.data.frame(res)
colnames(scer_DT) <- c("OR", "pval")
scer_DT$padj <- p.adjust(res[,2])
scer_DT$sig <- ifelse(scer_DT$padj < 0.05, "sig", "ns")
scer_DT$taille <- ifelse(scer_DT$padj < 0.05, 1.1, 1)
scer_DT$name <- rep("", nrow(scer_DT))
scer_DT$name[scer_DT$padj < 0.05] <- as.character(scer_kmer$kmer[scer_DT$padj < 0.05])
scer_DT$hjustify <- ifelse(scer_DT$padj < 0.05, scer_DT$OR/max(scer_DT$OR) , 0)


svg("figures/figure1/scer_DT_volcano.svg", height=5, width=5)

ggplot(scer_DT) + 
  geom_point(aes(x=OR, y=-log(pval), color=sig, size=taille), show.legend=F) + 
  geom_text(aes(x=OR, y=-log(pval), label=name), hjust=scer_DT$hjustify, nudge_y=1.5, size=5) + 
  labs(x="Odd's ratio", y="-log(pval)") + 
  scale_color_manual(values=c("#333333", "#bb2649")) + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.title = element_text(size=28)
  )

dev.off()



sig <- p.adjust(res[,2])<0.05
scer_kmer[sig,]
par(mfrow=c(2,1))
boxplot( scer_kmer$hi_young, scer_kmer$lo_young, scer_kmer$hi_old, scer_kmer$lo_old , ylim=c(0, 50), col=2)
boxplot( scer_kmer$hi_young[sig], scer_kmer$lo_young[sig], scer_kmer$hi_old[sig], scer_kmer$lo_old[sig], ylim=c(0, 100), col=2)


scer_kmerDF <- rbind( data.frame(class=rep("hiRD_0d", nrow(scer_kmer)), count=scer_kmer$hi_young, selection=rep("all DT", nrow(scer_kmer))),
                      data.frame(class=rep("loRD_0d", nrow(scer_kmer)), count=scer_kmer$lo_young, selection=rep("all DT", nrow(scer_kmer))),
                      data.frame(class=rep("hiRD_4d", nrow(scer_kmer)), count=scer_kmer$hi_old, selection=rep("all DT", nrow(scer_kmer))),
                      data.frame(class=rep("loRD_4d", nrow(scer_kmer)), count=scer_kmer$lo_old, selection=rep("all DT", nrow(scer_kmer))),
                      data.frame(class=rep("hiRD_0d", sum(sig)), count=scer_kmer$hi_young[sig], selection=rep("polyQ", sum(sig))),
                      data.frame(class=rep("loRD_0d", sum(sig)), count=scer_kmer$lo_young[sig], selection=rep("polyQ", sum(sig))),
                      data.frame(class=rep("hiRD_4d", sum(sig)), count=scer_kmer$hi_old[sig], selection=rep("polyQ", sum(sig))),
                      data.frame(class=rep("loRD_4d", sum(sig)), count=scer_kmer$lo_old[sig], selection=rep("polyQ", sum(sig)))
                      )
scer_kmerDF$class <- factor(scer_kmerDF$class, levels=c("hiRD_0d", "hiRD_4d", "loRD_0d", "loRD_4d"))


svg("figures/figure1/scer_DT_shifts.svg", height=5, width=5)

ggplot(scer_kmerDF) + 
  geom_boxplot(aes(x=class, y=count, fill=selection)) + 
  scale_fill_manual(values=c("#777777", "#bb2649")) + 
  scale_y_continuous(limits=c(0, 80)) +
  labs(x="", y="Count") + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.text.x = element_text(angle=35, hjust=1),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.title = element_blank(), 
    legend.position = c(0.2, 0.9)
  )

dev.off()



# --------------------- Figure 1F ---------------------

riboseq_YBR212W <- as.matrix(read.table("data/processed/riboseq_YBR212W.txt"))
riboseq_YMR070W <- as.matrix(read.table("data/processed/riboseq_YMR070W.txt"))


avrg_profile3 <- function(matIN){
  tmp <- matrix(NA, nrow=3, ncol=ncol(matIN))

  tmp[1,] <- colMeans(matIN[c(1:4),], na.rm=T)
  tmp[2,] <- colMeans(matIN[c(5:8),], na.rm=T)
  tmp[3,] <- colMeans(matIN[c(9:12),], na.rm=T)
  
  tmp
}

profile_YBR212W <- avrg_profile3(riboseq_YBR212W)
YBR212W_df <- as.data.frame(t(profile_YBR212W))
colnames(YBR212W_df) <- c("0d", "2d", "4d")
YBR212W_df$pos <- c(1: nrow(YBR212W_df) )
YBR212W_df.m <- melt(YBR212W_df, id='pos')

polydf <- data.frame(x=c(483, 490, 490, 483), y=c(0, 0, 5, 5), class=rep("polyQ", 4))


svg("figures/figure1/polyQ_YBR212W.svg", height=5, width=8)

ggplot(YBR212W_df.m) + 
  geom_polygon(data=polydf, inherit.aes=F, aes(x=x, y=y), show.legend= FALSE, alpha=0.2, fill="#bb2649") + 
  geom_line(aes(x=pos, y=value, col=variable, linetype=variable), size=1.2 ) + 
  scale_x_continuous(limits=c(450, 500)) + 
  scale_y_continuous(limits=c(0, 5)) + 
  #scale_color_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  scale_color_manual(values=c("#228b22", "#999999", "#222222")) + 
  labs(x="Position", y="Ribosome density", title="NGR1") + 
  theme_classic() + 
  theme(
    text = element_text(size=48),
    title = element_text(size=24),
    axis.title = element_text(size=40),
    legend.position = 'none'
  )

dev.off()


###

profile_YMR070W <- avrg_profile3(riboseq_YMR070W)
YMR070W_df <- as.data.frame(t(profile_YMR070W))
colnames(YMR070W_df) <- c("0d", "2d", "4d")
YMR070W_df$pos <- c(1: nrow(YMR070W_df) )
YMR070W_df.m <- melt(YMR070W_df, id='pos')

polydf1 <- data.frame(x=c(8, 11, 11, 8), y=c(0, 0, 5, 5), class=rep("Q1", 4))


svg("figures/figure1/polyQ_YMR070W.svg", height=5, width=8)

ggplot(YMR070W_df.m) + 
  geom_polygon(data=polydf1, inherit.aes=F, aes(x=x, y=y), show.legend= FALSE, alpha=0.2, fill="#bb2649") + 
  geom_line(aes(x=pos, y=value, col=variable, linetype=variable), size=1.2 ) + 
  scale_x_continuous(limits=c(1, 50)) + 
  scale_y_continuous(limits=c(0, 5)) + 
  scale_color_manual(values=c("#228b22", "#999999", "#222222")) + 
  labs(x="Position", y="Ribosome density", title="MOT3") + 
  theme_classic() + 
  theme(
    text = element_text(size=48),
    title = element_text(size=24),
    axis.title = element_text(size=40),
    legend.position = 'none'
  )

dev.off()