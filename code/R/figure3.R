# RiboAge | Figure 3
# S Pechmann (2025)


library(ggplot2)
library(reshape2)
library(cowplot)

setwd("~/M2/RiboAge")




# --------------------- Figure 3A ---------------------

scer_tRNA15_young <- scan('data/processed/scer2_tRNA_15_young_CC.txt')
scer_tRNA15_middle <- scan('data/processed/scer2_tRNA_15_middle_CC.txt')
scer_tRNA15_old <- scan('data/processed/scer2_tRNA_15_old_CC.txt')

scer_tRNA2_young <- scan('data/processed/scer2_tRNA_2_young_CC.txt')
scer_tRNA2_middle <- scan('data/processed/scer2_tRNA_2_middle_CC.txt')
scer_tRNA2_old <- scan('data/processed/scer2_tRNA_2_old_CC.txt')

scer_tRNA3_young <- scan('data/processed/scer2_tRNA_3_young_CC.txt')
scer_tRNA3_middle <- scan('data/processed/scer2_tRNA_3_middle_CC.txt')
scer_tRNA3_old <- scan('data/processed/scer2_tRNA_3_old_CC.txt')



scer_tRNA <- rbind(data.frame(age=0, opt=max(scer_tRNA15_young), t10=as.numeric(quantile(scer_tRNA15_young, 0.9)), fc='FC<1.5'),
                   data.frame(age=2, opt=max(scer_tRNA15_middle), t10=as.numeric(quantile(scer_tRNA15_middle, 0.9)), fc='FC<1.5'),
                   data.frame(age=4, opt=max(scer_tRNA15_old), t10=as.numeric(quantile(scer_tRNA15_old, 0.9)), fc='FC<1.5'),
                   data.frame(age=0, opt=max(scer_tRNA2_young), t10=as.numeric(quantile(scer_tRNA2_young, 0.9)), fc='FC<2'),
                   data.frame(age=2, opt=max(scer_tRNA2_middle), t10=as.numeric(quantile(scer_tRNA2_middle, 0.9)), fc='FC<2'),
                   data.frame(age=4, opt=max(scer_tRNA2_old), t10=as.numeric(quantile(scer_tRNA2_old, 0.9)), fc='FC<2'),
                   data.frame(age=0, opt=max(scer_tRNA3_young), t10=as.numeric(quantile(scer_tRNA3_young, 0.9)), fc='FC<3'),
                   data.frame(age=2, opt=max(scer_tRNA3_middle), t10=as.numeric(quantile(scer_tRNA3_middle, 0.9)), fc='FC<3'),
                   data.frame(age=4, opt=max(scer_tRNA3_old), t10=as.numeric(quantile(scer_tRNA3_old, 0.9)), fc='FC<3')
                  )
scer_tRNA$fc <- factor(scer_tRNA$fc, levels=c("FC<3", "FC<2", "FC<1.5"))



svg("figures/figure3/scer_tRNA_top10.svg", height=5, width=3.5)
ggplot(scer_tRNA) + 
  #geom_hline(aes(yintercept = opt[1]), size=0.5, color="#222222", alpha=0.5, linetype=5) +
  geom_line(aes(x=age, y=opt, color=fc), size=2) + 
  geom_point(aes(x=age, y=opt, color=fc), size=5) +
  geom_ribbon((aes(x=age, ymin=t10, ymax=opt, fill=fc)), alpha=0.5, show.legend=F) + 
  scale_color_manual(values=c("#186218", "#228b22", "#2cb42c")) +
  scale_fill_manual(values=c("#186218", "#228b22", "#2cb42c")) + 
  scale_x_continuous(breaks=c(0,2,4), labels=c(0,2,4)) +
  scale_y_continuous(limits=c(0.4,0.9)) + 
  labs(x="Age", y="Model fit (correlation)") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    legend.title = element_blank(),
    legend.position = c(0.7, 0.2)
    )

dev.off()







# --------------------- Figure 3B ---------------------

scer_tRNA15_young.tRNA <- read.table('data/processed/scer2_tRNA_15_young_tRNA.txt', header=T)
scer_tRNA15_middle.tRNA <- read.table('data/processed/scer2_tRNA_15_middle_tRNA.txt', header=T)
scer_tRNA15_old.tRNA <- read.table('data/processed/scer2_tRNA_15_old_tRNA.txt', header=T)

scer_tRNA2_young.tRNA <- read.table('data/processed/scer2_tRNA_2_young_tRNA.txt', header=T)
scer_tRNA2_middle.tRNA <- read.table('data/processed/scer2_tRNA_2_middle_tRNA.txt', header=T)
scer_tRNA2_old.tRNA <- read.table('data/processed/scer2_tRNA_2_old_tRNA.txt', header=T)

scer_tRNA3_young.tRNA <- read.table('data/processed/scer2_tRNA_3_young_tRNA.txt', header=T)
scer_tRNA3_middle.tRNA <- read.table('data/processed/scer2_tRNA_3_middle_tRNA.txt', header=T)
scer_tRNA3_old.tRNA <- read.table('data/processed/scer2_tRNA_3_old_tRNA.txt', header=T)



N <- 20

sortedmean <- function(trna, corr, ref, N){
  data <- trna[,sort(corr, decreasing=T, index.return=T)$ix[1:N]]
  m <- rowMeans(data/ref)
  m <- m[sort(m, decreasing=T, index.return=T)$ix]
  m
}


trnaDF <- read.csv("data/tRNA/yeast_tRNA_count.csv")

scer_trnaDF <- data.frame(ix = c(1:nrow(scer_tRNA15_young.tRNA)), 
                          young_fc15=sortedmean(scer_tRNA15_young.tRNA, scer_tRNA15_young, trnaDF$Abundance, 20),
                          young_fc2=sortedmean(scer_tRNA2_young.tRNA, scer_tRNA2_young, trnaDF$Abundance, 20),
                          young_fc3=sortedmean(scer_tRNA3_young.tRNA, scer_tRNA3_young, trnaDF$Abundance, 20),
                          middle_fc15=sortedmean(scer_tRNA15_middle.tRNA, scer_tRNA15_middle, trnaDF$Abundance, 20),
                          middle_fc2=sortedmean(scer_tRNA2_middle.tRNA, scer_tRNA2_middle, trnaDF$Abundance, 20),
                          middle_fc3=sortedmean(scer_tRNA3_middle.tRNA, scer_tRNA3_middle, trnaDF$Abundance, 20),
                          old_fc15=sortedmean(scer_tRNA15_old.tRNA, scer_tRNA15_old, trnaDF$Abundance, 20),
                          old_fc2=sortedmean(scer_tRNA2_old.tRNA, scer_tRNA2_old, trnaDF$Abundance, 20),
                          old_fc3=sortedmean(scer_tRNA3_old.tRNA, scer_tRNA3_old, trnaDF$Abundance, 20)
                          )
  

scer_trnaDF <- scer_trnaDF[trnaDF$Neigh_Nucl!="N",]   # without stop anticodons


svg("figures/figure3/scer_tRNA_fc.svg", height=5, width=2.2)

ggplot(scer_trnaDF) + 
  geom_hline(aes(yintercept = 1), size=0.7, color="#222222", alpha=0.5 ) +
  geom_hline(aes(yintercept = 2), size=0.5, color="#222222", alpha=0.5, linetype=2 ) +
  geom_hline(aes(yintercept = 1/2), size=0.5, color="#222222", alpha=0.5, linetype=2 ) +
  geom_line(aes(x=ix, y=young_fc2), size=1.2, col="#2cb42c") + 
  geom_line(aes(x=ix, y=middle_fc2), size=1.2, col="#228b22") + 
  geom_line(aes(x=ix, y=old_fc2), size=1.2, col="#186218") + 
  labs(y="FC", x="tRNAs") + 
  scale_y_log10(breaks=c(0.5, 1, 2), labels=c(0.5, 1, 2)) + 
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  )

dev.off()



# --------------------- Figure 3C ---------------------

acDF <- as.data.frame(read.table("data/tRNA/anticodon2.txt", header=T))
trnaDF <- trnaDF[trnaDF$Neigh_Nucl!="N",]   # without stop anticodons

res.YO <- data.frame(rd=c(), trna=c(), class=c())
res.YM <- data.frame(rd=c(), trna=c(), class=c())
for (i in 1:nrow(trnaDF)){
  current_ac <- as.character(trnaDF$Anticodon[i])
  current_co <- as.character(acDF$codon[which(acDF$anticodon == current_ac)])
  current_rd <- as.numeric(scer_rd[which(scer_rd$codon==current_co),c(3:5)])
  current_fc <- scer_trnaDF[i,]
  rel_rd.YO <- current_rd[3]/current_rd[1]
  rel_rd.YM <- current_rd[2]/current_rd[1]
  rel_trna.YO <- current_fc$old_fc2 / current_fc$young_fc2
  rel_trna.YM <- current_fc$middle_fc2 / current_fc$young_fc2
  res.YO <- rbind(res.YO, data.frame(rd=rel_rd.YO, trna=rel_trna.YO, class="4d / 0d"))
  res.YM <- rbind(res.YM, data.frame(rd=rel_rd.YM, trna=rel_trna.YM, class="2d / 0d"))
  rd_trna.DF <- rbind(res.YO, res.YM)
}


R1 <- cor(rd_trna.DF$rd[rd_trna.DF$class=="4d / 0d"], rd_trna.DF$trna[rd_trna.DF$class=="4d / 0d"])
R2 <- cor(rd_trna.DF$rd[rd_trna.DF$class=="2d / 0d"], rd_trna.DF$trna[rd_trna.DF$class=="2d / 0d"])



svg("figures/figure3/scer_rel_rd_tRNA.svg", height=5, width=4)

ggplot(rd_trna.DF, aes(x=rd, y=trna, color=class)) + 
  geom_smooth(method='glm', se=F) + 
  geom_point(aes(shape=class), size=2) + 
  labs(x="Relative RD", y="Relative [tRNA*]") + 
  scale_color_manual(values=c("#186218", "#2cb42c")) + 
  scale_x_continuous(breaks=c(0.8, 1, 1.2), labels=c(0.8, 1, 1.2)) + 
  scale_y_continuous(limits=c(0.83, 1.37)) + 
  scale_shape_manual(values=c(19, 15))+ 
  geom_text(x=0.75, y=0.85, hjust=0, label=paste("R =", round(R1, 2), sep=' '), col="#186218", size=5, stat="unique", show.legend=F) + 
  geom_text(x=1, y=0.85, hjust=0, label=paste("R =", round(R2, 2), sep=' '), col="#2cb42c", size=5, stat="unique", show.legend=F) + 
  
  theme_classic() + 
  theme(
    text= element_text(size=28), 
    legend.title = element_blank(), 
    legend.position = 'top'
  )

dev.off()









# --------------------- Figure 3D ---------------------

cele_tRNA15_young <- scan('data/processed/cele2_tRNA_15_young_CC.txt')
cele_tRNA15_middle <- scan('data/processed/cele2_tRNA_15_middle_CC.txt')
cele_tRNA15_old <- scan('data/processed/cele2_tRNA_15_old_CC.txt')

cele_tRNA2_young <- scan('data/processed/cele2_tRNA_2_young_CC.txt')
cele_tRNA2_middle <- scan('data/processed/cele2_tRNA_2_middle_CC.txt')
cele_tRNA2_old <- scan('data/processed/cele2_tRNA_2_old_CC.txt')

cele_tRNA3_young <- scan('data/processed/cele2_tRNA_3_young_CC.txt')
cele_tRNA3_middle <- scan('data/processed/cele2_tRNA_3_middle_CC.txt')
cele_tRNA3_old <- scan('data/processed/cele2_tRNA_3_old_CC.txt')

cele_tRNA <- rbind(data.frame(age=0, opt=max(cele_tRNA15_young), t10=as.numeric(quantile(cele_tRNA15_young, 0.9)), fc='FC<1.5'),
                   data.frame(age=6, opt=max(cele_tRNA15_middle), t10=as.numeric(quantile(cele_tRNA15_middle, 0.9)), fc='FC<1.5'),
                   data.frame(age=12, opt=max(cele_tRNA15_old), t10=as.numeric(quantile(cele_tRNA15_old, 0.9)), fc='FC<1.5'),
                   data.frame(age=0, opt=max(cele_tRNA2_young), t10=as.numeric(quantile(cele_tRNA2_young, 0.9)), fc='FC<2'),
                   data.frame(age=6, opt=max(cele_tRNA2_middle), t10=as.numeric(quantile(cele_tRNA2_middle, 0.9)), fc='FC<2'),
                   data.frame(age=12, opt=max(cele_tRNA2_old), t10=as.numeric(quantile(cele_tRNA2_old, 0.9)), fc='FC<2'),
                   data.frame(age=0, opt=max(cele_tRNA3_young), t10=as.numeric(quantile(cele_tRNA3_young, 0.9)), fc='FC<3'),
                   data.frame(age=6, opt=max(cele_tRNA3_middle), t10=as.numeric(quantile(cele_tRNA3_middle, 0.9)), fc='FC<3'),
                   data.frame(age=12, opt=max(cele_tRNA3_old), t10=as.numeric(quantile(cele_tRNA3_old, 0.9)), fc='FC<3')
)
cele_tRNA$fc <- factor(cele_tRNA$fc, levels=c("FC<3", "FC<2", "FC<1.5"))



svg("figures/figure3/cele_tRNA_top10.svg", height=5, width=3.5)
ggplot(cele_tRNA) + 
  #geom_hline(aes(yintercept = opt[1]), size=0.5, color="#222222", alpha=0.5, linetype=5) +
  geom_line(aes(x=age, y=opt, color=fc), size=2) + 
  geom_point(aes(x=age, y=opt, color=fc), size=5) +
  geom_ribbon((aes(x=age, ymin=t10, ymax=opt, fill=fc)), alpha=0.5, show.legend=F) + 
  scale_color_manual(values=c("#761d91", "#9826bb", "#b33dd7")) +
  scale_fill_manual(values=c("#761d91", "#9826bb", "#b33dd7")) + 
  scale_x_continuous(breaks=c(0,6,12), labels=c(1,6,12)) +
  scale_y_continuous(limits=c(0.4,0.92)) + 
  labs(x="Age", y="Model fit (correlation)") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    legend.title = element_blank(),
    legend.position = c(0.7, 0.2)
  )

dev.off()







# --------------------- Figure 3E ---------------------

cele_tRNA15_young.tRNA <- read.table('data/processed/cele2_tRNA_15_young_tRNA.txt', header=T)
cele_tRNA15_middle.tRNA <- read.table('data/processed/cele2_tRNA_15_middle_tRNA.txt', header=T)
cele_tRNA15_old.tRNA <- read.table('data/processed/cele2_tRNA_15_old_tRNA.txt', header=T)

cele_tRNA2_young.tRNA <- read.table('data/processed/cele2_tRNA_2_young_tRNA.txt', header=T)
cele_tRNA2_middle.tRNA <- read.table('data/processed/cele2_tRNA_2_middle_tRNA.txt', header=T)
cele_tRNA2_old.tRNA <- read.table('data/processed/cele2_tRNA_2_old_tRNA.txt', header=T)

cele_tRNA3_young.tRNA <- read.table('data/processed/cele2_tRNA_3_young_tRNA.txt', header=T)
cele_tRNA3_middle.tRNA <- read.table('data/processed/cele2_tRNA_3_middle_tRNA.txt', header=T)
cele_tRNA3_old.tRNA <- read.table('data/processed/cele2_tRNA_3_old_tRNA.txt', header=T)


N <- 20
sortedmean <- function(trna, corr, ref, N){
  data <- trna[,sort(corr, decreasing=T, index.return=T)$ix[1:N]]
  m <- rowMeans(data/ref)
  m <- m[sort(m, decreasing=T, index.return=T)$ix]
  m
}


trnaDF <- read.table("data/tRNA/cele_tRNA_count.txt", sep=" ", header=T)
trnaDF <- trnaDF[trnaDF$Abundance>0,]

cele_trnaDF <- data.frame(ix = c(1:nrow(cele_tRNA15_young.tRNA)), 
                          young_fc15=sortedmean(cele_tRNA15_young.tRNA, cele_tRNA15_young, trnaDF$Abundance, 20),
                          young_fc2=sortedmean(cele_tRNA2_young.tRNA, cele_tRNA2_young, trnaDF$Abundance, 20),
                          young_fc3=sortedmean(cele_tRNA3_young.tRNA, cele_tRNA3_young, trnaDF$Abundance, 20),
                          middle_fc15=sortedmean(cele_tRNA15_middle.tRNA, cele_tRNA15_middle, trnaDF$Abundance, 20),
                          middle_fc2=sortedmean(cele_tRNA2_middle.tRNA, cele_tRNA2_middle, trnaDF$Abundance, 20),
                          middle_fc3=sortedmean(cele_tRNA3_middle.tRNA, cele_tRNA3_middle, trnaDF$Abundance, 20),
                          old_fc15=sortedmean(cele_tRNA15_old.tRNA, cele_tRNA15_old, trnaDF$Abundance, 20),
                          old_fc2=sortedmean(cele_tRNA2_old.tRNA, cele_tRNA2_old, trnaDF$Abundance, 20),
                          old_fc3=sortedmean(cele_tRNA3_old.tRNA, cele_tRNA3_old, trnaDF$Abundance, 20)
)




svg("figures/figure3/cele_tRNA_fc.svg", height=5, width=2.2)

ggplot(cele_trnaDF) + 
  geom_hline(aes(yintercept = 1), size=0.7, color="#222222", alpha=0.5 ) +
  geom_hline(aes(yintercept = 2), size=0.5, color="#222222", alpha=0.5, linetype=2 ) +
  geom_hline(aes(yintercept = 1/2), size=0.5, color="#222222", alpha=0.5, linetype=2 ) +
  geom_line(aes(x=ix, y=young_fc2), size=1.2, col="#b33dd7") + 
  geom_line(aes(x=ix, y=middle_fc2), size=1.2, col="#9826bb") + 
  geom_line(aes(x=ix, y=old_fc2), size=1.2, col="#761d91") + 
  labs(y="FC", x="tRNAs") + 
  scale_y_log10(breaks=c(0.5, 1, 2), labels=c(0.5, 1, 2)) + 
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  )

dev.off()





acDF <- as.data.frame(read.table("data/tRNA/anticodon2.txt", header=T))

res.YO <- data.frame(rd=c(), trna=c(), class=c())
res.YM <- data.frame(rd=c(), trna=c(), class=c())
for (i in 1:nrow(trnaDF)){
  current_ac <- as.character(trnaDF$Anticodon[i])
  current_co <- as.character(acDF$codon[which(acDF$anticodon == current_ac)])
  current_rd <- as.numeric(cele_rd[which(cele_rd$codon==current_co),c(3:5)])
  current_fc <- cele_trnaDF[i,]
  rel_rd.YO <- current_rd[3]/current_rd[1]
  rel_rd.YM <- current_rd[2]/current_rd[1]
  rel_trna.YO <- current_fc$old_fc2 / current_fc$young_fc2
  rel_trna.YM <- current_fc$middle_fc2 / current_fc$young_fc2
  res.YO <- rbind(res.YO, data.frame(rd=rel_rd.YO, trna=rel_trna.YO, class="12d / 1d"))
  res.YM <- rbind(res.YM, data.frame(rd=rel_rd.YM, trna=rel_trna.YM, class="6d / 1d"))
  rd_trna.DF <- rbind(res.YO, res.YM)
}


R1 <- cor(rd_trna.DF$rd[rd_trna.DF$class=="12d / 1d"], rd_trna.DF$trna[rd_trna.DF$class=="12d / 1d"])
R2 <- cor(rd_trna.DF$rd[rd_trna.DF$class=="6d / 1d"], rd_trna.DF$trna[rd_trna.DF$class=="6d / 1d"])



svg("figures/figure3/cele_rel_rd_tRNA.svg", height=5, width=4)

ggplot(rd_trna.DF, aes(x=rd, y=trna, color=class)) + 
  geom_smooth(method='glm', se=F) + 
  geom_point(aes(shape=class), size=2) + 
  labs(x="Relative RD", y="Relative [tRNA*]") + 
  scale_color_manual(values=c("#761d91", "#b33dd7")) + 
  scale_shape_manual(values=c(19, 15))+ 
  geom_text(x=0.7, y=0.75, hjust=0, label=paste("R =", round(R1, 2), sep=' '), col="#761d91", size=5, stat="unique", show.legend=F) + 
  geom_text(x=0.9, y=0.75, hjust=0, label=paste("R =", round(R2, 2), sep=' '), col="#b33dd7", size=5, stat="unique", show.legend=F) + 
  #scale_x_continuous(breaks=c(0.8, 1, 1.2), labels=c(0.8, 1, 1.2)) + 
  scale_y_continuous(limits=c(0.75, 1.3)) + 
  theme_classic() + 
  theme(
    text= element_text(size=28), 
    legend.title = element_blank(), 
    legend.position = 'top'
  )

dev.off()





# --------------------- Figure 3F ---------------------

scer_abund <- read.table("data/processed/scer_translation_abundance.txt", header=T)
cele_abund <- read.table("data/processed/cele_translation_abundance.txt", header=T)

scer_abund$species <- rep("Scer", nrow(scer_abund))
cele_abund$species <- rep("Cele", nrow(cele_abund))

transl <- rbind( scer_abund[,c("species", "exp_young", "exp_middle", "exp_old")], cele_abund[,c("species", "exp_young", "exp_middle", "exp_old")] )
transl.m <- melt(transl, id="species")
transl.m$species <- factor(transl.m$species, levels=c("Scer", "Cele"))
transl.m$variable <- ifelse(transl.m$species=="Cele", paste(transl.m$variable, "C", sep=""), paste(transl.m$variable) )
transl.m$variable <- factor(transl.m$variable, levels=c("exp_young", "exp_middle", "exp_old", "exp_youngC", "exp_middleC", "exp_oldC"))


svg("figures/figure3/transl_abundance.svg", height=5, width=3.5)

ggplot(transl.m) + 
  geom_boxplot(aes(x=species, y=value, fill=variable)) + 
  scale_y_log10(limits=c(0.8, 10000)) + 
  scale_fill_manual(values=c("#2cb42c","#228b22", "#186218",  "#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="Expression") + 
  theme_classic() + 
  theme(
    text = element_text(size=24), 
    axis.title = element_text(size=28),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position = 'none'
  )


dev.off()


# --------------------- Figure 3G ---------------------


scer_abund <- read.table("data/processed/scer_translation_abundance.txt", header=T)

new_labels <- rep(NA, nrow(scer_abund))
new_labels[scer_abund$category=="Ribosome_biogenesis_factor"] <- "Ribogenesis"
new_labels[scer_abund$category=="Ribosome"] <- "Ribosome"
new_labels[scer_abund$category=="Translation_initiation"] <- "Initiation"
new_labels[scer_abund$category=="Translation_elongation"] <- "Elongation"
new_labels[scer_abund$category=="Translation_termination"] <- "Termination"
new_labels <- factor(new_labels, levels=c("Termination", "Elongation", "Initiation", "Ribogenesis", "Ribosome"))
scer_abund$category <- new_labels

norm_scer.middle <- median( scer_abund$exp_middle[scer_abund$category=="Ribosome"]/scer_abund$exp_young[scer_abund$category=="Ribosome"] )
norm_scer.old <- median( scer_abund$exp_old[scer_abund$category=="Ribosome"]/scer_abund$exp_young[scer_abund$category=="Ribosome"] )
scer_abund$expN_middle <- scer_abund$exp_middle / norm_scer.middle
scer_abund$expN_old <- scer_abund$exp_old / norm_scer.old

scer_abund_cat <- scer_abund[,c(5,2,7,8)]
colnames(scer_abund_cat) <- c('category', '0d', '2d', '4d')
scer_abund_cat.m <- melt(scer_abund_cat, id='category')
scer_abund_cat.m$variable <- factor(scer_abund_cat.m$variable, levels=c('4d', '2d', '0d'))


svg("figures/figure3/scer_abundance.svg", height=5, width=5)

ggplot(scer_abund_cat.m) + 
  geom_boxplot(aes(x=category, y=value, fill=variable), show.legend = T) + 
  labs(x="", y="Rel. expression") + 
  scale_fill_manual(values=c("#186218", "#228b22", "#2cb42c")) +
  scale_y_log10() + 
  coord_flip() + 
  #facet_grid(rows = vars(variable), scales="free_y") + 
  theme_classic() + 
  theme(
    text = element_text(size=24), 
    axis.title = element_text(size=28),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=20),
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(),
    legend.title = element_blank(), 
    legend.position = c(0.85, 0.15)
  )

dev.off()




scer_abund_elongation <- scer_abund[scer_abund$category=="Elongation",]


svg("figures/figure3/scer_elongation.svg", height=5, width=5)

ggplot(scer_abund_elongation ) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.3, color="#777777", linetype="twodash" ) +
  geom_point(aes(x=log10(exp_young), y=log10(expN_middle)), col="#2cb42c", size=2.5, shape=15, alpha=0.7) + 
  geom_point(aes(x=log10(exp_young), y=log10(expN_old)), col="#186218", size=2.5) + 
  scale_x_continuous(limits=c(0, 4.5), breaks=c(1, 2, 3,4), labels=c( 10, expression(10^2),  expression(10^3), expression(10^4))) + 
  scale_y_continuous(limits=c(0, 4.5), breaks=c(1, 2, 3,4), labels=c(10, expression(10^2),  expression(10^3), expression(10^4))) + 
  labs(x="Expression (0d)", y="Rel. expression (aged)", title="Translation elongation") + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    title = element_text(size=16),
    axis.text.y = element_text(size=24),
    axis.text.x = element_text(size=24, vjust=0),
    axis.title = element_text(size=30)
  )

dev.off()



# --------------------- Figure 3H ---------------------

cele_abund <- read.table("data/processed/cele_translation_abundance.txt", header=T)

new_labels <- rep(NA, nrow(cele_abund))
new_labels[cele_abund$category=="Ribosome_biogenesis_factor"] <- "Ribogenesis"
new_labels[cele_abund$category=="Ribosome"] <- "Ribosome"
new_labels[cele_abund$category=="Translation_initiation"] <- "Initiation"
new_labels[cele_abund$category=="Translation_elongation"] <- "Elongation"
new_labels[cele_abund$category=="Translation_termination"] <- "Termination"
new_labels <- factor(new_labels, levels=c("Termination", "Elongation", "Initiation", "Ribogenesis", "Ribosome"))
cele_abund$category <- new_labels

norm_cele.middle <- median( cele_abund$exp_middle[cele_abund$category=="Ribosome"]/cele_abund$exp_young[cele_abund$category=="Ribosome"] )
norm_cele.old <- median( cele_abund$exp_old[cele_abund$category=="Ribosome"]/cele_abund$exp_young[cele_abund$category=="Ribosome"] )
cele_abund$expN_middle <- cele_abund$exp_middle / norm_cele.middle
cele_abund$expN_old <- cele_abund$exp_old / norm_cele.old

cele_abund_cat <- cele_abund[,c(5,2,7,8)]
colnames(cele_abund_cat) <- c('category', '1d', '6d', '12d')
cele_abund_cat.m <- melt(cele_abund_cat, id='category')
cele_abund_cat.m$variable <- factor(cele_abund_cat.m$variable, levels=c('12d', '6d', '1d'))


svg("figures/figure3/cele_abundance.svg", height=5, width=5)

ggplot(cele_abund_cat.m) + 
  geom_boxplot(aes(x=category, y=value, fill=variable), show.legend = T) + 
  labs(x="", y="Rel. expression") + 
  scale_fill_manual(values=c("#761d91", "#9826bb", "#b33dd7")) +
  scale_y_log10(limits=c(0.9, 20000)) + 
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size=24), 
    axis.title = element_text(size=28),
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=20),
    axis.line.y = element_blank(), 
    axis.ticks.y = element_blank(),
    legend.title = element_blank(), 
    legend.position = c(0.85, 0.15)
  )

dev.off()



cele_abund <- read.table("data/processed/cele_translation_abundance.txt", header=T)

norm_cele.middle <- median( cele_abund$exp_middle[cele_abund$category=="Ribosome"]/cele_abund$exp_young[cele_abund$category=="Ribosome"] )
norm_cele.old <- median( cele_abund$exp_old[cele_abund$category=="Ribosome"]/cele_abund$exp_young[cele_abund$category=="Ribosome"] )
cele_abund$expN_middle <- cele_abund$exp_middle / norm_cele.middle
cele_abund$expN_old <- cele_abund$exp_old / norm_cele.old

cele_abund_elongation <- cele_abund[cele_abund$category=="Translation_elongation",]


svg("figures/figure3/cele_elongation.svg", height=5, width=5)

ggplot(cele_abund_elongation ) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.3, color="#777777", linetype="twodash" ) +
  geom_point(aes(x=log10(exp_young), y=log10(expN_middle)), col="#b33dd7", size=2.5, shape=15, alpha=0.7) + 
  geom_point(aes(x=log10(exp_young), y=log10(expN_old)), col="#761d91", size=2.5) + 
  scale_x_continuous(limits=c(-1.5, 4), breaks=c(0, 1, 2, 3,4), labels=c(1, 10, expression(10^2),  expression(10^3), expression(10^4))) + 
  scale_y_continuous(limits=c(-1.5, 4), breaks=c(0, 1, 2, 3,4), labels=c(1, 10, expression(10^2),  expression(10^3), expression(10^4))) + 
  labs(x="Expression (1d)", y="Rel. expression (aged)", title="Translation elongation") + 
  theme_classic() +
  theme(
    text = element_text(size=28), 
    title = element_text(size=16),
    axis.text.y = element_text(size=24),
    axis.text.x = element_text(size=24, vjust=0),
    axis.title = element_text(size=30)
  )

dev.off()
