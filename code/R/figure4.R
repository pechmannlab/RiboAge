# RiboAge | Figure 4
# S Pechmann (2025)


library(ggplot2)
library(reshape2)
library(cowplot)
library(uwot)
library(EnvStats)


setwd("~/M2/RiboAge")



# --------------------- Figure 4A ---------------------
scer_rd <- read.table("data/processed/scer_ribodensities.txt", header=T)
scer_codonusage = as.data.frame(read.table("data/processed/scer_codonusage.txt", header=T, sep='\t'))

scer_rd_cu.DF <- rbind(
  data.frame(rd = ( scer_rd$RD_old / scer_rd$RD_young ),
             cu = (scer_codonusage$old / scer_codonusage$young), 
             class = rep("4d / 0d", nrow(scer_rd)) ),
  data.frame(rd = ( scer_rd$RD_middle / scer_rd$RD_young ),
             cu = (scer_codonusage$middle / scer_codonusage$young), 
             class = rep("2d / 0d", nrow(scer_rd)) )
)

R1_scer = cor(scer_rd_cu.DF$rd[scer_rd_cu.DF$class=="4d / 0d"], scer_rd_cu.DF$cu[scer_rd_cu.DF$class=="4d / 0d"])
R2_scer = cor(scer_rd_cu.DF$rd[scer_rd_cu.DF$class=="2d / 0d"], scer_rd_cu.DF$cu[scer_rd_cu.DF$class=="2d / 0d"])

plot.scer <- ggplot(scer_rd_cu.DF, aes(x=rd, y=cu, color=class)) + 
  geom_smooth(method='glm', se=F) + 
  geom_point(aes(shape=class), size=2) + 
  labs(x="Relative RD", y="Relative CU", title="Scer") + 
  scale_color_manual(values=c("#186218", "#2cb42c")) + 
  geom_text(x=0.72, y=1.05, hjust=0, label=paste("R =", round(R1_scer, 2), sep=' '), col="#186218", size=5, stat="unique", show.legend=F) + 
  geom_text(x=0.72, y=1.0, hjust=0, label=paste("R =", round(R2_scer, 2), sep=' '), col="#2cb42c", size=5, stat="unique", show.legend=F) +   scale_shape_manual(values=c(19, 15))+ 
  theme_classic() + 
  theme(
    text = element_text(size=28) ,
    axis.title = element_text(size=28),
    title = element_text(size=16),
    legend.title = element_blank(), 
    legend.position = 'top'
  )




cele_rd <- read.table("data/processed/cele_ribodensities.txt", header=T)
cele_codonusage = as.data.frame(read.table("data/processed/cele_codonusage.txt", header=T, sep='\t'))

cele_rd_cu.DF <- rbind(data.frame(rd = (cele_rd$RD_old / cele_rd$RD_young),
                                  cu = (cele_codonusage$old / cele_codonusage$young), 
                                  class = rep("12d / 1d", nrow(cele_rd)) ),
                       data.frame(rd = (cele_rd$RD_middle / cele_rd$RD_young),
                                  cu = (cele_codonusage$middle / cele_codonusage$young), 
                                  class = rep("6d / 1d", nrow(cele_rd)) )
)


R1_cele = cor(cele_rd_cu.DF$rd[cele_rd_cu.DF$class=="12d / 1d"], cele_rd_cu.DF$cu[cele_rd_cu.DF$class=="12d / 1d"])
R2_cele = cor(cele_rd_cu.DF$rd[cele_rd_cu.DF$class=="6d / 1d"], cele_rd_cu.DF$cu[cele_rd_cu.DF$class=="6d / 1d"])

plot.cele <- ggplot(cele_rd_cu.DF, aes(x=rd, y=cu, color=class)) + 
  geom_smooth(method='glm', se=F) + 
  geom_point(aes(shape=class), size=2) + 
  labs(x="Relative RD", y="Relative CU", title="Cele") + 
  scale_color_manual(values=c("#761d91", "#b33dd7")) + 
  geom_text( x=0.65, y=1.14, hjust=0, label=paste("R =", round(R1_cele, 2), sep=' '), col="#761d91", size=5, stat="unique", show.legend=F) + 
  geom_text(x=0.65, y=1.1, hjust=0, label=paste("R =", round(R2_cele, 2), sep=' '), col="#b33dd7", size=5, stat="unique", show.legend=F) + 
  scale_shape_manual(values=c(19, 15))+ 
  theme_classic() + 
  theme(
    text = element_text(size=28) ,
    axis.title = element_text(size=28),
    title = element_text(size=16),
    legend.title = element_blank(), 
    legend.position = 'top'
  )


svg("figures/figure4/rel_rd_cu.svg", height=5, width=9)

plot_grid(plot.scer, plot.cele, ncol = 2, align = 'h')

dev.off()



# --------------------- Figure 4B ---------------------
scer_rd <- read.table("data/processed/scer_ribodensities.txt", header=T)
optDF <- read.table("data/processed/optimality_cu.txt", header=T)

p1 <- ggplot(optDF) + 
  geom_vline(aes(xintercept = 1), size=0.8, color="#222222", linetype=2) +
  geom_point(aes(x=CU_old/CU_young, y=tAI), color="#739bd0") + 
  labs(x="Relative CU (4d/0d)", y="tAI") + 
  theme_classic() + 
  theme(
    text = element_text(size=28)
  )


p2 <- ggplot(optDF) + 
  geom_vline(aes(xintercept = 1), size=0.8, color="#222222", linetype=2) +
  geom_point(aes(x=CU_old/CU_young, y=CSC), color="#0f4c81") + 
  labs(x="Relative CU (4d/0d)", y="CSC") + 
  theme_classic() + 
  theme(
    text = element_text(size=28)
  )


svg("figures/figure4/optimality_cu.svg", height=5, width=4)

plot_grid(p1, p2, ncol = 1, align = 'v')

dev.off()




scer_rd_tai.DF <- data.frame(rd = ( scer_rd$RD_old / scer_rd$RD_young ),
                            tAI = optDF$tAI, 
                            class = rep("4d / 0d", nrow(scer_rd)) )

plot.tAI <-  ggplot(scer_rd_tai.DF, aes(x=rd, y=tAI)) + 
  geom_smooth(method=lm, alpha=0.1, size=0.8, linetype=2, color="#777777", fill="#6667ab", se=T) + 
  geom_point(aes(shape=class), size=2, color="#739bd0") + 
  labs(x="Relative RD (4d/0d)", y="tAI") + 
  scale_x_continuous(breaks=c(0.8, 1, 1.2), labels=c(0.8, 1, 1.2)) + 
  #scale_color_manual(values=c("#186218", "#2cb42c")) + 
  theme_classic() + 
  theme(
    text = element_text(size=28) ,
    axis.title = element_text(size=28),
    title = element_text(size=16),
    legend.title = element_blank(), 
    legend.position = 'none'
  )


scer_rd_csc.DF <- data.frame(rd = ( scer_rd$RD_old / scer_rd$RD_young ),
             CSC = optDF$CSC, 
             class = rep("4d / 0d", nrow(scer_rd)) )

plot.CSC <- ggplot(scer_rd_csc.DF, aes(x=rd, y=CSC)) + 
  geom_smooth(method=lm, alpha=0.1, size=0.8, linetype=2, color="#777777", fill="#6667ab", se=T) + 
  geom_point(aes(shape=class), size=2, color="#0f4c81") + 
  labs(x="Relative RD (4d/0d)", y="CSC") + 
  scale_x_continuous(breaks=c(0.8, 1, 1.2), labels=c(0.8, 1, 1.2)) + 
  #scale_color_manual(values=c("#186218", "#2cb42c")) + 
  theme_classic() + 
  theme(
    text = element_text(size=28) ,
    axis.title = element_text(size=28),
    title = element_text(size=16),
    legend.title = element_blank(), 
    legend.position = 'none'
  )


svg("figures/figure4/optimality_RD.svg", height=5, width=4)

plot_grid(plot.tAI, plot.CSC, ncol = 1, align = 'v')

dev.off()



# --------------------- Figure 4C ---------------------

dwell <- read.table("data/SMoPT/codon_dwelltimes.txt", header=T)
smopt_young <- data.frame(smopt=dwell$young, RD=scer_rd$RD_young)

R_dwell = cor(scer_rd$RD_young, log(dwell$young) )

svg("figures/figure4/simDT_RD.svg", height=5, width=5)

ggplot(smopt_young, aes(x=RD, y=log(smopt) )) + 
  geom_smooth(method=lm, alpha=0.1, size=0.5, linetype=2, color="#777777", fill="#6667ab", se=F) + 
  geom_point(color="#bb2649", size=3) + 
  geom_text(x=0.6, y=1., hjust=0, label=paste("R =", round(R_dwell, 2), sep=' '), col="#bb2649", size=5, stat="unique", show.legend=F) + 
  labs(x="Experimental RD", y="Simulation dwelltimes (log)") + 
  theme_classic() + 
  theme(
    text = element_text(size=28)
  )

dev.off()




# --------------------- Figure 4D ---------------------

YAL003W <- read.table("data/SMoPT/YAL003W.txt", header=T)
YAL003W$collprob_young <- ifelse(YAL003W$elongation_young > 0, YAL003W$collision_young/YAL003W$elongation_young, 0)
YAL003W$collprob_old <- ifelse(YAL003W$elongation_old > 0, YAL003W$collision_old/YAL003W$elongation_old, 0)
YAL003W$idx <- c(1:nrow(YAL003W))
YAL003W$null <- rep(0, nrow(YAL003W))


p.young1 <- ggplot(YAL003W) + 
  geom_polygon(aes(x=idx, y=collprob_young ) , color="#739bd0", fill="#739bd0") + 
  labs(x="", y="", title="EFB1 - 0d") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.title = element_text(size=28),
    title = element_text(size=12),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
    )

p.old1 <- ggplot(YAL003W) + 
  geom_polygon(aes(x=idx, y=collprob_old) , color="#0f4c81", fill="#0f4c81" , size=1) + 
  labs(x="", y="", title="EFB1 - 4d") + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.title = element_text(size=28),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    title = element_text(size=12)
  )
  


YBR170C <- read.table("data/SMoPT/YBR170C.txt", header=T)
YBR170C$collprob_young <- ifelse(YBR170C$elongation_young > 0, YBR170C$collision_young/YBR170C$elongation_young, 0)
YBR170C$collprob_old <- ifelse(YBR170C$elongation_old > 0, YBR170C$collision_old/YBR170C$elongation_old, 0)
YBR170C$idx <- c(1:nrow(YBR170C))
YBR170C$null <- rep(0, nrow(YBR170C))


p.young2 <- ggplot(YBR170C) + 
  geom_polygon(aes(x=idx, y=collprob_young ) , color="#739bd0", fill="#739bd0") + 
  labs(x="", y="", title="NPL4 - 0d") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.title = element_text(size=28),
    title = element_text(size=12),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
  )

p.old2 <- ggplot(YBR170C) + 
  geom_polygon(aes(x=idx, y=collprob_old) , color="#0f4c81", fill="#0f4c81" , size=1) + 
  labs(x="", y="", title="NPL4 - 4d") + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.title = element_text(size=28),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    title = element_text(size=12)
  )


svg("figures/figure4/collision_examples.svg", height=5, width=3)

plot_grid(p.young1, p.old1, p.young2, p.old2, ncol = 1, align = 'v')

dev.off()



# ----------



collisions <- read.table('data/SMoPT/collisions_pergene.txt', header=T)
collisions$sig <- ifelse(collisions$p_adj < 0.05, "*", "ns")

svg("figures/figure4/collisions_volcano.svg", height=5, width=5)

ggplot(collisions, aes(x=log2(OR), y=-log(pval) )) + 
  geom_point(aes(color=sig), size=2, show.legend=F) + 
  scale_color_manual(values=c("#bb2649", "#777777")) + 
  labs(x="Odd's ratio (log2)", y="-log(p)") + 
  scale_y_continuous(limits=c(0, 80)) + 
  theme_classic() + 
  theme(
    text = element_text(size=48)
  )

dev.off()




svg("figures/figure4/collisions_hist.svg", height=5, width=5)

ggplot(collisions) + 
  geom_histogram(aes(x=log2(OR)), fill="#bb2649") + 
  labs(x="Odd's ratio (log2)", y="Count") + 
  theme_classic() + 
  theme(
    text = element_text(size=48)
  )

dev.off()




colldf <- read.table("data/SMoPT/collisions_counts.txt", header=T)
colnames(colldf) <- c("age", "rate", "cov.")
colldf.m <- melt(colldf, index='age')
colldf.m$age <- factor(colldf.m$age, levels=c('young', 'old'))


svg("figures/figure4/collisions_pct.svg", height=5, width=5)

ggplot(colldf.m) + 
  geom_col(aes(x=variable, y=value, fill=age), position=position_dodge2()) + 
  scale_fill_manual(values=c("#739bd0", "#0f4c81")) + 
  labs(x="", y="Percent") + 
  theme_classic() +
  theme(
    text = element_text(size=48), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

dev.off()




# --------------------- Figure 4E ---------------------

ccDF <- read.table("data/SMoPT/collisions_codoncounts.txt", header=T)
relCUdf <- data.frame(relCU=scer_codonusage$old/scer_codonusage$young,
                      relCOLL_site = (ccDF$old_site/sum(ccDF$old_site))/(ccDF$young_site/sum(ccDF$young_site)),
                      relCOLL_upstr= (ccDF$old_upstr/sum(ccDF$old_upstr))/(ccDF$young_upstr/sum(ccDF$young_upstr))                                                        
                      )

R_upstr <- cor(relCUdf$relCU, relCUdf$relCOLL_upstr)
R_site <- cor(relCUdf$relCU, relCUdf$relCOLL_site)


p.upstr <- ggplot(relCUdf, aes(x=relCU, y=relCOLL_upstr) ) + 
  geom_smooth(method=lm, alpha=0.1, size=0.8, linetype=2, color="#777777", fill="#777777", se=T) + 
  geom_point( color="#bb2649", size=2) + 
  geom_text(x=0.85, y=1.15, hjust=0, label=paste("R =", round(R_upstr, 2), sep=' '), col="#777777", size=5, stat="unique", show.legend=F) +   
  labs(x="Relative CU", y="Rel CU upstr.") + 
  theme_classic() + 
  theme(
    text = element_text(size=28)
  )


p.site <- ggplot(relCUdf, aes(x=relCU, y=relCOLL_site) ) + 
  geom_smooth(method=lm, alpha=0.1, size=0.8, linetype=2, color="#777777", fill="#777777", se=T) + 
  geom_point(color="#bb2649", size=2) + 
  geom_text(x=0.85, y=1.15, hjust=0, label=paste("R =", round(R_site, 2), sep=' '), col="#777777", size=5, stat="unique", show.legend=F) +  
  labs(x="Relative CU", y="Rel CU site") + 
  theme_classic() + 
  theme(
    text = element_text(size=28)
  )



rand.young <- as.matrix(read.table("data/SMoPT/rand_young.txt"))
rand.old <- as.matrix(read.table("data/SMoPT/rand_old.txt"))

res <- rep(NA, 1000)
for (i in 1:1000){
  y <- rand.young[,i]/sum(rand.young[,i])
  o <- rand.old[,i]/sum(rand.old[,i])
  rCU <- scer_codonusage$old/scer_codonusage$young
  rCOLL <- o/y
  res[i] <- cor(rCOLL, rCU)
}
randDF <- data.frame(rand=res)


p.rand <- ggplot(randDF) + 
  geom_histogram(aes(x=rand), alpha=0.5) + 
  geom_vline(aes(xintercept = R_upstr), size=2, color="#bb2649", alpha=0.5 ) +
  geom_vline(aes(xintercept = R_site), size=2, color="#bb2649", alpha=0.5, linetype=6 ) +
  geom_vline(aes(xintercept = quantile(res, 0.95)), size=1.5, color="black", alpha=0.5) +
  labs(x="Correlation", y="Count") + 
  theme_classic() + 
  theme(
    text = element_text(size=28)
  )





svg("figures/figure4/collisions_relCU.svg", height=5, width=4)

plot_grid(p.upstr, p.rand, ncol = 1, align = 'v')

dev.off()



# --------------------- Figure 4F ---------------------
total_young <- ccDF$young_site+ccDF$young_upstr
total_young <-total_young/sum(total_young)
total_old <- ccDF$old_site+ccDF$old_upstr
total_old <- total_old/sum(total_old)


scer_rd <- read.table("data/processed/scer_ribodensities.txt", header=T)
scer_codonusage = as.data.frame(read.table("data/processed/scer_codonusage.txt", header=T, sep='\t'))

scer_total <- data.frame(rd = ( scer_rd$RD_old / scer_rd$RD_young ),
             coll = (total_old / total_young )
              )


R_total = cor(scer_total$rd, scer_total$coll)

svg("figures/figure4/collisions_rRD_relCU.svg", height=5, width=5)

ggplot(scer_total, aes(x=rd, y=coll )) + 
  geom_smooth(method='glm', color="#186218", se=F) + 
  geom_point(color="#186218", size=2) + 
  labs(x="Relative RD", y="Relative Collision-CU") + 
  scale_color_manual(values=c("#186218", "#2cb42c")) + 
  geom_text(x=0.7, y=0.97, hjust=0, label=paste("R =", round(R_total, 2), sep=' '), col="#186218", size=5, stat="unique", show.legend=F) + 
  theme_classic() + 
  theme(
    text = element_text(size=28) ,
    axis.title = element_text(size=28),
    title = element_text(size=16),
    legend.title = element_blank(), 
    legend.position = 'top'
  )

dev.off()

# --------------------- Figure 4G ---------------------

optDF <- read.table("data/processed/optimality_cu.txt", header=T)


dwell.hypo <- as.matrix(read.table("data/SMoPT/dwell_hypo.txt"))
dwell.hypo1 <- as.matrix(read.table("data/SMoPT/dwell_hypo_CU1.txt"))


cor_self = rep(NA, 50)
cor_self1 = rep(NA, 50)
for (i in 3:50){
  cor_self[i]  <- cor(dwell$young, dwell.hypo[,i])
  cor_self1[i] <- cor(dwell$young, dwell.hypo1[,i])
}
plot(cor_self, t='l')
lines(cor_self1)

dwellh <- data.frame(pct=c(1:50), CU=cor_self, CU1=cor_self1)
dwellh.m  <- melt(dwellh, id='pct')



p.cor <- ggplot(dwellh.m, aes(x=pct, y=value)) + 
  geom_line(aes(color=variable), size=2) + 
  scale_color_manual(values=c("#bb2649", "#777777")) + 
  labs(x="Collision factor (%)", y="Correlation") + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    legend.position = 'none'
  )


boxdf = data.frame(c3=dwell.hypo[,3], 
                   c10=dwell.hypo[,10],
                   c20=dwell.hypo[,20],
                   c30=dwell.hypo[,30],
                   c40=dwell.hypo[,40],
                   c50=dwell.hypo[,50]
                   )
colnames(boxdf) <- c("3", "10", "20", "30", "40", "50")
boxdf.m <- melt(boxdf)


p.box <- ggplot(boxdf.m) + 
  geom_boxplot(aes(x=variable, y=value), fill="#bb2649") + 
  scale_y_continuous(limits=c(0, 2)) + 
  labs(x="Collision factor (%)", y="Sim. dwelltimes") + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )


svg("figures/figure4/sim_cor_box.svg", height=5, width=3.5)

plot_grid(p.cor, p.box, ncol = 1, align = 'v')

dev.off()


# --------------------- 

IC <- read.table("data/SMoPT/IC_hypo.txt", header=T)
IC1 <- read.table("data/SMoPT/IC_hypo_CU1.txt", header=T)

IC$relIC <- IC$IC / IC$IC[1]
IC$relIC1 <- IC1$IC/IC1$IC[1]

svg("figures/figure4/sim_entropy.svg", height=5, width=3.5)

ggplot(IC ) + 
  geom_line( aes(x=pct, y=relIC1), color="#777777", size=2) + 
  geom_line( aes(x=pct, y=relIC), color="#bb2649", size=2) + 
  #geom_point(size=2, shape=15) + 
  labs(x="Collision factor (%)", y="Rel. information content") + 
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
