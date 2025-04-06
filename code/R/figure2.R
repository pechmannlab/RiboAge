# RiboAge | Figure 2
# S Pechmann (2025)

library(ggplot2)
library(reshape2)
library(cowplot)
library(pROC)
library(uwot)

setwd("~/M2/RiboAge")




# --------------------- Figure 2A ---------------------
schematic



# --------------------- Figure 2B ---------------------
scer_rd <- read.table("data/processed/scer_ribodensities.txt", header=T)
scer_corr_young <- scan('data/processed/scer2_basic_young_CC.txt')
scer_DR_young <- as.matrix( read.table('data/processed/scer2_basic_young_DR.txt') )
scer_best_idx = sort(scer_corr_young, decreasing = T, index.return=T)$ix[1]
scer_best_CC = scer_corr_young[scer_best_idx]
scer_best_DR = scer_DR_young[,scer_best_idx]
scer_proof <- data.frame(TR=1/scer_rd$RD_young, PRED=scer_best_DR/max(scer_best_DR))

plot.scer <- ggplot(scer_proof, aes(x=TR, y=PRED)) + 
  geom_smooth(method=lm, alpha=0.1, size=0.8, linetype=2, color="#777777", fill="#228b22", se=T) + 
  geom_point(col="#228b22", size=3) + 
  geom_text(x=0.6, y=0.9, hjust=0, label=paste("R =", round(scer_best_CC, 2), sep=' '), size=7, stat="unique") + 
  scale_x_continuous(breaks=c(0.7, 1, 1.3, 1.6), labels=c(0.7, 1, 1.3, 1.6)) + 
  labs(x="TR = 1/RD", y="predicted TR", title="Scer") +
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.title = element_text(size=28),
    title = element_text(size=16)
  )


cele_rd <- read.table("data/processed/cele_ribodensities.txt", header=T)
cele_corr_young <- scan('data/processed/cele2_basic_young_CC.txt')
cele_DR_young <- as.matrix( read.table('data/processed/cele2_basic_young_DR.txt') )
cele_best_idx = sort(cele_corr_young, decreasing = T, index.return=T)$ix[1]
cele_best_CC = cele_corr_young[cele_best_idx]
cele_best_DR = cele_DR_young[,cele_best_idx]
cele_proof <- data.frame(TR=1/cele_rd$RD_young, PRED=cele_best_DR/max(cele_best_DR))

plot.cele <- ggplot(cele_proof, aes(x=TR, y=PRED)) + 
  geom_smooth(method=lm, alpha=0.1, size=0.8, linetype=2, color="#777777", fill="#9826bb", se=T) + 
  geom_point(col="#9826bb", size=3) + 
  geom_text(x=0.4, y=0.9, hjust=0, label=paste("R =", round(cele_best_CC, 2), sep=' '), size=7, stat="unique") + 
  labs(x="TR = 1/RD", y="predicted TR", title="Cele") +
  theme_classic() + 
  theme(
    text = element_text(size=28) ,
    axis.title = element_text(size=28),
    title = element_text(size=16)
  )


svg("figures/figure2/scatter_proof.svg", height=5, width=8)

plot_grid(plot.scer, plot.cele, ncol = 2, align = 'h')

dev.off()




# --------------------- Figure 2C ---------------------
scer_corr_young <- scan('data/processed/scer2_basic_young_CC.txt')
scer_corr_middle <- scan('data/processed/scer2_basic_middle_CC.txt')
scer_corr_old <- scan('data/processed/scer2_basic_old_CC.txt')

scer_cu_young <- scan('data/processed/scer3_cu_young_CC.txt')
scer_cu_middle <- scan('data/processed/scer3_cu_middle_CC.txt')
scer_cu_old <- scan('data/processed/scer3_cu_old_CC.txt')

scer_basic <- data.frame(age=c(0, 2, 4), 
                         opt=c(max(scer_corr_young), max(scer_corr_middle), max(scer_corr_old)),
                         t10=c(quantile(scer_corr_young, 0.9), quantile(scer_corr_middle, 0.9), quantile(scer_corr_old, 0.9)),
                         opt_cu=c(max(scer_cu_young), max(scer_cu_middle), max(scer_cu_old)),
                         t10_cu=c(quantile(scer_cu_young, 0.9), quantile(scer_cu_middle, 0.9), quantile(scer_cu_old, 0.9))
                        )



svg("figures/figure2/scer_basic_top10.svg", height=5, width=3.5)

ggplot(scer_basic) + 
  geom_hline(aes(yintercept = opt[1]), size=0.5, color="#222222", alpha=0.5, linetype=5) +
  geom_line(aes(x=age, y=opt), color="#228b22", size=2) + 
  geom_point(aes(x=age, y=opt), size=5, color="#228b22") +
  geom_ribbon((aes(x=age, ymin=t10, ymax=opt)), alpha=0.5, fill="#228b22") + 
  geom_line(aes(x=age, y=opt_cu), color="#777777", size=2) + 
  geom_point(aes(x=age, y=opt_cu), size=5, color="#777777") +
  geom_ribbon((aes(x=age, ymin=t10_cu, ymax=opt_cu)), alpha=0.5, fill="#777777") + 
  scale_x_continuous(breaks=c(0,2,4), labels=c(0,2,4)) +
  scale_y_continuous(limits=c(0.4,0.9)) + 
  labs(x="Age", y="Model fit (correlation)") + 
  theme_classic() + 
  theme(
    text = element_text(size=28)
  )

dev.off()




# --------------------- Figure 2D ---------------------
eng.young <- as.data.frame(read.table("data/processed/scer2_basic_young_E20.txt", header=T))
eng.young.m <- melt(eng.young)
eng.young.df <- data.frame(nt=eng.young.m[,1], energy=eng.young.m[,2], age=rep("0d", nrow(eng.young.m)))

eng.middle <- as.data.frame(read.table("data/processed/scer2_basic_middle_E20.txt", header=T))
eng.middle.m <- melt(eng.middle)
eng.middle.df <- data.frame(nt=eng.middle.m[,1], energy=eng.middle.m[,2], age=rep("2d", nrow(eng.middle.m)))

eng.old <- as.data.frame(read.table("data/processed/scer2_basic_old_E20.txt", header=T))
eng.old.m <- melt(eng.old)
eng.old.df <- data.frame(nt=eng.old.m[,1], energy=eng.old.m[,2], age=rep("4d", nrow(eng.old.m)))

eng.df <- rbind(eng.young.df, eng.middle.df, eng.old.df)




p.AA <- ggplot(eng.df[eng.df$nt=="AA",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="A") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

p.AC <- ggplot(eng.df[eng.df$nt=="AC",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),    
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p.AG <- ggplot(eng.df[eng.df$nt=="AG",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p.AT <- ggplot(eng.df[eng.df$nt=="AT",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p.CA <- ggplot(eng.df[eng.df$nt=="CA",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="C") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

p.CC <- ggplot(eng.df[eng.df$nt=="CC",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.CG <- ggplot(eng.df[eng.df$nt=="CG",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.CT <- ggplot(eng.df[eng.df$nt=="CT",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.GA <- ggplot(eng.df[eng.df$nt=="GA",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="G") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

p.GC <- ggplot(eng.df[eng.df$nt=="GC",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text = element_blank(), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p.GG <- ggplot(eng.df[eng.df$nt=="GG",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.GT <- ggplot(eng.df[eng.df$nt=="GT",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="")+
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.TA <- ggplot(eng.df[eng.df$nt=="TA",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="A", y="T") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )

p.TC <- ggplot(eng.df[eng.df$nt=="TC",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="C", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )

p.TG <- ggplot(eng.df[eng.df$nt=="TG",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="G", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )

p.TT <- ggplot(eng.df[eng.df$nt=="TT",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="T", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.y = element_blank(),    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )



svg("figures/figure2/scer_basic_energies.svg", height=7, width=8)
plot_grid(p.AA, p.AC, p.AG, p.AT, p.CA, p.CC, p.CG, p.CT, p.GA, p.GC, p.GG, p.GT, p.TA, p.TC, p.TG, p.TT, nrow=4, ncol = 4, align = 'hv')
dev.off()










# --------------------- Figure 2E ---------------------
reverse_rocdf <- function(roc){
  data <- data.frame(spec=roc$specificities, sens=roc$sensitivities)
  newspec <- c()
  newsens <- c()
  unique_spec <- unique(data$spec)
  for (i in 1:length(unique_spec)){
    sel <- data$spec == unique_spec[i]
    newspec <- c(newspec, rep(1-unique_spec[i], sum(sel)))
    newsens <- c(newsens, rev(data$sens[sel]) )
  }
  newdata <- data.frame(specificities=newspec, sensitivities=newsens)
}


parse_rocdf <- function(roc.data){
  roc_df <- data.frame(specificities=c(), sensitivities=c(), class=c(), weight=c(), feature=c() )
  for (i in (2:ncol(roc.data))){
    roc.1 <- roc(roc.data$label_0, roc.data[,i], direction="<")
    roc.1.rev <- reverse_rocdf(roc.1)
    roc.1.rev$class <- rep(paste("AUC",i-1, sep=''), length(roc.1.rev$specificities))
    roc.1.rev$weight <- rep(0.9, length(roc.1.rev$specificities))
    roc.1.rev$feature <- rep("all", length(roc.1.rev$specificities))
    roc_df <- rbind(roc_df, roc.1.rev)
    auc <- c(auc, roc.1$auc)
  }
  roc_df
}


compute_auc <- function(roc.data){
    auc <- c()
    for (i in (2:ncol(roc.data))){
    roc.1 <- roc(roc.data$label_0, roc.data[,i], direction="<")
    auc <- c(auc, roc.1$auc)
  }
  auc
}


roc.young <- read.table("data/processed/scer2_basic_young_ROC.txt", header=T)
roc.middle <- read.table("data/processed/scer2_basic_middle_ROC.txt", header=T)
roc.old <- read.table("data/processed/scer2_basic_old_ROC.txt", header=T)

roc_df_young <- parse_rocdf(roc.young)
roc_df_middle <- parse_rocdf(roc.middle)
roc_df_old <- parse_rocdf(roc.old)



p.roc.young <- ggplot(roc_df_young, aes(x=specificities, y=sensitivities, color=class)) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.5, linetype="twodash" ) +
  geom_line(aes(linetype=feature), size=0.5, show.legend=F) + 
  theme_classic() + 
  scale_colour_manual(values=rep("#2cb42c", length(unique(roc_df_young$class)))) +
  scale_x_continuous(breaks=c(0, 0.5, 1), labels=c(1, 0.5,  0)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1)) + 
  theme(
    axis.text = element_text(size=28),
    axis.title = element_text(size=28),
    title = element_text(size=16),
    text = element_text(size=28), 
    legend.position=c(0.85, 0.3), 
    legend.text = element_text(size=16),
    legend.title = element_blank()
  ) + 
  labs(x="Specificity", y="Sensitivity", title="0d")



p.roc.middle <- ggplot(roc_df_middle, aes(x=specificities, y=sensitivities, color=class)) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.5, linetype="twodash" ) +
  geom_line(aes(linetype=feature), size=0.5, show.legend=F) + 
  scale_colour_manual(values=rep("#228b22", length(unique(roc_df_middle$class)))) +
  theme_classic() + 
  scale_x_continuous(breaks=c(0, 0.5, 1), labels=c(1, 0.5,  0)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1)) + 
  theme(
    axis.text = element_text(size=28),
    axis.title = element_text(size=28),
    title = element_text(size=16),
    text = element_text(size=28), 
    legend.position=c(0.85, 0.3), 
    legend.text = element_text(size=16),
    legend.title = element_blank()
  ) + 
  labs(x="Specificity", y="Sensitivity", title="2d")


p.roc.old <-ggplot(roc_df_old, aes(x=specificities, y=sensitivities, color=class)) + 
  geom_abline(aes(intercept = 0, slope=1), size=0.5, linetype="twodash" ) +
  geom_line(aes(linetype=feature), size=0.5, show.legend=F) + 
  scale_colour_manual(values=rep("#186218", length(unique(roc_df_middle$class)))) +
  theme_classic() + 
  scale_x_continuous(breaks=c(0, 0.5, 1), labels=c(1, 0.5,  0)) + 
  scale_y_continuous(breaks=c(0, 0.5, 1), labels=c(0, 0.5, 1)) + 
  theme(
    axis.text = element_text(size=28),
    axis.title = element_text(size=28),
    title = element_text(size=16),
    text = element_text(size=28), 
    legend.position=c(0.85, 0.3), 
    legend.text = element_text(size=16),
    legend.title = element_blank()
  ) + 
  labs(x="Specificity", y="Sensitivity", title="4d")



auc_young <- compute_auc(roc.young)
auc_middle <- compute_auc(roc.middle)
auc_old <- compute_auc(roc.old)

aucDF <- data.frame(young=auc_young, middle=auc_middle, old=auc_old)
colnames(aucDF) <- c("0d", "2d", "4d")
aucDF.m <- melt(aucDF)



p.auc <- ggplot(aucDF.m) + 
  geom_boxplot(aes(x=variable, y=value, fill=variable)) +
  scale_fill_manual(values=c("#2cb42c", "#228b22", "#186218")) + 
  labs(x="", y="AUC") + 
  theme_classic() + 
  theme(
    text = element_text(size=28), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(), 
    legend.position = 'none'
  )



svg(file = "figures/figure2/scer_roc.svg", height = 7, width = 7)
plot_grid(p.roc.young, p.roc.middle, p.roc.old, p.auc, nrow=2, ncol = 2, align = 'hv')

dev.off()


# --------------------- Figure 2F ---------------------

eng.young <- as.data.frame(read.table("data/processed/scer2_basic_young_E20.txt", header=T))
eng.middle <- as.data.frame(read.table("data/processed/scer2_basic_middle_E20.txt", header=T))
eng.old <- as.data.frame(read.table("data/processed/sce2r_basic_old_E20.txt", header=T))
eng.mat <- rbind(eng.young, eng.middle, eng.old)
eng.mat <- as.matrix(eng.mat)



umap_data <- umap(eng.mat, metric='correlation', init="random", min_dist = 0.01, n_neighbors=10)
umap_data <- as.data.frame(umap_data)
colnames(umap_data) = c('umap1', 'umap2')


umap_data$subtype <- c(rep('0d', 20), rep('2d', 20), rep('4d', 20))
umap_data$subtype <- factor(umap_data$subtype, levels=c('0d', '2d', '4d') )
list_clusters <- unique(umap_data$subtype)



svg(file = "figures/figure2/scer_umap.svg", height = 5, width = 4.5)

ggplot(umap_data) + 
  geom_point(aes(x=umap1, y=umap2, color=subtype, shape=subtype), size=5) +
  labs(x="UMAP1", y="UMAP2") + 
  scale_color_manual(values=c("#2cb42c", "#228b22", "#186218")) +
  theme_classic()+ 
  theme(
    text = element_text(size=28), 
    legend.title = element_blank(), 
    legend.position = c(0.15, 0.8)
  )

dev.off()






# --------------------- Figure 2G ---------------------

cele_corr_young <- scan('data/processed/cele2_basic_young_CC.txt')
cele_corr_middle <- scan('data/processed/cele2_basic_middle_CC.txt')
cele_corr_old <- scan('data/processed/cele2_basic_old_CC.txt')

cele_cu_young <- scan('data/processed/cele3_cu_young_CC.txt')
cele_cu_middle <- scan('data/processed/cele3_cu_middle_CC.txt')
cele_cu_old <- scan('data/processed/cele3_cu_old_CC.txt')


cele_basic <- data.frame(age=c(0, 6, 12), 
                         opt=c(max(cele_corr_young), max(cele_corr_middle), max(cele_corr_old)),
                         t10=c(quantile(cele_corr_young, 0.9), quantile(cele_corr_middle, 0.9), quantile(cele_corr_old, 0.9)),
                         opt_cu=c(max(cele_cu_young), max(cele_cu_middle), max(cele_cu_old)),
                         t10_cu=c(quantile(cele_cu_young, 0.9), quantile(cele_cu_middle, 0.9), quantile(cele_cu_old, 0.9))
)



svg("figures/figure2/cele_basic_top10.svg", height=5, width=3.5)

ggplot(cele_basic) + 
  geom_hline(aes(yintercept = opt[1]), size=0.5, color="#222222", alpha=0.5, linetype=5) +
  geom_line(aes(x=age, y=opt), color="#9826bb", size=2) + 
  geom_point(aes(x=age, y=opt), size=5, color="#9826bb") +
  geom_ribbon((aes(x=age, ymin=t10, ymax=opt)), alpha=0.5, fill="#9826bb") + 
  geom_line(aes(x=age, y=opt_cu), color="#777777", size=2) + 
  geom_point(aes(x=age, y=opt_cu), size=5, color="#777777") +
  geom_ribbon((aes(x=age, ymin=t10_cu, ymax=opt_cu)), alpha=0.5, fill="#777777") + 
  scale_x_continuous(breaks=c(0,6,12), labels=c(1,6,12)) +
  scale_y_continuous(limits=c(0.4,0.9)) + 
  labs(x="Age", y="Model fit (correlation)") + 
  theme_classic() + 
  theme(
    text = element_text(size=28)
  )

dev.off()







# --------------------- Figure 2H ---------------------

eng.young <- as.data.frame(read.table("data/processed/cele2_basic_young_E20.txt", header=T))
eng.young.m <- melt(eng.young)
eng.young.df <- data.frame(nt=eng.young.m[,1], energy=eng.young.m[,2], age=rep("0d", nrow(eng.young.m)))

eng.middle <- as.data.frame(read.table("data/processed/cele2_basic_middle_E20.txt", header=T))
eng.middle.m <- melt(eng.middle)
eng.middle.df <- data.frame(nt=eng.middle.m[,1], energy=eng.middle.m[,2], age=rep("2d", nrow(eng.middle.m)))

eng.old <- as.data.frame(read.table("data/processed/cele2_basic_old_E20.txt", header=T))
eng.old.m <- melt(eng.old)
eng.old.df <- data.frame(nt=eng.old.m[,1], energy=eng.old.m[,2], age=rep("4d", nrow(eng.old.m)))

eng.df <- rbind(eng.young.df, eng.middle.df, eng.old.df)



p.AA <- ggplot(eng.df[eng.df$nt=="AA",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="A") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

p.AC <- ggplot(eng.df[eng.df$nt=="AC",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),    
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p.AG <- ggplot(eng.df[eng.df$nt=="AG",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p.AT <- ggplot(eng.df[eng.df$nt=="AT",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p.CA <- ggplot(eng.df[eng.df$nt=="CA",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="C") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

p.CC <- ggplot(eng.df[eng.df$nt=="CC",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.CG <- ggplot(eng.df[eng.df$nt=="CG",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.CT <- ggplot(eng.df[eng.df$nt=="CT",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="") + 
  theme_classic() +
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.GA <- ggplot(eng.df[eng.df$nt=="GA",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="G") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = 'none'
  )

p.GC <- ggplot(eng.df[eng.df$nt=="GC",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text = element_blank(), 
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )

p.GG <- ggplot(eng.df[eng.df$nt=="GG",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.GT <- ggplot(eng.df[eng.df$nt=="GT",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="", y="")+
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
  )


p.TA <- ggplot(eng.df[eng.df$nt=="TA",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="A", y="T") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )

p.TC <- ggplot(eng.df[eng.df$nt=="TC",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="C", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )

p.TG <- ggplot(eng.df[eng.df$nt=="TG",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="G", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )

p.TT <- ggplot(eng.df[eng.df$nt=="TT",]) + 
  geom_hline(aes(yintercept = 0), size=0.8, color="#222222") +
  geom_boxplot(aes(x=age, y=energy, fill=age)) + 
  scale_y_continuous(breaks=c(-6, -3, 0, 3), labels=c(-6, -3, 0, 3), limits=c(-7.5, 4)) + 
  scale_fill_manual(values=c("#b33dd7", "#9826bb", "#761d91")) + 
  labs(x="T", y="") + 
  theme_classic() + 
  theme(
    text = element_text(size=28),
    axis.text.y = element_blank(),    
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )



svg("figures/figure2/cele_basic_energies.svg", height=7, width=8)
plot_grid(p.AA, p.AC, p.AG, p.AT, p.CA, p.CC, p.CG, p.CT, p.GA, p.GC, p.GG, p.GT, p.TA, p.TC, p.TG, p.TT, nrow=4, ncol = 4, align = 'hv')
dev.off()






# --------------------- Figure 2J ---------------------

cele.eng.young <- as.data.frame(read.table("data/processed/cele2_basic_young_E20.txt", header=T))
cele.eng.middle <- as.data.frame(read.table("data/processed/cele2_basic_middle_E20.txt", header=T))
cele.eng.old <- as.data.frame(read.table("data/processed/cele2_basic_old_E20.txt", header=T))
cele.eng.mat <- rbind(cele.eng.young, cele.eng.middle, cele.eng.old)
cele.eng.mat <- as.matrix(cele.eng.mat)


umap_data <- umap(cele.eng.mat, metric='correlation', init="random", min_dist = 0.01, n_neighbors=10)
umap_data <- as.data.frame(umap_data)
colnames(umap_data) = c('umap1', 'umap2')


umap_data$subtype <- c(rep('0d', 20), rep('6d', 20), rep('12d', 20))
umap_data$subtype <- factor(umap_data$subtype, levels=c('0d', '6d', '12d') )
list_clusters <- unique(umap_data$subtype)





svg(file = "figures/figure2/cele_umap.svg", height = 5, width = 4.5)

ggplot(umap_data) + 
  geom_point(aes(x=umap1, y=umap2, color=subtype, shape=subtype), size=5) +
  labs(x="UMAP1", y="UMAP2") + 
  scale_color_manual(values=c("#b33dd7", "#9826bb", "#761d91")) +
  theme_classic()+ 
  theme(
    text = element_text(size=28), 
    legend.title = element_blank(), 
    legend.position = c(0.18, 0.18)
  )

dev.off()





# --------------------- Figure 2K ---------------------




scer.eng.young <- as.data.frame(read.table("data/processed/scer2_basic_young_E20.txt", header=T))
scer.eng.middle <- as.data.frame(read.table("data/processed/scer2_basic_middle_E20.txt", header=T))
scer.eng.old <- as.data.frame(read.table("data/processed/scer2_basic_old_E20.txt", header=T))


cele.eng.young <- as.data.frame(read.table("data/processed/cele2_basic_young_E20.txt", header=T))
cele.eng.middle <- as.data.frame(read.table("data/processed/cele2_basic_middle_E20.txt", header=T))
cele.eng.old <- as.data.frame(read.table("data/processed/cele2_basic_old_E20.txt", header=T))


data.young <- data.frame(scer=colMeans(scer.eng.young), cele=colMeans(cele.eng.young))
data.old<- data.frame(scer=colMeans(scer.eng.old), cele=colMeans(cele.eng.old))



p.young <- ggplot(data.young, aes(x=scer, y=cele)) + 
  geom_smooth(method=lm, alpha=0.1, size=0.8, linetype=2, color="#777777", fill="#6667ab", se=T) + 
  geom_point(color="#739bd0", size=3) + 
  labs(x="Scer (0d)", y="Cele (1d)") +
  geom_text(x=0, y=-2.3, hjust=0, label=paste("R =", round(cor(data.young)[1,2], 2), sep=' '), size=5, stat="unique") + 
  theme_classic() + 
  theme(
    text = element_text(size=24)
  )


p.old <- ggplot(data.old, aes(x=scer, y=cele)) + 
  geom_smooth(method=lm, alpha=0.1, size=0.8, linetype=2, color="#777777", fill="#6667ab", se=T) + 
  geom_point(color="#0f4c81", size=3) + 
  labs(x="Scer (4d)", y="Cele (12d)") +
  geom_text(x=0.6, y=-2.5, hjust=0, label=paste("R =", round(cor(data.old)[1,2], 2), sep=' '), size=5, stat="unique") + 
  theme_classic() + 
  theme(
    text = element_text(size=24)
  )




svg(file = "figures/figure2/correlations.svg", height = 5, width = 4)

plot_grid(p.young, p.old, ncol = 1, align = 'v')

dev.off()

