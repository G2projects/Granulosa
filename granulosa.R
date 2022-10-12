# R

library(limma)
library(affy)
library(affyPLM)
library(Biobase)
library(edgeR)
library(gplots)
library(ggplot2)
library(reshape2)
library(preprocessCore)
library(OptimalCutpoints)
library(VennDiagram)
library(huge)

if (FALSE) {
library(car)
library(plyr)
library(scatterplot3d)
library(mclust)
library(lme4)

library(igraph)
library(KEGGgraph)
library(lavaan)
library(semTools)

library(org.Hs.eg.db)
library(biomaRt)
library(clusterProfiler)
library(DOSE)
library(GOSemSim)

library(randomForest)
}

source("/home/fernando/euNet/euNet_core/rutils.R")


##### ADD STRING LEGEND !!!


#ctrl <- read.delim("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/data/proteomics_v2/clean/ctrl.txt",
#                   stringsAsFactors = FALSE)
#fsh <- read.delim("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/data/proteomics_v2/clean/fsh.txt",
#                  stringsAsFactors = FALSE)
#lh <- read.delim("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/data/proteomics_v2/clean/gonasi.txt",
#                 stringsAsFactors = FALSE)
#fshlh <- read.delim("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/data/proteomics_v2/clean/fsh_gonasi.txt",
#                    stringsAsFactors = FALSE)

core <- read.delim("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/data/proteomics_v2/scores.txt",
                   stringsAsFactors = FALSE)

x <- read.delim("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/data/proteomics_v2/proteins.txt",
                stringsAsFactors = FALSE)



table(x$control > 0)

table(x$fsh > 0)

table(x$lh > 0)

table(x$fshlh > 0)


plot(log2(x$nPeptides_control[x$control > 0 & x$nPeptides_control > 0] + 1),
     log2(x$control[x$control > 0 & x$nPeptides_control > 0] + 0.001))

cor(log2(x$nPeptides_control[x$control > 0 & x$nPeptides_control > 0] + 1),
    log2(x$control[x$control > 0 & x$nPeptides_control > 0] + 0.001))


x$log2w_control <- log2(x$control + 1)
x$log2w_fsh <- log2(x$fsh + 1)
x$log2w_lh <- log2(x$lh + 1)
x$log2w_fshlh <- log2(x$fshlh + 1)

x$log2n_control <- log2(x$nPeptides_control + 1)
x$log2n_fsh <- log2(x$nPeptides_fsh + 1)
x$log2n_lh <- log2(x$nPeptides_lh + 1)
x$log2n_fshlh <- log2(x$nPeptides_fshlh + 1)


cor.test(x$log2n_control, x$log2w_control)
cor.test(x$log2n_control, x$log2w_control)$p.value

pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Protein_count_score_correlation_control.pdf", width = 15, height = 10)
ggplot(x, aes(x = log2n_control, y = log2w_control, color = "darkblue")) +
geom_point(size = 4) +
geom_smooth(method = lm, se = TRUE, fill = "grey85", fullrange = TRUE,
            lty = 1, lwd = 0.8) +
scale_colour_manual(name = , values = c("darkblue")) +
#scale_shape_manual(values = c(17, 16)) +
theme_bw() +
theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 18),
  axis.title = element_text(size = 16, face = "bold"),
  #legend.key.size = unit(1, "cm"),
  #legend.text = element_text(size = 14),
  #legend.title = element_text(size = 16),
  legend.position = "none") +
annotate("text", x = 5, y = 3.5, size = 7,
         label = "Pearson correlation: 0.976 (0.973, 0.978)\np < 0.001",
         color = "darkblue",
         parse = FALSE) +
scale_x_continuous(limits = c(0, 8)) +
scale_y_continuous(limits = c(0, 12)) +
labs(x = "log2 peptides number", y = "log2 protein score")
dev.off()


cor.test(x$log2n_fsh, x$log2w_fsh)
cor.test(x$log2n_fsh, x$log2w_fsh)$p.value

pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Protein_count_score_correlation_FSH.pdf", width = 15, height = 10)
ggplot(x, aes(x = log2n_fsh, y = log2w_fsh, color = "darkblue")) +
geom_point(size = 4) +
geom_smooth(method = lm, se = TRUE, fill = "grey85", fullrange = TRUE,
            lty = 1, lwd = 0.8) +
scale_colour_manual(name = , values = c("darkblue")) +
#scale_shape_manual(values = c(17, 16)) +
theme_bw() +
theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 18),
  axis.title = element_text(size = 16, face = "bold"),
  #legend.key.size = unit(1, "cm"),
  #legend.text = element_text(size = 14),
  #legend.title = element_text(size = 16),
  legend.position = "none") +
annotate("text", x = 5, y = 3.5, size = 7,
         label = "Pearson correlation: 0.980 (0.978, 0.982)\np < 0.001",
         color = "darkblue",
         parse = FALSE) +
scale_x_continuous(limits = c(0, 8)) +
scale_y_continuous(limits = c(0, 12)) +
labs(x = "log2 peptides number", y = "log2 protein score")
dev.off()


cor.test(x$log2n_lh, x$log2w_lh)
cor.test(x$log2n_lh, x$log2w_lh)$p.value

pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Protein_count_score_correlation_hCG.pdf", width = 15, height = 10)
ggplot(x, aes(x = log2n_lh, y = log2w_lh, color = "darkblue")) +
geom_point(size = 4) +
geom_smooth(method = lm, se = TRUE, fill = "grey85", fullrange = TRUE,
            lty = 1, lwd = 0.8) +
scale_colour_manual(name = , values = c("darkblue")) +
#scale_shape_manual(values = c(17, 16)) +
theme_bw() +
theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 18),
  axis.title = element_text(size = 16, face = "bold"),
  #legend.key.size = unit(1, "cm"),
  #legend.text = element_text(size = 14),
  #legend.title = element_text(size = 16),
  legend.position = "none") +
annotate("text", x = 5, y = 3.5, size = 7,
         label = "Pearson correlation: 0.974 (0.971, 0.976)\np < 0.001",
         color = "darkblue",
         parse = FALSE) +
scale_x_continuous(limits = c(0, 8)) +
scale_y_continuous(limits = c(0, 12)) +
labs(x = "log2 peptides number", y = "log2 protein score")
dev.off()


cor.test(x$log2n_fshlh, x$log2w_fshlh)
cor.test(x$log2n_fshlh, x$log2w_fshlh)$p.value

pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Protein_count_score_correlation_FSH_and_hCG.pdf", width = 15, height = 10)
ggplot(x, aes(x = log2n_fshlh, y = log2w_fshlh, color = "darkblue")) +
geom_point(size = 4) +
geom_smooth(method = lm, se = TRUE, fill = "grey85", fullrange = TRUE,
            lty = 1, lwd = 0.8) +
scale_colour_manual(name = , values = c("darkblue")) +
#scale_shape_manual(values = c(17, 16)) +
theme_bw() +
theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 18),
  axis.title = element_text(size = 16, face = "bold"),
  #legend.key.size = unit(1, "cm"),
  #legend.text = element_text(size = 14),
  #legend.title = element_text(size = 16),
  legend.position = "none") +
annotate("text", x = 5, y = 3.5, size = 7,
         label = "Pearson correlation: 0.974 (0.972, 0.977)\np < 0.001",
         color = "darkblue",
         parse = FALSE) +
scale_x_continuous(limits = c(0, 8)) +
scale_y_continuous(limits = c(0, 12)) +
labs(x = "log2 peptides number", y = "log2 protein score")
dev.off()



if (FALSE) {
wilcox.test(log2(x$fsh + 0.001), log2(x$control + 0.001), conf.int = TRUE)
wilcox.test(log2(x$lh + 0.001), log2(x$control + 0.001), conf.int = TRUE)
wilcox.test(log2(x$fshlh + 0.001), log2(x$control + 0.001), conf.int = TRUE)

wilcox.test(log2(core$fsh + 0.001), log2(core$control + 0.001), conf.int = TRUE)
wilcox.test(log2(core$lh + 0.001), log2(core$control + 0.001), conf.int = TRUE)
wilcox.test(log2(core$fshlh + 0.001), log2(core$control + 0.001), conf.int = TRUE)


M <- data.frame(score = c(core$fsh, core$control), y = c(rep(1, nrow(core)), rep(0, nrow(core))))
#M <- data.frame(score = c(core$lh, core$control), y = c(rep(1, nrow(core)), rep(0, nrow(core))))
#M <- data.frame(score = c(core$fshlh, core$control), y = c(rep(1, nrow(core)), rep(0, nrow(core))))
#M <- data.frame(score = c(core$fsh, core$lh), y = c(rep(1, nrow(core)), rep(0, nrow(core))))
#M <- data.frame(score = c(log2(x$fsh + 1), log2(x$control + 1)), y = c(rep(1, nrow(x)), rep(0, nrow(x))))

optimal.cutpoints(X = "score", status = "y",
                  tag.healthy = 0,
                  methods = "SpEqualSe",
                  data = M)

#    FSH vs. Ctrl:   
#     LH vs. Ctrl: 
# FSH+LH vs. Ctrl: 
#    FSH vs.   LH: 
}


test.proteins <- c("P02751", "Q12805_5", "Q71U36", "Q15149", "Q09666",
                   "P46940", "P35579", "Q13813", "P08238", "Q13509",
                   "Q9Y490", "O75369", "P0DMV9", "O00468_6")


p <- matrix(nrow = nrow(x), ncol = 4)
for (i in 1:nrow(x)) {
	p[i, ] <- unlist(lapply(x[i, c(5, 15, 25, 35)], function(w) ifelse(w == 0, w + 1E-300, w)))
}



score <- x[, c(3, 13, 23, 33)]
rownames(score) <- x$protein
#flag <- 4
#score <- x[x$flag >= flag, c(3, 13, 23, 33)]
#z <- huge.npn(score)
#rownames(score) <- x$protein[x$flag >= flag]
dim(score)

#plot(density(log2(score[, 1] + 0.001)), col = "blue", lwd = 3)
#lines(density(log2(score[, 2] + 0.001)), col = "red3", lwd = 3)
#lines(density(log2(score[, 3] + 0.001)), col = "darkorange", lwd = 3)
#lines(density(log2(score[, 4] + 0.001)), col = "green3", lwd = 3)
#abline(v = -4.4, lty = 2)

#t <- -4.4
#t <- 2^t - 0.001
#score <- score[score$control >= t & score$fsh >= t & score$lh >= t & score$fshlh >= t,]
#dim(score)

score <- as.matrix(score)

#wilcox.test(log2(score[, 2] + 0.001), log2(score[, 1] + 0.001), conf.int = TRUE)

Q <- normalize.quantiles(score)
#Q <- round(Q)
rownames(Q) <- rownames(score)
colnames(Q) <- colnames(score)
score <- Q
rm(Q)
head(score)

eset <- ExpressionSet(assayData = score)

#pdf("Exprs_MAplot.pdf")
MAplot(eset, pairs = TRUE, plot.method = "smoothScatter",
             main = "M vs. A plot (expression set)",
             cex = 0.9)
#dev.off()

f <- factor(c("FSH", "OTHER", "OTHER", "OTHER"))

#pdf("MDS_MA.pdf")
score.mds <- plotMDS(exprs(eset), col = c("blue", "red")[f])
abline(v = 0, col = "gray40", lty = 2)
abline(h = 0, col = "gray40", lty = 2)
#dev.off()

#designMatrix <- model.matrix(~0+f)
#colnames(designMatrix) <- levels(f)

#plot(density(log2(score[, 1] + 0.001)), col = "blue", lwd = 3)
#lines(density(log2(score[, 2] + 0.001)), col = "red3", lwd = 3)
#lines(density(log2(score[, 3] + 0.001)), col = "darkorange", lwd = 3)
#lines(density(log2(score[, 4] + 0.001)), col = "green3", lwd = 3)

lbl <- rownames(score.mds$distance.matrix.squared)

y <- c("Control", "FSH", "hCG", "FSH+hCG")

Z <- data.frame(Sample = lbl,
                Phenotype = y,
                Dimension1 = score.mds$x/10,
                Dimension2 = score.mds$y/10)


pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Protein_score_MDS.pdf", width = 15, height = 10)
ggplot(Z, aes(x = Dimension1, y = Dimension2, color = Phenotype, shape = Phenotype)) +
#geom_point(shape = 16, size = 5, colour = "black") +
geom_point(shape = 16, size = 5.5,
           colour = c("darkblue", "darkred", "chocolate", "darkgreen")) +
geom_point(size = 4) +
#geom_smooth(method = lm, se = TRUE, fill = "grey85", fullrange = TRUE,
#            lty = 1, lwd = 0.8) +
scale_colour_manual(values = c("dodgerblue1", "red3", "green3", "darkgoldenrod1")) +
scale_shape_manual(values = c(16, 16, 16, 16)) +
theme_bw() +
theme(panel.border = element_blank(),
  panel.grid.major = element_line(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(size = 18),
  axis.title = element_text(size = 16, face = "bold"),
  legend.key.size = unit(1, "cm"),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16)) +
annotate("text", x = 10, y = 0.6, size = 7, color = "darkblue",
         label = "Control",
         parse = FALSE) +
annotate("text", x = -3.6, y = 3.3, size = 7, color = "darkred",
         label = "FSH",
         parse = FALSE) +
annotate("text", x = -2.7, y = -0.6, size = 7, color = "chocolate",
         label = "hCG",
         parse = FALSE) +
annotate("text", x = -2.3, y = -2.6, size = 7, color = "darkgreen",
         label = "Combined",
         parse = FALSE) +
geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
#facet_wrap(facets =  vars(Phenotype)) +
#scale_x_continuous(limits = c(-3, 4.5)) +
#scale_y_continuous(limits = c(-4, 3)) +
#geom_text(label = lbl, nudge_y = 0.1, show.legend = FALSE) +
#labs(x = "MDS dimension 1", y = "MDS dimension 2")
labs(x = "MDS dimension 1 (83.0%)", y = "MDS dimension 2 (14.0%)")
dev.off()



w.fsh <- log2((score[, 2] + 1)/(score[, 1] + 1))
w.lh <- log2((score[, 3] + 1)/(score[, 1] + 1))
w.fshlh <- log2((score[, 4] + 1)/(score[, 1] + 1))

# Formula
# =SE(O(B3<2, ASS(SEGNO(C3)+SEGNO(E3)+SEGNO(G3))<2),0,SE(C3>=1,1,SE(C3<=-1,-1,0)))


pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Protein_score_shift.pdf", width = 20, height = 10)
plot(density(w.fshlh), col = "green3", lwd = 3,
     main = "",
     xlab = "log2 fold change",
     cex.axis = 1.6,
     cex.lab = 1.4)
lines(density(w.fsh), col = "red3", lwd = 3)
lines(density(w.lh), col = "darkorange", lwd = 3)
abline(v = -1, lty = 2, lwd = 3)
abline(v = 1, lty = 2, lwd = 3)
legend("topright", fill = c("red3", "darkorange", "green3", "black"),
                     bg = "white",
legend = c("log2(FSH/Control)",
           "log2(hCG/Control)",
           "log2(Combined/Control)",
           "|log2 fold change| = 1"),
lty = c(1, 1, 1, 2),
cex = 1.6)
dev.off()



x.fsh <- rownames(score)[abs(w.fsh) >= 1]
x.lh <- rownames(score)[abs(w.lh) >= 1]
x.fshlh <- rownames(score)[abs(w.fshlh) >= 1]


chosen.fsh <- x.fsh[x.fsh %in% test.proteins]
chosen.fsh

chosen.lh <- x.lh[x.lh %in% test.proteins]
chosen.lh

chosen.fshlh <- x.fshlh[x.fshlh %in% test.proteins]
chosen.fshlh



shift <- read.delim("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/shift_analysis/Protein_shift_clean.txt", stringsAsFactors = FALSE)

x.fsh <- shift$protein[shift$shifted_FSH != 0]
x.lh <- shift$protein[shift$shifted_hCG != 0]
x.fshlh <- shift$protein[shift$shifted_combo != 0]



n1 <- length(x.fsh)
n2 <- length(x.lh)
n3 <- length(x.fshlh)
n12 <- length(x.fsh[x.fsh %in% x.lh])
n13 <- length(x.fsh[x.fsh %in% x.fshlh])
n23 <- length(x.lh[x.lh %in% x.fshlh])
n123 <- length(x.fsh[x.fsh %in% x.lh & x.fsh %in% x.fshlh])

pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Shifted_proteins_venn.pdf", width = 20, height = 20)
grid.newpage()
draw.triple.venn(area1 = n1,        # FSH
                 area2 = n2,        # LH
                 area3 = n3,        # FSH+LH
                 n12 = n12,
                 n13 = n13,
                 n23 = n23,
                 n123 = n123,
          		 col = c("red3", "gold", "royalblue"),
          		 fill = c("red3", "gold", "deepskyblue2"),
          		 alpha = c(0.1, 0.1, 0.1),
          		 fontfamily = "sans",
          		 cex = 6,
				 cat.fontfamily = "sans",
          		 cat.col = c("firebrick", "darkorange", "darkblue"),
          		 cat.cex = 6,
				 rotation = 1,
                 category = c("FSH", "hCG", "FSH+hCG"),
                 cat.dist = c(0.05, 0.05, 0.05))
dev.off()



x.fsh.up <- rownames(score)[w.fsh >= 1]
x.lh.up <- rownames(score)[w.lh >= 1]
x.fshlh.up <- rownames(score)[w.fshlh >= 1]

n1 <- length(x.fsh.up)
n2 <- length(x.lh.up)
n3 <- length(x.fshlh.up)
n12 <- length(x.fsh.up[x.fsh.up %in% x.lh.up])
n13 <- length(x.fsh.up[x.fsh.up %in% x.fshlh.up])
n23 <- length(x.lh.up[x.lh.up %in% x.fshlh.up])
n123 <- length(x.fsh.up[x.fsh.up %in% x.lh.up & x.fsh.up %in% x.fshlh.up])

pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Shifted_proteins_venn_UP.pdf", width = 20, height = 20)
grid.newpage()
draw.triple.venn(area1 = n1,        # FSH
                 area2 = n2,        # LH
                 area3 = n3,        # FSH+LH
                 n12 = n12,
                 n13 = n13,
                 n23 = n23,
                 n123 = n123,
          		 col = c("red3", "gold", "royalblue"),
          		 fill = c("red3", "gold", "deepskyblue2"),
          		 alpha = c(0.1, 0.1, 0.1),
          		 fontfamily = "sans",
          		 cex = 6,
				 cat.fontfamily = "sans",
          		 cat.col = c("firebrick", "darkorange", "darkblue"),
          		 cat.cex = 6,
				 rotation = 1,
                 category = c("FSH", "hCG", "FSH+hCG"),
                 cat.dist = c(0.05, 0.05, 0.05))
dev.off()



x.fsh.down <- rownames(score)[w.fsh <= -1]
x.lh.down <- rownames(score)[w.lh <= -1]
x.fshlh.down <- rownames(score)[w.fshlh <= -1]

n1 <- length(x.fsh.down)
n2 <- length(x.lh.down)
n3 <- length(x.fshlh.down)
n12 <- length(x.fsh.down[x.fsh.down %in% x.lh.down])
n13 <- length(x.fsh.down[x.fsh.down %in% x.fshlh.down])
n23 <- length(x.lh.down[x.lh.down %in% x.fshlh.down])
n123 <- length(x.fsh.down[x.fsh.down %in% x.lh.down & x.fsh.down %in% x.fshlh.down])

pdf("/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/Shifted_proteins_venn_DOWN.pdf", width = 20, height = 20)
grid.newpage()
draw.triple.venn(area1 = n1,        # FSH
                 area2 = n2,        # LH
                 area3 = n3,        # FSH+LH
                 n12 = n12,
                 n13 = n13,
                 n23 = n23,
                 n123 = n123,
          		 col = c("red3", "gold", "royalblue"),
          		 fill = c("red3", "gold", "deepskyblue2"),
          		 alpha = c(0.1, 0.1, 0.1),
          		 fontfamily = "sans",
          		 cex = 6,
				 cat.fontfamily = "sans",
          		 cat.col = c("firebrick", "darkorange", "darkblue"),
          		 cat.cex = 6,
				 rotation = 1,
                 category = c("FSH", "hCG", "FSH+hCG"),
                 cat.dist = c(0.05, 0.05, 0.05))
dev.off()



















































counts <- as.matrix(x[, c(6, 16, 26, 36)])
#counts <- as.matrix(x[, c(8, 18, 28, 38)])
rownames(counts) <- x$protein

# CTRL, FSH, LH, FSH+LH
grp <- c("X0", "X1", "X1", "X1")
#grp <- c("X0", "X1", "X0", "X1")
#grp <- c("X0", "X1", "X0", "X0")

head(counts)
dim(counts)


# Normalization
d <- DGEList(counts=counts, group=grp)
d <- calcNormFactors(d)
nrow(d)
d2 <- d[rowSums(cpm(d) > 1) > 2,]
nrow(d2)

# Maximize the NB conditional common likelihood to give the
# estimate of the common dispersion across all tags
d2 <- estimateCommonDisp(d2)

# Trended dispersion values calculated using empirical Bayes
d2 <- estimateTrendedDisp(d2)

# Estimate tagwise dispersion values calculated
# using weighted conditional maximum likelihood.
d2 <- estimateTagwiseDisp(d2)


if (FALSE) {
# Mean-variance relationship plot
#pdf("MeanVar.pdf")
jpeg("MeanVariance_plot.jpg", width = 15, height = 10, units = 'in', res = 300)
plotMeanVar(d2, show.raw.vars=TRUE, show.tagwise.vars=TRUE, xlab = 'Mean expression level (log10 scale)', ylab = 'Pooled gene-level variance (log10 scale)')
dev.off()

# Dispersion versus mean plot
#pdf("DispVsMean.pdf")
jpeg("DispersionMean_plot.jpg", width = 15, height = 10, units = 'in', res = 300)
plotBCV(d2, log = "y")
dev.off()

# multidimensional scaling plot
#pdf("MDS_RNAseq.pdf")
jpeg("MDS_plot.jpg", width = 15, height = 10, units = 'in', res = 300)
plotMDS(d2, col=c("blue", "red2")[factor(grp)])
abline(v = 0, col = "gray40", lty = 2)
abline(h = 0, col = "gray40", lty = 2)
dev.off()
}


DEG <- exactTest(d2, pair = c("X0", "X1"))

rstt <- topTags(DEG, n = nrow(d2))
dim(rstt$table[rstt$table$FDR <= 0.05 & rstt$table$logFC > 0,])
dim(rstt$table[rstt$table$FDR <= 0.05 & rstt$table$logFC < 0,])

den <- rstt$table
write.delim(den, "X1vsX0_DEPs.txt", rownames = TRUE)

den <- rstt$table[rstt$table$FDR <= 0.05,]

#jpeg("MAplot_MBvsCtrlF.jpg", width = 15, height = 10, units = 'in', res = 300)
plotSmear(DEG, de.tags = rownames(den), cex = 0.8, col = "gold", deCol = "red3")
abline(h = -0.585, col = "blue", lty = 2, lwd = 1.5)
abline(h = 0.585, col = "blue", lty = 2, lwd = 1.5)
abline(h = 0, lty = 3)
#dev.off()







counts <- as.matrix(x[, c(6, 16, 26, 36)])
rownames(counts) <- x$protein

# CTRL, FSH, LH, FSH+LH
grp <- c("CTRL", "FSH", "LH", "FSH_LH")
#grp <- c("X0", "X1", "X1", "X1")

head(counts)
dim(counts)


# Normalization
d <- DGEList(counts=counts, group=grp)
d <- calcNormFactors(d)
nrow(d)
d2 <- d[rowSums(cpm(d) > 1) > 2,]
nrow(d2)

# Maximize the NB conditional common likelihood to give the
# estimate of the common dispersion across all tags
d2 <- estimateCommonDisp(d2)

# Trended dispersion values calculated using empirical Bayes
d2 <- estimateTrendedDisp(d2)

# Estimate tagwise dispersion values calculated
# using weighted conditional maximum likelihood.
d2 <- estimateTagwiseDisp(d2)


#plotBCV(d2, log = "y")

#plotMDS(d2, col=c("blue", "red2")[factor(grp)])
#abline(v = 0, col = "gray40", lty = 2)
#abline(h = 0, col = "gray40", lty = 2)


bcv <- 0.35
#DEG <- exactTest(d2, pair = c("CTRL", "FSH"), dispersion = bcv^2)
#DEG <- exactTest(d2, pair = c("CTRL", "LH"), dispersion = bcv^2)
DEG <- exactTest(d2, pair = c("CTRL", "FSH_LH"), dispersion = bcv^2)
rstt <- topTags(DEG, n = nrow(d2))

dim(rstt$table[rstt$table$FDR <= 0.05 & rstt$table$logFC > 0,])
dim(rstt$table[rstt$table$FDR > 0.05 & rstt$table$FDR <= 0.2 & rstt$table$logFC > 0,])

dim(rstt$table[rstt$table$FDR <= 0.05 & rstt$table$logFC < 0,])
dim(rstt$table[rstt$table$FDR > 0.05 & rstt$table$FDR <= 0.2 & rstt$table$logFC < 0,])

#write.delim(rstt$table, "FSHvsCTRL_DEPs.txt", rownames = TRUE)
#write.delim(rstt$table, "LHvsCTRL_DEPs.txt", rownames = TRUE)
write.delim(rstt$table, "FSH-LHvsCTRL_DEPs.txt", rownames = TRUE)



write.delim(x.fsh.up[x.fsh.up %in% x.lh.up & x.fsh.up %in% x.fshlh.up],
            "/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/UP_protein_list_COMMON.txt")

write.delim(x.fsh.up[!(x.fsh.up %in% x.lh.up) & !(x.fsh.up %in% x.fshlh.up)],
            "/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/UP_protein_list_FSH.txt")

write.delim(x.lh.up[!(x.lh.up %in% x.fsh.up) & !(x.lh.up %in% x.fshlh.up)],
            "/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/UP_protein_list_hCG.txt")

write.delim(x.fshlh.up[!(x.fshlh.up %in% x.fsh.up) & !(x.fshlh.up %in% x.lh.up)],
            "/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/UP_protein_list_COMBINED.txt")


write.delim(x.fsh.down[x.fsh.down %in% x.lh.down & x.fsh.down %in% x.fshlh.down],
            "/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/DOWN_protein_list_COMMON.txt")

write.delim(x.fsh.down[!(x.fsh.down %in% x.lh.down) & !(x.fsh.down %in% x.fshlh.down)],
            "/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/DOWN_protein_list_FSH.txt")

write.delim(x.lh.down[!(x.lh.down %in% x.fsh.down) & !(x.lh.down %in% x.fshlh.down)],
            "/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/DOWN_protein_list_hCG.txt")

write.delim(x.fshlh.down[!(x.fshlh.down %in% x.fsh.down) & !(x.fshlh.down %in% x.lh.down)],
            "/home/fernando/GemelliBioinfoUnit/Granulosa_Milardi/DOWN_protein_list_COMBINED.txt")










