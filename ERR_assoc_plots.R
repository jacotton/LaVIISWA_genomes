library(ggplot2)
library(stringr)
library(wesanderson)
assoc_results <- read.table("plink.qassoc",header=TRUE)
#only keep big chromosomes:
fai <- read.delim("Schistosoma_mansoni_v7.fa.fai",header=FALSE)
big_chr <- fai[fai$V2 > 5000000,]
big_chr$V1 <- str_replace(big_chr$V1,"_[0-9]$","")

assoc_results_trim <- assoc_results[assoc_results$CHR %in% big_chr$V1,]

chr_sizes <- big_chr[,c(1,2,2,2)]
colnames(chr_sizes) <- c("CHR","LEN","START","TICK")
chr_sizes$START[1] <- 0
for (i in 2:nrow(chr_sizes)) { 
	partsum <- 0
	for(j in 1:(i-1)) { 
		partsum <- partsum +  chr_sizes$LEN[j]
	}
	chr_sizes$START[i] <- partsum
	chr_sizes$TICK[i] <- chr_sizes$START[i] + (chr_sizes$LEN[i] / 2 )
}
chr_sizes$SHORT_NAME <- str_replace(chr_sizes$CHR,"^SM_V7_","")
chr_sizes$TICK[1] <-chr_sizes$LEN[1] / 2
min(na.omit(assoc_results_trim$P))
min(p.adjust(na.omit(assoc_results_trim$P)))

assoc_results_trim2 <- merge(assoc_results_trim,chr_sizes,by="CHR")

#pretty plot
ggplot(assoc_results_trim2,aes(x=BP+START,y=-log(P,10),color=CHR)) + geom_point(size=2,alpha=0.5)  + theme_bw() + xlab("genomic position") + ylab("-log10(P)") + scale_color_manual(values=rep(wes_palette("Royal1",2),4)) 

g <- ggplot(assoc_results_trim2,aes(x=BP+START,y=-log(P,10),color=CHR)) + geom_point(size=1.5,alpha=0.5)  + theme_bw() + xlab("genomic position") + ylab("-log10(P)") + scale_color_manual(values=rep(wes_palette("Royal1",2),4)) + scale_x_continuous(breaks=chr_sizes$TICK,labels=chr_sizes$SHORT_NAME)
cairo_pdf(file="ERR_GWAS_Manhattan.pdf",onefile=FALSE)
g
dev.off()

cairo_png(file="ERR_GWAS_Manhattan.png")
g
dev.off()


#--- function of qqplots from https://gist.github.com/slowkow/9041570
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  obs <- -log10(sort(ps))
  obsn <- length(obs)
  df <- data.frame(
    observed = obs,
    expected = -log10(ppoints(obsn)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:obsn, shape2 = obsn:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:obsn, shape2 = obsn:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 1) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}


g2 <- gg_qqplot(assoc_results$P) + theme_bw()
cairo_pdf(file="ERR_GWAS_QQ.pdf",onefile=FALSE)
g2
dev.off()

png(file="ERR_GWAS_QQ.png")
g2 
dev.off()