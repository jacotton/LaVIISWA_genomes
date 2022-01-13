

setwd("~/Work/Schistosoma/JohnVianney")
library(ggplot2)
#multimodel inference
library(MuMIn)
library(knitr)
library(dplyr)
library(kableExtra)
data <- read.csv("fst_villages_3.csv")
#needed to make sure all models fit to only 'complete' data
options(na.action=na.fail)
#set -ve F_st to 1/10 of the smallest +ve value
data$pos_fst <- data$fst
data$pos_fst[data$pos_fst < 0] <- min(data$fst[data$fst > 0]) / 10
#transform fst to Nm Whitlock and McCauley, 1999, Heredity volume 82, pages 117–125. ITs not a good estimator of Nm (The population size scaled migration rate).

data$Nm <- 0.25 * (( 1 / data$pos_fst) - 1 )

#fit simple models of distance or 'island'

distance.model.fst.nolog <- lm(fst ~ linear_distance_km, data=data)
summary(distance.model.fst.nolog)
island.model.fst <- lm(fst ~ location, data=data)
summary(island.model.fst)
anova(island.model.fst)

#A typical gravity model models the number of migrants between populations in terms of the populations of two places and the distance between them:
#(e.g. http://ftp.iza.org/dp10329.pdf).
#for migration between $i$ and $j$ equal to $M_{ij}$ 
#$M_{ij} = G.\frac{P_i^a . P_j^b}{D_{ij}^Y}$
#where $P_i$ is the population of $i$, $D_{ij}$ is the distance between them, and $G$, $a$, $b$ and $Y$ are all estimated from the data.
#taking logs, this is estimated as a linear statistical model:


#or for a symmetric model - which I think is appropriate as $F_{st}$ is symmetric:
#$log(M_{ij}) = log(G) + a.log(P_i.P_j) - Y.log(D_{ij}) + \epsilon$

#which looks good - only 4 parameters.

gravity.model <- lm(log(Nm) ~ log(linear_distance_km) + log(pop1*pop2) ,data=data) 
summary(gravity.model)

full.model <- lm(log(Nm) ~ log(linear_distance_km)*location + log(pop1*pop2)*location ,data=data) 
summary(full.model)

#calculate AIC etc
d <- dredge(full.model)
kable(d,digits=2) %>% kable_styling(bootstrap_options = "condensed",font_size = 6)

#variable importance
good.models.99 <- get.models(d,cumsum(weight) <= .99)
kable(sw(good.models.99))

#So the interpretation here is that there is some weak evidence that distance and population size have some influence on $F_st$ the evidence that within- vs between- island is important is twice as strong.

## $F_{st}$ plot

#split character vectors
data$v1 <- sapply(strsplit(as.character(data$Villages),"_"),function(x) { x[1]})
data$v2 <- sapply(strsplit(as.character(data$Villages),"_"),function(x) { x[2]})
all_villages = unique(c(data$v1,data$v2))
#make sure these have the same levels
data$v1 <- factor(data$v1,levels=all_villages)
data$v2 <- factor(data$v2,levels=all_villages)
data$x <- factor(apply(cbind(data$v1,data$v2),1,function(x) all_villages[min(x)]),levels=all_villages)
data$y <- factor(apply(cbind(data$v1,data$v2),1,function(x) all_villages[max(x)]),levels=all_villages)


data$location2 <- c("within Damba island","between islands","between islands","within Damba island","within Koome island","between islands","within Damba island","between islands","within Damba island","between islands","between islands","between islands","within Koome island","between islands","within Damba island","between islands","between islands","within Koome island","between islands","between islands","between islands","within Koome island","between islands","between islands","between islands","between islands","between islands","between islands")

g2 <- ggplot(data,aes(x=x,y=y,color=fst,size=linear_distance_km,shape=location2) ) + geom_point() + scale_shape_manual(values=c(16,17,18),name="location") +
  theme_minimal() + scale_x_discrete(drop=FALSE)  +
  scale_y_discrete(drop=FALSE) + scale_size_continuous(range=c(2,8),name="linear distance (km)")  + 
  theme(axis.title.x = element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(angle=-90, hjust=0,vjust=0.5)) +
  scale_color_distiller(palette = "YlOrBr",name=bquote(F[st]))

cairo_pdf(file="S1figure_fstbyvillage_island.pdf",onefile=FALSE)
print(g2)
dev.off()
png(file="S1figure_fstbyvillage_island.png")
print(g2)
dev.off()


g <- ggplot(data,aes(y=fst,x=linear_distance_km,shape=location2,color=location2)) +
  geom_point(size=3) + scale_color_manual(name="",values=c("steelblue3","olivedrab3","green4")) + 
  theme_bw() + scale_shape_discrete(name="",) + 
  guides(shape = guide_legend(override.aes =list(color=c("steelblue3","olivedrab3","green4")))) + 
  ylab(bquote(F[st])) + xlab("linear distance (km)")


cairo_pdf(file="Figure4_distance_to_fst.pdf",onefile=FALSE)
print(g)
dev.off()
png(file="Figure4_disance_to_fst.png")
print(g)
dev.off()


