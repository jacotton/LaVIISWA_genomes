#---
#load some requirements
library(tidyr)
library(dplyr)
library(MCMCglmm)
library(colorspace)


#---
#Read and reshape some data
#----

#---
#including treatment arm
#
#intensive  : Busi, Kitosi, Katooke, Kisu
#standard: Kachanga, Kakeeka, Lugumba, Zingoola
village_2_arm <- cbind(c("Busi","Kitosi", "Katooke","Kisu","Kachanga","Kakeeka","Lugumba", "Zingoola"),c(rep("intensive",4),rep("standard",4)))
colnames(village_2_arm) <- c("village","arm")

data <- read.delim("Eggs_88 - Eggs_88.tsv")
#get the dataframe in a shape for analysis

new_data <- data %>% pivot_longer(c("kasepg_time1sample1","kasepg_time1sample2","kasepg_time1sample3","kasepg_time2sample1","kasepg_time2sample2","kasepg_time2sample3"),names_to = "tmnt",values_to="epg")
#not obvious what 'timepoint' means here
new_data$treatment[grep("time1",new_data$tmnt,fixed=TRUE)] <- 0
new_data$treatment[grep("time2",new_data$tmnt,fixed=TRUE)] <- 1
new_data$treatment <- as.numeric(new_data$treatment )
new_data$household <- paste0(new_data$vno,":",new_data$hhno)
data_for_model <- na.omit(new_data[c(1,8,9,22,23)])


colnames(data_for_model) <- c("individual","village","treatment","count","household")
data_for_model <- merge(data_for_model,village_2_arm,by.x="village")

#-----
#run the model
#-----
#don't have a random effect for site_name:treatment, hence marginalisation code below fails for this one.

#without "household" random effect
ranef_comb <-  formula(~us(1+treatment):individual + us(1):village + us(1):arm)
overdisp = formula(~idh(1):units)
prior2a <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = diag(c(1,1)), nu=0.002), 
                        G2=list(V = 1, nu=0.002),G3=list(V = 1, nu=0.002)))
full_form <- formula(count ~ treatment)
nitt <- 1000000; burnin <- 500000
full_comb <- MCMCglmm(fixed = full_form, random = ranef_comb,
                       rcov = overdisp, family = "poisson", data = data_for_model, prior = prior2a,
                       thin = 100, verbose = TRUE, nitt = nitt, burnin = burnin, pr=TRUE, pl=TRUE)

source("margeff.LaVIISWA.R")
source("ERR-dist-func.LaVIISWA.R")

#--
#get marginal posterior estimates for individuals, villages and treatment arms
out_full <- margeff(modlist = list(full_comb), data = data_for_model, nm = "treatment", marg = "individual")
out_village <- margeff(modlist = list(full_comb), data = data_for_model, nm = "treatment", marg = "village")
out_arm  <- margeff(modlist = list(full_comb), data = data_for_model, nm = "treatment", marg = "arm")
#---

my_err <- out_full$mn
my_err[my_err < 0] <- 0
write.table(cbind(out_full$id,out_full$id,my_err), file = "JV.ERR.phenotypes", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE)
#--------
#plotting full posteriors_of_stuff

arm_ERRdist <- ERRdist(modlist = list(full_comb), data = data_for_model, nm = "treatment", marg = 
"arm")
arm2 <- t(arm_ERRdist[,-ncol(arm_ERRdist)])
colnames(arm2) <- unique(data_for_model$arm)
arm3 <- pivot_longer(data.frame(arm2),cols=c("intensive","standard"),names_to="treatment",values_to="ERR")

vill_ERRdist <- ERRdist(modlist = list(full_comb), data = data_for_model, nm = "treatment", marg="village")
vill2 <- t(vill_ERRdist[,-ncol(vill_ERRdist)])
colnames(vill2) <- vill_ERRdist[,ncol(vill_ERRdist)]
vill3 <- pivot_longer(data.frame(vill2),cols=everything(),names_to="village",values_to="ERR")
vill3 <- merge(vill3,village_2_arm)
colnames(vill3)[3] <- "treatment"

library(wesanderson)
library(patchwork)
pal <-  wes_palette("Darjeeling1", n = 2)
g1 <- ggplot(arm3,aes(x=ERR,y=treatment,fill=treatment)) + geom_violin(linetype=0,alpha=0.75) + theme_bw() + xlim(0,1) + scale_fill_manual(values=pal) + geom_vline(xintercept=)
g2 <- ggplot(vill3,aes(x=ERR,y=village,fill=treatment)) + geom_violin(linetype=0,alpha=0.75) + theme_bw() + xlim(0,1) + scale_fill_manual(values=pal,guide=FALSE)
cairo_pdf(file="ERR_dist_by_village_and_arm.pdf",onefile=FALSE)
g1 / g2
dev.off()

library(scales)

## graphical parameters
smallpointsize <- 5
pointsize<-10
bigpointsize<-pointsize*1.5
fatlinesize<-1
thinlinesize <- 0.75
verythinlinesize <- 2
basefontsize<-20
smallfontsize <- 18
meanshape<-0
pointshape<-25
vjust<-1
hjust=-0.1
limits = c(0.01, 0.9999)
shade <- "grey60" 

max<-0.9999;min<-0.01
#just plotting posteriors for ERR for each individual

dfplot <- out_full
dfplot <- within(dfplot, {
  med[med>=max] <- max
  upr[upr>=max] <- max
  lwr[lwr>=max] <- max
  upr[upr<=min] <- min
  lwr[lwr<=min] <- min
  upr[upr<=min] <- min
  lwr[lwr<=min] <- min
})
#order them by ascending posterior mean
dfplot <- dfplot[order(dfplot$mn),]
# setting x-axis position to be 1:num of individuals
dfplot$x <- seq(1, nrow(dfplot))

extra <- unique(data_for_model[,c(2,1,6)])
colnames(extra)[1] <- "id"
dfplot_extra <- inner_join(dfplot,extra,by="id")
facet_labeller <- sub("  village","",dfplot_extra$village)
names(facet_labeller) <- dfplot_extra$site_name
dfplot_extra$x <- as.factor(dfplot_extra$x)

p <- ggplot(data=dfplot_extra, aes(y=mn, x=x), shape=46, col="black")+
  geom_linerange(aes(ymax=upr, ymin=lwr, x=x), size = 0.75, col="grey60")+
facet_grid(facets=.~village ,space="free_x",scales="free_x",drop=TRUE,labeller = labeller(site_name=facet_labeller))+  geom_point(size=1, col="black")+
  scale_y_continuous(labels=c(1, 10, 50 , 90, 99, 100),
                     breaks=c(min,0.1,0.5,0.9,0.99, max),
                      trans=logit_trans()) +
  geom_hline(aes(yintercept=.90),size = 0.75,linetype="dashed", col = "black")+
  theme_set(theme_bw(22))+
  xlab("Individual")+
  ylab("Estimated egg reduction rate (%)")+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=16,face="bold"),
        panel.background=element_blank(),strip.text.x = element_text( size = 7),axis.title.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())

##nasty grob manipulation code from https://stackoverflow.com/questions/41631806/change-facet-label-text-and-background-colour
g <- ggplot_gtable(ggplot_build(p))

striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) ) 
m <- merge(levels(as.factor(dfplot_extra$village)),village_2_arm,by.x="x",by.y="village")
#must be an easier way to add alpha to an RGB hex color!
r1 <- hex2RGB(pal[1])
r1.rgb <- rgb(r1@coords[1],r1@coords[2],r1@coords[3],alpha=0.75)
r2 <- hex2RGB(pal[2])
r2.rgb <- rgb(r2@coords[1],r2@coords[2],r2@coords[3],alpha=0.75)
m$col <- r2.rgb
m$col[m$arm == "intensive"] <- r1.rgb

fills <- m$col
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
cairo_pdf(file="individualERR_per_village.pdf",onefile=FALSE)
grid::grid.draw(g)
dev.off()
