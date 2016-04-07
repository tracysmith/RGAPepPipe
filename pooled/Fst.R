########my manhatten plot code
#####load data (on server)
dat <- read.table("C://Users/Mary/Documents/PepperellLab/Revisions/150507/SF_hiseq_Q20_mc6_manhatten.txt", header=T, sep='\t')
dat2 <- read.table("C://Users/Mary/Documents/PepperellLab/DiversityScales/Fst/150608_SF_Q20.fst.sort", header=F, sep='\t', na.strings="NA")

library(ggplot2)
library(reshape2)

#define colors based on fisher test
#bonferroni correction 0.05/14508
dat$col  <- with(dat, factor(ifelse(fet < 5.462638, 0, position)))


dat3 <- dat2[dat2$V4 > 5.462638,]

########theme set
theme_set(theme_bw(base_size = 10))

cbPalette <- c("grey", "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan",
               "red", "blue", "green", "violet", "cyan"
               )


p1 <- ggplot(dat, aes(x=position, y=fst)) +
  geom_point(aes(colour = col), size=3) +
  #facet_wrap(~key, nrow = 3) +
  #scale_colour_manual(values = c("grey", "red")) +
  scale_colour_manual(values = cbPalette) +
  ylab(expression(F[ST])) +
  xlab("Genomic Position (Mb)") +
  scale_x_continuous(breaks=c(0, 1e+06, 2e+06, 3e+06, 4e+06), labels=c("0", "1", "2", "3", "4")) +
  scale_y_continuous(limits=c(0,1)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        legend.title = element_blank()
  )
p1

p2 <- p1 + geom_text(data=subset(dat, fst > 0.25), aes(label=position),hjust=0,vjust=0) 
p2

vars <- colsplit(dat3$V10, ",", c("T1","T2","T3","T4","T5","T6","T7","T8","T9"))
vars$pos <- dat3$V5

p3 <- ggplot(vars) +
  geom_segment(aes(x="T1", y=T1, xend="T2", yend=T2)) +
  geom_segment(aes(x="T2", y=T2, xend="T3", yend=T3)) +
  geom_segment(aes(x="T3", y=T3, xend="T4", yend=T4))
p3


vars.m = melt(vars, id = "pos")

p4 <- ggplot(vars.m, aes(x=as.numeric(variable),y=(value))) + 
  geom_line(aes(group=pos, colour=as.factor(pos))) +
  scale_colour_manual(values = cbPalette) +
  ylab("Minor Allele Frequency") + 
  xlab("Timepoint") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        legend.title = element_blank())
p4  
  
  



