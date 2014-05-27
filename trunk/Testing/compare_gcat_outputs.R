#  Compare 2 GCAT output files
#  Yury V Bukhman, 2013-09-04

WORK.DIR = "/mnt/file.glbrc.org/ybukhman/Projects/B05_Yeast_HTS/B05.16 GCAT 3.2 testing Oct 2013/testing round 1"
INPUT.1 = "GCAT 2 output/output_gcat.fit_2013-10-03_17.56.55.txt"
INPUT.2 = "GCAT 3 output/output_gcat.fit_2013-10-03_18.08.26.txt"

setwd(WORK.DIR)

#  Read in the data
in1 = read.table(INPUT.1,header=T,sep="\t")
in2 = read.table(INPUT.2,header=T,sep="\t")

#  Compare values
same = signif(in1$spec.growth,2) == signif(in2$spec.growth,2) & signif(in1$tot.growth,2) == signif(in2$tot.growth,2) &  signif(in1$lag.time,2) == signif(in2$lag.time,2)
same[is.na(in1$spec.growth) & is.na(in2$spec.growth)] = T

#  Write out the output
#out = file(OUTPUT,open="w")
#if (all(same)) { writeLines("No significant differences identified",con=out) }
#close(out)

#  Manual exploration
#  Ascertain that results are quite different
summary(same)
head(same)
head(in1)
head(in2)
summary(in2$spec.growth/in1$spec.growth) # specific growth lower in GCAT 3
summary(in2$tot.growth/in1$tot.growth)
summary(in2$lag.time/in1$lag.time)  #  lag time also lower in GCAT 3

#  Plot specific growth rate of GCAT 3 vs. GCAT 2
x = in1$spec.growth
y = in2$spec.growth
plot(x,y,pch=16,main="Specific growth rate", xlab="old GCAT", ylab="new GCAT")
abline(0,1)
abline(lm(y~x),lty=2)
legend("bottomright",c("y=x line","regression line"),lty=c(1,2))

#  Investigate outliers
subset(in1, is.finite(in1$spec.growth) & (in1$spec.growth < in2$spec.growth - 0.003))
subset(in2, is.finite(in1$spec.growth) & (in1$spec.growth < in2$spec.growth - 0.003))

#  Clean up
rm(list=ls())