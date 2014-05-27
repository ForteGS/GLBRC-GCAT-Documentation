#test program; put test file in code dir and name file to match filename 

source("class.model.R")
source("class.well.R")
source("fit.model.R")
source("fitted.calculations.R")
source("GCAT.main.R")
source("misc.R")
source("normalize.and.transform.R")
source("plot.fit.R")
source("set.constants.R")
source("slope.analysis.R")
source("table2well.R")
source("table.output.R")
source("xlsx2well.R")
#time.format = "%Y-%m-%d %H:%M:%S"
time.input = "%m/%d/%Y %H:%M"
plate.nrow = 8
plate.ncol = 12

#Singleplate
#time.format = 1/3600
#filez = gcat.analysis.main("Revised_YPDAFEXglucoseTests_2-25-10.csv", out.dir = "/home/ndipiazza/Desktop/GCAT/R/stuff", single.plate = T, graphic.dir = "/home/ndipiazza/Desktop/GCAT/R/stuff", layout.file = NULL,  add.constant = 1 , blank.value = NULL, start.index = 2, growth.cutoff = .0455, points.to.remove = 0, remove.jumps = F, silent = F, verbose = F, return.fit = F, overview.jpgs = T)

time.input = time.format

#Multiplate
time.input = "%Y-%m-%d %H:%M:%S"
filez = gcat.analysis.main("august.csv", out.dir = "/home/ndipiazza/Desktop/GCAT/R/stuff", single.plate = F, graphic.dir = "/home/ndipiazza/Desktop/GCAT/R/stuff", layout.file = NULL,  add.constant = 1 , blank.value = NULL, start.index = 2, growth.cutoff = .0455, points.to.remove = 0, remove.jumps = F, silent = F, verbose = F, return.fit = F, overview.jpgs = T)

print(filez)
setwd("/home/ndipiazza/Desktop/GCAT/R")