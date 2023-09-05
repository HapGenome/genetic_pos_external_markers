## This project is to interpolate 2010 soy genetic position to external markers
# The information of external markers is collected Nicholas Dietz

library(readxl)
library(tidyverse)

# setup working directory
setwd("/Users/taozuo/Library/CloudStorage/OneDrive-Bayer/Project/2023/misc/externalMarkers")

# physical positions of the external markers
# collected by Nicholas Dietz and processed by Paul Blischat
phyPos = read.table("/Users/taozuo/Downloads/soysnp50k_wm82.a1_41317.marker_info.08072023.txt",header=T)

# old map soy 2010 map; These markers are used for FpString
# I did not include all 70l markers because 1) these 27k markers are in better quality; 2) 27k markers is enough for interpolcate

oldMap = read.table("/Users/taozuo/Downloads/Soybean__MON_2010__Map.txt",header=T)
table(oldMap$genMapChr,oldMap$scaffoldName)

# sort by physical pos and remove the ones with physPos == 0
oldMap = oldMap %>% arrange(scaffoldName,physPos) %>% filter(physPos != 0)
ggplot(oldMap,aes(physPos/10^6,genMapPos,col=scaffoldName)) + geom_point() + 
  theme_bw() + facet_wrap(vars(scaffoldName),scale = "free_x")

# remove another problematic marker, off position
oldMap = oldMap %>% filter(markerName != "NGMAX008343420")
summary(oldMap$genMapPos[2:nrow(oldMap)] - oldMap$genMapPos[1:(nrow(oldMap) - 1)])
table(oldMap$genMapPos[2:nrow(oldMap)] - oldMap$genMapPos[1:(nrow(oldMap) - 1)] < 0)
oldMap$type = "bayerMks"

# replace the chr names of new markers
phyPos$Scaffold = str_replace_all(phyPos$Scaffold,"Gm","Gm_W82_CR")
table(phyPos$Scaffold %in% oldMap$scaffoldName)
phyPos$genMapChr = NA
phyPos$genMapPos = NA
colnames(phyPos)[2] = "scaffoldName"
phyPos$type = "extMks"

# combine two data together
phyPos = bind_rows(phyPos,oldMap[,c(2:6,13)])
phyPos = phyPos %>% arrange(scaffoldName,physPos)

outRes = data.frame()
for(chr in unique(phyPos$scaffoldName)){
  tmp = phyPos %>% filter(scaffoldName == chr)
  tmp$genMapPos_2010 = interpolate(tmp$physPos,tmp$genMapPos)[,3]
  tmp$genMapChr_2010 = names(sort(table(tmp$genMapChr),decreasing = T))[1]
  
  # the first and last markers are not interpolated. need to add genetic positions
  naIndex = which(is.na(tmp$genMapPos_2010))
  minIndexBmk = min(which(tmp$type == "bayerMks"))
  maxIndexBmk = max(which(tmp$type == "bayerMks"))
  
  for(idx in naIndex){
    if(idx < minIndexBmk){
      tmp[idx,"genMapPos_2010"] = tmp[minIndexBmk,"genMapPos_2010"] - 0.1
    }else if(idx > minIndexBmk){
      tmp[idx,"genMapPos_2010"] = tmp[maxIndexBmk,"genMapPos_2010"] + 0.1
    }else{
      print(paste(chr,"Error",idx))
    }
  }
  
  outRes = bind_rows(outRes,tmp)
}


# save the data
write.csv(outRes %>% filter(type == "extMks"),"soysnp50k_wm82.a1_41317.csv",row.names=F)

# interpolate
interpolate <- function(x, y, tol=1e-5) {
  #x is a numeric vector with no missing values (e.g. physical positions on a chr)
  #y is a numeric vector the same length as x (e.g. genetic positions on a chr)
  #
  #the result is, after sorting y in the same order as x, interpolating the missing values of y 
  #based on the nearest non-missing values of y and their corresponding values of x
  if (!is.vector(x) || !is.vector(y) || length(x) != length(y) || !is.numeric(x) || !is.numeric(y) || any(is.na(x))) 
    stop("x and y need to be numeric vectors of the same length without any missing values in x")
  xsorted <- sort(x)
  xorder <- order(x)
  ysorted <- y[xorder]
  yna <- is.na(ysorted)
  if (sum(yna) == 0) return(cbind(x=xsorted,y.in=ysorted,y.out=ysorted))
  ynai <- which(yna)
  nyna <- which(!yna)
  starts <- sapply(ynai, function(z) max(nyna[nyna < z]))
  stops  <- sapply(ynai, function(z) min(nyna[nyna > z]))
  diffsy <- ysorted[stops] - ysorted[starts]
  diffsx <- xsorted[stops] - xsorted[starts]
  ratios <- ifelse(abs(diffsx) > tol, diffsy / diffsx, 0)
  interp <- sapply(1:sum(yna), function(z) ysorted[starts[z]] + ratios[z]*(xsorted[ynai[z]] - xsorted[starts[z]]))
  yout <- ysorted
  yout[ynai] <- interp
  #if(y[min(nyna)] == 0 & min(nyna) > 1) {yout[1:(min(nyna)-1)] <- 0} # the first overlap is 0, then assigned any position before the first match to 0
  out <- cbind(x=xsorted,y.in=ysorted,y.out=yout)
  out
}

