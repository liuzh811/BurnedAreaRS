setwd(".\\data")

library(rgdal)
library(raster)
linrary(maptools)

########################################################################################################
############### begin of reading into validation fire  and bead-extracted fire patches #################
#read tested regions
testr1 = readOGR(dsn=".",layer="testr1")
testr3 = readOGR(dsn=".",layer="testr3")

#read validation fire patches 
fire.sp <- readOGR(".", "FirePatch")

#slect fire within the fire polygon
fire.sp = rbind(fire.sp[testr1,],fire.sp[testr3,]) 

#read bead-extracted fire patches for tested region 1
testr1_1_r4 = readOGR(dsn=".",layer="testr1_1_r4_ndvi_90_2")
testr1_4_r4 = readOGR(dsn=".",layer="testr1_4_r4_ndvi")

#read bead-extracted fire patches for tested region 2
testr3_1_r4 = readOGR(dsn=".",layer="testr3_1_r4_ndvi_90_2")
testr3_4_r4 = readOGR(dsn=".",layer="testr3_4_r4_ndvi_90_2")

#bead-extracted fire patches into one files
polygons <- slot(testr1_1_r4[which(testr1_1_r4$dist_time == 2000169),2:7], "polygons")
polygons <- c(polygons, slot(testr1_4_r4[which(testr1_4_r4$dist_time> 2010),2:7], "polygons"))
polygons <- c(polygons,slot(testr3_1_r4[which(testr3_1_r4$dist_time == 2000169),2:7], "polygons"))
polygons <- c(polygons,slot(testr3_4_r4[which(testr3_4_r4$dist_time > 2010),2:7], "polygons"))

dat = testr1_1_r4@data[which(testr1_1_r4$dist_time == 2000169),2:7]
dat = rbind(dat, testr1_4_r4@data[which(testr1_4_r4$dist_time> 2010),2:7])
dat = rbind(dat, testr3_1_r4@data[which(testr3_1_r4$dist_time== 2000169),2:7])
dat = rbind(dat, testr3_4_r4@data[which(testr3_4_r4$dist_time> 2010),2:7])
                        
# rename IDs of Polygons
 for (i in 1:length(polygons)) {
   slot(polygons[[i]], "ID") <- paste(i)
 }

spatialPolygons <- SpatialPolygons(polygons)
spdf <- SpatialPolygonsDataFrame(spatialPolygons, data=data.frame(Id = 1:length(polygons)))
spdf@data = cbind(spdf@data,dat)
projection(spdf)<- projection(fire.sp)

#plot the results
plot(spdf)
polygonsLabel(spdf, labels = spdf@data$dist_time, method = "centroid", cex = 1, col = "green", pch=19)
plot(testr1,add= T)
plot(testr3,add= T)

############### end of reading into validation fire  and bead-extracted fire patches #################
########################################################################################################

########################################################################################################
###############                 begin of calulating detection rate                     #################
# for each reference burned patches, calculate 
# (1) whether detect or not 
# (2) the correct area mapped by the agrithom,
# (3) the percentage of area over/underestimated, 
# (4) the date difference

# find out corresponding polygons
dat.idx = list()
for (i in 1:length(fire.sp)){
  tmp = c()
  for (j in 1:length(spdf)){
      if (length(gIntersection (fire.sp[i,],spdf[j,]))>0){tmp = c(tmp, j)} 
  }
  dat.idx[[i]] <- tmp
}

#for first fire, it should be cut into study area
fire.sp.1=gIntersection (fire.sp[1,],testr1)

#doing some plot to see whether bead-extracted fire and validation corresponding to each other
plot(fire.sp[5,])
plot(spdf[dat.idx[[5]][3],], add = T, border = "red") #shule be 3

plot(fire.sp[7,])
plot(spdf[dat.idx[[7]][1],], add = T, border = "red") #shule be 1

plot(fire.sp[16,])
plot(spdf[dat.idx[[16]][2],], add = T, border = "red") #shule be 1 and 3

# doing some spatial overlay
area.validation = sapply(slot(fire.sp.1, "polygons"), slot, "area")
area.detection = sapply(slot(spdf[dat.idx[[1]][1],], "polygons"), slot, "area")
area.intersection = sapply(slot(gIntersection (fire.sp.1, spdf[dat.idx[[1]][1],]), "polygons"), slot, "area")
correct.rate = area.intersection/area.validation

# calculating detection rate
area.dif = sapply(slot(gDifference(spdf[dat.idx[[1]][1],],gIntersection(fire.sp.1, spdf[dat.idx[[1]][1],])), "polygons"), slot, "area")
omission.rate = (area.validation - area.detection)/area.validation ##
coomission.rate = area.dif/area.validation ##
date.dif = spdf[dat.idx[[1]][1],]$dist_time - fire.sp[1,]$startDOY

#get data together
dat = c(1, area.validation, correct.rate, omission.rate, coomission.rate, date.dif)
for (i in 2:length(fire.sp)){
  if(is.null(dat.idx[[i]])){
      dat = rbind(dat, c(i, 0, 0, 0, 0))
      } else 
{
  if(i ==5){
    
    area.validation = sapply(slot(fire.sp[i,], "polygons"), slot, "area")
    area.detection = sum(sapply(slot(spdf[dat.idx[[i]][3],], "polygons"), slot, "area"))
    area.intersection = sapply(slot(gIntersection (fire.sp[i,], spdf[dat.idx[[i]][3],]), "polygons"), slot, "area")
    area.dif = sapply(slot(gDifference(spdf[dat.idx[[i]][3],],gIntersection(fire.sp[i,], spdf[dat.idx[[i]][3],])), "polygons"), slot, "area")
    correct.rate = area.intersection/area.validation
    omission.rate = (area.validation - area.detection)/area.validation ##
    coomission.rate = area.dif/area.validation ##
    date.dif = spdf[dat.idx[[i]][3],]$dist_time - fire.sp[i,]$startDOY
    dat = rbind(dat, c(i, area.validation, correct.rate, omission.rate, coomission.rate, date.dif))
    
    
  } else {
    
    area.validation = sapply(slot(fire.sp[i,], "polygons"), slot, "area")
    area.detection = sum(sapply(slot(spdf[dat.idx[[i]],], "polygons"), slot, "area"))
    area.intersection = sapply(slot(gIntersection (fire.sp[i,], spdf[dat.idx[[i]],]), "polygons"), slot, "area")
    area.dif = sapply(slot(gDifference(spdf[dat.idx[[i]],],gIntersection(fire.sp[i,], spdf[dat.idx[[i]],])), "polygons"), slot, "area")
    
    correct.rate = area.intersection/area.validation
    omission.rate = (area.validation - area.detection)/area.validation ##
    coomission.rate = area.dif/area.validation ##
    date.dif = spdf[dat.idx[[i]][1],]$dist_time - fire.sp[i,]$startDOY
    dat = rbind(dat, c(i, area.validation, correct.rate, omission.rate, coomission.rate, date.dif))
    
 }
}
}

colnames(dat) <- c("ID", "area","correct.rate","omission.rate","coomission.rate","date.dif")
dat = data.frame(dat)
dat$area = dat$area/10000

###############                   end of calulating detection rate                     #################
########################################################################################################

########################################################################################################
###############                   begin of plotting detection rate                     #################
par(mfrow=c(2,2),mar=c(2,0,1,0),oma=c(6,0.5,0.5,0))

par(mar=c(0, 5, 0, 0.5)) 
plot(dat$area, dat$correct.rate, log = "x", xaxt="n",
     #xlab = "Fire Size (ha, in log scale)", 
     ylab = "Correct Rate",
     cex = 2, cex.lab = 1.5, cex.axis = 1.5)
text(120,0.9, "a)", cex = 2)

lines(lowess(dat$area, dat$correct.rate), col="black", lwd = 1.5) # lowess line (x,y) 

par(mar=c(0, 5, 0, 0.5)) 
plot(dat$area, dat$omission.rate, log = "x", xaxt="n",
     #xlab = "Fire Size (ha, in log scale)", 
     ylab = "Omission Rate",
     cex = 2, cex.lab = 1.5, cex.axis = 1.5)
text(120,0.65, "b)", cex = 2)
lines(lowess(dat$area, dat$omission.rate), col="black", lwd = 1.5) # lowess line (x,y)  

par(mar=c(0, 5, 0, 0.5)) 
plot(dat$area, dat$coomission.rate, log = "x", 
     #xlab = "Fire Size (ha, in log scale)", 
     ylab = "Commission Rate",
     cex = 2, cex.lab = 1.5, cex.axis = 1.5)
text(120,0.062, "c)", cex = 2)
lines(lowess(dat$area, dat$coomission.rate), col="black", lwd = 1.5) # lowess line (x,y) 


par(mar=c(0, 5, 0, 0.5)) 
plot(dat$area, dat$date.dif, log = "x", ylim = c(-30, 30), 
     #xlab = "Fire Size (ha, in log scale)", 
     ylab = "Date difference",
     cex = 2, cex.lab = 1.5, cex.axis = 1.5)
text(120,28, "d)", cex = 2)
dat2 = dat[which(dat$date.dif<100),]
lines(lowess(dat2$area, dat2$date.dif),col="black", lwd = 1.5) # lowess line (x,y) 

mtext("Fire Size (ha, in log scale)", side = 1, cex = 1.5,outer=TRUE,padj = 2.3)

###############                   end of plotting   detection rate                     #################
########################################################################################################

