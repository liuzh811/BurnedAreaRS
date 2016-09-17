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

########################################################################################################
### 4/21/2016 validate BEAD-derived polygon with MODIS fire products ###################################
### PART 1: MODIS burned area 
#downland modis burned area data to validate the results
#downland product == MCD45A1, version = 5.1 #MCD45A1.051

# load libraries
library(MODIS)
library(rgdal)
library(raster)

wkdir <- "D:/users/Zhihua/Landsat/MODIS_validation" # in rosa
setwd(wkdir)

#set spatial extent
dxal <- readOGR(dsn = "D:/users/Zhihua/Landsat/XinganImages/boundry", layer = "dxal_bj_proj_polygon")
proj.utm = projection(dxal)
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

####################### downlaod MODIS burned area data ##########################
#set MODIS path, see MODIS package for details
MODISoptions(localArcPath="D:/users/Zhihua/Landsat/MODIS_validation",
             outDirPath="D:/users/Zhihua/Landsat/MODIS_validation",
             gdalPath='c:/OSGeo4W64/bin')

#set date period for data
dates <- as.POSIXct( as.Date(c("1/1/2000","1/11/2015"),format = "%d/%m/%Y") )
dates2 <- transDate(dates[1],dates[2]) # Transform input dates from before

# getProduct() # list available MODIS products
#download MOD13Q1 data
runGdal(product="MCD45A1",  #VI/combined/Tile/500m/monthly
        begin=dates2$beginDOY, #start date
        end = dates2$endDOY,#end date
        tileH = 25,tileV = 3, #tile
        SDSstring = "11", #extract the first 3 layers
        extent = dxal, #crop to extent
        outProj=proj.utm, #reproject to UTM
        job = "MCD45A1") #download folder

# processing 
wkdir <- "F:/Rosa/Landsat/MODIS_validation"
setwd(wkdir)

#set spatial extent
dxal <- readOGR(dsn = "F:/Rosa/Landsat/XinganImages/boundry", layer = "dxal_bj_proj_polygon")
proj.utm = projection(dxal)
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

ba_date = preStack(path=".\\MCD45A1", pattern="*.burndate.tif$")
ba_quality = preStack(path=".\\MCD45A1", pattern="*.ba_qa.tif$")

#1. produce yearly burned area data
YearDOY = substr(ba_date, 20, 26)

ba.list = list()

for (yr in 2000:2015){
Year2010 = which(substr(YearDOY, 1, 4)==as.character(yr))

ba_date2010 = stack(ba_date[Year2010])
ba_quality2010 = stack(ba_quality[Year2010])

# remove bad data based on data quality
for (i in 1:nlayers(ba_date2010)) {
  
  ba_date2010[[i]][ba_quality2010[[i]] != 1] = 0
  
}

# to assess whether there are multiple burn or not within one year
Non.zero.length = function(x){length(which(x > 0))}
date2010_freq = calc(ba_date2010, Non.zero.length)
#freq(ba_date2010_1)

# pixels with multple date usually have the same date, but may contain some have different date, 
# use the first date to assign the date information

#get the location of the non-zero values
date2010_freq2 = date2010_freq > 0
date2010_freq2[date2010_freq2==0] = NA

Ex.pts.all.nonNA = function(x){
  proj.geo = projection(x)
  pts = rasterToPoints(x) #get raster coordinate, from left to right, top to bottom
  pts = data.frame(pts)
  pts <- SpatialPoints(coords = cbind(pts$x,pts$y),proj4string = CRS(proj.geo))
  pts.sp = SpatialPoints(coords = pts, proj4string = CRS(proj.geo))
  return(pts.sp)
}

#
date2010.sp = Ex.pts.all.nonNA(date2010_freq2)

#extract burn date values
date2010.date.df = extract(ba_date2010, date2010.sp)
#get the burn date for each pixels
date2010.date.df2 = apply(date2010.date.df, 1, function(x){x[which(x>1)[1]]})

#convert from point to raster
date2010.sp2 = data.frame(rasterToPoints(date2010_freq2)[,c(1,2)], date = date2010.date.df2)
coordinates(date2010.sp2) <- ~x+y
proj4string(date2010.sp2) = projection(date2010_freq)
date2010.sp2 = rasterize(date2010.sp2,date2010_freq2, field = "date")

date2010.sp2[date2010.sp2 > 273] = NA

ba.list[[yr - 1999]] = date2010.sp2

print(paste("Finish Extracting Year ", yr, " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

ba.list = stack(ba.list)
ba.list = crop(ba.list, dxal[3,])

#plot
names(ba.list) <- paste("Year", 2000:2015, sep = "")
ba.list2 = ba.list

#produce one raster, and get the burn date from the rasterStack
ba.grd = ba.list2[[1]]; ba.grd[] = 0
for (i in 1:nlayers(ba.list2)){
ba.grd[!is.na(ba.list2[[i]])] = ba.list2[[i]][!is.na(ba.list2[[i]])] + (i+1999)*1000
}

writeRaster(ba.grd,"F:/Rosa/Landsat/MODIS_validation/results/ba00-15.tif",format="GTiff", overwrite=TRUE) #write a raster stack files

#plot to see the figure
plot(date2010.sp2)
plot(dxal, add = T)

# read into manually-digitized fire polygon
fire.sp = readOGR(dsn="F:/Elias/XinganTM/GIS_data","FirePatch") 

# test PLOTING
#for test area 1
#read into test areas
require("rgeos")
r1 = readOGR(dsn=".\\BEAD_ploys",layer="testr1")

fire.sp1 = crop(fire.sp, r1)
ba.grd1 = crop(ba.grd, r1)

r1_firepoly1 = readOGR(dsn=".\\BEAD_ploys",layer="testr1_4_r4_ndvi")
r1_firepoly2 = readOGR(dsn=".\\BEAD_ploys",layer="testr1_4_r4_ndvi_90_2")

#for test area 3
#read into test areas
r3 = readOGR(dsn=".\\BEAD_ploys",layer="testr3")
fire.sp3 = crop(fire.sp, r3)
ba.grd3 = crop(ba.grd, r3)

r3_firepoly1 = readOGR(dsn=".\\BEAD_ploys",layer="testr3_1_r4_ndvi_90_2")
r3_firepoly2 = readOGR(dsn=".\\BEAD_ploys",layer="testr3_4_r4_ndvi_90_2")

# part 2: MODIS modis fire hotspot in this area to validate the results
#download data from http://modis-fire.umd.edu/pages/ActiveFire.php?target=GetData
# ftp://fuoco.geog.umd.edu/

#decompress the zip files
install.packages("R.utils")

library(R.utils)

# file.names = list.files(path = "./MCD14MLV5", pattern = "*.gz$")
file.names = list.files(path = "./MCD14MLV5", pattern = "*.asc$")

fire.df = c()

for (i in 1:length(file.names)){
#gunzip(paste("./MCD14MLV5/", file.names[i], sep = ""))

#file.tmp = substr(file.names[i], 1, nchar(file.names[i])-3)
file.tmp = data.frame(read.table(paste("./MCD14MLV5/", file.names[i], sep = ""), header = TRUE))

file.tmp = file.tmp[which(file.tmp$lon > 121 & file.tmp$lon < 127 & file.tmp$lat > 50 & file.tmp$lat < 53.5),]

fire.df = rbind(fire.df, file.tmp)

print(paste("Finishing for extracting ", i, " of ", length(file.names), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )

}

plot(density(fire.df$FRP), xlim = c(0,120))
plot(density(fire.df$conf), xlim = c(0,120))

#change into spatial database
library(sp)
fire.sp <- fire.df[which(fire.df$conf > 50),]
coordinates(fire.sp) <- ~lon+lat
projection(fire.sp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

#dxal <- readOGR(dsn = "D:/users/Zhihua/Landsat/XinganImages/boundry", layer = "dxal_bj_proj_polygon")
proj.utm = projection(dxal)
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "
dxal.geo = spTransform(dxal, CRS(proj.geo))

fire.sp.hz = fire.sp[dxal.geo[3,], ]

fire.sp.hz = spTransform(fire.sp.hz, CRS(proj.utm))

#crop into test area
fire.sp.hz1 = fire.sp.hz[r1, ]
fire.sp.hz3 = fire.sp.hz[r3, ]

######## plot final ######
#par(mfrow=c(2,1),mar=c(0,0,0,0))

png("region1-2.png",height = 3000, width = 3000, res = 300, units = "px")

plot(ba.grd1, legend=FALSE, axes=FALSE, box=FALSE)
plot(r1_firepoly2[!is.na(r1_firepoly2$dist_time),], border = "black", lwd = 2, add = TRUE)
plot(fire.sp1[1,], border = "black", lwd = 2, add = TRUE)
#plot(fire.sp1, border = "red", lwd = 2, add = TRUE)
plot(fire.sp.hz1, add = TRUE, pch = 20, col = "blue")

scalebar(10000, xy=c(505000, 5742500), type='bar', divs=4,below = "Meter")
#text(x=510000, y=5767000, "Region 1", cex = 1.5)
dev.off()

png("region2-2.png",height = 3000, width = 4500, res = 300, units = "px")

plot(ba.grd3, legend=FALSE, axes=FALSE, box=FALSE)
plot(r3_firepoly2[!is.na(r3_firepoly2$dist_time),], border = "black", lwd = 2, add = T)
plot(fire.sp3[1,], border = "black", lwd = 2, add = T)
#plot(fire.sp3, border = "red", lwd = 2, add = T)
plot(fire.sp.hz3, add = TRUE, pch = 20, col = "blue")

scalebar(10000, xy=c(498000, 5680000), type='bar', divs=4,below = "Meter")
#text(x=500000, y=5700000, "Region 2", cex = 1.5)
dev.off()


### 4/21/2016 validate BEAD-derived polygon with MODIS fire products ###################################
########################################################################################################
