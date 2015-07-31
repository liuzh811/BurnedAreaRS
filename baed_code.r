#3/9/2015
#by Zhihua Liu: liuzh811@gmail.com

# 1. this code aims at mapping burned area at high spatial and temporal resolution using 
#    Landsat and MODIS MOD13Q1 vegetation indices data
# 2. Landsat image should be normalized (Canty and Neilsen et al (2008)) to extract burned area
# 3. create a "data"  folder to store MODIS and Landsat data, as well as study area polygons
# 4. create a "result" folder to store fire polygon

############## section 1: read MODIS data and prepare functions
setwd("D:\\GitHub\\BurnedAreaRS")

# load libraries
library(rgdal)
library(raster)
library(dismo)
library(gbm)
library(rgeos)

#read test area polygon
testr1 <- readOGR(dsn = ".\\data", layer = "test_area") #read test area polygon file
proj.utm = projection(testr1)
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

#read MODIS NDVI data 
ref <- list.files(path=".\\data", pattern="*NDVI.tif$") #get NDVI data name
YearDOYmodis = substr(ref, nchar(ref[1])-28, nchar(ref[1])-22)  #get MODIS date
Yearmodis = as.numeric(substr(YearDOYmodis,1,4)) #get year
DOYmodis = as.numeric(substr(YearDOYmodis,5,7)) #get day of year
s <- stack(paste(".\\data", ref, sep = "")) #read NDVI data as raster stack

#read NDVI quality data info
ref.qa <- list.files(path=".\\data", pattern="*VI_Quality.tif$") #get NDVI data quality name
s.qa <- stack(paste(".\\data\\", ref.qa, sep = "")) #read NDVI quality data as raster stack

################# define two functions for future use  ######################
# function 1: find the MODIS VI good measurement 
QAfind.mt = function(m){ ##m is a matrix, row is points, and col is data\e
  nc = ncol(m);nr = nrow(m)
  #change value into bit
  qc.bit = apply(m,1,function(x){paste(rev(as.integer(intToBits(rev(x)))), collapse="")})
  bit.start = 31+32*(0:(nc-1))
  bit.end = 32+32*(0:(nc-1))
  #return a list, each element is a points
  qc.bit.band1 = lapply(qc.bit,function(x){strtoi(substring(x, bit.start, bit.end),2)}) 
  qc = unlist(qc.bit.band1)
  q.mt = matrix(qc, nrow = nr, ncol = nc, byrow = TRUE)
  return(q.mt)
}
# function 2: find start time for persisent NDVI change
findmaxtime4 = function(x, threshold1 = 3){ 
  x = as.numeric(x)
  x[which(is.na(x))]=2
  #fill values if there are only one or two missing values
  for(i in 1:(length(x)-4)){
      if (length(which(x[i:(i+4)]==2))<5){ #judge whether there are two many missing values
      if (x[i]==1 & x[i+1]==1 & x[i+2]==0 & x[i+3]==1 & x[i+4]==1) {x[i+2]=1}
      else if (x[i]==1 & x[i+1]==0 & x[i+2]==1 & x[i+3]==1 & x[i+4]==1) {x[i+1]=1}
    } #end of #judge whether there are two many missing values
    
  } #end of i
  
  #Find where the values change:
  diffs <- x[-1L] != x[-length(x)]
  #Get the indexes, and then get the difference in subsequent indexes:
  idx <- c(which(diffs), length(x))
  idx.1 = which(x[idx]==1) #find the location
  idx.0 = which(x[idx]==0) #find the location
  
  diff1 = diff(c(0, idx))[idx.1] #find numbers of consecutive 1
  max.dif = which(diff1>=threshold1)
  #max.dif.idx = rev(sort.int(max.dif, index.return=TRUE)$ix)
  #max.dif = max.dif[max.dif.idx]
   if(length(max.dif)>0){
    max.loc = idx[idx.1[max.dif]-1]+1
  } else {
    max.loc = "NA"
  }
  return(max.loc)
  }

############## section 2: read Landsat data and preprocessing
fn = c("lt2010229_p122r20_sub","le2011224_p122r20_sub") #Landsat file name
YearTM = as.numeric(substr(fn, 3,9)) #get year
testr1@data$area = sapply(slot(testr1, "polygons"), slot, "area")

# Read LANDSAT imageries at ENVI format and Create a bunch of list to store spectral indices
nbr.list = list() #NBR 
ndvi.list = list() #NDVI
band1.list = list()
band2.list = list()
band3.list = list()
band4.list = list()
band5.list = list()
band6.list = list()

tc1.list = list() #TC brightness
tc2.list = list() #TC greenness
tc3.list = list() #TC wetness
di.list = list() #disturbance index
cloud.list = list() #cloud mask
water.list = list() #water mask
for (i in 1:length(fn)){
  
  r = readGDAL(paste(".\\data\\", fn[i], sep = "")) #r is a spatialdataframe
  
  r1 <-raster(
    matrix(r@data@.Data[[1]], nrow = r@grid@cells.dim[2], ncol = r@grid@cells.dim[1], byrow = TRUE),
    xmn=r@bbox[1,1], xmx=r@bbox[1,2],
    ymn=r@bbox[2,1], ymx=r@bbox[2,2], 
    crs=CRS(proj.utm)
  )
  
  r2 <-raster(
    matrix(r@data@.Data[[2]], nrow = r@grid@cells.dim[2], ncol = r@grid@cells.dim[1], byrow = TRUE),
    xmn=r@bbox[1,1], xmx=r@bbox[1,2],
    ymn=r@bbox[2,1], ymx=r@bbox[2,2], 
    crs=CRS(proj.utm)
  )
  
  r3 <-raster(
    matrix(r@data@.Data[[3]], nrow = r@grid@cells.dim[2], ncol = r@grid@cells.dim[1], byrow = TRUE),
    xmn=r@bbox[1,1], xmx=r@bbox[1,2],
    ymn=r@bbox[2,1], ymx=r@bbox[2,2], 
    crs=CRS(proj.utm)
  )
  
  r4 <-raster(
    matrix(r@data@.Data[[4]], nrow = r@grid@cells.dim[2], ncol = r@grid@cells.dim[1], byrow = TRUE),
    xmn=r@bbox[1,1], xmx=r@bbox[1,2],
    ymn=r@bbox[2,1], ymx=r@bbox[2,2], 
    crs=CRS(proj.utm)
  )
  r5 <-raster(
    matrix(r@data@.Data[[5]], nrow = r@grid@cells.dim[2], ncol = r@grid@cells.dim[1], byrow = TRUE),
    xmn=r@bbox[1,1], xmx=r@bbox[1,2],
    ymn=r@bbox[2,1], ymx=r@bbox[2,2], 
    crs=CRS(proj.utm)
  )
  r6 <-raster(
    matrix(r@data@.Data[[6]], nrow = r@grid@cells.dim[2], ncol = r@grid@cells.dim[1], byrow = TRUE),
    xmn=r@bbox[1,1], xmx=r@bbox[1,2],
    ymn=r@bbox[2,1], ymx=r@bbox[2,2], 
    crs=CRS(proj.utm)
  )
  st = stack(r1,r2,r3,r4,r5,r6)
  st = crop(st, testr1)
  
  ## st[st < 0] <- NA #this process is very slow
  ## st[st > 9999] <- NA #this process is very slow
  #calculate NDVI
  ndvi = (st[[4]]-st[[3]])/(st[[4]]+st[[3]])
  
  #calculate NBR
  nbr = (st[[5]]-st[[6]])/(st[[5]]+st[[6]])
  
  #calculate TC index
  if(substr(fn[i],1,2) == "le"){
    brightness = 0.3561*st[[1]]+0.3972*st[[2]]+0.3904*st[[3]]+0.6966*st[[4]]+0.2286*st[[5]]+0.1596*st[[6]]  #brightness
    greenness = -0.3344*st[[1]]-0.3544*st[[2]]-0.4556*st[[3]]+0.6966*st[[4]]-0.0242*st[[5]]-0.2630*st[[6]]  #greenness
    wetness = 0.2626*st[[1]]+0.2141*st[[2]]+0.0926*st[[3]]+0.0656*st[[4]]-0.7629*st[[5]]-0.5388*st[[6]]  #wetness
  } else if (substr(fn[i],1,2) == "lt"){
    brightness = 0.3037*st[[1]]+0.2793*st[[2]]+0.4343*st[[3]]+0.5585*st[[4]]+0.5082*st[[5]]+0.1863*st[[6]]  #brightness
    greenness = -0.2848*st[[1]]-0.2435*st[[2]]-0.5436*st[[3]]+0.7246*st[[4]]+0.0840*st[[5]]-0.18*st[[6]]  #greenness
    wetness = 0.1509*st[[1]]+0.1793*st[[2]]+0.3299*st[[3]]+0.3406*st[[4]]-0.7112*st[[5]]-0.4572*st[[6]]  #wetness
  } else if (substr(fn[i],1,2) == "lc"){
    brightness = 0.3029*st[[1]]+0.2786*st[[2]]+0.4733*st[[3]]+0.5599*st[[4]]+0.508*st[[5]]+0.1872*st[[6]]  #brightness
    greenness = -0.2941*st[[1]]-0.243*st[[2]]-0.5424*st[[3]]+0.7276*st[[4]]+0.0713*st[[5]]-0.1608*st[[6]]  #greenness
    wetness = 0.1511*st[[1]]+0.1973*st[[2]]+0.3283*st[[3]]+0.3407*st[[4]]-0.7117*st[[5]]-0.4559*st[[6]]  #wetness
    }
  tc1 <- brightness
  tc2 <- greenness
  tc3 <- wetness
  
  #calculate disturbance index
  # find mature forest, 
  mat.for = ndvi > 0.8 & ndvi < 1
  if (length(unique(mat.for)) == 1 & unique(mat.for) == FALSE) {
    mat.for[] = -9999
    di = mat.for 
  } else {
    tc1.t = mat.for*tc1
    tc1.t[tc1.t==0] = NA
    tc1.t.mn = mean(getValues(tc1.t), na.rm = TRUE)
    tc1.t.sd = sd(getValues(tc1.t), na.rm = TRUE)
    tc1.t = (tc1-tc1.t.mn)/tc1.t.sd
    
    tc2.t = mat.for*tc2
    tc2.t[tc2.t==0] = NA
    tc2.t.mn = mean(getValues(tc2.t), na.rm = TRUE)
    tc2.t.sd = sd(getValues(tc2.t), na.rm = TRUE)
    tc2.t = (tc2-tc2.t.mn)/tc2.t.sd
    
    tc3.t = mat.for*tc3
    tc3.t[tc3.t==0] = NA
    tc3.t.mn = mean(getValues(tc3.t), na.rm = TRUE)
    tc3.t.sd = sd(getValues(tc3.t), na.rm = TRUE)
    tc3.t = (tc3 -tc3.t.mn)/tc3.t.sd
    
    di = tc1.t - tc2.t - tc3.t
  }
  
  #calculate NDSI to flag potential cloud
  ndsi = (st[[2]]-st[[5]])/(st[[2]]+st[[5]])
  MeanVis = (st[[1]]+st[[2]]+st[[3]])/3
  whiteness=abs((st[[1]]-MeanVis)/MeanVis)+abs((st[[2]]-MeanVis)/MeanVis)+abs((st[[3]]-MeanVis)/MeanVis) < 0.7
  b45 =  st[[4]]/st[[5]]>0.75
    
  cloud.t = st[[6]] > 3000 & ndsi<0.8 & ndvi < 0.8 & whiteness==1 & b45==1
  
  
  #water mask
  water = (ndvi<0.01 & st[[4]]<110)|(ndvi<0.1 & st[[4]]<50)
  
  nbr.list[[i]] = nbr
  ndvi.list[[i]] = ndvi
  band4.list[[i]] = st[[4]]
  band5.list[[i]] = st[[5]]
  tc1.list[[i]] = tc1
  tc2.list[[i]] = tc2
  tc3.list[[i]] = tc3
  di.list[[i]] = di
  cloud.list[[i]] = cloud.t
  water.list[[i]] = water
  
  print(paste("Finish Listing files ", i," of ", length(fn), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
  rm(list = c("r","r1","r2","r3","r4","r5","r6","brightness","greenness","wetness","st"))
}

#soomth spectral indices using a 3*3 window size  
nbr.list2 = list()
tc2.list2 = list()
tc3.list2 = list()
ndvi.list2 = list()
di.list2 = list()
b4.list2 = list()
b5.list2 = list()
cloud.list2 = list()
water.list2 = list()
for(i in 1:length(nbr.list)){
  if(is.null(di.list[[i]])==FALSE){
    nbr.t = focal(nbr.list[[i]], w= matrix(1/9, nc=3, nr=3), fun=median)
    nbr.list2[[i]]<-nbr.t
    
    tc2.t = focal(tc2.list[[i]], w= matrix(1/9, nc=3, nr=3), fun=median)
    tc2.list2[[i]]<-tc2.t
    
    tc3.t = focal(tc3.list[[i]], w= matrix(1/9, nc=3, nr=3), fun=median)
    tc3.list2[[i]]<-tc3.t
    
    ndvi.t = focal(ndvi.list[[i]], w= matrix(1/9, nc=3, nr=3), fun=median)
    ndvi.list2[[i]]<-ndvi.t
    
    di.t = focal(di.list[[i]], w= matrix(1/9, nc=3, nr=3), fun=median)
    di.list2[[i]]<-di.t
    
    
    b4.t = focal(band4.list[[i]], w= matrix(1/9, nc=3, nr=3), fun=median)
    b4.list2[[i]]<-b4.t
    
    b5.t = focal(band5.list[[i]], w= matrix(1/9, nc=3, nr=3), fun=median)
    b5.list2[[i]]<-b5.t
    
  }
  rm(list = c("nbr.t","ndvi.t","di.t","b4.t","b5.t"))
  print(paste("Finish calculating moving windows ", i," of ", length(di.list), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

cloud.list2 = cloud.list
water.list2 = water.list

#calculate spectral indices differece between consective time
nbr.list2.dif = list()
tc2.list2.dif = list()
tc3.list2.dif = list()
ndvi.list2.dif = list()
di.list2.dif = list()
b4.list2.dif = list()
b5.list2.dif = list()

for(i in 2:length(nbr.list2)){
  nbr.list2.dif[[i-1]] <- nbr.list2[[i]] - nbr.list2[[i-1]]
  tc2.list2.dif[[i-1]] <- tc2.list2[[i]] - tc2.list2[[i-1]]
  tc3.list2.dif[[i-1]] <- tc3.list2[[i]] - tc3.list2[[i-1]]
  ndvi.list2.dif[[i-1]] <- ndvi.list2[[i]] - ndvi.list2[[i-1]]
  di.list2.dif[[i-1]] <- di.list2[[i]] - di.list2[[i-1]]
  b4.list2.dif[[i-1]] <- b4.list2[[i]] - b4.list2[[i-1]]
  b5.list2.dif[[i-1]] <- b5.list2[[i]] - b5.list2[[i-1]]
  
  print(paste("Finish calculating difference ", i," of ", length(nbr.list2), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
}

############## section 3: extract burned polygons from LANDSAT data
set.seed(1000)

#set threshold for NBR, NDVI, DI
nbr.threshold = -0.05
ndvi.threshold = -0.05
di.threshold = 1
poly.list = list() #store burned patches, if more than 2 times of landsat data is available
for (i in 1:length(nbr.list2.dif)){
#3.1 find burned core pixels and unburned pixels
 ## select burned core pixels
	  burned.core = nbr.list2.dif[[i]] < nbr.threshold & ndvi.list2.dif[[i]] < ndvi.threshold & di.list2.dif[[i]] > di.threshold
	  burned.core[cloud.list2[[i]]==1 | cloud.list2[[i+1]]==1] = 0 #remove cloud influence
    
##select unburned pixels
      unburned = nbr.list2.dif[[i]] > 0 &  ndvi.list2.dif[[i]] > 0
      unburned[cloud.list2[[i]]==1 | cloud.list2[[i+1]]==1] = 0 #remove cloud influence
      
#3.2 using brt model to simulate burned probability
#3.2.1 preparing dataset   
      
      burned.core[burned.core == 0] = NA
      unburned[unburned == 0] = NA 
      r1.xy = rasterToPoints(burned.core) #get raster coordinate, from left to right, top to bottom
      r1.xy = data.frame(r1.xy)
      r1.xy.sp <- SpatialPoints(coords = cbind(r1.xy$x,r1.xy$y),proj4string = CRS(proj.utm))
      
      r1.xy2 = rasterToPoints(unburned) #get raster coordinate, from left to right, top to bottom
      r1.xy2 = data.frame(r1.xy2)
      r1.xy2.sp <- SpatialPoints(coords = cbind(r1.xy2$x,r1.xy2$y),proj4string = CRS(proj.utm))
      
#preparing environmental variables
      delta.ndvi = ndvi.list2.dif[[i]]
      delta.nbr = nbr.list2.dif[[i]]
      delta.di = di.list2.dif[[i]]
      delta.tc2 = tc2.list2.dif[[i]]
      delta.tc3 = tc3.list2.dif[[i]]
      delta.b4 = b4.list2.dif[[i]]
      delta.b5 = b5.list2.dif[[i]]
      predictors = stack(delta.ndvi, delta.nbr, delta.di, delta.tc2, delta.tc3, delta.b4,delta.b5)
      names(predictors) = c("ndvi", "nbr", "di", "tc2", "tc3", "b4","b5")
      
      r1.xy.sp.df = extract(predictors, r1.xy.sp)
      r1.xy.sp.df = data.frame(r1.xy.sp.df, burned = 1)
      
      r1.xy2.sp.df = extract(predictors, r1.xy2.sp)
      r1.xy2.sp.df = data.frame(r1.xy2.sp.df, burned = 0)
      
      dat.df = rbind(r1.xy.sp.df, r1.xy2.sp.df)
      
      dat.df1 = rbind(dat.df[sample(nrow(dat.df[dat.df$burned==1,]), 1000), ],dat.df[sample(nrow(dat.df[dat.df$burned==0,]), 2000), ]) 
#3.2.2 brt modeling
      brt1=gbm.step(data = dat.df1, gbm.x = 1:(ncol(dat.df1)-1),
                    gbm.y = ncol(dat.df1), 
                    family = "bernoulli",
                    tree.complexity = 3, 
                    n.trees = 25, learning.rate = 0.01, bag.fraction = 0.5)
  
  p <- predict(predictors, brt1, n.trees=brt1$gbm.call$best.trees, type="response")	
  
#find the threshold to convert p into a binary raster
  preds = predict(brt1, dat.df1, n.trees = brt1$gbm.call$best.trees, type="response")
  d <- cbind(dat.df1$burned, preds)
  pres <- d[d[,1]==1, 2]
  abse <- d[d[,1]==0, 2]
  e <- dismo::evaluate(p=pres, a=abse)
  p1 = p > threshold(e)$kappa 
      
#3.2.3 doing a 7 by 7 windows smoothing
      p2 = focal(p1, w= matrix(1/49, nc=7, nr=7))
      p2 = p2>=0.5
      p2[p2 == 0] = NA
      ##change to polygon
      p2.poly = rasterToPolygons(p2, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
      p2.poly.copy = p2.poly  #backup a copy of
      
#3.2.4 remove holes from the polygon
#doing some plot to see if there are hole in the polygons
      plot(p2.poly)
      for(ii in 1:length(p2.poly@polygons[[1]]@Polygons)){
        
        polygon(p2.poly@polygons[[1]]@Polygons[[ii]]@coords, col = ii)
        text(x = p2.poly@polygons[[1]]@Polygons[[ii]]@labpt[1], y = p2.poly@polygons[[1]]@Polygons[[ii]]@labpt[2],
             labels = p2.poly@polygons[[1]]@Polygons[[ii]]@hole)
        
      }
      
#remove hole in the polygons 
      j = 1
      n = length(p2.poly@polygons[[1]]@Polygons)
      repeat{
        #for(j in 1:length(patch.poly3.list2.sp@polygons[[1]]@Polygons)){
        if (p2.poly@polygons[[1]]@Polygons[[j]]@hole == "TRUE"){
          p2.poly@polygons[[1]]@Polygons[[j]] <- NULL
          n = n-1
          j = j-1
        } else {j = j+1}
        if (j>=n) break
      }
      
#break into different polygons
      p = 0
      patch.poly5.list2 = list()
      for(j in 1:length(p2.poly@polygons[[1]]@Polygons)){
        p = p+1
        patch.poly5.list2[[p]] = Polygons(list(Polygon(p2.poly@polygons[[1]]@Polygons[[j]]@coords)),ID=p)
        
      }
      
      polys <- SpatialPolygons(patch.poly5.list2, proj4string=CRS(proj.utm))
      dat.df = data.frame(value=1:length(patch.poly5.list2), row.names=1:length(patch.poly5.list2)) 
      polys = SpatialPolygonsDataFrame(polys, dat.df)
      
      polys@data$area = sapply(slot(polys, "polygons"), slot, "area")
      
# remove area smaller than 10 cells 
      polys = polys[which(polys@data$area>62500), ]  
     
# mark the potential cloud contaminated polygons   
      cloud.t = cloud.list2[[i]]==1 | cloud.list2[[i+1]]==1
      cloud.t[cloud.t==0] = NA
      
      polys@data$cloud = 0  #0 is outside cloud
      if(freq(cloud.t)[which(is.na(freq(cloud.t)[,1])),2] < ncell(p1)){
        cloud.poly = rasterToPolygons(cloud.t, n=4, na.rm=TRUE, digits=12, dissolve=TRUE)
        
        for (p2 in 1:length(polys)){
          tmp = polys[p2,]
          
          if(is.null(gIntersection(cloud.poly, tmp))){polys@data$cloud[p2] = 0}
          else if (gContains(cloud.poly, tmp) == TRUE) {polys@data$cloud[p2] = 2} # complete within cloud 
          else {polys@data$cloud[p2] = 1}  # on cloud boundary
        }
      
      }
  

  # plot results, 0-clear; 1-extract burned polygon is on the cloud boundary; 2-extract burned polygon is complete within the clouds
plot(polys)
polygonsLabel(polys, labels = polys@data$cloud, method = "centroid", cex = 1.5, col = "green", pch=19)
plot(cloud.poly, add = TRUE, col = "red")

  
#assign the time frame for disturbance from landsat data
      polys@data$StartDate =YearTM[i] 
      polys@data$EndDate =YearTM[i+1] 
      
      poly.list[[i]] = polys   
   
  print(paste("Finish extracting ", i," of ", length(nbr.list2.dif), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
  
}

#remove empty elements for the poly.list
poly.list2 = poly.list

if(length(which(sapply(poly.list,is.null),arr.ind=TRUE))>=1 ) {
  poly.list = poly.list[-which(sapply(poly.list,is.null),arr.ind=TRUE)]
}
#

##save the poly into a shapefile
for (i in 1:length(poly.list)){
    writeOGR(poly.list[[i]], ".\\result", "Fire_polys_without_date", driver="ESRI Shapefile", overwrite_layer=TRUE)
  
}

############## section 4: assign disturbance date from MODIS dataset
# 4.1 need to find persistent forest from landsat and use it to construct a model from MODIS data
# find persistent forest from landsat
set.seed(1001)	
ndvi2 = ndvi.list[[1]]>0.6 & ndvi.list[[1]] < 1 & ndvi.list[[2]] > 0.6 & ndvi.list[[2]] < 1 
ndvi2[ndvi2==0] = NA
ndvi2.pts = rasterToPoints(ndvi2) 
ndvi2.pts = data.frame(ndvi2.pts)

if (nrow(ndvi2.pts) > 5000){ndvi2.pts = ndvi2.pts[sample(1:nrow(ndvi2.pts), 5000), ]}

ndvi2.pts <- SpatialPoints(coords = cbind(ndvi2.pts$x,ndvi2.pts$y),proj4string = CRS(proj.utm))
ndvi2.df = data.frame(t(extract(s, ndvi2.pts)))

#read NDVI quality data info
ndvi2.qa.df = data.frame(t(extract(s.qa, ndvi2.pts)))
ndvi2.qa.df2 <-  QAfind.mt(ndvi2.qa.df)
bad.idx2 = which(ndvi2.qa.df2 > 0, arr.ind = TRUE)  #0 is good measurement, 1 is useable, but check other QA infos, 
ndvi2.df[bad.idx2] = NA #relpamce bad measruement with NA

ndvi2.df = data.frame(Year = Yearmodis, DOY = DOYmodis, value = apply(ndvi2.df,1,median,na.rm = T))

## second, fit a OLS model, and accout for seasonal fluctions
evi.ols = lm(value~ cos(2*3.14195*DOY/360)+sin(2*3.14195*DOY/360)+DOY, data = ndvi2.df)
df1.pred = predict(evi.ols, newdata = ndvi2.df,interval="predict", level = 0.75)
#matplot(df1.pred, type = "b", pch = c(1,2,2), col = c(1,2,2))
#matplot(df1.pred2, type = "b", pch = c(1,2,2), col = c(1,2,2))

# 4.2 assign disturbance date for each burned polygons
  firep <- polys
  firep@data$Fireid = 1:nrow(firep)
  firep@data$area = sapply(slot(firep, "polygons"), slot, "area")  #calculate area
  firep@data$dist_time = 0
  
  for(i in 1:length(firep@data$Fireid)){
    
    firep.tmp = firep[i,]
    
    #sampling data points based on based on burned area
    Rpts1 = spsample(firep.tmp, n=100, type='regular')
  
    Rpts1.v = extract(s,Rpts1) # Extract the EVI data for the available two layers from the generated stack
    Rpts1.qa.v = extract(s.qa,Rpts1) # Extract the EVI data quality informaiton 
    
   ##find only good measurement
    Rpts1.qa.v2 = QAfind.mt(Rpts1.qa.v)
    bad.idx2 = which(Rpts1.qa.v2 > 0, arr.ind = TRUE)  #0 is good measurement, 1 is useable, but check other QA infos, 
    Rpts1.v[bad.idx2] = NA #relpamce bad measruement with NA
    
    dat = data.frame(Burn = apply(Rpts1.v,2,median, na.rm = TRUE), df1.pred)
    
    # matplot(dat, col=c(1,2,3,3), lty = c(1,1,2,2), lwd = c(2,2,1,1),type = "b", ylab = "predicted y")
    # if it is significant bigger for more than a certain period of time, then, it should be fire, otherwise, not fire
    dat.dif = dat$lwr-dat$Burn

    dat.dif[which(dat.dif>0)] = 1
    dat.dif[which(dat.dif<=0)] = 0
    
    dist_time = as.numeric(YearDOYmodis[findmaxtime4(dat.dif, threshold = 3)])
    
    #if there are multiple date, then, check which one is in the landsat time frame
    dist_time_lt = which(dist_time>=firep@data$StartDate[i]&dist_time<=firep@data$EndDate[i])
    if(length(dist_time_lt)>0){
      
      firep@data$dist_time[i] = dist_time[dist_time_lt[1]]-8  ##assume the disturbance occur between the two date
      
    } else {firep@data$dist_time[i] = NA}
    
    
  print(paste("Finish find disturbance time ", i," of ", length(firep@data$Fireid), " at ", format(Sys.time(), "%a %b %d %X %Y"), sep = " ") )
    
  }
  
  
writeOGR(firep, ".\\result", "Fire_polys_with_date", driver="ESRI Shapefile", overwrite_layer=TRUE)
  
## plot results
plot(ndvi.list2.dif[[1]], col = colorRampPalette(c("red", "blue", "green"))(255))
plot(firep, add = TRUE, border = "red", lwd = 2)
polygonsLabel(firep, labels = firep@data$dist_time, method = "centroid", cex = 1.5, col = "green", pch=19)

