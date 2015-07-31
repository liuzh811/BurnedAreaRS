#7/31/2015 by liuzh811@gmail.com
#download MODIS VI data for BAED 

# load libraries
library(MODIS)
library(rgdal)
library(raster)

#set spatial extent
testr1 <- readOGR(dsn = ".\\data", layer = "test_area")
proj.utm = projection(testr1)
proj.geo = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0 "

#set MODIS path, see MODIS package for details
MODISoptions(localArcPath=".",
             outDirPath=".",
             gdalPath='c:/OSGeo4W64/bin')

#set date period for data
dates <- as.POSIXct( as.Date(c("1/5/2010","1/11/2011"),format = "%d/%m/%Y") )
dates2 <- transDate(dates[1],dates[2]) # Transform input dates from before

# getProduct() # list available MODIS products
#download MOD13Q1 data
runGdal(product="MOD13Q1",  #VI/combined/Tile/500m/monthly
        begin=dates2$beginDOY, #start date
        end = dates2$endDOY,#end date
        tileH = 25,tileV = 3, #tile
        SDSstring = "111", #extract the first 3 layers
        extent = testr1, #crop to extent
        outProj=proj.utm, #reproject to UTM
        job = "data") #download folder
