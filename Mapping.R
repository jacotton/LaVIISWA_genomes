library(tmap)
library(sf)
library(sp)

#Rivers
#Protected
#Lakes
#UgandaAdmin1
#Roads
#Urban
#Countryboudnars
#LVislands

tmap_mode("plot") 
tmap_options(check.and.fix=TRUE)

#from https://data.humdata.org/dataset/cod-ab-uga
Uganda_Admin1 <- read_sf("uga_admbnda_ubos_20200824_shp/uga_admbnda_adm1_ubos_20200824.shp")

#Country boundaries shape file from https://open.africa/dataset/africa-shapefiles
#version dated April 21, 2020
CountryBoundaries <- read_sf("afr_g2014_2013_0/afr_g2014_2013_0.shp")

#protected areas from:www.protectedplanet.net
ProtectedAreas <- read_sf("WDPA_June2017_UGA-shapefile-polygons/WDPA_June2017_UGA-shapefile-polygons.shp")

#islands from:
#https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/E0515D
LVislands <- read_sf("LVIslands_dataverse_files/LV_Islands_Polygon.shp")


#all below are from https://www.naturalearthdata.com

#download week of 27 July 2022. Version numbers in manuscript.
Urban <- read_sf("ne_10m_urban_areas/ne_10m_urban_areas.shp")

Rivers <- read_sf("ne_10m_rivers_lake_centerlines")

Roads <- read_sf("ne_10m_roads/ne_10m_roads.shp")

Lakes <- read_sf("ne_10m_lakes/ne_10m_lakes.shp")



plotlabs <- data.frame(name=c("Kampala","Lake Victoria","Jinja","Uganda"),long=c(32.5825,33,33.2026,33),lat=c(0.3476,-0.8,0.4479,2.5),size=c(6,12,6,12))
coordinates(plotlabs) <- c("long", "lat")
proj4string(plotlabs) <- CRS("+proj=longlat +datum=WGS84")

schools <- data.frame(name=c("Bwondha","Bugoto","Musubi","Kocoge"),long=c(33.56138,33.62837,33.6652,34.23175),lat=c(0.17775,0.32369,0.31105,0.77552),type=c("shoreline","shoreline","shoreline","inland"))
coordinates(schools) <- c("long", "lat")
proj4string(schools) <- CRS("+proj=longlat +datum=WGS84")

#coordinates of 'zoomed' map on largescale map.
x1 <- c(32.640075,32.640075,32.853028,32.853028,32.640075)
y1 <- c(-0.150663,0.053817,0.053817,-0.150663,-0.150663)

# assign the vertices to a `polygon` 
poly1 <- sp::Polygon(cbind(x1,y1))

# This step combines the last two together - making Polygons and then SpatialPolygons
bb2.Poly <- sp::SpatialPolygons(list(sp::Polygons(list(poly1),ID = "A")))


Island_labels <-data.frame(name=c("Damba\nisland","Koome island","Lugumba\nisland"),long=c(32.789374,32.75183,32.830776),lat=c(0.018277,-0.087252,-0.084675))
coordinates(Island_labels) <- c("long", "lat")
proj4string(Island_labels) <- CRS("+proj=longlat +datum=WGS84")

Villages <- data.frame(name=c("Lugumba","Kitosi","Busi","Zingoola","Kisu","Kakeeka","Kachanga","Katooke"),long=c(32.826376,32.800875,32.782971,32.714677,32.768072,32.764769,32.799936,32.803371),lat=c(-0.067780,-0.065548,-0.086709,-0.055394,0.013277,-0.039974,-0.008995,0.001041),type=c("standard","intensive","intensive","standard","intensive","standard","standard","intensive")) 

Village_text <- data.frame(name=c("Lugumba","Kitosi","Busi","Zingoola","Kisu","Kakeeka","Kachanga","Katooke"),long=c(32.826376,32.800875,32.782971,32.714677,32.768072,32.764769,32.799936,32.803371),lat=c(-0.062780,-0.060548,-0.081709,-0.050394,0.018277,-0.044974,-0.013995,0.006041),type=c("standard","intensive","intensive","standard","intensive","standard","standard","intensive")) 

coordinates(Villages) <- c("long", "lat")
proj4string(Villages) <- CRS("+proj=longlat +datum=WGS84")
coordinates(Village_text) <- c("long", "lat")
proj4string(Village_text) <- CRS("+proj=longlat +datum=WGS84")

#two panels were subsequently combined in illustrator
cairo_pdf(file="UgandaMap1.pdf",onefile=FALSE,width = 9, height = 7)

 tm_shape(Rivers,bbox=bbox1) + tm_lines(col="lightblue")  + tm_shape(ProtectedAreas) + tm_fill(col="darkolivegreen3") + tm_shape(Lakes) + tm_fill(col="lightblue")  + tm_shape(Uganda_Admin1) + tm_borders(lty = "dotted",col="grey70") + tm_shape(Roads) + tm_lines(col="lightgoldenrod1") + tm_shape(Urban) + tm_fill(col="goldenrod1") +  tm_shape(CountryBoundaries) + tm_borders(col="grey50") + tm_shape(plotlabs) + tm_text("name",size="size") + tm_shape(schools) + tm_dots(size=0.5,col="type",palette=c("deeppink","royalblue4")) + tm_layout(legend.outside = TRUE,bg.color="cornsilk") + tm_scale_bar(position=c("left", "bottom")) + tm_shape(bb2.Poly) + tm_borders()

dev.off()

cairo_pdf(file="UgandaMap2.pdf",onefile=FALSE,width = 9, height = 7)

 tm_shape(LVislands,bbox=bbox2) + tm_fill(col="cornsilk") + tm_layout(bg.color="lightblue") + tm_shape(Island_labels) + tm_text("name",size=0.75) + tm_shape(Villages) + tm_dots(size=0.25,col="type",palette=c("black","red")) + tm_shape(Village_text) + tm_text("name",size=1) + tm_scale_bar(position=c("left", "bottom"))

dev.off()
