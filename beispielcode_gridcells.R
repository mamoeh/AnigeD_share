library(sf)
library(sdcSpatial) # (nur für Beispiel-Datensatz)

data("dwellings") # Beispiel-Datensatz

# Erst muss man herausfinden, in welchem Koordinatenreferenzsystem (CRS) der
# Geocoder seine Daten ausspuckt. Im Beispieldatensatz sieht man in der Info, 
# dass es hier EPSG:28992 ist: https://epsg.io/28992
# Für Angaben in Längen-/Breitengraden dürfte es hingegen typischerweise so etwas
# wie WGS84 sein (o.ä.): https://epsg.io/4326

# 1. Koordinaten zu sf-Objekt
dwellings <- st_as_sf(dwellings, coords = c("x", "y"), crs = "EPSG:2899")

# Für die Grid-Aggregation wird im ESS einheitlich ETRS89-LAEA verwendet, das den
# EPSG-Code 3035 hat: https://epsg.io/3035

# 2. Koordinaten-Transformation
dwellings_neu <- st_transform(dwellings, crs = "EPSG:3035")
st_crs(dwellings_neu) # ist "ETRS89-extended / LAEA Europe"

# 3. sf-Objekt zurück zu x,y
dwellings_xy <- as.data.frame(st_coordinates(dwellings_neu))

# 4. Gitterzelle berechnen durch Runden
gz_size <- 100 # z.B. 100m x 100m

dwellings_xy$gz_x <- floor(dwellings_xy$X / gz_size) * gz_size
dwellings_xy$gz_y <- floor(dwellings_xy$Y / gz_size) * gz_size

# 5. Mittelpunkt bestimmen
dwellings_xy$mp_x <- dwellings_xy$gz_x + (gz_size / 2)
dwellings_xy$mp_y <- dwellings_xy$gz_y + (gz_size / 2)

plot(dwellings_xy$X, dwellings_xy$Y)       # nicht aggregiert
plot(dwellings_xy$mp_x, dwellings_xy$mp_y) # aggregiert

