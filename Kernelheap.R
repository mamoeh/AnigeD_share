library(Kernelheaping)
library(ptable)
library(sf)
library(raster)
library(ggplot2)
library(data.table)
library(dplyr)

# ----- Helper Functions -----

###### Look up the perturbed value by cell key from a pTable ################# # 
# x_ck : length-2 vector c(cell count, cell key) for one table cell
# pt   : pTable from a ptable object
############################################################################## #

get_perturbation_value <- function(x_ck, pt) {
  
  if(is.na(x_ck[1])) { 
    NA 
  } else {
    max_x <- max(pt$i)
    pt$v[pt$i == min(max_x, x_ck[1]) & data.table::between(x_ck[2], pt$p_int_lb, pt$p_int_ub)]
  }
}

###### #
#
###### #

get_jaccard <- function(var) {
  
  hs_org <- hots_spop[hots_spop$version == "original", var]
  hs_low <- hots_spop[hots_spop$version == "low", var]
  hs_mid <- hots_spop[hots_spop$version == "mid", var]
  hs_hig <- hots_spop[hots_spop$version == "high", var]
  
  jacc_low <- sum(hs_org & hs_low) / sum(hs_org | hs_low)
  jacc_mid <- sum(hs_org & hs_mid) / sum(hs_org | hs_mid)
  jacc_hig <- sum(hs_org & hs_hig) / sum(hs_org | hs_hig)
  
  c(jacc_low, jacc_mid, jacc_hig)
}

##### Turn raster version of hot spots into spatial polygons for plotting #### #
# var:      Variable for hot spot identification
# versions: vector of up to 4 of c("original", "low", "mid", "high")
############################################################################## #

hs_raster_to_sf <- function(var, versions) {
  
  r_list <- vector(mode = "list", length = length(versions))
  
  for(i in 1:length(versions)) {
    
    ras <- rasterFromXYZ(xyz = hots_spop[hots_spop$version == versions[i],
                                         c("lon", "lat", var)],
                         crs = st_crs(Berlin))
    
    ras[ras == 0] <- NA
    
    r_list[[i]] <- test <- st_as_sf(rasterToPolygons(ras, n = 4, digits = 1, dissolve = TRUE))
  }
  
  df <- bind_rows(r_list)
  st_crs(df) <- st_crs(Berlin)
  
  allvers <- c("original", "low", "mid", "high")
  df[, var] <- factor(versions, levels = allvers[match(versions, allvers)],
                      labels = ifelse(versions == "original", "original", 
                                      paste("CKM", allvers[match(versions, allvers)])))
  df
}  


# ----- Load Data -----

dat_path <- "kernelheaping/RendtelEtAl_material/kernelheaping/data/"


# ShapeFile
Berlin <- sf::st_read(file.path(dat_path, "shapefile/RBS_OD_LOR_2015_12.shp"))
sf::st_read(file.path(dat_path, "Open_map_berlin"))

# Inhabitant data
EWRDemog <- read.csv2(file.path(dat_path, "EWR201512E_Matrix.csv"))
EWRMigra <- read.csv2(file.path(dat_path, "EWRMIGRA201512H_Matrix.csv"),
                      colClasses = list(RAUMID = "character")) %>%
  mutate(RAUMID_int = as.integer(RAUMID))

EWR <- merge(EWRMigra, dplyr::select(EWRDemog, -c(ZEIT, BEZ, PGR, BZR, PLR, STADTRAUM)),
             by.x = "RAUMID_int", by.y = "RAUMID")

Berlin <- merge(Berlin, dplyr::select(EWR, -c(RAUMID_int, BEZ, PGR, BZR)), 
                by.x = "PLR", by.y = "RAUMID")

rm(EWR, EWRDemog, EWRMigra)


# ----- Prepare ptables -----

ptab_low <- ptable::create_cnt_ptable(D = 4,  V = 2.5, js = 2)
ptab_mid <- ptable::create_cnt_ptable(D = 8,  V = 3.5, js = 2)
ptab_hig <- ptable::create_cnt_ptable(D = 12, V = 5.0, js = 2) 

ptab_low@pTable$version <- "low"
ptab_mid@pTable$version <- "mid"
ptab_hig@pTable$version <- "high"

ptab_low@pTable$grad <- ptab_low@pTable$i / max(ptab_low@pTable$i)
ptab_mid@pTable$grad <- ptab_mid@pTable$i / max(ptab_mid@pTable$i)
ptab_hig@pTable$grad <- ptab_hig@pTable$i / max(ptab_hig@pTable$i)

# validate
ptab_low@empResults
ptab_mid@empResults
ptab_hig@empResults

# label
ptab_low@pTable$label <- FALSE
ptab_mid@pTable$label <- FALSE
ptab_hig@pTable$label <- FALSE
ptab_low@pTable$label[match(unique(ptab_low@pTable$i), ptab_low@pTable$i)] <- TRUE
ptab_mid@pTable$label[match(unique(ptab_mid@pTable$i), ptab_mid@pTable$i)] <- TRUE
ptab_hig@pTable$label[match(unique(ptab_hig@pTable$i), ptab_hig@pTable$i)] <- TRUE

pt <- rbind(ptab_low@pTable, ptab_mid@pTable, ptab_hig@pTable) %>%
  filter(i > 0)
pt$version <- factor(pt$version, levels = c("low", "mid", "high"))

# visualize ptabs
#ggplot(pt, aes(group = i, color = grad)) +
#  geom_line(aes(v, p), show.legend = FALSE) +
#  geom_point(aes(v, p), show.legend = FALSE) +
#  #geom_text(data = pt[pt$label, ], aes(v - 0.2, p + 0.01,label = i)) +
#  scale_colour_gradientn(colors = rev(viridis::magma(16))) +
#  facet_wrap(~version) +
#  theme_dark() +
#  xlab("noise") +
#  ylab("prob.")
#ggsave(filename = "ptabs.png", height = 800, width = 2400, units = "px")


# modify step
ptab_low <- modify_cnt_ptable(ptab_low, threshold = 0.2, seed = 45879571)
ptab_mid <- modify_cnt_ptable(ptab_mid, threshold = 0.2, seed = 98137615)
ptab_hig <- modify_cnt_ptable(ptab_hig, threshold = 0.2, seed = 26846183)


# ----- Create Cell Keys & Apply Perturbation -----

# Since we do not have access to microdata, we draw cell keys directly as a
# second-best option

Berlin_ck <- Berlin_prot_low <- Berlin_prot_mid <- Berlin_prot_hig <- 
  st_drop_geometry(Berlin)

set.seed(12164910)

# for migration background
Berlin_ck[, 12:21] <- runif(length(12:21) * nrow(Berlin_ck), 0, 1)

# for age groups
Berlin_ck[, 25:65] <- runif(length(25:65) * nrow(Berlin_ck), 0, 1)
Berlin_ck$E_E <- apply(Berlin_ck[, 25:65], MARGIN = 1, sum) 

# for sexes
Berlin_ck$E_EW <- runif(nrow(Berlin_ck), 0, 1)
Berlin_ck$E_EM <- Berlin_ck$E_E - Berlin_ck$E_EW
Berlin_ck$E_E  <- Berlin_ck$E_E  - floor(Berlin_ck$E_E)
Berlin_ck$E_EM <- Berlin_ck$E_EM - floor(Berlin_ck$E_EM)


# Apply noise

for(i in 12:65) {
  
  Berlin_prot_low[, i] <- Berlin_prot_low[, i] +
    apply(cbind(Berlin_prot_low[, i], Berlin_ck[, i]), 1, 
          get_perturbation_value, 
          pt = ptab_low@pTable)
  
  Berlin_prot_mid[, i] <- Berlin_prot_mid[, i] +
    apply(cbind(Berlin_prot_mid[, i], Berlin_ck[, i]), 1, 
          get_perturbation_value,
          pt = ptab_mid@pTable)
  
  Berlin_prot_hig[, i] <- Berlin_prot_hig[, i] +
    apply(cbind(Berlin_prot_hig[, i], Berlin_ck[, i]), 1, 
          get_perturbation_value,
          pt = ptab_hig@pTable)
}

# mark
versions <- c("original", "low", "mid", "high")

Berlin$version <- versions[1]
Berlin_prot_low$version <- versions[2]
Berlin_prot_mid$version <- versions[3]
Berlin_prot_hig$version <- versions[4]

# restore geometry
st_geometry(Berlin_prot_low) <- st_geometry(Berlin)
st_geometry(Berlin_prot_mid) <- st_geometry(Berlin)
st_geometry(Berlin_prot_hig) <- st_geometry(Berlin)

# collect in single data.frame
Berlin <- rbind(Berlin, Berlin_prot_low, Berlin_prot_mid, Berlin_prot_hig)
Berlin$version <- factor(Berlin$version, levels = versions)


rm(Berlin_prot_low, Berlin_prot_mid, Berlin_prot_hig, Berlin_ck)


#ggplot(Berlin[Berlin$version != "original", ]) +
#  geom_histogram(aes(x = Noise_E_E65U80, fill = version), 
#                 color = "black", binwidth = 1) +
#  facet_wrap(~version)


# ----- Vignette Part 1 -----

centroids <- st_coordinates(st_centroid(Berlin[Berlin$version == "original", ]))

# set common parameters
burn <- 5
samp <- 10
grsz <- 325

# data.frame to collect results in
dens_spop <- data.frame(lon = vector(mode = "numeric", length = 4*grsz^2),
                        lat = vector(mode = "numeric", length = 4*grsz^2),
                        version = rep(versions, each = grsz^2))
dens_spop$version <- factor(dens_spop$version, levels = versions)

# for all demographic variables ...
for(i in colnames(Berlin)[12:65]) {
  
  dens_spop[, i] <- 0
  
  # ... estimate four times (from orginal data, from CKM data low, mid, high)
  for(j in versions) {
    
    # estimate density
    est <- dshapebivr(data = cbind(centroids,
                                   st_drop_geometry(Berlin[Berlin$version == j, i])),
                      burnin = burn, samples = samp,
                      adaptive = FALSE,
                      shapefile = as_Spatial(Berlin[Berlin$version == j, ]),
                      gridsize = grsz, boundary = TRUE)
    
    # store estimated density
    dens_spop[dens_spop$version == j, i] <- as.vector(est$Mestimates$estimate)
    
    # set negative density to zero
    nd <- dens_spop[dens_spop$version == j, i] < 0
    dens_spop[dens_spop$version == j, i][nd] <- 0
    
    print(paste("Done:", i, j))
  }
}

# add grid point coordinates
grpts <- expand.grid(est$Mestimates$eval.points[[1]], est$Mestimates$eval.points[[2]])
dens_spop$lon <- rep(grpts[, 1], 4)
dens_spop$lat <- rep(grpts[, 2], 4)

#save(dens_spop, file = "dens_spop.RData")


## identify hotspots

hots_spop <- dens_spop
qtl <- 0.9

for(i in colnames(hots_spop)[4:57]) {
  
  for(j in versions) {
    
    dens <- dens_spop[dens_spop$version == j, i]
    
    hs_qtl <- quantile(dens[dens > 0], qtl)
    
    in_hs <- dens_spop[dens_spop$version == j, i] >= hs_qtl

    
    hots_spop[hots_spop$version == j, i][in_hs] <- 1
    hots_spop[hots_spop$version == j, i][!in_hs] <- 0
  }
}

## Calculate Jaccard similarities

jaccSim <- data.frame(variable = rep(colnames(Berlin)[12:65], 3),
                      version  = rep(c("low", "mid", "high"), each = 54))

jacc <- t(sapply(colnames(Berlin)[12:65], FUN = get_jaccard))
jaccSim$jaccSim <- round(c(jacc[, 1], jacc[, 2], jacc[, 3]), 2)


## Beispielkarten

ggplot() +
  geom_sf(data = Berlin) +
  geom_sf(data = hs_raster_to_sf(var = "E_E65U80", versions = c("original", "low", "high")) %>%
            arrange(desc(E_E65U80)),
          aes(color = E_E65U80, fill = E_E65U80), alpha = 0.3, lwd = 0.7) +
  theme_void()
#ggsave(filename = "Hotspots_E_E65U80.png", width = 2400, height = 1200, units = "px")

ggplot() +
  geom_sf(data = Berlin) +
  geom_sf(data = hs_raster_to_sf(var = "HK_EheJug", versions = c("original", "low", "high")) %>%
            arrange(desc(HK_EheJug)),
          aes(color = HK_EheJug, fill = HK_EheJug), alpha = 0.3, lwd = 0.7) +
  theme_void()
#ggsave(filename = "Hotspots_HK_EheJug.png", width = 2400, height = 1200, units = "px")



ggplot(dens_spop) +
  geom_raster(aes(lon, lat, fill = HK_EheJug)) +
  scale_fill_gradientn(colours = c("#FFFFFF", "#5c87c2", "#19224e")) +
  geom_sf(data = Berlin, fill = NA) +
  facet_wrap(~version) +
  theme_void()



# ----- Simulation Study -----
## adapted from: ?Kernelheaping::dbivr


# Mixed Normal Distributions

mu1 <- c(0, 0)
mu2 <- c(5, 3)
mu3 <- c(-4, 1)
SigmaB1 <- matrix(c(2, 0, 0, 2), 2, 2)
SigmaB2 <- matrix(c(1, 0, 0, 1), 2, 2)
SigmaB3 <- matrix(c(1, 0, 0, 3), 2, 2)
SigmaC1 <- matrix(c(4, 3, 3, 4), 2, 2)
SigmaC2 <- matrix(c(3, 0.5, 0.5, 1), 2, 2)
SigmaC3 <- matrix(c(5, 4, 4, 6), 2, 2)

mus <- rbind(mu1, mu2, mu3)
SigmasB <- rbind(SigmaB1, SigmaB2, SigmaB3)
SigmasC <- rbind(SigmaC1, SigmaC2, SigmaC3)
props <- c(1/3, 1/3, 1/3)

nsamp <- 500
set.seed(12345)
xtrueA <- rmvnorm.mixt(n = nsamp, mus = mu1, Sigmas = SigmaB2, props = 1)
xtrueB <- rmvnorm.mixt(n = nsamp, mus = mus, Sigmas = SigmasB, props = props)
xtrueC <- rmvnorm.mixt(n = nsamp, mus = mus, Sigmas = SigmasC, props = props)


# rounding

roundvalue <- c(0.75, 1.5, 2.25)

xroundA1 <- plyr::round_any(xtrueA, roundvalue[1])
xroundA2 <- plyr::round_any(xtrueA, roundvalue[2])
xroundA3 <- plyr::round_any(xtrueA, roundvalue[3])

xroundB1 <- plyr::round_any(xtrueB, roundvalue[1])
xroundB2 <- plyr::round_any(xtrueB, roundvalue[2])
xroundB3 <- plyr::round_any(xtrueB, roundvalue[3])

xroundC1 <- plyr::round_any(xtrueC, roundvalue[1])
xroundC2 <- plyr::round_any(xtrueC, roundvalue[2])
xroundC3 <- plyr::round_any(xtrueC, roundvalue[3])


pts <- rbind(xtrueA, xroundA1, xroundA2, xroundA3,
             xtrueB, xroundB1, xroundB2, xroundB3,
             xtrueC, xroundC1, xroundC2, xroundC3) %>% as.data.frame()
colnames(pts) <- c("x", "y")

pts$scenario <- as.factor(rep(c("A", "B", "C"), each = 4*nsamp))
pts$rounding <- as.factor(rep(rep(c("orig.", "r=0.75", "r=1.5", "r=2.25"), 
                                  each = nsamp), 3))
set.seed(59756160)
pts$rk <- runif(nrow(pts), 0, 1)

pts <- pts %>% group_by(x, y, scenario, rounding) %>%
  summarise(count = n(), ck = sum(rk))
pts$ck <- pts$ck - floor(pts$ck)

# perturb with a very small ptable
ptab_min <- create_cnt_ptable(D = 3, V = 2.5, js = 2, pstay = 0.5)
plot(ptab_min)

pts$count_ckm <- pts$count
pts$count_ckm[pts$rounding != "orig."] <- pts$count_ckm[pts$rounding != "orig."] +
  apply(pts[pts$rounding != "orig.", c("count", "ck")], 1, 
        get_perturbation_value, pt = ptab_min@pTable)


pts2 <- data.frame(x = c(pts[pts$scenario == "C" & pts$rounding == "orig.", ]$x,
                         rep(pts[pts$scenario == "C" & pts$rounding == "r=2.25", ]$x, 2)),
                   y = c(pts[pts$scenario == "C" & pts$rounding == "orig.", ]$y,
                         rep(pts[pts$scenario == "C" & pts$rounding == "r=2.25", ]$y, 2)),
                   count = c(pts[pts$scenario == "C" & pts$rounding == "orig.", ]$count,
                             pts[pts$scenario == "C" & pts$rounding == "r=2.25", ]$count,
                             pts[pts$scenario == "C" & pts$rounding == "r=2.25", ]$count_ckm),
                   type = rep(c("coordinates", "rounded", "rounded + CKM"), 
                            c(nrow(pts[pts$scenario == "C" & pts$rounding == "orig.", ]),
                              rep(nrow(pts[pts$scenario == "C" & pts$rounding == "r=2.25", ]), 2))))
pts2$type <- as.factor(pts2$type)

ggplot(pts2, aes(x, y)) +
  geom_point(aes(size = count), color = "#5c87c2", alpha = 0.4, show.legend = FALSE) +
  geom_text(aes(label = count)) +
  facet_wrap(~type)
#ggsave("concept.png", height = 800, width = 2400, units = "px")


# plot difference, estimate from CKM, do full sim., ...

# (...)





ggplot(pts) +
  geom_point(aes(x, y, size = count), color = "#5c87c2", show.legend = FALSE) +
  facet_grid(rounding ~ scenario) +
  theme_minimal()


#ggplot(pts[pts$rounding != "orig.", ], aes(count)) +
#  geom_histogram() +
#  facet_grid(rounding~scenario)




# estimating from rounded (heaped) data

estA1 <- dbivr(xroundA1, roundvalue = roundvalue[1], burnin = 5, samples = 10)
estA2 <- dbivr(xroundA2, roundvalue = roundvalue[2], burnin = 5, samples = 10)
estA3 <- dbivr(xroundA3, roundvalue = roundvalue[3], burnin = 5, samples = 10)

estB1 <- dbivr(xroundB1, roundvalue = roundvalue[1], burnin = 5, samples = 10)
estB2 <- dbivr(xroundB2, roundvalue = roundvalue[2], burnin = 5, samples = 10)
estB3 <- dbivr(xroundB3, roundvalue = roundvalue[3], burnin = 5, samples = 10)

estC1 <- dbivr(xroundC1, roundvalue = roundvalue[1], burnin = 5, samples = 10)
estC2 <- dbivr(xroundC2, roundvalue = roundvalue[2], burnin = 5, samples = 10)
estC3 <- dbivr(xroundC3, roundvalue = roundvalue[3], burnin = 5, samples = 10)

# estimating from true data

densA <- dmvnorm.mixt(x = expand.grid(estA1$Mestimates$eval.points[[1]],
                                      estA1$Mestimates$eval.points[[2]]),
                      mus = mu1, Sigmas = SigmaB2, props = 1)
densB <- dmvnorm.mixt(x = expand.grid(estB1$Mestimates$eval.points[[1]],
                                      estB1$Mestimates$eval.points[[2]]),
                      mus = mus, Sigmas = SigmasB, props = props)
densC <- dmvnorm.mixt(x = expand.grid(estC1$Mestimates$eval.points[[1]],
                                      estC1$Mestimates$eval.points[[2]]),
                      mus = mus, Sigmas = SigmasC, props = props)

d <- rbind(expand.grid(x = estA1$Mestimates$eval.points[[1]],
                       y = estA1$Mestimates$eval.points[[2]]),
           expand.grid(x = estA1$Mestimates$eval.points[[1]],
                       y = estA1$Mestimates$eval.points[[2]]),
           expand.grid(x = estA2$Mestimates$eval.points[[1]],
                       y = estA2$Mestimates$eval.points[[2]]),
           expand.grid(x = estA3$Mestimates$eval.points[[1]],
                       y = estA3$Mestimates$eval.points[[2]]),
           expand.grid(x = estB1$Mestimates$eval.points[[1]],
                       y = estB1$Mestimates$eval.points[[2]]),
           expand.grid(x = estB1$Mestimates$eval.points[[1]],
                       y = estB1$Mestimates$eval.points[[2]]),
           expand.grid(x = estB2$Mestimates$eval.points[[1]],
                       y = estB2$Mestimates$eval.points[[2]]),
           expand.grid(x = estB3$Mestimates$eval.points[[1]],
                       y = estB3$Mestimates$eval.points[[2]]),
           expand.grid(x = estC1$Mestimates$eval.points[[1]],
                       y = estC1$Mestimates$eval.points[[2]]),
           expand.grid(x = estC1$Mestimates$eval.points[[1]],
                       y = estC1$Mestimates$eval.points[[2]]),
           expand.grid(x = estC2$Mestimates$eval.points[[1]],
                       y = estC2$Mestimates$eval.points[[2]]),
           expand.grid(x = estC3$Mestimates$eval.points[[1]],
                       y = estC3$Mestimates$eval.points[[2]]))

d$density <- c(densA, 
               as.vector(estA1$Mestimates$estimate),
               as.vector(estA2$Mestimates$estimate),
               as.vector(estA3$Mestimates$estimate),
               densB, 
               as.vector(estB1$Mestimates$estimate),
               as.vector(estB2$Mestimates$estimate),
               as.vector(estB3$Mestimates$estimate),
               densC,
               as.vector(estC1$Mestimates$estimate),
               as.vector(estC2$Mestimates$estimate),
               as.vector(estC3$Mestimates$estimate))
d$density[d$density <= 0] <- 0

d$scenario <- as.factor(rep(c("A", "B", "C"), each = 4*4e+4))
d$rounding <- as.factor(rep(rep(c("orig.", "r=0.75", "r=1.5", "r=2.25"), each = 4e+4), 3))

ggplot(d[d$density > 0.01, ], aes(x = x, y = y, z = density)) +
  geom_contour_filled(color = "black", bins = 10, show.legend = FALSE) +
  facet_grid(rounding ~ scenario)


