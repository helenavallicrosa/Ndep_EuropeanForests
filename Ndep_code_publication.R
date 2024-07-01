#### Fig. 1 GAM's functional group ####
load("~/Desktop/WSL/Saturation/v4/data_v4.Rdata")

db <- subset(db, db$growth_pond < 0.083)
db <- subset(db, db$growth_pond > 0)
db$map_CRU <- db$map_CRU*10
db$total_5y <- db$total_5y/5
a <- plyr::count(db$plot_tree)
db <- merge(db, a, by.x="plot_tree", by.y="x")
db <- subset(db, db$freq > 1)
db_con <- subset(db, db$grp_tree_species == "conifers")
db_bro <- subset(db, db$grp_tree_species == "broadleaves")


db_con$plot_tree <- as.factor(db_con$plot_tree)
db_con$code_tree_species <- as.factor(db_con$code_tree_species)

library(mgcv)
# citation("mgcv")
library(gratia)
library(voxel)
# citation("voxel")
library(ggplot2)

db_con$plot <- as.factor(db_con$plot)

growth_con <- gam(growth_pond ~ s(total_5y,k=7, bs="ds") +
                    s(lat_decimal,k=5, bs="re"),
                  data=db_con, na.action=na.omit)

growth_bro <- gam(growth_pond ~ s(total_5y,k=6, bs="ds") +
                    s(lat_decimal,k=5, bs="re"),
                  data=db_bro, na.action=na.omit)

summary(growth_con)
summary(growth_bro)

con_plot <- plotGAM(growth_con, smooth.cov = "total_5y")+
  theme_classic()+
  xlab("N deposition (kgN/ha/yr)")+
  annotate(geom = 'text', x = 9, y = 0.004,
           label = paste("p. value < 0.0001 \n r2 = 0.034"), size=3) +
  xlim(c(0,50))+
  # ylim(c(0.003,0.022))+
  ylim(c(0.002,0.028))+
  ylab("Relative growth (cm/cm/yr)")+
  ggtitle("Conifers")
con_plot

hist(db_con$total_5y)

bro_plot <- plotGAM(growth_bro, smooth.cov = "total_5y")+
  theme_classic()+
  annotate(geom = 'text', x = 9, y = 0.004,
           label = paste("p. value < 0.0001 \n r2 = 0.019"), size=3) +
  xlab("N deposition (kgN/ha/yr)")+
  ylim(c(0.002,0.028))+
  xlim(c(0,50))+
  theme(axis.text.y = element_blank(),
        axis.title.y=element_blank())+
  ggtitle("Broadleaves")+
  ylab("Relative growth (cm/cm/yr)")
bro_plot

library(ggpubr)
tiff(filename="~/Desktop/WSL/Saturation/v4/con_bro_gam.tif", width=1000, height = 450, res= 150)
ggarrange(con_plot, bro_plot, ncol=2,widths = c(5.5,4.5))
dev.off()

#### Fig. 2 Functional group ####
load("~/Desktop/WSL/Saturation/v4/data_v4.Rdata")

db <- subset(db, db$growth_pond < 0.083)
db <- subset(db, db$growth_pond > 0)
db$map_CRU <- db$map_CRU*10
db$total_5y <- db$total_5y/5

a <- plyr::count(db$plot_tree)
db <- merge(db, a, by.x="plot_tree", by.y="x")
db <- subset(db, db$freq > 1)
db_con <- subset(db, db$grp_tree_species == "conifers")
db_bro <- subset(db, db$grp_tree_species == "broadleaves")

m_con <- lmer(growth_pond ~ mat_CRU+map_CRU+total_5y+survey_year+
                mat_CRU:total_5y+map_CRU:total_5y+
                (1|plot_tree), data=db_con)
summary(m_con)
hist(resid(m_con))

m_bro <- lmer(growth_pond ~ mat_CRU+map_CRU+total_5y+survey_year+
                mat_CRU:total_5y+map_CRU:total_5y+
                (1|plot_tree), data=db_bro)
summary(m_bro)

save(m_con, file="~/Desktop/WSL/Saturation/v4/m_con.Rdata")
save(m_bro, file="~/Desktop/WSL/Saturation/v4/m_bro.Rdata")

MuMIn::r.squaredGLMM(m_con)
pcon_temp <- plot_model(m_con, type = "pred", terms = c("total_5y", "mat_CRU [5,10,15]"),
                        colors=c("#ffb563","#f85e00", "#a41623"),se=T)+
  theme_classic() +
  ggtitle("Conifers")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0.2, 0, 0, 0.2),"cm"))+
  ylab("Relative growth (cm/cm/yr)")+
  ylim(c(0.006,0.021))+
  xlim(c(0,50))+
  annotate(geom = 'text', x = 22, y = 0.009,
           label = paste("p. value < 0.0001 \n r2 = 0.66"), size=3) +
  theme(legend.position = c(.63,.80),
        legend.box.background = element_rect(size=0.5),
        axis.text.x = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.title.x=element_blank())+
  guides(color=guide_legend(title="MAT ºC")) 
pcon_temp

pcon_map <- plot_model(m_con, type = "pred", terms = c("total_5y", "map_CRU [500,1000,1500]"),
                       colors=c("#a8dadc", "#457b9d", "#1d3557"), se=F)+
  theme_classic() +
  ggtitle("")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0),"inches"))+
  annotate(geom = 'text', x = 20, y = 0.009,
           label = paste("p. value < 0.0001 \n r2 = 0.66"), size=3) +
  ylim(c(0.006,0.021))+
  xlim(c(0,50))+
  theme(legend.position = c(.20,.80),
        legend.box.background = element_rect(size=0.5), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  guides(color=guide_legend(title="MAP (mm/y)")) 
pcon_map

MuMIn::r.squaredGLMM(m_bro)
pbro_temp <- plot_model(m_bro, type = "pred", terms = c("total_5y", "mat_CRU [5,10,15]"),
                        colors=c("#ffb563","#f85e00", "#a41623"),se=F)+
  theme_classic() +
  ggtitle("Broadleaved")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0.2, 0, 0.2, 0.2),"cm"))+
  annotate(geom = 'text', x = 10, y = 0.009,
           label = paste("p. value < 0.0001 \n r2=0.58"), size=3) +
  xlab("N deposition (kgN/ha/yr)") +
  ylab("Relative growth (cm/cm/yr)")+
  ylim(c(0.006,0.021))+
  theme(legend.position = "none")+
  guides(color=guide_legend(title="MAT ºC")) 
pbro_temp

pbro_map <- plot_model(m_bro_coex, type = "pred", terms = c("total_5y", "map_CRU [500,1000,1500]"),
                       colors=c("#a8dadc", "#457b9d", "#1d3557"))+
  theme_classic() +
  ggtitle("")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0, 0, 0.2, 0),"inches"))+
  annotate(geom = 'text', x = 34, y = 0.009,
           label = paste("p. value < 0.0001 \n r2=0.58"), size=3) +
  xlab("N deposition (kgN/ha/yr)") +
  ylim(c(0.007,0.021))+
  xlim(c(0,50))+
  theme(legend.position = "none",axis.text.y = element_text("none"),
        axis.title.y=element_blank())+
  guides(color=guide_legend(title="MAP (cm/yr)")) 
pbro_map

library(ggpubr)
tiff(filename="~/Desktop/WSL/Saturation/v4/Lmer_Con_bro_mat_map.tif", width=1100, height = 1100, res= 140)
ggarrange(pcon_temp, pcon_map, pbro_temp, pbro_map, ncol=2, nrow=2, align="h", widths = c(5.5,4.5),
          labels= c("a)","b)", "c)", "d)"))
dev.off()

#### Fig. 3 4 maps sensitivity ####
library(ncdf4)
library(raster)

load("~/Desktop/WSL/Saturation/v4/m_con.Rdata")
load("~/Desktop/WSL/Saturation/v4/m_bro.Rdata")

con_map <- stack("~/Desktop/WSL/Saturation/v4/vars_coniferous.gri")
con_map$MAP <- con_map$MAP*100
broad_map <- stack("~/Desktop/WSL/Saturation/v4/vars_broadleaved.gri")
broad_map$MAP <- broad_map$MAP*100

estim_con <- summary(m_con)
estim_con <- estim_con$coefficients
estim_bro <- summary(m_bro)
estim_bro <- estim_bro$coefficients

con_density <- raster("~/Desktop/WSL/Saturation/v4/con_density.tif")
bro_density <- raster("~/Desktop/WSL/Saturation/v4/bro_density.tif")

con_percent <- con_density/(con_density+bro_density)
bro_percent <- bro_density/(con_density+bro_density)

con_Ndep_sens_mat <- (estim_con[4,1]+(con_map$MAT*estim_con[6,1])+(mean(na.exclude(db_con$map_CRU))*estim_con[7,1]))
bro_Ndep_sens_mat <- (estim_bro[4,1]+(broad_map$MAT*estim_bro[6,1])+(mean(na.exclude(db_bro$map_CRU))*estim_bro[7,1]))

con_Ndep_sens_map <- (estim_con[4,1]+(con_map$MAP*estim_con[7,1])+(mean(na.exclude(db_con$mat_CRU))*estim_con[6,1]))
bro_Ndep_sens_map <- (estim_bro[4,1]+(broad_map$MAP*estim_bro[7,1])+(mean(na.exclude(db_bro$mat_CRU))*estim_bro[6,1]))

con_Ndep_sens_mat_df <- as.data.frame(con_Ndep_sens_mat, xy=T)
bro_Ndep_sens_mat_df <- as.data.frame(bro_Ndep_sens_mat, xy=T)
con_Ndep_sens_map_df <- as.data.frame(con_Ndep_sens_map, xy=T)
bro_Ndep_sens_map_df <- as.data.frame(bro_Ndep_sens_map, xy=T)

library(sf)
global_map_shp <- read_sf('~/Desktop/Batch maps/Global Maps/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')
options(scipen=999)

hist(con_Ndep_sens_mat_df$MAT)
hist(bro_Ndep_sens_mat_df$MAT)

con_Ndep_sens_mat_df$MAT[con_Ndep_sens_mat_df$MAT < -0.0003] <- -0.0003
con_Ndep_sens_mat_df$MAT[con_Ndep_sens_mat_df$MAT > 0.0003] <- 0.0003
bro_Ndep_sens_mat_df$MAT[bro_Ndep_sens_mat_df$MAT < -0.0003] <- -0.0003
bro_Ndep_sens_mat_df$MAT[bro_Ndep_sens_mat_df$MAT > 0.0003] <- 0.0003

library(ggplot2)
con_Ndep_sens_mat_2022 <- ggplot(con_Ndep_sens_mat_df) + 
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  geom_raster(aes(x, y, fill=MAT)) +
  scale_fill_gradient2(low = "#c1121f", high = "#669bbc", midpoint = 0, na.value = NA) + 
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1), 
        legend.title = element_text("Relative growth")) +
  ggtitle("Conifers: sensitivity to Ndep depending on MAT")+
  xlab("") +
  ylab("")
con_Ndep_sens_mat_2022

bro_Ndep_sens_mat_2022 <- ggplot(bro_Ndep_sens_mat_df) + 
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  geom_raster(aes(x, y, fill=MAT)) +
  scale_fill_gradient2(low = "#c1121f", high = "#669bbc", midpoint = 0, na.value = NA) + 
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1), 
        legend.title = element_text("Relative growth")) +
  ggtitle("Broadleaves: sensitivity to Ndep depending on MAT")+
  xlab("") +
  ylab("")
bro_Ndep_sens_mat_2022

hist(con_Ndep_sens_map_df$MAP)
hist(bro_Ndep_sens_map_df$MAP)

con_Ndep_sens_map_df$MAP[con_Ndep_sens_map_df$MAP < -0.0001] <- -0.0001
con_Ndep_sens_map_df$MAP[con_Ndep_sens_map_df$MAP > 0.0001] <- 0.0001
bro_Ndep_sens_map_df$MAP[bro_Ndep_sens_map_df$MAP < -0.0001] <- -0.0001
bro_Ndep_sens_map_df$MAP[bro_Ndep_sens_map_df$MAP > 0.0001] <- 0.0001

con_Ndep_sens_map_2022 <- ggplot(con_Ndep_sens_map_df) + 
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  geom_raster(aes(x, y, fill=MAP)) +
  scale_fill_gradient2(low = "#c1121f", high = "#669bbc", midpoint = 0, na.value = NA) + 
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1), 
        legend.title = element_text("Relative growth")) +
  ggtitle("Conifers: sensitivity to Ndep depending on MAP")+
  xlab("") +
  ylab("")
con_Ndep_sens_map_2022

bro_Ndep_sens_map_2022 <- ggplot(bro_Ndep_sens_map_df) + 
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  geom_raster(aes(x, y, fill=MAP)) +
  scale_fill_gradient2(low = "#c1121f", high = "#669bbc", midpoint = 0, na.value = NA) + 
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "bottom", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1), 
        legend.title = element_text("Relative growth")) +
  ggtitle("Broadleaves: sensitivity to Ndep depending on MAP")+
  xlab("") +
  ylab("")
bro_Ndep_sens_mat_2022

library(ggpubr)
tiff(filename="~/Desktop/WSL/Saturation/v4/Intermedis/4maps_int2.tif", width=1200, height = 1200, res= 120)
ggarrange(con_Ndep_sens_mat_2022, con_Ndep_sens_map_2022, bro_Ndep_sens_mat_2022,bro_Ndep_sens_map_2022,
          ncol=2, nrow=2, common.legend = F, legend="right", labels= c("a)","b)", "c)", "d)"))
dev.off()

#### Fig. 5 total Ndep sensitivity ####

total_temp_sens <- con_Ndep_sens_mat*con_percent + bro_Ndep_sens_mat*bro_percent
total_map_sens <- con_Ndep_sens_map*con_percent + bro_Ndep_sens_map*bro_percent

total <- total_map_sens+total_temp_sens
total_df <- as.data.frame(total, xy=T)
save(total_df, file="~/Desktop/WSL/Saturation/v4/Intermedis/total_Ndep_sensitivity.Rdata")
hist(total_df$layer)
total_df$layer[total_df$layer < -0.0002] <- -0.0002
total_df$layer[total_df$layer > 0.0002] <- 0.0002

plot(total)

total_plot <- ggplot(total_df) + 
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  # scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5, na.value = NA) +
  geom_raster(aes(x, y, fill=layer)) +
  scale_fill_gradient2(low = "#c1121f", high = "#669bbc", midpoint = 0, na.value = NA) + 
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "right", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1), 
        legend.title = element_text("Relative growth")) +
  # ggtitle("Areas above N saturation level in 2022")+
  ggtitle("N deposition effect on growth")+
  xlab("") +
  ylab("")
total_plot

tiff(filename="~/Desktop/WSL/Saturation/v4/total_map.tif", width=900, height = 800, res= 120)
total_plot
dev.off()

#### Fig 4. N-limitation ####
load("~/Desktop/WSL/Saturation/v4/data_v4.Rdata")

db <- subset(db, db$growth_pond < 0.083)
db <- subset(db, db$growth_pond > 0)
db$map_CRU <- db$map_CRU*10
db$total_5y <- db$total_5y/5
a <- plyr::count(db$plot_tree)
db <- merge(db, a, by.x="plot_tree", by.y="x")
db <- subset(db, db$freq > 1)
db_con <- subset(db, db$grp_tree_species == "conifers")
db_bro <- subset(db, db$grp_tree_species == "broadleaves")

db_con$group_mat[db_con$mat_CRU <= 5] <- "<5"
db_con$group_mat[db_con$mat_CRU > 5 & db_con$mat_CRU <= 10 ] <- "5-10"
db_con$group_mat[db_con$mat_CRU > 10 & db_con$mat_CRU <= 15 ] <- "10-15"
db_con$group_mat[db_con$mat_CRU > 15 ] <- ">15"
db_con <- subset(db_con, !is.na(db_con$group_mat))

db_bro$group_mat[db_bro$mat_CRU <= 5] <- "<5"
db_bro$group_mat[db_bro$mat_CRU > 5 & db_bro$mat_CRU <= 10 ] <- "5-10"
db_bro$group_mat[db_bro$mat_CRU > 10 & db_bro$mat_CRU <= 15 ] <- "10-15"
db_bro$group_mat[db_bro$mat_CRU > 15 ] <- ">15"
db_bro <- subset(db_bro, !is.na(db_bro$group_mat))

db_con$group_mat <- factor(db_con$group_mat, levels = c("<5", "5-10", "10-15", ">15"))
db_bro$group_mat <- factor(db_bro$group_mat, levels = c("<5", "5-10", "10-15", ">15"))

plyr::count(db_con$group_mat)
plyr::count(db_bro$group_mat)

Nlim_con <- ggplot(data=db_con, aes(x=group_mat, y=NP_limitation))+
  geom_boxplot(data=db_con, aes(x=group_mat, y=NP_limitation),fill="#fb8500")+
  theme_classic() +
  # stat_compare_means(comparisons=my_comparisons)+
  xlab("Mean annual temperature")+
  ylab("N limitation")+
  ggtitle("Conifers")
Nlim_con

Nlim_bro <- ggplot(data=db_bro, aes(x=group_mat, y=NP_limitation))+
  geom_boxplot(fill="#219ebc")+
  theme_classic() +
  # stat_compare_means(comparisons=my_comparisons, method = "t.test")+
  xlab("Mean annual temperature")+
  ylab("N limitation")+
  theme(axis.title.y=element_blank())+
  ggtitle("Broadleaves")

library(ggpubr)
tiff(filename="~/Desktop/WSL/Saturation/v4/Nlimit_con_bor.tif", width=1200, height = 500, res= 140)
ggarrange(Nlim_con, Nlim_bro)
dev.off()

#### SUPPLEMENTARY ####
#### Fig. S1 plots ####
library(sf)
load("~/Desktop/WSL/Saturation/v4/data_v4.Rdata")
global_map_shp <- read_sf('~/Desktop/Batch maps/Global Maps/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')

library(ggplot2)
map_plot_growth <- ggplot(db) + 
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray80") +
  # geom_point(aes(x=lon_decimal, y=lat_decimal, color=total_5y), shape=20, size=1, alpha=0.6)+
  geom_point(aes(x=lon_decimal, y=lat_decimal,
                 colour=factor(grp_tree_species)), size=0.8, alpha=0.5)+
  theme_minimal() +
  labs(colour="Leaf form") +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1)) +
  # labs(fill="Vegetation type", size="Sample size") +
  ylim(c(34,71)) +
  xlim(c(-30,35)) +
  ggtitle("ICP forest plots")+
  xlab("")+
  ylab("")
map_plot_growth

tiff(filename="~/Desktop/WSL/Saturation/v4/plot_map.tif", width=1000, height = 900, res= 140)
map_plot_growth
dev.off()

#### Fig. S2 distribution ####
con_density <- raster("~/Desktop/WSL/Saturation/v4/con_density.tif")
bro_density <- raster("~/Desktop/WSL/Saturation/v4/bro_density.tif")

con_percent <- con_density/(con_density+bro_density)
bro_percent <- bro_density/(con_density+bro_density)

con_density <- mask(con_density, con_map$MAT)
con_density_df <- as.data.frame(con_density, xy=T)
bro_density <- mask(bro_density, con_map$MAT)
bro_density_df <- as.data.frame(bro_density, xy=T)
con_percent <- mask(con_percent, con_map$MAT)
con_percent_df <- as.data.frame(con_percent, xy=T)
bro_percent <- mask(bro_percent, con_map$MAT)
bro_percent_df <- as.data.frame(bro_percent, xy=T)

hist(con_density_df$con_density)
hist(bro_density_df$bro_density)
con_density_df$con_density[con_density_df$con_density > 1500] <- 1500
bro_density_df$bro_density[bro_density_df$bro_density > 1500] <- 1500

con_density_plot <- ggplot(con_density_df) +
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  # scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5, na.value = NA) +
  geom_raster(aes(x, y, fill=con_density)) +
  scale_fill_gradient2(low = "#f2e8cf", high = "#386641", na.value = NA) +
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "right", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        legend.title = element_text("Relative growth")) +
  # ggtitle("Areas above N saturation level in 2022")+
  ggtitle("Conifer basal area")+
  xlab("") +
  ylab("")
con_density_plot

con_percent_plot <- ggplot(con_percent_df) +
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  # scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5, na.value = NA) +
  geom_raster(aes(x, y, fill=layer)) +
  scale_fill_gradient2(low = "#f2e8cf", high = "#386641", na.value = NA) +
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "right", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        legend.title = element_text("Relative growth")) +
  # ggtitle("Areas above N saturation level in 2022")+
  ggtitle("Conifer basal area")+
  xlab("") +
  ylab("")
con_percent_plot

bro_density_plot <- ggplot(bro_density_df) +
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  # scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5, na.value = NA) +
  geom_raster(aes(x, y, fill=bro_density)) +
  scale_fill_gradient2(low = "#f2e8cf", high = "#386641", na.value = NA) +
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "right", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        legend.title = element_text("Relative growth")) +
  # ggtitle("Areas above N saturation level in 2022")+
  ggtitle("Broadleaf basal area")+
  xlab("") +
  ylab("")
bro_density_plot

bro_percent_plot <- ggplot(bro_percent_df) +
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  # scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5, na.value = NA) +
  geom_raster(aes(x, y, fill=layer)) +
  scale_fill_gradient2(low = "#f2e8cf", high = "#386641", na.value = NA) +
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "right", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        legend.title = element_text("Relative growth")) +
  # ggtitle("Areas above N saturation level in 2022")+
  ggtitle("Broadleaf basal area")+
  xlab("") +
  ylab("")
bro_percent_plot

library(ggpubr)
tiff(filename="~/Desktop/WSL/Saturation/v4/SM/density_maps.tif", width=1600, height = 1600, res= 120)
ggarrange(con_density_plot, bro_density_plot, con_percent_plot, bro_percent_plot, ncol=2, nrow=2, common.legend = T, legend="right")
dev.off()

#### Fig. S3 mixed plots ####
load("~/Desktop/WSL/Saturation/v4/data_v4.Rdata")
db <- subset(db, db$growth_pond < 0.083)
db <- subset(db, db$growth_pond > 0)
db$map_CRU <- db$map_CRU*10
db$total_5y <- db$total_5y/5

a <- plyr::count(db$plot_tree)
db <- merge(db, a, by.x="plot_tree", by.y="x")
db <- subset(db, db$freq > 1)

####
coexisting <- db[,c(22,23)]
coex <- unique(coexisting)
a <- plyr::count(coex$plot)
coexisting_plots <- subset(a, a$freq > 1)
coexisting_plots$coexisting <- "YES"
coexisting_plots <- coexisting_plots[,-2]
standalone <- subset(a, a$freq < 2)
standalone$coexisting <- "NO"
standalone <- standalone[,-2]

coex <- rbind(coexisting_plots,standalone)
db <- merge(db,coex, by.x="plot", by.y="x")

db$bro_con_coex <- paste(db$grp_tree_species, db$coex, sep="_")
coex <- subset(db, db$coexisting == "YES")

coex_con <- subset(coex, coex$grp_tree_species == "conifers")
coex_bro <- subset(coex, coex$grp_tree_species == "broadleaves")

m_con_coex <- lmer(growth_pond ~ mat_CRU+map_CRU+total_5y+survey_year+
                     mat_CRU:total_5y+map_CRU:total_5y+
                     (1|plot_tree), data=coex_con)
summary(m_con_coex)

m_bro_coex <- lmer(growth_pond ~ mat_CRU+map_CRU+total_5y+survey_year+
                     mat_CRU:total_5y+map_CRU:total_5y+
                     (1|plot_tree), data=coex_bro)
summary(m_bro_coex)
MuMIn::r.squaredGLMM(m_con_coex)
pcon_temp <- plot_model(m_con_coex, type = "pred", terms = c("total_5y", "mat_CRU [5,10,15]"),
                        colors=c("#ffb563","#f85e00", "#a41623"),se=T)+
  theme_classic() +
  ggtitle("Conifers")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0.2, 0, 0, 0.2),"cm"))+
  ylab("Relative growth (cm/cm/yr)")+
  ylim(c(0.00,0.025))+
  xlim(c(0,50))+
  annotate(geom = 'text', x = 22, y = 0.009,
           label = paste("p. value < 0.0001 \n r2 = 0.69"), size=3) +
  theme(legend.position = c(.30,.80),
        legend.box.background = element_rect(size=0.5),
        axis.text.x = element_blank(),
        legend.key = element_rect(fill = "white"),
        axis.title.x=element_blank())+
  guides(color=guide_legend(title="MAT ºC")) 
pcon_temp

pcon_map <- plot_model(m_con_coex, type = "pred", terms = c("total_5y", "map_CRU [500,1000,1500]"),
                       colors=c("#a8dadc", "#457b9d", "#1d3557"), se=F)+
  theme_classic() +
  ggtitle("")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0),"inches"))+
  annotate(geom = 'text', x = 20, y = 0.009,
           label = paste("p. value < 0.0001 \n r2 = 0.69"), size=3) +
  ylim(c(0.00,0.025))+
  # xlim(c(0,50))+
  theme(legend.position = c(.20,.80),
        legend.box.background = element_rect(size=0.5), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  guides(color=guide_legend(title="MAP (mm/y)")) 
pcon_map

MuMIn::r.squaredGLMM(m_bro_coex)
summary(m_bro_coex)
pbro_temp <- plot_model(m_bro_coex, type = "pred", terms = c("total_5y", "mat_CRU [5,10,15]"),
                        colors=c("#ffb563","#f85e00", "#a41623"),se=F)+
  theme_classic() +
  ggtitle("Broadleaved")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0.2, 0, 0.2, 0.2),"cm"))+
  annotate(geom = 'text', x = 10, y = 0.009,
           label = paste("p. value < 0.0001 \n r2=0.67"), size=3) +
  xlab("N deposition (kgN/ha/yr)") +
  ylab("Relative growth (cm/cm/yr)")+
  ylim(c(0.00,0.025))+
  # xlim(c(0,230))+
  theme(legend.position = "none")+
  guides(color=guide_legend(title="MAT ºC")) 
pbro_temp

pbro_map <- plot_model(m_bro_coex, type = "pred", terms = c("total_5y", "map_CRU [500,1000,1500]"),
                       colors=c("#a8dadc", "#457b9d", "#1d3557"))+
  theme_classic() +
  ggtitle("")+
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.margin = unit(c(0, 0, 0.2, 0),"inches"))+
  annotate(geom = 'text', x = 34, y = 0.009,
           label = paste("p. value < 0.0001 \n r2=0.67"), size=3) +
  xlab("N deposition (kgN/ha/yr)") +
  ylim(c(0.00,0.025))+
  xlim(c(0,50))+
  theme(legend.position = "none",axis.text.y = element_text("none"),
        axis.title.y=element_blank())+
  guides(color=guide_legend(title="MAP (cm/yr)")) 
pbro_map

library(ggpubr)
tiff(filename="~/Desktop/WSL/Saturation/v4/SM/Lmer_Con_bro_mat_map_coexisting.tif", width=1100, height = 1100, res= 140)
ggarrange(pcon_temp, pcon_map, pbro_temp, pbro_map, ncol=2, nrow=2, align="h", widths = c(5.5,4.5),
          labels= c("a)","b)", "c)", "d)"))
dev.off()

#### Fig. S4 Correlation ####
load("~/Desktop/WSL/Saturation/v4/data_v4.Rdata")
library(corrplot)
colnames(db)
db_con <- subset(db, db$grp_tree_species == "conifers")
db_bro <- subset(db, db$grp_tree_species == "broadleaves")

db_pca <- na.exclude(db[,c(13,14,18,19,24,25)])

M1 <- cor(db_pca)
colnames(M1) <- c("MAT", "MAP", "Ndep", "Aridity", "Soil N", "N lim.")
rownames(M1) <- c("MAT", "MAP", "Ndep", "Aridity", "Soil N", "N lim.")

cor1 <- corrplot.mixed(M1, upper = 'number', lower = "square", tl.col = "black")

tiff(filename="~/Desktop/WSL/Saturation/v4/SM/correlation.tif", width=800, height = 750, res= 140)
cor1 <- corrplot.mixed(M1, upper = 'number', lower = "square", tl.col = "black")
dev.off()

#### Fig. S5 data distribution ####
load("~/Desktop/WSL/Saturation/v4/data_v4.Rdata")

mat <- ggplot(db, aes(x=mat_CRU, y=total_5y) ) +
  geom_hex(bins = 30) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()+
  xlab("MAT")+
  ylab("N deposition")
mat

map <- ggplot(db, aes(x=map_CRU, y=total_5y) ) +
  geom_hex(bins = 30) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()+
  xlab("MAP")+
  ylab("N deposition")
map

library(ggpubr)  
tiff(filename="~/Desktop/WSL/Saturation/v4/data_distribution.tif", width=800, height = 400, res= 120)
ggarrange(mat, map, ncol=2, align="h", common.legend = T, legend = "right")
dev.off()  

#### Fig. S6 area above saturation ####
library(rgdal)
library(ncdf4)
library(raster)
library(sf)
library(sp)
library(ggplot2)

EMEP_2022 <- function(dir){
  #create outputs
  EMEP_2022 <-vector("list", 5)
  for (i in 1:length(dir)) {
    emep <- nc_open(dir[i])
    lon <- ncvar_get(emep, "lon")
    lat <- t(ncvar_get(emep, "lat", verbose = F))
    fillvalue <- emep[["var"]][["DDEP_OXN_m2Grid"]][["missval"]]
    dry.ox <- ncvar_get(emep, "DDEP_OXN_m2Grid")
    dry.ox[dry.ox == fillvalue] <- NA
    wet.ox <- ncvar_get(emep, "WDEP_OXN")
    wet.ox[wet.ox == fillvalue] <- NA
    dry.red <- ncvar_get(emep, "DDEP_RDN_m2Grid")
    dry.red[dry.red == fillvalue] <- NA
    wet.red <- ncvar_get(emep, "WDEP_RDN")
    wet.red[wet.red == fillvalue] <- NA
    nc_close(emep) 
    
    #map obtention
    dry.ox_map <- flip(raster(t(dry.ox), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110")))
    wet.ox_map <- flip(raster(t(wet.ox), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110")))
    ox_map <- dry.ox_map+wet.ox_map
    dry.red_map <- flip(raster(t(dry.red), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110")))
    wet.red_map <- flip(raster(t(wet.red), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110")))
    red_map <- dry.red_map+wet.red_map
    total_map <- ox_map+red_map
    
    #map exportation prep.
    EMEP_2022[i] <- total_map
  }
  return(EMEP_2022)
}

mapFiles <- dir("~/Desktop/Batch maps/Global maps/N dep/EMEP2022/", full.names=TRUE)
EMEP_sat_2022 <- EMEP_2022(dir=mapFiles)
EMEP_2022_stack <- stack(EMEP_sat_2022)
Nsat_2022 <- calc(EMEP_2022_stack, sum)

Nsat_2022 <- crop(Nsat_2022, extent(-29.95, 42, 30.05, 81.95))
con_map <- stack("~/Desktop/WSL/Saturation/v4/vars_coniferous.gri")
plot(con_map$Ndep)
Nsat_2022 <- resample(Nsat_2022, con_map$MAT)
Nsat_2022 <- mask(Nsat_2022, con_map$MAT)


Nsat_2022_yr <- Nsat_2022/500
Nsat_2022_yr$layer[Nsat_2022_yr$layer < 30] <- NA
Nsat_2022_yr$layer[Nsat_2022_yr$layer > 30] <- 1
hist(getValues(Nsat_2022_yr))

N_30more <- as.data.frame(Nsat_2022_yr, xy=T)
global_map_shp <- read_sf('~/Desktop/Batch maps/Global Maps/ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp')

N_30more_map <- ggplot(N_30more) +
  geom_sf(data = global_map_shp$geometry, size=0.3, colour = "black", fill = "gray90") +
  # scale_fill_gradient2(low = "white", high = "red", midpoint = 0.5, na.value = NA) +
  geom_raster(aes(x, y, fill=layer)) +
  scale_fill_gradient(low = "red", high = "red", na.value = NA) +
  theme_minimal()+
  ylim(c(35,71)) +
  xlim(c(-25,42))+
  coord_sf(datum = st_crs("+proj=stere +lat_0=90 +lon_0=-32 +lat_ts=60 +x_0=8 +y_0=110"))+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1),
        legend.title = element_text("Relative growth")) +
  ggtitle("Areas above 30kg/ha/yr")+
  xlab("") +
  ylab("")
N_30more_map

tiff(filename="~/Desktop/WSL/Saturation/v4/SM/saturated.tif", width=800, height = 700, res= 140)
N_30more_map
dev.off()