#Brinkerhoff, et al. (2019) AMHG Analysis
#Written By: Craig Brinkerhoff
#Description: Runs all analysis and builds all figures necessary to repeat this study.
#Directions: requires this script, as well as the field-measurement dataset, amended
    #NHD ID join table, and the calibrated grain size dataset to all be stored in the same working directory


#set your working directory to your current working directory, and then run the script
setwd("C:\\Users\\cbrinkerhoff\\Box Sync\\Ongoing Projects\\hydraulic_geometry_project\\working\\final\\")


#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed and loaded
packages = c("TOSTER", "dataRetrieval", "stringr", "tidyverse", "dplyr", "tidyr", "stats", "ggplot2", "gridExtra", "cowplot", "RColorBrewer", "Metrics", 'lsr', 'grid', 'hydroGOF')
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

#set copy
med_NHD <- read.csv('C:NHD_join_table.csv')

#get Barber data (already filtered to 6+ stations and 20+ measurements)------------------------------------
all_rivers <- read.csv("field_measurements.csv")

#add reach slopes from NHD
all_rivers <- left_join(all_rivers, select(med_NHD, SLOPE, SOURCE_FEA), by = c("site_no" = "SOURCE_FEA"))

#intial cleaning of field measurements----------------------------------------------------------
all_rivers <- filter(all_rivers, is.finite(chan_width) ==1) %>%
  filter(is.finite(chan_velocity)==1) %>%
  filter(is.finite(chan_discharge)==1) %>%
  filter(is.finite(chan_area)==1) %>%
  filter(chan_width > 0) %>%
  filter(chan_velocity > 0) %>%
  filter(chan_discharge > 0) %>%
  filter(chan_area > 0) %>%
  filter(measured_rating_diff != 'Poor') %>%
  filter(measured_rating_diff != 'POOR') %>%
  filter(is.finite(SLOPE)==1)

#minimum 20 obs at-a-station
all_rivers <- group_by(all_rivers, site_no) %>%
  filter(n() >= 20)

#prepare NHD for downstream distances and s0 and dimminution coefficient
#SOURCE_FEAT is the gage ID snapped to the reaches, this file is every gage ID snapped to NHD med (with slopes and absolute sum and Measure, all needed)
med_NHD <- semi_join(med_NHD, all_rivers, by=c('SOURCE_FEA' = 'site_no'))
med_NHD$distDwnStrm <- med_NHD$ArbolateSu - ((med_NHD$Measure/100)*med_NHD$LENGTHKM) #measure is % from downstream end of reach that gage is located, ArboSu is total upstream river distance from downstream end of reach
all_rivers <- left_join(all_rivers, select(med_NHD, distDwnStrm, SOURCE_FEA), by = c("site_no" = "SOURCE_FEA"))

Slope0 <- group_by(all_rivers, river_name, site_no) %>%
  slice(1) %>%
  group_by(river_name) %>%
  slice(which.min(distDwnStrm))
Slope0 = Slope0 %>% rename('S0' = 'SLOPE')
all_rivers <- left_join(all_rivers, select(Slope0, S0, river_name), by = c("river_name" = "river_name"))

dimm <- group_by(all_rivers, river_name, site_no) %>%
  slice(1) %>%
  group_by(river_name) %>%
  summarise(dimm_coeff = ((n()*sum(distDwnStrm*log(SLOPE)))-(sum(distDwnStrm)*sum(log(SLOPE))))/
           ((n()*sum(distDwnStrm^2))-(sum(distDwnStrm))^2))
all_rivers <- left_join(all_rivers, select(dimm, dimm_coeff, river_name), by = c("river_name" = "river_name"))

#imperial to metric
all_rivers$chan_depth = all_rivers$chan_area/all_rivers$chan_width

all_rivers$chan_width = all_rivers$chan_width*0.305 #m
all_rivers$chan_depth = all_rivers$chan_depth*0.305 #m
all_rivers$chan_discharge = all_rivers$chan_discharge*0.028 #m3/s
all_rivers$chan_velocity = all_rivers$chan_velocity*0.305 #m/s

#calculate AHG
AHG_depth = group_by(all_rivers, site_no)  %>%
  do(model=lm(log(chan_depth) ~ log(chan_discharge), data = .)) %>%
  mutate(log_c = model$coefficients[1]) %>%
  mutate(c = exp(log_c)) %>%
  mutate(f = model$coefficients[2]) %>%
  mutate(ahg_depth_r2 = summary(model)$r.squared)
all_rivers <- left_join(all_rivers,AHG_depth, by="site_no")

AHG_width = group_by(all_rivers, site_no)  %>%
  do(model=lm(log(chan_width) ~ log(chan_discharge), data = .)) %>%
  mutate(log_a = model$coefficients[1]) %>%
  mutate(a = exp(log_a)) %>%
  mutate(b = model$coefficients[2]) %>%
  mutate(ahg_width_r2 = summary(model)$r.squared)
all_rivers <- left_join(all_rivers,AHG_width, by="site_no")

AHG_velocity = group_by(all_rivers, site_no)  %>%
  do(model=lm(log(chan_velocity) ~ log(chan_discharge), data = .)) %>%
  mutate(log_k = model$coefficients[1]) %>%
  mutate(k = exp(log_k)) %>%
  mutate(m = model$coefficients[2]) %>%
  mutate(ahg_velocity_r2 = summary(model)$r.squared)
all_rivers <- left_join(all_rivers,AHG_velocity, by="site_no")

#calculate AMHG
widthAMHG=group_by(all_rivers, river_name, site_no) %>%
  slice(1) %>%
  group_by(river_name) %>%
  do(model=lm(log_a~b,data= .)) %>%
  mutate(logQc_w=-model$coefficients[2]) %>%
  mutate(logwc=model$coefficients[1]) %>%
  mutate(wAMHGr2= summary(model)$r.squared)
all_rivers <- left_join(all_rivers, widthAMHG, by="river_name")

depthAMHG=group_by(all_rivers, river_name, site_no) %>%
  slice(1) %>%
  group_by(river_name) %>%
  do(model=lm(log_c~f,data= .)) %>%
  mutate(logQc_d=-model$coefficients[2]) %>%
  mutate(logdc=model$coefficients[1]) %>%
  mutate(dAMHGr2= summary(model)$r.squared)
all_rivers <- left_join(all_rivers, depthAMHG, by="river_name")

velocityAMHG=group_by(all_rivers, river_name, site_no) %>%
  slice(1) %>%
  group_by(river_name) %>%
  do(model=lm(log_k~m,data= .)) %>%
  mutate(logQc_v=-model$coefficients[2]) %>%
  mutate(logvc=model$coefficients[1]) %>%
  mutate(vAMHGr2= summary(model)$r.squared)
all_rivers <- left_join(all_rivers, velocityAMHG, by="river_name")

#for figure S2
yough = subset(all_rivers, river_name == 'youghiogheny r')
yough_plot <- ggplot(yough, aes(logQc_d, logdc)) +geom_smooth(aes(log(chan_discharge), log(chan_depth), group=site_no, color="#d95f02"), method="lm", se=FALSE, fullrange=TRUE)+geom_point(colour="black", size=5) + xlim(-10,15) + ggtitle("Youghiogheny River Depth AHG Rating Curves") + xlab("log(Q) m3/s") + ylab("log(Depth) m") + geom_vline(aes(xintercept=4.27), colour="#7570b3", size=1)+ geom_vline(aes(xintercept=2.43), colour="#1b9e77", size=1)

#for figure 1
bigsioux <- subset(all_rivers, river_name == 'big sioux r')
bigsioux_plot <- ggplot(bigsioux, aes(logQc_d, logdc)) +geom_smooth(aes(log(chan_discharge), log(chan_depth), group=site_no, color="#d95f02"), method="lm", se=FALSE, fullrange=TRUE)+geom_point(colour="black", size=5) + xlim(0,7) + ylim(-5,5) + ggtitle("AHG Rating Curves: Big Sioux River") + xlab(expression(log ~ Q ~ (m^3/s))) + ylab("log Depth (m)") + theme(legend.position = 'none')
bigsioux_plot2 <- ggplot(bigsioux, aes(log_c, f)) + geom_point(size=3) + geom_smooth(method='lm', se=FALSE, col='#7570b3') + annotate(geom = 'text', x=-1, y=.55, label=expression(r^2: ~ 0.29)) + xlab("log AHG coefficient c") + ylab('AHG exponent f') + ggtitle('AMHG: Big Sioux River')

ohio <- subset(all_rivers, river_name == 'ohio r')
ohio_plot <- ggplot(ohio, aes(logQc_d, logdc)) +geom_smooth(aes(log(chan_discharge), log(chan_depth), group=site_no, color="#d95f02"), method="lm", se=FALSE, fullrange=TRUE)+geom_point(colour="black", size=5) + xlim(-10,15) + ylim(-5,5) + ggtitle("AHG Rating Curves: Ohio RIver") + xlab(expression(log ~ Q ~ (m^3/s))) + ylab("log Depth (m)")+ theme(legend.position = 'none')
ohio_plot2 <- ggplot(ohio, aes(log_c, f)) + geom_point(size=3) + geom_smooth(method = 'lm', se=FALSE, col='#7570b3')  + annotate(geom = 'text', x=1, y=0.37, label=expression(r^2: ~ 0.98)) + xlab("log AHG coefficient c") + ylab('AHG exponent f') + ggtitle('AMHG: Ohio River')

#Calculate conductance coefficient
conductance = group_by(all_rivers, site_no)  %>%
  do(model=lm(log(chan_velocity) ~ log(chan_depth), data = .)) %>%
  mutate(log_conductanceK = model$coefficients[1]) %>%
  mutate(conductanceK = exp(log_conductanceK)) %>%
  mutate(conductanceP = model$coefficients[2]) %>%
  mutate(conductanceR2 = summary(model)$r.squared)
all_rivers <- left_join(all_rivers,conductance, by="site_no")

#Calibrate D from Ostendorf math
final_stationD <- read.csv("calibratedD.csv")
all_rivers <- merge(all_rivers, final_stationD, by.x="site_no", by.y='station')



#more parameters------------------------------------------------------------------
#calculate conductance K
all_rivers$K <- all_rivers$conductanceK / all_rivers$SLOPE^(1/2)

#outlier scrubber
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(0.25, 0.75), na.rm = na.rm, ...)
  H <- 1.5*IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

#bankfull hydraulics, defined as 2 yr return period of at-a-station flows
bankfullW = group_by(all_rivers, site_no) %>%
  mutate(rank = rank(chan_width, ties.method="first")) %>%
  mutate(n = n()) %>%
  mutate(desiredRank = round((n+1)/2, digits=0)) %>% #2 for a two year return period
  mutate(temp = which(rank == desiredRank)) %>%
  summarise(bankful_width = max(chan_width[temp]))
all_rivers <- left_join(all_rivers,bankfullW, by="site_no")

bankfullD <- group_by(all_rivers, site_no) %>%
  mutate(rank = rank(chan_depth, ties.method="first")) %>%
  mutate(n = n()) %>%
  mutate(desiredRank = round((n+1)/2, digits=0)) %>% #2 for a two year return period
  mutate(temp = which(rank == desiredRank)) %>%
  summarise(bankful_depth = max(chan_depth[temp]))
all_rivers <- left_join(all_rivers,bankfullD, by="site_no")

bankfullV <- group_by(all_rivers, site_no) %>%
  mutate(rank = rank(chan_velocity, ties.method="first")) %>%
  mutate(n = n()) %>%
  mutate(desiredRank = round((n+1)/2, digits=0)) %>% #2 for a two year return period
  mutate(temp = which(rank == desiredRank)) %>%
  summarise(bankful_vel = max(chan_velocity[temp]))
all_rivers <- left_join(all_rivers,bankfullV, by="site_no")

bankfullQ <- group_by(all_rivers, site_no) %>%
  mutate(rank = rank(chan_discharge, ties.method="first")) %>%
  mutate(n = n()) %>%
  mutate(desiredRank = round((n+1)/2, digits=0)) %>% #2 for a two year return period
  mutate(temp = which(rank == desiredRank)) %>%
  summarise(bankful_Q = max(chan_discharge[temp]))
all_rivers <- left_join(all_rivers,bankfullQ, by="site_no")

#calculate Sw for r shape parameter
r_by_group = group_by(all_rivers, site_no)  %>%
  do(model=lm(log(chan_width/bankful_width) ~ log(chan_depth/bankful_depth), data=.)) %>%
  mutate(Sw = model$coefficients[2]) %>% #SLOPE of function
  mutate(r_r2 = summary(model)$r.squared) #degree of fit of an r-channel
all_rivers <- left_join(all_rivers, r_by_group, by="site_no")

#Filter out rivers with few stations
all_rivers = group_by(all_rivers, site_no) %>%
  slice(1) %>%
  group_by(river_name) %>%
  filter(n() >= 6)

#calculate r, delta, and R
all_rivers$r <- (1/all_rivers$Sw) #from Dingman & Afshari, 2018
all_rivers$delta <- 1 + all_rivers$r + all_rivers$r*all_rivers$conductanceP
all_rivers$R <- (1+all_rivers$r)/all_rivers$r

#anlaytical expressions for AHG
all_rivers$aDingman <- all_rivers$bankful_width^((all_rivers$r+all_rivers$conductanceP*all_rivers$r)/all_rivers$delta)*(all_rivers$bankful_depth/all_rivers$R)^((1+all_rivers$conductanceP)/all_rivers$delta)*all_rivers$conductanceK^(1/all_rivers$delta)
all_rivers$cDingman <- all_rivers$bankful_width^(-all_rivers$r/all_rivers$delta)*(all_rivers$bankful_depth/all_rivers$R)^(1/all_rivers$delta)*(all_rivers$conductanceK)^(-all_rivers$r/all_rivers$delta)
all_rivers$kDingman <- all_rivers$bankful_width^((-all_rivers$r*all_rivers$conductanceP)/all_rivers$delta)*(all_rivers$bankful_depth/all_rivers$R)^(1/all_rivers$delta)*all_rivers$conductanceK^((1+all_rivers$r)/all_rivers$delta)

all_rivers$bDingman <- 1/all_rivers$delta
all_rivers$fDingman <- all_rivers$r/all_rivers$delta
all_rivers$mDingman <- (all_rivers$conductanceP*all_rivers$r)/all_rivers$delta

#inject congruent hydraulics
all_rivers$varyingwidth <- all_rivers$bankful_width^(all_rivers$r+(all_rivers$r*all_rivers$conductanceP))*(all_rivers$bankful_depth/all_rivers$R)^(-(1+all_rivers$conductanceP))*(all_rivers$conductanceK)^(-1)
all_rivers$varyingdepth <- all_rivers$bankful_width^(-1)*(all_rivers$bankful_depth/all_rivers$R)^(1/all_rivers$r)*(all_rivers$conductanceK)^(-1)
all_rivers$varyingvelocity <- all_rivers$bankful_width^(-1)*(all_rivers$bankful_depth/all_rivers$R)^(1/all_rivers$r)*(all_rivers$conductanceK)^((1+all_rivers$r)/(all_rivers$r*all_rivers$conductanceP))

all_rivers$congruentwidthTest <- exp(all_rivers$logwc)^(all_rivers$delta)*exp(all_rivers$logQc_w)^(-1)
all_rivers$congruentdepthTest <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_d)^(-1)
all_rivers$congruentvelocityTest <- exp(all_rivers$logvc)^(all_rivers$delta/(all_rivers$r*all_rivers$conductanceP))*exp(all_rivers$logQc_v)^(-1)

all_rivers <- filter(all_rivers, is.finite(congruentvelocityTest)==TRUE, is.finite(congruentwidthTest)==TRUE, is.finite(varyingvelocity)==TRUE, is.finite(varyingwidth)==TRUE, congruentwidthTest>0, congruentvelocityTest>0, varyingvelocity >0, varyingwidth>0)

subset20 <- filter(all_rivers, ahg_depth_r2 > 0.20, ahg_width_r2 > 0.20, ahg_velocity_r2 > 0.20)
subsetdepth30 <- filter(all_rivers, ahg_depth_r2 > 0.5, ahg_width_r2 > 0.5, ahg_velocity_r2 > 0.5)
subsetdepth85 <- filter(all_rivers, ahg_depth_r2 > 0.85, ahg_width_r2 > 0.85, ahg_velocity_r2 > 0.85)

#Run linear regressions for figure 2
lm1 <- data.frame(
  depthlm1 <- summary(lm(log10(remove_outliers(subset20$congruentdepthTest))~log10(remove_outliers(subset20$varyingdepth))))$r.squared,
  depthlm2 <- summary(lm(log10(remove_outliers(subsetdepth30$congruentdepthTest))~log10(remove_outliers(subsetdepth30$varyingdepth))))$r.squared,
  depthlm3 <- summary(lm(log10(remove_outliers(subsetdepth85$congruentdepthTest))~log10(remove_outliers(subsetdepth85$varyingdepth))))$r.squared,
  
  widthlm1 <- summary(lm(log10(remove_outliers(subset20$congruentwidthTest))~log10(remove_outliers(subset20$varyingwidth))))$r.squared,
  widthlm2 <- summary(lm(log10(remove_outliers(subsetdepth30$congruentwidthTest))~log10(remove_outliers(subsetdepth30$varyingwidth))))$r.squared,
  widthlm3 <- summary(lm(log10(remove_outliers(subsetdepth85$congruentwidthTest))~log10(remove_outliers(subsetdepth85$varyingwidth))))$r.squared,
  
  velocitylm1 <- summary(lm(log10(remove_outliers(subset20$congruentvelocityTest))~log10(remove_outliers(subset20$varyingvelocity))))$r.squared,
  velocitylm2 <- summary(lm(log10(remove_outliers(subsetdepth30$congruentvelocityTest))~log10(remove_outliers(subsetdepth30$varyingvelocity))))$r.squared,
  velocitylm3 <- summary(lm(log10(remove_outliers(subsetdepth85$congruentvelocityTest))~log10(remove_outliers(subsetdepth85$varyingvelocity))))$r.squared
)
lm1 <- t(lm1)

#Figure 2: plot congruent hydraulics in AHG
velocity1 <- ggplot(subset20, aes(log10(varyingvelocity), log10(congruentvelocityTest))) + geom_point(size = 4, color='#1b9e77') + geom_abline(size=1) + labs(title="Velocity") + annotate("text", x = -5, y = 10, label = expression(r^2: ~ 0.44)) + xlim(-10,15) + ylim(-10,15)
velocity2 <- ggplot(subsetdepth30, aes(log10(varyingvelocity), log10(congruentvelocityTest))) + geom_point(size = 4, color='#1b9e77') + geom_abline(size=1) + annotate("text", x = -5, y = 10, label = expression(r^2: ~ 0.49)) + xlim(-10,15) + ylim(-10,15)
velocity3 <- ggplot(subsetdepth85, aes(log10(varyingvelocity), log10(congruentvelocityTest))) + geom_point(size = 4, color='#1b9e77') + geom_abline(size=1) + annotate("text", x = -5, y = 10, label = expression(r^2: ~ 0.77)) + xlim(-10,15) + ylim(-10,15)
depth1 <- ggplot(subset20, aes(log10(varyingdepth), log10(congruentdepthTest))) + geom_point(size = 4, color='#d95f02') + geom_abline(size=1) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Depth")  + xlim(-2.5,2) + ylim(-2.5,2) + annotate("text", x = -1.5, y = 1, label = expression(r^2: ~ 0.28))
depth2 <- ggplot(subsetdepth30, aes(log10(varyingdepth), log10(congruentdepthTest))) + geom_point(size = 4, color='#d95f02') + geom_abline(size=1) + xlim(-2.5,2) + ylim(-2.5,2) + annotate("text", x = -1.5, y = 1, label = expression(r^2: ~ 0.32))
depth3 <- ggplot(subsetdepth85, aes(log10(varyingdepth), log10(congruentdepthTest))) + geom_point(size = 4, color='#d95f02') + geom_abline(size=1) + xlim(-2.5,2) + ylim(-2.5,2) + annotate("text", x = -1.5, y = 1, label = expression(r^2: ~ 0.51))
width1 <- ggplot(subset20, aes(log10(varyingwidth), log10(congruentwidthTest))) + geom_point(size = 4, color='#7570b3') + geom_abline(size=1)+ labs(title="Width")  + xlim(-300,300) + ylim(-300,300)
width2 <- ggplot(subsetdepth30, aes(log10(varyingwidth), log10(congruentwidthTest))) + geom_point(size = 4, color='#7570b3') + geom_abline(size=1) + xlim(-300,300) + ylim(-300,300)
width3 <- ggplot(subsetdepth85, aes(log10(varyingwidth), log10(congruentwidthTest))) + geom_point(size = 4, color='#7570b3') + geom_abline(size=1) + xlim(-300,300) + ylim(-300,300)

grid1 <- plot_grid(depth1+theme(legend.position="none", axis.title = element_blank()), 
                   width1+theme(legend.position="none", axis.title = element_blank()), 
                   velocity1+theme(legend.position="none", axis.title = element_blank()), 
                   depth2+theme(legend.position="none", axis.title = element_blank()), 
                   width2+theme(legend.position="none", axis.title = element_blank()), 
                   velocity2+theme(legend.position="none", axis.title = element_blank()), 
                   depth3+theme(legend.position="none", axis.title = element_blank()), 
                   width3+theme(legend.position="none", axis.title = element_blank()), 
                   velocity3+theme(legend.position="none", axis.title = element_blank()), ncol = 3)
plot1 <- plot_grid(grid1, rel_widths=c(3,0.7)) 
yTitle1 <- textGrob(expression(log[10]~River~Wide~Term), gp=gpar(fontface="bold", col="black", fontsize=18), rot=90)
xTitle1 <- textGrob(expression(log[10]~At~a~Station~Term), gp=gpar(fontface="bold", col="black", fontsize=18))

pdf("fig2.pdf", width = 9, height = 7) # Open a new pdf file
grid.arrange(arrangeGrob(plot1, left = yTitle1, bottom = xTitle1))
dev.off()

#Riv-wide models------------------------------------------------------------------------------
#Fb
all_rivers$bankful_Fb <- (all_rivers$bankful_vel)/(sqrt(all_rivers$bankful_depth*9.8))
all_rivers$chan_Fb <- (all_rivers$chan_velocity)/(sqrt(all_rivers$chan_depth*9.8))
fb_s <- group_by(all_rivers, site_no) %>%
  summarise(bankful_Fb = mean(bankful_Fb), fb_sSLOPE = mean(SLOPE))

model=lm(log(bankful_Fb)~log(fb_sSLOPE),data= fb_s)
fb_s$Fbx <- exp(model$coefficients[1])
fb_s$Fby <- model$coefficients[2]
fb_s$FbvsS_r2 <- summary(model)$r.squared

all_rivers <- left_join(all_rivers, fb_s, by='site_no')

all_rivers$BjerklieSLOPE <- (all_rivers$chan_Fb/1.07)^(1/(0.29))
all_rivers$BjerklieGlobalSLOPE <- (all_rivers$chan_Fb/2.85)^(1/(0.31))

#Einstein & Barbarossa
SLOPEvsD=group_by(all_rivers, river_name, site_no) %>%
  slice(1) %>%
  group_by(river_name) %>%
  do(model=lm(log(D)~log(SLOPE),data= .)) %>%
  mutate(capB=model$coefficients[2]) %>%
  mutate(logcapA=model$coefficients[1]) %>%
  mutate(capA=exp(logcapA))%>%
  mutate(SvsD_r2= summary(model)$r.squared)
all_rivers <- left_join(all_rivers, SLOPEvsD, by='river_name')

all_rivers$EinsteinSLOPE <- (all_rivers$D/all_rivers$capA)^(1/(all_rivers$capB))

#DHG
slopeDHG=group_by(all_rivers, river_name, site_no) %>%
  slice(1) %>%
  group_by(river_name) %>%
  do(model=lm(log(SLOPE)~log(bankful_Q),data= .)) %>%
  mutate(capY=model$coefficients[2]) %>%
  mutate(logcapX=model$coefficients[1]) %>%
  mutate(capX=exp(logcapX))%>%
  mutate(SDHG_r2= summary(model)$r.squared)
all_rivers <- left_join(all_rivers, slopeDHG, by='river_name')

all_rivers$dhgSLOPE <- all_rivers$capX*all_rivers$bankful_Q^all_rivers$capY

#Regime Theory
all_rivers$regimeD <- all_rivers$SLOPE*11*all_rivers$bankful_depth #(all_rivers$bankful_width*3.28/(0.93*(all_rivers$bankful_Q*35.3)^0.46))^(-0.15) #eq 10-67
all_rivers$RegimeSLOPE <- 0.44*all_rivers$regimeD^(1.15)*all_rivers$bankful_Q^(-0.46)

#graded profile
all_rivers$GradedSLOPE <- (all_rivers$S0*exp(all_rivers$dimm_coeff*all_rivers$distDwnStrm))

all_rivers$varyingwidth2 <- all_rivers$bankful_width^(all_rivers$r+(all_rivers$r*all_rivers$conductanceP))*(all_rivers$bankful_depth/all_rivers$R)^(-(1+all_rivers$conductanceP))*(all_rivers$K)^(-1)
all_rivers$varyingdepth2 <- all_rivers$bankful_width^(-1)*(all_rivers$bankful_depth/all_rivers$R)^(1/all_rivers$r)*(all_rivers$K)^(-1)
all_rivers$varyingvelocity2 <- all_rivers$bankful_width^(-1)*(all_rivers$bankful_depth/all_rivers$R)^(1/all_rivers$r)*(all_rivers$K)^((1+all_rivers$r)/(all_rivers$r*all_rivers$conductanceP))

#introduce river-wide models into expressions
all_rivers$congruentwidthEinstein <- exp(all_rivers$logwc)^(all_rivers$delta)*exp(all_rivers$logQc_w)^(-1)*all_rivers$EinsteinSLOPE^(1/2)
all_rivers$congruentwidthBjerklieGlobal <- exp(all_rivers$logwc)^(all_rivers$delta)*exp(all_rivers$logQc_w)^(-1)*all_rivers$BjerklieGlobalSLOPE^(1/2)
all_rivers$congruentwidthBjerklie <- exp(all_rivers$logwc)^(all_rivers$delta)*exp(all_rivers$logQc_w)^(-1)*all_rivers$BjerklieSLOPE^(1/2)
all_rivers$congruentwidthGraded <- exp(all_rivers$logwc)^(all_rivers$delta)*exp(all_rivers$logQc_w)^(-1)*all_rivers$GradedSLOPE^(1/(2))
all_rivers$congruentwidthRegime <- exp(all_rivers$logwc)^(all_rivers$delta)*exp(all_rivers$logQc_w)^(-1)*all_rivers$RegimeSLOPE^(1/2)
all_rivers$congruentwidthDHG <- exp(all_rivers$logwc)^(all_rivers$delta)*exp(all_rivers$logQc_w)^(-1)*all_rivers$dhgSLOPE^(1/2)

all_rivers$congruentdepthEinstein <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_d)^(-1)*all_rivers$EinsteinSLOPE^(1/2)
all_rivers$congruentdepthBjerklieGlobal <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_d)^(-1)*all_rivers$BjerklieGlobalSLOPE^(1/2)
all_rivers$congruentdepthBjerklie <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_d)^(-1)*all_rivers$BjerklieSLOPE^(1/2)
all_rivers$congruentdepthGraded <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_d)^(-1)*all_rivers$GradedSLOPE^(1/(2))
all_rivers$congruentdepthRegime <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_d)^(-1)*all_rivers$RegimeSLOPE^(1/2)
all_rivers$congruentdepthDHG <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_d)^(-1)*all_rivers$dhgSLOPE^(1/2)

all_rivers$congruentvelocityEinstein <- exp(all_rivers$logvc)^(all_rivers$delta/(all_rivers$r*all_rivers$conductanceP))*exp(all_rivers$logQc_v)^(-1)*all_rivers$EinsteinSLOPE^(-(1+all_rivers$r)/(2*all_rivers$r*all_rivers$conductanceP))
all_rivers$congruentvelocityBjerklieGlobal <- exp(all_rivers$logvc)^(all_rivers$delta/(all_rivers$r*all_rivers$conductanceP))*exp(all_rivers$logQc_v)^(-1)*all_rivers$BjerklieGlobalSLOPE^(-(1+all_rivers$r)/(2*all_rivers$r*all_rivers$conductanceP))
all_rivers$congruentvelocityBjerklie <- exp(all_rivers$logvc)^(all_rivers$delta/(all_rivers$r*all_rivers$conductanceP))*exp(all_rivers$logQc_v)^(-1)*all_rivers$BjerklieSLOPE^(-(1+all_rivers$r)/(2*all_rivers$r*all_rivers$conductanceP))
all_rivers$congruentvelocityGraded <- exp(all_rivers$logvc)^(all_rivers$delta/(all_rivers$r*all_rivers$conductanceP))*exp(all_rivers$logQc_v)^(-1)*all_rivers$GradedSLOPE^(-(1+all_rivers$r)/(2*all_rivers$r*all_rivers$conductanceP))
all_rivers$congruentvelocityRegime <- exp(all_rivers$logvc)^(all_rivers$delta/(all_rivers$r*all_rivers$conductanceP))*exp(all_rivers$logQc_v)^(-1)*all_rivers$RegimeSLOPE^(-(1+all_rivers$r)/(2*all_rivers$r*all_rivers$conductanceP))
all_rivers$congruentvelocityDHG <- exp(all_rivers$logvc)^(all_rivers$delta/(all_rivers$r*all_rivers$conductanceP))*exp(all_rivers$logQc_v)^(-1)*all_rivers$dhgSLOPE^(-(1+all_rivers$r)/(2*all_rivers$r*all_rivers$conductanceP))

#filter out infinites for lm regressions
all_rivers <- filter(all_rivers, is.finite(congruentwidthTest)==1) %>%
  filter(is.finite(congruentvelocityTest) == 1) %>%
  filter(congruentwidthTest >0) %>%
  filter(is.finite(varyingwidth)==1) %>%
  filter(varyingwidth > 0) %>%
  filter(is.finite(congruentvelocityBjerklie) ==1)

#must remove crazy outliers for t-test to run
velSubset <- filter(all_rivers, congruentvelocityEinstein < 1e+20) %>%
  filter(is.finite(congruentvelocityBjerklie) ==1) %>%
  filter(congruentvelocityBjerklieGlobal < 1e+20) %>%
  filter(congruentvelocityBjerklie > 0) %>%
  filter(congruentvelocityBjerklie < 1e+20) %>%
  filter(congruentvelocityGraded < 1e+20)

#tests on riv-wide models---------------------------------------------------------------
#modelFits
bjerklieglobalFits = group_by(all_rivers, river_name, site_no) %>%
  slice(1)%>%
  group_by(river_name) %>%
  mutate(model = (chan_Fb/2.85)^(1/0.31)) %>%
  do(lm=lm(log10(SLOPE)~log10(model), data= .), lm2=lm(log10(congruentdepthTest)~log10(varyingdepth),data= .), lm3=lm(log10(congruentwidthTest)~log10(varyingwidth),data= all_rivers), lm4=lm(log10(congruentvelocityTest)~log10(varyingvelocity),data= .)) %>%
  mutate(bjerklieglobalr2 = summary(lm)$r.squared) %>%
  mutate(congruentdepthr2 = summary(lm2)$r.squared) %>%
  mutate(congruentwidthr2 = summary(lm3)$r.squared) %>%
  mutate(congruentvelocityr2 = summary(lm4)$r.squared)

bjerklieFits = group_by(all_rivers, river_name, site_no) %>%
  slice(1)%>%
  group_by(river_name) %>%
  mutate(model = (chan_Fb/1.07)^(1/0.29)) %>%
  do(lm=lm(log10(SLOPE)~log10(model),data= ., na.action = na.omit)) %>%
  mutate(bjerklier2 = summary(lm)$r.squared)

regimeFits = group_by(all_rivers, river_name, site_no) %>%
  slice(1)%>%
  group_by(river_name) %>%
  mutate(model = 0.44*regimeD*0.0254^(1.15)*bankful_Q^(-0.46)) %>%
  do(lm=lm(log10(SLOPE)~log10(model),data= .), lm2=lm(log10(congruentdepthTest)~log10(varyingdepth),data= .), lm3=lm(log10(congruentwidthTest)~log10(varyingwidth),data= all_rivers), lm4=lm(log10(congruentvelocityTest)~log10(varyingvelocity),data= .)) %>%
  mutate(regimer2 = summary(lm)$r.squared) %>%
  mutate(congruentdepthr2 = summary(lm2)$r.squared) %>%
  mutate(congruentwidthr2 = summary(lm3)$r.squared) %>%
  mutate(congruentvelocityr2 = summary(lm4)$r.squared)

gradedFits = group_by(all_rivers, river_name, site_no) %>%
  slice(1)%>%
  group_by(river_name) %>%
  mutate(model = S0*exp(dimm_coeff*distDwnStrm)) %>%
  do(lm=lm(log10(SLOPE)~log10(model),data= .), lm2=lm(log10(congruentdepthTest)~log10(varyingdepth),data= .), lm3=lm(log10(congruentwidthTest)~log10(varyingwidth),data= all_rivers), lm4=lm(log10(congruentvelocityTest)~log10(varyingvelocity),data= .)) %>%
  mutate(gradedr2 = summary(lm)$r.squared) %>%
  mutate(congruentdepthr2 = summary(lm2)$r.squared) %>%
  mutate(congruentwidthr2 = summary(lm3)$r.squared) %>%
  mutate(congruentvelocityr2 = summary(lm4)$r.squared)

einsteinFits = group_by(all_rivers, river_name, site_no) %>%
  slice(1)%>%
  group_by(river_name) %>%
  mutate(model = (D/capA)^(1/capB)) %>%
  do(lm=lm(log10(SLOPE)~log10(model),data= .), lm2=lm(log10(congruentdepthTest)~log10(varyingdepth),data= .), lm3=lm(log10(congruentwidthTest)~log10(varyingwidth),data= .), lm4=lm(log10(congruentvelocityTest)~log10(varyingvelocity),data= .)) %>%
  mutate(einsteinr2 = summary(lm)$r.squared) %>%
  mutate(congruentdepthr2 = summary(lm2)$r.squared) %>%
  mutate(congruentwidthr2 = summary(lm3)$r.squared) %>%
  mutate(congruentvelocityr2 = summary(lm4)$r.squared)

dhgFits = group_by(all_rivers, river_name, site_no) %>%
  slice(1)%>%
  group_by(river_name) %>%
  mutate(model = capX*bankful_Q^capY) %>%
  do(lm=lm(log10(SLOPE)~log10(model),data= .), lm2=lm(log10(congruentdepthTest)~log10(varyingdepth),data= .), lm3=lm(log10(congruentwidthTest)~log10(varyingwidth),data= .), lm4=lm(log10(congruentvelocityTest)~log10(varyingvelocity),data= .)) %>%
  mutate(DHGr2 = summary(lm)$r.squared) %>%
  mutate(congruentdepthr2 = summary(lm2)$r.squared) %>%
  mutate(congruentwidthr2 = summary(lm3)$r.squared) %>%
  mutate(congruentvelocityr2 = summary(lm4)$r.squared)

#statistics
depthTestEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentdepthTest', i2='varyingdepth')), high_eqbound = 0.1, low_eqbound = -0.1)
depthEinsteinEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentdepthEinstein', i2='varyingdepth')), high_eqbound = 0.1, low_eqbound = -0.1)
depthGradedEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentdepthGraded', i2='varyingdepth')), high_eqbound = 0.1, low_eqbound = -0.1)
depthRegimeEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentdepthRegime', i2='varyingdepth')), high_eqbound = 0.1, low_eqbound = -0.1)
depthBjerklieEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentdepthBjerklie', i2='varyingdepth')), high_eqbound = 0.1, low_eqbound = -0.1)
depthBjerklieGlobalEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentdepthBjerklieGlobal', i2='varyingdepth')), high_eqbound = 0.1, low_eqbound = -0.1)
depthDHGEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentdepthDHG', i2='varyingdepth')), high_eqbound = 0.1, low_eqbound = -0.1)

widthTestEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentwidthTest', i2='varyingwidth')), high_eqbound = 0.1, low_eqbound = -0.1)
widthEinsteinEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentwidthEinstein', i2='varyingwidth')), high_eqbound = 0.1, low_eqbound = -0.1)
widthGradedEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentwidthGraded', i2='varyingwidth')), high_eqbound = 0.1, low_eqbound = -0.1)
widthRegimeEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentwidthRegime', i2='varyingwidth')), high_eqbound = 0.1, low_eqbound = -0.1)
widthBjerklieEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentwidthBjerklie', i2='varyingwidth')), high_eqbound = 0.1, low_eqbound = -0.1)
widthBjerklieGlobalEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentwidthBjerklieGlobal', i2='varyingwidth')), high_eqbound = 0.1, low_eqbound = -0.1)
widthDHGEqiv <- dataTOSTpaired(all_rivers, pairs=list(c(i1='congruentwidthDHG', i2='varyingwidth')), high_eqbound = 0.1, low_eqbound = -0.1)

velocityTestEqiv <- dataTOSTpaired(velSubset, pairs=list(c(i1='congruentvelocityTest', i2='varyingvelocity')), high_eqbound = 0.1, low_eqbound = -0.1)
velocityEinsteinEqiv <- dataTOSTpaired(velSubset, pairs=list(c(i1='congruentvelocityEinstein', i2='varyingvelocity')), high_eqbound = 0.1, low_eqbound = -0.1, plots = TRUE)
velocityGradedEqiv <- dataTOSTpaired(velSubset, pairs=list(c(i1='congruentvelocityGraded', i2='varyingvelocity')), high_eqbound = 0.1, low_eqbound = -0.1)
velocityRegimeEqiv <- dataTOSTpaired(velSubset, pairs=list(c(i1='congruentvelocityRegime', i2='varyingvelocity')), high_eqbound = 0.1, low_eqbound = -0.1)
velocityBjerklieEqiv <- dataTOSTpaired(velSubset, pairs=list(c(i1='congruentvelocityBjerklie', i2='varyingvelocity')), high_eqbound = 0.1, low_eqbound = -0.1)
velocityBjerklieGlobalEqiv <- dataTOSTpaired(velSubset, pairs=list(c(i1='congruentvelocityBjerklieGlobal', i2='varyingvelocity')), high_eqbound = 0.1, low_eqbound = -0.1)
velocityDHGEqiv <- dataTOSTpaired(velSubset, pairs=list(c(i1='congruentvelocityDHG', i2='varyingvelocity')), high_eqbound = 0.1, low_eqbound = -0.1)

lm <- data.frame(
  depthEinsteinlm <- summary(lm(log(all_rivers$congruentdepthEinstein)~log(all_rivers$varyingdepth2), data = all_rivers))$r.squared,
  depthGradedlm <- summary(lm(log(all_rivers$congruentdepthGraded)~log(all_rivers$varyingdepth2), data = all_rivers))$r.squared,
  depthRegimelm <- summary(lm(log(all_rivers$congruentdepthRegime)~log(all_rivers$varyingdepth2), data = all_rivers))$r.squared,
  depthBjerklielm <- summary(lm(log(all_rivers$congruentdepthBjerklie)~log(all_rivers$varyingdepth2), data = all_rivers))$r.squared,
  depthBjerklieGloballm <- summary(lm(log(all_rivers$congruentdepthBjerklieGlobal)~log(all_rivers$varyingdepth2), data = all_rivers))$r.squared,
  depthDHGlm <- summary(lm(log(all_rivers$congruentdepthDHG)~log(all_rivers$varyingdepth2), data = all_rivers))$r.squared,
  
  widthEinsteinlm <- summary(lm(log(all_rivers$congruentwidthEinstein)~log(all_rivers$varyingwidth2), data = all_rivers))$r.squared,
  widthGradedlm <- summary(lm(log(all_rivers$congruentwidthGraded)~log(all_rivers$varyingwidth2), data = all_rivers))$r.squared,
  widthRegimelm <- summary(lm(log(all_rivers$congruentwidthRegime)~log(all_rivers$varyingwidth2), data = all_rivers))$r.squared,
  widthBjerklielm <- summary(lm(log(all_rivers$congruentwidthBjerklie)~log(all_rivers$varyingwidth2), data = all_rivers))$r.squared,
  widthBjerklieGloballm <- summary(lm(log(all_rivers$congruentwidthBjerklieGlobal)~log(all_rivers$varyingwidth2), data = all_rivers))$r.squared,
  widthDHGlm <- summary(lm(log(all_rivers$congruentwidthDHG)~log(all_rivers$varyingwidth2), data = all_rivers))$r.squared,
  
  velocityEinsteinlm <- summary(lm(log(velSubset$congruentvelocityEinstein)~log(velSubset$varyingvelocity2), data = velSubset))$r.squared,
  velocityGradedlm <- summary(lm(log(velSubset$congruentvelocityGraded)~log(velSubset$varyingvelocity2), data = velSubset))$r.squared,
  velocityRegimelm <- summary(lm(log(velSubset$congruentvelocityRegime)~log(velSubset$varyingvelocity2), data = velSubset))$r.squared,
  velocityBjerklielm <- summary(lm(log(velSubset$congruentvelocityBjerklie)~log(velSubset$varyingvelocity), data = velSubset))$r.squared,
  velocityBjerklieGloballm <- summary(lm(log(velSubset$congruentvelocityBjerklieGlobal)~log(velSubset$varyingvelocity2), data = velSubset))$r.squared,
  velocityDHGlm <- summary(lm(log(velSubset$congruentvelocityDHG)~log(velSubset$varyingvelocity2), data = velSubset))$r.squared,
  
  Einsteinlm <- summary(lm(log(all_rivers$EinsteinSLOPE)~log(all_rivers$SLOPE), data = all_rivers))$r.squared,
  Gradedlm <- summary(lm(log(all_rivers$GradedSLOPE)~log(all_rivers$SLOPE), data = all_rivers))$r.squared,
  Regimelm <- summary(lm(log(all_rivers$RegimeSLOPE)~log(all_rivers$SLOPE), data = all_rivers))$r.squared,
  Bjerklielm <- summary(lm(log(all_rivers$BjerklieSLOPE)~log(all_rivers$SLOPE), data = all_rivers))$r.squared,
  BjerklieGloballm <- summary(lm(log(all_rivers$BjerklieGlobalSLOPE)~log(all_rivers$SLOPE), data = all_rivers))$r.squared,
  DHGlm <- summary(lm(log(all_rivers$dhgSLOPE)~log(all_rivers$SLOPE), data = all_rivers))$r.squared
)
lm <- t(lm)

#depth river-wide model plots-------------------------------------
einsteinPlot <- ggplot(all_rivers, aes(log10(varyingdepth2), log10(congruentdepthEinstein))) + geom_point(size = 3, color='#1b9e77') + geom_abline(size=1) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Depth Einstein") + annotate("text", x = -4, y = 0, label = "r2: 0.71")  + xlim(-6,0.5) + ylim(-6,0.5) #+ annotate("text", x = -4, y = -0.75, label = "nRMSE: 55.5%") #+ annotate("text", x = -1, y = -5, label = "9 Not Plotted")
bjerkliePlot <-ggplot(all_rivers, aes(log10(varyingdepth2), log10(congruentdepthBjerklie))) + geom_point(size = 3, color='#d95f02') + geom_abline(size=1)  + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Depth Bjerklie") + annotate("text", x = -4, y = 0, label = "r2: 0.12")  + xlim(-6,0.5) + ylim(-6,0.5) #+ annotate("text", x = -4, y = -0.75, label = "nRMSE: 115%") # + annotate("text", x = -1, y = -6, label = "12 Not Plotted")
bjerklieglobalPlot <- ggplot(all_rivers, aes(log10(varyingdepth2), log10(congruentdepthBjerklieGlobal))) + geom_point(size = 3, color='#7570b3') + geom_abline(size=1) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Depth Bjerklie Global") + annotate("text", x = -4, y = 0, label = "r2: 0.13")  + xlim(-6,0.5) + ylim(-6,0.5)# + annotate("text", x = -4, y = -0.75, label = "nRMSE: 111%") # + annotate("text", x = -1, y = -6, label = "16 Not Plotted")
gradedPlot <- ggplot(all_rivers, aes(log10(varyingdepth2), log10(congruentdepthGraded))) + geom_point(size = 3, color='#e7298a') + geom_abline(size=1) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Depth Graded Profile") + annotate("text", x = -4, y = 0, label = "r2: 0.35")  + xlim(-6,0.5) + ylim(-6,0.5)# + annotate("text", x = -4, y = -0.75, label = "nRMSE: 87.8%") #+ annotate("text", x = -1, y = -5, label = "10 Not Plotted")
regimePlot <- ggplot(all_rivers, aes(log10(varyingdepth2), log10(congruentdepthRegime))) + geom_point(size = 3, color='#66a61e') + geom_abline(size=1) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Depth Regime Theory") + annotate("text", x = -4, y = 0, label = "r2: 0.82")  + xlim(-6,0.5) + ylim(-6,0.5)# + annotate("text", x = -4, y = -0.75, label = "nRMSE: 50.3%") #+ annotate("text", x = -1, y = -5, label = "9 Not Plotted")
dhgPlot <- ggplot(all_rivers, aes(log10(varyingdepth2), log10(congruentdepthDHG))) + geom_point(size = 3, color='#e6ab02') + geom_abline(size=1) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Depth DHG") + annotate("text", x = -4, y = 0, label = "r2: 0.52")  + xlim(-6,0.5) + ylim(-6,0.5) #+ annotate("text", x = -4, y = -0.75, label = "nRMSE: 69.7%") # + annotate("text", x = -1, y = -5, label = "8 Not Plotted")

#velocity river-wide model plots
einsteinPlotvel <- ggplot(all_rivers, aes(log10(varyingvelocity2), log10(congruentvelocityEinstein))) + geom_point(size = 2, color= '#1b9e77') + geom_abline(size=1) + xlim(-100,100) + ylim(-100,100) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Velocity Einstein") + annotate("text", x = -50, y = 50, label = "r2: 0.98")# + annotate("text", x = -50, y = 25, label = "nRMSE: 21.0%") #+ annotate("text", x = 50, y = -50, label = "7 Not Plotted")
bjerkliePlotvel <-ggplot(all_rivers, aes(log10(varyingvelocity2), log10(congruentvelocityBjerklie))) + geom_point(size = 2, color= '#d95f02') + geom_abline(size=1) +xlim(-100,100) + ylim(-100,100) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Velocity Bjerklie") + annotate("text", x = -50, y = 50, label = "r2: 0.45") #+ annotate("text", x = -50, y = 25, label = "nRMSE: 50.3%") #+ annotate("text", x = 50, y = -50, label = "9 Not Plotted")
bjerklieglobalPlotvel <- ggplot(all_rivers, aes(log10(varyingvelocity2), log10(congruentvelocityBjerklieGlobal))) + geom_point(size = 2, color= '#7570b3') + geom_abline(size=1) +xlim(-100,100) + ylim(-100,100) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Velocity Bjerklie Global") + annotate("text", x = -50, y = 50, label = "r2: 0.81") #+ annotate("text", x = -50, y = 25, label = "nRMSE: 58.3%") #+ annotate("text", x = 50, y = -50, label = "10 Not Plotted")
gradedPlotvel <- ggplot(all_rivers, aes(log10(varyingvelocity2), log10(congruentvelocityGraded))) + geom_point(size = 2, color= '#e7298a') + geom_abline(size=1) + xlim(-100,100) + ylim(-100,100) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Velocity Graded Profile") + annotate("text", x = -50, y = 50, label = "r2: 0.91") #+ annotate("text", x = -50, y = 25, label = "nRMSE: 36.2%") #+ annotate("text", x = 50, y = -50, label = "10 Not Plotted")
regimePlotvel <- ggplot(all_rivers, aes(log10(varyingvelocity2), log10(congruentvelocityRegime))) + geom_point(size = 2, color= '#66a61e') + geom_abline(size=1) + xlim(-100,100) + ylim(-100,100) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Velocity Regime Theory") + annotate("text", x = -50, y = 50, label = "r2: 0.98") #+ annotate("text", x = -50, y = 25, label = "nRMSE: 25.0%")# + annotate("text", x = 50, y = -50, label = "8 Not Plotted")
dhgPlotvel <- ggplot(all_rivers, aes(log10(varyingvelocity2), log10(congruentvelocityDHG))) + geom_point(size = 2, color= '#e6ab02') + geom_abline(size=1) + xlim(-100,100) + ylim(-100,100) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Velocity DHG") + annotate("text", x = -50, y = 50, label = "r2: 0.93") #+ annotate("text", x = -50, y = 25, label = "nRMSE: 28.5%") #+ annotate("text", x = 50, y = -50, label = "9 Not Plotted")

#width river-wide model plots
einsteinPlotwidth <- ggplot(all_rivers, aes(log10(varyingwidth2), log10(congruentwidthEinstein))) + geom_point(size = 3, color = '#1b9e77') + geom_abline(size=1) + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Width Einstein") + annotate("text", x = -150, y = 200, label = "r2: 0.98") #+ annotate("text", x = -150, y = 120, label = "nRMSE: 14.6%") #+ annotate("text", x = 0, y = -50, label = "21 outliers")
bjerkliePlotwidth <-ggplot(all_rivers, aes(log10(varyingwidth2), log10(congruentwidthBjerklie))) + geom_point(size = 3, color = '#d95f02') + geom_abline(size=1)  + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Width Bjerklie") + annotate("text", x = -150, y = 200, label = "r2: 0.98") #+ annotate("text", x = -150, y = 120, label = "nRMSE: 14.6%") #+ annotate("text", x = 0, y = -50, label = "21 outliers")
bjerklieglobalPlotwidth <- ggplot(all_rivers, aes(log10(varyingwidth2), log10(congruentwidthBjerklieGlobal))) + geom_point(size = 3, color = '#7570b3') + geom_abline(size=1)  + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Width Bjerklie Global") + annotate("text", x = -150, y = 200, label = "r2: 0.98") #+ annotate("text", x = -150, y = 120, label = "nRMSE: 14.6%") #+ annotate("text", x = 0, y = -50, label = "21 outliers")
gradedPlotwidth <- ggplot(all_rivers, aes(log10(varyingwidth2), log10(congruentwidthGraded))) + geom_point(size = 3, color = '#e7298a') + geom_abline(size=1)  + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Width Graded Profile") + annotate("text", x = -150, y = 200, label = "r2: 0.98") #+ annotate("text", x = -150, y = 120, label = "nRMSE: 14.8%") #+ annotate("text", x = 0, y = -50, label = "21 outliers")
regimePlotwidth <- ggplot(all_rivers, aes(log10(varyingwidth2), log10(congruentwidthRegime))) + geom_point(size = 3, color = '#66a61e') + geom_abline(size=1)  + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Width Regime Theory") + annotate("text", x = -150, y = 200, label = "r2: 0.98") #+ annotate("text", x = -150, y = 120, label = "nRMSE: 14.7%")# + annotate("text", x = 0, y = -50, label = "21 outliers")
dhgPlotwidth <- ggplot(all_rivers, aes(log10(varyingwidth2), log10(congruentwidthDHG))) + geom_point(size = 3, color = '#e6ab02') + geom_abline(size=1)  + ylab("log10(Congruent Hydraulics)") + xlab("log10(At-a-Station)") + labs(title="Width DHG") + annotate("text", x = -150, y = 200, label = "r2: 0.98") #+ annotate("text", x = -150, y = 120, label = "nRMSE: 14.7%")# + annotate("text", x = 0, y = -50, label = "21 outliers")

#Figure S1
gridwidth <- plot_grid(einsteinPlotwidth+theme(legend.position="none", axis.title = element_blank()), 
                       bjerkliePlotwidth+theme(legend.position="none", axis.title = element_blank()), 
                       bjerklieglobalPlotwidth+theme(legend.position="none", axis.title = element_blank()), 
                       gradedPlotwidth+theme(legend.position="none", axis.title = element_blank()), 
                       regimePlotwidth+theme(legend.position="none", axis.title = element_blank()), 
                       dhgPlotwidth+theme(legend.position="none", axis.title = element_blank()),
                       ncol = 3)
plotwidth <- plot_grid(gridwidth, rel_widths = c(3,0.7))
yTitlewidth <- textGrob(expression(log[10]~River~Wide~Term), gp=gpar(fontface="bold", col="black", fontsize=18), rot=90)
xTitlewidth <- textGrob(expression(log[10]~At~a~Station~Term), gp=gpar(fontface="bold", col="black", fontsize=18))

pdf("figS1.pdf", width = 9, height = 7) # Open a new pdf file
grid.arrange(arrangeGrob(plotwidth, left = yTitlewidth, bottom = xTitlewidth))
dev.off()

#figure 3---------------------------------------------------
gridCombo <- plot_grid(gradedPlot+theme(legend.position="none", axis.title = element_blank()), 
                       gradedPlotvel+theme(legend.position="none", axis.title = element_blank()),
                       regimePlot+theme(legend.position="none", axis.title = element_blank()),
                       regimePlotvel+theme(legend.position="none", axis.title = element_blank()),
                       einsteinPlot+theme(legend.position="none", axis.title = element_blank()), 
                       einsteinPlotvel+theme(legend.position="none", axis.title = element_blank()),
                       bjerkliePlot+theme(legend.position="none", axis.title = element_blank()),
                       bjerkliePlotvel+theme(legend.position="none", axis.title = element_blank()), 
                       bjerklieglobalPlot+theme(legend.position="none", axis.title = element_blank()),
                       bjerklieglobalPlotvel+theme(legend.position="none", axis.title = element_blank()), 
                       dhgPlot+theme(legend.position="none", axis.title = element_blank()), 
                       dhgPlotvel+theme(legend.position="none", axis.title = element_blank()),
                        ncol = 2)
plotCombo <- plot_grid(gridCombo, align='v', axis=c('tblr'))
yTitleCombo <- textGrob(expression(log[10]~River~Wide~Term), gp=gpar(fontface="bold", col="black", fontsize=18), rot=90)
xTitleCombo <- textGrob(expression(log[10]~At~a~Station~Term), gp=gpar(fontface="bold", col="black", fontsize=18))

pdf("fig3.pdf", width = 6, height = 10) # Open a new pdf file
grid.arrange(arrangeGrob(gridCombo, left = yTitleCombo, bottom = xTitleCombo))
dev.off()

#figure4-----------------------------------------
cols <- c(1:6, 13:18)
lmPlot <- as.data.frame(lm[cols,])
colnames(lmPlot)[colnames(lmPlot)=="lm[cols, ]"] <- "recoveryFit"
lmPlot$modelFit <- lm[19:24,]
sig <- c(1:6)
lmPlot$Hydraulic <- ifelse(seq_len(nrow(lmPlot)) %in% sig, "Depth", "Velocity")
lmPlot$Model <- c('Einstein', 'Graded Profile', 'Regime Theory', 'Bjerklie', 'Bjerklie Global', 'DHG', 'Einstein', 'Graded Profile', 'Regime Theory', 'Bjerklie', 'Bjerklie Global', 'DHG')

depthFit <- lm(recoveryFit~modelFit, data = filter(lmPlot, Hydraulic == 'Depth'))
velocityFit <- lm(log(recoveryFit)~log(modelFit), data = filter(lmPlot, Hydraulic == 'Velocity'))
new <- data.frame(modelFit = lmPlot$modelFit)
lmPlot$modeled <- ifelse(lmPlot$Hydraulic == 'Depth', lmPlot$modelFit*0.70687+0.13616, exp(0.03934)*lmPlot$modelFit^0.17922)

library(ggrepel)

fig4 <-ggplot(lmPlot, aes(modelFit, recoveryFit, color = Hydraulic))  + 
  geom_point(size = 3)+
  scale_color_manual(values = c("Depth" = '#1b9e77','Velocity' = '#7570b3')) +
  ylab(expression(Degree ~ of ~ AMHG-AHG ~ Reconciliation ~ (r^2))) + 
  xlab(expression(River ~ Wide ~ Model ~ Strength ~ (r^2))) + 
  xlim(0,1)+
  ylim(0,1) +
  geom_smooth(method='lm', se=FALSE) +
  geom_label_repel(label=lmPlot$Model, show.legend = FALSE) +
  theme(legend.title = element_blank())
fig4
ggsave("fig4.pdf", width = 7, height = 7)

#for figure S2---------------------------------------------------------------------
all_rivers$congruentdepthTestDEPTH <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_d)^(-1)
all_rivers$congruentdepthTestWIDTH <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_w)^(-1)
all_rivers$congruentdepthTestVEL <- exp(all_rivers$logdc)^(all_rivers$delta/all_rivers$r)*exp(all_rivers$logQc_v)^(-1)

depthDEPTH <- ggplot(all_rivers, aes(log10(varyingdepth), log10(congruentdepthTestDEPTH))) + geom_point(size = 2, color='#d95f02') + geom_abline(size=1) + ylab("log10(River-Wide Term)") + xlab("log10(At-a-Station Term)") + labs(title="Depth AMHG Using Qcd") + xlim(-4,2) + ylim(-4,2)
depthWIDTH <- ggplot(all_rivers, aes(log10(varyingdepth), log10(congruentdepthTestWIDTH))) + geom_point(size = 2, color='#d95f02') + geom_abline(size=1) + ylab("log10(River-Wide Term)") + xlab("log10(At-a-Station Term)") + labs(title="Depth AMHG Using Qcw") + xlim(-4,2) + ylim(-4,2)
depthVEL <- ggplot(all_rivers, aes(log10(varyingdepth), log10(congruentdepthTestVEL))) + geom_point(size = 2, color='#d95f02') + geom_abline(size=1) + ylab("log10(River-Wide Term)") + xlab("log10(At-a-Station Term)") + labs(title="Depth AMHG Using Qcv")  + xlim(-4,2) + ylim(-4,2)

gridcongruent <- plot_grid(depthDEPTH+theme(legend.position="none", axis.title = element_blank()), 
                       depthWIDTH+theme(legend.position="none", axis.title = element_blank()), 
                       depthVEL+theme(legend.position="none", axis.title = element_blank()),
                       yough_plot+theme(legend.position="none", axis.title = element_blank()))
plotcongruent <- plot_grid(gridcongruent, rel_widths = c(3,0.7))
yTitlecongruent <- textGrob(expression(log[10]~River~Wide~Term), gp=gpar(fontface="bold", col="black", fontsize=18), rot=90)
xTitlecongruent <- textGrob(expression(log[10]~At~a~Station~Term), gp=gpar(fontface="bold", col="black", fontsize=18))

pdf("figS2.pdf", width = 9, height = 7) # Open a new pdf file
grid.arrange(arrangeGrob(plotcongruent, left = yTitlecongruent, bottom = xTitlecongruent))
dev.off()

#figure 1--------------------------------------------------------------------------------------------------
gridGW15 <- plot_grid(bigsioux_plot,
                      bigsioux_plot2,
                      ohio_plot,
                      ohio_plot2)
plotGW15 <- plot_grid(gridGW15)


grid.arrange(plotGW15)
save_plot("fig1.pdf", plotGW15, ncol=2, base_height = 7, base_width = 5)