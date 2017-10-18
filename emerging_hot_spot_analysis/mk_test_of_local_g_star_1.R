#install.packages('trend')
library('trend')
library('spdep')
#install.packages('rgeos')
library(maptools)
library('reshape2')
library("rgdal")

setwd("~/Documents/Dropbox_1/TACC/emerging_hot_spot_analysis")

dat=read.csv('subset_minimum_tarrant_1.csv')
head(dat)

dat$Date=as.Date(dat[,3], format = "%m/%d/%y")
head(dat)

dat$month=as.numeric(format(dat$Date, "%m"))
head(dat)

# raster grid
res=0.03
xmn=min(dat[,1])-res;ymn=min(dat[,2])-res;xmx=max(dat[,1])+res;ymx=max(dat[,2])+res
library(raster)
r <- raster(xmn=xmn, ymn=ymn, xmx=xmx, ymx=ymx, res=res)
r[] <- 0

## caculate z-value of localG for each cell on the grid over 12 months
month = sort(unique(dat$month))
my.list <- vector('list', length(month))
col_name=vector()

for (i in month){
    print(i)
    one_month_dat=dat[dat$month==i,c("X","Y")]
    
    sp <- SpatialPoints(coords = one_month_dat)
    tab <- table(cellFromXY(r, sp))
    r[as.numeric(names(tab))] <- tab
    
    plot(r)
    points(sp, pch=20)
    
    # x, y is the center of each grid
    d <- data.frame(coordinates(r), count=r[])
    
    xycoords <- cbind(d$x, d$y)
    # lower distance bound=0; upper distance bound =30
    nb3 <- dnearneigh(xycoords, 0, 3*res)
    # localG statist  
    G3 <- localG(d$count, nb2listw(nb3, style="B"))
    col_name=c(col_name,sprintf("locaG_month__%s", i))
    print(col_name)
    my.list[[i]] <- data.frame(as.vector(G3))
}

z_value_DF <- do.call('cbind', my.list)

colnames(z_value_DF)=col_name

#### calculate mk.test statistic
n = dim(z_value_DF)[1]
my.list <- vector('list', n)

for (i in 1:n){
  mk=mk.test(as.vector(t(z_value_DF[i,])))
  mk_z_score=mk$statistic
  mk_p_value=mk$p.value
  my.list[[i]] <- data.frame(mk_z_score,mk_p_value)
}
mk_DF <- do.call('rbind', my.list)

df_xy_mk=cbind(coordinates(r), mk_DF)


###  Plot mk.test z score, the 323 row of mk.test z score are plotted top-left to top-right and 
###  then moving down, ending at bottom-right cell
r <- rasterFromXYZ(df_xy_mk[,c(1,2,3)])
plot(r,col=cm.colors(255))


#### whole dataset
data_all=cbind(df_xy_mk, z_value_DF)
data_all_order_by_mk_z_score=data_all[order(data_all$mk_z_score),]
head(data_all_order_by_mk_z_score)




