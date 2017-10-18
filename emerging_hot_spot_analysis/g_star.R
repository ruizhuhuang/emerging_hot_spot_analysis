#install.packages('trend')
library('trend')
library('spdep')
# load up area shape file:
install.packages('rgeos')
library(maptools)

setwd("~/Documents/Dropbox_1/TACC/emerging_hot_spot_analysis")

dat=read.csv('subset_minimum_tarrant_1.csv')
head(dat)

dat$Date=as.Date(dat[,3], format = "%m/%d/%y")
head(dat)

dat$month=as.numeric(format(dat$Date, "%m"))
head(dat)

library('reshape2')

# Use this or some other method to add a column of row indices.
dat$row <- with(dat, ave(month==month, month, FUN = seq_along))

m <- melt(dat[,c("X","Y","row","month")], id.vars = c("row", "month"))

data_cube <- acast(m, row ~ variable ~ month)


# raster grid
res=0.03
xmn=min(dat[,1])-res;ymn=min(dat[,2])-res;xmx=max(dat[,1])+res;ymx=max(dat[,2])+res
library(raster)
r <- raster(xmn=xmn, ymn=ymn, xmx=xmx, ymx=ymx, res=res)
r[] <- 0

## caculate z-value of localG for each cell on the grid over 12 months
month = dim(data_cube)[3]
my.list <- vector('list', month)

for (z in 1:month){
    one_month_dat=na.omit(data_cube[ , ,z])
    
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
    col_name=sprintf("z_month_%s", z)
    my.list[[z]] <- data.frame(as.vector(G3))
}

z_value_DF <- do.call('cbind', my.list)

col_name=paste(replicate(month,"locaG_month_"), seq(1,month),sep='')
colnames(z_value_DF)=col_name

#### calculate mk.test statistic
n = dim(z_value_DF)[1]
my.list <- vector('list', n)

for (z in 1:n){
  mk=mk.test(as.vector(t(z_value_DF[z,])))
  mk_z_score=mk$statistic
  mk_p_value=mk$p.value
  my.list[[z]] <- data.frame(mk_z_score,mk_p_value)
}
mk_DF <- do.call('rbind', my.list)

df_xy_mk=cbind(coordinates(r), mk_DF)
x=sort(unique(df_xy_mk[,1]))
y=sort(unique(df_xy_mk[,2]))
z = df_xy_mk[,3]

brks <- seq(-5,5,1)
cm.col <- cm.colors(length(brks)-1)
# x, y have to be increasing unqie x and y axis value, z is plotted from top left to right
image(x, y, t(matrix(z, nrow=17, ncol=19, byrow=TRUE)),
      breaks=brks, col=cm.col, asp=1)

#### whole dataset
data_all=cbind(df_xy_mk, z_value_DF)
data_all_order_by_mk_z_score=data_all[order(data_all$mk_z_score),]
head(data_all_order_by_mk_z_score)

### 
r <- rasterFromXYZ(df_xy_mk[,c(1,2,3)])
plot(r,col=rev(heat.colors(255)))
