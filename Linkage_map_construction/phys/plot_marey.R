library(scales)

for(lg in 1:5) {
	print(lg)
	png(paste0("marey", lg, ".png"), width=1024, height=1024)
	map=read.table(paste0("order", lg, ".mapped"))
	map_tmp = map[map$V1==paste0("LG", lg),]

	
	map_tmp$V4 = -min(map_tmp$V4) + map_tmp$V4
	map_tmp$V3 = -min(map_tmp$V3) + map_tmp$V3
	ymax = max(max(map_tmp$V4), max(map_tmp$V3))
	plot(map_tmp$V2, map_tmp$V3, xlab="Position (Mb)",ylab="Recombination Distance (cM)",xaxt="n", main=paste0("LG", lg), col=rgb(1.0,0,0,0.1), pch=20, cex=1.0, ylim=c(0,ymax))
	points(map_tmp$V2, map_tmp$V4, col=rgb(0,0,1.0,0.1), pch=20, cex=1.0)

#	agp <- read.table(paste0("chr", lg, ".agp"))
#	agp=agp[agp$V1==paste0("LG",lg) & agp$V5=="W",c(1:3)]

#	segments(agp$V2,2*c(1:2),agp$V3,2*c(1:2),col=c("red","blue"),lwd=2,lend=1)
	axis(1,seq(0,500000000,5000000),seq(0,500,5))
#	segments(agp$V2,0,agp$V2,ymax,lwd=0.5,col="darkgray")
	dev.off()
}

