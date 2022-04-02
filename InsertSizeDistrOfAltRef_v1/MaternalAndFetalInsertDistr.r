#!/usr/bin/Rscript

# 输入;
myArg <- commandArgs(FALSE)
BinDir = substring(myArg[grep('--file=',myArg)],8)
ScriptName = basename(BinDir)
BinDir = dirname(BinDir)
ArgsId = grep('--args',myArg)
ArgsId = ArgsId + 1
myArg = myArg[ArgsId:length(myArg)]
if(length(myArg) != 7 && length(myArg) != 9 && length(myArg) != 10 && length(myArg) != 11)
{
	stop("The Arguments' number was not 7 9 10 or 11.
	
	Example: Rscript ",ScriptName," Maternal.xls Fetal.xls Image.pdf Title SubTitle XLable YLable Lagend4Data1 Legend4Data2 PlotType PlotSize
	
	本脚本用于展示母本和子代cfDNA片段长度分布。
	
	")
}
MaternalFile = myArg[1]
FetalFile = myArg[2]
ImgFile = myArg[3]
Title = myArg[4]
Sub = myArg[5]
XLab = myArg[6]
YLab = myArg[7]
Legend4Data1 = "Maternal cfDNA"
Legend4Data2 = "Fetal cfDNA"
if(length(myArg) >= 9)
{
	Legend4Data1 = myArg[8]
	Legend4Data2 = myArg[9]
}
PlotType = "l"
if(length(myArg) >= 10)
{
	PlotType = myArg[10]
}
PlotSize = 1
if(length(myArg) >= 11)
{
	PlotSize = as.numeric(myArg[11])
}
PlostWidth = 1.5

# 读入数据;
MaternalData = read.table(MaternalFile,header = F)
FetalData = read.table(FetalFile,header = F)
# 获得X和Y的峰值;
MinX = min(MaternalData$V1,FetalData$V1)
MaxX = max(MaternalData$V1,FetalData$V1)
MaxY = max(MaternalData$V2,FetalData$V2)
Ratio = 1.1
MaxY = MaxY * Ratio

# 保存PDF;
pdf(ImgFile,width = 5,height = 4,family="GB1")
# 图片布局、边界等，par对ggplot无效;
par(mar=c(2.6,2.6,3.1,1.1),mgp=c(1.5,1,0),oma=c(1,1,1,1))
# 作图
plot(x = MaternalData$V1,y = MaternalData$V2,col = "blue",type = PlotType,cex = PlotSize,lwd = PlostWidth,xlim = c(MinX,MaxX),ylim = c(0,MaxY),xlab = "",ylab = "",axes = FALSE,bty = "n",xaxs = "i", yaxs = "i")
par(new=TRUE)
plot(x = FetalData$V1,y = FetalData$V2,col = "red",type = PlotType,cex = PlotSize,lwd = PlostWidth,xlim = c(MinX,MaxX),ylim = c(0,MaxY),xlab = "",ylab = "",axes = FALSE,bty = "n",xaxs = "i", yaxs = "i")
mtext(Title,side = 3,line = 2,cex = 1.5,font = 2,outer=FALSE)
mtext(Sub,side = 3,line = 1,cex = 0.75,font = 1,outer=FALSE)
title(xlab = XLab,ylab = YLab,cex.lab = 1.1,font.lab = 2)
axis(side = 1, lwd = 2,lwd.ticks = 2,padj = -0.9,tck = 0.02)
axis(side = 2, lwd = 2,lwd.ticks = 2,padj = 0.9,tck = 0.02)
legend("topright",legend = c(Legend4Data1,Legend4Data2),col = c("blue","red"),lty = 1,lwd = 2,bty = "n",cex = 0.7)

# 标线;
Range = length(MaternalData$V1) - 4
for (Id in 5:Range) {
	if (MaternalData$V2[Id] > MaternalData$V2[Id - 1] && MaternalData$V2[Id - 1] > MaternalData$V2[Id - 2] && MaternalData$V2[Id - 2] > MaternalData$V2[Id - 3] && MaternalData$V2[Id - 3] > MaternalData$V2[Id - 4] && MaternalData$V2[Id] > MaternalData$V2[Id + 1] && MaternalData$V2[Id + 1] > MaternalData$V2[Id + 2] && MaternalData$V2[Id + 2] > MaternalData$V2[Id + 3] && MaternalData$V2[Id + 3] > MaternalData$V2[Id + 4]) {
		lines(c(MaternalData$V1[Id],MaternalData$V1[Id]),c(0,MaternalData$V2[Id]),col = "#696969",lty = 2, lwd = 1)
		text(x = MaternalData$V1[Id],y = MaternalData$V2[Id],MaternalData$V1[Id],cex = 0.7,adj = c(0,0))
	}
}
Range = length(FetalData$V1) - 5
for (Id in 6:Range) {
	LeftOn = 0
	if (FetalData$V2[Id] > FetalData$V2[Id - 1]) {
		LeftOn = LeftOn + 1
	}
	if (FetalData$V2[Id - 1] > FetalData$V2[Id - 2] || FetalData$V2[Id] > FetalData$V2[Id - 2]) {
		LeftOn = LeftOn + 1
	}
	if (FetalData$V2[Id - 2] > FetalData$V2[Id - 3] || FetalData$V2[Id] > FetalData$V2[Id - 3]) {
		LeftOn = LeftOn + 1
	}
	if (FetalData$V2[Id - 3] > FetalData$V2[Id - 4] || FetalData$V2[Id] > FetalData$V2[Id - 4]) {
		LeftOn = LeftOn + 1
	}
	if (FetalData$V2[Id - 4] > FetalData$V2[Id - 5] || FetalData$V2[Id] > FetalData$V2[Id - 5]) {
		LeftOn = LeftOn + 1
	}
	
	RightOn = 0
	if (FetalData$V2[Id] > FetalData$V2[Id + 1]) {
		RightOn = RightOn + 1
	}
	if (FetalData$V2[Id + 1] > FetalData$V2[Id + 2] || FetalData$V2[Id] > FetalData$V2[Id + 2]) {
		RightOn = RightOn + 1
	}
	if (FetalData$V2[Id + 2] > FetalData$V2[Id + 3] || FetalData$V2[Id] > FetalData$V2[Id + 3]) {
		RightOn = RightOn + 1
	}
	if (FetalData$V2[Id + 3] > FetalData$V2[Id + 4] || FetalData$V2[Id] > FetalData$V2[Id + 4]) {
		RightOn = RightOn + 1
	}
	if (FetalData$V2[Id + 4] > FetalData$V2[Id + 5] || FetalData$V2[Id] > FetalData$V2[Id + 5]) {
		RightOn = RightOn + 1
	}
	
	if (LeftOn >= 5 && RightOn >= 5) {
		lines(c(FetalData$V1[Id],FetalData$V1[Id]),c(0,FetalData$V2[Id]),col = "#696969",lty = 2, lwd = 1)
		text(x = FetalData$V1[Id],y = FetalData$V2[Id],FetalData$V1[Id],cex = 0.7,font = 2,adj = c(0.2,0))
	}
}

# 保存;
dev.off()