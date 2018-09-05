#facets 后续 segment 和 拷贝数， ploidy， purity 在R 里完成。  可以衔接上，代码如下：
#以下为  读取18个人的 snp file， 生成 pdf图片和 segment 文件， segment 返回了 log2 ratio， 绝对拷贝数目等。 具体 每个segment 加基因 可以 用你们之前的代码。   
args<-commandArgs(T)
#library(argparser)
#p <- add_argument(p, "--od", help='输出目录')
#p <- add_argument(p, "--path",help='输入文件路径')
#p <- add_argument(p, "--name", help='输入文件名称')
#-----------传递参数
#argv <- parse_args(p)
library(facets)
path<-args[1]
purity<-matrix(data = NA, nrow = 1, ncol = 6, byrow = FALSE,
               dimnames = NULL)
rcmat<-readSnpMatrix(paste(args[1],args[2],".snp",sep=""))
xx<-preProcSample(rcmat,cval=100,ndepth=15,hetscale=FALSE,het.thresh=0.2,ndepthmax=1500)
oo<-procSample(xx,cval=100)
fit=emcncf(oo)
purity[1,1]<-args[2]
purity[1,2]<-fit$purity
purity[1,3]<-fit$ploidy
purity[1,4]<-oo$dipLogR
purity[1,5]<-toString(oo$flags[1])
purity[1,6]<-toString(oo$flags[2])
file<-paste(args[3], args[2],".pdf",sep="")
pdf(file)
plotSample(x=oo,emfit=fit)
dev.off()
#file <- paste(args[3], args[2],".seg",sep="")
file<-gsub(".pdf",".seg",file)
write.csv(fit$cncf,file,quote=F,row.names = F)


#write.table(purity,file="purity.txt",quote=F,row.names=F,sep="\t",append=T,col.names=F)

