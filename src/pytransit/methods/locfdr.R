require(locfdr)
require(GGMridge)

args = commandArgs(trailingOnly=TRUE)


data = read.table(args[1],sep='\t',head=T)

#lim=10
z = data$Z.score.of.concentration.dependence
#z[z > lim] = lim
#z[z < -lim] = -lim


NT = 2
mod = locfdr(z,nulltype=NT, df=20)
data$locfdr = mod$fdr

temp = data
temp = temp[order(temp$locfdr),]

EFDR = NULL
for (i in 1:nrow(temp))
{
  efdr = mean(temp$locfdr[1:i])
  EFDR = rbind(EFDR,efdr)
}
temp$efdr = EFDR
data$efdr = temp[rownames(data),"efdr"]

#cat(sprintf("RESULT: %s %d %d %d", args[1],length(data[(data$efdr < 0.05)& (z < 0),1]), length(data[(data$efdr < 0.05)& (z > 0),1]),length(data[(data$efdr < 0.05),1])))



write.table(data,args[1],sep="\t",row.names=F,quote=F)