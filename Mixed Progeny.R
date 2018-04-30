library("VariantAnnotation") #load the package
library("vcfR")
library("memgene")

### Setw
setwd("/Users/kbuttons/Documents/Research/P01/Data/Core C/")


### read in file 
Cross_SNPs <- read.vcfR("Cross.HQ.filt.recal.vcf")
colnames(Cross_SNPs@gt)
Cross_SNPs@gt[,c("Sample_ART-1","Sample_NF54_GFP_LUC")]

### Put these two columns into vectors 
Parent1_GT <- rep(-1, length(Cross_SNPs@gt[,"Sample_ART-1"]))
Parent2_GT <- rep(-1, length(Cross_SNPs@gt[,"Sample_NF54_GFP_LUC"]))

### loop through the entire vector and insert the 0/0 or 0/1 into vectors
for(i in 1:length(Cross_SNPs@gt[,"Sample_ART-1"])){
  Parent1_GT[i] <- unlist(strsplit(Cross_SNPs@gt[i,"Sample_ART-1"],":"))[1]
  Parent2_GT[i] <- unlist(strsplit(Cross_SNPs@gt[i,"Sample_NF54_GFP_LUC"],":"))[1]
}

### SNPs where Parents 1 and 2 are homozygotes	
Parental_SNPs <- which((Parent1_GT=="0/0"|Parent1_GT=="1/1"|Parent1_GT=="0/1"|Parent1_GT=="1/0")&(Parent2_GT=="0/0"|Parent2_GT=="1/1"|Parent2_GT=="1/0"|Parent2_GT=="0/1"))


### Store parental genotypes at SNPs where parents meet this above condition
SNPmatrix <- matrix(-1,ncol=(length(colnames(Cross_SNPs@gt))+6),nrow=length(Cross_SNPs@gt[,"Sample_ART-1"]))
for(i in 1:length(Cross_SNPs@gt[,"Sample_ART-1"])){
  SNPmatrix[i,1:7] <- Cross_SNPs@fix[i,1:7]
  for(j in 2:length(colnames(Cross_SNPs@gt))){
    SNPmatrix[i,j+6] <- unlist(strsplit(Cross_SNPs@gt[i,j],"[:]"))[1]
  }
}

colnames(SNPmatrix) <- c(colnames(Cross_SNPs@fix[,1:7]),colnames(Cross_SNPs@gt[,-1]))

### subset SNPs based on homozygosity of parents
Parental_homozySNPs <- SNPmatrix[Parental_SNPs,]

### Sort SNPs based on location
Parental_homozySNPs <- Parental_homozySNPs[order(Parental_homozySNPs[,"CHROM"],as.numeric(Parental_homozySNPs[,"POS"])),]

### Count heterozygous calls in matrix
het_counts <- matrix(0,nrow=length(colnames(Cross_SNPs@gt[,-1])),ncol=1)
for(j in 2:length(colnames(Cross_SNPs@gt))){
  het_counts[j-1,1] <- length(which(Parental_homozySNPs[,j+6]=="0/1"|Parental_homozySNPs[,j+6]=="1/0"))/length(which(Parental_homozySNPs[,j+6]=="0/1"|Parental_homozySNPs[,j+6]=="1/0"|Parental_homozySNPs[,j+6]=="1/1"|Parental_homozySNPs[,j+6]=="0/0"))
}

rownames(het_counts) <- colnames(Cross_SNPs@gt[,-1])

### histogram of heterozygous call counts
hist(het_counts,breaks=100,xlab="Proportion of heterozygous SNPs per progeny",ylab="Counts")
rownames(het_counts)[which(het_counts>700)]

### Sliding window size of average cross over to look for above average number of mutations
### From histogram, take 130 as our average heterozygous error rate for this 
window_size = 50
sliding_window_hetz <- matrix(0,nrow=dim(Parental_homozySNPs)[1]/window_size,ncol=length(colnames(Cross_SNPs@gt[,-1])))

for(i in 1:dim(sliding_window_hetz)[1]){
  for(j in 2:length(colnames(Cross_SNPs@gt))){
    sliding_window_hetz[i,j-1] <- length(which(Parental_homozySNPs[(((i-1)*window_size)+1):(i*window_size),j+6]=="0/1"|Parental_homozySNPs[(((i-1)*window_size)+1):(i*window_size),j+6]=="1/0"))/window_size
  }
}

colnames(sliding_window_hetz) <- colnames(Cross_SNPs@gt[,-1])

hist(sliding_window_hetz,breaks=30)
boxplot(sliding_window_hetz,ylab="Number of SNPs in 30 SNP sliding window",xlab="Samples")

### HeatMap of genome scan

library(RColorBrewer)
library(gplots)
library(vegan)
library(grDevices)

## Make vector of colors for values below threshold
rc1 <- colorpanel(8,"yellow", "yellow")
## Make vector of colors for values above threshold
rc2 <- colorpanel(8,"lightblue", "darkblue")
rampcols <- c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
#rampcols[c(-nHalf, nHalf+1)] <- rgb(t(col2rgb("green")), maxColorValue=256) 
rb1 <- seq(0, (2*1/window_size), length.out=8)
rb2 <- seq((2*1/window_size), 1, length.out=8)[-1]
rampbreaks <- c(-4,rb1,rb2,4)


mat <- sliding_window_hetz

png("HetScanMap.png", height=20, width=24, units="in", res=220)
par(oma=c(0,0,0,0))

heatmap.2(mat, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
          trace="none", density="none", scale="none",col = rampcols, 
          breaks = rampbreaks, cexRow=1,cexCol=1,na.rm=TRUE,na.color="grey90",
          key=FALSE,margins=c(5, 15), srtCol=360, labCol=NULL, labRow=NULL)

legend(0,1,c("background error","above background heterzygousity","high heterozygosity"),col=c("yellow","lightblue","darkblue"),pch=c(15,15,15),bty="n",y.intersp=1,cex=1.5)

dev.off()


### Check which columns have the maximum heterozygous counts for a window above expected
max_het <- rep(0,length(colnames(Cross_SNPs@gt[,-1])))
for(i in 1:length(colnames(Cross_SNPs@gt[,-1]))){
  max_het[i] <- max(sliding_window_hetz[,i])
}

colnames(sliding_window_hetz)[which(max_het>15)]

### Sliding window size of average cross over based on BP location of SNP
### From histogram, take 130 as our average heterozygous error rate for this 
window_size <- 8*30000   ### from Su et al Science 1999, 15 
tot_genome <- 23000000
sliding_window_hetz <- matrix(0,nrow=(tot_genome/window_size),ncol=length(colnames(Cross_SNPs@gt[,-1])))

### Find SNPs that most closely match window size by chromosome
chrom <- unique(Parental_homozySNPs[,"CHROM"])
chrom_length <- c(643292,947102,1060087,1204112,1343552,1418244,1501717,1419563,1541723,1687655,2038337,2271478,2895605,3281971)  ### lengths for 3D7 pulled from http://www.genome.jp/kegg-bin/show_organism?org=pfa
a <- c(1)
chrom_count <- 1

while(chrom_count <= 14){
  if(a[length(a)] + window_size < chrom_length[chrom_count]){
    a[[length(a)+1]] <- a[length(a)] + window_size
    print("incrementing a within chrom")
  }
  if(a[length(a)] + window_size >= chrom_length[chrom_count]){
    a[[length(a)+1]] <- chrom_length[chrom_count]
    a[[length(a)+1]] <- 1
    chrom_count <- chrom_count + 1
    print("incrementing chrom")
  }
}


for(i in 1:dim(sliding_window_hetz)[1]){
  for(j in 2:length(colnames(Cross_SNPs@gt))){
    sliding_window_hetz[i,j-1] <- length(which(Parental_homozySNPs[(((i-1)*window_size)+1):(i*window_size),j+6]=="0/1"|Parental_homozySNPs[(((i-1)*window_size)+1):(i*window_size),j+6]=="1/0"))
  }
}

colnames(sliding_window_hetz) <- colnames(Cross_SNPs@gt[,-1])

hist(sliding_window_hetz,breaks=30)
boxplot(sliding_window_hetz,ylab="Number of SNPs in 30 SNP sliding window",xlab="Samples")


### Check which columns have the maximum heterozygous counts for a window above expected
max_het <- rep(0,length(colnames(Cross_SNPs@gt[,-1])))
for(i in 1:length(colnames(Cross_SNPs@gt[,-1]))){
  max_het[i] <- max(sliding_window_hetz[,i])
}

colnames(sliding_window_hetz)[which(max_het>15)]

