# Libraries:
library(plotrix)
library(igraph)
library(quantsmooth)

arc = function(x0, x1, y, xr, yr, col, lwd) {
  x = (x0 + x1)/2  # Center of arc
  xr = x - x0 	 # x-radius of arc
  
  x_points = seq(x - xr, x + xr, length.out = 200)
  y_points = y + yr * sqrt( 1  -  ( (x_points-x) / xr )^2 )
  
  lines(x_points,
        y_points,
        col = col,
        lwd = lwd
        )
}

# INPUTS:
# when running from command line:
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
#if (length(args)!=4) {
#  stop("Arguments: DATA_DIR(path), donor, normal, tumor", call.=FALSE)
#}
#path <- args[1] #input path
#donor <- args[2] #Sample name
#nSample <- args[3] #Normal Sample name
#tSample <- args[4] #Tumor Sample name
rescue = "yes"
saveClusteredMutations = "no"



###########################################
samples=c("I-H-106917-T2-1-D1-2","I-H-106917-T2-2-D1-2","I-H-106917-T2-3-D1-2","I-H-106917-T2-4-D1-2","I-H-130718-T1-1-D1-2","I-H-130718-T1-10-D1-2","I-H-130718-T1-11-D1-1","I-H-130718-T1-12-D1-1","I-H-130718-T1-2-D1-2","I-H-130718-T1-4-D1-2","I-H-130718-T1-6-D1-2","I-H-130718-T1-9-D1-2","I-H-130719-T1-2-D1-2","I-H-130719-T1-4-D1-2","I-H-130719-T1-5-D1-2","I-H-130719-T1-6-D1-2","I-H-130720-T1-2-D1-2","I-H-130720-T1-3-D1-2","I-H-130720-T1-4-D1-2","I-H-130720-T1-5-D1-2","I-H-130720-T1-8-D1-2","I-H-130720-T1-9-D1-1")
for(sample in samples) {
#sample="I-H-106917-T2-1-D1-2"
tSample=sample
cna=paste0("/Users/yellapav/Desktop/p220_2019/genome_plots/input/",sample,"_subclones.txt")
sv=paste0("/Users/yellapav/Desktop/p220_2019/genome_plots/input/",sample,"_sv.txt")
snv=paste0("/Users/yellapav/Desktop/p220_2019/genome_plots/input/",sample,"_vcf.txt")
purity=paste0("/Users/yellapav/Desktop/p220_2019/genome_plots/input/",sample,"_rho.txt")
#caveman needs CHROM and POSITION

# INFO PURITY SAMPLE:
info <- read.table(purity, sep="\t", header=T, stringsAsFactors = F)

# SV:
brass <- read.table(sv, sep="\t", header=T, stringsAsFactors = F)
brass$POS1=brass$START1
brass$POS2=brass$START2
tSample=as.character(head(brass$SAMPLE,n=1))
#######brass <- brass[brass$CUSTOM_FILTER != "OUT",]


# CNA (Battenberg curated output):
battenberg <- read.table(cna, sep = "\t", header = T, stringsAsFactors = F)
battenberg$nMaj2_A[is.na(battenberg$nMaj2_A)] = 0
battenberg$nMin2_A[is.na(battenberg$nMin2_A)] = 0
battenberg$frac2_A[is.na(battenberg$frac2_A)] = 0
battenberg$total = battenberg$nMaj1_A*battenberg$frac1_A + battenberg$nMin1_A*battenberg$frac1_A + battenberg$nMaj2_A*battenberg$frac2_A + battenberg$nMin2_A*battenberg$frac2_A
battenberg$minor = battenberg$nMin1_A*battenberg$frac1_A + battenberg$nMin2_A*battenberg$frac2_A

# SNV (combined of caveman, mutect2, strelka, muse):
caveman = read.table(snv, sep = "\t", header = T, stringsAsFactors = F)
caveman = caveman[,1:4]

chr_list = 1:24; names(chr_list) = c(1:22,"X","Y")
caveman$CHROM = chr_list[caveman$CHROM]
caveman = caveman[order(caveman$CHROM, caveman$POSITION),]
caveman$intermutdist = c(0, caveman$POSITION[2:nrow(caveman)] - caveman$POSITION[1:(nrow(caveman)-1)] + 1)
caveman$intermutdist[caveman$intermutdist < 0] = NA
caveman$CHROM[caveman$CHROM == 23] = "X"
caveman$CHROM[caveman$CHROM == 24] = "Y"

caveman$muttype = paste(caveman$REF, caveman$ALT, sep=">")
muttype_list = sort(unique(caveman$muttype))
colvec = c("lightpink1","darkolivegreen3","gray75","lightskyblue","black","red","red","black","lightskyblue","gray75","darkolivegreen3", "lightpink1")
names(colvec) = muttype_list

if(saveClusteredMutations == "yes"){
  clustermuts = which(intermutdist<1000)
  clustermuts = sort(unique(c(clustermuts,clustermuts+1)))
  write.table(caveman[clustermuts,], file="clustered_subs.5cols", col.names=T, row.names=F, sep="\t", quote=F)
  write.table(caveman[-clustermuts,], file="unclustered_subs.5cols", col.names=T, row.names=F, sep="\t", quote=F)
}


# CANCER CONSENSUS GENES::   
gene = read.table("~/Downloads/INFO/cancer_gene_census_hg19_COSMICv85.csv", sep=",", stringsAsFactors = F, header = T)
gene$chr <- unlist(lapply(gene$Genome.Location, function(x) strsplit(x, split = ":")[[1]][1]))
gene$start <- as.numeric(unlist(lapply(gene$Genome.Location, function(x) ifelse(strsplit(x, split = ":")[[1]][2] != "-", strsplit(strsplit(x, split = ":")[[1]][2], "-")[[1]][1], NA))))
gene$end <- as.numeric(unlist(lapply(gene$Genome.Location, function(x) ifelse(strsplit(x, split = ":")[[1]][2] != "-", strsplit(strsplit(x, split = ":")[[1]][2], "-")[[1]][2], NA))))
gene <- gene[!is.na(gene$start), ]

# CHROMOSOME LENGTHS:
chr_lens = read.table("~/Downloads/INFO/hg19.chrom_sizes", header=F, sep="\t", row.names=1, colClasses = c("character", "numeric"))
temp = rownames(chr_lens)
chr_lens = chr_lens[,1]
names(chr_lens) = temp
chr_lens <- chr_lens[names(chr_lens) %in% c(1:22, "X", "Y")]


# LINEAR GENOME PLOTS:
#pdf(paste0(path, donor, "/VC/", nSample, "_vs_", tSample, "/final/combine/svs/", nSample, "_vs_", tSample, "_genome_plot_lineal.pdf"), height=7, width=15)
pdf(paste0("/Users/yellapav/Desktop/p220_2019/genome_plots/plots/",sample,"_genome_plot_tripleSV.pdf"), height=7, width=15)
yrange = c(0, ceiling(max(battenberg$total))+1)
if(yrange[2] == 3){ yrange[2] <- yrange[2]+1 }
yrange_size = yrange[2] - yrange[1]
ylim = c(yrange[1], yrange[2] + 1.4*yrange_size)


# GLOBAL PROFILE:
layout(matrix(c(1,2),ncol=1), widths=c(7,7), heights=c(1,2), TRUE) 
par(mar=c(0,4,4,8)+0.1, xpd=F, cex.axis=0.75, cex.lab=0.75)

xlim = c(0, sum(chr_lens))

# SNV: 
# SNVS - Rainfall plot:
cumSumsChr <- cumsum(chr_lens)
caveman_global <- caveman
for(i in 2:22){
  caveman_global$POSITION[caveman_global$CHROM == i] <- caveman_global$POSITION[caveman_global$CHROM == i] + cumSumsChr[i-1]
}
caveman_global$POSITION[caveman_global$CHROM == "X"] <- caveman_global$POSITION[caveman_global$CHROM == "X"] + cumSumsChr[22]
caveman_global$POSITION[caveman_global$CHROM == "Y"] <- caveman_global$POSITION[caveman_global$CHROM == "Y"] + cumSumsChr["X"]


# Density plot: plot(density(caveman_global$POSITION, 100))
plot(1, type="n", xlab="", xlim=xlim, ylim=c(1, max(caveman$intermutdist,na.rm=T)), log="y", 
     main=paste0(tSample, " (Predicted purity: ", info$Predicted_purity, " | Predicted ploidy: ", info$Predicted_ploidy, "): Global profile"), 
     ylab="Intermutation distance", las=2, bty ="n", xaxt="n", yaxs = "i", xaxs = "i" )
par(new=TRUE, xpd=T)
points(x = caveman_global$POSITION, y = caveman_global$intermutdist, pch=20, cex=0.5, col=colvec[caveman_global$muttype])

legend("right", bty = "n", inset=c(-0.0975,0), ncol=2,
       legend=names(colvec), 
       pch=20, col=colvec, cex=0.7)

# CHR separation
segments(x0 = cumsum(chr_lens), x1=cumsum(chr_lens), y0=1, y1=max(caveman$intermutdist,na.rm=T), col="gray20", lwd = 1) #ceiling(max(ascat$tot))+4
mtext(c(1:22, "X", "Y"), side = 3, 0, at = cumsum(chr_lens)-(cumsum(chr_lens)-c(0, cumsum(chr_lens)[1:23]))/2, adj = 0.5, cex=1)


# CNA:
par(mar=c(4,4,0,8)+0.1, xpd=F)
plot(1, type="n", xlab="", xlim=xlim, ylim=ylim, ylab="", las=2,bty ="n", xaxt="n", yaxt="n", yaxs = "i", xaxs = "i")
par(new=TRUE, xpd=T)
rect(par("usr")[1],-.5,par("usr")[2],0.5,col = "gray95", border = FALSE, lwd = 0.001)

col<- rep(c("gray99","gray95"), yrange[2])
for(z in seq(1,yrange[2]-1,1)) { 
  rect(par("usr")[1],z-0.5,par("usr")[2],z+0.5,col = col[z], border = FALSE, lwd = 0.001)
}

segments(xlim[1], seq(yrange[1], yrange[2]-1,1), 
         xlim[2], seq(yrange[1],yrange[2]-1,1),
         lwd=0.5, col="gray75", lty=1)

axis(2, line=0, at = seq(0, yrange[2]-1, 1),  col.axis = "black",
     labels=seq(0, yrange[2]-1, 1), lwd.ticks =1, cex.axis = 0.75, las=2)
mtext("Copy Number", 2, +2.5, at = median(c(yrange[1], yrange[2]-1)), adj = 0.5, cex=0.75)


for(i in c(1:24)){
  if(i == 1){
    paintCytobands(i, pos = c(0, -0.55), units = "hg19", width = 0.2, cex.leg = 0.1, 
                   bands = "major", legend=F, length.out=chr_lens[i], bleach=0.3) 
  }
  else{
    if(i == 23){chrom = "X"}
    else if(i == 24){chrom = "Y"}
    else{ chrom = i}
    paintCytobands(chrom, pos = c(cumSumsChr[i-1], -0.55), units = "hg19", width = 0.2, cex.leg = 0.1, 
                   bands = "major", legend=F, length.out=chr_lens[i], bleach=0.3) 
  }
}


# BATTENBERG segments:
battenberGlobal <- battenberg#[battenberg$chr %in% 1:22, ]
for(i in 2:22){
  battenberGlobal$startpos[battenberGlobal$chr == i] = battenberGlobal$startpos[battenberGlobal$chr == i] + cumSumsChr[i-1]
  battenberGlobal$endpos[battenberGlobal$chr == i] = battenberGlobal$endpos[battenberGlobal$chr == i] + cumSumsChr[i-1]
}
battenberGlobal$startpos[battenberGlobal$chr == "X"] = battenberGlobal$startpos[battenberGlobal$chr == "X"] + cumSumsChr[22]
battenberGlobal$endpos[battenberGlobal$chr == "X"] = battenberGlobal$endpos[battenberGlobal$chr == "X"] + cumSumsChr[22]

battenberGlobal$startpos[battenberGlobal$chr == "Y"] = battenberGlobal$startpos[battenberGlobal$chr == "Y"] + cumSumsChr["X"]
battenberGlobal$endpos[battenberGlobal$chr == "Y"] = battenberGlobal$endpos[battenberGlobal$chr == "Y"] + cumSumsChr["X"]



segments(battenberGlobal$startpos, battenberGlobal$minor, 
         battenberGlobal$endpos, battenberGlobal$minor,
         lwd=3, col="darkslategray", lty=1)

segments(battenberGlobal$startpos, battenberGlobal$total, 
         battenberGlobal$endpos, battenberGlobal$total,
         lwd=3, col="goldenrod2", lty=1)


# CHR separation
segments(x0 = cumsum(chr_lens), x1=cumsum(chr_lens), y0=-1, y1=yrange[2], col="gray20", lwd = 1) #ceiling(max(ascat$tot))+4


# SVs:
# First dotted lines:
segments(
  x0 = c(1, 1), x1 = c(xlim[2], xlim[2]),
  y0 = c(yrange[2] + 0.3*yrange_size, yrange[2] + 0.75*yrange_size),
  lty = 3
)
segments(x0 = c(1, 1), x1 = c(xlim[2], xlim[2]),y0=yrange[2] + 0.5, lty=1) 

# BRASS output:
brassGlobal <- brass#[brass$chr1 %in% 1:22 & brass$chr2 %in% 1:22, ]
for(i in 2:22){
  brassGlobal$POS1[brassGlobal$CHR1 == i] = brassGlobal$POS1[brassGlobal$CHR1 == i] + cumSumsChr[i-1]
  brassGlobal$POS2[brassGlobal$CHR2 == i] = brassGlobal$POS2[brassGlobal$CHR2 == i] + cumSumsChr[i-1]
}
brassGlobal$POS1[brassGlobal$CHR1 == "X"] = brassGlobal$POS1[brassGlobal$CHR1 == "X"] + cumSumsChr[22]
brassGlobal$POS2[brassGlobal$CHR1 == "X"] = brassGlobal$POS2[brassGlobal$CHR1 == "X"] + cumSumsChr[22]

brassGlobal$POS1[brassGlobal$CHR1 == "Y"] = brassGlobal$POS1[brassGlobal$CHR1 == "Y"] + cumSumsChr["X"]
brassGlobal$POS2[brassGlobal$CHR1 == "Y"] = brassGlobal$POS2[brassGlobal$CHR1 == "Y"] + cumSumsChr["X"]


if(nrow(brassGlobal)>0) {
  for (j in (1:nrow(brassGlobal))) { 
    #for (j in (1:10)) {
    ltyBrass = 1 #ifelse(brassGlobal$assembly_score[j] != "_", 1, 2) # solid line brass2
    if(brassGlobal$CHR1[j] == brassGlobal$CHR2[j]){
      if(brassGlobal$STRAND1[j] == "+" & brassGlobal$STRAND2[j] == "+"){  # DELETION
        segments(brassGlobal$POS1[j], 0, brassGlobal$POS1[j], yrange[2] + 0.3*yrange_size, col="firebrick2", lwd=1, lty=ltyBrass)
        segments(brassGlobal$POS2[j], 0, brassGlobal$POS2[j], yrange[2] + 0.3*yrange_size, col="firebrick2", lwd=1, lty=ltyBrass)
        arc(x0  = brassGlobal$POS1[j],
            x1  = brassGlobal$POS2[j],
            y   = yrange[2]+.3*yrange_size, 
            yr  = ifelse(brassGlobal$STRAND2[j] == "+", 1, -1) * 0.2 * yrange_size,
            col = "firebrick2",
            lwd = 1
        )
      }
      else if(brassGlobal$STRAND1[j] == "-" & brassGlobal$STRAND2[j] == "-") { # tandem-duplication
        segments(brassGlobal$POS1[j], 0, brassGlobal$POS1[j], yrange[2] + 0.3*yrange_size, col="dodgerblue1", lwd=1, lty=ltyBrass)
        segments(brassGlobal$POS2[j], 0, brassGlobal$POS2[j], yrange[2] + 0.3*yrange_size, col="dodgerblue1", lwd=1, lty=ltyBrass)
        arc(
          x0  = brassGlobal$POS1[j],
          x1  = brassGlobal$POS2[j],
          y   = yrange[2]+.3*yrange_size, 
          yr  = ifelse(brassGlobal$STRAND2[j] == "+", 1, -1) * 0.2 * yrange_size,
          col = "dodgerblue1",
          lwd = 1
        )
      }
      else{ # Inversion
        inversionCol <- "palegreen3"
        segments(brassGlobal$POS1[j], 0, brassGlobal$POS1[j], yrange[2] + 0.75*yrange_size, col=inversionCol, lwd=1, lty=ltyBrass)
        segments(brassGlobal$POS2[j], 0, brassGlobal$POS2[j], yrange[2] + 0.75*yrange_size, col=inversionCol, lwd=1, lty=ltyBrass)
        arc(
          x0  = brassGlobal$POS1[j],
          x1  = brassGlobal$POS2[j],
          y   = yrange[2]+.75*yrange_size, 
          yr  = ifelse(brassGlobal$STRAND2[j] == "+", 1, -1) * 0.2 * yrange_size,
          col = inversionCol,
          lwd = 1
        )
      }
    }
    
    else { # inter-chromosomal translocation
      segments(brassGlobal$POS1[j], 0, brassGlobal$POS1[j], yrange[2] + 0.75*yrange_size, col="gray20", lwd=1, lty=ltyBrass)
      segments(brassGlobal$POS2[j], 0, brassGlobal$POS2[j], yrange[2] + 0.75*yrange_size, col="gray20", lwd=1, lty=ltyBrass)
      arc(
        x0  = brassGlobal$POS1[j],
        x1  = brassGlobal$POS2[j],
        y   = yrange[2]+.75*yrange_size, 
        yr  = ifelse(brassGlobal$STRAND2[j] == "+", 1, -1) * 0.2 * yrange_size,
        col = "gray20",
        lwd = 1
      )
    }
  }
}


# add legend:
legend("right", bty = "n", inset=c(-0.105,0),
       legend=c("Total CN", "Minor CN", 
                "Deletion", "Tandem-dup", "Inversion", "Translocation"), 
       lty=1, lwd=c(4,4, rep(1,5)), col=c("goldenrod2", "darkslategray", "firebrick2", "dodgerblue1",  "palegreen3", "gray20"), cex=0.75)




# CHROMS SEPARATED:
for(i in c(1:22, "X", "Y")) {
  
  xlim <- c(0, chr_lens[names(chr_lens) == as.character(i)])
  layout(matrix(c(1,2),ncol=1), widths=c(7,7), heights=c(1,2), TRUE) 
  par(mar=c(0,4,4,8)+0.1, xpd=F, cex.axis=0.75, cex.lab=0.75)
  
  # SNVS - Rainfall plot:
  caveman_chr = caveman[caveman$CHROM == i, ]
  plot(1, type="n", xlab="", xlim=xlim, ylim=c(1, max(caveman$intermutdist,na.rm=T)), log="y", 
       main=paste0(tSample, ": chromosome ", i), 
       ylab="Intermutation distance", las=2, bty ="n", xaxt="n", yaxs = "i", xaxs = "i" )
  par(new=TRUE, xpd=T)
  points(x = caveman_chr$POSITION, y = caveman_chr$intermutdist, pch=20, cex=0.5, col=colvec[caveman_chr$muttype])
  
  legend("right", bty = "n", inset=c(-0.0925,0), ncol=2,
         legend=names(colvec), 
         pch=20, col=colvec, cex=0.7)

  # CNA:
  par(mar=c(10,4,0,8)+0.1, xpd=F)
  battenberg_sample_filter_chr <- battenberg[battenberg$chr == i, ]
  
  plot(1, type="n", xlab="", xlim=xlim, ylim=ylim, ylab="", las=2,bty ="n", xaxt="n", yaxt="n", yaxs = "i", xaxs = "i")
  par(new=TRUE, xpd=T)
  
  rect(par("usr")[1],-.5,par("usr")[2],0.5,col = "gray95", border = FALSE, lwd = 0.001)
  
  col= rep(c("gray99","gray95"), yrange[2])
  for(z in seq(1,yrange[2]-1,1)) { 
    rect(par("usr")[1],z-0.5,par("usr")[2],z+0.5,col = col[z], border = FALSE, lwd = 0.001)
  }
  
  axis(1, pos = -.85, at = pretty(xlim)[pretty(xlim) <= xlim[2]], labels = pretty(xlim)[pretty(xlim) <= xlim[2]] / 1e6)
  #mtext(paste0("chr", i, " (Mb)"), 1, 1, at = median(yrange), adj = 1, cex=0.5)
  
  axis(2, line=0, at = seq(0, yrange[2]-1, 1),  col.axis = "black",
       labels=seq(0, yrange[2]-1, 1), lwd.ticks =1, cex.axis = 0.75, las=2)
  mtext("Copy Number", 2, +2.5, at = median(c(yrange[1], yrange[2]-1)), adj = 0.5, cex=0.75)
  
  
  paintCytobands(i, pos = c(0, -0.55), units = "hg19", width = 0.3, cex.leg = 0.1, 
                 bands = "major", legend=F, length.out = xlim[2])
  
  segments(xlim[1], seq(yrange[1], yrange[2]-1,1), 
           xlim[2], seq(yrange[1],yrange[2]-1,1),
           lwd=0.5, col="gray75", lty=1)
  
  
  # BATTENBERG segments:
  segments(battenberg_sample_filter_chr$startpos, battenberg_sample_filter_chr$minor, 
           battenberg_sample_filter_chr$endpos, battenberg_sample_filter_chr$minor,
           lwd=4, col="darkslategray", lty=1)
  
  segments(battenberg_sample_filter_chr$startpos, battenberg_sample_filter_chr$total, 
           battenberg_sample_filter_chr$endpos, battenberg_sample_filter_chr$total,
           lwd=4, col="goldenrod2", lty=1)
  
  
  # SVs:
  # First dotted lines:
  segments(
    x0 = c(1, 1), x1 = c(xlim[2], xlim[2]),
    y0 = c(yrange[2] + 0.3*yrange_size, yrange[2] + 0.75*yrange_size),
    lty = 3
  )
  segments(x0 = c(1, 1), x1 = c(xlim[2], xlim[2]),y0=yrange[2] + 0.5, lty=1)    
    
    
  # Colors used in the other script:
  #td_col         = "darkorange4" => head to tail > tandem duplication
  #del_col        = "darkslateblue"  => tail to head > deletion
  #inter_chrs_col = "darkorchid1"
  #tail_tail_col  = "cadetblue4" => inversion
  #head_head_col  = "chartreuse2" => inversion
    
  # BRASS output:
  brass_sample_chr <- brass[brass$CHR1 == i | brass$CHR2 ==i, ]  
  if(nrow(brass_sample_chr)>0) {
    for (j in (1:nrow(brass_sample_chr))) { 
      ltyBrass <- 1 # ifelse(brass_sample_chr$assembly_score[j] != "_", 1, 2) # solid line brass2
      if(brass_sample_chr$CHR1[j] == brass_sample_chr$CHR2[j]){
        if(brass_sample_chr$STRAND1[j] == "+" & brass_sample_chr$STRAND2[j] == "+"){  # DELETION
          segments(brass_sample_chr$POS1[j], 0, brass_sample_chr$POS1[j], yrange[2] + 0.3*yrange_size, col="firebrick2", lwd=1, lty=ltyBrass)
          segments(brass_sample_chr$POS2[j], 0, brass_sample_chr$POS2[j], yrange[2] + 0.3*yrange_size, col="firebrick2", lwd=1, lty=ltyBrass)
          arc(x0  = brass_sample_chr$POS1[j],
              x1  = brass_sample_chr$POS2[j],
              y   = yrange[2]+.3*yrange_size, 
              yr  = ifelse(brass_sample_chr$STRAND2[j] == "+", 1, -1) * 0.2 * yrange_size,
              col = "firebrick2",
              lwd = 1
          )
        }
        else if(brass_sample_chr$STRAND1[j] == "-" & brass_sample_chr$STRAND2[j] == "-") { # tandem-duplication
          segments(brass_sample_chr$POS1[j], 0, brass_sample_chr$POS1[j], yrange[2] + 0.3*yrange_size, col="dodgerblue1", lwd=1, lty=ltyBrass)
          segments(brass_sample_chr$POS2[j], 0, brass_sample_chr$POS2[j], yrange[2] + 0.3*yrange_size, col="dodgerblue1", lwd=1, lty=ltyBrass)
          arc(
            x0  = brass_sample_chr$POS1[j],
            x1  = brass_sample_chr$POS2[j],
            y   = yrange[2]+.3*yrange_size, 
            yr  = ifelse(brass_sample_chr$STRAND2[j] == "+", 1, -1) * 0.2 * yrange_size,
            col = "dodgerblue1",
            lwd = 1
          )
        }
        else{ # Inversion
          inversionCol <- "palegreen3"
          segments(brass_sample_chr$POS1[j], 0, brass_sample_chr$POS1[j], yrange[2] + 0.75*yrange_size, col=inversionCol, lwd=1, lty=ltyBrass)
          segments(brass_sample_chr$POS2[j], 0, brass_sample_chr$POS2[j], yrange[2] + 0.75*yrange_size, col=inversionCol, lwd=1, lty=ltyBrass)
          arc(
            x0  = brass_sample_chr$POS1[j],
            x1  = brass_sample_chr$POS2[j],
            y   = yrange[2]+.75*yrange_size, 
            yr  = ifelse(brass_sample_chr$STRAND2[j] == "+", 1, -1) * 0.2 * yrange_size,
            col = inversionCol,
            lwd = 1
          )
        }
      }
      
      else { # inter-chromosomal translocation
        #low end in chr
        if(brass_sample_chr$CHR1[j] == i) {
          segments(brass_sample_chr$POS1[j], 0, brass_sample_chr$POS1[j], yrange[2] + yrange_size, col="gray20", lwd=1, lty=ltyBrass)
          segments(
            x0 = brass_sample_chr$POS1[j],
            y0 = yrange[2] + yrange_size,
            x1 = brass_sample_chr$POS1[j] + ifelse(brass_sample_chr$STRAND1[j] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            y1 = ifelse(brass_sample_chr$STRAND1[j] == "+", yrange[2] + 1.1 * yrange_size, yrange[2] + 1.3 * yrange_size),
            xpd = NA, lty=ltyBrass
          )
          text(
            as.character(brass_sample_chr$CHR2[j]),
            x = brass_sample_chr$POS1[j] + ifelse(brass_sample_chr$STRAND1[j] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            y = yrange[2] + ifelse(brass_sample_chr$STRAND1[j] == "+", 1.2, 1.4) * yrange_size,
            xpd = NA
          )
        }
        # high end in chr
        else {
          segments(brass_sample_chr$POS2[j], 0, brass_sample_chr$POS2[j],
                   yrange[2] + yrange_size, col="gray20", lwd=1, lty=ltyBrass)
          segments(
            x0 = brass_sample_chr$POS2[j],
            y0 = yrange[2] + yrange_size,
            x1 =  brass_sample_chr$POS2[j] + ifelse(brass_sample_chr$STRAND2[j] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            y1 = ifelse(brass_sample_chr$STRAND2[j] == "+", yrange[2] + 1.1 * yrange_size, yrange[2] + 1.3 * yrange_size),
            xpd = NA,
            lty=ltyBrass
          )
          text(
            as.character(brass_sample_chr$CHR1[j]),
            x = brass_sample_chr$POS2[j] + ifelse(brass_sample_chr$STRAND2[j] == "+", -1, +1) * (par("usr")[2] - par("usr")[1])/50,
            y = yrange[2] + ifelse(brass_sample_chr$STRAND2[j] == "+", 1.2, 1.4) * yrange_size,
            xpd = NA
          )
        }
      }
    }
  }
  
  # add legend:
  legend("right", bty = "n", inset=c(-0.105,0),
         legend=c("Total CN", "Minor CN", 
                  "Deletion", "Tandem-dup", "Inversion", "Translocation"), 
         lty=1, lwd=c(4,4, rep(1,5)), col=c("goldenrod2", "darkslategray", "firebrick2", "dodgerblue1", "palegreen3", "gray20"), cex=0.75)
  
  
  # CANCER GENES:
  gene_chr= gene[gene$chr == i & gene$Hallmark == "Yes",] # intentar millorarho
  if(nrow(gene_chr) > 0){
    for(z in (1:nrow(gene_chr))) {
      text(gene_chr$start[z], -3.2*ylim[2]/10, gene_chr$Gene.Symbol[z], srt = 60, cex=1, col="brown3", adj=1)
      segments(gene_chr$start[z], -0.85, gene_chr$start[z], -3*ylim[2]/10, lty = 2)
      points(gene_chr$start[z], -3*ylim[2]/10, type = "p", col="black", pch=21, bg="red", cex=1)
    }
  }
}

dev.off()
}

