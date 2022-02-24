
source("~/MS/ForSubmissionFunctions.R")
load("../data/hCSE.rda")
load("../data/tcgaNTMatched.rda")
load("../data/tcgaStage.rda")
load("../data/tcgaGrade.rda")
load("../data/HPA2022Jan9.RData")
load("../data//TCGAentrezInsClass202111.RData")


# Figure 1A. Transcriptomic landscape of hCSE.

reshCSE <- list()
reshCSE$hCSE <- hCSE
reshCSE$tcgaNTMatched <- tcgaNTMatched

pwGenes <- reshCSE$hCSE$EntrezID[1:115]
TarDEG <- reshCSE$tcgaNTMatched
ref = T
no = 1 #
# TarDEG <- reshCSE$tcgaStageDEG; ref=F; no=2
# TarDEG <- reshCSE$tcgaGradeDEG; ref=F; no=3

names(TarDEG$deg.list)
name.inx0 <- names(TarDEG$pval.list[[1]])
tmp <- lapply(TarDEG$pval.list, function(x)
  x[name.inx0])
compare.pvals = data.frame(tmp)


name.inx <- names(TarDEG$logFC.list[[1]])
tmp2 <- lapply(TarDEG$logFC.list, function(x)
  x[name.inx])
compare.logFC = data.frame(tmp2)

tcgaStudyAbbreviations.pre$Abbreviation
vec.pre
if (nchar(names(compare.pvals)[1]) > 5) {
  names(compare.pvals) <- toTCGAAcronym(names(compare.pvals))
  names(compare.logFC) <- toTCGAAcronym(names(compare.logFC))
}
names(compare.logFC)

compare.pvals.tmp <- compare.pvals[pwGenes, ]
compare.pvals.tmp[compare.pvals.tmp > .05] <- 0
compare.pvals.tmp[compare.pvals.tmp > 0] <- 1
compare.logFC2 <- compare.logFC[pwGenes, ] * compare.pvals.tmp
rownames(compare.logFC2) <-
  mapg(rownames(compare.logFC2))
dim(compare.logFC2)

ignoreNonSig = T
mat = compare.logFC[pwGenes, ]
if (ignoreNonSig)
  mat[compare.pvals[pwGenes, ] > .05] <- 0
mat[mat > 1] = 1
mat[mat < (-1)] = -1


mat2 <- ifelse(compare.logFC[pwGenes, ] > 0, 1,-1)
mat2[compare.pvals[pwGenes, ] > .05] <- 0

if (ref) {
  tmp1 = sapply(as.data.frame(t(mat2)), function(x) {
    tmp = length(which(x == 0))
    ifelse(abs(sum(sign(x))) != (length(x) - tmp), sum(sign(x)), NA)
  }, simplify = T)
  tmp1.sort = sort(tmp1[!is.na(tmp1)], decreasing = T)
  
  tmp2 = sapply(as.data.frame(t(mat2)), function(x) {
    tmp = length(which(x == 0))
    ifelse(abs(sum(sign(x))) == (length(x) - tmp), sum(sign(x)), NA)
  }, simplify = T)
  tmp2.sort = sort(tmp2[!is.na(tmp2)], decreasing = T)
  rowOrder = c(sp1 <-
                 names(tmp2.sort[tmp2.sort > 0]),
               sp2 <-
                 names(tmp1.sort),
               sp3 <- names(tmp2.sort[tmp2.sort < 0]))
  mapg(sp1)
  mapg(sp3)
  row_split = c(rep(1, length(sp1)), rep(2, length(sp2)), rep(3, length(sp3)))
}

mat = mat[rowOrder, ]
rowOrder.ref = rowOrder
colOrder = c(
  "LIHC",
  "CHOL",
  "LUAD",
  "LUSC",
  "BLCA",
  "UCEC",
  "PRAD",
  "READ",
  "COAD",
  "HNSC",
  "BRCA",
  "THCA",
  "KIRP",
  "KIRC",
  "KICH"
)
mat.pval = compare.pvals[rowOrder, ]
rownames(mat) <-
  paste0(hCSE[rowOrder, ][, 3], " (", hCSE[rowOrder, ][, 2], ") ")
row_names.col <- hCSE[rowOrder, ][, 5]
row_names.col
column_names.col <- "black"

if (no == 2) {
  mat = mat[, c(6, 17, 9, 10, 1, 13, 8, 15, 2, 11, 5, 7, 16, 18, 14, 3, 4, 12)]
  cluster_columns = F
  r.tmp = .042
} else if (no == 3) {
  mat = mat[, c(3, 10, 6, 9, 4, 5, 1, 8, 7, 2)]
  cluster_columns = F
  r.tmp = .06
} else {
  mat = mat[, colOrder]
  cluster_columns = F
  r.tmp = .05
}

library(pheatmap)
pheatmap(mat, cluster_rows = F, cluster_cols = F, fontsize_row = 8)





# Figure 1B. Enrichment study of hCSE
library(org.Hs.eg.db)

# Compute FDR adjusted DEGs
TarDEG <- reshCSE$tcgaNTMatched
pwGenes <- rownames(reshCSE$hCSE)[1:115]
# pwGenes <- rownames(reshCSE$hCSE)[116:225]

name.inx0 <- names(TarDEG$pval.list[[1]])
tmp <- lapply(TarDEG$pval.list, function(x)
  x[name.inx0])
compare.pvals = data.frame(tmp)
dim(compare.pvals)

name.inx <- names(TarDEG$logFC.list[[1]])
tmp2 <- lapply(TarDEG$logFC.list, function(x)
  x[name.inx])
compare.logFC = data.frame(tmp2)

names(compare.pvals) <- toTCGAAcronym(names(compare.pvals))
names(compare.logFC) <- toTCGAAcronym(names(compare.logFC))

keys = keys(org.Hs.eg.db, keytype = "ENZYME")
enzymeGenes <-
  mapIds(
    org.Hs.eg.db,
    keys = keys,
    column = "ENTREZID",
    keytype = "ENZYME",
    multiVals = "list"
  )
enzymeGenes  <-
  unique(do.call("c", enzymeGenes))
length(enzymeGenes) # 2230

op <- par()
par(mfcol = c(2, 5))
seed <- getOption("mySeedV")
set.seed(seed)
nsim = 9999


ran.res.wg <- ran.res.enzyme <- list()
inx = intersect(rownames(TarDEG[[1]][[1]][[1]]), pwGenes)
temp.func.deg = function(inx) {
  compare.pvals.inx <- compare.pvals[inx,]
  
  compare.pvals.inx[compare.pvals.inx > .05] <- 0
  length(compare.pvals.inx[compare.pvals.inx > 0])
}

stat = temp.func.deg(inx)
for (j in 1:nsim) {
  inx.wg = intersect(rownames(TarDEG[[1]][[1]][[1]]), sample(rownames(TarDEG[[1]][[1]][[1]]), length(pwGenes)))
  inx.enzyme = intersect(rownames(TarDEG[[1]][[1]][[1]]), sample(enzymeGenes, length(pwGenes)))
  ran.res.wg[[j]] <- temp.func.deg(inx.wg)
  ran.res.enzyme[[j]] <- temp.func.deg(inx.enzyme)
}

for (k in c("ran.res.wg", "ran.res.enzyme")) {
  pval = (length(which(do.call(
    "c", eval(as.symbol(k))
  ) >= stat)) + 1) / (nsim + 1)
  pval
  hist(do.call("c", eval(as.symbol(k))))
  abline(v = stat, lty = 2, lwd = 0.3)
}


# Compute kappa & icc
ran.res.wg <- ran.res.enzyme <- list()
inx = intersect(rownames(TarDEG[[1]][[1]][[1]]), pwGenes)
temp.func.icc = function(inx) {
  compare.logFC.inx <-  compare.logFC[inx,]
  
  compare.logFC.inx1 <- compare.logFC.inx
  compare.logFC.inx1[compare.pvals[inx,] > .05] <- 0
  Kappa(compare.logFC.inx1)$value
  
}
stat = temp.func.icc(inx)

for (j in 1:nsim) {
  inx.wg = intersect(rownames(TarDEG[[1]][[1]][[1]]), sample(rownames(TarDEG[[1]][[1]][[1]]), length(pwGenes))) #
  inx.enzyme = intersect(rownames(TarDEG[[1]][[1]][[1]]), sample(enzymeGenes, length(pwGenes)))
  ran.res.wg[[j]] <- temp.func.icc(inx.wg)
  ran.res.enzyme[[j]] <- temp.func.icc(inx.enzyme)
}


for (k in c("ran.res.wg", "ran.res.enzyme")) {
  pval = (length(which(do.call(
    "c", eval(as.symbol(k))
  ) >= stat)) + 1) / (nsim + 1)
  pval
  hist(do.call("c", eval(as.symbol(k))))
  abline(v = stat, lty = 2, lwd = 0.3)
}


# Druggable(FDA)
HPA <-
  read.delim2('https://www.proteinatlas.org/search?format=tsv')
names(HPA)
tar = unique(.mapid(rownames(TarDEG[[1]][[1]][[1]])))
duplicated.inx
HPA.tar = HPA[which(HPA$Gene %in% tar), c('Gene', colnames(HPA)[grep('Pathology.prognostics', colnames(HPA))])]
HPA.tar <-
  HPA.tar[-(which(HPA.tar$Gene %in% HPA.tar$Gene[which(duplicated(HPA.tar$Gene))])), ]
rownames(HPA.tar) <-
  HPA.tar[, 1]
HPA.tar <- HPA.tar[,-1]
peep(HPA.tar)
HPA.tar.df = data.frame(lapply(HPA.tar, function(x) {
  if (length(which(x == "")) != 0)
    x[which(x == "")] <- 0
  x[grep("unprognostic \\(", x)] <- 0
  x[which(x != "0")] <-
    1
  x
}))
rownames(HPA.tar.df) = rownames(HPA.tar)
colnames(HPA.tar.df) <-
  gsub("Pathology.prognostics...", "", colnames(HPA.tar.df))
length(which(HPA.tar.df == "1"))
length(which(HPA.tar.df == "0"))
FDA.tar = HPA[grep("FDA approved drug targets", HPA[, 8]), "Gene"]

potent.tar = HPA[grep("Potential drug targets", HPA[, 8]), "Gene"]



bg.fdr = F
ran.res.wg <- ran.res.enzyme <- list()
FDA.tar.entrez = mapg(FDA.tar)
if (bg.fdr) {
  inx = intersect(FDA.tar.entrez, pwGenes)
  stat = temp.func.deg(inx)
} else {
  stat = length(intersect(FDA.tar.entrez, pwGenes))
}
for (j in 1:nsim) {
  if (bg.fdr) {
    inx.wg = intersect(FDA.tar.entrez, sample(rownames(TarDEG[[1]][[1]][[1]]), length(pwGenes)))
    inx.enzyme = intersect(FDA.tar.entrez, sample(enzymeGenes, length(pwGenes)))
    ran.res.wg[[j]] <- temp.func.deg(inx.wg)
    ran.res.enzyme[[j]] <- temp.func.deg(inx.enzyme)
  } else {
    ran.res.wg[[j]]  = length(intersect(FDA.tar.entrez, sample(
      rownames(TarDEG[[1]][[1]][[1]]), length(pwGenes)
    )))
    ran.res.enzyme[[j]] = length(intersect(FDA.tar.entrez, sample(enzymeGenes, length(pwGenes))))
  }
}
for (k in c("ran.res.wg", "ran.res.enzyme")) {
  pval = (length(which(do.call(
    "c", eval(as.symbol(k))
  ) >= stat)) + 1) / (nsim + 1)
  pval
  hist(do.call("c", eval(as.symbol(k))))
  abline(v = stat, lty = 2, lwd = 0.3)
}



# Druggable(potential)
bg.fdr = F
ran.res.wg <- ran.res.enzyme <- list()
potent.tar.entrez = mapg(potent.tar)
if (bg.fdr) {
  inx = intersect(potent.tar.entrez, pwGenes)
  stat = temp.func.deg(inx)
} else {
  stat = length(intersect(potent.tar.entrez, pwGenes))
}
for (j in 1:nsim) {
  if (bg.fdr) {
    inx.wg = intersect(potent.tar.entrez, sample(rownames(TarDEG[[1]][[1]][[1]]), length(pwGenes)))
    inx.enzyme = intersect(potent.tar.entrez, sample(enzymeGenes, length(pwGenes)))
    ran.res.wg[[j]] <- temp.func.deg(inx.wg)
    ran.res.enzyme[[j]] <- temp.func.deg(inx.enzyme)
  } else {
    ran.res.wg[[j]] = length(intersect(potent.tar.entrez, sample(
      rownames(TarDEG[[1]][[1]][[1]]), length(pwGenes)
    )))
    ran.res.enzyme[[j]] = length(intersect(potent.tar.entrez, sample(enzymeGenes, length(pwGenes))))
  }
}
for (k in c("ran.res.wg", "ran.res.enzyme")) {
  pval = (length(which(do.call(
    "c", eval(as.symbol(k))
  ) >= stat)) + 1) / (nsim + 1)
  pval
  hist(do.call("c", eval(as.symbol(k))))
  abline(v = stat, lty = 2, lwd = 0.3)
}

# Survival genes
pwGenes2 <- mapg(pwGenes)
enzymeGenes2 <- mapg(enzymeGenes)
bgGene2 = rownames(HPA.tar.df)

inx = intersect(bgGene2, pwGenes2)
temp.func.surv = function(inx) {
  if (length(inx) == 0)
    next
  length(which(HPA.tar.df[inx, ] == "1"))
}
stat <- temp.func.surv(inx)

ran.res.wg <- ran.res.enzyme <- list()
for (j in 1:nsim) {
  inx.wg = intersect(bgGene2, sample(bgGene2, length(pwGenes2)))
  inx.enzyme = intersect(bgGene2, sample(enzymeGenes2, length(pwGenes2)))
  ran.res.wg[[j]] <- temp.func.surv(inx.wg)
  ran.res.enzyme[[j]] <- temp.func.surv(inx.enzyme)
}
for (k in c("ran.res.wg", "ran.res.enzyme")) {
  pval = (length(which(do.call(
    "c", eval(as.symbol(k))
  ) >= stat)) + 1) / (nsim + 1)
  pval
  hist(do.call("c", eval(as.symbol(k))))
  abline(v = stat, lty = 2, lwd = 0.3)
}



# Figure 1C. Penalized regression - leave one cancer types out cross validation 
load("../data/tcgaNTMatched.rda")
expClList = tcgaNTMatched$expClList
names.tmp <- names(expClList) <-  toTCGAAcronym(names(expClList))
expClList <-
  lapply(expClList, function(x) {
    list(x = x$exp, y = as.vector(x$cl[, 1]))
  })
d.merged <- expClList[colOrder]

seed <- getOption("mySeedV")
set.seed(seed)
d.merged$LIHC$x <-
  d.merged$LIHC$x[c("5831", "6241",     "32", "1036", "84803", "3158", "34", "5105"),]
d.merged <- summ.gene(d.merged, filter = NULL, "gmatch")
key.l <- label.l <- rep(list(c("C", "E")), 15)
s <- .tarExtract(d.merged, key.l, label.l)
names(s)
devCohorts = rep(1, 15)
penalty.factor = rep(1, dim(s[[1]]$x)[1])

alphas.seq = seq(0.01, 1, length = 10)

alphas.seq = 1
res <-
  biPDSmlv2(
    s,
    CVal = "LOSOCV",
    GlobalOp = "none",
    devCohorts = devCohorts,
    verbosePlotTable = T,
    alphas.seq = alphas.seq,
    penalty.factor = penalty.factor
  )


res$coefDf
res$cv.plot.list[[3]]
res$cv.plot.list[[4]]
res$cv.plot.list[[6]]
res$table.cv
res$table.val
pdf("./res.pdf", width = 5, height = 4)
res$cv.plot.list
res$val.plot.list
dev.off()

autoplot(sscurves.cv.pre)

pw.in = list(hCSENonZeroCoef = c("5831", "6241", "32", "1036", "84803", "3158", "34", "5105"))
names(s)[c(3, 2, 6, 8, 9, 13)]

for (i in c(3, 2, 6, 8, 9, 13)) {
  PDS.res <-
    PDS(s[[i]]$x[, which(s[[i]]$y[, 1] == "C")], s[[i]]$x[, which(s[[i]]$y[, 1] ==
                                                                    "E")], pw.in)
  .inhousePCplot(
    PDS.res,
    names(PDS.res$curves),
    normals = which(s[[i]]$y[, 1] == "C"),
    PDF = F,
    preSetPar = F,
    plot.type = "2d",
    mainTitle = names(s)[i],
    cloud.cex = .3,
    ticktype = "detailed"
  )
}






# Figure 2. Multiomics interation
target.pwGenes.inx = 1:115
pwGenes = reshCSE$hCSE$EntrezID[target.pwGenes.inx]
pwGenesHugo = reshCSE$hCSE$pwGenesHugo[target.pwGenes.inx]
length(pwGenesHugo)
hCSE = reshCSE$hCSE[target.pwGenes.inx,]
colOrder = c(
  "LIHC",
  "CHOL",
  "LUAD",
  "LUSC",
  "BLCA",
  "UCEC",
  "PRAD",
  "READ",
  "COAD",
  "HNSC",
  "BRCA",
  "THCA",
  "KIRP",
  "KIRC",
  "KICH"
)
paste(colOrder, toTCGAAcronym(colOrder))

print(unique(unlist(strsplit(HPA[, c("Protein.class")], ', '))))

names(reshCSE$tcgaNTMatchedDEGFiltered$expClList) = toTCGAAcronym(names(reshCSE$tcgaNTMatchedDEGFiltered$expClList))

names(reshCSE$tcgaNTMatchedDEGFiltered$deg.list) = toTCGAAcronym(names(reshCSE$tcgaNTMatchedDEGFiltered$deg.list))


ix = colOrder[1]
ix
TT.tar = ix
TT.tar

if (!exists(paste0(TT.tar, "Ins"))) {
  load(sprintf("~/DB/Firehose/Ins/%sIns.RData", TT.tar))
}  # <<=======
TT <- eval(as.symbol(paste0(TT.tar, "Ins")))

reservoir.omics = list()
subtype.contrast = T
subtype.origin = "NATvsPT"

if (subtype.origin == "NATvsPT") {
  Ix = intersect(colnames(TT@x),
                 colnames(reshCSE$tcgaNTMatchedDEGFiltered$expClList[[TT.tar]]$exp))
  cl = data.frame(sub = ifelse(Ix %in% Ix[grep("-11", Ix)], "NAT", "PT"))
  rownames(cl) = Ix
} else {
  cl = data.frame(sub = cl[, 1])
  rownames(cl) = rownames(cl.pre)
  
}
cl
sub1 = levels(factor(cl[, "sub"]))[1]
sub1
sub2 = levels(factor(cl[, "sub"]))[2]
sub2
reservoir.omics$subtype = cl


require(maftools)
maf <-
  TT@meta$fireData@Mutation[which(TT@meta$fireData@Mutation$Hugo_Symbol %in% pwGenesHugo), ]
colnames(maf)
dim(maf)
maf$Matched_Norm_Sample_Barcode

if (subtype.contrast) {
  maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, 1, 15)
  maf.list = list()
  for (i in levels(factor(cl[, 1]))) {
    maf.list[[i]] = maf[which(maf$Tumor_Sample_Barcode %in% rownames(cl[which(cl[, 1] ==
                                                                                i), , drop = F])),]
  }
}
names(maf.list)
if (subtype.contrast) {
  maf.list = c(list(maf = maf), maf.list)
} else {
  maf.list = list(maf = maf)
}
maf.summary.list = list()
for (i in names(maf.list)) {
  maf2 = maf.list[[i]][, c("Hugo_Symbol", "Variant_Classification")]
  if (length(which(maf2$Variant_Classification %in% "Silent")) != 0)
    maf2 = maf2[-which(maf2$Variant_Classification %in% "Silent"), ]
  
  n.row = length(h <- unique(maf2$Hugo_Symbol))
  n.col = length(v <- unique(maf2$Variant_Classification))
  maf.summary = matrix(rep(0, n.row * n.col), n.row, n.col)
  colnames(maf.summary) = v
  rownames(maf.summary) = h
  
  maf.split = split(maf2$Variant_Classification, maf2$Hugo_Symbol)
  for (j in rownames(maf.summary)) {
    maf.summary[j, rownames(table(maf.split[[j]]))] <-
      table(maf.split[[j]])
  }
  maf.summary = t(maf.summary)
  maf.summary = maf.summary[names(sort(rowSums(maf.summary), decreasing = T)), , drop = F]
  maf.summary = maf.summary[, names(sort(colSums(maf.summary), decreasing = T)), drop = F]
  maf.summary
  dim(maf.summary)
  t(maf.summary)
  maf.summary.list[[i]] = maf.summary
}
maf.summary.list


maf.tar.subtype = levels(factor(cl[, 1]))[2]
if (subtype.contrast) {
  inx = intersect(colnames(maf.summary.list[[2]]), colnames(maf.summary.list[[3]]))
  if (length(inx) != 0) {
    res.maf = maf.summary.list[[maf.tar.subtype]]
    res.maf = res.maf[,-which(colnames(res.maf) %in% inx)]
  } else {
    res.maf = maf.summary.list[[maf.tar.subtype]]
  }
} else {
  res.maf = t(maf.summary.list[[1]])
}
res.maf

reservoir.omics$maf <- t(res.maf)
mafCol = adjustcolor(.col("mutation")(19), 1)
names(mafCol) = c(
  "Frame_Shift_Del",
  "Missense_Mutation",
  "In_Frame_Del",
  "Splice_Site",
  "Multi_Hit",
  "Frame_Shift_Ins",
  "Complex_Event",
  "In_Frame_Del",
  "Nonsense_Mutation",
  "IGR",
  "Silent",
  "RNA",
  "Intron",
  "Nonstop_Mutation",
  "ITD",
  "In_Frame_Ins",
  "Translation_Start_Site",
  "Amp",
  "Del"
)

# SCNA
p.adj.method = "BH"
tcga.mat <- TT@x
rownames(tcga.mat) <- .mapid(rownames(tcga.mat))
dim(tcga.mat)
tcga.mat[1:3, 1:3]
x1 <- tcga.mat

inx0 = intersect(rownames(x1), pwGenesHugo)
x1 <- x1[inx0,]
x2 <- TT@meta$fireData@GISTIC@ThresholdedByGene
dim(x2)
colnames(x2) <- gsub("\\.", "-", substr(colnames(x2), 1, 15))
rownames(x2) <- x2[, 1]
x2 <- x2[, -c(1:3)]
inx <-
  intersect(rownames(x1), rownames(x2))
inx2 <- intersect(colnames(x1), colnames(x2))
x1 <- x1[inx, inx2]
dim(x1)
x2 <- x2[inx, inx2]
dim(x2)
x2 <-
  as.matrix(sapply(x2, function(x)
    as.numeric(x)))
rownames(x2) <- inx



if (subtype.contrast) {
  inx3 = intersect(rownames(cl), colnames(x2))
  inx3
  x1 = x1[, inx3, drop = F]
  x2 = x2[, inx3, drop = F]
  tmp.cl <- cl[inx3, , drop = F]
  if (length(levels(factor(tmp.cl[, "sub"]))) != 1) {
    scna.chi.pval <- list()
    for (i in 1:dim(x2)[1]) {
      scna.tmp <-
        wilcox.test(x2[i, tmp.cl[, "sub"] == sub1], x2[i, tmp.cl[, "sub"] == sub2])
      scna.chi.pval[[i]] <- scna.tmp$p.value
    }
    CNV.stat = data.frame(Negative.log10.adj.pval = -log10(p.adjust(unlist(scna.chi.pval), method = p.adj.method)))
  } else {
    scna.oneSampleWilcox.pval <- list()
    for (i in 1:dim(x2)[1]) {
      scna.oneSampleWilcox.pval[[i]] <-
        wilcox.test(x2[i, ], mu = 0)$p.value
    }
    CNV.stat = data.frame(Negative.log10.adj.pval = -log10(p.adjust(
      unlist(scna.oneSampleWilcox.pval), method = p.adj.method
    )))
    rownames(CNV.stat) <- rownames(x2)
  }
  
  direction = rowMeans(x2[, tmp.cl$sub == sub2]) - rowMeans(x2[, tmp.cl$sub ==
                                                                 sub1])
  if (is.na(all(direction))) {
    CNV.stat$Direction = ifelse(rowMeans(x2) > 0, "Amplication", "Deletion")
  } else {
    CNV.stat$Direction = ifelse(direction > 0, "Amplication", "Deletion")
  }
} else {
  CNV.stat = NULL
}

results.cnv <- data.frame()
for (i in rownames(x1)) {
  results.cnv[i, "ES"] <-  cor.test(x1[i,],  x2[i,])$estimate[[1]]
  results.cnv[i, "Pval"] <-  cor.test(x1[i,],  x2[i,])$p.value
}
if (subtype.contrast) {
  if (!is.null(CNV.stat)) {
    results.cnv = cbind(results.cnv, CNV.stat)
  } else {
    results.cnv$Negative.log10.adj.pval = 0
    results.cnvl$Direction = 0
  }
  reservoir.omics$CNV <- results.cnv
}
results.cnv.sorted = sortDf(results.cnv, T, list(0, 0.05), list(">", "<"))

# Methylation ====
tcga.mat <- TT@x
rownames(tcga.mat) <- .mapid(rownames(tcga.mat))
dim(tcga.mat)
tcga.mat[1:3, 1:3]

x1 <- tcga.mat

inx0 = intersect(rownames(x1), pwGenesHugo)
x1 <- x1[inx0,]
x2 <- TT@meta$fireData@Methylation
dim(x2)
x2[1:3, 1:3]
colnames(x2)  <- gsub("\\.", "-", colnames(x2))
inx <-
  intersect(rownames(x1), rownames(x2))
inx2 <- intersect(colnames(x1), colnames(x2))
x1 <- x1[inx, inx2]
dim(x1)
x2 <- x2[inx, inx2]
dim(x2)
x2 <-
  as.matrix(sapply(x2, function(x)
    as.numeric(x)))
rownames(x2) <- inx
str(x1)
str(x2)

if (subtype.contrast) {
  inx3 = intersect(rownames(cl), colnames(x2))
  x1 = x1[, inx3, drop = F]
  x2 = x2[, inx3, drop = F]
  tmp.cl <- cl[inx3, , drop = F]
  dim(x1)
  dim(x2)
  dim(tmp.cl)
  if (length(levels(factor(tmp.cl[, "sub"]))) > 1) {
    methyl.stat <-
      data.frame(Negative.log10.adj.pval = -log10(
        uni.stat(x2, tmp.cl[, "sub", drop = F], stat = "ano", ano.method = "wilcox")$sub$stat[, 2]
      ))
    direction = rowMeans(x2[, tmp.cl$sub == sub2]) - rowMeans(x2[, tmp.cl$sub ==
                                                                   sub1])
    methyl.stat$Direction = ifelse(direction > 0, "Hyper", "Hypo")
    methyl.stat
  } else {
    methyl.stat <- NULL
  }
}

results.methyl <- data.frame()
for (i in rownames(x1)) {
  results.methyl[i, "ES"] <-
    cor.test(x1[i,],  x2[i,])$estimate[[1]]
  results.methyl[i, "Pval"] <-
    cor.test(x1[i,],  x2[i,])$p.value
}
if (subtype.contrast) {
  if (!is.null(methyl.stat)) {
    results.methyl = cbind(results.methyl, methyl.stat)
  } else {
    results.methyl$Negative.log10.adj.pval = 0
    results.methyl$Direction = 0
  }
  reservoir.omics$methyl <- results.methyl
}

results.methyl.sorted = sortDf(results.methyl, F, list(-0.2, 0.05), list("<", "<"))


# miRNA
require(MatrixEQTL)
exp <- TT@x
tar.miRNA <- TT@meta$fireData@miRNASeqGene
colnames(tar.miRNA) <- substr(colnames(tar.miRNA), 1, 15)
inx = intersect(colnames(exp), colnames(tar.miRNA))
x1 = exp[pwGenes, inx]
x2 = tar.miRNA[, inx]

x2 = sd.filter(x2, min.std = 0.2, zero.filter = round(dim(x2)[2] * .2, 0))
x3 = log2(x2 + 1)

expr <- SlicedData$new()
miRNA <- SlicedData$new()

expr$CreateFromMatrix(x1)
miRNA$CreateFromMatrix(x3)
show(expr)
show(miRNA)
res = Matrix_eQTL_engine(miRNA,
                         expr,
                         output_file_name = "",
                         pvOutputThreshold = 1)
res$all$eqtls$snps
res$all$eqtls$pvalue
res$all$eqtls$FDR
res1 = res$all$eqtls[, c(6, 5, 1:4)]
res1 = sortDf(res1, F, c(0, 0.05), c("<", "<"))

if (dim(res1)[1] > 200) {
  res1 = res1[1:200,]
  reservoir.omics$meta$miRNAQTLBetathreshold = res1[200, 1]
} else {
  reservoir.omics$meta$miRNAQTLBetathreshold = 0
}
reservoir.omics$miRNA <- res1

if (subtype.contrast) {
  inx3 = intersect(rownames(cl), colnames(x2))
  x1 = x1[, inx3, drop = F]
  x2 = x2[, inx3, drop = F]
  tmp.cl <- cl[inx3, , drop = F]
  dim(x1)
  dim(x2)
  dim(tmp.cl)
  miRNA.stat <-
    uni.stat(x2,
             tmp.cl[, "sub", drop = F],
             stat = "ano",
             ano.method = "moderateT",
             count = T)$sub$topTable[, c(1, 5)]
  try({
    tmp.chr = mapG(rownames(miRNA.stat))[[1]][, c(1:3, 8)]
  }, silent = T)
  try({
    miRNA.stat = cbind(miRNA.stat, tmp.chr[match(rownames(miRNA.stat), tmp.chr$Name), ])[, c(3:5, 1:2)]
  }, silent = T)
} else {
  miRNA.stat = NULL
}
reservoir.omics$miRNA.stat <- miRNA.stat

# mRNA
if (subtype.origin == "NATvsPT") {
  if (subtype.contrast) {
    res.stat = na.omit(sortDf(reshCSE$tcgaNTMatched$deg.list[[TT.tar]]$tT.filter[pwGenes, c(1, 5, 2:3)]))
    rownames(res.stat) <- mapg(rownames(res.stat))
    reservoir.omics$subtypeDE <- res.stat
  }
} else {
  if (subtype.contrast) {
    exp <- TT@x[pwGenes,]
    exp[1:3, 1:3]
    dim(exp)
    inx3 = intersect(rownames(cl), colnames(exp))
    inx3
    exp <-
      exp[, inx3, drop = F]
    tmp.cl <- cl[inx3, , drop = F]
    dim(tmp.cl)
    #  log2(2^exp)[1:3, 1:3]==exp[1:3, 1:3]
    res.stat <-
      uni.stat((2 ^ exp - 1),
               tmp.cl[, "sub", drop = F],
               plot = T,
               stat = "ano",
               ano.method = "moderateT",
               count = T
      )
    dat = res.stat$plots$data[res.stat$plots$data$adj.P.Val < sort(res.stat$plots$data$adj.P.Val, decreasing = F)[5],]
    res.stat$plots + geom_text_repel(
      data = dat,
      aes(label = hCSE[rownames(dat), "Pathway"]),
      show.legend = F,
      size = 3,
      arrow = arrow(length = unit(0.02, "npc"))
    )
    res.stat = res.stat[[1]]$topTable[, c(1, 5, 2:4, 6)][rownames(exp),]
    res.stat = sortDf(res.stat)
    rownames(res.stat) <- mapg(rownames(res.stat))
    reservoir.omics$subtypeDE <- res.stat
  }
}

reservoir.omics$meta$labels <-
  list(pwGenes = pwGenes,
       pwGenesHugo = pwGenesHugo,
       hCSE = hCSE)
mirEvidance = "predicted"
val0 = mapG(
  reservoir.omics$meta$labels[[1]],
  F,
  expand = 'validated.miRNA',
  mirEvidance = mirEvidance,
  miRNACutoffPercentage = 30
)
val = sapply(val0$mature_mirna_id, function(x)
  strsplit(x, " \\| "))
names(val) = val0$SYMBOL
length(as.character(unlist(val)))


# eQTmiRNA
res1 <- reservoir.omics$miRNA
dim(res1)
tmp.chr = mapG(res1$snps)[[1]][, c(1:3, 8)]
tmp.chr2 = cbind(res1, tmp.chr[match(res1$snps, tmp.chr$Name), ])
bed.gene = mapG(unique(res1$gene), F, NULL)[, c(5:7, 4)]

bed.gene$gene2 = rownames(bed.gene)
tmp.chr3 = cbind(tmp.chr2, bed.gene[match(tmp.chr2$gene, bed.gene$gene2), ])
tmp.chr3$gene == tmp.chr3$gene2

val2 <- val
names(val2) = mapg(names(val))
inx.intersect = intersect(unique(res1$gene), names(val2))
inx.intersect
if (length(inx.intersect) != 0) {
  val2 <- val2[inx.intersect]
  
  res3 <- res2 <- list()
  for (i in 1:length(val2)) {
    res2[[i]] <-
      intersect(res1[which(res1$gene %in% names(val2)[i]), "snps"], gsub("-3p|-5p", "", tolower(unlist(val2[[i]]))))
    res3[[i]]  <-
      res1[which(res1$gene %in% names(val2)[i]),][which(res1[which(res1$gene %in% names(val2)[i]), "snps"] %in% gsub("-3p|-5p", "", tolower(unlist(val2[[i]])))),]
  }
}
res2
res3
val3 <- do.call("rbind", res3)
tmp.chr4 = tmp.chr3[rownames(val3),]
valAndeQTmiRNA = na.omit(tmp.chr4[, c(1, 7:14)])
valAndeQTmiRNA


if (!is.null(reservoir.omics$miRNA.stat)) {
  miRNA.df.pre = reservoir.omics$miRNA.stat[reservoir.omics$miRNA.stat$adj.P.Val <
                                              0.05, , drop = F]
} else {
  miRNA.df.pre <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(miRNA.df.pre) <-
    c(
      "ChromosomeStart",
      "start",
      "End",
      "logFC",
      "adj.P.Val",
      "eqtl.beta",
      "forCircosLink"
    )
}
if (!is.null(reservoir.omics$miRNA.stat) &&
    dim(miRNA.df.pre)[1] != 0) {
  miRNA.df.pre$eqtl.beta = NA
  forCircosLink = intersect(rownames(miRNA.df.pre), reservoir.omics$miRNA$snps) # <<<<<<<<<<
  miRNA.df.pre$forCircosLink = ifelse(rownames(miRNA.df.pre) %in% forCircosLink, 1, 0)
}

if (dim(miRNA.df.pre)[1] != 0 &&  dim(valAndeQTmiRNA)[1] != 0) {
  tmp.val2 = valAndeQTmiRNA[!duplicated(valAndeQTmiRNA$Name) &
                              !(valAndeQTmiRNA$Name %in% forCircosLink),]
} else {
  tmp.val2 = valAndeQTmiRNA[!duplicated(valAndeQTmiRNA$Name),]
}

miRNA.eqtl = data.frame(matrix(rep(NA, dim(tmp.val2)[1] * dim(miRNA.df.pre)[2]), dim(tmp.val2)[1], dim(miRNA.df.pre)[2]))
colnames(miRNA.eqtl) = colnames(miRNA.df.pre)
rownames(miRNA.eqtl) = tmp.val2$Name
miRNA.eqtl[, c(1:3, dim(miRNA.df.pre)[2] - 1)] <-
  tmp.val2[, c(2:4, 1)] 
miRNA.df.pre = rbind(miRNA.df.pre, miRNA.eqtl)
tail(miRNA.df.pre)

valeQTmiRNADEmiRNA = valAndeQTmiRNA[which(valAndeQTmiRNA$Name %in% forCircosLink),]
valeQTmiRNADEmiRNA = valeQTmiRNADEmiRNA[which(miRNA.df.pre[as.character(valeQTmiRNADEmiRNA$Name),]$logFC *reservoir.omics$subtypeDE[as.character(valeQTmiRNADEmiRNA$SYMBOL),]$logFC < 0), ]



# Figure 4-6. Survival analysis - subtype identification
names(TCGA.pancan@pheno)
names(TCGA.pancan@meta$PancanClinical)
TCGA.pancan@meta$PancanClinical2[, 1:10]
names(TCGA.pancan@meta$TCGAClinical)   # <<========
tcgaCl.tmp <- TCGA.pancan@meta$PancanClinical2
rownames(tcgaCl.tmp) <-
  TCGA.pancan@meta$PancanClinical2$bcr_patient_barcode
inx = intersect(pwGenes, rownames(TCGA.pancan@x))

res <- list()
canInx = c("pancreatic adenocarcinoma",
           "adrenocortical cancer",
           "kidney clear cell carcinoma")
i = canInx[1]
i
cl <-
  cbind(TCGA.pancan@pheno[inx2 <-
                            rownames(TCGA.pancan@pheno)[intersect(
                              which(TCGA.pancan@pheno$sample_type %in% "Primary Tumor"),
                              which(TCGA.pancan@pheno$X_primary_disease %in% i)
                            )], names(TCGA.pancan@pheno)[c(3, 4, 19, 20)]] , data.frame(tcgaCl.tmp)[substr(inx2, 1, 12), c("pathologic_stage"), drop =
                                                                                                      F])
colnames(cl) <- c("time", "status",  "age", "sex", "stage")

exp = TCGA.pancan@x[inx, inx2]
dim(exp)
rownames(exp) <- mapg(rownames(exp))
dim(exp)
dim(cl)
exp = exp[-which(
  rownames(exp) %in% c(
    "GFPT1",
    "GCH1",
    "GCLC",
    "GCLM",
    "PYCR1",
    "PYCR2",
    "PYCR3",
    "RRM1",
    "RRM2",
    "RRM2B",
    "GALK2",
    "SPTLC1",
    "SPTLC3",
    "PLPP4",
    "LPIN1",
    "LPIN3",
    "PLPP5",
    "LPIN2",
    "MOGAT3",
    "GPAT3",
    "GPAT4",
    "TH"
  )
),]   # <<<<<<<<<<<<<<<

pam.exp = consensusP(exp, maxK = 6)

cl$sex <- ifelse(cl$sex == "", NA, cl$sex)
cl$stage <- factor(cl$stage)
cl$Subtype <- pam.exp$res[[pam.exp$optK]]$consensusClass
cl$Subtype <- factor(cl$Subtype)
levels(cl$Subtype) <-  paste0("CPS_", 1:length(levels(cl$Subtype)))
require(survival)
bb <-
  cbind(cl, t(exp))
bb
bb$sex <-
  as.factor(bb$sex)
bb$stage <- as.factor(bb$stage)
dim(bb)
fit <-
  survival::survfit(survival::Surv(time, status) ~ Subtype, data = bb)
survminer::ggsurvplot(
  fit,
  data = bb,
  pval = T,
  risk.table = "nrisk_cumevents",
  surv.median.line = "hv",
  tables.height = .25,
  conf.int = F,
  size = 0.5,
  pval.size = 4,
  risk.table.fontsize = 3.5
)


# Dichotomizing
cl$Subtype <-
  ifelse(cl$Subtype %in% c("CPS_2", "CPS_6"), "hCSE_CPS_2", "hCSE_CPS_1")
# cl$Subtype <- ifelse(cl$Subtype %in% c( "CPS_2"), "hCSE_CPS_2", "hCSE_CPS_1") # adrenocortical cancer
bb <-
  cbind(cl, t(exp))
bb
bb$sex <- as.factor(bb$sex)
bb$stage <- as.factor(bb$stage)
dim(bb)
rownames(bb)
fit <-
  survival::survfit(survival::Surv(time, status) ~ Subtype, data = bb)
surv_diff <-
  survival::survdiff(survival::Surv(time, status) ~ Subtype, data = bb)
survminer::ggsurvplot(
  fit,
  data = bb,
  pval = T,
  risk.table = "nrisk_cumevents",
  surv.median.line = "hv",
  tables.height = .25,
  conf.int = F,
  size = 0.5,
  pval.size = 4,
  risk.table.fontsize = 3.5
)
coxph(Surv(time, status) ~ ., bb[, c(1:6)])










