## ----tidy=TRUE-----------------------------------------------------------
rlib <- "../Rlib"
.libPaths(rlib)

## ----tidy=FALSE, size="small"--------------------------------------------
require( restrSiteUtils )
rseq <- c("ACTAGT", "CCTAGG", "GCTAGC")
makeRestrDataPackage( tempdir(), rseq, "hg18",
                     "BSgenome.Hsapiens.UCSC.hg18",
                     "6-CUTTER",
                     installLib=TRUE)

## ----read-data, tidy=FALSE, size="small"---------------------------------
sites.df <- read.csv("sites.df.csv")

sites.gr <-
    with(sites.df,
         GRanges(seqnames,
                 IRanges( start=start, width=2),
                 strand=strand, Sequence=Sequence,
                 BID=BID))

colnames(sites.df)


## ----get-dist, tidy=FALSE, size="small"----------------------------------
require( restrEnz.Hsapiens.UCSC.hg18.RENZ.6.CUTTER )
site.dist <- distRL( sites.gr, width.GR )

## ----tidy=FALSE, size="small"--------------------------------------------
sites.reverse <- sites.gr
strand( sites.reverse ) <-
    ifelse( strand( sites.gr ) == "+", "-", "+")
site.dist.reverse <- distRL( sites.reverse, width.GR )

## ----summarize, tidy=FALSE, size="small"---------------------------------
summary( site.dist )
summary( site.dist.reverse )

## ----mrcs, tidy=FALSE, size="small"--------------------------------------
mrc.gr <- sampleDist(3, site.dist, width.GR,
                     parentNames=sites.gr$Sequence)

mrc.dist <- distRL( mrc.gr, width.GR)

mrc.dist.matrix <-
    do.call(rbind, tapply(mrc.dist, mrc.gr$parentNames, c))

## ----tidy=FALSE, size="small"--------------------------------------------
 
## Be sure rows are in same order as sites.gr$Sequence or the
## following check will not work:

all(length(sites.gr$Sequence),
    sites.gr$Sequence==rownames(mrc.dist.matrix))

## Each MRC is found to have the same distance as its parent:

all(site.dist==mrc.dist.matrix)


## ----combos, tidy=FALSE, size="small"------------------------------------

mcols(mrc.gr) <- DataFrame(Sequence=mrc.gr$parentNames)

mcols(sites.gr)$BID <- NULL

combo.gr <- c(sites.gr,mrc.gr)

## use `type' to distinguish MRCs and real sites

mcols(combo.gr)$type <-
    rep(c("insertion","match"),c(length(sites.gr),length(mrc.gr)))


## ----get.flank, tidy=FALSE, size="small"---------------------------------

combo.flank.20 <-
    shift( flank( combo.gr, 10L, start=TRUE, both=TRUE),
          ifelse( strand( combo.gr ) == "+", 1L, -1L))

seq.flank.20 <-
      getSeq( Hsapiens, combo.flank.20, as.character=TRUE)

t( sapply( 1:20, function(x) {
    y <- substring( seq.flank.20, x, x)
    tab <- prop.table( table( y, combo.gr$type ), 2)
    round( log( tab[,1] / tab[,2], 2 ), 2 )}))


## ----results="asis", tidy=FALSE------------------------------------------
toLatex(sessionInfo())

