plot_MA = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot", pch=20) {

    plot(logCounts, logFoldChange, col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);

}


plot_Volcano = function(logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Volcano plot", pch=20) {

   plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);

}


plot_MA_and_Volcano = function(logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot") {

    def.par = par(no.readonly = TRUE) # save default, for resetting...

    gridlayout = matrix(c(1:4),nrow=2,ncol=2, byrow=TRUE);
    layout(gridlayout, widths=c(1,1,1,1), heights=c(1,1,1,1)) 

    plot_MA(logCounts, logFoldChange, FDR);
    plot_Volcano(logFoldChange, FDR);

    # draw again, but use a smaller dot for data points
    plot_MA(logCounts, logFoldChange, FDR, pch='.');
    plot_Volcano(logFoldChange, FDR, pch='.');
    

    par(def.par)   
        
    
}

## GO
run_GO_enrichment <- function(path, contrast, functional_info){
  setwd(paste0(path,  '/intermediate/DE_F_culorum/Gly/', contrast))
  Analysed_genes <- read.table(paste0(path, '/intermediate/DE_F_culorum/Gly/', contrast, '/logFC-CPM-FDR.txt'), header = T, sep ='\t')
  Analysed_genes <- merge(Analysed_genes, functional_info, by = 'gene_id')
  Analysed_genes <- Analysed_genes[!is.na(Analysed_genes$GO_terms),]
  write.table(Analysed_genes[,c('gene_id')], paste0(path, '/intermediate/DE_F_culorum/Gly/', contrast, '/logFC-CPM-FDR_GO_Terms.txt'), col.names = F, sep ='\t', quote = F, row.names = F)
  
  #select DE egnes that have GO terms
  DE_edgeR_results <- read.table(paste0(path, '/intermediate/DE_F_culorum/Gly/', contrast, '/DE-edgeR-results.txt'), header = T, sep ='\t')
  DE_edgeR_results <- merge(DE_edgeR_results, functional_info, by = 'gene_id')
  DE_edgeR_results <- DE_edgeR_results[!is.na(DE_edgeR_results$GO_terms),]
  if (nrow(DE_edgeR_results)>=1){
  write.table(DE_edgeR_results[,c('gene_id')], paste0(path, '/intermediate/DE_F_culorum/Gly/', contrast, '/DE-edgeR-results_GO_Terms.txt'), col.names = F, sep ='\t', quote = F, row.names = F)
  #Run goatools
  writeLines(text = 'ln -s ~/SMS_5783_21_F_culmorum_RNA_Seq_annotation/data/meta_data/annotation_F_culorum/annotation_gene_id_GO.txt
ln -s /Users/nimra236/git/goatools/go-basic.obo
ln -s /Users/nimra236/git/goatools/goslim_generic.obo
#conda activate goatools
~/git/goatools/scripts/find_enrichment.py DE-edgeR-results_GO_Terms.txt logFC-CPM-FDR_GO_Terms.txt annotation_gene_id_GO.txt --outfile=goea_results.xlsx,goea_results.tsv --pvalcalc fisher --pval 0.05', con = 'goatools.sh') 
  cat(paste0('cd ', path,  '/intermediate/DE_F_culorum/Gly/', contrast,'\n', 'sh goatools.sh'))
  }else{
    print('There are no genes with GO term to perform GO enrichment analysis')
  }
}

##KEGG
library(KEGGprofile)
library(clusterProfiler)
run_KEGG <- function(DE_edgeR_results){
  DE_edgeR_results <- DE_edgeR_results[order(DE_edgeR_results$logFC, decreasing = T),]
  rownames(DE_edgeR_results) <- DE_edgeR_results$gene_id
  kegg = bitr_kegg(geneID = unique(DE_edgeR_results$Entrez_ID), fromType = "ncbi-geneid",
                   toType = c("kegg"),
                   organism = 'fox')
  kk <- enrichKEGG(gene         = kegg[,2],
                   organism     = 'fox',
                   pvalueCutoff = 0.05)
  ordered_DE_edgeR_results_logFC <- DE_edgeR_results$logFC
  names(ordered_DE_edgeR_results_logFC) <- DE_edgeR_results$gene_id
  # I removed permutation from gseKEGG                 nPerm        = 1000,
  kk2 <- gseKEGG(geneList     = ordered_DE_edgeR_results_logFC,
                 organism     = 'fox',
                 minGSSize    = 2,
                 pvalueCutoff = 1,
                 verbose      = FALSE)
  kk_list <- list(kk,kk2)
  names(kk_list) <- c('Overrepresentation', 'Set_enrichment_analyses')
  return(kk_list)
}
##heatmap.3
### pulled from here, and then tweaked slightly: http://www.biostars.org/p/18211/
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.1,
                      #cexRow = 0.2 + 1/log10(max(nr,2)),
                      #cexCol = 0.2 + 1/log10(max(nc,2)),
        cexRow = 0.2,
        cexCol = 0.2,                  

        scaleRangeMin,
        scaleRangeMax,


    cex.main = 1,
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
 
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }



    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }

    retval <- list()


    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)

    dendrogram <- match.arg(dendrogram)

    trace <- match.arg(trace)

    density.info <- match.arg(density.info)

    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")

    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")

    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE

    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE

    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")

    nr <- di[1]
    nc <- di[2]

    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    #print(paste("nr:", nr, "nc:", nc, "cexCol:", cexCol, "cexRow:", cexRow))
    #stop("debug")



    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")

    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))

    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"

            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }

    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"

            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
 
   if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
 
   if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }

    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()

    x <- x[rowInd, colInd]  # rearrange matrix according to dendrograms
    x.unscaled <- x

    cellnote <- cellnote[rowInd, colInd]  # also rearrange the cellnotes

    # get labels 
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]


    ## do scaling of matrix according to Z-scores
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }

    # number of breaks
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }

    # set breakpoints
    if (length(breaks) == 1) {
        if (missing(scaleRangeMin))
            scaleRangeMin = min(x, na.rm=na.rm)

        if (missing(scaleRangeMax))
            scaleRangeMax = max(x, na.rm=na.rm)


        if (!symbreaks) {
            breaks <- seq(scaleRangeMin, scaleRangeMax, length=breaks);
        } else {
            #extreme <- max(abs(x), na.rm = TRUE)
            extreme = max(abs(c(scaleRangeMin,scaleRangeMax)), na.rm=na.rm)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }

    nbr <- length(breaks)
    ncol <- length(breaks) - 1

    if (class(col) == "function")
        col <- col(ncol)

    min.breaks <- min(breaks)
    max.breaks <- max(breaks)

    # adjust for out-of-range given break settings
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks

    # layout height
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)

    # layout width
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)

    # define the layout
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
 
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || ncol(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of ncol(x) ", nc, " columns")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
            side_height = min(side.height.fraction*nrow(ColSideColors), 1);
            lhei=c(lhei[1], side_height, lhei[2])
        }
 
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || nrow(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of nrow(x) ", nr, " rows.  It currently has ", nrow(RowSideColors), " rows.")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
            side_width = min(side.height.fraction*ncol(RowSideColors), 1);
			lwid <- c(lwid[1], side_width, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
 
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    

	op <- par(no.readonly = TRUE)
    on.exit(par(op))
 
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
 
	###########################################
	## Draw the colorbars for the annotations:
	###########################################	

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[rowInd, , drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            # print(rsc)
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            #print("RSC: ", rsc)
            #print(rsc.colors)    
            image(1:nrow(rsc), 1:ncol(rsc), rsc, col = as.vector(rsc.colors), axes = FALSE, xlab="", ylab="")
		
			# add labels
            if (length(colnames(RowSideColors)) > 0) {  
                #axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
				#axis(1, 0:(nrow(rsc)-1), colnames(RowSideColors), las = 2, tick = T) # ncol because transposed
            	axis(1, 1:ncol(RowSideColors), labels=colnames(RowSideColors), las=2, cex.axis=0.5, tick=F, xlab="", ylab="")

			}
        }
    }
    


    if (!missing(ColSideColors)) {
 
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[, colInd, drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            #print(csc)
            image(1:nrow(t(csc)), 1:ncol(t(csc)), t(csc), col = as.vector(csc.colors), axes = FALSE, xlab="", ylab="")

			# add labels
            if (length(rownames(ColSideColors)) > 0) {
                #axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
				axis(2, 1:(nrow(ColSideColors)), labels=rownames(ColSideColors), las = 2, tick = FALSE, cex.axis=0.5)
            }
        }
    }
 


    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    
	# draw the central heatmap
	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
	
	# store the matrix drawn
	retval$carpet <- x
    
	# store the dendrograms
	if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    
	# store the breaks
	retval$breaks <- breaks
    
	# store the colormap used
	retval$col <- col
	
	# specially color in the na values	
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", col = na.color, add = TRUE)
    }

	# X-axis column labels
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)

	# X-axis title
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)

	# Y-axis row labeling
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)

	# Y-axis title
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
 
   	if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    

	min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    
	# column trace
	if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }

	# row trace
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }

	# add cell labels
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), col = notecol, cex = notecex)

	###########################
	## Plot the row dendrogram
	###########################

    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()

	#############################
	## Plot the column dendrogram
	#############################

    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()

    if (!is.null(main))
        title(main, cex.main=cex.main) #cex.main = 1.5 * op[["cex.main"]])


	############################
	## Add the Color Chart
	############################

    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(c(x,breaks), na.rm = TRUE)
            max.raw <- max(c(x,breaks), na.rm = TRUE)
        }
 
        message('for plotting:: min.raw: ', min.raw, ' max.raw: ', max.raw);
        
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()

    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], high = retval$breaks[-1], color = retval$col)

    invisible(retval)
}



colorpanel = function (n, low, mid, high) 
{
    if (missing(mid) || missing(high)) {
        low <- col2rgb(low)
        if (missing(high)) 
            high <- col2rgb(mid)
        else high <- col2rgb(high)
        red <- seq(low[1, 1], high[1, 1], length = n)/255
        green <- seq(low[3, 1], high[3, 1], length = n)/255
        blue <- seq(low[2, 1], high[2, 1], length = n)/255
    }
    else {
        isodd <- odd(n)
        if (isodd) {
            n <- n + 1
        }
        low <- col2rgb(low)
        mid <- col2rgb(mid)
        high <- col2rgb(high)
        lower <- floor(n/2)
        upper <- n - lower
        red <- c(seq(low[1, 1], mid[1, 1], length = lower), seq(mid[1, 
            1], high[1, 1], length = upper))/255
        green <- c(seq(low[3, 1], mid[3, 1], length = lower), 
            seq(mid[3, 1], high[3, 1], length = upper))/255
        blue <- c(seq(low[2, 1], mid[2, 1], length = lower), 
            seq(mid[2, 1], high[2, 1], length = upper))/255
        if (isodd) {
            red <- red[-(lower + 1)]
            green <- green[-(lower + 1)]
            blue <- blue[-(lower + 1)]
        }
    }
    rgb(red, blue, green)
}


greenred = function (n)  {
    colorpanel(n, "green", "black", "red")
}

odd = function (x) {
    x%%2 == 1
}

even = function (x) {
    x%%2 == 0
}


sample_matrix_to_color_assignments = function(sampleAnnotationsMatrix, colors) {

	if (missing(colors))
		colors = rainbow(nrow(sampleAnnotationsMatrix))

	nsamples = nrow(sampleAnnotationsMatrix);

	if (length(colors) < nrow(sampleAnnotationsMatrix))
		stop("Error, only ", length(colors), " colors specified, but have ", nsamples, " samples");

	for (i in 1:nrow(sampleAnnotationsMatrix)) {
		c = colors[i]
		sampleAnnotationsMatrix[i,] = sapply(sampleAnnotationsMatrix[i,], function(x) ifelse( x, as.character(c), 'white'))
	}

	return(sampleAnnotationsMatrix);

}

##edgeR
library(edgeR)
sample_matrix_to_color_assignments = function(sampleAnnotationsMatrix, colors) {

	if (missing(colors))
		colors = rainbow(nrow(sampleAnnotationsMatrix))

	nsamples = nrow(sampleAnnotationsMatrix);

	if (length(colors) < nrow(sampleAnnotationsMatrix))
		stop("Error, only ", length(colors), " colors specified, but have ", nsamples, " samples");

	for (i in 1:nrow(sampleAnnotationsMatrix)) {
		c = colors[i]
		sampleAnnotationsMatrix[i,] = sapply(sampleAnnotationsMatrix[i,], function(x) ifelse( x, as.character(c), 'white'))
	}

	return(sampleAnnotationsMatrix);

}
sample_correlation_heatmap <- function(data.TMM.matrix){
      ##Plot samples  correlation
    #set group/condition with replicates in a matrix
    samples_data<-cbind.data.frame(group,colnames(data.TMM.matrix))
    sample_types = as.character(unique(samples_data[,1]))
    rep_names = as.character(samples_data[,2])
    data_r=matrix(0,nrow(data.TMM.matrix),length(rep_names),byrow = T)
    for(n in 1:length(rep_names))
    {
      data_r[,n]<-data.TMM.matrix[,rep_names[n]]
    }
    
    colnames(data_r)<-rep_names
    data_r = data_r[, colnames(data_r) %in% samples_data[,2], drop=F ] #select samples from matrix based on group comparisons
    nsamples = length(sample_types)
    sample_colors = rainbow(nsamples)
    names(sample_colors) = sample_types
    #divide replicates to group/condition
    sample_type_list = list()
    for (i in 1:nsamples) {
      samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
      sample_type_list[[sample_types[i]]] = as.vector(samples_want)
    }
    sample_factoring = colnames(data_r)
    
    for (i in 1:nsamples) {
      sample_type = sample_types[i]
      replicates_want = sample_type_list[[sample_type]]
      sample_factoring[ colnames(data_r) %in% replicates_want ] = sample_type
    }
    # reorder according to sample type.
    tmp_sample_reordering = order(sample_factoring)
    data_r = data_r[,tmp_sample_reordering,drop=F]
    sample_factoring = sample_factoring[tmp_sample_reordering]
    initial_matrix = data_r # store before doing various data transformations
#    data_r = log2(data_r+1)
    sample_factoring = colnames(data_r)
    sampleAnnotations = matrix(ncol=ncol(data_r),nrow=nsamples)
    for (i in 1:nsamples) {
      sampleAnnotations[i,] = colnames(data_r) %in% sample_type_list[[sample_types[i]]]
    }
    sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
    #set the same color for replicates of each group, separately
    sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
    rownames(sampleAnnotations) = as.vector(sample_types)
    colnames(sampleAnnotations) = colnames(data.TMM.matrix)
    data_r = as.matrix(data_r) # convert to matrix
    #generate correlation between samples
    sample_cor = cor(data_r, method='pearson', use='pairwise.complete.obs')
    sample_dist = dist(t(data_r), method='euclidean')
    hc_samples = hclust(sample_dist, method='complete')
    sample_cor_for_plot = sample_cor
    pdf("Sample_correlation_matrix.pdf")
    heatmap.3(sample_cor_for_plot, dendrogram='both', Rowv=as.dendrogram(hc_samples), 
              Colv=as.dendrogram(hc_samples), col = greenred(75), scale='none', symm=TRUE, key=TRUE,
              density.info='none', trace='none', symkey=FALSE, symbreaks=F, margins=c(10,10), cexCol=1, 
              cexRow=1, cex.main=0.75, main=paste("sample correlation matrix of DE genes"),
              ColSideColors=sampleAnnotations, RowSideColors=t(sampleAnnotations))
    dev.off()
}
run_edge_R <- function(data.TMM, A, B, annotation, group, nsample){
  if(nsample == ''){
    nsample = 3
  }
  no_fdr <- 0
  ndata<-data.TMM[,c(A,B)] #Selecting columns
  y <- DGEList(counts=ndata,group=group) #Creating edgeR object 
  keep<-rowSums(cpm(y) > 1) >= nsample ##At least four of the samples should have cpm above 1. CPM is one of the metrics commonly used to measure the expression
  y<-y[keep,] 
  #y <- calcNormFactors(y) to normalize the data but it may not be needed for your data
  y <- estimateCommonDisp(y) #Estimate Common Negative Binomial Dispersion by Conditional Maximum Likelihood. Please read edgeR manual
  y <- estimateTagwiseDisp(y) # Estimate Empirical Bayes Tagwise Dispersion Values. Please read edgeR manual

  # Plot BCV/Dispersion
  png("plotBCV.png", width=1500,height=1000,res=125)
  plot(BCV)	
  dev.off()
	
  et = exactTest(y) #Pairwise comparison ()
  tTags = topTags(et,n=NULL) #list tested genes
  FC.CPM.FDR<-cbind(row.names(tTags$table),tTags$table) #name the rows. 
  colnames(FC.CPM.FDR)[1]<-"gene_id" #Rename the column to use for adding annotation info frmo annotation matrix/table/data.frame
  write.table(y$samples,"library-size-norm.factors.txt",row.names=F,quote=F,sep="\t") #Write normalization factors. 
  
  png("plot_MA.png",width=1500,height=1000,res=125)
  plot_MA(et$table$logCPM,et$table$logFC,et$table$PValue)
  dev.off()
  png("Volcano.png",width=1000,height=1000,res=125)
  plot_Volcano(et$table$logFC, et$table$PValue) 
  dev.off()
  write.table(FC.CPM.FDR,"logFC-CPM-FDR.txt",row.names=F,quote=F,sep="\t")
  
  png("plot_MA-FDR.png",width=1500,height=1000,res=125)
  plot_MA(FC.CPM.FDR$logCPM,FC.CPM.FDR$logFC,FC.CPM.FDR$FDR)
  dev.off()
  png("Volcano-FDR.png",width=1000,height=1000,res=125)
  plot_Volcano(FC.CPM.FDR$logFC, FC.CPM.FDR$FDR) 
  dev.off()

  
  
  ##To get DE genes
  library(reshape2)
  de <- decideTestsDGE(et, p=0.05,adjust.method="fdr",lfc=1) #Select genes with FDR < 0.05 and log fold change >=1 or <=-1
  DE_summary<-summary(de)
  if(sum(DE_summary[c(1,3)]) == 0){
    de <- decideTestsDGE(et, p=0.05,adjust.method="none",lfc=1) #Select genes with FDR < 0.05 and log fold change >=1 or <=-1
    DE_summary<-summary(de)
    print(' PLEASE NOTE THAT WE SWITCHED OFF THE MULTIPLE TESTING DUE TO THE FACT THAT AFTER CORRECTION NO GENE WAS SELECTED BY FDR <0.01')
    no_fdr <- 1
  }
  write.table(DE_summary,"DE-summary-gene-count.txt",sep="\t",quote=F,row.names=T,col.names=F)
  isDE<-as.logical(de)
  DEnames<-as.matrix((rownames(y)[isDE]))
  colnames(DEnames)<-"gene_id"
  tmp <- merge(DEnames,annotation,by="gene_id")#,all.x=T)
  FinalDE<-merge(tmp,FC.CPM.FDR,by="gene_id")#,all.x=T)
  
  if(sum(DE_summary[c(1,3)])>=2){
    ##heatmap of DE and sample correlation
    library(cluster)
    library(gplots)
    
    FinalDE<-merge(FinalDE,data.TMM[,c(1,c(A,B))],by="gene_id",all.x=T)
    #Write significant or DE genes with stats and original normalized expression
    write.table(FinalDE,"DE-edgeR-results.txt",row.names=F,sep="\t",quote=F,col.names=T)
    
    ##Sample correlation heatmap
    data.TMM.matrix = FinalDE[,c(7:ncol(FinalDE))] # remove the gene column since its now the rowname value
    rownames(data.TMM.matrix) = FinalDE$gene_id # set rownames to gene identifiers
    data.TMM.matrix = as.matrix(data.TMM.matrix) # convert to matrix
    sample_correlation_heatmap(data.TMM.matrix = data.TMM.matrix)
    rm(data.TMM.matrix)
    
    ##Heatmap of DEG
    data.TMM.matrix = FinalDE[,c(7:ncol(FinalDE))] # remove the gene column since its now the rowname value
    rownames(data.TMM.matrix) = FinalDE$gene_id # set rownames to gene identifiers
    data.TMM.matrix = as.matrix(data.TMM.matrix) # convert to matrix
    
    cr = cor(data.TMM.matrix, method='spearman')
    data.TMM.matrix = log2(data.TMM.matrix+1)
    data.TMM.matrix = t(scale(t(data.TMM.matrix), scale=F)) # center rows, mean substracted
    head(data.TMM.matrix)
    gene_dist = dist(data.TMM.matrix, method='euclidean')
    hc_genes = hclust(gene_dist, method='complete')
    hc_samples = hclust(as.dist(1-cr), method="complete") # cluster conditions
    if(length(unlist(as.hclust(hc_genes)[3])) < 6){
      gene_partition_assignments <- cutree(as.hclust(hc_genes), k=length(unlist(as.hclust(hc_genes)[3])))
    }else{
      gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6)
    }
    partition_colors = rainbow(length(unique(gene_partition_assignments)), start=0.4, end=0.95)
    gene_colors = partition_colors[gene_partition_assignments]
    quantBrks = quantile(data.TMM.matrix, c(0.03, 0.97))
    myheatcol = bluered(75)
    pdf("DE-edgeR-results-heatmap_with_genes.pdf", height=10, width= 9)
    cex.main = 0.7
    heatmap.2(data.TMM.matrix, dendrogram='row', Rowv=as.dendrogram(hc_genes), Colv=F, 
              col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,2.5,1.5), margins=c(7,7), breaks=seq(quantBrks[1], quantBrks[2], length=76),main="log2 (normalized expression + 1)")
    dev.off()
    ##Boxplot per gene in FinalDE
    return(FinalDE)
  }
  if(sum(DE_summary[c(1,3)]) == 1){
    FinalDE<-merge(FinalDE[(FinalDE$PValue<=0.05),],data.TMM[,c(1,c(A,B))],by="gene_id",all.x=T)
    write.table(FinalDE,"DE-edgeR-results.txt",row.names=F,sep="\t",quote=F,col.names=T)
    data.TMM.matrix = FinalDE[,c(7:ncol(FinalDE))] # remove the gene column since its now the rowname value
    rownames(data.TMM.matrix) = FinalDE$gene_id # set rownames to gene identifiers
    data.TMM.matrix = as.matrix(data.TMM.matrix) # convert to matrix
    data.TMM.matrix <- cbind.data.frame(Expression = t(data.TMM.matrix), Group = group)
    return(FinalDE)
  }
  if(sum(DE_summary[c(1,3)]) == 0){
    print('There are either no genes significant with FDR/P-value < 0.05 or less than two genes')
  }
}
