#@ This file contains all the functions to 
#@ get the differential table

TableDiff_print <- function(input,BaseContrast,resDiff, info = NULL) 
{
  VarInt = input$VarInt
  dds = resDiff$dds
  counts = resDiff$counts
  target = resDiff$target
  group = ""
  if(length(VarInt)>1) for(i in VarInt){ group = paste(group,target[,i], sep="-") }
  else group = target[,VarInt]
  conds = unique(group)
  
  result = list()
  alpha = as.numeric(input$AlphaVal)
  cooksCutoff = ifelse(input$CooksCutOff!='Auto',ifelse(input$CooksCutOff!=Inf,input$CutOffVal,Inf),TRUE)
  result[[input$ContrastList_table]] <- results(dds,contrast=BaseContrast[,input$ContrastList_table],pAdjustMethod=input$AdjMeth,
                                                cooksCutoff=cooksCutoff,
                                                independentFiltering=input$IndFiltering,alpha=alpha,parallel = TRUE)
  
  #names(result) <- gsub("_", " ", names(result))
  if (is.null(info))  info <- data.frame(Id = rownames(result[[1]]))
  else names(info)[1] <- "Id"
  if (any(duplicated(info[, 1]))) stop("Duplicated IDs in the annotations")
  counts <- data.frame(Id = rownames(counts(dds)), counts(dds), round(counts(dds, normalized = TRUE)))
  colnames(counts) <- c("Id", colnames(counts(dds)), paste0("norm.", colnames(counts(dds))))
  bm <- data.frame(Id = rownames(result[[1]]), baseMean = round(result[[1]][,"baseMean"], 2))
  
  base <- merge(merge(info, counts, by = "Id", all.y = TRUE), bm, by = "Id")
  tmp <- base[, paste("norm", colnames(counts(dds)), sep = ".")]
  for (cond in conds) {
    base[, cond] <- round(apply(as.data.frame(tmp[, group == cond]), 1, mean), 0)
  }
  complete.complete <- base
  complete <- vector("list", length(result))
  names(complete) <- names(result)
  for (name in names(result)) {
    complete.name <- base
    
    ### ??????????
    conds.supp <- setdiff(conds, gsub("\\(|\\)", "", unlist(strsplit(name, " vs "))))
    if (length(conds.supp) > 0) {
      complete.name <- complete.name[, -which(names(complete.name) %in% conds.supp)]
      samples.supp <- colnames(counts(dds))[group %in% conds.supp]
      col.supp <- c(samples.supp, paste("norm", samples.supp, sep = "."))
      complete.name <- complete.name[, -which(names(complete.name) %in% col.supp)]
    }
    ### ??????????
    
    res.name <- data.frame(Id = rownames(result[[name]]), 
                           FC = round(2^(result[[name]][, "log2FoldChange"]), 3),
                           log2FoldChange = round(result[[name]][, "log2FoldChange"], 3),
                           pvalue = result[[name]][, "pvalue"], 
                           padj = result[[name]][, "padj"])
    complete.name <- merge(complete.name, res.name, by = "Id")
    mcols.add <- data.frame(Id = rownames(counts(dds)), 
                            dispGeneEst = mcols(dds)$dispGeneEst, dispFit = mcols(dds)$dispFit, 
                            dispMAP = mcols(dds)$dispMAP, dispersion = mcols(dds)$dispersion, 
                            betaConv = mcols(dds)$betaConv, maxCooks = mcols(dds)$maxCooks)
    if (is.null(cooksCutoff)){
      m <- nrow(attr(dds, "modelMatrix"))
      p <- ncol(attr(dds, "modelMatrix"))
      cooksCutoff <- qf(0.99, p, m - p)
    }
    mcols.add$outlier <- ifelse(mcols(dds)$maxCooks > cooksCutoff, "Yes", "No")
    complete.name <- merge(complete.name, mcols.add, by = "Id")
    complete[[name]] <- complete.name
    complete.name=complete.name[order(complete.name$padj),]
    significant.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv), ]
    up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange>=0.0), ]
    ## useless order
    #up.name <- up.name[order(up.name$padj), ]
    down.name <- complete.name[which(complete.name$padj<=alpha & complete.name$betaConv & complete.name$log2FoldChange<=0.0), ]
    #down.name <- down.name[order(down.name$padj), ]
    
    name <- gsub(" ", "", name)
    keep <- c("Id","baseMean","FC","log2FoldChange","padj")
    #       complete.complete[, paste(name, keep, sep = ".")] <- complete.name[, keep]
  }
  #return(list(complete=complete.name,up=up.name,down=down.name))
  return(list(significant=significant.name[,keep],complete=complete.name[,keep],up=up.name[,keep],down=down.name[,keep]))
}
