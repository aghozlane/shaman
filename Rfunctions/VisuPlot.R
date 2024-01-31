#@ This file contains all the functions for the
#@ visualisation plots of SHAMAN


##                       ##
##        Barplot ####
##                       ##

Plot_Visu_Barplot <-
  function(input, resDiff, colors, others, plot = "classic", community = NULL)
  {
    req(colors)
    target = resDiff$target
    #Sort the palette for the barchart
    rgb <- grDevices::col2rgb(colors)
    tsp <- TSP::as.ATSP(dist(t(rgb)))
    sol <- TSP::solve_TSP(tsp, control = list(repetitions = 1e3))
    colors <- colors[sol]
    ## Get Input for BarPlot
    VarInt = input$VisuVarInt
    ind_taxo = input$selectTaxoPlot
    nbKept = length(ind_taxo)
    if(is.null(community)){
      others <- others[!(others %in% ind_taxo)]
      other_tmp_combined = GetDataToPlot(input, resDiff, VarInt, others)
      others_counts_tmp_combined = other_tmp_combined$counts
      OthersSamplesNames = other_tmp_combined$namesCounts
      if (!is.null(OthersSamplesNames))
        OthersnamesTax = colnames(others_counts_tmp_combined)
      else if (!is.null(others))
        OthersnamesTax = others
    }
    else{
      
    }
    dataNull = data.frame(x = c(0, 0), y = c(1, 2))
    
    all_tmp_combined <-
      GetDataToPlot(input, resDiff, VarInt, c(ind_taxo, others))
    all_counts_tmp_combined <- all_tmp_combined$counts
    allSamplesNames <- all_tmp_combined$namesCounts
    allnamesTax = c(ind_taxo, others)
    plotd3 = nvd3Plot(
      x ~ y ,
      data = dataNull,
      type = "multiBarChart",
      id = 'barplotTaxoNyll',
      height = ifelse(
        is.na(input$heightVisu),
        800,
        as.character(input$heightVisu)
      ),
      width = if (input$modifwidthVisu) {
        as.character(input$widthVisu)
      }
    )
    gg = NULL
    
    if (length(VarInt) > 0 &&
        !is.null(all_counts_tmp_combined) &&
        nrow(all_counts_tmp_combined) > 0)
    {
      val = NULL
      #Get the modalities
      for (i in 1:length(VarInt))
      {
        ## Replace "-" by "."
        target[, VarInt[i]] =  gsub("-", ".", target[, VarInt[i]])
        
        Tinput = paste("input$", "ModVisu", VarInt[i], sep = "")
        expr = parse(text = Tinput)
        ## All the modalities for all the var of interest
        val = c(val, eval(expr))
      }
      
      ## Create the others data frame for the plot function
      all_dataBarPlot_mat = c()
      all_tmp_mat = matrix(0, ncol = 3, nrow = length(allnamesTax))
      all_tmp_counts = c()
      
      ## Taxo
      all_tmp_mat[1:length(allnamesTax), 1] = allnamesTax
      for (i in 1:(nrow(all_counts_tmp_combined)))
      {
        ## Counts
        
        tmpProp =  all_counts_tmp_combined[i, ]
        if (input$CountsOrProp == "prop")
        {
          tmpProp = round(tmpProp / sum(tmpProp), 3)
          tmpProp = as.numeric(tmpProp / sum(tmpProp) * 100)
        }
        all_tmp_counts = c(all_tmp_counts, tmpProp)
        
        ## Meta data
        all_tmp_mat[1:length(allnamesTax), 3] = as.character(rep(allSamplesNames[i], length(allnamesTax)))
        
        ## Combined the sample
        all_dataBarPlot_mat = rbind(all_dataBarPlot_mat, all_tmp_mat)
      }
      
      all_dataBarPlot_mat = as.data.frame(all_dataBarPlot_mat)
      colnames(all_dataBarPlot_mat) = c("Taxonomy", "Proportions", "AllVar")
      all_dataBarPlot_mat[, 2] = all_tmp_counts
      
      main_dataBarPlot_mat <- all_dataBarPlot_mat %>%
        dplyr::group_by(AllVar) %>%
        dplyr::slice(1:nbKept)
      
      others_dataBarPlot_mat <-  all_dataBarPlot_mat %>%
        dplyr::group_by(AllVar) %>%
        dplyr::slice((nbKept + 1):ncol(all_counts_tmp_combined)) %>%
        dplyr::summarise(Proportions = sum(Proportions, na.rm = TRUE)) %>%
        dplyr::mutate(Taxonomy = 'Others')
      
      others_dataBarPlot_mat <-
        as.data.frame(others_dataBarPlot_mat)
      main_dataBarPlot_mat <- as.data.frame(main_dataBarPlot_mat)
      if(length(VarInt) == 1){
        main_dataBarPlot_mat$AllVar <- factor(main_dataBarPlot_mat$AllVar, levels = val)
        
        main_dataBarPlot_mat <- rbind(main_dataBarPlot_mat, others_dataBarPlot_mat)
        
        main_dataBarPlot_mat <- main_dataBarPlot_mat %>%
          dplyr::arrange(AllVar)
      }
      else{
        main_dataBarPlot_mat$AllVar <- factor(main_dataBarPlot_mat$AllVar, levels = unique(others_dataBarPlot_mat$AllVar))
        main_dataBarPlot_mat <- rbind(main_dataBarPlot_mat, others_dataBarPlot_mat)
      }
      main_dataBarPlot_mat <- as.data.frame(main_dataBarPlot_mat)
      if (input$SensPlotVisu == "Vertical")
        Sens = "multiBarChart"
      else
        Sens = "multiBarHorizontalChart"
      XRotate = input$rotateXLabel
      colors = rev(rep(colors, ceiling(nbKept / length(colors))))
      plotd3 <-
        nvd3Plot(
          Proportions ~ AllVar |
            Taxonomy,
          data = main_dataBarPlot_mat,
          type = Sens,
          id = 'barplotTaxo',
          height = ifelse(
            is.na(input$heightVisu),
            toString(800),
            toString(as.character(input$heightVisu))
          ),
          width = if (input$modifwidthVisu) {
            as.character(input$widthVisu)
          }
        )
      plotd3$chart(stacked = TRUE, color = colors)
      if (input$SensPlotVisu == "Vertical") {
        plotd3$chart(reduceXTicks = FALSE)
        plotd3$xAxis(rotateLabels = XRotate)
      }
      
      
      
      ##                                 ##
      ## Same plot in ggplot2 for export
      ##                                 ##
      
      # tax.colors=rep(c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78", "#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5"
      #                  ,"#8c564b","#c49c94","#e377c2","#f7b6d2","#7f7f7f", "#c7c7c7","#bcbd22","#dbdb8d","#17becf","#9edae5"),ceiling(nbKept/20))
      
      # main_dataBarPlot_mat$Taxonomy = factor(main_dataBarPlot_mat$Taxonomy,
      #                                        levels = c(allnamesTax[!(allnamesTax %in% others)], 'Others'))
      # main_dataBarPlot_mat$AllVar = factor(main_dataBarPlot_mat$AllVar, levels = rev(unique(main_dataBarPlot_mat$AllVar)))
      
      # main_dataBarPlot_mat <- main_dataBarPlot_mat %>%
      #   dplyr::arrange(AllVar)
      
      #colors <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(main_dataBarPlot_mat$Taxonomy)))
      main_dataBarPlot_mat$Taxonomy <- factor(main_dataBarPlot_mat$Taxonomy, levels = unique(main_dataBarPlot_mat$Taxonomy))
      #colors <- rev(rep(colors, length.out = length(unique(main_dataBarPlot_mat$Taxonomy))))
      colors <- rev(colors[1:length(unique(main_dataBarPlot_mat$Taxonomy))])
      gg = ggplot(main_dataBarPlot_mat,
                  aes(
                    x = AllVar,
                    y = Proportions,
                    fill = forcats::fct_rev(Taxonomy)
                  )) 
      if (input$CountsOrProp == "counts" &&
          input$positionBarPlot == "fill")
        gg = gg + geom_bar(position=position_stack(reverse = TRUE), stat='identity')
      else
        gg = gg + geom_bar(position = input$positionBarPlot, stat = 'identity')
      gg = gg + theme_bw() + scale_fill_manual(values = colors)
      gg = gg + theme(
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = XRotate),  
        axis.text.y = element_text(size = input$sizeBarplotTitle),
        axis.title.x = element_text(size = input$sizeBarplotTitle * 1.3),
        axis.title.y = element_text(size = input$sizeBarplotTitle * 1.3),
        legend.title = element_text(size = input$sizeBarplotLegend *1.3),
        legend.text = element_text(size = input$sizeBarplotLegend)) +
        guides(fill = guide_legend(input$TaxoSelect))
      
      if (input$CountsOrProp == "prop")
        gg = gg + labs(y = "Relative abundance (%)", x = "")
      if (input$CountsOrProp == "counts" &&
          input$NormOrRaw == "norm")
        gg = gg + labs(y = "Normalized counts", x = "")
      else if (input$CountsOrProp == "counts" &&
               input$NormOrRaw != "norm")
        gg = gg + labs(y = "Raw counts", x = "")
      if (input$SensPlotVisu == "Horizontal")
        gg = gg + coord_flip()
      
    }
    return(list(plotd3 = plotd3, gg = gg))
  }




##                       ##
##          HEATMAP ####
##                       ##
Plot_Visu_Heatmap <- function(input, resDiff, export = FALSE) {
  VarInt = input$VisuVarInt
  ind_taxo = input$selectTaxoPlot
  sortrow = input$SortHeatRow
  sortcol = input$SortHeatColumn
  
  counts_tmp_combined = GetDataToPlot(input, resDiff, VarInt, ind_taxo)$counts
  
  if (!is.null(counts_tmp_combined) &&
      nrow(counts_tmp_combined) > 0)
  {
    ## Transform to log2
    counts_tmp_combined = log2(GetDataToPlot(input, resDiff, VarInt, ind_taxo)$counts +
                                 1)
    
    col <- switch(
      input$colors,
      "green-blue" = colorRampPalette(brewer.pal(9, "GnBu"))(256),
      "blue-white-red" = colorRampPalette(rev(brewer.pal(9, "RdBu")))(200),
      "purple-white-orange" = colorRampPalette(rev(brewer.pal(9, "PuOr")))(200),
      "red-yellow-green" = colorRampPalette(rev(brewer.pal(9, "RdYlGn")))(200)
    )
    
    ## Transpose matrix if Horizontal
    if (input$SensPlotVisu == "Horizontal")
      counts_tmp_combined = t(as.matrix(counts_tmp_combined))
    if (nrow(counts_tmp_combined) == 1)
      sortrow = "No"
    if (ncol(counts_tmp_combined) == 1)
      sortcol = "No"
    if (!export) {
      plot = d3heatmap(
        counts_tmp_combined,
        dendrogram = "none",
        Rowv = (sortrow == "Yes"),
        Colv = (sortcol == "Yes"),
        na.rm = TRUE,
        width = ifelse(input$modifwidthVisu, input$widthVisu, "100%"),
        height = input$heightVisu,
        show_grid = FALSE,
        colors = col,
        scale = input$scaleHeatmap,
        cexRow = as.numeric(input$LabelSizeHeatmap),
        margins = c(12, 30),
        cexCol = as.numeric(input$LabelSizeHeatmap),
        offsetCol = input$LabelColOffsetHeatmap,
        offsetRow = input$LabelRowOffsetHeatmap
      )
      
    }
    if (export) {
      dendrogram = "none"
      if (input$SortHeatColumn == "Yes" &&
          input$SortHeatRow == "Yes")
        dendrogram = "both"
      else if (input$SortHeatColumn == "Yes")
        dendrogram = "column"
      else if (input$SortHeatRow == "Yes")
        dendrogram = "row"
      plot = heatmap.2(
        counts_tmp_combined,
        dendrogram = dendrogram,
        Rowv = (sortrow == "Yes"),
        Colv = (sortcol == "Yes"),
        na.rm = TRUE,
        density.info = "none",
        margins = c(
          as.numeric(input$lowerMargin),
          as.numeric(input$rightMargin)
        ),
        trace = "none",
        srtCol = 45,
        col = col,
        scale = input$scaleHeatmap,
        cexRow = input$LabelSizeHeatmap,
        cexCol = input$LabelSizeHeatmap,
        offsetCol = input$LabelColOffsetHeatmap,
        offsetRow = input$LabelRowOffsetHeatmap,
        symm = FALSE,
        symkey = FALSE,
        symbreaks = FALSE
      )
    }
    return(plot)
  }
}



##                       ##
##          BOXPLOTS ####
##                       ##
Plot_Visu_Boxplot <-
  function(input, resDiff, alpha = 0.7, dataDiff, colors) {
    req(colors)
    gg = NULL
    
    ## Get Input for BoxPlot
    VarInt = input$VisuVarInt
    ind_taxo = input$selectTaxoPlot
    
    tmp_merge = GetDataToPlot(input, resDiff, VarInt, ind_taxo, aggregate =
                                FALSE)
    counts_tmp_combined = tmp_merge$counts
    levelsMod = tmp_merge$levelsMod
    nbKept = length(ind_taxo)
    
    if (!is.null(counts_tmp_combined) &&
        nrow(counts_tmp_combined) > 0 && !is.null(levelsMod))
    {
      if (input$typeDataBox == "Relative")
        counts_tmp_combined = tmp_merge$prop_all
      else
        counts_tmp_combined = log2(counts_tmp_combined + 1)
      if (nbKept == 1)
        colnames(counts_tmp_combined) = ind_taxo
      
      ## Create the data frame for the plot function
      dataBarPlot_mat = c()
      tmp_mat = matrix(0, ncol = 3, nrow = nbKept)
      tmp_counts = c()
      
      for (i in 1:(nrow(counts_tmp_combined)))
      {
        ## Taxo
        tmp_mat[1:nbKept, 1] = colnames(counts_tmp_combined)
        
        ## Counts
        tmpProp = counts_tmp_combined[i, ]
        tmp_counts = c(tmp_counts, tmpProp)
        
        ## Meta data
        tmp_mat[1:nbKept, 3] = as.character(rownames(counts_tmp_combined)[i], nbKept)
        
        ## Conbined the sample
        dataBarPlot_mat <- rbind(dataBarPlot_mat, tmp_mat)
        #print(dataBarPlot_mat)
      }
      
      
      dataBarPlot_mat = as.data.frame(dataBarPlot_mat)
      
      colnames(dataBarPlot_mat) = c("Taxonomy", "Value", "Samples")
      dataBarPlot_mat[, 2] = tmp_counts
      
      if (!is.null(dataDiff)) {
        dataDiff$complete$pvalue_adjusted[is.na(dataDiff$complete$pvalue_adjusted)] = 0
        labels <-
          sprintf(
            "%s\np-value: %.3f",
            dataDiff$complete$Id,
            dataDiff$complete$pvalue_adjusted
          )
        dataBarPlot_mat$TaxonomyPvalue = factor(
          dataBarPlot_mat$Taxonomy,
          levels = as.vector(dataDiff$complete$Id),
          labels = labels
        )
      }
      if (is.null(input$BoxColorBy) ||
          length(VarInt) <= 1) {
        dataBarPlot_mat$Colors = dataBarPlot_mat$Samples
      }
      if (!is.null(input$BoxColorBy) && length(VarInt) > 1)
      {
        tmp = strsplit(as.character(dataBarPlot_mat$Samples), "-")
        ind = which(VarInt %in% input$BoxColorBy)
        dataBarPlot_mat$Colors = rapply(tmp, function(x)
          paste(x[ind], collapse = "-"), how = "unlist")
      }
      dataBarPlot_mat$Samples = factor(dataBarPlot_mat$Samples, levels = levelsMod)
      dataBarPlot_mat$Colors = factor(dataBarPlot_mat$Colors, levels = levelsMod)
      #column that is used to print the p-value under the title of each graph
      gg = ggplot(dataBarPlot_mat, aes(x = Samples, y = Value, fill = Colors))  + geom_boxplot(alpha =
                                                                                                 alpha) + theme_bw()
      gg = gg + theme(
        axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        strip.text = element_text(size = 20, hjust = 0.5)
      )
      gg = gg + ylab(paste(input$typeDataBox, "abundance")) + scale_fill_manual(values = colors) + guides(fill =
                                                                                                            FALSE)
      if (input$CheckAddPointsBox)
        gg = gg + geom_point(position = position_jitterdodge(dodge.width = 0.9))
      if (input$SensPlotVisu == "Horizontal")
        gg = gg + coord_flip()
      if (nbKept > 1)
        gg = gg + facet_wrap( ~ Taxonomy, scales = input$ScaleBoxplot)
      if (!is.null(dataDiff) &&
          !is.null(input$SelectSpecifTaxo) &&
          (length(input$ContrastList_table_Visu) == 1) &&
          input$SelectSpecifTaxo == "Diff")
        gg = gg + facet_wrap( ~ TaxonomyPvalue, scales = input$ScaleBoxplot)
      req(gg)
    }
    
    return(gg)
  }

##                       ##
##          KRONA ####
##                       ##
Plot_Visu_Krona <- function(input, resDiff, CT_OTU, taxo_table) {
  res = NULL
  ## Get Input for Krona
  VarInt = input$VisuVarInt
  ind_taxo = input$selectTaxoPlot
  
  ## Removed column with only 1 modality
  ind = which(apply(
    taxo_table,
    2,
    FUN = function(x)
      length(unique(x[!is.na(x)]))
  ) == 1)
  if (length(ind) > 0)
    taxo_table = taxo_table[, -ind]
  
  #print(counts_tmp_combined)
  if (nrow(CT_OTU) > 0 &&
      !is.null(CT_OTU) &&
      nrow(taxo_table) > 0 && !is.null(taxo_table))
  {
    tmp = CreateTableTree(input, resDiff, CT_OTU, taxo_table, VarInt)
    
    if (nrow(tmp$counts) > 0 && !is.null(tmp$counts))
    {
      merge_dat = c()
      for (cond in tmp$levelsMod) {
        merge_dat = rbind(merge_dat, cbind(cond, merge(
          round(tmp$counts[cond, ]), taxo_table, by = "row.names"
        )))
      }
      # Reorder columns
      merge_dat = merge_dat[, c(3, 1, 4:dim(merge_dat)[2], 2)]
      # Remove zero counts
      # Required for Krona
      res = merge_dat[merge_dat[, 1] > 0, ]
    }
  }
  return(res)
}

##                       ##
##      Phylo PLOT ####
##                       ##
Plot_Visu_Phylotree = function(input,
                               resDiff,
                               CT_OTU,
                               taxo_table,
                               treeseq) {
  res = NULL
  VarInt = input$VisuVarInt
  ind_taxo = input$selectTaxoPlot
  
  ## Removed column with only 1 modality
  ind = which(apply(
    taxo_table,
    2,
    FUN = function(x)
      length(unique(x[!is.na(x)]))
  ) == 1)
  if (length(ind) > 0)
    taxo_table = taxo_table[, -ind]
  
  #print(counts_tmp_combined)
  if (nrow(CT_OTU) > 0 &&
      !is.null(CT_OTU) &&
      nrow(taxo_table) > 0 && !is.null(taxo_table))
  {
    tmp = CreateTableTree(input, resDiff, CT_OTU, taxo_table, VarInt)
    
    if (nrow(tmp$counts) > 0 && !is.null(tmp$counts))
    {
      #merge_dat = c()
      #for(cond in tmp$levelsMod){
      #  merge_dat = rbind(merge_dat, cbind(cond, merge(round(tmp$counts[cond,]), taxo_table, by="row.names")))
      #}
      # Reorder columns
      #merge_dat=merge_dat[, c(3,1,4:dim(merge_dat)[2],2)]
      # Remove zero counts
      # Required for Krona
      #res = merge_dat[merge_dat[,1]>0,]
      counts = round(t(tmp$counts))
      data = as.data.frame(rbind(c("#name", colnames(counts)), cbind(rownames(counts), counts)))
      res = PhyloTreeMetaR(treeseq, data)
    }
  }
  #if(input$TransDataPhyloTree =="log2") data = cbind(target,log2(t(counts)+1),div)
  #else if(input$TransDataPhyloTree =="none") data = cbind(target,t(counts),div)
  #print(counts)
  return(res)
}


##                       ##
##      SCATTER PLOT ####
##                       ##
Plot_Visu_Scatterplot <-
  function(input,
           resDiff,
           export = FALSE,
           lmEst = FALSE,
           CorEst = FALSE) {
    plot = NULL
    regCoef = NULL
    Rsq = NULL
    cor.est = NULL
    cor.pvalue = NULL
    div = NULL
    dds = resDiff$dds
    
    if (input$NormOrRaw == "norm")
    {
      counts = as.data.frame(round(counts(dds, normalized = TRUE)))
    }
    else
    {
      counts = as.data.frame(round(counts(dds, normalized = FALSE)))
    }
    
    target = as.data.frame(resDiff$target)
    
    ## Get the diversity values
    tmp_div = Plot_Visu_Diversity(input, resDiff, ForScatter = TRUE)$dataDiv
    
    if (!is.null(tmp_div)) {
      div = cbind(
        round(tmp_div$value[tmp_div$diversity == "Alpha"], 3),
        round(tmp_div$value[tmp_div$diversity == "Shannon"], 3),
        round(tmp_div$value[tmp_div$diversity == "Inv.Simpson"], 3),
        round(tmp_div$value[tmp_div$diversity == "Simpson"], 3)
      )
      colnames(div) = c("Alpha div", "Shannon div", "Inv.Simpson div", "Simpson div")
    }
    if (input$TransDataScatter == "log2")
      data = cbind(target, log2(t(counts) + 1), div)
    else if (input$TransDataScatter == "none")
      data = cbind(target, t(counts), div)
    
    ## Get Input for ScatterPlot
    Xvar = input$Xscatter
    Yvar = input$Yscatter
    ColBy = input$ColorBy
    PchBy = input$PchBy
    PointSize = input$PointSize
    NamesData = colnames(data)
    x_var = NULL
    y_var = NULL
    col_var = NULL
    symbol_var = NULL
    size_var = NULL
    
    if (!is.null(Xvar)) {
      if (Xvar %in% NamesData) {
        x_var = data[, Xvar]
      }
    }
    if (!is.null(Yvar)) {
      if (Yvar %in% NamesData) {
        y_var = data[, Yvar]
      }
    }
    if (!is.null(ColBy)) {
      if (ColBy %in% NamesData) {
        col_var = data[, ColBy]
      }
    }
    if (!is.null(PchBy)) {
      if (PchBy %in% NamesData) {
        symbol_var = data[, PchBy]
      }
    }
    if (!is.null(PointSize)) {
      if (PointSize %in% NamesData) {
        size_var = data[, PointSize]
      }
    }
    
    if (!export &&
        !input$AddRegScatter &&
        !lmEst && !CorEst && !is.null(x_var) && !is.null(y_var)) {
      plot = scatterD3(
        x = x_var,
        y = y_var,
        lab = rownames(data),
        xlab = Xvar,
        ylab = Yvar,
        col_var = col_var,
        col_lab = ColBy,
        symbol_var = symbol_var,
        symbol_lab = PchBy,
        size_var = size_var,
        size_lab = PointSize,
        key_var = rownames(data),
        height = input$heightVisu,
        width = if (input$modifwidthVisu) {
          input$widthVisu
        },
        point_opacity = 0.7,
        labels_size = input$SizeLabelScatter,
        transitions = TRUE
      )
      return(plot)
    }
    
    if (export || input$AddRegScatter) {
      if (!lmEst && !CorEst) {
        col_var = if (ColBy == "None" || is.null(ColBy))
          1
        else
          data[, ColBy]
        symbol_var = if (PchBy == "None" ||
                         is.null(PchBy))
          factor(rep(1, nrow(data)))
        else
          data[, PchBy]
        size_var = if (PointSize == "None" ||
                       is.null(PointSize))
          1
        else
          data[, PointSize]
        
        plot = ggplot(data, aes(x = x_var, y = y_var)) + geom_point(aes(
          color = col_var,
          size = size_var,
          shape = symbol_var
        ),
        alpha = 0.7) + theme_bw()
        if (input$AddRegScatter)
          plot = plot + geom_smooth(method = "lm")
        if (input$SizeLabelScatter != 0)
          plot = plot + geom_text(
            aes(
              label = rownames(data),
              color = col_var,
              size = as.numeric(input$SizeLabelScatter) / 10
            ),
            vjust = 0,
            nudge_y = 0.05
          )
        plot = plot + xlab(Xvar) + ylab(Yvar)
        
        return(plot)
      }
    }
    if (lmEst && !CorEst)
    {
      res = lm(y_var ~ x_var)
      sumRes = summary(res)
      regCoef = sumRes$coefficients
      rownames(regCoef) = c("Intercept", Xvar)
      Rsq = sumRes$r.squared
      return(list(regCoef = regCoef, Rsq = Rsq))
    }
    if (CorEst)
    {
      #typesTarget = sapply(target,class)
      #print(typesTarget)
      #numInd = (typesTarget=="numeric")[1:ncol(target)]
      #print(numInd)
      typesTarget = sapply(target, is.numeric)
      numInd = which(typesTarget[2:ncol(target)])
      if (any(numInd))
        data = cbind(target[, numInd], log2(t(counts) + 1), div)
      if (!any(numInd))
        data = cbind(log2(t(counts) + 1), div)
      #typesdata = sapply(data,is.numeric)
      cor.est = round(cor(
        dplyr::select_if(as.data.frame(data), is.numeric),
        method = input$CorMeth
      ), 3)
      #cor.est = round(cor(data[,numInd],method = input$CorMeth),3)
      #cor.pvalue = cor.test(data,method = input$CorMeth)
      return(list(cor.est = cor.est))
    }
  }



##                       ##
##      Diversity ####
##                       ##
Plot_Visu_Diversity <-
  function(input,
           resDiff,
           colors = NULL,
           ForScatter = FALSE,
           alpha_transparency = 0.8) {
    gg = NULL
    dataTmp = NULL
    dataDiv = NULL
    dds = resDiff$dds
    counts = round(counts(dds, normalized = TRUE))
    
    ## Get Input for the plot
    if (!ForScatter)
    {
      VarInt = input$VisuVarInt
      VarIntBoxDiv = input$VarBoxDiv
      VarIntDivCol = input$VarDivCol
      ind_taxo = rownames(counts)
      tmp = GetDataToPlot(input,
                          resDiff,
                          VarInt,
                          ind_taxo,
                          aggregate = FALSE,
                          rarefy = TRUE)
      counts_tmp_combined = tmp$counts
      targetInt = tmp$targetInt
      levelsMod = tmp$levelsMod
    }
    if (ForScatter)
    {
      counts_tmp_combined = t(counts)
      targetInt = resDiff$target
      targetInt$AllVar = targetInt[, 1]
      levelsMod = NULL
    }
    if (nrow(counts_tmp_combined) > 0 &&
        !is.null(counts_tmp_combined) && !is.null(targetInt))
    {
      cond = table(targetInt$AllVar)
      cond = cond[cond != 0]
      sqrt.nb = sqrt(cond)
      
      #save(counts_tmp_combined,targetInt,file = "testDiv.RData")
      alpha <-
        tapply(TaxoNumber(counts_tmp_combined), targetInt$AllVar, mean)
      alpha_selected = !is.na(alpha)
      alpha = alpha[alpha_selected]
      alpha_sd = tapply(TaxoNumber(counts_tmp_combined), targetInt$AllVar, sd)
      alpha_sd = alpha_sd[alpha_selected]
      ci.alpha.down = pmax(alpha - 1.96 * alpha_sd / sqrt.nb, 0)
      ci.alpha.up = alpha + 1.96 * alpha_sd / sqrt.nb
      
      #print(alpha)#mean per condition
      #print(TaxoNumber(counts_tmp_combined)) #mean per species per condition
      
      shan <-
        tapply(
          vegan::diversity(counts_tmp_combined, index = "shannon"),
          targetInt$AllVar,
          mean
        )
      shan_selected = !is.na(shan)
      shan = shan[shan_selected]
      shan_sd = tapply(vegan::diversity(counts_tmp_combined, index = "shannon"),
                       targetInt$AllVar,
                       sd)
      shan_sd = shan_sd[shan_selected]
      ci.shan.down = pmax(shan - 1.96 * shan_sd / sqrt.nb, 0)
      ci.shan.up = shan + 1.96 * shan_sd / sqrt.nb
      
      simpson <-
        tapply(
          vegan::diversity(counts_tmp_combined, index = "simpson"),
          targetInt$AllVar,
          mean
        )
      simpson_selected = !is.na(simpson)
      simpson = simpson[simpson_selected]
      simpson_sd = tapply(vegan::diversity(counts_tmp_combined, index = "simpson"),
                          targetInt$AllVar,
                          sd)
      simpson_sd = simpson_sd[simpson_selected]
      ci.simpson.down = pmax(simpson - 1.96 * simpson_sd / sqrt.nb, 0)
      ci.simpson.up = simpson + 1.96 * simpson_sd / sqrt.nb
      
      invsimpson <-
        tapply(
          vegan::diversity(counts_tmp_combined, index = "invsimpson"),
          targetInt$AllVar,
          mean
        )
      invsimpson_selected = !is.na(invsimpson)
      invsimpson = invsimpson[invsimpson_selected]
      invsimpson_sd = tapply(
        vegan::diversity(counts_tmp_combined, index = "invsimpson"),
        targetInt$AllVar,
        sd
      )
      invsimpson_sd = invsimpson_sd[invsimpson_selected]
      ci.invsimpson.down = pmax(invsimpson - 1.96 * invsimpson_sd / sqrt.nb, 0)
      ci.invsimpson.up = invsimpson + 1.96 * invsimpson_sd / sqrt.nb
      
      braycurtis = tapply(colSums(as.matrix(
        vegan::vegdist(counts_tmp_combined, method = "bray")
      )), targetInt$AllVar, mean)
      braycurtis_selected <- !is.na(braycurtis)
      braycurtis = braycurtis[braycurtis_selected]
      braycurtis_sd = tapply(colSums(as.matrix(
        vegan::vegdist(counts_tmp_combined, method = "bray")
      )),  targetInt$AllVar, sd)
      braycurtis_sd = braycurtis_sd[braycurtis_selected]
      ci.braycurtis.down = pmax(braycurtis - 1.96 * braycurtis_sd / sqrt.nb, 0)
      ci.braycurtis.up = braycurtis + 1.96 * braycurtis_sd / sqrt.nb
      
      gamma <- TaxoNumber(counts_tmp_combined, targetInt$AllVar)
      gamma = gamma[!is.na(gamma)]
      beta = gamma / alpha - 1
      nb = length(alpha)
      dataTmp = data.frame(
        value = c(alpha, beta,
                  braycurtis,
                  gamma, shan, simpson, invsimpson),
        ci.down = c(
          ci.alpha.down,
          beta,
          ci.braycurtis.down,
          gamma,
          ci.shan.down,
          ci.simpson.down,
          ci.invsimpson.down
        ),
        ci.up = c(
          ci.alpha.up,
          beta,
          ci.braycurtis.up,
          gamma,
          ci.shan.up,
          ci.simpson.up,
          ci.invsimpson.up
        ),
        diversity = c(
          rep("Alpha", nb),
          rep("Beta", nb),
          rep("Bray-Curtis", nb),
          rep("Gamma", nb),
          rep("Shannon", nb),
          rep("Simpson", nb),
          rep("Inv.Simpson", nb)
        ),
        Var = as.character(rep(names(alpha), 7))
      )
      richness <- data.frame(
        Var = names(TaxoNumber(counts_tmp_combined)),
        valAlpha = TaxoNumber(counts_tmp_combined),
        valBray = colSums(as.matrix(
          vegan::vegdist(counts_tmp_combined, method = "bray")
        )),
        valShannon = vegan::diversity(counts_tmp_combined, index = "shannon"),
        valSimpson = vegan::diversity(counts_tmp_combined, index = "simpson"),
        valInvSimpson = vegan::diversity(counts_tmp_combined, index = "invsimpson")
      )
      dataDiv = dataTmp
      if (!ForScatter)
      {
        dataTmp = dataTmp[dataTmp$diversity %in% input$WhichDiv, ]
        
        ## Order of the modalities
        dataTmp$Var = factor(dataTmp$Var, levels = levelsMod)
        richness$Var = factor(richness$Var, levels = levelsMod)
        req(dataTmp$Var)
        tmp.mat = matrix(unlist((lapply(
          as.matrix(as.character(dataTmp$Var)), strsplit, "-"
        ))), ncol = length(VarInt), byrow = T)
        tmp.level = matrix(unlist((lapply(
          as.matrix(as.character(levelsMod)), strsplit, "-"
        ))), ncol = length(VarInt), byrow = T)
        richness.mat = matrix(unlist((lapply(
          as.matrix(as.character(richness$Var)), strsplit, "-"
        ))), ncol = length(VarInt), byrow = T)
        
        indVar = VarInt %in% VarIntBoxDiv
        if (length(which(indVar)) >= 1) {
          richness$VarCol = paste0(sub("-.*", "", richness$Var))
          if (length(which(indVar)) >= 2) {
            tmp.levelX = apply(tmp.level[, which(indVar)], 1, paste, collapse = "-")
            dataTmp$VarX = factor(apply(tmp.mat[, which(indVar)], 1, paste, collapse = "-"),
                                  levels = unique(tmp.levelX))
            richness$VarCol = factor(apply(richness.mat[, which(indVar)], 1, paste, collapse = "-"),
                                     levels = unique(tmp.levelX))
          }
          if (length(which(indVar)) == 1) {
            tmp.levelX = tmp.level[, which(indVar)]
            dataTmp$VarX = factor(tmp.mat[, which(indVar)], levels = unique(tmp.levelX))
            richness$VarCol = factor(richness.mat[, which(indVar)], levels = unique(tmp.levelX))
          }
        }
        
        if (is.null(VarIntBoxDiv)) {
          dataTmp$VarX = tmp.mat[, 1]
          richness$VarCol = richness.mat[, 1]
        }
        dataTmp$VarCol = dataTmp$VarX
        if (length(which(!indVar)) >= 1) {
          if (length(which(!indVar)) >= 2) {
            tmp.levelCol = apply(tmp.level[, which(!indVar)], 1, paste, collapse = "-")
            dataTmp$VarCol = factor(apply(tmp.mat[, which(!indVar)], 1, paste, collapse = "-"),
                                    levels = unique(tmp.levelCol))
          }
          if (length(which(!indVar)) == 1) {
            tmp.levelCol = tmp.level[, which(!indVar)]
            dataTmp$VarCol = factor(tmp.mat[, which(!indVar)], levels = unique(tmp.levelCol))
          }
        }
        
        if (input$DiversityPlots == 2) {
          dataDiv <- dataTmp
          gg = ggplot(dataTmp, aes(
            x = Var,
            y = value,
            fill = VarX
          ))
          gg = gg + theme_bw() + theme(
            axis.text.x = element_text(
              angle = 90,
              hjust = 1,
              vjust = 0.5
            ),
            legend.title = element_blank()
          )
          gg = gg + geom_col(
            width = 0.5,
            position = position_dodge(width = 0.5),
            alpha = alpha_transparency
          )
          if (input$SensPlotVisu == "Horizontal")
            gg = gg + coord_flip() + facet_wrap( ~ diversity, scales = "fixed")
          if (input$SensPlotVisu == "Vertical")
            gg = gg + facet_wrap( ~ diversity, scales = input$DivScale)
          gg = gg + xlab(paste(VarIntBoxDiv, collapse = "-")) + ylab("Diversity")
          gg = gg + scale_fill_manual(values = colors[1:length(unique(dataTmp[, 6]))]) + scale_color_manual(values = colors[1:length(unique(dataTmp[, 6]))])
        }
        else{
          if (input$DiversityPlots == 1) {
            dataDiv <- richness
            if (input$WhichDiv == 'Alpha') {
              gg = ggpubr::ggboxplot(
                richness,
                x = "Var",
                y = "valAlpha",
                merge = TRUE,
                color = "VarCol",
                add = c("jitter"),
                width = 0.75
              )
              y_range = c(min(richness$valAlpha), max(richness$valAlpha))
            }
            if (input$WhichDiv == 'Bray-Curtis') {
              gg = ggpubr::ggboxplot(
                richness,
                x = "Var",
                y = "valBray",
                merge = TRUE,
                color = "VarCol",
                add = c("jitter"),
                width = 0.75
              )
              y_range = c(min(richness$valBray), max(richness$valBray))
            }
            else if (input$WhichDiv == 'Shannon') {
              gg = ggpubr::ggboxplot(
                richness,
                x = "Var",
                y = "valShannon",
                merge = TRUE,
                color = "VarCol",
                add = c("jitter"),
                width = 0.75
              )
              y_range = c(min(richness$valShannon),
                          max(richness$valShannon))
            }
            else if (input$WhichDiv == 'Simpson') {
              gg = ggpubr::ggboxplot(
                richness,
                x = "Var",
                y = "valSimpson",
                merge = TRUE,
                color = "VarCol",
                add = c("jitter"),
                width = 0.75
              )
              y_range = c(min(richness$valSimpson),
                          max(richness$valSimpson))
            }
            else if (input$WhichDiv == 'Inv.Simpson') {
              gg = ggpubr::ggboxplot(
                richness,
                x = "Var",
                y = "valInvSimpson",
                merge = TRUE,
                color = "VarCol",
                add = c("jitter"),
                width = 0.75
              )
              y_range = c(min(richness$valInvSimpson),
                          max(richness$valInvSimpson))
            }
            gg = gg + stat_summary(
              fun.y = mean,
              geom = "errorbar",
              aes(ymax = ..y.., ymin = ..y..),
              width = .75,
              linetype = "dashed"
            )
          }
          gg = gg + theme_bw() + theme(
            axis.text.x = element_text(
              angle = 90,
              hjust = 1,
              vjust = 0.5
            ),
            legend.title = element_blank()
          )
          if (input$SensPlotVisu == "Horizontal")
            gg = gg + coord_flip() + geom_jitter(color = "black",
                                                 size = 0.4,
                                                 alpha = 0.9)
          if (input$SensPlotVisu == "Vertical")
            gg = gg + geom_jitter(color = "black",
                                  size = 0.4,
                                  alpha = 0.9)
          req(gg)
          req(richness)
          gg = gg +
            ggpubr::stat_compare_means(
              comparisons = lapply(input$SelectizePairs, function(element)
                strsplit(element, " VS ")[[1]]),
              method = "wilcox.test",
              paired = FALSE,
              p.adjust.method = input$AdjMeth,
              size = input$sizeDiversityTests) +
            ggpubr::stat_compare_means(
              label.x.npc = "center",
              label.y = y_range[1] * 0.9,
              paired = FALSE,
              size = input$sizeDiversityTests
            ) + coord_cartesian()
          gg = gg + labs(
            x = paste(VarIntBoxDiv, collapse = "-"),
            y = paste0("Richness (", input$WhichDiv, ")")
          ) +
            scale_fill_manual(values = colors[1:length(unique(richness[, 7]))]) + scale_color_manual(values = colors[1:length(unique(richness[, 7]))])
        }
        
      }
      ## Get interactivity
      #ff = ggplotly(gg)
      gg <- gg +
        theme(axis.text.x = element_text(size = input$sizeDiversityTitle),  
              axis.text.y = element_text(size = input$sizeDiversityTitle),   
              axis.title.x = element_text(size = input$sizeDiversityTitle * 1.3),  
              axis.title.y = element_text(size = input$sizeDiversityTitle * 1.3),
              legend.text = element_text(size = input$sizeDiversityLegend)) +
        guides(fill = guide_legend(title="Groups"))
    }
    return(list(plot = gg, dataDiv = dataDiv))
    
  }


##                       ##
##       RAREFACTION ####
##                       ##
Plot_Visu_Rarefaction <-
  function(input, resDiff, xlim, ylim, ylab = "Species") {
    PlotRare = NULL
    dds = resDiff$dds
    
    ## Taxo in columns
    
    if (input$NormOrRaw == "norm")
    {
      counts = t(round(counts(dds, normalized = TRUE)))
    }
    else
    {
      counts = t(round(counts(dds, normalized = FALSE)))
    }
    
    if (nrow(counts) > 0 && !is.null(counts))
    {
      max <- max(rowSums(counts))
      raremax <- min(rowSums(counts))
      #PlotRare = rarefaction_curve(counts, step = 10, taxo = "Taxonomy level")
      options(warn = -1)
      PlotRare = rarecurve(
        counts,
        step = max(1, ceiling(max / 60)),
        sample = raremax,
        col = "blue",
        cex = 0.9,
        xlim = xlim,
        ylim = ylim,
        ylab = ylab
      )
      options(warn = 0)
    }
    return(PlotRare)
  }




##                                      ##
##
##          Useful functions ####
##
##                                      ##

## Get the non-zero taxo by sample
TaxoNumber <-  function (x, groups, mar = 1)
{
  if (!missing(groups))
  {
    if (length(groups) == 1)
      groups = rep(groups, nrow(x))
    x = aggregate(x, list(groups), max)
    rownames(x) = x[, 1]
    x = x[,-1]
  }
  if (length(dim(x)) > 1)
    res = apply(x > 0, mar, sum)
  else
    res = sum(x > 0)
  return(res)
}


## Modified version of expand.grid
expand.grid2.list <- function(listInput)
{
  n = length(listInput)
  if (is.list(listInput) && n > 1)
  {
    l1 = listInput[[1]]
    l2 = listInput[[2]]
    res = c()
    
    for (i in l1) {
      for (j in l2) {
        res = rbind(res, paste(i, j, sep = "-"))
      }
    }
    listInput[[1]] = res
    listInput = listInput[-2]
    if (length(listInput) > 1 &&
        is.list(listInput))
      res = expand.grid2.list(listInput)
  }
  else
    res = listInput
  return(res)
}


## Put the data in the right format to be plot
GetDataToPlot <-
  function(input,
           resDiff,
           VarInt,
           ind_taxo,
           sec_variable = NULL,
           aggregate = TRUE,
           rarefy = FALSE)
  {
    sec_variable_added_to_VarInt <- FALSE
    if (!is.null(sec_variable)) {
      if (!is.element(sec_variable, VarInt)) {
        VarInt <- c(VarInt, sec_variable)
        sec_variable_added_to_VarInt <- TRUE
      }
    }
    
    dds = resDiff$dds
    val = c()
    list.val = list()
    if (input$NormOrRaw == "norm")
    {
      counts = as.data.frame(round(counts(dds, normalized = TRUE)))
    }
    if (input$NormOrRaw == "vst") {
      counts = as.data.frame(SummarizedExperiment::assay(varianceStabilizingTransformation(dds)))
    }
    else
    {
      counts = as.data.frame(round(counts(dds, normalized = FALSE)))
    }
    if (rarefy) {
      set.seed(1234)
      counts = t(rrarefy(t(counts), min(colSums(counts))))
    }
    
    target = resDiff$target
    counts_tmp_combined = NULL
    prop_tmp_combined = NULL
    targetInt = NULL
    namesCounts = NULL
    levelsMod = NULL
    prop_all = NULL
    ## Select a subset within the taxonomy level (default is the 12 most abundant)
    nbKept = length(ind_taxo)
    Taxonomy = rownames(counts)
    
    if (length(VarInt) > 0 && nbKept > 0)
    {
      ## Get the modalities to keep
      for (i in 1:length(VarInt))
      {
        ## Replace "-" by "."
        target[, VarInt[i]] =  gsub("-", ".", target[, VarInt[i]])
        
        Tinput = paste("input$", "ModVisu", VarInt[i], sep = "")
        expr = parse(text = Tinput)
        ## All the modalities for all the var of interest
        val = c(val, eval(expr))
        if (sec_variable_added_to_VarInt) {
          mod <- eval(expr)
          if (is.null(mod)) {
            mod <- as.character(unique(as.factor(target[, VarInt[i]])))
          }
          list.val[[i]] = mod
        }
        else{
          list.val[[i]] = eval(expr)
        }
      }
      if (!is.null(val) && !is.null(list.val))
      {
        ## Create the variable to plot
        targetInt = as.data.frame(target[, VarInt])
        rownames(targetInt) = rownames(target)
        ## Combining the Varint
        if (length(VarInt) > 1) {
          targetInt$AllVar = apply(targetInt, 1, paste, collapse = "-")
          targetInt$AllVar = factor(targetInt$AllVar, levels =  expand.grid2.list(list.val))
        }
        if (length(VarInt) <= 1) {
          targetInt$AllVar = target[, VarInt]
          targetInt$AllVar = factor(targetInt$AllVar, levels = val)
        }
        colnames(targetInt) = c(VarInt, "AllVar")
        
        ## Keep only the selected modalities
        ind_kept = which(!is.na(targetInt$AllVar))
        targetInt = targetInt[ind_kept, ]
        
        levelsMod = levels(targetInt$AllVar)
        
        ## Create the counts matrix only for the selected subset
        counts_tmp = counts[Taxonomy %in% ind_taxo, ]
        counts_tmp = counts_tmp[, colnames(counts_tmp) %in% rownames(targetInt)]
        
        ## Proportions over all the taxonomies
        ## Proportion verified
        prop_all = t(counts) / rowSums(t(counts))
        prop_all = as.data.frame(prop_all[, Taxonomy %in% ind_taxo])
        prop_all = as.matrix(prop_all[rownames(prop_all) %in% rownames(targetInt), ])
        rownames(prop_all) = targetInt$AllVar
        
        ## Be careful transposition !
        if (aggregate &&
            nrow(counts_tmp) > 0 && nrow(targetInt) > 0)
        {
          counts_tmp_combined = aggregate(t(counts_tmp), by = list(targetInt$AllVar), mean)
          rownames(counts_tmp_combined) = counts_tmp_combined$Group.1
          namesCounts = counts_tmp_combined$Group.1
          counts_tmp_combined = as.matrix(counts_tmp_combined[, -1])
        }
        else if (!aggregate &&
                 nrow(counts_tmp) > 0 && nrow(targetInt) > 0)
        {
          ## Proportion verified
          counts_tmp_combined = t(counts_tmp)
          prop_tmp_combined = counts_tmp_combined / colSums(counts_tmp)
          rownames(counts_tmp_combined) = targetInt$AllVar
          namesCounts = targetInt$AllVar
          rownames(prop_tmp_combined) = targetInt$AllVar
        }
        
        ## Ordering the counts
        if (!is.null(counts_tmp_combined))
        {
          MeanCounts = apply(counts_tmp_combined, 2, mean)
          ord = order(MeanCounts, decreasing = TRUE)
          if (ncol(counts_tmp_combined) > 1)
            counts_tmp_combined = counts_tmp_combined[, ord]
          if (!aggregate) {
            if (ncol(counts_tmp_combined) > 1)
              prop_tmp_combined = prop_tmp_combined[, ord]
          }
          if (ncol(prop_all) > 1)
            prop_all = prop_all[, ord]
        }
      }
    }
    return(
      list(
        counts = counts_tmp_combined,
        targetInt = targetInt,
        prop = prop_tmp_combined,
        namesCounts = namesCounts,
        levelsMod = levelsMod,
        prop_all = prop_all
      )
    )
    
  }


## Create the data table for the tree representation
CreateTableTree <-
  function(input,
           resDiff,
           CT_Norm_OTU,
           taxo_table,
           VarInt,
           ind_taxo = rownames(CT_Norm_OTU))
  {
    dds = resDiff$dds
    val = c()
    list.val = list()
    counts = CT_Norm_OTU
    target = resDiff$target
    counts_tmp_combined = NULL
    prop_tmp_combined = NULL
    targetInt = NULL
    namesCounts = NULL
    levelsMod = NULL
    ## Select a subset within the taxonomy level (default is the 12 most abundant)
    nbKept = length(ind_taxo)
    Taxonomy = rownames(counts)
    
    if (length(VarInt) > 0 && nbKept > 0)
    {
      ## Get the modalities to keep
      for (i in 1:length(VarInt))
      {
        ## Replace "-" by "."
        target[, VarInt[i]] =  gsub("-", ".", target[, VarInt[i]])
        
        Tinput = paste("input$", "ModVisu", VarInt[i], sep = "")
        expr = parse(text = Tinput)
        ## All the modalities for all the var of interest
        val = c(val, eval(expr))
        list.val[[i]] = eval(expr)
      }
      if (!is.null(val) && !is.null(list.val))
      {
        # save(val,list.val,VarInt,target,Taxonomy,counts,ind_taxo,file="testTree.RData")
        
        ## Create the variable to plot
        targetInt = as.data.frame(target[, VarInt])
        rownames(targetInt) = rownames(target)
        ## Combining the Varint
        if (length(VarInt) > 1) {
          targetInt$AllVar = apply(targetInt, 1, paste, collapse = "-")
          targetInt$AllVar = factor(targetInt$AllVar, levels =  expand.grid2.list(list.val))
        }
        else{
          targetInt$AllVar = target[, VarInt]
          targetInt$AllVar = factor(targetInt$AllVar, levels = val)
        }
        colnames(targetInt) = c(VarInt, "AllVar")
        
        ## Keep only the selected modalities
        ind_kept = which(!is.na(targetInt$AllVar))
        targetInt = targetInt[ind_kept, ]
        
        levelsMod = levels(targetInt$AllVar)
        
        ## Create the counts matrix only for the selected subset
        counts_tmp = counts[Taxonomy %in% ind_taxo, ]
        counts_tmp = counts_tmp[, colnames(counts_tmp) %in% rownames(targetInt)]
        
        ## Be careful transposition !
        # Group per condition
        if (nrow(counts_tmp) > 0 && nrow(targetInt) > 0)
        {
          counts_tmp_combined = aggregate(t(counts_tmp), by = list(targetInt$AllVar), mean)
          namesCounts = counts_tmp_combined$Group.1
          rownames(counts_tmp_combined) = namesCounts
          
          counts_tmp_combined = as.matrix(counts_tmp_combined[, -1])
        }
        
        
        ## Ordering the counts
        #       if(!is.null(counts_tmp_combined))
        #       {
        #         MeanCounts = apply(counts_tmp_combined,2,mean)
        #         ord = order(MeanCounts,decreasing=TRUE)
        #         counts_tmp_combined = as.matrix(counts_tmp_combined[,ord])
        #
        #       }
        
      }
    }
    
    return(
      list(
        counts = counts_tmp_combined,
        targetInt = targetInt,
        namesCounts = namesCounts,
        levelsMod = levelsMod
      )
    )
    
    
  }




##                       ##
##        Tree
##                       ##

## The count matrix must be given at the leaf level.

Plot_Visu_Tree <- function(input, resDiff, CT_Norm_OTU, taxo_table)
{
  res = NULL
  ## Get Input for BarPlot
  VarInt = input$VisuVarInt
  ind_taxo = input$selectTaxoPlot
  nodeFind = input$TaxoTree
  ## Removed column with only 1 modality
  ind = which(apply(
    taxo_table,
    2,
    FUN = function(x)
      length(unique(x[!is.na(x)]))
  ) == 1)
  if (length(ind) > 0)
    taxo_table = taxo_table[, -ind]
  # tmp_combined = GetDataToPlot(input,resDiff,VarInt,ind_taxo,CT_Norm_OTU=CT_Norm_OTU)
  
  if (nrow(CT_Norm_OTU) > 0 &&
      !is.null(CT_Norm_OTU) &&
      nrow(taxo_table) > 0 && !is.null(taxo_table))
  {
    tmp = CreateTableTree(input, resDiff, CT_Norm_OTU, taxo_table, VarInt)
    
    if (nrow(tmp$counts) > 0 && !is.null(tmp$counts))
    {
      #save(tmp,taxo_table,nodeFind,file="testTree.RData")
      merge_dat = merge(cbind(taxo_table, rownames(taxo_table)), round(t(tmp$counts)), by =
                          "row.names")
      colnames(merge_dat)[1] = "Category"
      colnames(merge_dat)[dim(taxo_table)[2] + 2] = "OTU"
      merge_dat[is.na(merge_dat)] = ""
      levels <- c("Category", colnames(taxo_table), "OTU")
      conditions <- rownames(tmp$counts)
      #nodeFind = input$TaxoTree
      nodeFind = input$selectTaxoPlot
      if (length(input$selectTaxoPlot) == 0)
        nodeFind = NULL
      #save(merge_dat,conditions,levels,nodeFind,file="abuntree.RDATA")
      res = treeWeightD3(
        merge_dat,
        conditions,
        levels,
        nodeFind = nodeFind,
        height = input$heightVisu + 10,
        width = if (input$modifwidthVisu) {
          input$widthVisu
        },
        input$isColorRandom,
        input$colorLevel
      )
    }
  }
  return(res)
}


##                       ##
##      NETWORK
##                       ##

#Function that computes a permutation test and return an adjacency matrix
#input: mat, max iteration number
#output: adjacency matrix

# compute_pcor <- function(input, mat, n_iter = 100) {
#   
#   n_rows <- nrow(mat)
#   correlation_storage <- vector("list", length = n_rows)
#   for (i in 1:n_rows) {
#     correlation_storage[[i]] <- numeric(n_iter * (n_rows - 1))
#   }
#   adjacency_matrix <- matrix(0, nrow = n_rows, ncol = n_rows)
#   correlation_matrix <- matrix(NA, nrow = n_rows, ncol = n_rows)
#   
#   # Loop through each pair of rows
#   for (i in 1:(n_rows - 1)) {
#     for (j in (i + 1):n_rows) {
#       # Initialize vectors to store correlations between rows i and j for each shuffle
#       correlations_i <- numeric(n_iter)
#       correlations_j <- numeric(n_iter)
#       
#       # Shuffle and compute correlations n_iter times
#       for (shuffle in 1:n_iter) {
#         # Shuffle only the i-th and j-th row
#         mat_shuffled_i <- mat
#         mat_shuffled_j <- mat
#         mat_shuffled_i[i, ] <- sample(mat[i, ])
#         mat_shuffled_j[j, ] <- sample(mat[j, ])
#         
#         # Compute correlation of the shuffled i-th row with row j and vice versa
#         correlations_i[shuffle] <- cor(mat_shuffled_i[i, ], mat[j, ])
#         correlations_j[shuffle] <- cor(mat_shuffled_j[j, ], mat[i, ])
#       }
#       
#       # Store the shuffled correlations
#       correlation_storage[[i]][((j - 1) * n_iter + 1):(j * n_iter)] <- correlations_i
#       correlation_storage[[j]][((i - 1) * n_iter + 1):(i * n_iter)] <- correlations_j
#       
#       # Compute the observed correlation once for each unique pair of rows
#       observed_correlation <- cor(mat[i, ], mat[j, ])
#       correlation_matrix[i, j] <- observed_correlation
#       correlation_matrix[j, i] <- observed_correlation  # symmetry
#       
#       # Calculate p-values and fill the adjacency matrix using the observed_correlation
#       p_value_i <- min(sum(correlations_i < observed_correlation), sum(correlations_i > observed_correlation)) / n_iter
#       p_value_j <- min(sum(correlations_j < observed_correlation), sum(correlations_j > observed_correlation)) / n_iter
#       
#       # Determine significance based on the p-value and the sign of the observed correlation
#       if ((p_value_i <= 0.025 || p_value_i >= 0.975) && observed_correlation > 0) {
#         adjacency_matrix[i, j] <- 1
#         adjacency_matrix[j, i] <- 1
#       } else if ((p_value_i <= 0.025 || p_value_i >= 0.975) && observed_correlation < 0) {
#         adjacency_matrix[i, j] <- -1
#         adjacency_matrix[j, i] <- -1
#       }
#     }
#   }
#   
#   # Return the adjacency matrix and the correlation matrix
#   return(list(adjacency_matrix = adjacency_matrix, correlation_matrix = correlation_matrix))
# }


# Register the parallel backend
# numCores <- detectCores() - 1
# doParallel::registerDoParallel(cores=numCores)
# # Function to compute partial correlations
# compute_partial_correlations <- function(mat, n_iter, method, threshold) {
#   n_rows <- nrow(mat)
#   corr_storage <- vector("list", length = n_rows)
#   adjacency <- matrix(0, nrow = n_rows, ncol = n_rows)
#   corr_matrix <- matrix(NA, nrow = n_rows, ncol = n_rows)
#   
#   # Perform parallel computation
#   results <- foreach(i = 1:(n_rows - 1), .combine = rbind, .multicombine = TRUE, .packages = c("stats")) %dopar% {
#     local_corrs <- numeric(n_iter)
#     for (shuffle in 1:n_iter) {
#       mat_shuffled <- mat
#       mat_shuffled[i, ] <- sample(mat[i, ])
#       local_corrs[shuffle] <- cor(mat_shuffled[i, ], mat[j, ], method = method)
#     }
#     
#     # Compute the observed correlation for the i-th row
#     observed_corr <- sapply((i + 1):n_rows, function(j) cor(mat[i, ], mat[j, ], method = method))
#     p_values <- sapply((i + 1):n_rows, function(j) {
#       local_pval <- min(sum(local_corrs < observed_corr[j-i]), sum(local_corrs > observed_corr[j-i])) / n_iter
#       local_pval
#     })
#     
#     cbind(i, (i + 1):n_rows, observed_corr, p_values)
#   }
#   
#   # Process results
#   for (result in results) {
#     i <- result[, 1]
#     j <- result[, 2]
#     observed_corr <- result[, 3]
#     p_value <- result[, 4]
#     
#     # Fill the adjacency matrix
#     adjacency[i, j] <- ifelse(p_value <= threshold/2 || p_value >= 1 - threshold/2, sign(observed_corr), 0)
#     corr_matrix[i, j] <- observed_corr
#   }
#   
#   # Ensure symmetry
#   adjacency <- adjacency + t(adjacency)
#   corr_matrix <- corr_matrix + t(corr_matrix)
#   
#   list(adjacency = adjacency, corr_matrix = corr_matrix)
# }
# 
# # Stop the parallel cluster when you are done
# doParallel::stopImplicitCluster()

Plot_network <-
  function(input,
           resDiff,
           compute_pcor,
           availableTaxo,
           ind_taxo,
           qualiVariable,
           dataInput,
           export = FALSE,
           colors = NULL) {
    plot = NULL
    dataVN = NULL
    cluster = NULL
    cluster_df = NULL
    voronoi = NULL
    pcor_mat= NULL
    community_barplot = NULL
    VarInt = input$VisuVarInt
    
    set.seed(08)
    if (isolate(input$colorCorr == 'corr')) {
      sec_variable = isolate(input$sec_variable)
    }
    else{
      sec_variable = NULL
    }
    
    data <-
      GetDataToPlot(
        input,
        resDiff,
        VarInt,
        availableTaxo,
        sec_variable = sec_variable,
        aggregate = FALSE
      )
    if (!is.null(data) && !is.null(data$targetInt)) {
      counts_tmp_combined <- data$counts
      dataVariables <- as.matrix(data$targetInt)
      if (isolate(input$colorCorr == 'corr') &&
          isolate(qualiVariable()) &&
          !is.null(dataVariables)) {
        dataVariables[, sec_variable] <-
          sapply(dataVariables[, sec_variable], function(x)
            if (is.element(x, isolate(input$values1))) {
              1
            } else{
              0
            })
      }
      
      if (!is.null(counts_tmp_combined)) {
        #################################
        ###### partial correlation ######
        #################################
        
        countsMatrix <- as.matrix(counts_tmp_combined)
        # ppcor <-
        #   ppcor::pcor(countsMatrix, method = input$pcorrMethod)
        # pcor <- ppcor$estimate
        # rownames(pcor) <- colnames(countsMatrix)
        # colnames(pcor) <- colnames(countsMatrix)
        #
        # # Step 1: Calculate the variance-covariance matrix
        # covMatrix <- cov(countsMatrix)
        # 
        # # Step 2: Compute the inverse
        # invCovMatrix <- MASS::ginv(covMatrix)
        # 
        # 
        # # Step 3: Calculate partial correlations
        # pcorMatrix <- -invCovMatrix / sqrt(outer(diag(invCovMatrix), diag(invCovMatrix)))
        # 
        # # Set diagonal to zero (self-correlation is not meaningful)
        # diag(pcorMatrix) <- 0
        # 
        # # Step 4: Calculate p-values
        # t_values <- pcorMatrix * sqrt((n - 2) / (1 - pcorMatrix^2))
        # p_values <- 2 * pt(-abs(t_values), df = n - 2)
        # 
        
        ################################
        ####### adjacency matrix #######
        ################################
        
        
        #if (input$colorCorr == "pcorr") {
          #   adjacency <- matrix(apply(p_values, c(1, 2), function(p_val, pcor_val) {
          #     if (!is.na(p_val) && p_val <= 0.05) {
          #       return(1) 
          #     } else {
          #       return(0)
          #     }
          #   }, pcor_val = pcor), nrow = nrow(p_values), ncol = ncol(p_values))
          
          # adjacency <- matrix(apply(pcor, c(1, 2), function(x) {
          #   # Check for NA values in x and ensure input$pcorrThreshold is not NA
          #   if (!is.na(x) && !is.na(input$pcorrThreshold) && abs(x) > as.numeric(input$pcorrThreshold)) {
          #     return(x)
          #   } else {
          #     return(0)
          #   }
          # }), nrow = nrow(pcor), ncol = ncol(pcor))
          
          permutation <- compute_pcor
          adjacency <- permutation$adjacency
          rownames(adjacency) <- colnames(countsMatrix)
          colnames(adjacency) <- colnames(countsMatrix)
          observed_correlation <- permutation$corr_matrix
        #}
        # else{
        #   adjacency <-
        #     matrix(mapply(function(a, b) {
        #       mapply(function(x, y) {
        #         x * y
        #       }, x = a, y = b)
        #     }, a = pval_bool, b = cor_sgn), nrow = n)
        #   rownames(adjacency) <- colnames(countsMatrix)
        #   colnames(adjacency) <- colnames(countsMatrix)
        # }
        # ### Remove rows and columns with only NA        # this way, elements with the same count in all sample (often 0 in this case) will not appear
        # adjacency <- adjacency[apply(adjacency, 1, function(y) !all(is.na(y))),]
        # adjacency <- t(adjacency)
        # adjacency <- adjacency[apply(adjacency, 1, function(y) !all(is.na(y))),]
        # adjacency <- t(adjacency)
        
        ### Replace NA by zeros (ie "no correlation")     # this way, those elements will appear as single nodes
        adjacency[is.na(adjacency)] <- 0
        #adjacency <- abs(adjacency)
        req(adjacency)
        adjacency <- adjacency[, ind_taxo]
        adjacency <- adjacency[ind_taxo, ]
        
        igraphGraph <-
          graph_from_adjacency_matrix(adjacency,
                                      diag = FALSE,
                                      mode = "upper",
                                      weighted = TRUE) # mode = "upper" for adjusted p-value, mode = "lower" for p-value not adjusted
        
        #################################
        ######## nodes and edges ########
        #################################
        
        
        list_to_label <- isolate(input$ToLabelNetwork)
        req(igraphGraph)
        dataVN <- toVisNetworkData(igraphGraph)
        dataVN$nodes$title <- paste0("<b>", dataVN$nodes$id, "</b>")
        dataVN$nodes$label <- dataVN$nodes$id
        # dataVN$nodes$label <- sapply(dataVN$nodes$id, function(x)if(is.element(x, list_to_label)){x}else{""})
        dataVN$nodes$cor <- 0
        dataVN$edges$pcor <- 0
        #Find the linked nodes and set them a color according to their partial correlation
        if (isolate(input$colorCorr == 'pcorr')) {
          
          # Extract unique nodes from the edges
          nodesPcor <- unique(c(dataVN$edges$from, dataVN$edges$to))
          
          # Create a mapping from node names to their indices based on the column names of 'pcor'
          node_mapping <- match(nodesPcor, colnames(countsMatrix))
          
          # Initialize a vector to store partial correlation values
          pcor_values <- numeric(nrow(dataVN$edges))
          
          # Loop over the edges to get the existing partial correlation values
          for (i in seq_len(nrow(dataVN$edges))) {
            from_index <- node_mapping[which(nodesPcor == dataVN$edges$from[i])]
            to_index <- node_mapping[which(nodesPcor == dataVN$edges$to[i])]
            
            # Store the partial correlation value
            pcor_values[i] <- observed_correlation[from_index, to_index]
          }
          
          # Store the correlation values back into the 'pcor' column
          dataVN$edges$pcor <- pcor_values
          
          # Assign colors based on the correlation: Red for positive, Blue for negative
          dataVN$edges$color <- ifelse(dataVN$edges$pcor >= 0, 
                                       toString(isolate(input$edgeColorPositive)),
                                       toString(isolate(input$edgeColorNegative)))
          
          
          
          # dataVN$edges$color <-
          #   ifelse(dataVN$edges$pcor > 0,
          #          toString(isolate(input$edgeColorPositive)),
          #          toString(isolate(input$edgeColorNegative)))
        }
        else{
          dataVN$edges$color <- "#000000"
          dataVN$edges$pcor <- 0
          dataVN$edges$weights <- 0
        }
        
        dataVN$edges$weight <- dataVN$edges$pcor
        if (!is.null(sec_variable) &&
            (isolate(input$colorCorr) == "corr")) {
          cor <-
            sapply(dataVN$nodes$id, function(x)
              cor(as.numeric(dataVariables[, sec_variable]), countsMatrix[, x]))
          dataVN$nodes$cor <- round(cor, digits = 5)
          scale <-
            if (isolate(input$scaleFree)) {
              max(c(max(dataVN$nodes$cor), -min(dataVN$nodes$cor)))
            } else{
              1
            }
          dataVN$nodes$color.background <-
            sapply(dataVN$nodes$cor, function(x)
              colorRampPalette(rev(brewer.pal(
                9, isolate(input$colorPalette)
              )))(200)[round(x, digits = 2) * 100 / scale + 100])
          dataVN$nodes$color.highlight.background <-
            dataVN$nodes$color.background
        }
        palette <- isolate(input$colorsdiagVisuPlot)
        dataVN$edges$weight <- abs(dataVN$edges$weight)
        #################################
        ########### clustering ##########
        #################################
        
        
        if (input$pcorrClustAlgo == "louvain") {
          if (input$clusterWeightsParam == "Yes" && input$NormOrRaw != "vst")
            cluster <-
              igraph::cluster_louvain(igraphGraph, weights = NULL)
          else
            cluster <-
              igraph::cluster_louvain(igraphGraph, weights = NA)
        }
        
        
        else if (input$pcorrClustAlgo == "leiden") {
          if (input$clusterWeightsParam == "Yes" && input$NormOrRaw != "vst")
            cluster <-
              igraph::cluster_leiden(
                igraphGraph,
                weights = NULL,
                objective_function = "modularity",
                n_iterations = -1
              )
          else
            cluster <-
              igraph::cluster_leiden(
                igraphGraph,
                weights = NA,
                objective_function = "modularity",
                n_iterations = -1
              )
        }
        
        
        else if (input$pcorrClustAlgo == "walktrap") {
          if (input$clusterWeightsParam == "Yes" && input$NormOrRaw != "vst")
            cluster <-
              igraph::cluster_walktrap(igraphGraph, weights = NULL)
          else
            cluster <-
              igraph::cluster_walktrap(igraphGraph, weights = NA)
        }
        
        #Let's compute the output to shape the nodes and edges
        cluster_df <-
          data.frame(as.list(membership(cluster)), check.names = FALSE)
        cluster_df <- as.data.frame(t(cluster_df))
        rownames(cluster_df) <- rownames(dataVN$nodes)
        cluster_df$id <- rownames(cluster_df)
        
        dataVN$nodes <-
          dplyr::left_join(dataVN$nodes, cluster_df, by = "id")
        
        names(dataVN$nodes)[names(dataVN$nodes) == "V1"] <-
          "community"
        
        if (!is.na(isolate(input$SelectTaxoNetwork)) && !is.na(isolate(input$SelectTaxoNetwork)) != "Linked") {
          # Nodes that are not in edges
          nodes_not_in_edges <- dataVN$nodes$id[!(dataVN$nodes$id %in% c(dataVN$edges$from, dataVN$edges$to))]
          
          # Check if there are any such nodes and if the vector is not NA
          if (!is.na(nodes_not_in_edges) && length(nodes_not_in_edges) > 0) {
            dataVN$nodes[!(dataVN$nodes$id %in% c(dataVN$edges$from, dataVN$edges$to)),]$community <- 0
            dataVN$nodes$community <- as.integer(factor(dataVN$nodes$community))
          }
        }
        if (isolate(input$colorCorr != 'corr'))
        {
          dataVN$nodes$color.background <-
            colors[1 + (dataVN$nodes$community - 1) %% length(colors)]
          dataVN$nodes$color.highlight.background <-
            dataVN$nodes$color.background
        }
        
        
        dataVN$nodes$color.border <- isolate(input$colorBorder)
        dataVN$nodes$color.highlight.border <-
          isolate(input$colorHighlightBorder)
        
        dataVN$nodes <-
          dataVN$nodes[order(dataVN$nodes$community),]
        dataVN$edges[which(!dataVN$edges$from %in% dataVN$nodes),]
        dataVN$edges[which(!dataVN$edges$to %in% dataVN$nodes),]
        
        
        plot <-
          visNetwork(nodes = dataVN$nodes, edges = dataVN$edges)
        plot <-
          visIgraphLayout(plot,
                          layout = "layout_nicely",
                          physics = FALSE,
                          smooth = FALSE)
        plot <- visLegend(plot, 
                          useGroups = T)
        plot <- visNodes(plot, size = 20)
        plot <-
          visEdges(plot,
                   width = 1,
                   color = list(color = dataVN$edges$color))
        plot <- visGroups(plot, groupname = "community")
        
        plot <-
          visOptions(
            plot,
            width = if (isolate(input$modifwidthVisu)) {
              isolate(input$widthVisu)
            },
            height = isolate(input$heightVisu),
            autoResize = FALSE,
            highlightNearest = list(
              enabled = TRUE,
              degree = 1,
              hover = FALSE
            )
          )
        if (isolate(input$colorCorr) != 'corr')
          plot <-
          visOptions(
            plot,
            selectedBy = list(
              variable = "community",
              style = 'width: 200px; height: 26px; background: #f8f8f8; color: black; font: caption; border:none; outline:none;',
              values = c(
                min(dataVN$nodes$community):max(dataVN$nodes$community)
              )
            )
          )
        plot <- visLayout(plot, randomSeed = 22)
        
        # #Get a new graph from computed nodes and edges
        # graph_igraph <-
        #   graph_from_data_frame(dataVN$edges,
        #                         directed = FALSE,
        #                         vertices = dataVN$nodes)
        # V(graph_igraph)$label.cex = 0.75
        # membership <- cluster$membership
        # 
        # #Find only lonely nodes in the membership list (they have a unique integer) then set it to 0
        # membership[!(membership %in% membership[duplicated(membership)])] <-
        #   0
        # cluster$membership <-
        #   as_membership(as.integer(factor(membership)))
        # 
        # V(igraphGraph)$community <- membership(cluster)
        # if (isFALSE(input$showLabelNetwork)) {
        #   if (input$SelectTaxoNetwork == "Linked") {
        #     voronoi <-
        #       plot(
        #         cluster,
        #         graph_igraph,
        #         mark.groups = communities(cluster),
        #         mark.border = "black",
        #         #vertex.color = V(igraphGraph)$community ,
        #         edge.color = dataVN$edges$color,
        #         vertex.size = input$nodeSizeNetwork / 2,
        #         edge.arrow.size = input$linkWidth / 10,
        #         vertex.label = NA,
        #         main = paste0(
        #           "Community detection using ",
        #           input$pcorrClustAlgo,
        #           " algorithm for ",
        #           input$TaxoSelect
        #         ),
        #         layout = igraph::layout_nicely(graph_igraph)
        #       )
        #   }
        #   else{
        #     voronoi <-
        #       plot(
        #         graph_igraph,
        #         mark.groups = communities(cluster),
        #         mark.border = "black",
        #         #vertex.color = dataVN$nodes$color.background,
        #         edge.color = dataVN$edges$color,
        #         vertex.size = input$nodeSizeNetwork / 2,
        #         edge.arrow.size = input$linkWidth / 10,
        #         vertex.label = NA,
        #         main = paste0(
        #           "Community detection using ",
        #           input$pcorrClustAlgo,
        #           " algorithm for ",
        #           input$TaxoSelect
        #         ),
        #         layout = igraph::layout_nicely(graph_igraph)
        #       )
        #   }
        # }
        # else{
        #   if (input$SelectTaxoNetwork == "Linked") {
        #     voronoi <-
        #       plot(
        #         cluster,
        #         graph_igraph,
        #         mark.groups = communities(cluster),
        #         mark.border = "black",
        #         #vertex.color = dataVN$nodes$color.background,
        #         edge.color = dataVN$edges$color,
        #         vertex.label.color = "black",
        #         vertex.label.dist = 1.5,
        #         vertex.size = input$nodeSizeNetwork / 2,
        #         edge.arrow.size = input$linkWidth / 10,
        #         vertex.label = V(graph_igraph)$label,
        #         main = paste0(
        #           "Community detection using Louvain algorithm for ",
        #           input$TaxoSelect
        #         ),
        #         layout = igraph::layout_nicely(graph_igraph)
        #       )
        #   }
        #   else{
        #     voronoi <-
        #       plot(
        #         graph_igraph,
        #         mark.groups = communities(cluster),
        #         mark.border = "black",
        #         vertex.color = dataVN$nodes$color.background,
        #         edge.color = dataVN$edges$color,
        #         vertex.label.color = "black",
        #         vertex.label.dist = 1.5,
        #         vertex.size = input$nodeSizeNetwork / 2,
        #         edge.arrow.size = input$linkWidth / 10,
        #         vertex.label = V(graph_igraph)$label,
        #         main = paste0(
        #           "Community detection using Louvain algorithm for ",
        #           input$TaxoSelect
        #         ),
        #         layout = igraph::layout_nicely(graph_igraph)
        #       )
        #   }
        # }
        
        #################################
        ###### taxonomy data table ######
        #################################
        
        taxo <- data.frame(dataInput$data$taxo)
        colnames(cluster_df) <- c("Community", as.character(input$TaxoSelect))
        taxo <- data.frame(RowName = rownames(taxo), taxo)
        # Perform data transformation
        taxo <- dplyr::inner_join(taxo, cluster_df, by = as.character(input$TaxoSelect)) %>%
          dplyr::mutate(across(where(is.factor), as.character)) %>% # Convert factors to characters
          tidyr::replace_na(setNames(lapply(taxo, function(x) "Unknown"), names(taxo))) %>%
          dplyr::arrange(Community)
        
        # Reassign the column 'RowName' as row names
        rownames(taxo) <- taxo$RowName
        # Optionally, remove the 'RowName' column if it's no longer needed
        taxo <- taxo %>% dplyr::select(-RowName)
        # Convert back to data frame if needed
        resTaxo <- as.data.frame(taxo)
        
        
      }
      pcor_mat = pcor
    }
    req(pcor_mat)
    return(list(
      plot = plot,
      data = dataVN,
      #voronoi = voronoi,
      pcor_mat = pcor, 
      taxo = resTaxo,
      pcor = pcor
    ))
  }


Plot_network_sunburst <- function(input,
                                  resDiff,
                                  compute_pcor,
                                  availableTaxo,
                                  ind_taxo,
                                  qualiVariable, 
                                  dataInput,
                                  colors){
  res = isolate(Plot_network(input,
                             resDiff,
                             compute_pcor,
                             availableTaxo,
                             ind_taxo,
                             qualiVariable,
                             dataInput,
                             colors = colors)$taxo)
  req(input$CommunitySunburst)
  req(res)
  # Filter for the first community
  res_community <- res %>%
    dplyr::filter(Community == as.integer(input$CommunitySunburst)) %>% 
    dplyr::select(-Community)
  
  # Replace hyphens in taxonomy names to avoid issues with sunburstR paths
  res_community <- res_community %>%
    dplyr::mutate(across(everything(), ~ gsub("-", "_", .)))
  
  # Create the path string for sunburstR from the taxonomy columns
  taxonomy_cols <- setdiff(names(res), "Community")
  res_community$path <- do.call(paste, c(res_community[taxonomy_cols], sep = "-"))
  
  # Count the occurrences of each unique path
  res_counts <- res_community %>%
    dplyr::count(path) %>%
    dplyr::rename(count = n)
  
  # Create the sunburst plot data frame
  sunburst_data <- data.frame(paths = res_counts$path, counts = res_counts$count)
  
  # Generate the sunburst plot for the first community
  p <- sunburstR::sunburst(sunburst_data, count = TRUE, legend = 0, colors = colors, height = ifelse(
    is.na(input$heightVisu),
    toString(paste0(800, "px")),
    toString(paste0(as.character(input$heightVisu), "px"))
  ),
  width = if (input$modifwidthVisu) {
    as.character(input$widthVisu)
  })
  
  # Return the sunburst plot for the first community
  return(p)
  
}
