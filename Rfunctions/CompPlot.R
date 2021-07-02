#@ This file contains all the functions for the 
#@ comparison plots of SHAMAN


##############################
##       HEATMAP
##############################
Plot_Visu_Heatmap_FC <- function(input, BaseContrast, resDiff, ContrastListDebounce, SelectTaxoPlotCompDebounce, export=FALSE){
  
  res = NULL
  #SelContrast = input$ContrastList_table_FC
  SelContrast = ContrastListDebounce()
  selcontrast_matrix = as.matrix(BaseContrast[,SelContrast])
  colnames(selcontrast_matrix) = SelContrast
  log2FC = Get_log2FC_padj(input,selcontrast_matrix,resDiff, info = NULL)$log2FC
  
  if(!is.null(log2FC) && length(SelContrast)>=2)
  { 
    cont = which(colnames(log2FC)%in%SelContrast)
    log2FC = log2FC[,SelContrast]
    ###
    ind_taxo = SelectTaxoPlotCompDebounce()
    #ind_taxo = input$selectTaxoPlotComp
    ind = rownames(log2FC)%in%ind_taxo
    log2FC = as.matrix(log2FC[ind,])
    
    if(input$SortHeatComp =="Selection") tmp_ord = match(ind_taxo, rownames(log2FC))
    if(input$SortHeatComp =="Names") tmp_ord = order(rownames(log2FC))
    if(input$SortHeatComp =="Values") tmp_ord = order(log2FC[,1])
    
    if(input$SortHeatComp !="Auto") log2FC = log2FC[tmp_ord,]
    
    col1 <- c(colorRampPalette(c("royalblue4","royalblue3","royalblue2","royalblue1","white"))(n = 100),colorRampPalette(c("white",  "firebrick1", "firebrick2", "firebrick3", "firebrick4"))(n = 100))
    breaks <- c(seq(min(log2FC,-0.01,na.rm = TRUE), 0,length=100),seq(0.01,max(log2FC,0.02,na.rm = TRUE),length=100))
    colorFunc <- col_bin(col1, bins = rescale(breaks))
    ## Transpose matrix if Horizontal
    if(input$SensPlotVisuComp=="Horizontal") log2FC = t(as.matrix(log2FC))
    
    if(!export && nrow(log2FC)>0) res = d3heatmap(log2FC, dendrogram = "none", Rowv = (input$SortHeatComp =="Auto"), Colv = FALSE, na.rm = TRUE, height = input$heightVisuComp, show_grid = FALSE, colors = colorFunc, scale = input$scaleHeatmapComp,cexRow = input$LabelSizeHeatmapComp,cexCol =input$LabelSizeHeatmapComp, offsetCol=input$LabelColOffsetHeatmapComp,offsetRow=input$LabelRowOffsetHeatmapComp)
    if(export && nrow(log2FC)>0) heatmap.2(log2FC, dendrogram = "none", Rowv = (input$SortHeatComp =="Auto"), Colv = FALSE, na.rm = TRUE, margins=c(input$lowerMarginComp,input$rightMarginComp), density.info="none", trace="none", col = col1, scale = input$scaleHeatmapComp,cexRow = input$LabelSizeHeatmapComp,cexCol =input$LabelSizeHeatmapComp, 
                                           offsetCol=input$LabelColOffsetHeatmapComp,offsetRow=input$LabelRowOffsetHeatmapComp,symm=FALSE,symkey=TRUE,symbreaks=TRUE)
  }
  return(res)
}

##############################
##      P VALUE DENSITY PLOT
##############################
Plot_pValue_Density <- function(input, BaseContrast, resDiff, ContrastListDebounce, alphaVal, InputpValueDensityfocus){
  res = NULL
  #SelContrast = input$ContrastList_table_FC
  SelContrast = ContrastListDebounce()
  if(length(SelContrast)>=1){
  selcontrast_matrix = as.matrix(BaseContrast[,SelContrast])
  colnames(selcontrast_matrix) = SelContrast
  
  padj = data.frame(na.omit(Get_log2FC_padj(input,selcontrast_matrix,resDiff)$padj))
  
  data <- data.frame()
  for (cont in SelContrast){
    data_cont <- data.frame(padj[,cont])
    data_cont$contrast <- cont
    colnames(data_cont) <- c("padj","contrast")
    data <- rbind(data, data_cont)
  }
  data$contrast <- factor(data$contrast, levels = SelContrast)
  p <- ggplot(data, aes(x = padj, color = contrast, fill = contrast)) + theme_minimal()
  InputpValueDensityfocus() # for reactivity
  if(isolate(input$adaptBWtoFocus)){p <- p + geom_density(alpha = input$fillOpacity, size = input$lineWidth, adjust = as.numeric(alphaVal), n = 2 ^ (input$numberPoints))}else{p <- p + geom_density(alpha = input$fillOpacity, size = input$lineWidth, adjust = 1, n = 2 ^ (input$numberPoints))}
  if(isolate(input$focusUnderThreshold)){p <- p + coord_cartesian(xlim = c(0,as.numeric(alphaVal)))}
  p <- p + theme(axis.title = element_text(size = input$FontSizepValueDensity), 
                 axis.text = element_text(size = input$FontSizepValueDensity),
                 legend.title = element_text(size = input$FontSizepValueDensity),
                 legend.text = element_text(size = input$FontSizepValueDensity))
  res = p
  }
  return(res)
}

##############################
##      LOGIT PLOT
##############################
Plot_Comp_Logit <- function(input, BaseContrast, resDiff, SelectTaxoPlotCompDebounce, export = FALSE)
  {
  plot = NULL
  if(input$Contrast1 != input$Contrast2 && input$Contrast1 != "..." && input$Contrast2 != "..."){
    SelContrast = c(input$Contrast1, input$Contrast2)
    selcontrast_matrix = as.matrix(BaseContrast[,SelContrast])
    colnames(selcontrast_matrix) = SelContrast
    padj = na.omit(Get_log2FC_padj(input,selcontrast_matrix,resDiff)$padj)
   
    #data to plot
    x <- padj[,input$Contrast1]
    logit_x <- log(x/(1-x))
    y <- padj[,input$Contrast2]
    logit_y <- log(y/(1-y))
    
    #classification for colours
    alpha_logit <- log(as.numeric(input$AlphaVal)/(1-as.numeric(input$AlphaVal)))
    x_signif <- logit_x < alpha_logit
    y_signif <- logit_y < alpha_logit
    matrix <- cbind(x_signif, y_signif)
    color_var <-
      apply(matrix, 1, function(point)
        if (point[[1]]) {
          if (point[[2]]) {"Significant for both contrasts"}
          else {"Significant for contrast 1"}
        } 
        else{
          if (point[[2]]) {"Significant for contrast 2"}
          else{"Not significant"}
          }
      )
    
    #labels
    names <- row.names(padj)
    ###
    points_to_label <- SelectTaxoPlotCompDebounce()
    #points_to_label <- input$selectTaxoPlotComp
    labels <-
      sapply(names, function(name)
        if (is.element(name, points_to_label)) {
          points_to_label[match(name, points_to_label)]
        } else{
          ""
        })
    
    if(!export){
      plot <- scatterD3(x = logit_x, 
                        y = logit_y,
                        xlab = paste("logit p value", input$Contrast1),
                        ylab = paste("logit p value", input$Contrast2),
                        hover_opacity = 1,
                        tooltip_text = row.names(padj),
                        lines = if (input$showSignifThresholdsLogitPlot) {
                          if(input$showDiagonal){
                            data.frame(
                              slope = c(Inf, 0, 1),
                              intercept = c(alpha_logit,alpha_logit,0),
                              stroke = "black",
                              stroke_width = input$linesWidthLogitPlot,
                              stroke_dasharray = 5
                            )
                          }
                          else{
                          data.frame(
                            slope = c(Inf, 0),
                            intercept = c(alpha_logit,alpha_logit),
                            stroke = "black",
                            stroke_width = input$linesWidthLogitPlot,
                            stroke_dasharray = 5
                          )}
                        }
                        else{if (input$showDiagonal){
                          data.frame(
                            slope = c(1),
                            intercept = c(0),
                            stroke = "black",
                            stroke_width = input$linesWidthLogitPlot,
                            stroke_dasharray = 5
                          )}
                        else{NULL}},
                        fixed = input$fixed11,
                        col_var = color_var,
                        col_lab = if(input$showLegendLogit){"Legend"}else{NA}
                        ,
                        colors = c(
                          "Not significant" = input$colour01,
                          "Significant for contrast 1" = input$colour02,
                          "Significant for contrast 2" = input$colour03,
                          "Significant for both contrasts" = input$colour04
                        ),
                        lab = labels,
                        point_opacity = input$pointOpacityLogit,
                        point_size = input$pointSizeLogit,
                        axes_font_size = paste(input$axisFontSizeLogit, "0%", sep = ""),
                        labels_size = input$labelsSizeLogit,
                        legend_font_size = paste(input$legendFontSizeLogit, "0%", sep = ""), 
                        legend_width = input$legendWidthLogit,
                        # dom_id_reset_zoom = "scatterD3-reset-zoomLogit",
                        # dom_id_svg_export = "scatterD3-svg-exportLogit",
                        # menu = FALSE,
                        # disable_wheel = TRUE
                        )
      }
    
    if(export){
      data <- data.frame(logit_x, logit_y, color_var)
      plot <- ggplot(data, aes(x = logit_x, y = logit_y)) + theme_minimal()
      plot <- plot + geom_point(aes(color = color_var), size = (input$pointSizeLogit %/% 30), alpha = input$pointOpacityLogit)
      plot <- plot + scale_x_continuous(paste("logit p value", input$Contrast1))
      plot <- plot + scale_y_continuous(paste("logit p value", input$Contrast2))
      plot <- plot + theme(axis.title = element_text(size = input$axisFontSizeLogit), 
                           axis.text = element_text(size = input$axisFontSizeLogit),
                           legend.title = element_text(size = input$legendFontSizeLogit),
                           legend.text = element_text(size = input$legendFontSizeLogit))
      plot <- plot + geom_text(aes(label = labels, color = color_var), size = input$labelsSizeLogit %/% 3, hjust = 1, vjust = 1, show.legend = FALSE)
      plot <- plot + scale_color_manual(name = "Legend title", values = c("Not significant" = input$colour01, "Significant for contrast 1" = input$colour02,
                                                                          "Significant for contrast 2" = input$colour03, "Significant for both contrasts" = input$colour04))
      if (input$fixed11) {
        plot <- plot + coord_fixed()
      }
      if (input$showSignifThresholdsLogitPlot) {
        plot <- plot + geom_hline(yintercept = alpha_logit, linetype = "dashed", size = input$linesWidthLogitPlot/2)
        plot <- plot + geom_vline(xintercept = alpha_logit, linetype = "dashed", size = input$linesWidthLogitPlot/2)
      }
      if(input$showDiagonal){
        plot <- plot + geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = input$linesWidthLogitPlot/2)
      }
  }
  return(plot)
}}


##############################
##      VENN DIAGRAM
##############################
Plot_Visu_Venn <- function(input,BaseContrast,resDiff, ContrastListVennDebounce, export=FALSE){
  res = NULL
  contrasts_without_diff = NULL
  #SelContrast = input$ContrastList_table_FCVenn
  SelContrast = ContrastListVennDebounce()
  if(length(SelContrast)>=2 & length(SelContrast)<=4){
    gotData = GetData_venn(input,SelContrast,BaseContrast,resDiff)
    
    data2 = gotData$df.tot
    if(length(SelContrast) - length(colnames(data2[sapply(data2,function(x) all(is.na(x)))])) < 2){contrasts_without_diff = colnames(data2[sapply(data2,function(x) all(is.na(x)))])
    }
    else{
  data = gotData$res
  res = venn_tooltip(d3vennR(data=data))
  }}
  return(list(res = res, contrasts_without_diff=contrasts_without_diff))
}

##############################
##      UPSET
##############################
Plot_UpSet <- function(input,BaseContrast, resDiff, ContrastListDebounce, export=FALSE){
  plot = NULL
  df = NULL
  contrasts_without_diff = NULL
  #SelContrast = input$ContrastList_table_FC
  SelContrast = ContrastListDebounce()
  if(length(SelContrast)>=2){
    gotData = GetData_venn(input,SelContrast,BaseContrast,resDiff)
    data = gotData$df.tot
    
    if(length(SelContrast) - length(colnames(data[sapply(data,function(x) all(is.na(x)))])) < 2){contrasts_without_diff = colnames(data[sapply(data,function(x) all(is.na(x)))])
    }
    else{
    listInput <- list()
    n <- ncol(data)
    for(i in 1:n){
      new <- list(na.omit(data[,i]))
      listInput <- c(listInput, new)
    }
    names(listInput) <- colnames(data)
    
    # for plot
    plot <- upset(fromList(listInput), nsets = n, nintersects = NA, set_size.show = input$showNumbers, 
                  show.numbers = if(input$showNumbers){"yes"}else{"no"}, sets = SelContrast, keep.order = !(input$orderBySize),
                  point.size = input$pointSizeUpSet, line.size = 1.5, order.by = input$orderByUpset, 
                  mainbar.y.label = "Number of differential features in common", sets.x.label = "Number of differential features", 
                  text.scale = c(2, 1.6, 1.6, 1.6, 2, 1.5), decreasing = input$decreasingUpset)
    
    # for table
    combos <- Reduce(c,lapply(rev(1:length(listInput)), function(x) combn(rev(colnames(data)),x,simplify=FALSE) ))
    
    listSignif <- Reduce(function(a,b) {unique(c(as.character(a),as.character(b)))}, listInput)
    n_col = length(combos)
    n_row = length(listSignif)
    df <- matrix(rep(listSignif, times = n_col), ncol = n_col)
    lst <- list()
    for(j in 1:n_col){for(i in 1:n_row){group_contrasts <- unlist(combos[j])
                                        if (is.element(df[i,j], lst)){df[i,j] <- NA}
                                        else{if (!is.element(df[i,j], Reduce(intersect, listInput[group_contrasts]))){df[i,j] <- NA}
                                             else{lst <- c(lst,df[i,j])}
                                            }
                                        }
                      }
    colnames(df) <- lapply(combos, function(lstCont) paste(lstCont, collapse = " & "))

    df = as.data.frame(apply(df,2,Go_data_top))
    maxRow = max(apply(df,2,FUN=function(x) length(which(!is.na(x)))))
    df = df[1:max(maxRow,1),]
    df = df[,which(apply(!is.na(df),2,any))]
    }}
  if(export) res = list(plot=plot,table=df, contrasts_without_diff=contrasts_without_diff)
  else res = list(table=df, contrasts_without_diff=contrasts_without_diff)
  return(res)
}


##############################
##   Contrasts comparison
##############################
Plot_MultipleVenn <- function(input,BaseContrast, resDiff, ContrastListDebounce){
  plot = NULL
  contrasts_without_diff = NULL
  SelContrast = ContrastListDebounce()
  if(length(SelContrast)>=2){
  gotData = GetData_venn(input,SelContrast,BaseContrast,resDiff)
  
  data2 = gotData$df.tot
  if(length(SelContrast) - length(colnames(data2[sapply(data2,function(x) all(is.na(x)))])) < 2){contrasts_without_diff = colnames(data2[sapply(data2,function(x) all(is.na(x)))])
  }
  else{
  
  data = gotData$res_multiple_venn

  plot <- ggplot(data,aes(x=x,y=y)) + geom_point() + theme_bw() + 
          geom_label_repel(aes(label=name), size = input$labelSizemultipleVenn) + xlim(c(0,1)) + ylim(c(0,1)) + 
          xlab(bquote(Contrast1 *intersect(Contrast2) ~ "/ Contrast1")) + 
          ylab(bquote(Contrast1 *intersect(Contrast2) ~ "/ Contrast2"))
  plot <- plot + theme(axis.title = element_text(size = input$FontSizeMultipleVenn), 
                 axis.text = element_text(size = input$FontSizeMultipleVenn - 2))
  }}
  return(list(plot = plot, contrasts_without_diff=contrasts_without_diff))
}


##############################################################
##
##          Useful functions
##
##############################################################

## Get the non NA data at the top of the dataframe
Go_data_top <- function(vect)
{
  n = length(vect)
  tmp = rep(NA,n)
  ind = which(!is.na(vect))
  n1 = length(ind)
  if(n1>0)    tmp[1:n1] =  vect[ind]
  return(tmp)
}


## Get the data to plot the heatmap with the FC
## and for the logit plot
Get_log2FC_padj <-function(input,BaseContrast,resDiff, info = NULL)
{
  log2FC = NULL
  padj = NULL
  dds = resDiff$dds
  counts = resDiff$counts
  target = resDiff$target
  SelContrast = colnames(BaseContrast)
  nbCont = length(SelContrast)
  result = list()
  alpha = as.numeric(input$AlphaVal)
  cooksCutoff = ifelse(input$CooksCutOff!='Auto',ifelse(input$CooksCutOff!=Inf,input$CutOffVal,Inf),TRUE)
  
  if(nbCont>=1)
  {
    for(i in 1:nbCont)
    { 
      cont = as.character(SelContrast[i])
      result[[cont]] <- results(dds,contrast=BaseContrast[,cont],pAdjustMethod=input$AdjMeth,
                                cooksCutoff=cooksCutoff,
                                independentFiltering=input$IndFiltering,alpha=alpha)
    }
    log2FC = as.matrix(round(result[[SelContrast[1]]][, "log2FoldChange"], 3))
    padj = as.matrix((result[[SelContrast[1]]][, "padj"]))
    if(nbCont>1)
    {
      for(i in 2:nbCont)
      {
        log2FC = cbind(log2FC,round(result[[SelContrast[i]]][, "log2FoldChange"], 3))
        padj = cbind(padj,(result[[SelContrast[i]]][, "padj"]))
      }
    }  
    colnames(log2FC) = names(result)
    colnames(padj) = names(result)
    
    rownames(log2FC) = rownames(result[[SelContrast[1]]])
    rownames(padj) = rownames(result[[SelContrast[1]]])
  }
  return(list(log2FC=as.data.frame(log2FC),padj=padj))
}


## Add tooltips on venn digramm
venn_tooltip <- function(venn){
  venn$x$tasks[length(venn$x$tasks)+1] <- list(
    htmlwidgets::JS('
                    function(){
                    var div = d3.select(this);
                    
                    // add a tooltip
                    var tooltip = d3.select("body").append("div")
                    .attr("class", "venntooltip")
                    .style("position", "absolute")
                    .style("text-align", "center")
                    .style("width", 128)
                    .style("height", 16)
                    .style("background", "#333")
                    .style("color","#ddd")
                    .style("padding","2px")
                    .style("border","0px")
                    .style("border-radius","8px")
                    .style("opacity",0);
                    
                    div.selectAll("path")
                    .style("stroke-opacity", 0)
                    .style("stroke", "#fff")
                    .style("stroke-width", 0)
                    
                    // add listeners to all the groups to display tooltip on mousover
                    div.selectAll("g")
                    .on("mouseover", function(d, i) {
                    
                    // sort all the areas relative to the current item
                    venn.sortAreas(div, d);
                    
                    // Display a tooltip with the current size
                    tooltip.transition().duration(400).style("opacity", .9);
                    tooltip.text(d.size);
                    
                    // highlight the current path
                    var selection = d3.select(this).transition("tooltip").duration(400);
                    selection.select("path")
                    .style("stroke-width", 3)
                    .style("fill-opacity", d.sets.length == 1 ? .4 : .1)
                    .style("stroke-opacity", 1);
                    })
                    
                    .on("mousemove", function() {
                    tooltip.style("left", (d3.event.pageX) + "px")
                    .style("top", (d3.event.pageY - 28) + "px");
                    })
                    
                    .on("mouseout", function(d, i) {
                    tooltip.transition().duration(50).style("opacity", 0);
                    var selection = d3.select(this).transition("tooltip").duration(400);
                    selection.select("path")
                    .style("stroke-width", 0)
                    .style("fill-opacity", d.sets.length == 1 ? .25 : .0)
                    .style("stroke-opacity", 0);
                    });
                    }
                    ')
    )
  return(venn)
  }


## Transform the data for the venn diagram and the "contrasts comparison"
GetData_venn <-function(input,SelContrast,BaseContrast,resDiff)
{
  res = list()
  res_multiple_venn = NULL
  df.tot = NULL
  VarInt = input$VarInt
  dds = resDiff$dds
  counts = resDiff$counts
  target = resDiff$target
  nbCont = length(SelContrast)
  result = list()
  alpha = as.numeric(input$AlphaVal)
  cooksCutoff = ifelse(input$CooksCutOff!='Auto',ifelse(input$CooksCutOff!=Inf,input$CutOffVal,Inf),TRUE)
  
  if(nbCont>=2)
  {
    for(i in 1:nbCont)
    { 
      cont = as.character(SelContrast[i])
      result[[cont]] <- results(dds,contrast=BaseContrast[,cont],pAdjustMethod=input$AdjMeth,
                                cooksCutoff=cooksCutoff,
                                independentFiltering=input$IndFiltering,alpha=alpha)
    }
    padj = round(result[[SelContrast[1]]][, "padj"], 3)
    
    # save(result,padj,SelContrast,file = "test1.RData")
    df = as.matrix(rownames(result[[1]]))
    if(length(which(padj>alpha))>0) df[which(padj>alpha),]=NA
    if(any(is.na(padj))) df[which(is.na(padj)),]=NA

    if(nbCont>1)
    {
      for(i in 2:nbCont)
      {
        padj = round(result[[SelContrast[i]]][, "padj"], 3)
        df.tmp = as.matrix(rownames(result[[i]]))
        if(length(which(padj>alpha))>0) df.tmp[which(padj>alpha),]=NA 
        if(any(is.na(padj))) df.tmp[which(is.na(padj)),]=NA  
        df = cbind(df,df.tmp)
      }
      colnames(df) = SelContrast
      df = as.data.frame(df)
    }
    ## Keep the entire dataframe
    df.tot = as.data.frame(apply(df,2,Go_data_top))
    maxRow = max(apply(as.data.frame(apply(df,2,Go_data_top)),2,FUN=function(x) length(which(!is.na(x)))))
    df.tot = df.tot[1:max(maxRow,1),]
    ## Remove col with only NA
    df = df[,which(apply(!is.na(df),2,any))]
    
    ncont = ncol(as.data.frame(df))
    names.df = names(df)
    rownames_multiple_venn = c()
    cmp=1
    if(ncont>1 && !is.null(ncont))
    {
      for(i in 1:(ncont))
      {
        for(j in i:ncont)
          # For Venn diagram (until intersection of 4 contrasts)
        { if(i==j) {res[[cmp]] = list(sets=list(names.df[i]),size= length(which(!is.na(intersect(df[,i],df[,i])))))
                    cmp=cmp+1}
          if(i!=j) {res[[cmp]] = list(sets=list(names.df[i],names.df[j]),size= length(which(!is.na(intersect(df[,i],df[,j])))))
                    cmp=cmp+1
                    for(k in j:ncont){
                      if(k!=j){
                        res[[cmp]] = list(sets=list(names.df[i],names.df[j],names.df[k]), size = length(which(!is.na(intersect(intersect(df[,i],df[,j]),df[,k])))))
                        cmp=cmp+1
                        for(l in k:ncont){
                          if(l!=k){
                            res[[cmp]] = list(sets=list(names.df[i],names.df[j],names.df[k],names.df[l]), size = length(which(!is.na(intersect(intersect(intersect(df[,i],df[,j]),df[,k]),df[,l])))))
                            cmp=cmp+1}
                        }
                        }
                    }
          }
          # For 'Contrasts comparison'
          if(i!=j) {
            if (is.null(res_multiple_venn))
              {res_multiple_venn <- data.frame(name = c(rownames_multiple_venn, paste(names.df[i],names.df[j], sep = " vs ")), 
                                                   x = c(length(which(!is.na(intersect(df[,i],df[,j]))))/length(which(!is.na(intersect(df[,i],df[,i]))))), 
                                                   y = c(length(which(!is.na(intersect(df[,i],df[,j]))))/length(which(!is.na(intersect(df[,j],df[,j]))))),
                                                   stringsAsFactors=FALSE)}
            else{new_row <- data.frame(name = c(rownames_multiple_venn, paste(names.df[i],names.df[j], sep = " vs ")),
                         x = c(length(which(!is.na(intersect(df[,i],df[,j]))))/length(which(!is.na(intersect(df[,i],df[,i]))))),
                         y = c(length(which(!is.na(intersect(df[,i],df[,j]))))/length(which(!is.na(intersect(df[,j],df[,j]))))))
                res_multiple_venn <- rbind(res_multiple_venn, new_row)}
          }
        }
      }
    }
  }
  return(list(res=res,df.tot=df.tot,res_multiple_venn=res_multiple_venn))
}

