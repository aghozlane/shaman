#@ This file contains all the functions for the 
#@ comparison plots of SHAMAN


##############################
##       HEATMAP
##############################
Plot_Visu_Heatmap_FC <- function(input,BaseContrast,resDiff,export=FALSE){
  
  res = NULL
  SelContrast = input$ContrastList_table_FC
  selcontrast_matrix = as.matrix(BaseContrast[,SelContrast])
  colnames(selcontrast_matrix) = SelContrast
  log2FC = Get_log2FC_padj(input,selcontrast_matrix,resDiff, info = NULL)$log2FC
  
  
  if(!is.null(log2FC) && length(SelContrast)>=2)
  { 
    cont = which(colnames(log2FC)%in%SelContrast)
    log2FC = log2FC[,SelContrast] 
    ind_taxo = input$selectTaxoPlotComp
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
##      VENN DIAGRAM
##############################
Plot_Visu_Venn <- function(input,BaseContrast,resDiff,export=FALSE){
  
  res = NULL
  SelContrast = input$ContrastList_table_FC
  data = GetData_venn(input,SelContrast,BaseContrast,resDiff)$res
  res = venn_tooltip(d3vennR(data=data))
  
  return(res)
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
    padj = as.matrix(round(result[[SelContrast[1]]][, "padj"], 3))
    if(nbCont>1)
    {
      for(i in 2:nbCont)
      {
        log2FC = cbind(log2FC,round(result[[SelContrast[i]]][, "log2FoldChange"], 3))
        padj = cbind(padj,round(result[[SelContrast[i]]][, "padj"], 7))
      }
      colnames(log2FC) = names(result)
      colnames(padj) = names(result)
    }
    
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


## Transform the data for the venn diagram
GetData_venn <-function(input,SelContrast,BaseContrast,resDiff)
{
  res = list()
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
    cmp=1
    if(ncont>1 && !is.null(ncont))
    {
      for(i in 1:(ncont))
      {
        for(j in i:ncont)
        {
          if(i!=j) res[[cmp]] = list(sets=list(names.df[i],names.df[j]),size= length(which(!is.na(intersect(df[,i],df[,j])))))
          if(i==j) res[[cmp]] = list(sets=list(names.df[i]),size= length(which(!is.na(intersect(df[,i],df[,i])))))
          cmp=cmp+1
        }
      }
    }
    
  }
  return(list(res=res,df.tot=df.tot))
}

