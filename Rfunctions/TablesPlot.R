#@ This file contains all the functions for the
#@ plots in tab "Tables" of SHAMAN


###########################
##        Bar chart
###########################
Bar_Chart_Tables <- function(input, data, export = FALSE) {
  data_table <- data.table::data.table(data$complete)
  data_table <-
    data_table[order(eval(parse(text = input$ColumnOrder)), decreasing = input$Decreasing)]
  
  data_table <-
    data_table[data_table$Id %in% input$selectTaxoPlotBarChart,]
  
  #set colors
  x <- data_table[["log2FoldChange"]]
  y <- data_table[["pvalue_adjusted"]]
  y_log <- -1 * log10(y)
  y_signif <- -log10(as.numeric(input$AlphaVal))
  x_signif <-
    if (!input$setXsignifThreshold) {
      0
    } else{
      input$XsignifThreshold
    }
  x_up <- x > x_signif
  x_down <- x < -x_signif
  y_sig <- y_log > y_signif
  matrix <- cbind(x_up, x_down, y_sig)
  Group <- 
    apply(matrix, 1, function(point)
    if ((anyNA(point)) |
        (!point[[3]])) {
      "Not significant"
    } else{
      if (point[[1]]) {
        "Up"
      } else{
        if (point[[2]]) {
          "Down"
        } else{
          "Not significant"
        }
      }
    })
  colors <-
    apply(matrix, 1, function(point)
      if ((anyNA(point)) |
          (!point[[3]])) {
        input$colour2
      } else{
        if (point[[1]]) {
          input$colour3
        } else{
          if (point[[2]]) {
            input$colour1
          } else{
            input$colour2
          }
        }
      })
  data_table[["color"]] <- colors
  data_table[["Group"]] <- Group
  
  if(!export){
    plot <- amBarplot(
      x = "Id",
      y = "log2FoldChange",
      data = data_table,
      horiz = TRUE,
      ylab = "log 2 Fold Change",
      creditsPosition = "bottom-right"
    )
    plot <- setProperties(plot, fontSize = input$fontSize)
  }
  
  else{
    data <- data_table
    data$Id <- factor(data$Id, levels = rev(data$Id))
    plot <- ggplot(data, aes(x = Id, y = log2FoldChange)) + theme_minimal()
    plot <- plot + geom_bar(stat = "identity", aes(color = Group, fill = Group), width = 0.7)
    plot <- plot + coord_flip()
    plot <- plot + scale_color_manual(values = c("Down" = input$colour1,"Not significant" = input$colour2,"Up" = input$colour3))
    plot <- plot + scale_fill_manual(values = c("Down" = input$colour1,"Not significant" = input$colour2,"Up" = input$colour3))
    plot <- plot + labs(x = "", y = "log 2 Fold Change")
    plot <- plot + theme(axis.title = element_text(size = input$fontSize),
                         axis.text = element_text(size = input$fontSize),
                         legend.title = element_text(size = input$fontSize),
                         legend.text = element_text(size = input$fontSize))
    if(input$removeLabelsBarChart){plot <- plot + theme(axis.ticks.y = element_blank(),
                                                        axis.text.y = element_blank())}
  }
  return(plot)
}

###########################
##        Volcano plot
###########################
Volcano_Plot <- function(input, data, export = FALSE) {

  x <- data$complete[["log2FoldChange"]]
  y <- data$complete[["pvalue_adjusted"]]
  y_log <- -1 * log10(y)
  
  #set axes limits, horizontally centering the chart
  xmin <- -1 * max(-1 * min(x[!is.na(x)]), max(x[!is.na(x)]))
  xmax <- -xmin
  deltaX <- xmax - xmin
  ymin <- min(y_log[!is.na(y_log)])
  ymax <- max(y_log[!is.na(y_log)])
  deltaY <- ymax - ymin
  margin <- 0.1
  xlimits <- c(xmin - margin * deltaX, xmax + margin * deltaX)
  ylimits <- c(ymin - margin * deltaY / 10, ymax + margin * deltaY)
  
  #significance thresholds
  y_signif <- -log10(as.numeric(input$AlphaVal))
  x_signif <-
    if (!input$setXsignifThreshold) {
      0
    } else{
      input$XsignifThreshold
    }
  
  #comparison to thresholds
  x_up <- x > x_signif
  x_down <- x < -x_signif
  y_sig <- y_log > y_signif
  matrix <- cbind(x_up, x_down, y_sig)
  color_var <-
    apply(matrix, 1, function(point)
      if ((anyNA(point)) |
          (!point[[3]])) {
        "Not significant"
      } else{
        if (point[[1]]) {
          "Up"
        } else{
          if (point[[2]]) {
            "Down"
          } else{
            "Not significant"
          }
        }
      })
  
  names <- data$complete[["Id"]]
  points_to_label <- input$selectTaxoLabelVolcano
  labels <-
    sapply(names, function(name)
      if (is.element(name, points_to_label)) {
        points_to_label[match(name, points_to_label)]
      } else{
        ""
      })
  if(!export){
  plot <-
    scatterD3(
      x = x,
      y = y_log,
      xlab = "Log 2 Fold Change",
      ylab = "- Log 10 p value adjusted",
      hover_opacity = 1,
      tooltip_text = names,
      xlim = xlimits,
      ylim = ylimits,
      lines = if (input$showSignifThresholds) {
                  data.frame(
                            slope = c(Inf, Inf, 0),
                            intercept = c(-x_signif, x_signif, y_signif),
                            stroke = "black",
                            stroke_width = input$signifThresholdsWidth,
                            stroke_dasharray = 5
                            )
                  } 
              else{NULL},
      col_var = color_var,
      col_lab = NA, #if(input$showLegend){"Legend title"}else{NA}
      colors = c(
        "Down" = input$colour1,
        "Not significant" = input$colour2,
        "Up" = input$colour3
      ),
      lab = labels,
      point_opacity = input$pointOpacity,
      point_size = input$pointSize,
      axes_font_size = paste(input$axisFontSize, "0%", sep = ""),
      labels_size = input$labelsSize,
      # legend_font_size = paste(input$legendFontSize, "%", sep = ""), 
      # legend_width = input$legendWidth,
      dom_id_reset_zoom = "scatterD3-reset-zoomVolcano",
      dom_id_svg_export = "scatterD3-svg-exportVolcano",
      menu = FALSE,
      height = input$heightVolcanoPlot,
      width=if(input$modifwidthVolcano){input$widthVolcanoPlot},
      disable_wheel = FALSE
    )
  return(plot)
  }
  
  if(export){
    data <- data.frame(x,y_log, color_var)
    plot <- ggplot(data, aes(x = x, y = y_log)) + theme_minimal()
    plot <- plot + geom_point(aes(color = color_var), size = (input$pointSize %/% 30), alpha = input$pointOpacity)
    plot <- plot + scale_x_continuous("Log 2 Fold Change", limits = xlimits)
    plot <- plot + scale_y_continuous("- Log 10 p value adjusted", limits = ylimits)
    plot <- plot + theme(axis.title = element_text(size = input$axisFontSize), 
                         axis.text = element_text(size = input$axisFontSize),
                         legend.title = element_text(size = input$axisFontSize),
                         legend.text = element_text(size = input$axisFontSize))
    plot <- plot + geom_text(aes(label = labels, color = color_var), size = input$labelsSize %/% 3, hjust = 1, vjust = 1, show.legend = FALSE)
    plot <- plot + scale_color_manual(name = "Group", values = c("Down" = input$colour1,"Not significant" = input$colour2,"Up" = input$colour3))
    if(input$showSignifThresholds){
      plot <- plot + geom_hline(yintercept = y_signif, linetype = "dashed", size = input$signifThresholdsWidth/2)
      plot <- plot + geom_vline(xintercept = x_signif, linetype = "dashed", size = input$signifThresholdsWidth/2) + geom_vline(xintercept = -x_signif, linetype = "dashed", size = input$signifThresholdsWidth/2)
    }
    return(plot)
  }
}