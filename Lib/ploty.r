
library(dplyr)
library(magrittr)
library(cowplot)

library(Seurat)

SingleExIPlot2 <- function(
  data,
  idents,
  split = NULL,
  type = 'violin',
  sort = FALSE,
  y.max = NULL,
  adjust = 1,
  pt.size = 0,
  cols = NULL,
  log = FALSE
) {
  set.seed(seed = 42)
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
  }
  feature <- colnames(x = data)
  data$ident <- idents
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(
        x = tapply(
          X = data[, feature],
          INDEX = data$ident,
          FUN = mean
        ),
        decreasing = grepl(pattern = paste0('^', tolower(x = sort)), x = 'decreasing')
      )))
    )
  }
  if (log) {
    noise <- rnorm(n = length(x = data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- rnorm(n = length(x = data[, feature])) / 100000
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else{
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- ifelse(test = log, yes = 'Log Expression Level', no = 'Expression Level')
  y.max <- max(data[, feature])
  if (is.null(x = split) || type != 'violin') {
    vln.geom <- geom_violin
    fill <- 'ident'
  } else {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- 'split'
  }
  switch(
    EXPR = type,
    'violin' = {
      x <- 'ident'
      y <- paste0("`", feature, "`")
      xlab <- 'Identity'
      ylab <- axis.label
      geom <- list(
        vln.geom(scale = 'width', adjust = adjust, trim = TRUE),
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      jitter <- geom_jitter(height = 0, size = pt.size)
      log.scale <- scale_y_log10()
      axis.scale <- ylim
    },
    'ridge' = {
      x <- paste0("`", feature, "`")
      y <- 'ident'
      xlab <- axis.label
      ylab <- 'Identity'
      geom <- list(
        geom_density_ridges(scale = 4),
        theme_ridges(),
        scale_y_discrete(expand = c(0.01, 0)),
        scale_x_continuous(expand = c(0, 0))
      )
      jitter <- geom_jitter(width = 0, size = pt.size)
      log.scale <- scale_x_log10()
      axis.scale <- function(...) {
        invisible(x = NULL)
      }
    },
    stop("Unknown plot type: ", type)
  )
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]
  ) +
    labs(x = xlab, y = ylab, title = feature, fill = NULL) +
    theme_cowplot()
  plot <- do.call(what = '+', args = list(plot, geom))
  plot <- plot + if (log) {
    log.scale
  } else {
    axis.scale(min(data[, feature]), y.max)
  }
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      } else {
        unlist(x = lapply(
          X = idents,
          FUN = function(pattern, x) {
            x.mod <- gsub(
              pattern = paste0(pattern, '.'),
              replacement = paste0(pattern, ': '),
              x = x,
              fixed = TRUE
            )
            x.keep <- grep(pattern = ': ', x = x.mod, fixed = TRUE)
            x.return <- x.mod[x.keep]
            names(x = x.return) <- x[x.keep]
            return(x.return)
          },
          x = unique(x = as.vector(x = data$split))
        ))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    } else {
      labels <- unique(x = as.vector(x = data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  
  return(plot)
}






myvenn=function(res,fout){
  library(VennDiagram)
  library(grid)

  dnames=names(res)
  nc=length(dnames)
  cols=c('red','blue','yellow','green','grey')
  aaa=venn.diagram(res,fout,category.names =dnames,filename = NULL,fill=cols[1:nc])
  jpeg(fout);
  grid.draw(aaa);
  dev.off();

}






figExp=function(emb,colr,size=0.2,alpha=0.4){
  
  # colr=p2$counts[cname,'CCL20']
  #size=0.01
  #alpha=0.5
  
  
  embedding=emb
  
  
  plot.df <- tibble::rownames_to_column(as.data.frame(embedding), 
                                        "CellName")
  colnames(plot.df)[2:3] <- c("x", "y")
  geom_point_w <- ggplot2::geom_point
  
  
  colors=colr
  
  color.range <- range(colors)
  
  gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x = x, y = y)) 
  
  gg <- gg + geom_point_w(col='#d8d0d0', 
                          alpha = alpha, size = size)
  
  #gg <- gg + ggplot2::scale_color_gradient2(low = "#0000ff",  mid = "#d8d0d0", high = "#ff0000", limits = color.range)
  
  rownames(plot.df)=plot.df[,1]
  nname=names(colr[colr>0])
  
  plot.df2=plot.df[nname,]
  
  #gg=gg+geom_point(data=plot.df2, aes(x = x, y = y),alpha = 0.4, size = 0.01,aes(col = colors[nname]))
  
  gg=gg+geom_point(data=plot.df2, aes(x = x, y = y,colour = colors[nname]),alpha = alpha, size = size)
  
  gg <- gg + ggplot2::scale_color_gradient2(low = "#0000ff",  mid = "#d8d0d0", high = "#ff0000", limits = color.range)
  
  gg <- gg + ggplot2::theme(legend.position = "none")
  gg <- gg + ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                            axis.text = ggplot2::element_blank())
  gg <- gg + ggplot2::theme(axis.title = ggplot2::element_blank())
  
  
  gg <- gg + theme(   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  
  return(gg)
  
}

# t=p2$counts[cname,gene]
#   a=figExp(emb,t,size=0.01)
#  a=a+ theme(plot.title = element_text(size=12,hjust = 0.5))+ggtitle(gene)
