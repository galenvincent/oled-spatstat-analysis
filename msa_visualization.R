# MSA Visualization of Results
library(rapt)
library(dplyr)
library(plotly)
library(RColorBrewer)

# Upload APT data
pos <- readPOS('Z:/Galen/Stephan\ APT\ Analysis/R34_06365-v01-CentralCubeExclPole.pos')
rng.simple <- readRRNG('Z:/Galen/Stephan\ APT\ Analysis/6365_Alalloy_SimpleRngs.RRNG')

ms <- createSpec(pos, res = 0.05)
ms.log <- transformIntensity(ms, method = 'log10')

#prettyPlot(ms.log, rng = rng.full, xlim = c(0, 75), main= 'Full RRNG')
#prettyPlot(ms.log, rng = rng.simple, xlim = c(0, 75), main = 'Simple RRNG')

pos.r <- rngPOS(pos, rng.simple)
cnts <- rngCount(pos.r, rng.simple)
cnts.totals <- group_by(cnts, name) %>%
  summarize(counts = sum(counts), fraction = sum(fraction))
cnts.totals

frac <- sum(cnts.totals$fraction[2:3])

win.r <- box3(c(min(pos.r$x), max(pos.r$x)), 
              c(min(pos.r$y), max(pos.r$y)), 
              c(min(pos.r$z), max(pos.r$z)))

# Upload MSA results
load('Z:/Galen/MSA/msa.sweep.andrew.RData')

# Set msa parameters used
nmins <- c(5, 10, 15, 20)
dmaxs <- seq(1, 2, 0.1)
params <- expand.grid(nmins, dmaxs)

# Create pp3 of different elements
pp3.full <- createSpat(pos.r, win = win.r)
marks(pp3.full) <- pos.r$mark

pp3.df <- data.frame(x = pp3.full$data$x, y = pp3.full$data$y, z = pp3.full$data$z, mark = pp3.full$data$marks)
pp3.df.zn <- filter(pp3.df, mark == 'Zn1')
pp3.df.mg <- filter(pp3.df, mark == 'Mg1')
pp3.df.al <- filter(pp3.df, mark == 'Al1')
pp3.df.aggregates <- filter(pp3.df, mark == 'Mg1' | mark == 'Zn1' )
pp3.aggregates <- subset(pp3.full, marks == 'Mg1' | marks == 'Zn1' )


dmax.to.plot <- c(1.1, 1.4, 1.7, 2.0)
nmin.to.plot <- c(5, 10, 15, 20)
inds.to.do <- which(round(params[,2],1) %in% dmax.to.plot)
n <- length(inds.to.do)

data.list <- list()
cnt <- 1
for(i in inds.to.do){
  if(is.na(msa.sweep[[i]])){
    data.list[[cnt]] <- data.frame(x = c(0), y = c(0), z = c(0), mark = '', col.mark = '')
  }else{
    data.list[[cnt]] <- pp3.df[unlist(msa.sweep[[i]]$A),]
    nums <- sample(1:length(msa.sweep[[i]]$A), length(msa.sweep[[i]]$A))
    data.list[[cnt]]$col.mark <- rep(nums, sapply(msa.sweep[[i]]$A, length))
  }
  cnt <- cnt + 1
}

dropdowns <- list()
dropdowns[[1]] <- list(label = 'Zn+Mg',
                       method = 'update',
                       args = list(list(visible = rep(FALSE, n+1))))
dropdowns[[1]]$args[[1]]$visible[1] <- TRUE
for(i in 1:n){
  dropdowns[[i+1]] <- list(label = paste('dmax: ', toString(params[inds.to.do[i],2]), ', Nmin: ', toString(params[inds.to.do[i],1]), sep = ''),
                         method = 'update',
                         args = list(list(visible = rep(FALSE, n+1))))
  #Show background
  dropdowns[[i+1]]$args[[1]]$visible[1] <- TRUE
  
  #Highlight clusters
  dropdowns[[i+1]]$args[[1]]$visible[i+1] <- TRUE
}

p <- plot_ly()

p <- add_trace(p, type = 'scatter3d', mode = 'markers', data = pp3.df.aggregates,
               x = ~x, y = ~y, z = ~z,
               marker = list(opacity = 0.1, 
                             size = 3,
                             color = 'gray'))
for(i in 1:n){
  p <- add_trace(p, type = 'scatter3d', mode = 'markers', data = data.list[[i]],
                   x = ~x, y = ~y, z = ~z,
                   marker = list(opacity = 0.3, 
                                 size = 3,
                                 color = ~col.mark, 
                                 colorscale = "Jet"))
}

p <- layout(p, 
            title = 'Interactive MSA Plot',
            scene = list(xaxis = list(range = c(-20.2, 20.2)),
                         yaxis = list(range = c(-20.2, 20.2)),
                         zaxis = list(range = c(-20.2, 20.2)),
                         aspectratio = list(x = 1, y = 1, z = 1)),
            showlegend = FALSE,
            updatemenus = list(
              list(
                #xanchor = 'right',
                #yanchor = 'top',
                #x = 1, 
                #y = 1,
                active = 0,
                showactive = TRUE,
                buttons = dropdowns
              )))

p


library(ks)
imap3d <- kde(coords(pp3.aggregates), gridsize = 100, xmin = c(-20.2, -20.2, -20.2), xmax = c(20.2, 20.2, 20.2))
plot(imap3d, cont = c(5), alphavec = c(0.25))
plot3d.pp3(pp3.aggregates)


