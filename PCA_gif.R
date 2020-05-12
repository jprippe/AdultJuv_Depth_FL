require(vegan)

# Generate you PCA dataset (you'll need to change the grouping variables)
pp0=capscale(ibsMat~1)
pp0.admix <- data.frame(id=bams, age=substr(bams,4,4), habitat=substr(bams,3,3), admix=cluster.admix, pp0$CA$u)

# 3D PCA in plot.ly -------------------------------------------------------

library(plotly)

fig <- plot_ly(pp0.admix, x = ~MDS1, y = ~MDS2, z = ~MDS3, mode = 'markers', color = ~habitat, symbol = ~age, 
               colors = c('skyblue1', 'palegreen1', 'plum1'), marker = list(size=4)) %>% 
  add_markers() %>% 
  layout(scene = list(xaxis = list(title = 'MDS1'),
                      yaxis = list(title = 'MDS2'),
                      zaxis = list(title = 'MDS3'),
                      camera = list(eye = list(x = 1, y = 1, z = 1))))
p <- plotly_build(fig)
p$x$data[[2]]$marker$symbol <- 'diamond'
p

# Output a rotational image series to create a gif
for(i in seq(0,6.3,by=0.1)){
  outfile <- paste("PCA",round(i,digits=2), sep = "_")
  cam.zoom = 2
  ver.angle = 0
  fig <- plot_ly(pp0.admix, x = ~MDS1, y = ~MDS2, z = ~MDS3, 
                 mode = 'markers', color = ~habitat, symbol = ~age, 
                 colors = c('skyblue1', 'palegreen1', 'plum1'), marker = list(size=4)) %>% 
    add_markers() %>% 
    layout(scene = list(xaxis = list(title = 'MDS1'),
                        yaxis = list(title = 'MDS2'),
                        zaxis = list(title = 'MDS3'),
                        camera = list(eye = list(x = cos(i)*cam.zoom,y = sin(i)*cam.zoom, z=0.2),
                                      center = list(x = 0, y = 0, z = 0))))
  p <- plotly_build(fig)
  p$x$data[[2]]$marker$symbol <- 'diamond'
  p
  
  cat("Now rendering iteration:", i,"\n")
  orca(fig, paste(outfile,"png", sep="."), width = 1200, height = 1050)
}
