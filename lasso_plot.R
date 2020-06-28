#自定义颜色
mycol <- rep(c("#223D6C","#D20A13","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)


plotCoef_plus <- function (beta, norm, lambda, df, dev, label = FALSE, legend = FALSE, xvar = c("norm", 
                                                                                                "lambda", "dev"), xlab = iname, ylab = "Coefficients", ...) 
{
  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = "Log Lambda"
    approx.f = 0
  }, dev = {
    index = dev
    iname = "Fraction Deviance Explained"
    approx.f = 1
  })
  dotlist = list(...)
  type = dotlist$type
  
  if (legend){
    #在右侧留出画图例的地方
    par(xpd = T, mar = par()$mar + c(0,0,0,6))
  }
  
  #修改bty，换个更好看的边框，还可以改成，o / n / 7 / l / c / u / ]
  if (is.null(type)) 
    matplot(index, t(beta), lty = 1, lwd = 2,
            xlab = xlab, ylab = ylab, 
            xlim = c(0, xmax), #设置x轴最大值
            col = mycol,#线的颜色
            type = "l", cex.lab=1.2, cex.axis=1,
            bty="n", ...)#不画右边框
  else matplot(index, t(beta), lty = 1, lwd = 2,
               xlab = xlab, ylab = ylab, 
               xlim = c(0, xmax), 
               col = mycol,
               type = "l", cex.lab=1.2, cex.axis=1,
               bty="n", ...)
  atdf = pretty(index)
  prettydf = approx(x = index, y = df, xout = atdf, rule = 2, 
                    method = "constant", f = approx.f)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA)
  
  ## 这边很重要，显示基因名称
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[gene_found, ncol(beta)]
    
    #原函数打印序号，修改为打印基因名
    text(xpos, ypos, gene_found,
         cex = 0.8, #基因名字体大小
         #基因名的颜色跟线一样
         col = mycol,
         #如果你不想要彩色的字，就用下面这行
         #col = "black",
         pos = pos)
  }
  if (legend) {
    #画图例
    legend("topright",
           inset=c(-0.12,0),#图例画到图外面
           legend = rownames(myexpr), #图例文字
           col = mycol, #图例线的颜色，与文字对应
           lwd = 3, #图例中线的粗细
           cex = 1, #图例字体大小
           bty = "n") #不显示图例边框
  }
  par(xpd=FALSE)
}

plot.glmnet_plus <- function (x, xvar = c("norm", "lambda", "dev"), label = FALSE, legend = FALSE,
                              ...) 
{
  xvar = match.arg(xvar)
  plotCoef_plus(x$beta, lambda = x$lambda, df = x$df, dev = x$dev.ratio, 
                label = label, legend = legend, xvar = xvar, ...)
}

