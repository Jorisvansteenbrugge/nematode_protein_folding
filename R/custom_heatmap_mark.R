
scale.legend <- function(input,bg.col,col.scale,x.lab,start.zero,lwd.nu,cex.nu,xlab.line,ylab.line,line.col){
                          if(missing(input)){                                   stop("no input defined")}
                          if(missing(bg.col)){                                  bg.col <- "white"}
                          if(missing(col.scale)){                               col.scale <- c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF","#F7F7F7","#E6F5D0","#B8E186","#7FBC41","#4D9221","#276419")}
                          if(missing(x.lab)){                                   x.lab <- ""}
                          if(missing(start.zero)){                              start.zero <- TRUE}
                          if(missing(lwd.nu)){                                  lwd.nu <- 2}
                          if(missing(cex.nu)){                                  cex.nu <- 1}
                          if(missing(xlab.line)){                               xlab.line <- 2.5}
                          if(missing(ylab.line)){                               ylab.line <- 2.5}
                          if(missing(line.col)){                                line.col <- "white"}

                          ###Remove NA

                          ###Make the color scale
                          colors.use <- color.scale(col.scale,input)

                          ###first plotting layer
                          if(start.zero==TRUE){x.range <- c(0,max(input,na.rm=T))}
                          if(start.zero==FALSE){x.range <- c(min(input,na.rm=T),max(input,na.rm=T))}

                          x.val <- axis.output(x.range[1],x.range[2])[[1]]
                              x.val <- x.val[c(1:which(x.val > max(x.range))[1])]
                          y.val <- hist(input,breaks=x.val,plot=FALSE)$counts

                          y.range <- c(0,max(y.val))
                          x.range <- c(min(x.val),max(x.val))

                          plot(x=x.val,y=c(y.val,y.val[length(y.val)]),type="s",xaxs="i",axes=F,xlab="",ylab="",xlim=x.range,ylim=y.range)
                          col.grad <- seq(x.range[1],x.range[2],length.out=length(colors.use)+1)
                          for(i in 2:length(col.grad)){
                              rect(xleft=col.grad[i-1],xright=col.grad[i]+10000,ybottom=-40^6,ytop=40^6,border=NA,col=colors.use[i-1])
                          }
                          points(x=x.val,y=c(y.val,y.val[length(y.val)]),type="s",col=line.col,lwd=lwd.nu,lend=1)

                          ###close it up, maxe axes etc
                          box()
                          fixed.axis <- axis.output(x.range[1],x.range[2])
                          axis(1,fixed.axis[[1]],labels=F,tcl=-0.2)
                          if(start.zero){axis(1,c(0,fixed.axis[[2]]),c(0,fixed.axis[[2]]),mgp=c(3,0.7,0),cex.axis=cex.nu,tcl=-0.5,las=1)}
                          if(start.zero==FALSE){axis(1,fixed.axis[[2]],fixed.axis[[2]],mgp=c(3,0.7,0),cex.axis=cex.nu,tcl=-0.5,las=1)}
                          mtext(side=1,text=x.lab,cex=cex.nu+0.5,line=xlab.line)
                          fixed.axis <- axis.output(y.range[1],y.range[2])
                          axis(2,fixed.axis[[1]],labels=F,tcl=-0.2)
                          axis(2,fixed.axis[[2]],fixed.axis[[2]],mgp=c(3,0.7,0),cex.axis=cex.nu,tcl=-0.5,las=1)
                          mtext(side=2,"Counts",cex=cex.nu+0.5,line=ylab.line)
                         }


axis.output <- function(min.nu,max.nu){

                        output <- NULL; output <- as.list(output)
                        for(i in sort(c(10^c(seq(-2,8,by=1)),2*10^c(seq(-2,8,by=1)),4*10^c(seq(-2,8,by=1)),6*10^c(seq(-2,8,by=1)),8*10^c(seq(-2,8,by=1))),decreasing=T)){
                            ###determine the scale of the differences
                            if(max.nu-min.nu < i*2){
                                small.axis <- seq(floor(min.nu/i)*i,ceiling(max.nu/i)*i*2,by=i/20)
                                large.axis <- seq(floor(min.nu/i)*i,ceiling(max.nu/i)*i*2,by=i/4)[-1]
                            }
                        }
                        output[[1]] <- small.axis
                        output[[2]] <- large.axis
                        return(output)
}
