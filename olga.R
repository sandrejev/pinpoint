library(ddply)
library(reshape2)
library(ggplot2)
library(gridExtra)

benchmar = function()
{
    files = list.files("benchmark2/BY/", pattern=".*\\.tab$", full.names=T)
    
    data = data.frame()
    for(file in files)
    {
        d = read.delim(file, stringsAsFactors=FALSE)
        d$day = as.numeric(gsub(".*(\\d+)day.*", "\\1", basename(file)))
        d$format = as.numeric(gsub("(\\d+)_.*", "\\1", basename(file)))
        d.wide = dcast(d, row ~ col, value.var="s.area") > 100
        d.wide2 = matrix(0, nrow=nrow(d.wide)+2, ncol=ncol(d.wide)+2)
        d.wide2[2:c(nrow(d.wide)+1), 2:ncol(d.wide)] = d.wide[,-1]
        d.wide = distmap(d.wide2, metric="manhattan")[2:c(nrow(d.wide)+1), 2:ncol(d.wide)]
        d = merge(d, melt(d.wide,value.name="edge.dist"), by.x=c("row", "col"), by.y=c("Var1", "Var2"))    
        d$row_96 = ceiling(d$row / (format.dim(d$format[1], 1) / 8))
        d$col_96 = ceiling(d$col / (format.dim(d$format[1], 2) / 12))
        d$edge.distn = d$edge.dist / (format.dim(d$format[1], 1) / 8)
        
        data = rbind(data, d)
    }
    
    
    data.sum = ddply(data, .(row_96, col_96, day), summarize,
                     saturation.b.mean2=mean(saturation.b.mean),
                     value.b.mean2=mean(value.b.mean),
                     hue.b.mean2=mean(value.b.mean),
                     red.b.mean2=mean(red.b.mean),
                     green.b.mean2=mean(green.b.mean),
                     blue.b.mean2=mean(blue.b.mean)
    )
    
    pdf("edge.effect5.pdf", paper="a4r", height=842, width=595)
    
    grid.arrange(
        ggplot(subset(data, day %in% 1:2 & edge.dist>0)) + geom_density(aes(x=saturation.b.mean, fill=factor(edge.dist == 1)), alpha=0.3) +
            facet_grid(format ~ day, scales="free") + theme_classic(base_size = 20) + scale_fill_discrete(name="Distance == 1"),
        ggplot(subset(data, day %in% 3:4 & edge.dist>0)) + geom_density(aes(x=saturation.b.mean, fill=factor(edge.dist == 1)), alpha=0.3) +
            facet_grid(format ~ day, scales="free") + theme_classic(base_size = 20) + scale_fill_discrete(name="Distance == 1")
    )
    
    grid.arrange(
        ggplot(subset(data, day %in% 1:2 & edge.dist>1)) + geom_density(aes(x=saturation.b.mean, fill=factor(edge.dist == 2)), alpha=0.3) +
            facet_grid(format ~ day, scales="free") + theme_classic(base_size = 20) + scale_fill_discrete(name="Distance == 2"),
        ggplot(subset(data, day %in% 3:4 & edge.dist>1)) + geom_density(aes(x=saturation.b.mean, fill=factor(edge.dist == 2)), alpha=0.3) +
            facet_grid(format ~ day, scales="free") + theme_classic(base_size = 20) + scale_fill_discrete(name="Distance == 2")
    )
    
    dev.off()
    
    
    pdf("edge.effect.avg.pdf", paper="a4r", height=842, width=595)
    
    data2 = ddply(data, .(day, format, row_96, col_96, edge.dist), summarize, saturation.b.mean=mean(saturation.b.mean))    
    grid.arrange(
        ggplot(subset(data2, day %in% 1:2 & edge.dist>0)) + geom_density(aes(x=saturation.b.mean, fill=factor(edge.dist == 1)), alpha=0.3) +
            facet_grid(format ~ day, scales="free") + theme_classic(base_size = 20) + scale_fill_discrete(name="Distance == 1"),
        ggplot(subset(data2, day %in% 3:4 & edge.dist>0)) + geom_density(aes(x=saturation.b.mean, fill=factor(edge.dist == 1)), alpha=0.3) +
            facet_grid(format ~ day, scales="free") + theme_classic(base_size = 20) + scale_fill_discrete(name="Distance == 1")
    )
    
    grid.arrange(
        ggplot(subset(data2, day %in% 1:2 & edge.dist>1)) + geom_density(aes(x=saturation.b.mean, fill=factor(edge.dist == 2)), alpha=0.3) +
            facet_grid(format ~ day, scales="free") + theme_classic(base_size = 20) + scale_fill_discrete(name="Distance == 2"),
        ggplot(subset(data2, day %in% 3:4 & edge.dist>1)) + geom_density(aes(x=saturation.b.mean, fill=factor(edge.dist == 2)), alpha=0.3) +
            facet_grid(format ~ day, scales="free") + theme_classic(base_size = 20) + scale_fill_discrete(name="Distance == 2")
    )
    
    dev.off()
    
    pdf("correlation.pdf", paper="a4r", height=842, width=595)
    
    data2.384_96 = merge(subset(data2, format==384), subset(data2, format==96), by=c("row_96", "col_96", "day"), suffixes=c("_384", "_96"))
    ggplot(data2.384_96, aes(saturation.b.mean_384, saturation.b.mean_96, group=factor(day))) +
        geom_text(aes(color=edge.dist_384, label=day)) + 
        geom_smooth(method="lm") + theme_classic()
    
    data2.384_96.cor = ddply(data2.384_96, .(day), summarize, correlation=cor(saturation.b.mean_384, saturation.b.mean_96, method="spearman"), x=min(saturation.b.mean_384), y=max(saturation.b.mean_96))
    ggplot(data2.384_96, aes(saturation.b.mean_384, saturation.b.mean_96, group=factor(day))) +
        geom_text(aes(color=edge.dist_384, label=day)) + 
        geom_smooth(method="lm") + theme_classic() + 
        geom_text(aes(x, y, label=paste("S =", gsub("^.", "", round(correlation, 2)))), data=data2.384_96.cor) +
        facet_wrap( ~ day, scales="free")
    
    dev.off()
}