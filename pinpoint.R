library(EBImage)
library(reshape2)
library(plyr)
source("transformations.R")
source("visualizations.R")


find.edges = function(img.bw)
{
    #img.bw_median = EBImage::medianFilter(img.bw, 3)
    img.bw_median = img.bw
    
    f1 = matrix(c(-1,-2,-1, 0,0,0, 1,2,1), nrow=3)
    f1r = matrix(c(1,2,1, 0,0,0, -1,-2,-1), nrow=3)
    f2 = matrix(c(-1,0,1, -2,0,2, -1,0,1), nrow=3)
    f2r = matrix(c(1,0,-1, 2,0,-2, 1,0,-1), nrow=3)
    g1 = filter2(img.bw_median+0.5, f1)
    g1r = filter2(img.bw_median+0.5, f1)
    g2 = filter2(img.bw_median+0.5, f2)
    g2r = filter2(img.bw_median+0.5, f2)
    img.bw_sobel = sqrt(((g1+g1r)/2)^2 + ((g2+g2r)/2)^2)
    
    list(img=img.bw_sobel)
}

find.peaks = function(mask, format, r, margin=1, p=F, trace=1)
{
    # The row/collumn format is inverse of x/y format used by EBImage
    margin.i = setdiff(1:2, margin)
    dim_exp = format.dim(format, margin.i)
        
    library(zoo)
    library(quantmod)
    pe = apply(mask>0, margin, mean, na.rm=T)
    pe.rle = rle(pe)
    pe.rle = pe.rle$length[pe.rle$values > .95]
    pe.rle = pe.rle[pe.rle > r*1.5]
    s = ifelse(length(pe.rle) > dim_exp/2, median(pe.rle), r*1.5)
    s = pmax(s, r*1.5)
    
    pe.mean = zoo::rollmean(pe,  round.odd(s), na.pad=T, fill=c(pe[1], NA, last(pe)))
    
    px = quantmod::findPeaks(pe.mean)
    px = px[pe.mean[px] > 0.05]
    px.orig = px
    
    px.diff = diff(px)
    

    .trace("Estimating peaks difference", trace, 1)
    px.diff = px.diff[px.diff >= (median(px.diff) - round(r/2)) & px.diff <= (median(px.diff) + round(r/2))]
    px.diff_exp = mean(px.diff)

    
    stopifnot(!is.na(px.diff_exp))
    
    px_exp = round(seq(0, by=px.diff_exp, length.out=dim_exp))
    px = px[c(1, which(diff(px) > px.diff_exp/4) + 1)]
        
    #if(length(px) < dim_exp) 
    #{
    #    plot(pe, type="l", col="#000000")
    #    lines(pe.mean, type="l", col="#FF000066")
    #    abline(v=px.orig, col="grey")
    #    abline(v=px, col="red")
    #    stop(paste0("Didn't find enough candidate peaks for grid detection (", length(px), " < ", dim_exp, ")"))
    #}
        
    peak.select = function(px_m, px_exp) {
        px_exp.tmp = px_exp
        px_fit = rep(NA, length(px_exp))
        for(i in 1:length(px_exp)) {
            x = px_exp[i]
            fi = which.min(abs(x - px_m))
            px_fit[i] = px_m[fi]
            px_m[fi] = NA
        }
        
        px_fit
    }
    
    px.fit = list()
    for(px.start in px)
        px.fit[[length(px.fit)+1]] = optim(px.start, function(z1.offset, z1.px_exp, z1.px) {
            res = sum(sapply(z1.px, function(z2, z2.offset, z2.px_exp) { min(abs(z2 - z2.px_exp - z2.offset)) }, z2.px_exp=z1.px_exp, z2.offset=z1.offset)^2)
            res
        }, method="SANN", z1.px_exp=px_exp, z1.px=px, lower=0, upper=max(px))
    offset = px.fit[[which.min(sapply(px.fit, function(z) min(z$value)))]]$par 
    
    px.fit_cor = sapply(1:length(px_exp), function(i, ipem, ioffset, ipx.diff_exp, ipx_exp) {
        x = ipx_exp[i] + ioffset
        x.local = (x-floor(ipx.diff_exp/2)+1):(x+floor(ipx.diff_exp/2)-1)
        x.middle = ceiling(length(x.local)/2)
        opt.ctrl = list(maxit = 20000, trace=1000)
        f.opt = function(z, pem) -pem[x.local][round(z)]
        x.opt = x.local[1] + optim(x.middle, f.opt, pem=ipem, method="Brent", control=opt.ctrl, lower=1, upper=length(x.local))$par
        x.opt = ifelse(ipem[x] > ipem[x.opt], x, x.opt)
        
        ipem.rle = rle(ipem)
        ipem.rle.x = c(0, cumsum(ipem.rle$lengths))
        peak.start = ipem.rle.x[last(which(x.opt >= ipem.rle.x))]
        peak.end = ipem.rle.x[first(which(x.opt <= ipem.rle.x))]
        
        x.opt = round((peak.start + peak.end)/2)
        x.opt = ifelse(peak.end >= length(ipem)*0.99, ipem.rle.x[length(ipem.rle.x)-1], x.opt)
        x.opt = ifelse(peak.start <= length(ipem)*0.01, ipem.rle.x[2], x.opt)
        if(all(ipem[x.local] < 0.05)) x.opt = median(x.local)
        
        return(x.opt)
    }, ipem=pe, ioffset=offset, ipx.diff_exp=px.diff_exp, ipx_exp=px_exp)
    
    px.fit_cor2 = round(px.fit_cor)
    
    if(p)
    {
        plot(pe, type="l", main=paste0("Identified ", dim_exp, " peaks in ", margin.i, " margin"))
        lines(pe.mean, col="#FF000066")
        arrows(px.fit_cor2, 1.2, px.fit_cor2, pe[px.fit_cor2]+0.01, angle=rep(15, length(px.fit)), col="#00FF0066", lwd=2)
        segments(px.fit_cor2, pe[px.fit_cor2]-0.01, px.fit_cor2, 0, col="#CCCCCCDD", lwd=0.5)
    }
    
    sort(px.fit_cor2)
}


find.plate = function(img.bw, p=F)
{
    library(PET)
    
    img.bw.dmax = min(dim(img.bw))
    img.bw.square = resize(img.bw, w=img.bw.dmax, h=img.bw.dmax)
    img.bw.square.edges = find.edges(img.bw.square)$img
    sinogram = PET::radon(img.bw.square.edges)$rData
    sinogram.h = hist(sinogram, breaks=2^10, plot=F)
    sinogram.m = length(sinogram)-cumsum(c(0, sinogram.h$counts))
    sinogram.thr = sinogram.h$breaks[which(sinogram.m < img.bw.dmax/4)[1]]
    img.plate = PET::iradon(sinogram > sinogram.thr, img.bw.dmax, img.bw.dmax)$irData
    img.plate = resize(img.plate, w=dim(img.bw)[1], h=dim(img.bw)[2])
    img.plate.edge = bwlabel(erode(.inverse(img.plate > 1), makeBrush(round.odd(img.bw.dmax/200), "box")))
    img.plate.edge.ftrs = computeFeatures.shape(img.plate.edge)
    img.plate.edge.biggest = which.max(img.plate.edge.ftrs[,"s.area"])
    img.plate.edge[!(img.plate.edge %in% img.plate.edge.biggest)] = 0
    img.plate.edge = img.plate.edge > 0
    
    which.frame = function(mask.bulk, d, FUN)
    {
        FUN.s = FUN
        if(FUN=="max") {
            FUN = function(x) ifelse(sum(x > 0), max(which(x > 0)), NA)
        } else {
            FUN = function(x) ifelse(sum(x > 0), min(which(x > 0)), NA)
        }
        
        d.inv = setdiff(1:2, d)
        ndim = dim(mask.bulk)[d]
        ndim.inv = dim(mask.bulk)[d.inv]
        res = data.frame(V1=1:ndim, V2=apply(mask.bulk, d, FUN))
        if(FUN.s == "max") { 
            f = res$V2 > ndim.inv*0.75 
        } else {
            f = res$V2 < ndim.inv*0.25 
        }
        res = res[complete.cases(res) & f,]
        res = median(res$V2, na.rm=T)
        coef = 0.0075 + (1-dim(mask.bulk)[d]/max(dim(mask.bulk))) * 0.05
        res = res + coef * ifelse(FUN.s == "max", -1, 1) * dim(mask.bulk)[d] 
            
        
        res
    }
    
    f.bottom = which.frame(img.plate.edge, 2, "max")
    f.bottom = ifelse(is.na(f.bottom), ncol(img.plate.edge), f.bottom)
    f.top = which.frame(img.plate.edge, 2, "min")
    f.top = ifelse(is.na(f.top), 0, f.top)
    f.right = which.frame(img.plate.edge, 1, "max")
    f.right = ifelse(is.na(f.right), nrow(img.plate.edge), f.right)
    f.left = which.frame(img.plate.edge, 1, "min")
    f.left = ifelse(is.na(f.left), 0, f.left)
    
    img.plate.edge2 = .blank(img.plate.edge)
    img.plate.edge2[f.top:f.bottom,f.left:f.right] = 1
    
    if(p)
    {
        display(normalize(img.bw + 0.5*img.plate.edge +0.2*img.plate.edge2))
    }
    
    img.plate.edge2
}

.find.spaces = function(pe, px, r, p=F)
{
    pe.rle = rle(pe)
    pe.rle = pe.rle$length[pe.rle$values > .95]
    pe.rle = pe.rle[pe.rle > r*1.5]
    s = ifelse(length(pe.rle) > (length(px)+1)/2, median(pe.rle), r*1.5)
    s = pmax(s, r*1.5)
    pe.mean = zoo::rollmean(pe,  round.odd(s), na.pad=T, fill=c(pe[1], NA, last(pe)))
    
    sx.init = zoo::rollmean(px, 2)
    sx.init = round(c(first(sx.init) - mean(diff(sx.init)), sx.init, last(sx.init) + mean(diff(sx.init))))
    sx.init = data.frame(x=sx.init, left=c(0, px+round(r/4)), right=c(px-round(r/4), length(pe)))
    
    sx = apply(sx.init, 1, function(z, zpem) {
        z.ctrl = list(maxit=20000, trace=1000)
        f.opt = function(zz, zzpem) -zzpem[round(zz)]
        z.opt = optim(z["x"], f.opt, zzpem=zpem, method="Brent", control=z.ctrl, lower=z["left"], upper=z["right"])$par
        
        zpem.lsize = 3
        zpem.local = zpem[(z.opt-zpem.lsize):(z.opt+zpem.lsize)]
        z.opt = z.opt + (which.max(zpem.local) - zpem.lsize - 1)
        
        zpem.rle = rle(zpem)
        zpem.rle.x = c(0, cumsum(zpem.rle$lengths))
        peak.start = zpem.rle.x[last(which(z.opt >= zpem.rle.x))]
        peak.end = zpem.rle.x[first(which(z.opt <= zpem.rle.x))]
        z.opt = round((peak.start + peak.end)/2)
        z.opt = ifelse(peak.end >= length(zpem)*0.99, zpem.rle.x[length(zpem.rle.x)-1]+1, z.opt)
        z.opt = ifelse(peak.start <= length(zpem)*0.01, zpem.rle.x[2], z.opt)
        
    }, zpem=pe.mean)
    
    if(p)
    {
        plot(pe, type="l")
        lines(pe.mean, col="#FF000066")
        for(i in 1:nrow(sx.init)) polygon(x=c(sx.init[i, 2:3], sx.init[i, 3:2]), y=c(0, 0, 2, 2), border=NA, col="#FF000022")
        abline(v=sx.init$x, col="red")    
        abline(v=sx, col="green")
    }
    
    sx
}

img.threshold = function(img.bw, format, r)
{
    fdim = format.dim(format)
    
    img.thr = EBImage::thresh(img.bw, r*1.5, r*1.5, 0.005)
    s.radius = round.odd(.005*min(dim(img.bw)))
    img.thr = EBImage::opening(img.thr, EBImage::makeBrush(s.radius, shape="disc"))
    img.thr = EBImage::bwlabel(img.thr)

    img.thr = fillHull(img.thr)
    mask = EBImage::opening(img.thr > 0, EBImage::makeBrush(round.odd(r), shape="disc"))
    mask = EBImage::bwlabel(mask)
    mask[mask %in% unique(c(mask[,c(1, ncol(mask))], mask[c(1, nrow(mask)), ]))] = 0
    
    mask
}

find.grid = function(mask, format, r, trace=1, p=F)
{
    fdim = format.dim(format)
    
    # Find grid
    cy = find.peaks(mask, format, r, margin=2, p=p, trace=trace)
    cx = find.peaks(mask, format, r, margin=1, p=p, trace=trace)

    if(p) {
        .trace("Display detected centers", trace, 1)
        mask.plot = mask
        mask.plot[cx,] = 1
        mask.plot[,cy] = 1
        EBImage::display(mask.plot)
    }
    
    mask.inv = .inverse(mask)
    mask.plot = mask.inv
    mask.inv[c(cx-1,cx,cx+1),] = 0
    mask.inv[,c(cy-1,cy,cy+1)] = 0
    mask.inv[pmin(nrow(mask), max(cx) + round(mean(diff(cx))/2)):nrow(mask),] = 1
    mask.inv[1:pmax(1, min(cx) - round(mean(diff(cx))/2)),] = 1
    mask.inv[,pmin(ncol(mask), max(cy) + round(mean(diff(cy))/2)):ncol(mask)] = 1
    mask.inv[,1:pmax(1, min(cy) - round(mean(diff(cy))/2))] = 1
    mask.inv[pmin(nrow(mask), max(cx) + round(mean(diff(cx)))):nrow(mask),] = 0
    mask.inv[1:pmax(1, min(cx) - round(mean(diff(cx)))),] = 0
    mask.inv[,pmin(ncol(mask), max(cy) + round(mean(diff(cy)))):ncol(mask)] = 0   
    mask.inv[,1:pmax(1, min(cy) - round(mean(diff(cy))))] = 0   
    sx = .find.spaces(rowMeans(mask.inv), cx, r, p=T)
    sy = .find.spaces(colMeans(mask.inv), cy, r, p=T)
        
    if(p)
    {
        mask.plot[sx, ] = 0
        mask.plot[, sy] = 0
        EBImage::display(mask.plot)
    }
       
    mask.sqr = .blank(mask, 1)
    mask.sqr[sx, ] = 0
    mask.sqr[, sy] = 0
    mask.sqr.area = .blank(mask, 0)
    mask.sqr.area[min(sx):max(sx),min(sy):max(sy)] = 1
    mask.sqr = mask.sqr * mask.sqr.area
    mask.sqr = EBImage::bwlabel(mask.sqr > 0)
    
    if(p)
    {  
        displayc(mask.sqr, img.bw)
    }
    
    mask.sqr.max = max(mask.sqr)
    if(mask.sqr.max != prod(fdim)) {
        stop(paste0("Number of detected grid cells is incorrect (found: ", mask.sqr.max, ", expected: ", prod(fdim), ")"))
    }
        
    # Plot results
    ret.ind = expand.grid(col=1:length(cx), row=1:length(cy))
    ret.df = expand.grid(x=cx, y=cy) 
    ret = cbind(ret.ind, ret.df)
    ret$label = 1:nrow(ret)
    
    .trace("Find coordinates of pin squares locations", trace, 1)    
    cx.left = sx[sapply(cx, FUN=function(z, zsx) last(which(zsx < z)), zsx=sx)]
    cx.right = sx[sapply(cx, FUN=function(z, zsx) which(zsx > z)[1], zsx=sx)]
    cy.top = sy[sapply(cy, FUN=function(z, zsy) last(which(zsy < z)), zsy=sy)]
    cy.bottom = sy[sapply(cy, FUN=function(z, zsy) which(zsy > z)[1], zsy=sy)]
    ret$xl = cx.left[match(ret$x, cx)] + 1
    ret$xr = cx.right[match(ret$x, cx)]
    ret$yt = cy.top[match(ret$y, cy)] + 1
    ret$yb = cy.bottom[match(ret$y, cy)]
    
    mask.sqr[sx[2:length(sx)],] = mask.sqr[sx[2:length(sx)]-1,]
    mask.sqr[,sy[2:length(sy)]] = mask.sqr[,sy[2:length(sy)]-1]

    seeds = .blank(mask.sqr, 0)
    seeds[ret$x, ret$y] = mask.sqr[ret$x, ret$y]
    
    spaces = .blank(mask.sqr, 0)
    spaces[sx, sy] = 1
    
    list(data=ret, mask=mask.sqr, seeds=seeds, spaces=spaces)
}


find.radius = function(img, format, p=F)
{
    fdim = format.dim(format)
    
    img.thr = EBImage::thresh(img, 20, 20, 0.005)
    mask = EBImage::opening(img.thr, EBImage::makeBrush(5, shape="disc"))
    mask = EBImage::bwlabel(mask)
    mask.ftrs = EBImage::computeFeatures.shape(mask)
    bad.objects = which(mask.ftrs[,"s.area"] < prod(dim(img)[1:2])/400^2 | mask.ftrs[,"s.radius.min"]*1.5 < mask.ftrs[,"s.radius.max"])
    mask[mask %in% bad.objects] = 0
    
    mask = EBImage::fillHull(mask)
    mask.ftrs = EBImage::computeFeatures.shape(mask)
    
    order.thr = sort(mask.ftrs[,"s.area"])[pmax(nrow(mask.ftrs)-4, 1)]
    order.thr = c(order.thr/10, order.thr)
    mask.ftrs = mask.ftrs[mask.ftrs[,"s.area"] > order.thr[1] & mask.ftrs[,"s.area"] < order.thr[2],]
    objects.order = order(mask.ftrs[,"s.area"], decreasing=T)
    objects = rownames(mask.ftrs)[objects.order]
    
    r.median = median(mask.ftrs[objects, "s.radius.mean"])
    
    if(p)
    {
        mask2 = mask
        mask2[!(mask2 %in% objects)] = 0
        mask2 = channel(mask2, "rgb")
        mask2 = drawCircle(mask2, nrow(mask2)%/%2, ncol(mask2)%/%2, radius=round(r.median), col="#FF9466", fill=T)
        display(mask2)
    }
    
    r.median
}

find.intensity = function(grid.data, img.bw, mode="center", format=max(grid.data$data$col) * max(grid.data$data$row), p=F)
{
    fdim = format.dim(format)
    if(mode == "center")
    {
        mr = 0.1
        minmax = seq(-mr, mr, length.out=15)
        lmax.grid = expand.grid(x=minmax, y=minmax)
        lmax.grid = lmax.grid[sqrt(lmax.grid$x^2 + lmax.grid$y^2) <= mr,]
        lmax.grid = rbind(
            plyr::ddply(lmax.grid, "x", function(z) data.frame(y=min(z["y"]))),
            plyr::ddply(lmax.grid, "x", function(z) data.frame(y=max(z["y"])))
        )
        
        nr = nrow(img.bw)
        lmax = apply(lmax.grid, 1, function(z, d, znr) {
            s = pmin(d$yb-d$yt, d$xr-d$xl)
            ceiling(d$y + z[2]*s)*znr + 
            ceiling(d$x + z[1]*s)
        }, d=grid.data$data, znr=nr)
        lmax.vals = matrix(as.vector(img.bw[as.numeric(lmax)]), nrow=prod(fdim))
        lmax.mean = apply(lmax.vals, 1, function(z) {
            z.dens = density(z)
            z.dens$x[which.max(z.dens$y)]
            #mean(quantile(z, c(.25, .75)))
        })
            
        if(p)
        {
            img.plot = channel(img.bw, "rgb")
            img.plot[rep(as.numeric(lmax), each = 3)] = 1
            img.plot[as.numeric(lmax)] = 0
            img.plot[as.numeric(lmax)] = 0
            display(img.plot)
            
            plot(density(as.numeric(lmax.vals[1,])), col="#00000033")
            for(i in 1:nrow(lmax)) {
                lines(density(as.numeric(lmax.vals[i,])), col="#00000033")
            }
        }
        
        
        return(lmax.mean)
    }
    
    if(mode == "background")
    {
        img.bw
        
        q = .5
        bg = with(grid.data$data, cbind(
            tl = apply(cbind(
                img.bw[yt*nrow(img.bw) + xl], 
                img.bw[(yt+1)*nrow(img.bw) + xl], 
                img.bw[yt*nrow(img.bw) + (xl + 1)], 
                img.bw[(yt+1)*nrow(img.bw) + (xl + 1)]
            ), 1, quantile, probs=q), 
            
            tr = apply(cbind(
                img.bw[yt*nrow(img.bw) + xr],
                img.bw[(yt+1)*nrow(img.bw) + xr],
                img.bw[yt*nrow(img.bw) + (xr - 1)],
                img.bw[(yt+1)*nrow(img.bw) + (xr - 1)]
            ), 1, quantile, probs=q),  
            
            bl = apply(cbind(
                img.bw[yb*nrow(img.bw) + xl], 
                img.bw[(yb-1)*nrow(img.bw) + xl], 
                img.bw[yb*nrow(img.bw) + (xl + 1)],
                img.bw[(yb-1)*nrow(img.bw) + (xl + 1)]
            ), 1, quantile, probs=q), 
            
            br = apply(cbind(
                img.bw[yb*nrow(img.bw) + xr],
                img.bw[(yb-1)*nrow(img.bw) + xr],
                img.bw[yb*nrow(img.bw) + (xr - 1)],
                img.bw[(yb-1)*nrow(img.bw) + (xl - 1)]
            ), 1, quantile, probs=q)
        ))
        
        return(bg)
    }
}



gradient.bg = function(grid.data, img.bw, img.bw.stack, format)
{
    pfdim = dim(img.bw.stack)[3]
    frame.dim = prod(dim(img.bw.stack)[1:2])
    grad.x = row(img.bw.stack[,,1])/nrow(img.bw.stack)
    grad.y = col(img.bw.stack[,,1])/ncol(img.bw.stack)
    grad.br = grad.x*grad.y
    grad.bl = (1-grad.x)*grad.y
    grad.tl = (1-grad.x)*(1-grad.y)
    grad.tr = grad.x*(1-grad.y)
    bg = find.intensity(grid.data, img.bw, "background", format=format)
    
    img.bw.stack_bgbr = img.bw.stack_bgbl = img.bw.stack_bgtr = img.bw.stack_bgtl = img.bw.stack
    img.bw.stack_bgtl[T] = rep(as.numeric(grad.tl), pfdim) * rep(bg[,"tl"], each=frame.dim)
    img.bw.stack_bgtr[T] = rep(as.numeric(grad.tr), pfdim) * rep(bg[,"tr"], each=frame.dim)
    img.bw.stack_bgbl[T] = rep(as.numeric(grad.bl), pfdim) * rep(bg[,"bl"], each=frame.dim)
    img.bw.stack_bgbr[T] = rep(as.numeric(grad.br), pfdim) * rep(bg[,"br"], each=frame.dim)

    img.bw.stack_bg = img.bw.stack_bgtl + img.bw.stack_bgtr + img.bw.stack_bgbl + img.bw.stack_bgbr
    
    img.bw.stack_bg
}

.make.crown.fft.cache <<- list()
make.crown.fft = function(d, idim, thickness=3)
{
    cache.name = paste(d, idim[1], idim[2], thickness)
    cache.copy = .make.crown.fft.cache
    
    wfilter = NULL
    if(!(cache.name %in% names(cache.copy))) {
        wf = makeCrownBrush(d, thickness)
        cf = dim(wf)%/%2
        cx = idim%/%2
        
        wfilter = matrix(0.0, nr=idim[1], nc=idim[2])
        wfilter[(cx[1]-cf[1]):(cx[1]+cf[1]),(cx[2]-cf[2]):(cx[2]+cf[2])] = wf
        wfilter = fft(wfilter)
        cache.copy[[cache.name]] = wfilter
        .make.crown.fft.cache <<- cache.copy
    } else {
        wfilter = cache.copy[[cache.name]]
    }
    
    wfilter
}

find.circles = function(img.edges, format, diameter.bounds=NA)
{
    fdim = format.dim(format)
    s = min(dim(img.edges)[1:2]/rev(fdim))
    
    thickness = 3
    
    diameters = c()
    if(all(is.na(diameter.bounds))) {
        diameters = seq(s, s/10, -thickness)
    } else {
        diameter.bounds = sort(diameter.bounds, decreasing=T)
        diameters = seq(diameter.bounds[1], diameter.bounds[2], -thickness)
    }
    
    diameters = unique(sapply(diameters, round.odd))
    diameters = c(diameters[diameters > 10], 10)
    
    img.edges.fft = fft(img.edges)
    dx = dim(img.edges)
    cx = dx%/%2
    
    .trace("Trying multiple radii values", 1,1)
    img.cicles = Image(0, dim=c(dim(img.edges)[1:2], length(diameters)))
    for (i in 1:length(diameters)) {
        wfilter = make.crown.fft(diameters[i], dim(img.edges), thickness)
        cx = dim(img.edges)%/%2
        
        index1 = c(cx[1]:dx[1], 1:(cx[1]-1))
        index2 = c(cx[2]:dx[2], 1:(cx[2]-1))
        img.cicles[,,i] = Re(fft(img.edges.fft*wfilter, inverse=T)/prod(dim(img.edges)))[index1,index2]
    }
    .trace("Trying multiple radii values [DONE]", 1,1)
    
    img.cicles.sig = which(img.cicles > 0.5, arr.ind=T)
    img.cicles.sig = img.cicles.sig[order(img.cicles.sig[,3], decreasing = T),]
    
    img.centers = .blank(img.edges, 0)
    img.centers[img.cicles.sig[,1], img.cicles.sig[,2]] = diameters[img.cicles.sig[,3]]
    img.centers.lab = bwlabel(closing(img.centers > 0, makeBrush(round.odd(thickness), "disc")))
    img.centers.mf = computeFeatures.moment(img.centers.lab)
    img.centers.bf = computeFeatures.basic(img.centers.lab, img.centers)

    img.circles = .blank(img.edges, 0)
    for(i in 1:nrow(img.centers.mf))
        img.circles = drawCircle(img.circles, x=img.centers.mf[i,"m.cx"], img.centers.mf[i,"m.cy"], radius=img.centers.bf[i,"b.mean"]/2, i, fill=T)
    
    img.circles
}

find.pins = function(img.bw, format, r, ..., img.edges=NULL, img.circles=NULL, p=F, trace=1)
{
    fdim = format.dim(format)
    s = round.odd(.005*min(dim(img.bw)))
    
    .trace("Find grid", trace, 1)
    mask = img.threshold(img.bw, format, r)
    ret.grid = find.grid(mask, r=r, format=format, p=p)
    ret.grid$mask.flat = as.numeric(ret.grid$mask[ret.grid$mask > 0])
    
    .trace("Create stacked images (bw/class/seed)", trace, 1)
    img.bw.stack = stackObjects2(ret.grid$mask, img.bw)
    frame.dim = prod(dim(img.bw.stack)[1:2])
    img.bw.stack_class = EBImage::Image(rep(1:prod(fdim), each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_bg = gradient.bg(ret.grid, img.bw, img.bw.stack, format=format)
    mask.stack = stackObjects2(ret.grid$mask, mask)
    
    if(is.null(img.circles)) {
        .trace("Estimate perimiters for Sobel thresholding", trace, 1)
        img.bw.stack_bgmask = (abs(img.bw.stack - img.bw.stack_bg) > 0.1) | (mask.stack > 0)
        img.bw.stack_bgmask[is.na(img.bw.stack)] = 0
        img.bw.stack_bgmask = closing(img.bw.stack_bgmask, makeBrush(round.odd(2*s), "disc"))
        img.bw.stack_bgmask = fillHull(img.bw.stack_bgmask)
        img.bw.stack_bgmask = (img.bw.stack_bgmask > 0) * img.bw.stack_class
        img.bw.stack_bgmask[is.na(img.bw.stack)] = NA
        img.bw.stack_bgmask.tile = tile2(img.bw.stack_bgmask, ret.grid$mask)
    } else {
        img.bw.stack_bgmask.tile = img.circles > 0
    }
    
    img.bw.stack_bgmask.tile = propagate(img.bw, ret.grid$seeds, mask=img.bw.stack_bgmask.tile > 0)
    img.bw.stack_bgmask.tile.ftrs = computeFeatures.shape(img.bw.stack_bgmask.tile)
    
    img.edges.mask = NULL
    r.max = pmin(mean(ret.grid$data$xr - ret.grid$data$xl), mean(ret.grid$data$yb-ret.grid$data$yt))/2
    if((r.max-r)/r.max/2 > 0.1)
    {
        perimeters = img.bw.stack_bgmask.tile.ftrs[as.character(1:prod(fdim)),"s.perimeter"]
        
        .trace("Sobel edge detection", trace, 1)
        if(is.null(img.edges))
            img.edges = find.edges(medianFilter(img.bw, 3))$img
        
        .trace("Sobel thresholding", trace, 1)
        img.edges.wells = as.numeric(img.edges[ret.grid$mask > 0])
        img.edges.label = aggregate(img.edges.wells, by=list(label=ret.grid$mask.flat), function(z) z)
        img.edges.stack = img.bw.stack
        img.edges.stack[!is.na(img.edges.stack)] = unlist(img.edges.label$x)
        img.edges.thr = sapply(1:prod(fdim), function(z, z.img, z.per) {
            z.v =  z.img$x[[which(z.img$label==z)]]
            z.h = hist(z.v, breaks=2^10, plot=F)
            z.m = length(z.v)-cumsum(c(0, z.h$counts))
            z.h$breaks[which(z.m < 3*z.per[z])[1]]
        }, z.img=img.edges.label, z.per=perimeters)
        img.edges.stack.thr = EBImage::Image(rep(img.edges.thr, each=frame.dim), dim=dim(img.edges.stack))
        img.edges.stack.mask = (img.edges.stack > img.edges.stack.thr)
        
        .trace("Unnecessary shit", trace, 1)
        img.edges.mask = tile2(img.edges.stack.mask, ret.grid$mask)
        if(!is.null(img.circles)) img.edges.mask = img.edges.mask | (img.bw.stack_bgmask.tile > 0)
        img.edges.mask = fillHull(img.edges.mask > 0)
        img.edges.mask = propagate(img.bw, ret.grid$seeds, mask=img.edges.mask)
        img.edges.mask = closing(img.edges.mask > 0, makeBrush(round.odd(r), "disc"))
        img.edges.mask = propagate(img.bw, ret.grid$seeds, mask=img.edges.mask)
    } else {
        img.edges.mask = img.bw.stack_bgmask.tile
    }
    
    if(p)
    {
        displayc(img.edges.mask, img.bw)
    }
    
    list(mask.pins=img.edges.mask, mask.sqr=ret.grid$mask, data=ret.grid$data)
}


parse.file = function()
{
    files = list.files("benchmark2/BY/", pattern=".*\\.JPG$", full.names=T)
    
    img.big = EBImage::readImage(files[7])
    img = EBImage::resize(img.big, nrow(img.big)/3, ncol(img.big)/3)
    plate.mask = find.plate(rgb2hsv(img)[,,3], p=T)
    
    for(file in files[1:4])
    {
        .trace(paste0(which(file==files), " / ", length(files), ": ", file), 1, 1)
        format = as.numeric(gsub("(\\d+)_.*", "\\1", basename(file)))
        img.big = EBImage::readImage(file)
        img = EBImage::resize(img.big, nrow(img.big)/3, ncol(img.big)/3)
        img = crop(plate.mask, img)
        img.hsv = rgb2hsv(img)
        img.bw = img.hsv[,,3]

        img.edges = find.edges(medianFilter(img.bw, 3))$img[,,1]
        
        r = NULL; img.circles.3=NULL;
        if(format == 1536)
        {
            r = find.radius(img.bw, format, p=F)
        } else {
            img.circles = find.circles(img.edges, format) #, diameter.bounds=c(r/2, r*4))
            r = median(computeFeatures.shape(bwlabel(img.circles))[,"s.radius.mean"])
            img.circles.3 = erode(img.circles > 0, makeBrush(3, "disc"))
        }
        
        ret = find.pins(img.bw, format, img.edges=img.edges, img.circles=img.circles.3, r=r, p=T)
        displayc(ret$mask.pins, img)
        
        .trace("Calculate color features", 1, 1)
        img.ftrs = EBImage::computeFeatures.shape(ret$mask.pins)
        img.ftrs = cbind(img.ftrs, EBImage::computeFeatures.moment(ret$mask.pins))
        img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img, c("red", "green", "blue")))
        #img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img.cmyk, c("cyan", "magenta", "yellow", "black")))
        img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img.hsv, c("hue", "saturation", "value")))
        img.ftrs$label = as.numeric(rownames(img.ftrs))
        
        color.ftrs = colnames(img.ftrs)[grepl("red|green|blue|cyan|magenta|yellow|black", colnames(img.ftrs)) & grepl("mean", colnames(img.ftrs))]
        img.pca = princomp(img.ftrs[,color.ftrs])
        .trace(paste0("First compoment explains ", round(100*(img.pca$sdev[1]^2 / sum(img.pca$sdev^2)), 1), "% variance"), 1,1)
        img.ftrs$pca1 = img.pca$scores[,1]
        
        .trace("Writing results", 1,1)
        png(paste0(file, "_features.png"), width=nrow(ret$mask.pins), height=ncol(ret$mask.pins))
        displayf(ret$mask.pins, img.ftrs, variable="pca1", method="raster")
        dev.off()
        
        png(paste0(file, "_contour.png"), width=nrow(ret$mask.pins), height=ncol(ret$mask.pins))
        displayc(ret$mask.pins, img, method="raster")
        dev.off()
        
        data.export = merge(ret$data, img.ftrs, by="label", all.x=T)
        write.table(data.export, file=paste0(file, "_data.tab"), quote=F, row.names=F, sep="\t", na="")
    }
}

