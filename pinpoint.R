rgb2cmyk = function(img2)
{
    img2 = img2[,,1:3]
    img2.1 = img2[,,1]
    img2.2 = img2[,,2]
    img2.3 = img2[,,3]
    
    K = img2.1
    K[img2.2 > K] = img2.2[img2.2 > K]
    K[img2.3 > K] = img2.3[img2.3 > K]

    K = 1 - K
    C = (1-img2.1-K) / (1-K)
    M = (1-img2.2-K) / (1-K)
    Y = (1-img2.3-K) / (1-K)
    
    
    library(abind)
    abind(C, M, Y, K, along=3)
}

computeFeatures.basic2 = function(mask2, img2, layers="rgb")
{
    res = NULL
    if(length(dim(img)) == 2)
    {   
        res = as.data.frame(EBImage::computeFeatures.basic(mask2, img2))
    } else {
        if(length(layers) == 1) layers.v = unlist(strsplit(layers, ""))
        else layers.v = layers
        
        for(i in 1:dim(img)[3])    
        {
            res.tmp = as.data.frame(EBImage::computeFeatures.basic(mask2, img2[,,i]))
            colnames(res.tmp) = paste0(layers.v[i], ".", colnames(res.tmp))
            
            if(is.null(res)) res = res.tmp
            else res = cbind(res, res.tmp)
        }
    }   
    
    res
}

rgb2hsv = function(v)
{
    a = grDevices::rgb2hsv(as.numeric(v[,,1]), as.numeric(v[,,2]), as.numeric(v[,,3]), maxColorValue=1)
    v.hsv = v
    v.hsv[,,1] = a[1,]
    v.hsv[,,2] = a[2,]
    v.hsv[,,3] = a[3,]
    EBImage::colorMode(v.hsv) = EBImage::Grayscale
    
    v.hsv
}



displayc = function(mask, img, method="browser")
{
    if(length(dim(img)) == 2 || dim(img)[3] == 1)
    {
        img = EBImage::rgbImage(img, img, img)
    }

    EBImage::display(EBImage::paintObjects(mask, img, col="#FF000066"), method=method)
}

displayb = function(mask, img, method="browser")
{
    img.na = img
    img.na[mask==0] = NA
    img.na = normalize(img.na)
    img.na[is.na(img.na)] = 0
    EBImage::display(img.na, method=method)
}

displayf = function(mask2, ftrs, variable="variable", method="browser")
{
    if(is.null(ftrs$label)) 
        ftrs$label = as.numeric(rownames(ftrs))
    
    mask2.data = ftrs[match(mask2[mask2!=0], ftrs$label), variable]
    
    mask2[mask2!=0] = mask2.data
    mask3 = EBImage::Image(abind(1/3+mask2/2, 1/4+mask2/2, 1/4+mask2/2, along=3), colormode=EBImage::Color)
    mask3[,,1] = NA
    mask3[mask2==0] = NA 
    
    mask3 = EBImage::normalize(mask3)
    mask3[is.na(mask3)] = 0
    EBImage::display(mask3, method=method)
}

format.dim = function(format, dim=c(1,2))
{
    d = format
    if(length(format) == 1)
    {
        formats = list()
        formats[["1700"]] = c(34, 50)
        formats[["1536"]] = c(32, 48)
        formats[["884"]] = c(26, 34)
        formats[["768"]] = c(24, 32)
        formats[["468"]] = c(18, 26)
        formats[["384"]] = c(16, 24)
        formats[["140"]] = c(10, 14)
        formats[["96"]] = c(8, 12)
        d = formats[[as.character(format)]]
    }
    
    d[dim]
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
        
    if(length(px) < dim_exp) 
    {
        plot(pe, type="l", col="#000000")
        lines(pe.mean, type="l", col="#FF000066")
        abline(v=px.orig, col="grey")
        abline(v=px, col="red")
        stop(paste0("Didn't find enough candidate peaks for grid detection (", length(px), " < ", dim_exp, ")"))
    }
    
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
    
    
    
    fit.offsets = 1:round(length(pe) / 4)
    fit.candidates = sapply(fit.offsets, function(o) px - o)
    px.fit = apply(fit.candidates, 2, peak.select, px_exp=px_exp)
    px.fit_err = apply(px.fit, 2, function(x) sum((px_exp-x)^2))
    
    offset = fit.offsets[which.min(px.fit_err)] 
    
    plot(pe, type="l")
    lines(pe.mean, type="l", col="grey")
    abline(v=px_exp + offset)
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


find.plate = function(img.bw, r, p=F)
{
    s = 0.01*min(dim(img.bw))
    img.thr = EBImage::thresh(img.bw,  2*s, 2*s, 0.005)
    img.thr = EBImage::opening(img.thr, EBImage::makeBrush(round.odd(s), shape="disc"))
    mask = EBImage::bwlabel(img.thr)
    mask.orig = mask
    
    mask.ftrs = as.data.frame(EBImage::computeFeatures.shape(mask))
    
    # remove frame
    bulk.objects = which(mask.ftrs$s.radius.max > r*4)
    mask.bulk = mask
    mask.bulk[!(mask %in% bulk.objects)] = 0
    
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
        res = res + 0.0125 * ifelse(FUN.s == "max", -1, 1) * dim(mask.bulk)[d]
        
        res
    }
    
    f.bottom = which.frame(mask.bulk, 2, "max")
    f.bottom = ifelse(is.na(f.bottom), ncol(mask.bulk), f.bottom)
    f.top = which.frame(mask.bulk, 2, "min")
    f.top = ifelse(is.na(f.top), 0, f.top)
    f.right = which.frame(mask.bulk, 1, "max")
    f.right = ifelse(is.na(f.right), nrow(mask.bulk), f.right)
    f.left = which.frame(mask.bulk, 1, "min")
    f.left = ifelse(is.na(f.left), 0, f.left)

    mask2 = mask
    mask2[T] = 0
    mask2[(f.top+r):(f.bottom-r),(f.left+r):(f.right-r)] = 1
    
    if(p)
    {
        mask.plot = channel(mask.orig, "rgb")
        mask.plot[mask.plot %in% bulk.objects] = rep(c(1, 0.58, 0.4), each=sum(mask.plot %in% bulk.objects)/3)
        mask.plot = mask.plot + 0.5*channel(mask2, "rgb")    
        display(mask.plot)
    }
    
    mask2
}

round.odd = function (x) 
{
    x = floor(x)
    if (x%%2 == 0) 
        x = x + 1
    x
}
    
.inverse = function(mask)
{
    mask.inv = mask
    mask.inv[T] = 0
    mask.inv[mask == 0] = 1
    mask.inv
}

.trace = function(text, level, req)
{
    if(level >= req)
        writeLines(text)
}

.find.spaces = function(pe, px, p=F)
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

find.grid = function(img.bw, format, r=NULL, trace=1, p=F, ..., plate.mask=NULL)
{
    if(is.null(r)) r = find.radius(img.bw, format, p=p)
    
    fdim = format.dim(format)
    
    .trace("Convert to black and white", trace, 1)
    img.thr = EBImage::thresh(img.bw, r*1.5, r*1.5, 0.005)
    s.radius = round.odd(.005*min(dim(img.bw)))
    img.thr = EBImage::opening(img.thr, EBImage::makeBrush(s.radius, shape="disc"))
    img.thr = EBImage::bwlabel(img.thr)
    if(!is.null(plate.mask))
    {
        .trace("Remove objects outside plate pins area", trace, 1)
        frame.objects = unique(as.vector(img.thr * (plate.mask == 0)))
        img.thr[img.thr %in% frame.objects] = 0
        img.thr = reenumerate(img.thr)
    }
    
    img.thr = fillHull(img.thr)
    mask = EBImage::opening(img.thr > 0, EBImage::makeBrush(round.odd(r), shape="disc"))
    mask = EBImage::bwlabel(mask)
    
    
    # Find grid
    .trace("Find pins centers", trace, 1)
    cy = find.peaks(mask, format, r, margin=2, p=p, trace=trace)
    cx = find.peaks(mask, format, r, margin=1, p=p, trace=trace)


    if(p) {
        .trace("Display detected centers", trace, 1)
        mask.plot = mask
        mask.plot[cx,] = 1
        mask.plot[,cy] = 1
        EBImage::display(mask.plot)
    }
    
    
    
    .trace("Find spaces between pins", trace, 1)
    mask.inv = .inverse(mask)
    mask.inv[(max(cx) + round(mean(diff(cx))/2)):nrow(mask),] = 1
    mask.inv[1:(min(cx) - round(mean(diff(cx))/2)),] = 1
    mask.inv[,(max(cy) + round(mean(diff(cy))/2)):ncol(mask)] = 1
    mask.inv[,1:(min(cy) - round(mean(diff(cy))/2))] = 1 
    mask.inv[(max(cx) + round(mean(diff(cx)))):nrow(mask),] = 0
    mask.inv[1:(min(cx) - round(mean(diff(cx)))),] = 0
    mask.inv[,(max(cy) + round(mean(diff(cy)))):ncol(mask)] = 0   
    mask.inv[,1:(min(cy) - round(mean(diff(cy))))] = 0   
    sx = .find.spaces(rowMeans(mask.inv), cx, p=T)
    sy = .find.spaces(colMeans(mask.inv), cy, p=T)
         
    if(p)
    {
        .trace("Display detected spaces", trace, 1)
        mask.plot = mask.inv
        mask.plot[sx, ] = 0
        mask.plot[, sy] = 0
        EBImage::display(mask.plot)
    }
    
    .trace("Find squares for around pins", trace, 1)    
    mask.sqr = mask
    mask.sqr[T] = 1
    mask.sqr[sx, ] = 0
    mask.sqr[, sy] = 0
    mask.sqr.area = mask.sqr
    mask.sqr.area[T] = 0
    mask.sqr.area[min(sx):max(sx),min(sy):max(sy)] = 1
    mask.sqr = mask.sqr * mask.sqr.area
    mask.sqr = EBImage::bwlabel(mask.sqr > 0)
    
    if(p)
    {
        .trace("Plot squares for around pins", trace, 1)    
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

    
    list(data=ret, mask=mask.sqr)
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
    
    objects.order = order(mask.ftrs[,"s.area"], decreasing=T)
    objects = rownames(mask.ftrs)[objects.order]
    objects = objects[5:pmin(prod(fdim)/2, nrow(mask.ftrs)/3)]
    
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

find.intensity = function(grid.data, img.bw, mode="center", format=max(grid.data$col) * max(grid.data$row), p=F)
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

find.pins = function(img.bw, format, r=NULL, p=F, ..., plate.mask=NULL)
{
    fdim = format.dim(format)
    
    if(is.null(r)) r = find.radius(img.bw, format, p=p)
    
    ret.grid = find.grid(img.bw, r=r, format=format, p=p, plate.mask=plate.mask)
    
    f = makeBrush(round.odd(r/2), shape='disc', step=FALSE)
    f = f/sum(f)
    img.bw.low = filter2(img.bw, f)
    lmax = find.intensity(ret, img.bw.low, "center", format=format)
    bg = find.intensity(ret, img.bw.low, "background", format=format)
    bg.avg = rowMeans(bg)

    img.bw.stack = EBImage::stackObjects(ret.grid$mask, img.bw, bg.col=-1)
    img.bw.stack[img.bw.stack < 0] = NA
    img.bw.stack.orig = img.bw.stack
    frame.dim = prod(dim(img.bw.stack)[1:2])
    
    grad.x = row(img.bw.stack[,,1])/nrow(img.bw.stack)
    grad.y = col(img.bw.stack[,,1])/ncol(img.bw.stack)
    grad.br = grad.x*grad.y
    grad.bl = (1-grad.x)*grad.y
    grad.tl = (1-grad.x)*(1-grad.y)
    grad.tr = grad.x*(1-grad.y)
    img.bw.stack_bgtl = EBImage::Image(rep(as.numeric(grad.tl), prod(fdim)) * rep(bg[1:prod(fdim),"tl"], each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_bgtr = EBImage::Image(rep(as.numeric(grad.tr), prod(fdim)) * rep(bg[1:prod(fdim),"tr"], each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_bgbl = EBImage::Image(rep(as.numeric(grad.bl), prod(fdim)) * rep(bg[1:prod(fdim),"bl"], each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_bgbr = EBImage::Image(rep(as.numeric(grad.br), prod(fdim)) * rep(bg[1:prod(fdim),"br"], each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_bg = img.bw.stack_bgtl + img.bw.stack_bgtr + img.bw.stack_bgbl + img.bw.stack_bgbr
    img.bw.stack[is.na(img.bw.stack)] = img.bw.stack_bg[is.na(img.bw.stack)]
    
    img.bw.stack_max = EBImage::Image(rep(lmax[1:prod(fdim)], each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_class = EBImage::Image(rep(1:prod(fdim), each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_hat = EBImage::selfcomplementaryTopHatGreyScale(img.bw.stack, EBImage::makeBrush(round.odd(r/2)))
    img.bg.stack = EBImage::Image(rep(bg.avg > lmax, each=frame.dim), dim=dim(img.bw.stack))
    
    
    f = makeBrush(round.odd(r/1.5), shape='disc', step=FALSE)
    img.bw.stack_low1.5 = filter2(img.bw.stack, f/sum(f))
    f = makeBrush(round.odd(r/2), shape='disc', step=FALSE)
    img.bw.stack_low2 = filter2(img.bw.stack, f/sum(f))
    f = makeBrush(round.odd(r/3), shape='disc', step=FALSE)
    img.bw.stack_low3 = filter2(img.bw.stack, f/sum(f))
    img.bw.stack_loc = abs(img.bw.stack - img.bw.stack_low3) > 0.05
    img.bw.stack_loc = closing(img.bw.stack_loc, makeBrush(round.odd(r/1.5)))
    img.bw.stack_loc = fillHull(img.bw.stack_loc)
    
    display(tile(img.bw.stack_max, 12))
    display(tile(img.bw.stack_bg, 12))
    display(tile(img.bw.stack, 12))
    display(tile(img.bw.stack_loc, 12))
    
    mask.stack.a = !is.na(img.bw.stack.orig)
    mask.stack.b = img.bw.stack > img.bw.stack_bg + 0.1*(img.bw.stack_max - img.bw.stack_bg)
    mask.stack.c = abs(img.bw.stack_max - img.bw.stack_bg) / img.bw.stack_bg > 1/4

    mask.stack = mask.stack.a & img.bw.stack_loc
    mask.stack = EBImage::bwlabel(mask.stack)
    mask.stack = EBImage::fillHull(mask.stack)
    mask.stack = EBImage::opening(mask.stack > 0, EBImage::makeBrush(round.odd(r/2)))
    mask.stack = (mask.stack > 0) * img.bw.stack_class
    
    mask = ret.grid$mask
    mask[T] = 0    
    mask.sqr.vector = which(ret.grid$mask>0)[order(as.numeric(ret.grid$mask[ret.grid$mask>0]))]
    mask[mask.sqr.vector] = mask.stack[!is.na(img.bw.stack.orig)]
    
    seeds = mask
    seeds[T] = 0
    seeds[ret$data$x, ret$data$y] = mask[ret$data$x, ret$data$y]
    mask.prop = EBImage::propagate(img.bw, seeds, mask)
        
    if(p)
    {
        displayc(EBImage::tile(mask.stack, fg.col="#000000", nx=fdim[2]), EBImage::tile(img.bw.stack, fg.col=gray(mean(bg)), nx=fdim[2]))
        displayc(mask.prop, img.bw)
    }
    
    list(mask.pins=mask.prop, mask.sqr=ret.grid$mask, data=ret$data)
}

autoconfig.color = function()
{
    files = list.files("benchmark2/p28_X200_96/", pattern=".*\\.JPG$", full.names=T)
    images.val = data.frame()
    for(file in files)
    {
        writeLines(paste0(which(file==files), " / ", length(files), ": ", file))
        img.big = EBImage::readImage(file)
        img.sample = img.big[round(0.4*nrow(img.big)):round(0.6*nrow(img.big)), round(0.4*ncol(img.big)):round(0.6*ncol(img.big)),]
        img.val = cbind(as.numeric(img.hsv[,,1]), as.numeric(img.hsv[,,2]), as.numeric(img.hsv[,,3]), as.numeric(img[,,1]), as.numeric(img[,,2]), as.numeric(img[,,3]), as.numeric(img.cmyk[,,1]), as.numeric(img.cmyk[,,2]), as.numeric(img.cmyk[,,3]), as.numeric(img.cmyk[,,3]))
        images.val = rbind(images.val, img.val)
    }
    
    images.val.s = images.val[sample.int(nrow(images.val), 10000),]
    colnames(images.val.s) = c("h", "s", "v", "r", "g", "b", "c", "m", "y", "k")
    img.val.pca = princomp(images.val.s)
    writeLines(paste0("First compoment explains ", round(100*(img.val.pca$sdev[1]^2 / sum(img.val.pca$sdev^2)), 1), "% variance"))
    pc = 9
    img.pca = 
        img.val.pca$loadings[1,pc] * img.hsv[,,1] + img.val.pca$loadings[2,pc] * img.hsv[,,2] + img.val.pca$loadings[3,pc] * img.hsv[,,3] +
        img.val.pca$loadings[4,pc] * img[,,1] + img.val.pca$loadings[5,pc] * img[,,2] + img.val.pca$loadings[6,pc] * img[,,3] +
        img.val.pca$loadings[7,pc] * img.cmyk[,,1] + img.val.pca$loadings[8,pc] * img.cmyk[,,2] + img.val.pca$loadings[9,pc] * img.cmyk[,,3] + img.val.pca$loadings[10,pc] * img.cmyk[,,4]

    img.pca = 2*img.hsv[,,2]
    display(normalize(img.pca))
    a = 180:(nrow(img)-160)
    b = 110:(ncol(img)-110)
    
   # img.bw = 
    img.cmyk[a,b,3]

    par(mfrow=c(3,4))
    plot(density(img.cmyk[a,b,1]), main="c", col="cyan")
    lines(density(as.numeric(img[a,b,1])), main="r", col="red")
    
    lines(density(img.cmyk[a,b,2]), main="m", col="magenta")
    lines(density(img.cmyk[a,b,4]), main="k", col="black")
    lines(density(as.numeric(img[a,b,1])), main="r", col="red")
    lines(density(as.numeric(img[a,b,2])), main="g", col="green")
    lines(density(as.numeric(img[a,b,3])), main="b", col="blue")
    
    
    
    #lines(density(img.cmyk[a,b,3]), main="y", col="yellow")
    #lines(density(img.hsv[a,b,1]), main="h", col="#AA33FF")
    #lines(density(img.hsv[a,b,2]), main="s", col="#CC66FF")
    #lines(density(img.hsv[a,b,3]), main="v", col="#EE99FF")
    par(mfrow=c(1,1))
    
    display(normalize(img[a,b,1]))
    
    
    display(normalize(img.cmyk[a,b,4] * img.hsv[a,b,3]))
    
    - 0.5*img.cmyk[,,1] - 0.3*img.cmyk[,,3] - 0.9*img.cmyk[,,4]
    
}


parse.file = function()
{
    files = list.files("benchmark2/p28_X200_96/", pattern=".*\\.JPG$", full.names=T)
    for(file in files)
    {
        writeLines(paste0(which(file==files), " / ", length(files), ": ", file))
        #file = "benchmark/1536_4day_2X.JPG"
        
        format = as.numeric(gsub("(\\d+)_.*", "\\1", basename(file)))
        # Read file
        img.big = EBImage::readImage(file)
        img = EBImage::resize(img.big, nrow(img.big)/3, ncol(img.big)/3)
        img.hsv = rgb2hsv(img)
        img.cmyk = rgb2cmyk(img)
        #img.bw = EBImage::channel(img, "grey")
        img.bw = img.hsv[,,3]
        
        r = find.radius(img.bw, format, p = T)
        plate.mask = find.plate(img.bw, r, T)
        ret = find.pins(img.bw, format, r=r, p=T, plate.mask=plate.mask)
        displayc(ret$mask.pins, img)
        
        img.ftrs = EBImage::computeFeatures.shape(ret$mask.pins)
        img.ftrs = cbind(img.ftrs, EBImage::computeFeatures.moment(ret$mask.pins))
        img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img, c("red", "green", "blue")))
        img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img.cmyk, c("cyan", "magenta", "yellow", "black")))
        img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img.hsv, c("hue", "saturation", "value")))
        img.ftrs$label = as.numeric(rownames(img.ftrs))
        
        color.ftrs = colnames(img.ftrs)[grepl("red|green|blue|cyan|magenta|yellow|black", colnames(img.ftrs)) & grepl("mean", colnames(img.ftrs))]
        head(img.ftrs[,color.ftrs])
        img.pca = princomp(img.ftrs[,color.ftrs])
        writeLines(paste0("First compoment explains ", round(100*(img.pca$sdev[1]^2 / sum(img.pca$sdev^2)), 1), "% variance"))
        img.ftrs$pca1 = img.pca$scores[,1]
        
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