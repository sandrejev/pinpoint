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

    #K = 1 - apply(img2, 1:2, max)
    
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
        res = as.data.frame(computeFeatures.basic(mask2, img2))
    } else {
        if(length(layers) == 1) layers.v = unlist(strsplit(layers, ""))
        else layers.v = layers
        
        for(i in 1:dim(img)[3])    
        {
            res.tmp = as.data.frame(computeFeatures.basic(mask2, img2[,,i]))
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
    colorMode(v.hsv) = Grayscale
    
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
    mask3 = Image(abind(1/3+mask2/2, 1/4+mask2/2, 1/4+mask2/2, along=3), colormode=Color)
    mask3[,,1] = NA
    mask3[mask2==0] = NA 
    
    mask3 = normalize(mask3)
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


find.peaks = function(mask, format, r, margin=1, p=F)
{
    # The row/collumn format is inverse of x/y format used by EBImage
    margin.i = setdiff(1:2, margin)
    dim_exp = format.dim(format, margin.i)
    
    library(zoo)
    library(quantmod)
    pe = apply(mask>0, margin, mean, na.rm=T)
    pe.mean = zoo::rollmean(pe, round.odd(r*1.5))
    px = quantmod::findPeaks(pe.mean)
    px.orig = px
    
    px.diff = diff(px)
    px.diff = px.diff[px.diff >= (median(px.diff) - round(r/2)) & px.diff <= (median(px.diff) + round(r/2))]
    px.diff_exp = mean(px.diff)
    px_exp = round(seq(0, by=px.diff_exp, length.out=dim_exp))
    
    px = px[c(1, which(diff(px) > px.diff_exp/4) + 1)]

    if(length(px) < dim_exp) 
    {
        plot(pe.mean, type="l")
        abline(v=px.orig, col="grey")
        abline(v=px, col="red")
        stop("Didn't find enough candidate peaks for grid detection")
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
    
    px.fit_cor = sapply(1:length(px_exp), function(i, ipem, ioffset, ipx.diff_exp, ipx_exp) {
        x = ipx_exp[i] + ioffset
        x.local = (x-floor(ipx.diff_exp/2)+1):(x+floor(ipx.diff_exp/2)-1)
        x.middle = ceiling(length(x.local)/2)
        opt.ctrl = list(maxit = 20000, trace=1000)
        f.opt = function(z, pem) -pem[x.local][round(z)]
        x.opt = x.local[1] + optim(x.middle, f.opt, pem=ipem, method="Brent", control=opt.ctrl, lower=1, upper=length(x.local))$par
        x.opt = ifelse(ipem[x] > ipem[x.opt], x, x.opt)
        
        return(x.opt)
    }, ipem=pe, ioffset=offset, ipx.diff_exp=px.diff_exp, ipx_exp=px_exp)
    
    px.fit_cor2 = sapply(1:length(px.fit_cor), function(i) {
        x = px.fit_cor[i]
        x.local = (x-round(r/3)):(x+round(r/3))
        x + which.max(pe[x.local]) - which(x.local==x)
    })
    
    # Move last and first peak tighter together. This happens sometimes in case
    # of tie result
    pe.rle = with(rle(pe), lengths[c(1, length(lengths))])
    lst = px.fit_cor2[length(px.fit_cor2)]
    px.fit_cor2[length(px.fit_cor2)] = ifelse(lst > length(pe) - pe.rle[2] + 1, length(pe) - pe.rle[2] + 1, lst)
    px.fit_cor2[1] = ifelse(px.fit_cor2[1] < pe.rle[1], pe.rle[1], px.fit_cor2[1])
    
    if(p)
    {
        plot(pe, type="l", main=paste0("Identified ", dim_exp, " peaks in ", margin.i, " margin"))
        lines(pe.mean, col="#FF000066")
        arrows(px.fit_cor2, 1.2, px.fit_cor2, pe[px.fit_cor2]+0.01, angle=rep(15, length(px.fit)), col="#00FF0066", lwd=2)
        segments(px.fit_cor2, pe[px.fit_cor2]-0.01, px.fit_cor2, 0, col="#CCCCCCDD", lwd=0.5)
    }
    
    sort(px.fit_cor2)
}


find.plate = function(img.bw, r=find.radius(img.bw), ..., mask=NULL)
{
    if(is.null(mask)) {
        img.thr = EBImage::thresh(img.bw, r*2, r*2, 0.005)
        img.thr = EBImage::opening(img.thr, EBImage::makeBrush(round.odd(r), shape="disc"))
        mask = EBImage::bwlabel(img.thr)
    }
    
    mask.ftrs = as.data.frame(EBImage::computeFeatures.shape(mask))
    
    # remove frame
    bulk.objects = which(mask.ftrs$s.area > 16*pi*(r^2))
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
        
        res
    }
    
    f.bottom = which.frame(mask.bulk, 2, "max")
    f.top = which.frame(mask.bulk, 2, "min")
    f.right = which.frame(mask.bulk, 1, "max")
    f.left = which.frame(mask.bulk, 1, "min")

    
    mask2 = mask
    mask2[T] = 0
    mask2[(f.top+r):(f.bottom-r),(f.left+r):(f.right-r)] = 1
    
    mask2
}

round.odd = function (x) 
{
    x = floor(x)
    if (x%%2 == 0) 
        x = x + 1
    x
}
    

#find.borders = function(x) {
#    x.diff = diff(x)
#    x.min = min(x) - round(max(x.diff))
#    x.max = max(x) + round(max(x.diff))
#    round(c(x.min, zoo::rollmean(x, 2), x.max))
#}

.find.spaces = function(pe0, r, x)
{
    r2 = round(max(diff(x)))
    optmin = function(pe) { 
        pe = pe > 0; 
        round(optim(.75*r2, function(x) pe[round(x)], method="Brent", lower=1, upper=length(pe))$par) 
    }    
    
    m1 = min(x); m2 = max(x);
    c(m1-optmin(rev(pe0[(m1-r2):m1])), x, m2 - optmin(pe0[m2:(m2+r2)]))
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

find.grid = function(img.bw, format, r=find.radius(img.bw), trace=1, p=F)
{
    fdim = format.dim(format)
    
    .trace("Convert to black and white", trace, 1)
    img.thr = EBImage::thresh(img.bw, r*2, r*2, 0.005)
    img.thr = EBImage::opening(img.thr, EBImage::makeBrush(round.odd(r), shape="disc"))
    mask = EBImage::bwlabel(img.thr)
    
    # Remove frame
    .trace("Detect and remove frame around pins", trace, 1)
    plate.mask = find.plate(mask=mask, r=r)
    frame.objects = unique(as.vector(mask * (plate.mask == 0)))
    mask[mask %in% frame.objects] = 0

    # Find grid
    .trace("Find pins centers", trace, 1)
    cy = find.peaks(mask, format, r, margin=2, p=p)
    cx = find.peaks(mask, format, r, margin=1, p=p)
    cy = round(cy)
    cx = round(cx)

    if(p) {
        .trace("Display detected centers", trace, 1)
        mask.plot = mask
        mask.plot[cx,] = 1
        mask.plot[,cy] = 1
        EBImage::display(mask.plot)
    }
    
    .trace("Find spaces between pins", trace, 1)
    mask.inv = .inverse(mask)
    mask.inv[(max(cx) + 2*r):nrow(mask),] = 1
    mask.inv[1:(min(cx) - 2*r),] = 1
    mask.inv[,(max(cy) + 2*r):ncol(mask)] = 1
    mask.inv[,1:(min(cy) - 2*r)] = 1    
    sy = find.peaks(mask.inv, fdim+1, r, margin=2, p=p)
    sx = find.peaks(mask.inv, fdim+1, r, margin=1, p=p)
    #sy[1] = 2*sy[2] - sy[3]; sx[1] = 2*sx[2] - sx[3]
    #sy[length(sy)] = 2*sy[length(sy)-1] - sy[length(sy)-2]
    #sx[length(sx)] = 2*sx[length(sx)-1] - sx[length(sx)-2]
    sx = round(sx)
    sy = round(sy)
    
    if(p)
    {
        .trace("Display detected spaces", trace, 1)
        mask.plot = mask.inv
        mask.plot[sx, ] = 0
        mask.plot[, sy] = 0
        EBImage::display(mask.inv)
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

    #mask.test = mask.sqr
    #d = sample(0:1536, 1537)
    #for(j in 1:1536) {
    #    f = mask.test==j
    #    mask.test[f] = d[j]
    #}
    #display(mask.test)
    #display(mask.test/1536)
    
    list(data=ret, mask=mask.sqr)
}


find.radius = function(img)
{
    img.thr = EBImage::thresh(img, 20, 20, 0.005)
    mask = EBImage::opening(img.thr, EBImage::makeBrush(5, shape="disc"))
    mask = EBImage::bwlabel(mask)
    mask.ftrs = EBImage::computeFeatures.shape(mask)
    bad.objects = which(mask.ftrs[,"s.area"] < 10 | mask.ftrs[,"s.radius.min"]*1.5 < mask.ftrs[,"s.radius.max"])
    mask[mask %in% bad.objects] = 0
    median(mask.ftrs[,"s.radius.mean"], na.rm=T)
}

find.intensity = function(grid.data, img.bw, mode="center", format=max(grid.data$col) * max(grid.data$row))
{
    fdim = format.dim(format)
    if(mode == "center")
    {
        lmax.grid = expand.grid(x=seq(1/3, 2/3, length.out=11), y=seq(1/3, 2/3, length.out=11))
        lmax = apply(lmax.grid, 1, function(z, d) {
            ceiling(z[2]*d$yt + (1-z[2])*d$yb)*nrow(img.bw) + ceiling(z[1]*d$xl + (1-z[1])*d$xr)
        }, d=grid.data$data)
        lmax = matrix(img.bw[lmax], nrow=prod(fdim))
        lmax = apply(lmax, 1, quantile, 0.95)
        return(lmax)
    }
    
    if(mode == "background")
    {
        bg = with(grid.data$data, cbind(img.bw[yt*nrow(img.bw) + xl], img.bw[yt*nrow(img.bw) + xr],
            img.bw[yb*nrow(img.bw) + xl], img.bw[yb*nrow(img.bw) + xr]))
        bg = apply(bg, 1, min)
        return(bg)
    }
}

find.pins = function(img.bw, format, p=F)
{
    r=find.radius(img.bw)
    
    fdim = format.dim(format)
    ret = find.grid(img.bw, format=format, p=p)
    
    lmax = find.intensity(ret, img.bw, "center", format=format)
    bg = find.intensity(ret, img.bw, "background", format=format)

    img.bw.stack = EBImage::stackObjects(ret$mask, img.bw, bg.col=-1)
    img.bw.stack.orig = img.bw.stack
    frame.dim = prod(dim(img.bw.stack)[1:2])
    
    img.bw.stack_bg = EBImage::Image(rep(bg[1:prod(fdim)], each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack[img.bw.stack < 0] = img.bw.stack_bg[img.bw.stack < 0]
    
    img.bw.stack_max = EBImage::Image(rep(lmax[1:prod(fdim)], each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_class = EBImage::Image(rep(1:prod(fdim), each=frame.dim), dim=dim(img.bw.stack))
    img.bw.stack_hat = EBImage::selfcomplementaryTopHatGreyScale(img.bw.stack, EBImage::makeBrush(round.odd(r/2)))
    
    mask.stack.a = img.bw.stack > img.bw.stack_bg + 0.1*(img.bw.stack_max - img.bw.stack_bg)
    mask.stack.b = (img.bw.stack_max - img.bw.stack_bg) / img.bw.stack_bg > 1/4
    mask.stack = mask.stack.a & mask.stack.b
    mask.stack = EBImage::bwlabel(mask.stack)
    mask.stack = EBImage::fillHull(mask.stack)
    mask.stack = EBImage::opening(mask.stack > 0, EBImage::makeBrush(round.odd(r/2)))
    mask.stack = (mask.stack > 0) * img.bw.stack_class
    
    mask = ret$mask
    mask[T] = 0    
    mask.sqr.vector = which(ret$mask>0)[order(as.numeric(ret$mask[ret$mask>0]))]
    mask[mask.sqr.vector] = mask.stack[img.bw.stack.orig>=0]
    
    seeds = mask
    seeds[T] = 0
    seeds[ret$data$x, ret$data$y] = mask[ret$data$x, ret$data$y]
    mask.prop = propagate(img.bw, seeds, mask)
        
    if(p)
    {
        displayc(EBImage::tile(mask.stack, fg.col="#000000", nx=fdim[2]), EBImage::tile(img.bw.stack, fg.col=gray(mean(bg)), nx=fdim[2]))
        displayc(mask.prop, img.bw)
        EBImage::display(EBImage::tile(mask.stack.b, 48))
    }
    
    list(mask.pins=mask.prop, mask.sqr=ret$mask, data=ret$data)
}



parse.file = function()
{
    files = list.files("benchmark", pattern=".*\\.JPG$", full.names=T)
    for(file in files)
    {
        writeLines(paste0(which(file %in% files), " / ", length(files), ": ", file))
        #file = "benchmark/1536_4day.JPG"
        
        format = as.numeric(gsub("(\\d+)_.*", "\\1", basename(file)))
        # Read file
        img.big = EBImage::readImage(file)
        img = EBImage::resize(img.big, nrow(img.big)/3, ncol(img.big)/3)
        img.bw = EBImage::channel(img, "grey")
        img.hsv = rgb2hsv(img)
        img.cmyk = rgb2cmyk(img)
        
        ret = find.pins(img.bw, format, p=T)
        
        img.ftrs = computeFeatures.shape(ret$mask.pins)
        img.ftrs = cbind(img.ftrs, computeFeatures.moment(ret$mask.pins))
        img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img, c("red", "green", "blue")))
        img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img.cmyk, c("cyan", "magenta", "yellow", "black")))
        img.ftrs = cbind(img.ftrs, computeFeatures.basic2(ret$mask.pins, img.hsv, c("hue", "saturation", "value")))
        img.ftrs$label = as.numeric(rownames(img.ftrs))
        
        color.ftrs = colnames(img.ftrs)[grepl("red|green|blue|cyan|magenta|yellow|black", colnames(img.ftrs)) & grepl("mean", colnames(img.ftrs))]
        head(img.ftrs[,color.ftrs])
        img.pca = princomp(img.ftrs[,color.ftrs])
        writeLines(paste0("First compoment explains ", round(100*(img.pca$sdev[1]^2 / sum(img.pca$sdev^2)), 1), "% variance"))
        img.ftrs$variable = img.pca$scores[,1]
        
        displayf(ret$mask.pins, img.ftrs)
        displayc(ret$mask.pins, img)
    }
}