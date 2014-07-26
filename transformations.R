library(EBImage)
library(reshape2)
library(plyr)

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

.inverse = function(mask)
{
    mask.inv = mask
    mask.inv[T] = 0
    mask.inv[mask == 0] = 1
    mask.inv
}


blur = function(img, size)
{
    f = makeBrush(round.odd(size), shape='disc', step=FALSE)
    f = f/sum(f)
    filter2(img, f)
}

stackObjects2 = function(mask, img)
{
    img.stack = EBImage::stackObjects(mask, img, bg.col=-1)
    img.stack[img.stack < 0] = NA
    
    img.stack
}

tile2 = function(img.stack, mask.sqr)
{
    mask = mask.sqr
    mask[T] = 0    
    mask.sqr.vector = which(mask.sqr>0)[order(as.numeric(mask.sqr[mask.sqr>0]))]
    mask[mask.sqr.vector] = img.stack[!is.na(img.stack)]
    
    return(mask)
}


.corner.objects = function(mask)
{
    mask.corner = img.bw.stack
    mask.corner[T] = 1
    mask.corner[is.na(img.bw.stack)] = 0
    mask.corner = !(erode(mask.corner, kern=makeBrush(3, "box")) > 0)
}

remove.grain = function(mask, grain.size, bg.size=1, grain.fill=1e-3)
{
    f = matrix(0, nrow=grain.size+bg.size*2, ncol=grain.size+bg.size*2)
    f[1:bg.size,] = f[(nrow(f)-bg.size+1):nrow(f),] = f[,(ncol(f)-bg.size+1):ncol(f)] = f[,1:bg.size] = 1
    f = f/sum(f)
    
    grain.fill = ifelse(grain.fill > 1, grain.fill/sum(f>0), grain.fill) 
    mask.ret = mask & (filter2(mask > 0, f) >= grain.fill)
    mask.ret
}

express.objects = function(mask, grain.size, bg.size, grain.fill=1e-3)
{
    f2 = makeBrush(grain.size, "gaussian")
    f2.con = f2[1,ceiling(nrow(f2)/2)]
    f2[f2 >= f2.con] = 1
    f2[f2 < f2.con] = 0

    f3 = makeBrush(grain.size + bg.size*2, "gaussian")
    f3.con = f3[1,ceiling(nrow(f3)/2)]
    f3[f3 >= f3.con] = 1
    f3[f3 < f3.con] = 0
    f3[(1+bg.size):(nrow(f3)-bg.size),(1+bg.size):(nrow(f3)-bg.size)] = f3[(1+bg.size):(nrow(f3)-bg.size),(1+bg.size):(nrow(f3)-bg.size)] - f2
    
    
    mask.ret = filter2(mask > 0, f2/sum(f2)) * filter2(mask > 0, f3/sum(f3))
    display(tile(mask.ret, 12), "filter")
    mask.ret = mask.ret > grain.fill
    display(tile(mask.ret, 12), "bin")
    mask.ret = fillHull(mask.ret)
    display(tile(mask.ret, 12), "fill")
    mask.ret = closing(mask.ret, kern=makeBrush(round.odd(r/4), "disc"))
    display(tile(mask.ret, 12), "close")
    
    mask.ret
}
