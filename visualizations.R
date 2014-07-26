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
