round.odd = function (x) 
{
    x = floor(x)
    if (x%%2 == 0) 
        x = x + 1
    x
}


.trace = function(text, level, req)
{
    if(level >= req)
        writeLines(paste0(format(Sys.time(), "%X"), ": ", text))
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