########################
## fmcsBatch Function ##
########################
.packageName <- 'fmcsR'

fmcsBatch <-
function(querySdf, sdfset, al = 0, au = 0, bl = 0, bu = 0, 
			matching.mode = "static",timeout=0) {

    if(class(querySdf)=="SDFset") querySdf <- querySdf[[1]]
    al = as.integer(al) 
    au = as.integer(au) 
    bl = as.integer(bl) 
    bu = as.integer(bu) 
    
    if (matching.mode == "static") {
        matching.int = 0
    } else if (matching.mode == "aromatic") {
        matching.int = 1
    } else if (matching.mode == "ring") {
        matching.int = 2
    } else {
        stop("matching.mode needs to be static, aromatic, or ring.")
    }
    
    matching.int = as.integer(matching.int)
    
    s1 <- as(querySdf, "character") 
    s1 <- paste(s1, collapse="\n")
    
    targetSize <- numeric(length(sdfset))
    mcsSize <- numeric(length(sdfset))
    
    for (i in 1:length(sdfset)) {
        # print(i)
        s2 <- as(sdfset[[i]], "character")
        s2 <- paste(s2, collapse="\n")
        
        result_data <- 
        .C('fmcs_R_wrap', s1, s2, al, au, bl, bu, 
            matching.int, as.integer(0), as.integer(timeout), 
            sdfOne="", sdfTwo="",
            sdf1Size = "", sdf2Size = "", mcsSize = "", PACKAGE = 'fmcsR')
            
        querySize <- as.integer(result_data$sdf1Size)
        targetSize[i] <- as.integer(result_data$sdf2Size)
        mcsSize[i] <- as.integer(result_data$mcsSize)
        
    }
    minsize <- ifelse(querySize < targetSize, querySize, targetSize)
    tanimoto <- mcsSize/(querySize + targetSize - mcsSize)
    overlap <- mcsSize/minsize
    searchMA <- cbind(Query_Size = querySize, Target_Size = targetSize, MCS_Size = mcsSize, Tanimoto_Coefficient = tanimoto, Overlap_Coefficient = overlap) 
    rownames(searchMA) <- cid(sdfset)
    return(searchMA)
}

