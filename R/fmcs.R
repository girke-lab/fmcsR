###################
## fmcs Function ##
###################
.packageName <- 'fmcsR'

fmcs <-
function(sdf1, sdf2, al = 0, au = 0, bl = 0, bu = 0, matching.mode = "static", fast = FALSE) {
    if(class(sdf1)=="SDF") sdf1_name <- "CMP1" 
    if(class(sdf1)=="SDFset") {
        sdf1_name <- cid(sdf1) 
	sdf1 <- sdf1[[1]]
    }
    if(class(sdf2)=="SDF") sdf2_name <- "CMP2" 
    if(class(sdf2)=="SDFset") {
        sdf2_name <- cid(sdf2) 
	sdf2 <- sdf2[[1]]
    }
    
    s1 = as(sdf1, "character") 
    s1 = paste(s1, collapse="\n")
    s2 = as(sdf2, "character")
    s2 = paste(s2, collapse="\n")
    
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
    
    running.mode = 1
    if (fast) {
        running.mode = 0
    }
 
    matching.int = as.integer(matching.int)
    running.mode = as.integer(running.mode)
 
    result_data = 
        .C('fmcs_R_wrap', s1, s2, al, au, bl, bu, 
            matching.int, running.mode, 0, 
            idxOne="", idxTwo="",
            sdf1Size = "", sdf2Size = "", mcsSize = "", PACKAGE = 'fmcsR')
            
    querySize = as.integer(result_data$sdf1Size)
    targetSize = as.integer(result_data$sdf2Size)
    mcs = as.integer(result_data$mcsSize)
   
    stats <- c(Query_Size=querySize, Target_Size=targetSize, MCS_Size=mcs, Tanimoto_Coefficient=(mcs/(querySize + targetSize - mcs)), Overlap_Coefficient=(mcs/min(c(querySize, targetSize)))) 
    if (fast) {
        return(stats)
    }

    if (result_data$idxOne != "") {
    
        idx1 <- strsplit(result_data$idxOne, "\n")[[1]]
        idx2 <- strsplit(result_data$idxTwo, "\n")[[1]]
        
        sdfList1 <- list()
        sdfList2 <- list()
        
        for (idx in seq(along=idx1)) {
            idxVector <- strsplit(idx1[idx], " ")[[1]]
            idxVector <- as.integer(idxVector)
            idxVector <- as.vector(idxVector)
            sdfList1[[idx]] <- idxVector
        }
        listNames <- c()
        for (i in 1:length(sdfList1)) {
            listNames <- c(listNames, paste(sdf1_name, "_fmcs_", sprintf("%d", i), sep=""))
        }
        names(sdfList1) <- listNames

        for (idx in seq(along=idx2)) {
            idxVector <- strsplit(idx2[idx], " ")[[1]]
            idxVector <- as.integer(idxVector)
            idxVector <- as.vector(idxVector)
            sdfList2[[idx]] <- idxVector
        }
        listNames <- c()
        for (i in 1:length(sdfList2)) {
            listNames <- c(listNames, paste(sdf2_name, "_fmcs_", sprintf("%d", i), sep=""))
        }
        names(sdfList2) <- listNames
	sdf1 <- as(sdf1, "SDFset"); cid(sdf1) <- sdf1_name
	sdf2 <- as(sdf2, "SDFset"); cid(sdf2) <- sdf2_name
	return(new("MCS", stats=stats, mcs1=list(query=sdf1, mcs1=sdfList1), mcs2=list(target=sdf2, mcs2=sdfList2)))
            
    }
    
}

##########################################
## Class and Method Definitions for MCS ##
##########################################
## Define MCS class
setClass("MCS", representation(stats="numeric", mcs1="list", mcs2="list"))

## Methods to return components of MCS
setGeneric(name="stats", def=function(x) standardGeneric("stats"))
setMethod(f="stats", signature="MCS", definition=function(x) {return(x@stats)}) 
setGeneric(name="mcs1", def=function(x) standardGeneric("mcs1"))
setMethod(f="mcs1", signature="MCS", definition=function(x) {return(x@mcs1)}) 
setGeneric(name="mcs2", def=function(x) standardGeneric("mcs2"))
setMethod(f="mcs2", signature="MCS", definition=function(x) {return(x@mcs2)}) 

## Behavior of "[[" operator for MCS to return components wiht mcs[["stats"]], mcs[["mcs1"]], etc.
setMethod(f="[[", signature="MCS", definition=function(x, i, ..., drop) {
	if(i=="stats") return(x@stats)                 
	if(i=="mcs1") return(x@mcs1)                 
	if(i=="mcs2") return(x@mcs2)                 
})

## Constructor method
## List to MCS with: as(mylist, "MCS")
setAs(from="list", to="MCS", 
	def=function(from) {
		new("MCS", stats=from$stats, mcs1=from$mcs1, mcs2=from$mcs2)
})

## Define print behavior for MCS
setMethod(f="show", signature="MCS",                
   definition=function(object) {
         cat("An instance of ", "\"", class(object), "\" ", "\n", sep="")
         cat(c("", "Number of MCSs:", length(object@mcs1[[2]]), "\n", 
                   paste(cid(object@mcs1[[1]]), ":", sep=""), object@stats[1], "atoms", "\n", 
                   paste(cid(object@mcs2[[1]]), ":", sep=""), object@stats[2], "atoms", "\n", 
                   "MCS:", object@stats[3], "atoms", "\n", 
                   "Tanimoto Coefficient:", round(object@stats[4], 5), "\n",
                   "Overlap Coefficient:", round(object@stats[5], 5), sep="", "\n")
	 )
})

#################################
## Return MCS Object as SDFset ##
#################################
## Helper function to use atomsubset() on MCS objects in order to 
## obtain SDFset objects for their results.
mcs2sdfset <- function(x, ...) {
	sdfmcs1 <- SDFset() 
	sdfmcs2 <- SDFset()
	for(i in seq(along=mcs1(x)[[2]])) {
		 tmpsdfmcs1 <- atomsubset(mcs1(x)[[1]][[1]], atomrows=mcs1(x)[[2]][[i]], ...)	
		 sdfmcs1 <- suppressWarnings(c(sdfmcs1, tmpsdfmcs1))
		 tmpsdfmcs2 <- atomsubset(mcs2(x)[[1]][[1]], atomrows=mcs2(x)[[2]][[i]], ...)	
		 sdfmcs2 <- suppressWarnings(c(sdfmcs2, tmpsdfmcs2))
	}
	cid(sdfmcs1) <- names(mcs1(x)[[2]])	
	cid(sdfmcs2) <- names(mcs2(x)[[2]])	
	return(list(query=sdfmcs1, target=sdfmcs2))
}
## Usage:
# mcs2sdfset(x=mcs, type="new")
# mcs2sdfset(x=mcs, type="old")[[1]][[1]]
# plot(mcs2sdfset(x=mcs, type="new")[[1]][1]) # Plotting works only for type="new" 

#######################
## Plot MCS Function ##
#######################
plotMCS <- function(x, mcs=1, print=FALSE, ...) {
	if(class(x)!="MCS") stop("Input needs to be of MCS class")
	par(mfrow=c(length(mcs), 2))
	for(i in mcs) {
		plot(mcs1(x)[[1]][[1]], colbonds=mcs1(x)[[2]][[i]], main=cid(mcs1(x)[[1]]), ...)
		plot(mcs2(x)[[1]][[1]], colbonds=mcs2(x)[[2]][[i]], main=cid(mcs2(x)[[1]]), ...)
	}
}



