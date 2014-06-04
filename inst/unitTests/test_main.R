

library(fmcsR)


test_a.fmcs <- function(){

	data(sdfsample)
	mcs = fmcs(sdfsample[1],sdfsample[2])
	#print("class:")
	#print(class(mcs))

	s=stats(mcs)
	#print(s)

	checkEquals(mcs1(mcs)$mcs1[[1]], c(8,17,18,16,7,11,19,20,24,29,1))
	checkEquals(mcs2(mcs)$mcs2[[1]], c(4,10,12,15,8,18,9,11,14,17,1))
	checkEquals(s["Query_Size"],33,checkNames=FALSE )
	checkEquals(s["Target_Size"],26,checkNames=FALSE )
	checkEquals(s["MCS_Size"],11,checkNames=FALSE )
	checkEquals(s["Tanimoto_Coefficient"], 0.2291667,
					checkNames=FALSE,tolerance=0.00001)
	checkEquals(s["Overlap_Coefficient"], 0.4230769 ,
					checkNames=FALSE,tolerance=0.00001)


	mcs = fmcs(sdfsample[1],sdfsample[2],timeout= 20)
	#print(stats(mcs))

}

test_b.fmcsBatch <- function(){

	#DEACTIVATED("disabled")
	data(sdfsample)
	t1=system.time(results <<-
		fmcsBatch(sdfsample,sdfsample,numParallel=1))[[3]]
	checkEquals(ncol(results),5)
	checkEquals(nrow(results),length(sdfsample))
	checkEqualsNumeric(results[5,],c(33,23,14,0.3333333,0.6086957),tolerance=0.00001)
	timeout=20
	t2 = system.time(fmcsBatch(sdfsample,sdfsample,timeout=timeout))[[3]]
	message("first run time: ",t1," second: ",t2)
	checkTrue(t2 < t1)
}
