failureRegions <- function(experiments,parameterspace,fail) {
	if (class(experiments) != "poi") stop("experiments has to be an object of class poi")
	if (class(parameterspace) != "poi") stop("parameterspace has to be an object of class poi")
	fail<-as.logical(fail)
    solution <- .Call(rPorta::rConeReductionINT,getNumerator(experiments),getDenominator(experiments),getNumerator(parameterspace),getDenominator(parameterspace),as.logical(fail));	
	cones<-list()
	for (i in c(1:length(solution[[3]]))) {
		cones[[i]]<-new("ieq",num=solution[[3]][[i]],den=solution[[4]][[i]])
	}	
	conesources<-index(experiments,fail==1,)
	return(new("paramspacerestriction",parameterspace=parameterspace,feasiblepoints=solution[[1]],excludingcone=solution[[2]],cones=cones,conesources=conesources))
}

.gcd <- function(a,b) ifelse (b==0, a, .gcd(b, a %% b))

.d2q <- function(v) {
	maxValue=.Machine$integer.max
	if (is.data.frame(v)) v<-as.matrix(v)
	if (max(v)>maxValue) stop(paste("value",max(v),"exceeds the allowed maximal value"))
	precision<-min(10^floor(log10(abs(maxValue))),10^floor(log10(abs(maxValue/max(v)))))
	num<-round(v*precision)
	den<-num
	den[]<-precision
	div<-.gcd(num,den)
	return(list(num/div,den/div))
}

as.poi <- function(v) {
	if (class(v)=="poi") return(v)
	if (is.null(dim(v))) {
		nrow=1
		ncol=length(v)
	} else {
		nrow=dim(v)[1]
		ncol=dim(v)[2]
	}		
	if (is.character(v)) {
		temp<-as.integer(unlist(strsplit(v,"/")))
		temp<-list(matrix(temp[((1:(length(temp)/2))*2)-1],nrow=nrow,ncol=ncol),matrix(temp[(1:(length(temp)/2))*2],nrow=nrow,ncol=ncol))		
	} else {
		temp<-.d2q(v)
	}	
	return(new("poi",num=matrix(temp[[1]],nrow=nrow,ncol=ncol),den=matrix(temp[[2]],nrow=nrow,ncol=ncol)))		
}

as.ieq <- function(v,sign=NULL) {
	if (class(v)=="ieq") return(v)
	v<-as.poi(v)
	if (dim(v@num)[2]==1) stop("More values are needed")
	if (is.null(sign)) sign=sign=rep(1,dim(v@num)[1])
	return(new("ieq",num=v@num,den=v@den,sign=sign))		
}

as.ieqFile <- function(v,sign=NULL) {
	if (class(v)=="ieqFile") return(v)
	v<-as.ieq(v,sign)
	return(new("ieqFile",inequalities=v))		
}

as.poiFile <- function(v,hull=TRUE) {
	if (class(v)=="poiFile") return(v)
	v<-as.poi(v)
	if (hull) return(new("poiFile",convex_hull=v)) else return(new("poiFile",convex_cone=v))				
}


read.portaFile <- function(file) {
		if (file.exists(file)) {
			file2<-NULL
			options <- "-R";
			solution <- .Call(rPorta::rPortaInterface,file,options,file2);		
			return (.convertReturn(solution));				
		} else {
			stop("read.portaFile needs a valid path to a porta file")
		}
}

example.ieq <- function() {
	return(new("ieqFile",
				valid=new("poi",num=matrix(c(3,3,0,2,3),nrow=1,ncol=5)),
				lower_bounds=c(0,1,2,2,2),
				upper_bounds=c(2,2,2,5,5),
				elimination_order=c(2,0,1,0,3),
				inequalities=new("ieq",
					num=t(matrix(c(27,-28,0,57,-37,0,
								 0,0,0,-1,1,1,
								 0,1,0,0,-2,-3,
								 0,0,-1,0,0,0,
								 0,-2,0,0,1,0,
								 0,-4,0,0,-1,-1),ncol=6,nrow=6)),
					den=t(matrix(c(rep(1,30),
								 1,15,1,1,15,1),ncol=6,nrow=6)),
					sign=c(0,0,-1,-1,-1,1))))
}			

example.poi <- function() {
	return(new("poiFile",
				convex_hull=new("poi",num=t(matrix(c(3,3,0,
													 5,1,0,
													 1,5,0),nrow=3,ncol=3)),
									  den=t(matrix(c(1,1,1,
													 3,1,1,
													 1,2,1),nrow=3,ncol=3))),
				convex_cone=new("poi",num=matrix(c(0,0,2),nrow=1,ncol=3),
									  den=matrix(c(1,1,3),nrow=1,ncol=3))
				))
}	

.portaCall <- function (fileObject,options,fileObject2=NULL){
		if (is.na(max(fileObject))) stop("Not enough data to run PORTA.")
		if (max(fileObject)>.Machine$integer.max) stop(paste("PORTA does not accept integers bigger than",.Machine$integer.max))
		filename <- paste(getwd(),"/porta",sep="");
		filename2<-NULL
		writeToFile(fileObject,file=filename);
		if (!is.null(fileObject2)) {
			if (is.na(max(fileObject2))) stop("Not enough data to run PORTA.")
			if (max(fileObject2)>.Machine$integer.max) stop(paste("PORTA does not accept integers bigger than",.Machine$integer.max))		
			writeToFile(fileObject2,file=filename); 
			if (class(fileObject2) == "ieqFile"){
				filename2 <- paste(filename,".ieq",sep="");
			} else if (class(fileObject2) == "poiFile"){
				filename2 <- paste(filename,".poi",sep="");
			}
		}				
	if (class(fileObject) == "ieqFile"){
		filename <- paste(filename,".ieq",sep="");
	} else if (class(fileObject) == "poiFile"){
		filename <- paste(filename,".poi",sep="");
	}
	solution <- .Call(rPorta::rPortaInterface,filename,options,filename2);		
	if (file.exists(filename)) file.remove(filename)
	if ((!is.null(filename2))&&file.exists(filename2)) file.remove(filename2)
	return (solution)
}

.convertReturn <- function(ret) {
	if (ret[[1]]==0) {
		valid=new("poi")
		validity=matrix()		
		if (length(ret[[2]])>0) valid=new("poi",num=ret[[2]],den=ret[[3]])
		if (length(ret[[8]])>0) validity=ret[[8]]			
		fileObj<-new("ieqFile",
				valid=valid, 
				inequalities=new("ieq",
					num=rbind(ret[[4]],ret[[6]]),
					den=rbind(ret[[5]],ret[[7]]),
					sign=c(rep(0,dim(ret[[4]])[1]),rep(-1,dim(ret[[6]])[1]))),
				strong_validity=validity)
	} else if (ret[[1]]==1) {
		convex=new("poi")
		cone=new("poi")
		validity=matrix()
		if (length(ret[[2]])>0) cone=new("poi",num=ret[[2]],den=ret[[3]])
		if (length(ret[[4]])>0) convex=new("poi",num=ret[[4]],den=ret[[5]])	
		if (length(ret[[6]])>0) validity=ret[[6]]	
		fileObj<-new("poiFile",	convex_hull=convex,	convex_cone=cone,strong_validity=validity)		
	} else if (ret[[1]]==2) {
		fileObj<-ret[[2]]	
		if (length(ret[[3]])>0) {
			cat("equations :\n");
			show(new("ieq",num=ret[[3]],den=ret[[4]],sign=rep(0,dim(ret[[4]])[1])))
		}
	} else if (ret[[1]]==3){
		fileObj <- c()
		for (i in 2:length(ret)){
			ineq_convex=new("poi")
			ineq_cone=new("poi")
			ineq_validity=matrix()
			ineq_cone_valid <- FALSE
			ineq_conv_valid <- FALSE
			if (length(ret[[i]][[1]])>0){
				if (length(ret[[i]][[1]][[1]])>0){
					ineq_cone=new("poi",num=ret[[i]][[1]][[1]],den=ret[[i]][[1]][[2]])
					ineq_cone_valid <- TRUE
				}
				if (length(ret[[i]][[1]][[3]])>0){
					ineq_convex=new("poi",num=ret[[i]][[1]][[3]],den=ret[[i]][[1]][[4]])	
					ineq_conv_valid <- TRUE
				}
			}
			fO <- "no points"
			if (ineq_conv_valid || ineq_cone_valid){
				fO<-new("poiFile",	convex_hull=ineq_convex,convex_cone=ineq_cone,strong_validity=ineq_validity)
			}
			fileObj <- c(fileObj,fO)
			eq_convex=new("poi")
			eq_cone=new("poi")
			eq_validity=matrix()
			eq_cone_valid <- FALSE
			eq_conv_valid <- FALSE
			if (length(ret[[i]][[2]])>0){
				if (length(ret[[i]][[2]][[1]])>0){
					eq_cone=new("poi",num=ret[[i]][[2]][[1]],den=ret[[i]][[2]][[2]])
					eq_cone_valid <- TRUE
				}
				if (length(ret[[i]][[2]][[3]])>0){
					eq_convex=new("poi",num=ret[[i]][[2]][[3]],den=ret[[i]][[2]][[4]])	
					eq_conv_valid<-TRUE
				}
			}
			fO <- "no points"
			if (eq_conv_valid || eq_cone_valid){
				fO<-new("poiFile",	convex_hull=eq_convex,convex_cone=eq_cone,strong_validity=eq_validity)
			}
			fileObj<-c(fileObj,fO)
		}
		fileObj <- as.list(fileObj) 
	}		
	return(fileObj)
}
.gridi <- function(llim=-1, ulim=1, dis=0.5, fac=2)  
{ axe <- seq(llim,ulim,dis)
	n <- length(axe) 
	fac2 <- rep(axe,rep(n,n))
	fac1 <- rep(axe,n)
	gri <- data.frame(fac1,fac2)
	if(fac >= 3)
	{ for(i in 3:fac)
		{  fac3 <- rep(axe, rep(dim(gri)[1],n))
			gri3 <- gri
			for(j in 2:n) gri3 <- rbind(gri3, gri)  
			gri <- data.frame(gri3,fac3)
		} 
	}
	return(gri)
}

.exampleTraf <- function() {
	return(traf(example.ieq(),opt_elim=TRUE));
}	

.exampleASOP <- function() {
	candidates <- .gridi(llim=0, ulim=1, dis=0.2, fac=3)  

	design1 <- data.frame(fac1=c(0.8,0.4,0.2,0.2,1.0,0.2,0.8,0.8,0.2,0.2,0.8,0.2,0.6,0.0,0.8,0.2,0.8,0.4,0.4,0.8),
		fac2=c(0.2,0.2,0.2,0.8,0.2,0.2,0.8,0.8,0.2,0.8,0.8,0.2,0.0,0.6,0.6,0.8,1.0,0.4,0.8,0.2),
		fac3=c(0.8,0.0,1.0,0.8,0.2,0.4,0.8,0.2,0.2,0.2,0.0,0.6,0.2,0.2,0.6,0.6,0.4,0.2,1.0,1.0))

	fail <- c(0,1,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,1,0,0)
	return(failureRegions(as.poi(design1),as.poi(candidates),fail));
}