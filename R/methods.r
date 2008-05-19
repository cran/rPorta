setGeneric("convertToString",function(object) standardGeneric("convertToString"))
setGeneric("getNumerator",function(object) standardGeneric("getNumerator"))
setGeneric("getDenominator",function(object) standardGeneric("getDenominator"))
setGeneric("index",function(x, ...) standardGeneric("index"))
setGeneric("getFeasiblePoints",function(object) standardGeneric("getFeasiblePoints"))
setGeneric("getUnfeasiblePoints",function(object) standardGeneric("getUnfeasiblePoints"))
setGeneric("writeToFile",function(object,file="porta") standardGeneric("writeToFile"))
setGeneric("traf",function(object, opt_elim=FALSE,chernikov_rule_off=FALSE,validity_table_out=FALSE,long_arithmetic=FALSE) standardGeneric("traf"))
setGeneric("portsort",function(object) standardGeneric("portsort"))
setGeneric("fmel",function(object,chernikov_rule_off=FALSE,long_arithmetic=FALSE) standardGeneric("fmel"))
setGeneric("vint",function(object) standardGeneric("vint"))
setGeneric("iespo",function(object,poiObject,validity_table_out=FALSE) standardGeneric("iespo"))
setGeneric("posie",function(object,poiObject) standardGeneric("posie"))
setGeneric("fctp",function(object,poiObject) standardGeneric("fctp"))
setGeneric("getParamspaceInfo",function(object) standardGeneric("getParamspaceInfo"))


setMethod("writeToFile","portaFile",
	function(object,file="porta") {
		if (class(object)=="portaFile") stop("Class portaFile is an abstract class.")
		suffix<-"ieq"
		if (class(object)=="poiFile") suffix<-"poi"
		targetString <- paste(file,suffix,sep=".")
		textfile=file(targetString,open="w")
		cat(convertToString(object),file=textfile)
		close(textfile)		
	}	
)

setMethod("show","portaFile",
	function(object) {
		if (class(object)=="portaFile") stop("Class portaFile is an abstract class.")
		cat(convertToString(object))
		if ( (!any(is.na(object@strong_validity))) && (length(object@strong_validity)>0)  ) {
			cat("\nstrong validity table (rows=points, columns=inequalities):\n")
			show(object@strong_validity)
		}				
	}	
)


setMethod("convertToString","poi",
	function(object) {
		if (length(object@den)==0) {
			conversion=.d2q(object@num)
			object@num=conversion[[1]]
			object@den=conversion[[2]]
		}	
		str<-""
		if (any(dim(as.matrix(object@num))!=dim(as.matrix(object@den)))) stop("numerators and denominators have different sizes")
		for (i in c(1:dim(as.matrix(object@num))[1])) {
			for (j in c(1:dim(as.matrix(object@num))[2])) {
				num=as.matrix(object@num)[i,j]
				den=as.matrix(object@den)[i,j]
				if (den<0) {
						num<--num
						den<--den
				}		
				if ((den==1)||(num==0)) str<-paste(str,num," ",sep="") 	
				else str<-paste(str,num,"/",den," ",sep="")
			}
			str<-paste(str,"\n")	 	
		}	
		return(str)			
	}
)

setMethod("getNumerator","poi",
	function(object) {
		if (length(object@den)==0) {
			conversion=.d2q(object@num)
			object@num=conversion[[1]]
		}		
		return(object@num)
	}
)

setMethod("getDenominator","poi",
	function(object) {
		if (length(object@den)==0) {
			conversion=.d2q(object@num)
			object@den=conversion[[2]]
		}		
		return(object@den)
	}
)

setMethod("index","poi",
	function(x, ...) {
		num=x@num[...]
		if (!is.matrix(num)) num<-as.matrix(t(num))
		if (length(x@den)==0) {
			return(new("poi",num=num,den=x@den))
		} else {	
			den=x@den[...]
			if (!is.matrix(den)) den<-as.matrix(t(den))
			return(new("poi",num=num,den=den))
		}	
	}
)		

setMethod("show","poi",
	function(object) {
		cat(convertToString(object))
	}	
)

setMethod("as.matrix","poi",
	function(x) {
	  return(x@num/x@den)
	}	
)

setMethod("as.character","poi",
	function(x) {
	  paste(x@num,x@den,sep="/")
	}	
)

setMethod("max","poi",
	function(x,na.rm=TRUE) {
	  if (all(is.na(x@num))) numMax=NA else numMax=max(abs(x@num),na.rm=TRUE)
  	  if (all(is.na(x@den))) denMax=NA else denMax=max(abs(x@den),na.rm=TRUE)
  	  if (all(is.na(c(numMax,denMax)))) poiMax=NA else poiMax=max(numMax,denMax,na.rm=TRUE)
	  return(poiMax)
	}	
)


setMethod("convertToString","poiFile",
	function(object) {
		str<-""
		dimension<-0
		hull<-FALSE
		cone<-FALSE
		if (length(object@convex_hull@num)>0) {		
			hull<-TRUE
			dimension<-dim(as.matrix(object@convex_hull@num))[2]
		}
		if (length(object@convex_cone@num)>0) {		
			cone<-TRUE
			dimension2<-dim(as.matrix(object@convex_cone@num))[2]
			if ((dimension>0)&&(dimension!=dimension2)) stop ("convex hull and convex cone have different dimensions")
			dimension<-dimension2
		}	
		if (dimension==0) {
			cat("no points")
		} else {		
			str<-paste(str,"DIM =",dimension,"\n \n")
			if (hull) {		
				str<-paste(str,"CONV_SECTION\n",convertToString(object@convex_hull),"\n")
			}	
			if (cone) {		
				str<-paste(str,"CONE_SECTION\n",convertToString(object@convex_cone),"\n")
			}	
			str<-paste(str,"END\n")
		}
		return(str)	
	}
)

setMethod("max","poiFile",
	function(x,na.rm=TRUE) {
  	  if (all(is.na(c(max(x@convex_hull),max(x@convex_cone))))) poiMax=NA else poiMax=max(max(x@convex_hull),max(x@convex_cone),na.rm=TRUE)
	  return(poiMax)
	}	
)

setMethod("convertToString","ieq",
		function(object) {
			if (length(object@den)==0) {
				conversion=.d2q(object@num)
				object@num=conversion[[1]]
				object@den=conversion[[2]]
			}
			str<-""
			if (any(dim(as.matrix(object@num))!=dim(as.matrix(object@den)))) stop("numerators and denominators have different sizes")
			for (i in c(1:dim(as.matrix(object@num))[1])) {
				sign<-" <= "
				if (length(object@sign)>0) {
					if (length(object@sign)!=dim(as.matrix(object@num))[1]) stop("length of sign differs from the number of rows")
					if (object@sign[i]==0) sign<-" == "
					if (object@sign[i]==1) sign<-" >= "				
				}
				str<-paste(str,"(",i,") ",sep="")		
				for (j in c(1:(dim(as.matrix(object@num))[2]-2))) {
					num=as.matrix(object@num)[i,j]
					den=as.matrix(object@den)[i,j]
					if (j>1) {
						num=abs(num)
						den=abs(den)
					} else {
						if (den<0) {
							num<--num
							den<--den
						}
					}	
					op=" + "
					if ((as.matrix(object@den)[i,j+1]*as.matrix(object@num)[i,j+1])<0) op=" - "
					if ((den==1)||(num==0)) value=num else value=paste(num,den,sep="/")
					str<-paste(str,value,"x",j,op,sep="")		     	
				}
				num=abs(as.matrix(object@num)[i,j+1])
				den=abs(as.matrix(object@den)[i,j+1])
				if ((den==1)||(num==0)) value=num else value=paste(num,den,sep="/")
				str<-paste(str,value,"x",j+1,sign,sep="");
				num=as.matrix(object@num)[i,j+2]
				den=as.matrix(object@den)[i,j+2]
				if (den<0) {
					num<--num
					den<--den
				}    		
				if ((den==1)||(num==0)) value=num else value=paste(num,den,sep="/")
				str<-paste(str,value, "\n",sep="");
			}	
			return(str)			
		}
		)

setMethod("show","ieq",
	function(object) {
		cat(convertToString(object))
	}	
)

setMethod("convertToString","ieqFile",
	function(object) {
		str<-""
		dimension<-0
		inequalities<-FALSE
		if (length(object@inequalities@num)>0) {		
			inequalities<-TRUE
			dimension<-dim(as.matrix(object@inequalities@num))[2]-1
		}
		if (dimension==0) {
			cat("no inequalities")
		} else {		
			str<-paste(str,"DIM =",dimension,"\n \n")
			if (length(object@valid@num)>0) {
				if(length(object@valid@num)!=dimension) stop("length of vector valid does not match problem dimension")
				str<-paste(str,"VALID\n",convertToString(object@valid),"\n")
			}
			if ((length(object@lower_bounds)>0)||(length(object@upper_bounds)>0)) {
				if(length(object@lower_bounds)!=dimension) stop("length of vector lower_bounds does not match problem dimension")
				if(length(object@upper_bounds)!=dimension) stop("length of vector upper_bounds does not match problem dimension")
				str<-paste(str,"LOWER_BOUNDS\n",gsub(",","",toString(object@lower_bounds)),"\n")
				str<-paste(str,"UPPER_BOUNDS\n",gsub(",","",toString(object@upper_bounds)),"\n \n")
	 		}		
			if (length(object@elimination_order)>0) {
				if(length(object@elimination_order)!=dimension) stop("length of vector elimination_order does not match problem dimension")
				str<-paste(str,"ELIMINATION_ORDER\n",gsub(",","",toString(object@elimination_order)),"\n \n")
			} 		
			if (inequalities) {		
				str<-paste(str,"INEQUALITIES_SECTION\n",convertToString(object@inequalities),"\n")
			}			
			str<-paste(str,"END\n")
		}	
		return(str)	
	}
)

setMethod("max","ieqFile",
	function(x,na.rm=TRUE) {
  	  if (all(is.na(c(max(x@valid),max(x@inequalities))))) ieqMax=NA else ieqMax=max(max(x@valid),max(x@inequalities),na.rm=TRUE)
	  return(ieqMax)
	}	
)

setMethod("show","paramspacerestriction",
	function(object) {
		denominator<-FALSE
		if (length(object@parameterspace@den)>0) denominator<-TRUE
		cat("Feasible Points:\n")
		if (denominator) {
			submatrix<-index(object@parameterspace,object@feasiblepoints,)
			feasible<-matrix(as.character(submatrix),nrow=dim(submatrix@num)[1],ncol=dim(submatrix@num)[2],dimnames=dimnames(submatrix@num))
			feasible<-data.frame(feasible)
			show(feasible)				
		} else {
			show(object@parameterspace@num[object@feasiblepoints,])	
		}
		cat("\nUnfeasible Points:\n")
		if (denominator) {
			submatrix<-index(object@parameterspace,-object@feasiblepoints,)
			unfeasible<-matrix(as.character(submatrix),nrow=dim(submatrix@num)[1],ncol=dim(submatrix@num)[2],dimnames=dimnames(submatrix@num))
			unfeasible<-cbind(unfeasible,object@excludingcone)
			dimnames(unfeasible)[[2]][length(dimnames(unfeasible)[[2]])]<-"PCC"
			unfeasible<-data.frame(unfeasible)			
			show(unfeasible)				
		} else {
			unfeasible<-cbind(object@parameterspace@num[-object@feasiblepoints,],object@excludingcone)
			dimnames(unfeasible)[[2]][length(dimnames(unfeasible)[[2]])]<-"PCC"
			show(unfeasible)	
		}
		cat("\nCones:\n\n")
		for (i in c(1:length(object@cones))) {
			cat("PCC for point ")
			print(index(object@conesources,i,))
			show(object@cones[[i]])
			cat("\n")
		}	
	}
)

setMethod("getUnfeasiblePoints","paramspacerestriction",
	function(object) {
		if (length(object@parameterspace@den)>0) {
			return(index(object@parameterspace,-object@feasiblepoints,))		
		} else {
			return(object@parameterspace@num[-object@feasiblepoints,])
		}		
	}
)

setMethod("getFeasiblePoints","paramspacerestriction",
	function(object) {
		if (length(object@parameterspace@den)>0) {
			return(index(object@parameterspace,object@feasiblepoints,))		
		} else {
			return(object@parameterspace@num[object@feasiblepoints,])
		}		
	}
)

setMethod("getParamspaceInfo","paramspacerestriction",
	function(object) {
		paramspace<-matrix(as.character(object@parameterspace),nrow=dim(object@parameterspace@num)[1],ncol=dim(object@parameterspace@num)[2],dimnames=dimnames(object@parameterspace@num))
		feasible<-rep(0,dim(object@parameterspace@num)[1])
		feasible[-object@feasiblepoints]=object@excludingcone
		sources<-vector("character",length(object@cones))
		for (i in c(1:length(object@cones))) {
			sources[i]<-convertToString(index(object@conesources,i,))
			sources[i]<-substr(sources[i],1,nchar(sources[i])-3)
		}
		sources<-c(sources,"-")
		temp<-feasible
		temp[temp==0]<-length(object@cones)+1
		sources<-sources[temp] 
		return(data.frame(paramspace,PCC=feasible,PCC_source=sources))
	}
)

setMethod("traf","portaFile",
	function(object, opt_elim=FALSE,chernikov_rule_off=FALSE,validity_table_out=FALSE,long_arithmetic=FALSE){
		options <- "-T";
		if (opt_elim) options<-paste(options,"o",sep="");
		if (chernikov_rule_off) options<-paste(options,"c",sep="");
		if (validity_table_out) options<-paste(options,"v",sep="");
		if (long_arithmetic) options<-paste(options,"l",sep="");		
		return (.convertReturn(.portaCall(object,options)));
	}
)

setMethod("portsort","portaFile",
	function(object){
		options <- "-S";
		return (.convertReturn(.portaCall(object,options)));
	}
)	

setMethod("dim","poiFile",
		function(x) {
			options <- "-D";
			return (.convertReturn(.portaCall(x,options)));			
		}
)

setMethod("fmel","ieqFile",
	function(object,chernikov_rule_off=FALSE,long_arithmetic=FALSE){
		options <- "-F";
		if (chernikov_rule_off) options<-paste(options,"c",sep="");		
		if (long_arithmetic) options<-paste(options,"l",sep="");				
		return (.convertReturn(.portaCall(object,options)));
	}
)	

setMethod("vint","ieqFile",
	function(object){
		options <- "-V";
		return (.convertReturn(.portaCall(object,options,object)));
	}
)	

setMethod("iespo","ieqFile",
	function(object,poiObject,validity_table_out=FALSE){
	 	if (!(class(poiObject) == "poiFile")) stop("needs an object of class poiFile")
		options <- "-I";
		if (validity_table_out) options<-paste(options,"v",sep="");
		return (.convertReturn(.portaCall(object,options,poiObject)));
	}
)	

setMethod("posie","ieqFile",
	function(object,poiObject){
	 	if (!(class(poiObject) == "poiFile")) stop("needs an object of class poiFile")
		options <- "-P";
		return (.convertReturn(.portaCall(object,options,poiObject)));
	}
)	

setMethod("fctp","ieqFile",
	function(object,poiObject){
	 	if (!(class(poiObject) == "poiFile")) stop("needs an object of class poiFile")
		options <- "-C";
		return (.convertReturn(.portaCall(object,options,poiObject)));
	}
)	