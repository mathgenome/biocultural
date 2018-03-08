#R functions for biocultural diversity
#Data format:
#Column 1: Name of the cultural unit
#Column 2: Count associated to cultural unit
#Column 3 ... Species and their counts
#Generate data frame with conditional frequencies of species/cultures
#d is the original data frame
spIc<-function(d){
	x<-d[-c(1,2)]/rowSums(d[-c(1,2)]) #Conditional frequencies
	x<-cbind(d[c(1,2)],x) #Complete data frame
	#Get frequencies of cultural units
	x[2]<-x[2]/sum(x[2])
	x
}


#Generate a frequency matrix culturexspecies
#From spIc object

toMatrixBC<-function(s){
	my.f<-function(x)s[2]*x #Function for product with culture frequencies
	matF<-as.data.frame(apply(s[-c(1:2)],2,my.f))
	matF<-as.matrix(matF)
	dimnames(matF)<-NULL
	matF
}

#Parameter estimation

bc.info<-function(dat){
	datF<-spIc(dat)
	matF<-toMatrixBC(datF)
	#Entropy function
	MyLog2p<-function(x){if(x==0) 0 else x*log(x,2)} 
	entropy<-function(x){-sum(sapply(x,MyLog2p))}
	#Vectors needed
		fi<-rowSums(matF)  #fi.
		fj<-colSums(matF)  #f.j
	#Matrices needed
		fj.i<-datF[-c(1,2)]  #fj|i
		#fi|j
		my.f<-function(x)x/fj
		fi.j<-t(apply(matF,1,my.f))
	#Specificities
	my.s<-function(j){entropy(fi)-entropy(fi.j[,j])}
	specificities<-sapply(c(1:length(fj)),my.s)
	#Relative specificities
	r.specificities<-specificities/log(dim(datF)[1],2)
	#Specialization of cultural groups
	specializations<-as.matrix(fj.i)%*%specificities
	specializations<-as.vector(specializations)
	#Relative specializations
	r.specializations<-specializations/log(dim(datF)[1],2)
	s.diversity<-apply(datF[-c(1,2)],1,entropy)
	#Mutual information
	MutInfo<-fj%*%specificities
	MutInfo<-MutInfo[1,1]
	#Biocultural complexity
	b.c<-2^MutInfo
	#Relative B
	r.bc<-b.c/dim(datF)[1]
	#Final object
	obj<-list("cu"=as.character(dat[,1]),"species"=names(dat)[-c(1,2)],"specificities"=specificities,"r.specificities"=r.specificities,"s.diversity"=s.diversity,"specializations"=specializations,"r.specializations"=r.specializations,"mutualInfo"=MutInfo,"BC"=b.c,"r.BC"=r.bc)
	obj
	
}


