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

#######
#Bootstrap functions

library(reshape2)
dat.long<-function(dat){
	melt(dat,id.vars=c(1,2),value.name="count")
} #Convert to long format

#datL is data in long format

get.row<-function(datL,n){
	csum<-cumsum(datL$count)
	c(1:dim(datL)[1])[n<=csum][1]
} #Get the row in the long format to which the n-th specimen belongs

boot.sample<-function(datL){
	sam<-sample(1:sum(datL$count),replace=T)
	sapply(sam,function(x)get.row(datL,x))
} #Get the row for a random sample of specimens with replacement

library(tidyr)

get.new.data<-function(datL){
	x<-boot.sample(datL)
	frame<-as.data.frame(table(x))
	ndat<-datL
	ndat$count[as.numeric(as.vector(frame$x))]<-frame$Freq
	spread(ndat,variable,count)
} #Get a new data frame with random sampling with replacement

boot.data<-function(data){
	datL<-dat.long(data)
	get.new.data(datL)
} #Get a new data frame with random sampling with replacement from data in wide format

#Bootstrap n times

boot.n<-function(data,n){
	f<-vector("list",n) #Specificities
	fr<-vector("list",n) #Relative Specificities
	a<-vector("list",n) #Specialization
	ar<-vector("list",n) #Relative specialization
	s<-vector("list",n) #Shannon diversity
	b<-NA #Biocultural complexity
	br<-NA #Relative biocultural complexity
	length(b)<-n
	length(br)<-n
	for(i in 1:n){
	x<-bc.info(boot.data(data))
	f[[i]]<-x$specificities
	fr[[i]]<-x$r.specificities
	a[[i]]<-x$specializations
	ar[[i]]<-x$r.specializations
	s[[i]]<-x$s.diversity
	b[i]<-x$BC
	br[i]<-x$r.BC}
	list(f,fr,a,ar,s,b,br)
} #Object list with: specificities, r specificities, specializations, r specializations, Shannon diversity, BC and r BC from n bootstrap resamplings.

#Get standard errors of estimates
#From object created with boot.n

estimate.se<-function(bn){
	f<-lapply(bn[[1]],as.data.frame)
	f<-as.data.frame(f)
	sf<-apply(f,1,sd) #Standard errors for specificities
	mf<-apply(f,1,mean) #Bootstrap means for specificities
	fr<-lapply(bn[[2]],as.data.frame)
	fr<-as.data.frame(fr)
	sfr<-apply(fr,1,sd) #Standard errors for r specificities
	mfr<-apply(fr,1,mean) #Bootstrap means for r specificities
	a<-lapply(bn[[3]],as.data.frame)
	a<-as.data.frame(a)
	sa<-apply(a,1,sd) #Standard errors for specializations
	ma<-apply(a,1,mean) #Bootstrap means for specializations
	ar<-lapply(bn[[4]],as.data.frame)
	ar<-as.data.frame(ar)
	sar<-apply(ar,1,sd) #Standard errors for r specializations
	mar<-apply(ar,1,mean) #Bootstrap means for r specializations
	s<-lapply(bn[[5]],as.data.frame)
	s<-as.data.frame(s)
	ss<-apply(s,1,sd) #Standard errors for Shannon entropy
	ms<-apply(s,1,mean) #Bootstrap means for Shannon entropy
	sb<-sd(bn[[6]])
	sbr<-sd(bn[[7]])
	clb<-quantile(bn[[6]],probs=c(.025,.975))
	clbr<-quantile(bn[[7]],probs=c(.025,.975))
	meanb<-mean(bn[[6]])
	meanbr<-mean(bn[[7]])
list("m_specif"=mf,"sd_specif"=sf,"m_rspecif"=mfr,"sd_rspecif"=sfr,"m_specia"=ma,"sd_specia"=sa,"m_rspecia"=mar,"sd_rspecia"=sar,"m_sdiv"=ms,"sd_sdiv"=ss,"m_BC"=meanb,"sd_BC"=sb,"m_rBC"=meanbr,"sd_rBC"=sbr, "cl_BC"=clb,"cl_rBC"=clbr)
} #Object with means and standard errors of sspecificities, r specificities, specializations, r specializations, Shannon diversity, BC and r BC. The last ones are the 95%CI for BC and rBC.

########
#Generate complete tables
#dat is the data table, and n is the number of bootstrap resamplings
#Bias correction is made as follows:
#bias = bootstrap mean - estimate
#bcorr.estimate=estimate-bias
#bcorr.estimate=estimate-(bootstrap mean - estimate)
#bcorr.estimate=2*estimate-bootstrap mean
#See page 44 Manly book.. Randomization, Bootstrap...
bc.tables<-function(data,n){
	mybc<-bc.info(data)
	mybn<-boot.n(data,n)
	myse<-estimate.se(mybn)
	tabSpecies<-data.frame("Species"=mybc$species,"Specificity"=mybc$specificities,"BCorSpec"=2*mybc$specificities-myse$m_specif,"SE.S"=myse$sd_specif, "R.Specificity"=mybc$r.specificities,"BCorRSpec"=2*mybc$r.specificities-myse$m_rspecif,"SE.RS"=myse$sd_rspecif)
	tabCultures<-data.frame("Culture"=mybc$cu,"Specialization"=mybc$specializations,"BCorSpecia"=2*mybc$specializations-myse$m_specia,"SE.S"=myse$sd_specia,"R.Specialization"=mybc$r.specializations,"BCorRSpecia"=2*mybc$r.specializations-myse$m_rspecia,"SE.RS"=myse$sd_rspecia,"S.Diversity"=mybc$s.diversity,"BCorSDiver"=2*mybc$s.diversity-myse$m_sdiv,"SE.SD"=myse$sd_sdiv)
	tabBC<-data.frame("BC"=mybc$BC,"BCorBC"=2*mybc$BC-myse$m_BC,"SE.BC"=myse$sd_BC,"R.BC"=mybc$r.BC,"BCorRBC"=2*mybc$r.BC-myse$m_rBC,"SE.RBC"=myse$sd_rBC)
	#First type of percentile ci (p 50 Manly)
	taBC_CI_1<-data.frame("BC_CL"=myse$cl_BC[1],"BC_CH"=myse$cl_BC[2],"RBC_CL"=myse$cl_rBC[1],"RBC_CH"=myse$cl_rBC[2])
	#Second type of percentile ci (p 50 Manly). Better for skewed distributions
	taBC_CI_2<-data.frame("BC_CL"=2*mybc$BC-myse$cl_BC[2],"BC_CH"=2*mybc$BC-myse$cl_BC[1],"RBC_CL"=2*mybc$r.BC-myse$cl_rBC[2],"RBC_CH"=2*mybc$r.BC-myse$cl_rBC[1])
	row.names(taBC_CI_1)<-c("1")
	row.names(taBC_CI_2)<-c("1")
list(
tabSpecies=tabSpecies,tabCultures=tabCultures,tabBC=tabBC,
taBC_CI_1=taBC_CI_1,taBC_CI_2=taBC_CI_2
)
}
