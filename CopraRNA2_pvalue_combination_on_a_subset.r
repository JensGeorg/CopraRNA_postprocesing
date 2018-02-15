
mnum=40000		
ooi<-"NC_000913"
mRNA_only<-TRUE

da<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv") ## edit prw changed file name 
options <- read.table("CopraRNA_option_file.txt", sep=":") 
root<-as.numeric(as.character(options[14,2]))



en<-grep("Annotation", colnames(da))
da2<-da[,3:(en-1)]
namegenomes<-colnames(da2)
genomes<-list()

myfun<-function(da2){
	int_table<-gsub("^([^\\|]*\\|[^\\|]*\\|)", "",da2)
	int_table<-gsub("\\|.*","",int_table)
	int_table
}


int_table<-apply(da2, 2, myfun)


int_table<-matrix(as.numeric(int_table),nrow(int_table),ncol(int_table))
colnames(int_table)<-colnames(da2)
rownames(int_table)<-seq(1,nrow(int_table))

weight<-read.csv("zscore.weight", header=F, sep=";")

weight<-weight[ match(colnames(da2),toupper(weight[,1])),]


if(root==0){
	weight[,2]<-1
}

if(root>0){
	weight[,2]<-weight[,2]^(1/root)
}


qtrans<-function(dat2){
		if(is(dat2)[1]=="numeric"){
			dat2<-matrix(dat2,1,length(dat2))
		}
        for(i in 1:ncol(dat2)){
                dat2[,i]<-qnorm(dat2[,i], lower.tail=FALSE)
        }
        dat2
}


if(ooi_close_analysis==TRUE){

mnum<-40000
if(chooseclose==TRUE){
	close_name_vector<-c()
}

if(auto_close==TRUE){
	command<-paste("clustalo -i ", "input_sRNA.fa", " --distmat-out=distmatout.txt --full --percent-id --output-order=input-order --force --max-hmm-iterations=-1", sep="")
	system(command)
	temp<-read.delim("distmatout.txt",sep="",header=F, , skip=1)
	unlink("distmatout.txt")
	na<-temp[,1]
	temp<-temp[,2:ncol(temp)]
	colnames(temp)<-na
	rownames(temp)<-na
	ooi_pos<-grep(ooi, colnames(temp))
	ooi_col<-temp[,ooi_pos]
	names(ooi_col)<-na
	ooi_col<-sort(ooi_col,decreasing = TRUE)
	num<-floor(length(ooi_col)*close)
	ooi_col<-ooi_col[1:4]
}

load("conservation_table.Rdata")



con_table<-conservation_table[[3]]
con_table_sub<-conservation_table[[4]]

p_table<-con_table
p_table[]<-NA

mRNA_test<-function(x){
	out<-NA
	if(is.na(x[1])==F){
		if(length(x)==4){
			if(x[2]=="TRUE"){
				out<-as.numeric(x[1])
			}
		}
	}
	out
}
mRNA_and_sRNA_test<-function(x){
	out<-NA
	if(is.na(x[1])==F){
	if(length(x)==4){
		if(x[2]=="TRUE"){
			if(x[3]=="TRUE" | x[4]=="TRUE"){
				out<-as.numeric(x[1])
			}
		}
	}
	}
	out
}


for(i in 1:nrow(p_table)){
	
	if(mRNA_only==TRUE){
		con_temp<-con_table[i,]
		con_temp<-strsplit(con_temp, "\\|")
		con_temp_sub<-con_table_sub[i,]
		con_temp_sub<-strsplit(con_temp_sub, "\\|")
		con_temp<-unlist(lapply(con_temp,mRNA_test))
		con_temp_sub<-unlist(lapply(con_temp_sub,mRNA_test))
		na<-which(is.na(con_temp))
		if(length(na)>0){
			con_temp[na]<-con_temp_sub[na]
		}
		
		p_table[i,]<-con_temp
	}
	
	if(mRNA_only==FALSE){
		con_temp<-con_table[i,]
		con_temp<-strsplit(con_temp, "\\|")
		con_temp_sub<-con_table_sub[i,]
		con_temp_sub<-strsplit(con_temp_sub, "\\|")
		con_temp<-unlist(lapply(con_temp,mRNA_and_sRNA_test))
		con_temp_sub<-unlist(lapply(con_temp_sub,mRNA_and_sRNA_test))
		na<-which(is.na(con_temp))
		if(length(na)>0){
			con_temp[na]<-con_temp_sub[na]
		}
		
		p_table[i,]<-con_temp
	}
}
nan<-colnames(p_table)
p_table<-matrix(as.numeric(p_table),nrow(p_table),ncol(p_table))
colnames(p_table)<-nan


close_pos<-match(names(ooi_col),colnames(p_table))
p_table<-p_table[,close_pos]

posvec<-rep(NA,10^7)
vari<-rep(NA,10^7)
count<-0
for(i in 1:nrow(p_table)){
	temp<-order(p_table[i,], na.last=NA)
	
	if(length(temp)>2){              
		for(j in 1:(length(temp)-2)){
			count<-count+1
			vari[count]<-paste(sort(temp[1:(j+2)]), collapse="_")
			posvec[count]<-i
		}
	}
}

vari<-na.omit(vari)
posvec<-na.omit(posvec)
dups<-which(duplicated(vari))
first_occurence<-match(vari[dups], vari)
posvec2<-posvec

for(i in 1:length(dups)){
	posvec2[first_occurence[i]]<-paste(posvec2[first_occurence[i]],posvec[dups[i]], sep="_")

}
posvec2<-posvec2[-dups]
vari<-unique(vari)


spa<-floor(length(vari)/mnum)
print(length(vari))
rest<-length(vari)-(mnum*spa)
spal<-rep(mnum,spa)
if(rest>0){
spa<-spa+1
spal<-c(spal,rest)
}

ooilist<-rep((list(rep(NA,nrow(int_table)))), spa)

na<-paste("_",vari,"_",sep="")
count<-1
ooi_pos<-grep(ooi, colnames(p_table))
for(jj in 1:spa){

out<-rep((list(rep(NA,nrow(int_table)))), spal[jj])

for(i in count:(count+spal[jj]-1)){
	
	temp1<-as.numeric(strsplit(vari[i],"_")[[1]])
	temp<-na.omit(int_table[,temp1])
	wtemp<-weight[temp1,2]
	position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	print(i)
	dat2<-qtrans(temp)
    a<-rowSums(t(t(dat2)*wtemp))
    b<-sum(wtemp)^2
    d<-sum(wtemp^2)
	l<-length(a)
    k<-1/l
	y<-k * seq(1,l)
    dat3<-data.frame(y,a,d,b)
    
	rhotemp<-coef(nls(y~sort(pnorm(a/sqrt((1-rho)*d+rho*b), lower.tail=F)), data=dat3, start=list(rho=0), control=list(minFactor = 1/128, tol = 1e-05, warnOnly = TRUE)))
	
	temp<-na.omit(p_table[position,temp1])
	wtemp<-weight[temp1,2]
	position<-as.numeric(strsplit(posvec2[i],"_")[[1]])
	
	dat2<-qtrans(temp)
    a<-rowSums(t(t(dat2)*wtemp))
    b<-sum(wtemp)^2
    d<-sum(wtemp^2)
	
	out[[i-count+1]][position]<-pnorm(a/sqrt((1-rhotemp)*d+rhotemp*b),lower.tail=F)
	
	
}
	names(out)<-na[count:(count+spal[jj]-1)]
	save(out, file=paste("ooi_ooi_cons_full_table_",jj,".Rdata",sep=""))
	temp<-grep(paste("_",ooi_pos,"_",sep=""), na[count:(count+spal[jj]-1)])
	if(length(temp)>0){
	ooilist[[jj]]<-do.call(pmin, c(out[temp],list(na.rm=T)))
	}

count<-count+spal[jj]

}
################################

out_ooi<-do.call(pmin, c(ooilist,list(na.rm=T)))

out_ooi_fdr<-p.adjust(out_ooi, method="BH")
out_ooi<-cbind(out_ooi_fdr,out_ooi)
colnames(out_ooi)<-c("fdr","p-value")
out_ooi<-cbind(out_ooi, da[,3:ncol(da)])

initial_sorting<-seq(1,nrow(p_table))
out_ooi<-cbind(out_ooi, initial_sorting)
out_ooi<-out_ooi[order(as.numeric(out_ooi[,2])),]

write.table(out_ooi, file="ooi_ooi_cons_close_orgs.csv",sep=",", quote=F, row.names=F)


}