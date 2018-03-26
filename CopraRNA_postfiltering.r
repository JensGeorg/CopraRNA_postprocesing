

# post filter copraRNA results based on "conservation classification"


load("heatmap_data.Rdata")
out_table2<-res[[1]]


ooi<-read.csv("CopraRNA2_final_all_ooi.csv", sep=",")
ooi_cons<-read.csv("CopraRNA2_final_all_ooi_ooiconsensus.csv",sep=",")
bal<-read.csv("CopraRNA2_final_all_balanced.csv", sep=",")
bal_cons<-read.csv("CopraRNA2_final_all_balanced_consensus.csv", sep=",")

dat<-list(ooi,ooi_cons,bal,bal_cons)
neg<-which(out_table2[,1]>4)
for(i in 1:length(dat)){
	dat[[i]]<-dat[[i]][order(dat[[i]][,"initial_sorting"]),]
	if(length(neg)>0){
		dat[[i]]<-dat[[i]][-neg,]
		dat[[i]]<-dat[[i]][order(dat[[i]][,"p.value"]),]
	}
}



# remove duplicate clusters

for(i in 1:length(dat)){
	temp<-dat[[i]][,3]
	temp<-gsub("\\(.*","",temp)
	temp<-which(duplicated(temp))
	if(length(temp)>0){
		dat[[i]]<-dat[[i]][-temp,]
	}
}


for(i in 1:length(dat)){
	nam<-paste(nam1[i],"_filtered.csv",sep="")
	write.table( dat[[i]],file=nam ,sep=",", col.names=F)


}
