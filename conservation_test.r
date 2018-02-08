
# script by Jens Georg

# R --slave -f  ../conservation_test.r --args thres=0.35 thres2=0.5 thres_overlap=0.5




args <- commandArgs(trailingOnly = TRUE)


thres<-0.32 # IntaRNA p-value threshold for positives
thres2<- 0.51 # IntaRNA p-value threshold for unsures
thres_overlap<-0.45 # Threshold for overlap with consensus sites in mRNA and sRNA



if (length(args)>0) {
 for(i in 1:length(args)){
	temp<-strsplit(args[i],"=")
	temp<-temp[[1]]
	temp1<-temp[1]
	temp2<-temp[2]
	assign(as.character(temp1),temp2)
 }
}


thres<-as.numeric(thres)
thres2<-as.numeric(thres2)
thres_overlap<-as.numeric(thres_overlap)

# The sRNA/target interaction in a organism is considered as conserved 
# if the respective IntaRNA p-value is <= a given threshold (default: 0.35)
# and the interaction sites in the mRNA and sRNA matches the consensus. 
# conserved: IntaRNA p-value <= threshold + mRNA and sRNA match consensus (opt or subopt)
# conserved2: IntaRNA p-value <= threshold + only mRNA matches the consensus(opt or subopt)
# unsure: IntaRNA p-value <= threshold but no consensus match in mRNA
# unsure2: IntaRNA p-value > threshold but <= threshold 2 and mRNA matches consensus
# not conserved: IntaRNA p-value > threshold 2 but mRNA matches consensus
# not_conserved2: IntaRNA p-value > threshold + no consensus match in mRNA




load("conservation_table.Rdata") # for p_values
con_table<-conservation_table[[1]]
con_table_sub<-conservation_table[[2]]

load("interaction_positions.Rdata") # for predicted interaction sites 
mRNA_sites<-interaction_positions[[1]]
mRNA_sites_sub<-interaction_positions[[2]]
sRNA_sites<-interaction_positions[[3]]
sRNA_sites_sub<-interaction_positions[[4]]

load("consensus_positions.Rdata") # positions of consensus interaction sites
consensus_mRNA<-consensus_both[[1]]
consensus_sRNA<-consensus_both[[2]]


out_table<-con_table
out_table[]<-NA

p_table<-out_table			# relevant p-value, default optimal p-value besides optimal prediction does not cover the consensus in mRNA and supoptimal prediction does and suboptimal p-vale < thres
mRNA_table<-out_table		# predicted interaction matches the mRNA consensus (TRUE/FALSE)
sRNA_table<-out_table		# predicted interaction matches the sRNA consensus (TRUE/FALSE)
opt_sub_table<-out_table	# information of optimal or suboptimal predictions are used ("opt", "sub")


overlap<-function(int, cons, ov_thres=0.45){ # int= vector with 2 elements, start + end of predicted site; cons= vector with 2 elements, start + end of consensus site
	out<-FALSE
	site<-seq(int[1],int[2])
	con<-seq(cons[1],cons[2])
	len_site<-length(site)
	len_cons<-length(con)
	inter<-intersect(con,site)
	over<-length(inter)/min(len_site,len_cons)
	if(over>=ov_thres){
		out<-TRUE
	}
	out
}

# cons_test<-function(p, p_sub, msites, msites_sub, ssites, ssites_sub, con_mRNA, con_sRNA, thres, thres2, ov_thres){  # classifies the targets as conserved, unsure or not conserved
	
	# out<-rep(NA, length(p))
	# out_p<-rep(NA, length(p))
	# out_mRNA<-rep(NA, length(p))
	# out_sRNA<-rep(NA, length(p))
	# out_opt_sub<-rep(NA, length(p))
	
	# for(i in 1:length(p)){
		# #print(i)
		# if(is.na(p[[i]][1])==F){	# if p[[i]][1] is NA then domclust detected no homolog in this organism
			# opt<-p[[i]][1]			# optimal p_value
			# subopt<-p_sub[[i]][1]	# sub optimal p_value
			# m_sub<-FALSE
			# s_sub<-FALSE
			# s2_sub<-FALSE
			# if(is.na(subopt)==F){
				# m_sub<-overlap(as.numeric(msites_sub[[i]]),con_mRNA ,ov_thres)
				# s_sub<-overlap(as.numeric(ssites_sub[[i]]),con_sRNA[[1]] ,ov_thres)
				
				# if(is.null(con_sRNA[[2]])==F){
					# s2_sub<-overlap(as.numeric(ssites_sub[[i]]),con_sRNA[[2]] ,ov_thres)
				# }
			# }
			# if(is.na(subopt)==T){
				# subopt<-1
			# }
			
			# # test for overlap with consensus site 
			# m<-overlap(as.numeric(msites[[i]]),con_mRNA ,ov_thres)
			# s<-overlap(as.numeric(ssites[[i]]),con_sRNA[[1]] ,ov_thres)
			# s2<-FALSE
			# if(is.null(con_sRNA[[2]])==F){
				# s2<-overlap(as.numeric(ssites[[i]]),con_sRNA[[2]] ,ov_thres)
			# }
			
			
			# # classification based on p_value and predicted sites
			
			# # conserved: IntaRNA p-value <= threshold + mRNA and sRNA match consensus opt
			# if((opt<=thres & m==T & s==T) | (opt<=thres & m==T & s2==T)){
				# out[i]<-"conserved"
				# out_p[i]<-opt
				# out_mRNA[i]<-TRUE
				# out_sRNA[i]<-TRUE
				# out_opt_sub[i]<-"opt"
			# }
			# # conserved: IntaRNA p-value <= threshold + mRNA and sRNA match consensus subopt
			# if( (opt>thres & subopt<=thres & m_sub==T & s_sub==T) | (opt>thres & subopt<=thres & m_sub==T & s2_sub==T) | (m==F & subopt<=thres & m_sub==T & s_sub==T)){
				# out[i]<-"conserved"
				# out_p[i]<-subopt
				# out_mRNA[i]<-TRUE
				# out_sRNA[i]<-TRUE
				# out_opt_sub[i]<-"sub"
			# }
			
			
			# # conserved2: IntaRNA p-value <= threshold + only mRNA matches the consensus opt
			# if((opt<=thres & m==T & s==F & s2==F) ){
				# out[i]<-"conserved2"
				# out_p[i]<-opt
				# out_mRNA[i]<-TRUE
				# out_sRNA[i]<-FALSE
				# out_opt_sub[i]<-"opt"
			# }
			# # conserved2: IntaRNA p-value <= threshold + only mRNA matches the consensus subopt
			# if( (m==F & subopt<=thres & m_sub==T & s_sub==F & s2_sub==F) | (opt>thres & subopt<=thres & m_sub==T & s_sub==F & s2_sub==F) ){
				# out[i]<-"conserved2"
				# out_p[i]<-subopt
				# out_mRNA[i]<-TRUE
				# out_sRNA[i]<-FALSE
				# out_opt_sub[i]<-"sub"
			# }
			
			
			
			# # unsure2: IntaRNA p-value <= threshold but no consensus match in mRNA opt --> false site
			# if((opt<=thres & m==F & m_sub==F )  ){
				# out[i]<-"unsure2"
				# out_p[i]<-opt
				# out_mRNA[i]<-FALSE
				# out_sRNA[i]<-FALSE
				# if(s==T | s2==T){
					# out_sRNA[i]<-TRUE
				# }	
				# out_opt_sub[i]<-"opt"				
			# }
			
			
			# # unsure: IntaRNA p-value > threshold but <= threshold 2 and mRNA matches consensus opt --> too high p-value
			# if((opt>thres & opt<=thres2 & m==T)){
				# out[i]<-"unsure"
				# out_p[i]<-opt
				# out_mRNA[i]<-TRUE
				# out_sRNA[i]<-FALSE
				# if(s==T | s2==T){
					# out_sRNA[i]<-TRUE
				# }				
				# out_opt_sub[i]<-"opt"				
			# }
			
			# # unsure: IntaRNA p-value > threshold but <= threshold 2 and mRNA matches consensus subopt --> too high p-value
			# if((subopt>thres & subopt<=thres2 & m_sub==T & m==F ) | (subopt>thres & subopt<=thres2 & m_sub==T & opt>thres2) ){
				# out[i]<-"unsure"
				# out_p[i]<-subopt
				# out_mRNA[i]<-TRUE
				# out_sRNA[i]<-FALSE
				# if(s_sub==T | s2_sub==T){
					# out_sRNA[i]<-TRUE
				# }				
				# out_opt_sub[i]<-"sub"				
			# }
			
			
			# # not conserved: IntaRNA p-value > threshold 2 but mRNA matches consensus opt --> much too high p-value
			# if(opt>thres2  & m==T ){
				# out[i]<-"not_conserved"
				# out_p[i]<-opt
				# out_mRNA[i]<-TRUE
				# out_sRNA[i]<-FALSE
				# if(s==T | s2==T){
					# out_sRNA[i]<-TRUE
				# }				
				# out_opt_sub[i]<-"opt"
				

			# }
			
			
			
			# # not conserved: IntaRNA p-value > threshold 2 but mRNA matches consensus subopt --> much too high p-value
			# if( (opt>thres2 & m_sub==T & m==F)| (subopt>thres2 & m_sub==T & m==F & opt>thres)){
				# out[i]<-"not_conserved"
				# out_p[i]<-subopt
				# out_mRNA[i]<-TRUE
				# out_sRNA[i]<-FALSE
				# if(s_sub==T | s2_sub==T){
					# out_sRNA[i]<-TRUE
				# }				
				# out_opt_sub[i]<-"sub"				
			# }
			
			
			
			# # not_conserved2: IntaRNA p-value > threshold + no consensus match in mRNA
			# if((opt>thres  & m==F & m_sub==F)){
				# out[i]<-"not_conserved2"
				# out_p[i]<-opt
				# out_mRNA[i]<-FALSE
				# out_sRNA[i]<-FALSE
				# if(s==T | s2==T ){
					# out_sRNA[i]<-TRUE
				# }				
				# out_opt_sub[i]<-"opt"
			# }
			
		# }
	# }
	# out<-list(out,out_p,out_mRNA,out_sRNA,out_opt_sub)
# }

cons_test<-function(p, p_sub, msites, msites_sub, ssites, ssites_sub, con_mRNA, con_sRNA, thres, thres2, ov_thres){  # classifies the targets as conserved, unsure or not conserved
	
	out<-rep(NA, length(p))
	out_p<-rep(NA, length(p))
	out_mRNA<-rep(NA, length(p))
	out_sRNA<-rep(NA, length(p))
	out_opt_sub<-rep(NA, length(p))
	
	for(i in 1:length(p)){
		#print(i)
		if(is.na(p[[i]][1])==F){	# if p[[i]][1] is NA then domclust detected no homolog in this organism
			opt<-p[[i]][1]			# optimal p_value
			subopt<-p_sub[[i]][1]	# sub optimal p_value
			m_sub<-FALSE
			s_sub<-FALSE
			s2_sub<-FALSE
			if(is.na(subopt)==F){
				m_sub<-overlap(as.numeric(msites_sub[[i]]),con_mRNA ,ov_thres)
				s_sub<-overlap(as.numeric(ssites_sub[[i]]),con_sRNA[[1]] ,ov_thres)
				
				if(is.null(con_sRNA[[2]])==F){
					s2_sub<-overlap(as.numeric(ssites_sub[[i]]),con_sRNA[[2]] ,ov_thres)
				}
			}
			if(is.na(subopt)==T){
				subopt<-1
			}
			
			# test for overlap with consensus site 
			m<-overlap(as.numeric(msites[[i]]),con_mRNA ,ov_thres)
			s<-overlap(as.numeric(ssites[[i]]),con_sRNA[[1]] ,ov_thres)
			s2<-FALSE
			if(is.null(con_sRNA[[2]])==F){
				s2<-overlap(as.numeric(ssites[[i]]),con_sRNA[[2]] ,ov_thres)
			}
			
			
			# classification based first on consensus site in mRNA, second on p-value
						
			# optimal mRNA hit positive
			if(m==T){						
				if(opt<=thres){
					out[i]<-"conserved2"
					out_p[i]<-opt
					out_mRNA[i]<-TRUE
					out_sRNA[i]<-FALSE
					if(s==T | s2==T){
						out_sRNA[i]<-TRUE
						out[i]<-"conserved"
					}				
					out_opt_sub[i]<-"opt"	
				}
				if(opt>thres & opt<=thres2){
					out[i]<-"unsure"
					out_p[i]<-opt
					out_mRNA[i]<-TRUE
					out_sRNA[i]<-FALSE
					if(s==T | s2==T){
						out_sRNA[i]<-TRUE
					}				
					out_opt_sub[i]<-"opt"	
				}
				if(opt>thres2 ){
					out[i]<-"not_conserved"
					out_p[i]<-opt
					out_mRNA[i]<-TRUE
					out_sRNA[i]<-FALSE
					if(s==T | s2==T){
						out_sRNA[i]<-TRUE
					}				
					out_opt_sub[i]<-"opt"	
				}
			}	
			
			# optimal mRNA hit negative, suboptimal hit positive			
			if(m==F & m_sub==T){				
				if(subopt<=thres){
					out[i]<-"conserved2"
					out_p[i]<-subopt
					out_mRNA[i]<-TRUE
					out_sRNA[i]<-FALSE
					if(s_sub==T | s2_sub==T){
						out_sRNA[i]<-TRUE
						out[i]<-"conserved"
					}				
					out_opt_sub[i]<-"sub"	
				}
				if(subopt>thres & subopt<=thres2){
					out[i]<-"unsure"
					out_p[i]<-subopt
					out_mRNA[i]<-TRUE
					out_sRNA[i]<-FALSE
					if(s_sub==T | s2_sub==T){
						out_sRNA[i]<-TRUE
					}				
					out_opt_sub[i]<-"sub"	
				}
				if(subopt>thres2){
					out[i]<-"not_conserved"
					out_p[i]<-subopt
					out_mRNA[i]<-TRUE
					out_sRNA[i]<-FALSE
					if(s_sub==T | s2_sub==T){
						out_sRNA[i]<-TRUE
					}				
					out_opt_sub[i]<-"sub"	
				}
			}
				
				
			# both optimal mRNA hit and suboptimal hit negative
				
			if(m==F & m_sub==F){
				if(opt<=thres){
					out[i]<-"unsure2"
					out_p[i]<-opt
					out_mRNA[i]<-FALSE
					out_sRNA[i]<-FALSE
					if(s==T | s2==T){
						out_sRNA[i]<-TRUE						
					}				
					out_opt_sub[i]<-"opt"	
				}
				if(opt>thres){
					out[i]<-"not_conserved2"
					out_p[i]<-opt
					out_mRNA[i]<-FALSE
					out_sRNA[i]<-FALSE
					if(s==T | s2==T){
						out_sRNA[i]<-TRUE						
					}				
					out_opt_sub[i]<-"opt"	
				}
			}			
		}
	}
	out<-list(out,out_p,out_mRNA,out_sRNA,out_opt_sub)
}				
			


for(i in 1:nrow(con_table)){
	#print(i)
	p_temp<-con_table[i,]
	p_temp<-strsplit(p_temp, "\\|")
	p_temp_sub<-con_table_sub[i,]
	p_temp_sub<-strsplit(p_temp_sub, "\\|")
	
	msite_temp<-strsplit(mRNA_sites[i,], "\\|")
	msite_temp_sub<-strsplit(mRNA_sites_sub[i,], "\\|")
	ssite_temp<-strsplit(sRNA_sites[i,], "\\|")
	ssite_temp_sub<-strsplit(sRNA_sites_sub[i,], "\\|")
	
	cons_m<-consensus_mRNA[[i]]
	cons_s<-consensus_sRNA[[i]]
	
	temp<-cons_test(p_temp, p_temp_sub, msite_temp,msite_temp_sub, ssite_temp, ssite_temp_sub,cons_m, cons_s, thres, thres2,thres_overlap)
	out_table[i,]<-temp[[1]]
	p_table[i,]<-temp[[2]]	
	mRNA_table[i,]<-temp[[3]]
	sRNA_table[i,]<-temp[[4]]
	opt_sub_table[i,]<-temp[[5]]
}

out_table2<-matrix(,nrow(out_table),ncol(out_table))
out_table2[which(out_table=="conserved")]<-1
out_table2[which(out_table=="conserved2")]<-2
out_table2[which(out_table=="unsure")]<-3
out_table2[which(out_table=="unsure2")]<-4
out_table2[which(out_table=="not_conserved")]<-5
out_table2[which(out_table=="not_conserved2")]<-6
out_table2[which(is.na(out_table)==T)]<-7


res<-list(out_table2,p_table,mRNA_table,sRNA_table,opt_sub_table,out_table)
save(res,file="heatmap_data.Rdata")


evo_table<-read.csv("CopraRNA2_final_all_evo.csv", sep="\t")  #path needs to be adapted to clean output

evo<-out_table2[evo_table[,"initial_sorting"],]

evo2<-matrix(,nrow(evo),4)
colnames(evo2)<-c("%_conserved","%_unsure","%_not_conserved","%_no_homolog")

for(i in 1:nrow(evo)){
	evo2[i,1]<-length(which(evo[i,]==1 | evo[i,]==2))/ncol(evo)
	evo2[i,2]<-length(which(evo[i,]==3 | evo[i,]==4))/ncol(evo)
	evo2[i,3]<-length(which(evo[i,]==5 | evo[i,]==6))/ncol(evo)
	evo2[i,4]<-length(which(evo[i,]==7))/ncol(evo)
}


evo_table<-cbind(evo2,evo_table[,1],evo_table[,5:ncol(evo_table)])

write.table(evo_table, file="CopraRNA2_final_all_evo.csv", sep="\t", row.names=F)


