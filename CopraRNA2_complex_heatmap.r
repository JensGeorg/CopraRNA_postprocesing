sel="genelist.txt" # case sensitve list of locus tags
inputfile="CopraRNA2_final_all_ooi_ooiconsensus.csv"
num=25
consensus="overall" # alternative: "ooi"
select=F
clustering="ribosomal" # alternative: "ribosomal", "default"
coprarna_reference_file="CopraRNA_available_organisms.txt"
int_p_thres=0.35
int_p_thres2=0.5
prefix="sRNA"



require(ComplexHeatmap)
require(ape) 
require(phangorn)
require(seqinr)
#selected_heatmap<-function(int_p_thres=0.35,int_p_thres2=0.5, sel="genelist.txt", select=T, clustering="ribosomal", consensus="overall", inputfile="CopraRNA2_final_all_ooi.csv", num=25,coprarna_reference_file="CopraRNA_available_organisms.txt", prefix="sRNA"){

clus=TRUE
load("heatmap_data.Rdata")
out_table2<-res[[1]]
p_table<-res[[2]]
mRNA_table<-res[[3]]
sRNA_table<-res[[4]]
opt_sub_table<-res[[5]]



evo_analysis<-read.csv(inputfile,sep=",", header=T) 
selection<-evo_analysis[1:num,"initial_sorting"]



evo_analysis<-read.csv("CopraRNA2_prep_anno_addhomologs_padj_amountsamp.csv",sep=",", header=T) 


as.hclust.phylo2 <- function(x, ...)
{
    if (!is.ultrametric(x)) stop("the tree is not ultrametric")
    if (!is.binary.phylo(x)) stop("the tree is not binary")
    if (!is.rooted(x)) stop("the tree is not rooted")
    n <- length(x$tip.label)
    x$node.label <- NULL # by Jinlong Zhang (2010-12-15)
    bt <- branching.times(x)
    N <- n - 1L

    x <- reorder(x, "postorder")
    m <- matrix(x$edge[, 2], N, 2, byrow = TRUE)
    anc <- x$edge[c(TRUE, FALSE), 1]
    bt <- bt[as.character(anc)] # 1st, reorder
    ## 2nd, sort keeping the root branching time in last (in case of
    ## rounding error if there zero-lengthed branches nead the root)
    bt <- c(sort(bt[-N]), bt[N])
    o <- match(names(bt), anc)
    m <- m[o, ]

    ## first renumber the tips:
    TIPS <- m <= n
    m[TIPS] <- -m[TIPS]

    ## then renumber the nodes:
    oldnodes <- as.numeric(names(bt))[-N]
    m[match(oldnodes, m)] <- 1:(N - 1)

    names(bt) <- NULL
    obj <- list(merge = m, height = 2*bt, order = 1:n, labels = x$tip.label,
                call = match.call(), method = "unknown")
    class(obj) <- "hclust"
    obj
}

mafft<-function(filename="ncrna.fa", outname="ncrna_aligned.fa", mode="accurate"){
	if(mode=="accurate"){
		command<-paste("mafft --maxiterate 1000 --localpair --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	
	if(mode=="fast"){
		command<-paste("mafft --retree 2 --maxiterate 0 --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	if(mode=="very_fast"){
		command<-paste("mafft --retree 1 --maxiterate 0 --quiet --inputorder ", filename, " > ", outname, sep="" )
	}
	system(command)
	#fas<-read.fasta(outname)
	#fas
} 

find_sRNA<-function(x, y){
	e<-grep("Annotation", colnames(y))-1
	temp<-evo_analysis[,3:e]
	genes<-apply(temp, 1, paste, collapse="")
	x<-as.matrix(x)
	out<-unlist(apply((x), 1, find1, genes))
	out
}



find1<-function(x, genes){
	temp<-grep(x, genes, ignore.case=F)
	temp
}

if(select==T){
	genelist<-read.csv(sel, sep="\t", header=F)[,1]
	selection<-unique(find_sRNA(genelist, evo_analysis))
}




ee<-grep("Annotation", colnames(evo_analysis))-1
colnames(out_table2)<-colnames(evo_analysis)[3:ee]
copref<-read.delim(coprarna_reference_file, sep="\t", header=T,comment.char = "#")	
nam<-c()
for(i in 1:ncol(out_table2)){
	tnam<-grep(gsub("\\..*","",colnames(out_table2)[i]),copref[,1])
	nam<-c(nam,as.character(copref[tnam,2]))
	
	
}

nam2<-c()
for(i in 1:length(nam)){
	temp<-substr(nam[i],1,3)
	temp2<-strsplit(nam[i],"_")[[1]]
	temp<-paste(temp,"_",temp2[2], sep="")
	if(length(temp2)>2){
		temp<-paste(temp, temp2[length(temp2)], sep="_")
	}
	nam2<-c(nam2,temp)
}
nam2<-paste(nam2, colnames(out_table2), sep="_")



gene_anno<-c()
evo_analysis<-as.matrix(evo_analysis)
for(i in 1:nrow(evo_analysis)){

	e<-grep("Annotation", colnames(evo_analysis))-1
	temp<-evo_analysis[i,3:e]
	temp<-temp[which(temp!="")][1]
	locus_tag<-gsub("\\(.*","",as.character(temp))
	genename<-gsub(".*\\(","",as.character(temp))
	genename<-gsub("\\|.*", "", genename)
	genename<-gsub("\\/", "", genename)
	na<-paste(genename, locus_tag,sep="_")
	gene_anno<-c(gene_anno, na)
}



l<-2
len<-length(which(duplicated(gene_anno)))
while(len>0){
	#print(l)
	gene_anno[which((duplicated(gene_anno)))]<-paste(gene_anno[which((duplicated(gene_anno)))],l,sep="_")
	l<-l+1
	len<-length(which(duplicated(gene_anno)))
}
row.names(out_table2)<-gene_anno


if(clustering=="sRNA"){
command<-paste("clustalo -i ", "input_sRNA.fa", " --distmat-out=distmatout1.txt --full --use-kimura --output-order=input-order --force --max-hmm-iterations=-1", sep="")
system(command)
temp<-read.delim("distmatout1.txt",sep="",header=F, , skip=1)
unlink("distmatout1.txt")
na<-temp[,1]
temp<-temp[,2:ncol(temp)]
colnames(temp)<-na
rownames(temp)<-na
dis<-as.dist(temp) 
clus<-(hclust(dis,method="average"))

ord<-clus$label

ord2<-match(ord, gsub("\\..*","",colnames(out_table2)))




clus$label<-nam2
}
if(clustering=="ribosomal"){

	# rRNA<-matrix(,nrow(coor2),2)

		# con <- dbConnect(SQLite(), dbname="/home/jens/jensSicherung/Paper und Drafts/GLASSgo/SILVA_132_SSUParc_tax_silva.fasta/rRNA.db", ":memory:")
		# # #rRNA_data <- dbReadTable(con, 'rRNA_data')
		# # ids<-rRNA_data[,1]
		# # dbWriteTable(con, "IDs",ids)
		# #
		# for(j in 1:nrow(coor2)){
		# #for(i in 1:2){
			# temp<-gsub("\\..*","",coor2[j,1])
			# temp<-paste("SELECT sequence   FROM rRNA_data WHERE ID == '", temp,"'",sep="")
			# #print(temp)
			# temp<-dbGetQuery(con,temp )
			# #print(temp)
			# rRNA[j,2]<-as.character(temp[1,1])
			# #dbClearResult(temp)
			# rRNA[j,1]<-coor2[j,"fin"]
		# }

		# dbDisconnect(con)

		# fasta<-c()
		# for(j in 1:nrow(rRNA)){
			# fasta<-c(fasta,paste(">", rRNA[j,1], sep=""),gsub("u","t",rRNA[j,2]))

		

		# write.table(fasta, file="16s.fasta", row.names=F, col.names=F, quote=F)


		mafft(filename="16s_sequences.fa")
		
	
	tempf<-read.fasta("ncrna_aligned.fa")
	write.fasta(tempf, file.out="ncrna_aligned.fa", names=names(tempf), nbchar=100000)
	dat<-read.phyDat("ncrna_aligned.fa", format="fasta", type="DNA")
	dm <- dist.ml(dat, model="F81")
	treeNJ <- NJ(dm)
	fitJC = pml(treeNJ, data=dat)
	fit2<-chronos(fitJC$tree)
	clus<-as.hclust.phylo2(midpoint(fit2))


# d<-read.delim("distmat.out", sep="\t", skip=6, header=FALSE)
 # dd<-d[2:nrow(d),2:(ncol(d)-2)]
 # rownames(dd)<-gsub(" .*","",d[2:nrow(d),ncol(d)])
# colnames(dd)<-gsub(" .*","",d[2:nrow(d),ncol(d)])
 # dis<-as.dist(t(dd))
 
 
 # clus<-(hclust(dis,method="average"))

ord<-clus$label

ord2<-match(ord, gsub("\\..*","",colnames(out_table2)))




clus$label<-nam2
 
}
out_table2<-out_table2[selection,]
out_table2<-out_table2[,ord2]
p_table<-p_table[selection,]
mRNA_table<-mRNA_table[selection,]
sRNA_table<-sRNA_table[selection,]
opt_sub_table<-opt_sub_table[selection,]
p_table<-p_table[,ord2]
mRNA_table<-mRNA_table[,ord2]
sRNA_table<-sRNA_table[,ord2]
opt_sub_table<-opt_sub_table[,ord2]
colnames(out_table2)<-nam2
#require(pheatmap)

my_palette<-colorRampPalette(c("olivedrab2","olivedrab3","goldenrod1","orange","orangered2","orangered3","white"))(n=7)

nam<-paste(prefix,"conservation_heatmap.pdf", sep="_" )


 
col_fun = circlize::colorRamp2(c(-1,-0.9,0, int_p_thres,int_p_thres+0.1, 0.5,1), c("white","darkolivegreen1","darkolivegreen1","olivedrab3", "orange", "orangered2","orangered2"), space = "RGB")
#my_palette<-colorRampPalette(c("lightblue1","lightblue2","lightyellow1","lightyellow2","gray88","gray94","white"))(n=7)
 
lab<-c(	paste("p<", int_p_thres, " m+sRNA site cons.", sep=""),
		paste("p<", int_p_thres, " mRNA site cons.", sep=""),
		paste(int_p_thres,">p<", int_p_thres2, " mRNA site cons.", sep=""),
		paste("p<", int_p_thres, " mRNA site NOT cons.", sep=""),
		paste("p>", int_p_thres2, " mRNA site cons.", sep=""),
		paste("p>", int_p_thres, " mRNA site NOT cons.", sep=""),
		"no homolog")		


	

ha = HeatmapAnnotation(	                       
							intarna = rep(-1,ncol(out_table2)),#seq(0,1,by=0.1), 
							col = list(
										intarna = circlize::colorRamp2(c(-1,-0.9,0, int_p_thres, 0.5,1), c("white","olivedrab2","olivedrab2", "orange", "orangered2","orangered2"))),
										which=c("row"),
										width=unit(0, "npc"),
										gap = unit(0, "npc"),
										annotation_width = 0,
										annotation_legend_param = list(
																		
																		intarna=list(title="IntaRNA p-value\n(circular inlay)",at=seq(0,1,by=0.1),labels_gp = gpar(cex=0.5), title_gp = gpar(cex=0.8))
																		)		
							
							)							

pdf(nam, width = 2.3+ncol(out_table2)*0.2, height = 3.4+nrow(out_table2)*0.2,useDingbats=F) ## edit prw // onefile=FALSE fix for two page pdf

a<-Heatmap(	out_table2,
			
			col=my_palette,
			#top_annotation_height =  unit(0, "mm"),
			cluster_rows = F, 
			cluster_columns =clus, 
			column_title=paste(consensus, "_consensus_based"),
			show_heatmap_legend = T,
			heatmap_legend_param = list(title = "Target prediction",at = seq(1,7), labels = lab,labels_gp = gpar(cex=0.5), gp=list(col=my_palette),color_bar = "discrete"),
			#rect_gp = gpar(col = "darkgrey", lty = 1, lwd = 2),
			column_names_gp = gpar( rot = 30, cex=0.7),
			row_names_gp = gpar( rot = 0, cex=0.75),
			rect_gp=gpar(col="grey",lwd=0.5,lty=1),
		
				
			
			cell_fun=function(j, i, x, y, width=width, height=width, fill){
				s = min(unit.c(convertWidth(width, "cm"), convertHeight(height, "cm")))
				
				#rect_gp = gpar(col = "lightgrey", lty = 1, lwd = 0.15)
				if(out_table2[i,j]!=7){
					#print(c(width,height))
					grid.circle(x = x-0.25*unit.c(width), y = y+0.25*unit.c(height), r = max(unit.c(width),unit.c(height))*0.18,  
						gp = gpar(fill = col_fun((as.numeric(p_table[i, j]))), col = "grey",lwd=0.3))
						
					#
					if(opt_sub_table[i, j]=="sub"){
						grid.text(opt_sub_table[i, j], x = x, y = y-0.3*unit.c(height), gp = gpar(cex=0.3, col="gray34"))	
					}
					if(mRNA_table[i, j]==TRUE & sRNA_table[i, j]==TRUE){
						#grid.text("+/+", x = x, y = y-0.1*unit.c(height), gp = gpar(cex=0.55, col="gray19"))
					
					}
					if(mRNA_table[i, j]==F & sRNA_table[i, j]==F){
						grid.text("-/-", x = x, y = y-0.1*unit.c(height), gp = gpar(cex=0.55, col="gray34"))
					
					}
					if(mRNA_table[i, j]==T & sRNA_table[i, j]==F){
						grid.text("+/-", x = x, y = y-0.1*unit.c(height), gp = gpar(cex=0.55, col="gray34"))
					
					}
					if(mRNA_table[i, j]==F & sRNA_table[i, j]==T){
						grid.text("-/-", x = x, y = y-0.1*unit.c(height), gp = gpar(cex=0.55, col="gray34"))
					
					}
				}
				}
			
			)
add_heatmap(a,ha)
#draw(a)
dev.off()





#}
#selected_heatmap(select=select, clustering=clustering, consensus=consensus, inputfile=inputfile, sel=sel, num=num, coprarna_reference_file=coprarna_reference_file,int_p_thres=int_p_thres, int_p_thres2=int_p_thres2,prefix=prefix)

