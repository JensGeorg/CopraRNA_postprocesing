
sel="genelist.txt" # case sensitve list of locus tags
inputfile="CopraRNA2_final_all_ooi.csv"
num=25
consensus="overall" # alternative: "ooi"
select=TRUE
clustering="ribosomal" # alternative: "ribosomal", "default"
coprarna_reference_file="CopraRNA_available_organisms.txt"
int_p_thres=0.35
int_p_thres2=0.5
prefix="sRNA"



require(ComplexHeatmap)
require(ape) 
require(phangorn)

selected_heatmap<-function(int_p_thres=0.35,int_p_thres2=0.5, sel="genelist.txt", select=T, clustering="ribosomal", consensus="overall", inputfile="CopraRNA2_final_all_ooi.csv", num=25,coprarna_reference_file="CopraRNA_available_organisms.txt", prefix="sRNA"){

clus=TRUE
load("heatmap_data.Rdata")
out_table2<-res[[1]]
p_table<-res[[2]]
mRNA_table<-res[[3]]
sRNA_table<-res[[4]]
opt_sub_table<-res[[5]]



evo_analysis<-read.csv(inputfile,sep=",", header=T) 
selection<-evo_analysis[1:num,"initial_sorting"]



evo_analysis<-evo_analysis[order(evo_analysis[,"initial_sorting"]),]

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
d<-read.delim("distmat.out", sep="\t", skip=6, header=FALSE)
 dd<-d[2:nrow(d),2:(ncol(d)-2)]
 rownames(dd)<-gsub(" .*","",d[2:nrow(d),ncol(d)])
colnames(dd)<-gsub(" .*","",d[2:nrow(d),ncol(d)])
 dis<-as.dist(t(dd))
 
 
 clus<-(hclust(dis,method="average"))

ord<-clus$label

ord2<-match(ord, gsub("\\..*","",colnames(out_table2)))

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


clus$label<-nam2
 
}
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

pdf(nam,paper = "a4r", width = 0, height = 0,useDingbats=F) ## edit prw // onefile=FALSE fix for two page pdf

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
					grid.circle(x = x-0.25*unit.c(width), y = y+0.25*unit.c(height), r = unit.c(width)*0.3, 
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





}
selected_heatmap(select=select, clustering=clustering, consensus=consensus, inputfile=inputfile, sel=sel, num=num, coprarna_reference_file=coprarna_reference_file,int_p_thres=int_p_thres, int_p_thres2=int_p_thres2,prefix=prefix)








cor_mat = cor(mat)
od = hclust(dist(cor_mat))$order
cor_mat = cor_mat[od, od]
nm = rownames(cor_mat)
col_fun = circlize::colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
# `col = col_fun` here is used to generate the legend
Heatmap(cor_mat, name = "correlation", col = col_fun, rect_gp = gpar(type = "none"), 
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
        if(i == j) {
            grid.text(nm[i], x = x, y = y)
        } else if(i > j) {
            grid.circle(x = x, y = y, r = abs(cor_mat[i, j])/2 * min(unit.c(width, height)), 
                gp = gpar(fill = col_fun(cor_mat[i, j]), col = NA))
        } else {
            grid.text(sprintf("%.1f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 8))
        }
    }, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, show_column_names = FALSE)








