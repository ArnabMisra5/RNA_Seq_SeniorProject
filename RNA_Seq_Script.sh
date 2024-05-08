#!/bin/bash
input1=$(ls /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/Sample_*/*.fastq.gz)
input_files=(${input1[@]})
#echo $input_files[@]}
/home/sheetalshetty/Downloads/FastQC/fastqc --extract --outdir=/media/sheetalshetty/NewDrive1/Arnab_PNPLA7/RawData ${input1}
#md5sum
#Generate checksum (unique value for a files contents)

Fastqc -o /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/PreAlignmentQC -t 16 /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/RawData/*/*.fastq.gz
#fastqc run

STAR --runThreadN 32 --runMode genomeGenerate --genomeDir /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/STAR_Alignment_Files \
--genomeFastaFiles /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/Refrences.HG38/Homo_sapiens.GRCh38.dna.primary_assembly.fasta \
--sjdbGTFfile /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/Refrences.HG38/gencode.v45.chr_patch_hapl_scaff.basic.annotation.gtf --sjdbOverhang 99
#initial star test run

rawsFolder=/media/sheetalshetty/NewDrive1/Arnab_PNPLA7/RawData
myFolderPaths=$(ls $rawsFolder)
myFolders=($myFolderPaths[@])
for name in ${myFolders[29]}; do
	echo $rawsFolder/$name
	fastq1=$rawsFolder/$name/*_R1*.fast.gz
	fastq2=$rawsFolder/$name/*_R2*.fast.gz
	echo $fastq1
	echo $fastq2
	mkdir -p /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/STAR_Outputs/ReMapped/$name
	STAR --runThreadN 32 --genomeDir /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/STAR_Alignment_Files \
	--readFilesIn $fastq1 $fastq2 /
	--readFilesCOmmand gunzip -c --outFileNamePrefix /media/sheetalshetty/NewDrive!/Arnab_PNPLA7/STAR_Outputs/ReMapped/$name/ --outSAMtype BAM unsorted
done
#STAR Alignmentcode

htseq-count --format=bam /media/sheetlashetty/NewDrive1/Arnab_PNPLA7/STAR_Outputs/*.bam /media/sheetalshetty/NewDrive1/Arnab_PNPLA7/Refrences.HG38/Homo_sapiens.GRCh38.111.gtf > AllHTSeq.tsv 
#htseq-count and gff/gtf files

setwd("~/Kruer's Lab CP Research")
pasCts <- read.csv("AllHtseqCount.tsv", sep="\t",row.names = 1) #locate a specific file within a package directory and assign the full path to that file to the variable pasCts
colnames(pasCts) <- c("Sample_PNPLA-742-3_174_019","Sample_PNPLA-742-3_174_020","Sample_PNPLA-742-3_174_021","Sample_GMO8399-GLU-5_144_049","Sample_GMO8399-GLU-4_156_037","Sample_GMO8399-GLU-3_168_025","Sample_GMO8399-GLU-2_180_013","Sample_GMO8339-GLU-1_192_001","Sample_GMO8339-BSA-3_108_085","Sample_GMO8339-BSA-2_120_073","Sample_GMO8339-BSA-1_258_319","Sample_GMO8398-BSA-3_190_003","Sample_GMO8398-BSA-2_107_086","Sample_GMO8398-BSA-1_119_074","Sample_GMO2978-BSA-2_189_004","Sample_GMO2978-BSA-1_106_087","Sample_GMO8447-GLU-3_117_076","Sample_GMO8447-GLU-2_129_064","Sample_GMO8447-GLU-1_141_052","Sample_GMO5565-GLU-2_152_041","Sample_GMO5565-GLU-1_164_029","Sample_GMO3440-GLU-3_176_017","Sample_GMO3440-GLU-2_188_005","Sample_GMO3440-GLU-1_105_088","F-731-6","F-731-5","F-731-3","F-731-2","F-731-1","F731-4") #renaming column names
pasAnno= read.csv("CountsDataReplicated.tsv", sep="\t",row.names = 1) #retrieves the full path to the file "CountsData" within the "Kruer's Lab CP Research" package's "extdata" directory and stores it in the variable pasAnno.
pasAnno$condition= factor(pasAnno$condition) # converts the values in the "condition" column of the coldata data frame into a factor and updates the data frame with this conversion.
pasAnno$type = factor(pasAnno$type) #
pasAnno$sample = factor(pasAnno$sample)

all(rownames(pasAnno) %in% colnames(pasCts)) #checks if all the row names of coldata are present as column names in cts. If all row names are present in the column names, it will return TRUE; otherwise, it will return FALSE


library("DESeq2") #opens DESeq2 package
dds = DESeqDataSetFromMatrix(countData = pasCts, colData = pasAnno, design = ~ condition) #creates a DESeqDataSet object from count data (cts) and sample metadata (coldata). The formula ~ condition specifies the experimental design, where condition is a categorical variable indicating different experimental conditions or groups.
dds$condition = relevel(dds$condition, ref = "unaffected") #sets the reference level of the condition factor to "untreated".
sample_names <- c("Sample_PNPLA7", "Sample_GMO8399", "Sample_GMO8339", "Sample_GMO8398", "Sample_GMO2978", "Sample_GMO8447", "Sample_GMO5565", "Sample_GMO3440", "F-731")
run_identifiers <- paste0("run", seq_len(nrow(dds)))
dds$run <- paste0("run",sample_names)
run_identifiers <- paste0("run", 1:length(colnames(dds)))
dds$run <- run_identifiers
ddscoll= collapseReplicates(dds, dds$sample, dds$run)
ddscoll2 = DESeq(ddscoll) #collapsed data
resCollapsed <- results(ddscoll2) #collapsed data results

dfDESeq_Ouputs_Collapsed = data.frame(resCollapsed)
write.csv(dfDESeq_Ouputs_Collapsed, file = "DESeq_Outputs_Collapsed.csv")
DESeq_Outputs_Collapsed = read.csv("DeSeq_Outputs_Collapsed.csv")

# Extract p-values from DESeq2 resultsp
p_values <- res$pvalue

# Sort p-values in ascending order
sorted_p_values <- sort(p_values)
ResOrdered <-res[order(res$padj), ]
write.table(ResOrdered, file = "P_adj.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

#code for PCA graph
library(ggplot2)
vsd <- vst(DESeq_Outputs_Collapsed, blind=FALSE)
plotPCA(vsd, intgroup=c("condition", "Cell.Line"))
pcaData <- plotPCA(vsd, intgroup=c("condition", "Cell.Line"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=Cell.Line)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + 
  ggtitle("PCA Plot by Condtion and Cell Line")

#code for log2FC(y) vs padj(x)

DESeq_Outputs_Collapsed$minuslog10padj = -(log10(DESeq_Outputs_Collapsed$padj))
#code for boxplots 
library(ggplot2) 
p = ggplot(data=DESeq_Outputs_Collapsed , aes(x = DESeq_Outputs_Collapsed$log2FoldChange, y = DESeq_Outputs_Collapsed$minuslog10padj, color=true_val)) + geom_point() #scatterplot of log2FoldChange vs padj
p = p + labs(x="log2FoldChange", y="-log10(padj)", title="-log10(padj) vs log2FoldChange Volcano Plot") + scale_color_gradient(low="blue", high="red") #labels axis 
which(DESeq_Outputs_Collapsed$minuslog10padj >= -log10(0.05), arr.ind=TRUE) #use if not ordered
x= which(DESeq_Outputs_Collapsed$minuslog10padj >= -log10(0.05), arr.ind=TRUE) # passing to variable
DESeq_Outputs_Collapsed$true_val=0 #creating a column where are values are 0
DESeq_Outputs_Collapsed$true_val[x]=1#assigning data with 1 value to whichever entries in which log10 value is above vutoff value



ggplot(DESeq_Outputs_Collapsed, aes(x="padj", y="log2FoldChange")) + geom_boxplot(outlier.colour="red", outlier.shape=8,outlier.size=4) # Change outlier, color, shape and size
p + stat_summary(fun.y=mean, geom="point", shape=23, size=4)# Box plot with mean points
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)# Box plot with dot plot
p + geom_jitter(shape=16, position=position_jitter(0.2))# Box plot with jittered points # 0.2 : degree of jitter in x direction
p<-ggplot(DESeq_Outputs, aes(x=dose, y=len, color=dose)) + geom_boxplot()# Change box plot line colors by groups
p


#code for biomart annotations
library("biomaRt") #opens biomaRt package
searchDatasets(mart = ensembl, pattern = "hsapiens")  #generate ensembl object
ensembl = useEnsembl(biomart = "genes")  #generate ensembl object
ensembl = useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)  #generate ensembl object
annotation <- getBM(filters = 'ensembl_gene_id', attributes=c('hgnc_symbol','gene_biotype','entrezgene_id','ensembl_gene_id'), values = DESeq_Outputs$genes, mart = ensembl) #running annotations
New_Row_Order <- match(DESeq_Outputs$genes, annotation$ensembl_gene_id) #moves our row order of our second data frame column to be in the row order of the first data frame column
Annotation_Rearranged <- annotation[New_Row_Order, , drop=FALSE] #create new Dataframe with new order
DESeq_Outputs$HGNC=Annotation_Rearranged$hgnc_symbol #upending the newly ordered column from the second dataframe to the first original data frame
DESeq_Outputs$Entrez=Annotation_Rearranged$entrezgene_id

library(biomaRt)

# Specify the Ensembl dataset
ensembl <- useMart("ensembl")

# Now, you can run your getBM() function with the selected dataset
annotation <- getBM(filters = 'ensembl_gene_id', 
                    attributes = c('hgnc_symbol', 'gene_biotype', 'entrezgene_id', 'ensembl_gene_id'), 
                    values = DESeq_Outputs_Collapsed$genes, 
                    mart = ensembl)

#GSEA Plots

#KEGG Bar Graph
library(ggplot2)# Load required libraries
df_GSEA_KEGG <- read.csv("GSEA_result_KEGG.tsv", sep="\t")# dataframe creation
df_GSEA_KEGG$Description <- factor(df_GSEA_KEGG$Description, levels = df_GSEA_KEGG$Description[order(-df_GSEA_KEGG$p.adjust)]) # Reorder the ID factor variable based on p.adj values in decreasing order
ggplot(df_GSEA_KEGG, aes(x = p.adjust, y = Description)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "GSEA KEGG Results Bar Graph", x = "p.adjust", y = "KEGG Term") +
  theme_minimal()# Bar plot


#REACTOME Bar Graph
library(ggplot2)# Load required libraries
df_GSEA_REACTOME_15 <- read.csv("GSEA_result_REACTOME15.csv")# dataframe creation
df_GSEA_REACTOME_15$setSize <- factor(df_GSEA_REACTOME_15$setSize, levels = df_GSEA_REACTOME_15$setSize[order(-df_GSEA_KEGG_15$p.adjust)])
ggplot(df_GSEA_REACTOME_15, aes(x = p.adjust, y = setSize)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  labs(title = "GSEA REACTOME Results Bar Graph", x = "p.adjust", y = "REACTOME Term") +
  theme_minimal()# Bar plot


#KEGG Expanding dot plot
library(ggplot2)
df_GSEA_KEGG <- read.csv("GSEA_result_KEGG.tsv", sep="\t")# dataframe creation
ggplot(df_GSEA_KEGG, aes(x = enrichmentScore, y = Description, size = p.adjust)) +
  geom_point(color = "purple") +
  scale_size_continuous(range = c(3, 10)) +  # Adjust the range of dot sizes according to your preference
  labs(title = "GSEA KEGG Results",
       x = "Enrichment Score",
       y = "KEGG Term",
       size = "Adjusted P value") +
  theme_minimal()

#REACTOME Expanding Dot Plot
library(ggplot2)
df_GSEA_REACTOME_15 <- read.csv("GSEA_result_REACTOME15.csv")
ggplot(df_GSEA_REACTOME_15, aes(x = enrichmentScore, y = setSize, size = p.adjust)) +
  geom_point(color = "darkmagenta") +
  scale_size_continuous(range = c(3, 10)) +  # Adjust the range of dot sizes according to your preference
  labs(title = "GSEA REACTOME Results",
       x = "Enrichment Score",
       y = "REACTOME Term",
       size = "Adjusted P value") +
  theme_minimal()

#GO Expanding dot plot
library(ggplot2)
df_GSEA_GO_15 <- read.csv("GSEA_results_GO15.csv")# dataframe creation
ggplot(df_GSEA_GO_15, aes(x = enrichmentScore, y = setSize, size = p.adjust)) +
  geom_point(color = "hotpink") +
  scale_size_continuous(range = c(3, 10)) +  # Adjust the range of dot sizes according to your preference
  labs(title = "GSEA GO Results",
       x = "Enrichment Score",
       y = "GO Term",
       size = "Adjusted P value") +
  theme_minimal()

#GO Bar Graph
library(ggplot2)# Load required libraries
df_GSEA_GO_15 <- read.csv("GSEA_results_GO15.csv")# dataframe creation
df_GSEA_GO_15$setSize <- factor(df_GSEA_GO_15$setSize, levels = df_GSEA_GO_15$setSize[order(-df_GSEA_GO_15$p.adjust)])
ggplot(df_GSEA_REACTOME_15, aes(x = p.adjust, y = setSize)) +
  geom_bar(stat = "identity", fill = "turquoise") +
  labs(title = "GSEA GO Results", x = "p.adjust", y = "GO Term") +
  theme_minimal()# Bar plot

DESeq_Outputs = read.csv("DESeq_Outputs.csv") #to read DESeq_Outputs back in

library(clusterProfiler)# Using gene ontology
goGO <- enrichGO(gene = DESeq_Outputs$Entrez[1:16], OrgDb = "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "fdr", pvalueCutoff = 0.05) 
write.csv(goGO, file = "Overrepresentation_GO.csv", row.names = FALSE)
goKEGG <- enrichKEGG(gene=DESeq_Outputs$Entrez[1:16], organism="hsa", keyType = "kegg", pvalueCutoff = 0.05) #Using KEGG database go
write.csv(goKEGG, file = "Overrepresentation_KEGG.csv", row.names = FALSE)
goREACTOME <- enrichPC(gene=DESeq_Outputs$HGNC[1:16], source='reactome', keyType = "hgnc", pvalueCutoff = 0.05) # Using REACTOME database
write.csv(goREACTOME, file = "Overrepresentation_REACTOME.csv", row.names = FALSE)

cluster_summary = data.frame(goGO)
clusterKEGG_summary = data.frame(goKEGG)
clusterREACTOME_summary = data.frame(goREACTOME)


#0.1 p value testing
goGO1MF <- enrichGO(gene = DESeq_Outputs$Entrez[1:55], OrgDb = "org.Hs.eg.db", ont = "MF", pAdjustMethod = "fdr", pvalueCutoff = 0.1) 
goGO1BP <- enrichGO(gene = DESeq_Outputs$Entrez[1:55], OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.1) 
goGO1CC <- enrichGO(gene = DESeq_Outputs$Entrez[1:55], OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.1) 
goKEGG1 <- enrichKEGG(gene=DESeq_Outputs$Entrez[1:55], organism="hsa", keyType = "kegg", pvalueCutoff = 0.1) #Using KEGG database go
goREACTOME1 <- enrichPC(gene=DESeq_Outputs$HGNC[1:55], source='reactome', keyType = "hgnc", pvalueCutoff = 0.1) # Using REACTOME database

cluster_summary1MF = data.frame(goGO1MF)
cluster_summary1BP = data.frame(goGO1BP)
cluster_summary1CC = data.frame(goGO1CC)
clusterKEGG_summary1 = data.frame(goKEGG1)
clusterREACTOME_summary1 = data.frame(goREACTOME1)

DESeq_Outputs$L2FC=sign(DESeq_Outputs$log2FoldChange)
DESeq_Outputs$minuslog10padj = -(log10(DESeq_Outputs$padj))
DESeq_Outputs$RankedL2FC = DESeq_Outputs$minuslog10padj/DESeq_Outputs$L2FC
y = DESeq_Outputs[,c("genes", "RankedL2FC")]
write.table(y,file="PNPLA7_DifferentialExpression.rnk",quote=F,sep="/t",row.names=F)


                            
GSEA_prep = read.csv("DESeq_Outputs.csv", header=TRUE)
gene_list_original <- GSEA_prep$log2FoldChange
names(gene_list_original) <- GSEA_prep$genes
gene_list<-na.omit(gene_list_original)
gene_list = sort(gene_list, decreasing = TRUE)
library(org.Hs.eg.db)
gseGO <- gseGO(geneList = gene_list,
                        ont ="ALL",
               keyType = "ENSEMBL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "none")
write.table(gseGO@result, file = "GSEA_results.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

# Convert gene IDs for gseKEGG function

ids<-bitr(names(gene_list_original), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)# We will lose some genes here because not all IDs will be converted
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
gseKEGG = GSEA_prep[GSEA_prep$genes %in% dedup_ids$ENSEMBL,]# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
gseKEGG$entrez = dedup_ids$ENTREZID# Create a new column in df2 with the corresponding ENTREZ IDs
kegg_gene_list <- gseKEGG$log2FoldChange# Create a vector of the gene unuiverse
names(kegg_gene_list) <- gseKEGG$entrez# Name vector with ENTREZ ids
kegg_gene_list<-na.omit(kegg_gene_list)# omit any NA values 
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)# sort the list in decreasing order (required for clusterProfiler)
KEGG_gsea <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
write.table(KEGG_gsea@result, file = "GSEA_result_KEGG.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

library(ReactomePA)
REACTOME_gsea <- gsePathway(kegg_gene_list, 
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = FALSE)
write.table(REACTOME_gsea@result, file = "GSEA_result_REACTOME.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

read.csv("GSEA_result_KEGG.tsv", sep="\t")
