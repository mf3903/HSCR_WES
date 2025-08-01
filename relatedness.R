#use Ming's modified algorithm derived from yling

gt=read.table('gt.GT.FORMAT.simple', sep = '\t', header = T, na.strings = c('./.', '.'), stringsAsFactors = F)

#check the uniq gt of all individuals
unique(as.vector(as.matrix(subset(gt, select=-c(CHROM,POS)))))

#create  dataframe of pairwise ID
indID=colnames(gt)[3:ncol(gt)]
SampleID1=rep(indID,each=length(indID))
SampleID2 = rep(indID, length(indID))
samplepair = data.frame(cbind(SampleID1, SampleID2))
samplepair$relation = ifelse(samplepair$SampleID1==samplepair$SampleID2, 'self', NA)

#load similarity function (derived from yling's script on hpc)
genotypes=subset(gt, select=-c(CHROM,POS))
individuals=samplepair
ifconditional='no'

  geno<-genotypes
  ind<-individuals
  ifc<-ifconditional # if "no" uses all common markers if "yes" calculates conditional similarity
  
 
  #initialise all possible gt
  ind$"0_2"=NA
  ind$"0_1"=NA
  ind$"0_0"=NA
  ind$"1_2"=NA
  ind$"1_1"=NA
  ind$"1_0"=NA
  ind$"2_2"=NA
  ind$"2_1"=NA
  ind$"2_0"=NA
  
  #initialize pairwise similarity parameters
  ind$n_x_tot=NA
  ind$n_y_tot=NA
  ind$n_xy_tot=NA
  ind$S=NA


 
  for (i in 1:nrow(ind)){
    # grep sample IDs from genotype matrix
    n1<-grep(paste0(ind[i,1],"$"), colnames(geno))
    n2<-grep(paste0(ind[i,2],"$"), colnames(geno))
    
    # display all possibilities of geno combination in each pair of relatives
    ind$"0_2"[i]<-length(which(geno[,n1]==0 & geno[,n2]==2))
    ind$"0_1"[i]<-length(which(geno[,n1]==0 & geno[,n2]==1))
    ind$"0_0"[i]<-length(which(geno[,n1]==0 & geno[,n2]==0))
    ind$"1_2"[i]<-length(which(geno[,n1]==1 & geno[,n2]==2))
    ind$"1_1"[i]<-length(which(geno[,n1]==1 & geno[,n2]==1))
    ind$"1_0"[i]<-length(which(geno[,n1]==1 & geno[,n2]==0))
    ind$"2_2"[i]<-length(which(geno[,n1]==2 & geno[,n2]==2))
    ind$"2_1"[i]<-length(which(geno[,n1]==2 & geno[,n2]==1))
    ind$"2_0"[i]<-length(which(geno[,n1]==2 & geno[,n2]==0))   

    # numerator of similarity
    ind$n_x_tot[i]=sum(ind$"0_0"[i], ind$"2_2"[i], ind$"0_2"[i],  ind$"2_0"[i], ind$"0_1"[i], ind$"2_1"[i], 2*ind$"1_0"[i], 2*ind$"1_2"[i], 2*ind$"1_1"[i], na.rm = T)
    ind$n_y_tot[i] = sum(ind$"0_0"[i], ind$"2_2"[i], ind$"0_2"[i],  ind$"2_0"[i], 2*ind$"0_1"[i], 2*ind$"2_1"[i], ind$"1_0"[i], ind$"1_2"[i], 2*ind$"1_1"[i], na.rm = T)
    ind$n_xy_tot[i] = sum(ind$"0_0"[i], ind$"2_2"[i], ind$"0_2"[i],  ind$"2_0"[i], ind$"0_1"[i], ind$"2_1"[i], ind$"1_0"[i], ind$"1_2"[i], 2*ind$"1_1"[i], na.rm = T)
    
    ind$S[i] = 0.5*ind$n_xy_tot[i] * (1/ind$n_x_tot[i] + 1/ind$n_y_tot[i])
      }



# S = 0 whenever nx, ny or nxy is zero
ind$S[ind$n_x_tot==0]=0
ind$S[ind$n_y_tot==0]=0
ind$S[ind$n_xy_tot==0]=0


write.table(ind, "ind_pairwise_similarity.txt", row.names = F, quote = F, sep = "\t")
