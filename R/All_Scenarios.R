#' #library
## ------------------------------------------------------------------------------------------------------------
options(java.parameters = '-Xmx8000m')
library(dplyr)
library(xlsx)
library(readxl)
library(rJava)
library("viper")
library(gridExtra)
library(ggplot2)
library(grid)
library(pheatmap)
library("car")
library(RColorBrewer)
library(gdata)
library(readr)
library(vioplot)
#-----read data-------------------------
## RC proteins-->RC_c, RC_N---------------------------------
RC<- read.delim("./proteins_rc.txt", header= TRUE, sep="\t")
RC_C <- select(RC, starts_with("C"))
RC_C[is.na(RC_C)]=0

RC_N <- select(RC, starts_with("N"))
RC_N[is.na(RC_N)]=0

## RC peptides-->RC_peptides_uniq-----------------------------------
RC_peptides2= read.delim("./peptides_rc.txt")
RC_peptides <- read.delim("./peptides_rc.txt")
rownames(RC_peptides) <- make.names(RC_peptides[,1],unique=T)
RC_peptides=RC_peptides[,5:ncol(RC_peptides)]
RC_peptides[is.na(RC_peptides)]=0

tmpname= unique(RC_peptides2[,1])
RC_peptides_uniq=matrix(data=0,nrow = length(tmpname),ncol = ncol(RC_peptides))
rownames(RC_peptides_uniq) <- tmpname
colnames(RC_peptides_uniq) <- colnames(RC_peptides)
for (j in 1:ncol(RC_peptides))
{
  for(i in 1:nrow(RC_peptides))
  {
    splitm=rownames(RC_peptides[i,])
    splits=strsplit(splitm,".",fixed=TRUE)[[1]] # remove the decimal points in protein names
    tmpscore=RC_peptides[i,j]
    if(as.numeric(tmpscore)>0)
    {
      target=which((rownames(RC_peptides_uniq) %in% splits ))
      RC_peptides_uniq[target,j]=1
    }
  }
}
## CCLE dataset-->RC_tpu,Liver_tpu------
sce_e_get_ccle <- function(data_path, tissue_type)
{
  CCLE <- as.data.frame(read_excel(data_path, sheet = "Normalized Protein Expression", col_names =T))
  tissue_ccle <- select(CCLE, c("Gene_Symbol", "Uniprot_Acc",contains(tissue_type)))
  if (tissue_type=="KIDNEY"){
    n_col = 13
  }
  else{
    n_col = 16
  }
  
  #cancer proteins
  tissue_ccle[tissue_ccle[,]=="NA"]=NA
  #for proteins appeared in all samples
  tumor_prot=c()
  for(i in 1:nrow(tissue_ccle))
  {
    l=length(which(is.na(tissue_ccle[i,3:n_col]) %in% "TRUE"))
    if(l<1) #value in all samples should not be zero
    {
      tumor_prot=rbind(tumor_prot,tissue_ccle[i,2])
    }
  }
  
  for (i in 1:nrow(tumor_prot)) {
    tumor_prot[i] = strsplit(tumor_prot[i],"-")[[1]][1]
  }
  tumor_prot_uniq = unique(tumor_prot)
  length(tumor_prot_uniq)
  tpu = tumor_prot_uniq
  return (tpu)
}
RC_tpu <- sce_e_get_ccle("./ccle.xlsx","KIDNEY")
Liver_tpu <- sce_e_get_ccle("./ccle.xlsx","LIVER")


## hela siha DDA proteins-->Hela_ddaProteins,Siha_ddaProteins---------------------------------------------------------------
Hela_ddaProteins <- read.csv("./proteins_hela_dda.csv")
Siha_ddaProteins <- read.csv("./proteins_siha_dda.csv")

colnames(Hela_ddaProteins)=c("Proteins","811","812","813")
rownames(Hela_ddaProteins)=Hela_ddaProteins[,1]
Hela_ddaProteins=subset(Hela_ddaProteins, select=-Proteins)

colnames(Siha_ddaProteins)=c("Proteins","811","812","813")
rownames(Siha_ddaProteins)=Siha_ddaProteins[,1]
Siha_ddaProteins=subset(Siha_ddaProteins, select=-Proteins)


## hela siha DDA peptides--> heladdapepuniq, sihaddapepuniq---------------------------------------------------------------
Hela_ddapep <- read.csv("./peptides_hela_dda.csv",header=TRUE)

heladdapep=c()
heladdapep1=c()
heladdapep2=c()
heladdapep3=c()
#Filter peptides
for (rowIndex in 1:nrow(Hela_ddapep))
{
  s = Hela_ddapep[rowIndex, 18] # rowIndex, 18 or 14
  splitt = strsplit(s, ":", fixed = TRUE)
  splitm=unlist(splitt)
  for(i in 1:length(splitm))
  {
    splits=strsplit(splitm[i],"|",fixed=TRUE)[1]
    splitS=splits[[1]]
    if (length(splitS) >= 1) {
      subS = splitS[1]
      if(!is.na(subS))
      {
        if(!(subS %in% heladdapep))
          heladdapep=rbind(heladdapep,subS)
        if (!(subS %in% heladdapep1) && Hela_ddapep[rowIndex,1]=="1") { # rowIndex, 1 or 11
          heladdapep1 = rbind(heladdapep1, subS)
        }
        if (!(subS %in% heladdapep2) && Hela_ddapep[rowIndex,1]=="2") {
          heladdapep2 = rbind(heladdapep2, subS)
        }
        if (!(subS %in% heladdapep3) && Hela_ddapep[rowIndex,1]=="3") {
          heladdapep3 = rbind(heladdapep3, subS)
        }
      }
    }
  }
}
heladdapepuniq=matrix(data=0,nrow = length(heladdapep),ncol = 3)
rownames(heladdapepuniq)=heladdapep

for(i in 1:nrow(heladdapep))
{
  tmp=heladdapep[i]
  if(tmp %in% heladdapep1) heladdapepuniq[i,1]=1
  if(tmp %in% heladdapep2) heladdapepuniq[i,2]=1
  if(tmp %in% heladdapep3) heladdapepuniq[i,3]=1
}
heladdapepuniq=as.data.frame(heladdapepuniq)
colnames(heladdapepuniq)=c("sample1","sample2","sample3")
--

Siha_ddapep <- read.csv("./peptides_siha_dda.csv",header=TRUE)
sihaddapep=c()
sihaddapep1=c()
sihaddapep2=c()
sihaddapep3=c()
#Filter peptides
for (rowIndex in 1:nrow(Siha_ddapep))
{
  s = Siha_ddapep[rowIndex, 18]
  splitt = strsplit(s, ":", fixed = TRUE)
  splitm=unlist(splitt)
  for(i in 1:length(splitm))
  {
    splits=strsplit(splitm[i],"|",fixed=TRUE)[1]
    splitS=splits[[1]]
    if (length(splitS) >= 1) {
      subS = splitS[1]
      if(!is.na(subS))
      {
        if(!(subS %in% sihaddapep))
          sihaddapep=rbind(sihaddapep,subS)
        if (!(subS %in% sihaddapep1) && Siha_ddapep[rowIndex,1]=="1") {
          sihaddapep1 = rbind(sihaddapep1, subS)
        }
        if (!(subS %in% sihaddapep2) && Siha_ddapep[rowIndex,1]=="2") {
          sihaddapep2 = rbind(sihaddapep2, subS)
        }
        if (!(subS %in% sihaddapep3) && Siha_ddapep[rowIndex,1]=="3") {
          sihaddapep3 = rbind(sihaddapep3, subS)
        }
      }
    }
  }
}
sihaddapepuniq=matrix(data=0,nrow = length(sihaddapep),ncol = 3)
rownames(sihaddapepuniq)=sihaddapep

for(i in 1:nrow(sihaddapep))
{
  tmp=sihaddapep[i]
  if(tmp %in% sihaddapep1) sihaddapepuniq[i,1]=1
  if(tmp %in% sihaddapep2) sihaddapepuniq[i,2]=1
  if(tmp %in% sihaddapep3) sihaddapepuniq[i,3]=1
}
sihaddapepuniq=as.data.frame(sihaddapepuniq)
colnames(sihaddapepuniq)=c("sample1","sample2","sample3")


## hela siha DIA proteins-->hela_dia_prot_list1, hela_dia_prot_list2, siha_dia_prot_list1, siha_dia_prot_list2, siha_dia_prot_list3-----
##DIA
Hela_diaProteins <- read.csv("./proteins_hela_dia.csv")
Siha_diaProteins <- read.csv("./proteins_siha_dia.csv")

colnames(Hela_diaProteins)=c("Proteins","first","second")

hela_dia_prot_list1 <- as.data.frame(unique(Hela_diaProteins[which(as.numeric(Hela_diaProteins$first)>0),1]))
colnames(hela_dia_prot_list1)='proteins'

hela_dia_prot_list2 <- as.data.frame(unique(Hela_diaProteins[which(as.numeric(Hela_diaProteins$second)>0),1]))
colnames(hela_dia_prot_list2)='proteins'

colnames(Siha_diaProteins)=c("Proteins","first","second","third")

siha_dia_prot_list1 <- as.data.frame(unique(Siha_diaProteins[which(as.numeric(Siha_diaProteins$first)>0),1]))
colnames(siha_dia_prot_list1)='proteins'

siha_dia_prot_list2 <- as.data.frame(unique(Siha_diaProteins[which(as.numeric(Siha_diaProteins$second)>0),1]))
colnames(siha_dia_prot_list2)='proteins'

siha_dia_prot_list3 <- as.data.frame(unique(Siha_diaProteins[which(as.numeric(Siha_diaProteins$third)>0),1]))
colnames(siha_dia_prot_list3)='proteins'

## LC proteins-->LC_T------
#LC_T
#protein list
LC <- as.data.frame(read_excel("./proteins_lc.xlsx",sheet = "c_proteins"))
rownames(LC) <- LC[,1]
LC <- LC[,-1]

LC_T <- select(LC, contains("T"))
LC_T[LC_T[,]==1]=0
LC_T <- as.data.frame(LC_T)

## HCC proteins-->HCC------
#HCC
for(i in 1:5){
  name=paste0("Rep", as.character(i))
  test <- as.data.frame(read_excel("./proteins_hcc.xls",sheet = paste0(name)))[c(1,3:8)]
  colnames(test)=c("Protein", test[2,2:7])
  test=test[3:nrow(test),]
  test = test %>%
    dplyr::group_by(Protein) %>%
    dplyr::mutate(Protein=substring(Protein, 4, 9))
  assign(paste0("Rep", i), test)
  assign(paste0(name, "_prot"), unique(test[,1]))
}

df1 <- merge(x=Rep1[,1:4], y=Rep2[,1:4], by="Protein",all=TRUE, suffixes = c("_P1", "_P2"))
df2 <- merge(x=Rep3[,1:4], y=Rep4[,1:4], by="Protein",all=TRUE, suffixes = c("_P3", "_P4"))
df3 <- merge(x=df1, y=df2, by="Protein",all=TRUE)
df4 <- unique(merge(x=df3, y=Rep5[,1:4], by="Protein",all=TRUE))
HCC <- df4[,2:16]
colnames(HCC)=c("P1_HCC1", "P1_HCC2", "P1_HCC3", "P2_HCC1", "P2_HCC2", "P2_HCC3","P3_HCC1", "P3_HCC2", "P3_HCC3", "P4_HCC1", "P4_HCC2", "P4_HCC3", "P5_HCC1", "P5_HCC2", "P5_HCC3")
rownames(HCC)=make.names(df4[,1], unique = T)
HCC[is.na(HCC)]=0
HCC=HCC[,colnames(HCC)[order(colnames(HCC))]]
HCC <- data.frame(sapply(HCC, function(x) as.numeric(as.character(x))))
rownames(HCC)=make.names(df4[,1], unique = T)



## complexes-->complex_vector2018, complex_vector2022------------------------------------------------------------------------------------------------------------
human_complexes_2018 <- read.delim("./human_complexes_2018.txt")
human_complexes_2018=human_complexes_2018[human_complexes_2018[,3]=="Human",]
complex_vector2018 <- strsplit(as.vector(human_complexes_2018[,6]), ';')
complex_vector2018 <- setNames(complex_vector2018, human_complexes_2018[,1])

human_complexes_2022 <- read.delim("./human_complexes_2022.txt")
human_complexes_2022=human_complexes_2022[human_complexes_2022[,3]=="Human",]
complex_vector2022 <- strsplit(as.vector(human_complexes_2022[,6]), ';')
complex_vector2022 <- setNames(complex_vector2022, human_complexes_2022[,1])

year = '2018'
if(year=='2018')
{
  complex_vector<- complex_vector2018}else{
    complex_vector<- complex_vector2022
  }
## FCS filtering-->fcs_normal_mat,fcs_cancer_mat,protprob_normal_mat,protprob_cancer_mat,FCS_sig_complex_n,FCS_sig_complex_c--------------
# fcs results, rows: proteins, columns: normal samples, values: probability of proteins being in samples
fcs_normal_mat <- read.table("./normal_protein_probabilities.txt", sep="\t", header=T)
fcs_cancer_mat <- read.table("./cancer_protein_probabilities.txt", sep="\t", header=T)

# protrec results, rows: proteins, columns: normal samples, values: probability of proteins being in samples
normal_mat <- read.table("./from_protprob_normal_protein_probabilities.txt", sep="\t", header=T)
cancer_mat <- read.table("./from_protprob_cancer_protein_probabilities.txt", sep="\t", header=T)

#the proteins associated with significant complexes in FCS
# rows: samples, columns: complexes, values: 1 - (probability of complexes present in samples)
FCS_sig_complex_n <- read.table("./FCS_output_RC_N.txt", header=T, check.names=F)
FCS_sig_complex_c <- read.table("./FCS_output_RC_C.txt", header=T, check.names=F)

#-----Protrec methods:------------------------------------------------------------------------------------------------------------
PROTREC_cplx_prob <- function(data, complex, fdr, threshold)
{
  output_mat <- c()
  for (x in 1:ncol(data))
  {
    p_val <- c()
    for (i in 1:length(complex))
    {
      size_complex <- length(complex[[i]])
      
      if(size_complex < threshold)
      {
        p_val <- append(p_val, 1)
      }
      else
      {
        detected_prot <- rownames(data)[as.numeric(data[,x])!=0]
        prob_complex_exists <- (length(intersect(detected_prot, complex[[i]]))/size_complex) * (1 - fdr)
        prob_complex_exists_pval <- (1 - prob_complex_exists)
        
        p_val <- append(p_val, prob_complex_exists_pval)
      }
    }
    output_mat <- rbind(output_mat, p_val)
  }
  colnames(output_mat) <- names(complex)
  rownames(output_mat) <- colnames(data)
  return(output_mat)
}

PROTREC_protprob <- function(cplx, p, prot, fdr,mode)
{
  new_mat <- c()
  for (i in 1:length(cplx))
  { 
    protein_probabilities <- c() # protein_probabilities is a list containing the conditional probabilities of the constituent proteins in this complex
    size_complex <- length(cplx[[i]]) # number of proteins that make up this complex
    complex_names_vec <- rep(names(cplx)[i], size_complex)
    complex_prot_names_vec <- cplx[[i]]
    
    for(j in 1:length(cplx[[i]])) 
    {
      if (cplx[[i]][j] %in% prot) #this protein is found in the reported list then its probability is dependent on whether the complex is formed or not
      { 
        protein_probabilities <- append(protein_probabilities, (p[names(p) %in% names(cplx)[i]]) + (1 - fdr)*(1 - p[names(p) %in% names(cplx)[i]])) 
      }
      else #this is an additional protein and it takes on the probability of the complex is formed only
      {
        protein_probabilities <- append(protein_probabilities, p[names(p) %in% names(cplx)[i]]) #simply the probability the complex is  formed
      }
    }
    
    temp_mat <- cbind(complex_prot_names_vec, complex_names_vec, round(protein_probabilities, 4))
    new_mat <- rbind(new_mat, temp_mat)
  }
  
  prot_list <- unique(sort(unlist(cplx))) # total number of unique proteins in the complex database = 3674
  output <- c()
  for (i in 1:length(prot_list))
  {
    #max
    if(mode=='max'){output<- rbind(output,cbind(prot_list[i], max(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]),3]))))}
    #min
    else if(mode=='min'){output<- rbind(output,cbind(prot_list[i], min(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]),3]))))}
    #median
    else if(mode=='median'){output<- rbind(output,cbind(prot_list[i], median(unlist(as.numeric(new_mat[which(new_mat[,1]%in%prot_list[i]),3])))))}
    #mean
    else if(mode=='mean'){output<- rbind(output,cbind(prot_list[i], mean(unlist(as.numeric(new_mat[which(new_mat[,1]%in%prot_list[i]),3])))))}
  }
  return(output)
}

#-----complex_prob------------------------------------------------------------
#RC_N
rc_nprotrec <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, 5)
rc_nprotrec10 <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, 10)
rc_nprotrec20 <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, 20)
rc_nprotrec30 <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, 30)
rc_nprotrec40 <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, 40)
rc_nprotrec50 <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, 50)
rc_nprotrec200 <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, 200)

#RC_C
rc_cprotrec <- PROTREC_cplx_prob(RC_C, complex_vector, 0.01, 5)
rc_cprotrec10 <- PROTREC_cplx_prob(RC_C, complex_vector, 0.01, 10)
rc_cprotrec20 <- PROTREC_cplx_prob(RC_C, complex_vector, 0.01, 20)
rc_cprotrec30 <- PROTREC_cplx_prob(RC_C, complex_vector, 0.01, 30)
rc_cprotrec40 <- PROTREC_cplx_prob(RC_C, complex_vector, 0.01, 40)
rc_cprotrec50 <- PROTREC_cplx_prob(RC_C, complex_vector, 0.01, 50)
rc_cprotrec200 <- PROTREC_cplx_prob(RC_C, complex_vector, 0.01, 200)

#hela
hela_ddaprotrec <- PROTREC_cplx_prob(Hela_ddaProteins, complex_vector, 0.01, 5)
hela_ddaprotrec10 <- PROTREC_cplx_prob(Hela_ddaProteins, complex_vector, 0.01, 10)
hela_ddaprotrec20 <- PROTREC_cplx_prob(Hela_ddaProteins, complex_vector, 0.01, 20)
hela_ddaprotrec30 <- PROTREC_cplx_prob(Hela_ddaProteins, complex_vector, 0.01, 30)
hela_ddaprotrec40 <- PROTREC_cplx_prob(Hela_ddaProteins, complex_vector, 0.01, 40)
hela_ddaprotrec50 <- PROTREC_cplx_prob(Hela_ddaProteins, complex_vector, 0.01, 50)
hela_ddaprotrec200 <- PROTREC_cplx_prob(Hela_ddaProteins, complex_vector, 0.01, 200)

#siha
siha_ddaprotrec <- PROTREC_cplx_prob(Siha_ddaProteins, complex_vector, 0.01, 5)
siha_ddaprotrec10 <- PROTREC_cplx_prob(Siha_ddaProteins, complex_vector, 0.01, 10)
siha_ddaprotrec20 <- PROTREC_cplx_prob(Siha_ddaProteins, complex_vector, 0.01, 20)
siha_ddaprotrec30 <- PROTREC_cplx_prob(Siha_ddaProteins, complex_vector, 0.01, 30)
siha_ddaprotrec40 <- PROTREC_cplx_prob(Siha_ddaProteins, complex_vector, 0.01, 40)
siha_ddaprotrec50 <- PROTREC_cplx_prob(Siha_ddaProteins, complex_vector, 0.01, 50)
siha_ddaprotrec200 <- PROTREC_cplx_prob(Siha_ddaProteins, complex_vector, 0.01, 200)

#LC
#LC_T
lc_tprotrec <- PROTREC_cplx_prob(LC_T, complex_vector, 0.01, 5)
lc_tprotrec10 <- PROTREC_cplx_prob(LC_T, complex_vector, 0.01, 10)
lc_tprotrec20 <- PROTREC_cplx_prob(LC_T, complex_vector, 0.01, 20)
lc_tprotrec30 <- PROTREC_cplx_prob(LC_T, complex_vector, 0.01, 30)
lc_tprotrec40 <- PROTREC_cplx_prob(LC_T, complex_vector, 0.01, 40)
lc_tprotrec50 <- PROTREC_cplx_prob(LC_T, complex_vector, 0.01, 50)
lc_tprotrec200 <- PROTREC_cplx_prob(LC_T, complex_vector, 0.01, 200)

#HCC
hcc_protrec <- PROTREC_cplx_prob(HCC, complex_vector, 0.01, 5)
hcc_protrec10 <- PROTREC_cplx_prob(HCC, complex_vector, 0.01, 10)
hcc_protrec20 <- PROTREC_cplx_prob(HCC, complex_vector, 0.01, 20)
hcc_protrec30 <- PROTREC_cplx_prob(HCC, complex_vector, 0.01, 30)
hcc_protrec40 <- PROTREC_cplx_prob(HCC, complex_vector, 0.01, 40)
hcc_protrec50 <- PROTREC_cplx_prob(HCC, complex_vector, 0.01, 50)
hcc_protrec200 <- PROTREC_cplx_prob(HCC, complex_vector, 0.01, 200)




#-----common functions-----------
## sce acde recovery function------------------------------------------------------------------------------------------------------------
pairwise_recovery_protrec <- function(prot_predict_list, original_prot_list, check_prot_list, complex_vec, protrecscoreset)
{
  #extract proteins with significant probabilities >= protrecscoreset (p-value- BAyes conversion)
  
  
  prot_predict_list <- prot_predict_list[which(as.numeric(prot_predict_list[,2]) >= protrecscoreset),]
  
  
  add_proteins <- setdiff(prot_predict_list[,1], original_prot_list) #this is the set of recovered proteins associated with the significant complexes
  
  obs_intersect <- round(length(intersect(add_proteins, check_prot_list))/length( add_proteins), 3)
  
  theoretical_vec <- c()
  
  for (i in 1:1000)
  {
    S_prime <- sample(unlist(complex_vec), length(add_proteins)) 
    add_proteins_prime <- setdiff(S_prime, original_prot_list)
    theoretical_intersect <- round(length(intersect(add_proteins_prime, check_prot_list))/length(add_proteins_prime), 3)
    theoretical_vec <- append(theoretical_vec, theoretical_intersect)
  }  
  
  p_value <- sum(as.numeric(theoretical_vec >= obs_intersect))/length(theoretical_vec)
  
  
  return(c(obs_intersect, p_value, length(add_proteins), length(intersect(add_proteins, check_prot_list)), add_proteins,"a",intersect(add_proteins, check_prot_list)))
}
#-----sce a-------------------------
## sce a RC ---------------------------------------
get_result_sce_a_rc <- function(rc_nprotrec,rc_cprotrec,protrecscoreset,mode){
  RC_recov_protrec_NT1NT2 <-c()
  for (i in c(1:6))  #loop for NT1 -> NT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_nprotrec[i,], rownames(RC_N),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector,protrecscoreset)
  RC_recov_protrec_NT1NT2<- rbind(RC_recov_protrec_NT1NT2,val[1:4])
  print(i)
  }
  
  RC_recov_protrec_NT2NT1 <- c()
  for (j in c(7:12))  #loop for NT2 -> NT2_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_nprotrec[j,], rownames(RC_N),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_N)[(RC_N[,j])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,j-6])!=0],complex_vector,protrecscoreset)
  RC_recov_protrec_NT2NT1<- rbind(RC_recov_protrec_NT2NT1,val[1:4])
  print(j)
  }  
  
  RC_recov_protrec_CT1CT2 <- c()
  for (k in c(1:6))  #loop for CT1 -> CT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[k,], rownames(RC_C),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_C)[(RC_C[,k])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,k+18])!=0],complex_vector,protrecscoreset)
  RC_recov_protrec_CT1CT2<- rbind(RC_recov_protrec_CT1CT2,val[1:4])
  print(k)
  }  
  
  RC_recov_protrec_CT2CT1 <- c()
  for (l in c(7:12))  #loop for CT2 -> CT2_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[l,], rownames(RC_C),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_C)[(RC_C[,l])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,l])!=0],complex_vector,protrecscoreset)
  RC_recov_protrec_CT2CT1<- rbind(RC_recov_protrec_CT2CT1,val[1:4])   
  print(l)
  }  
  
  
  RC_recov_cross_protrec<- data.frame(RC_recov_protrec_NT1NT2,RC_recov_protrec_NT2NT1,RC_recov_protrec_CT1CT2,RC_recov_protrec_CT2CT1)  
  
  colnames(RC_recov_cross_protrec)<- c(rep("NT1->NT2",4),rep("NT2->NT1",4),rep("CT1->CT2",4),rep("CT2->CT1",4))
  
  return(RC_recov_cross_protrec)
}

pathh="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_july/corum2022/RC/"

for(protrecscoreset in c(0.4,0.5,0.6,0.7,0.8,0.95))
{
  for(mode in c("max","min","median","mean"))
  {
    result = get_result_sce_a_rc(rc_nprotrec,rc_cprotrec,protrecscoreset,mode)
    write.csv(result,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_5_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",5,",mode))
    
    result1 = get_result_sce_a_rc(rc_nprotrec10,rc_cprotrec10,protrecscoreset,mode)
    write.csv(result1,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_10_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",10,",mode))
    
    result2 = get_result_sce_a_rc(rc_nprotrec20,rc_cprotrec20,protrecscoreset,mode)
    write.csv(result2,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_20_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",20,",mode))
    
    result3 = get_result_sce_a_rc(rc_nprotrec30,rc_cprotrec30,protrecscoreset,mode)
    write.csv(result3,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_30_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",30,",mode))
    
    result4 = get_result_sce_a_rc(rc_nprotrec40,rc_cprotrec40,protrecscoreset,mode)
    write.csv(result4,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_40_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",40,",mode))
    
    result5 = get_result_sce_a_rc(rc_nprotrec50,rc_cprotrec50,protrecscoreset,mode)
    write.csv(result5,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_50_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",50,",mode))
    
    result6 = get_result_sce_a_rc(rc_nprotrec200,rc_cprotrec200,protrecscoreset,mode)
    write.csv(result6,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_200_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",200,",mode))
  }
}


## sce a hela siha-------------------------------------------
get_result_sce_a_hela <- function(hela_ddaprotrec,protrecscoreset,mode){
  Hela_ddacrossrecovprotrec1 <-c()
  PROTREC_heladdaprot1 <- data.frame(PROTREC_protprob(complex_vector, 1-hela_ddaprotrec[1,], rownames(Hela_ddaProteins),0.01,mode))
  rownames(PROTREC_heladdaprot1)<- PROTREC_heladdaprot1[,1]
  PROTREC_heladdaprot1[,2]<-as.numeric(as.character(PROTREC_heladdaprot1[,2]))
  PROTREC_heladdaprot1<- data.frame(PROTREC_heladdaprot1)
  
  PROTREC_heladdaprot2 <- data.frame(PROTREC_protprob(complex_vector, 1-hela_ddaprotrec[2,], rownames(Hela_ddaProteins),0.01,mode))
  rownames(PROTREC_heladdaprot2)<- PROTREC_heladdaprot2[,1]
  PROTREC_heladdaprot2[,2]<-as.numeric(as.character(PROTREC_heladdaprot2[,2]))
  PROTREC_heladdaprot2<- data.frame(PROTREC_heladdaprot2)
  
  PROTREC_heladdaprot3 <- data.frame(PROTREC_protprob(complex_vector, 1-hela_ddaprotrec[3,], rownames(Hela_ddaProteins),0.01,mode))
  rownames(PROTREC_heladdaprot3)<- PROTREC_heladdaprot3[,1]
  PROTREC_heladdaprot3[,2]<-as.numeric(as.character(PROTREC_heladdaprot3[,2]))
  PROTREC_heladdaprot3<- data.frame(PROTREC_heladdaprot3)
  
  heladdaprotrec811812val<- pairwise_recovery_protrec(PROTREC_heladdaprot1,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,1]!=0)],heladdapep2[,1],complex_vector,protrecscoreset)
  heladdaprotrec811813val<- pairwise_recovery_protrec(PROTREC_heladdaprot1,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,1]!=0)],heladdapep3[,1],complex_vector,protrecscoreset)
  Hela_ddacrossrecovprotrec1<- rbind(Hela_ddacrossrecovprotrec1,heladdaprotrec811812val[1:4])
  Hela_ddacrossrecovprotrec1<- rbind(Hela_ddacrossrecovprotrec1,heladdaprotrec811813val[1:4])
  heladdaprotrec812813val<- pairwise_recovery_protrec(PROTREC_heladdaprot2,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,2]!=0)],heladdapep3[,1],complex_vector,protrecscoreset)
  heladdaprotrec812811val<- pairwise_recovery_protrec(PROTREC_heladdaprot2,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,2]!=0)],heladdapep1[,1],complex_vector,protrecscoreset)
  Hela_ddacrossrecovprotrec1<- rbind(Hela_ddacrossrecovprotrec1,heladdaprotrec812811val[1:4])
  Hela_ddacrossrecovprotrec1<- rbind(Hela_ddacrossrecovprotrec1,heladdaprotrec812813val[1:4])
  heladdaprotrec813811val<- pairwise_recovery_protrec(PROTREC_heladdaprot3,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,3]!=0)],heladdapep1[,1],complex_vector,protrecscoreset)
  heladdaprotrec813812val<- pairwise_recovery_protrec(PROTREC_heladdaprot3,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,3]!=0)],heladdapep2[,1],complex_vector,protrecscoreset)
  Hela_ddacrossrecovprotrec1<- rbind(Hela_ddacrossrecovprotrec1,heladdaprotrec813811val[1:4])
  Hela_ddacrossrecovprotrec1<- rbind(Hela_ddacrossrecovprotrec1,heladdaprotrec813812val[1:4])
  
  colnames(Hela_ddacrossrecovprotrec1)=c(rep("Protrec",4))
  rownames(Hela_ddacrossrecovprotrec1)=c("1--2","1--3","2--1","2--3","3--1","3--2")
  return(Hela_ddacrossrecovprotrec1)
}

##siha
get_result_sce_a_siha <- function(siha_ddaprotrec,protrecscoreset,mode){
  Siha_ddacrossrecovprotrec1 <-c()
  PROTREC_sihaddaprot1 <- data.frame(PROTREC_protprob(complex_vector, 1-siha_ddaprotrec[1,], rownames(Siha_ddaProteins),0.01,mode))
  rownames(PROTREC_sihaddaprot1)<- PROTREC_sihaddaprot1[,1]
  PROTREC_sihaddaprot1[,2]<-as.numeric(as.character(PROTREC_sihaddaprot1[,2]))
  PROTREC_sihaddaprot1<- data.frame(PROTREC_sihaddaprot1)
  
  PROTREC_sihaddaprot2 <- data.frame(PROTREC_protprob(complex_vector, 1-siha_ddaprotrec[2,], rownames(Siha_ddaProteins),0.01,mode))
  rownames(PROTREC_sihaddaprot2)<- PROTREC_sihaddaprot2[,1]
  PROTREC_sihaddaprot2[,2]<-as.numeric(as.character(PROTREC_sihaddaprot2[,2]))
  PROTREC_sihaddaprot2<- data.frame(PROTREC_sihaddaprot2)
  
  PROTREC_sihaddaprot3 <- data.frame(PROTREC_protprob(complex_vector, 1-siha_ddaprotrec[3,], rownames(Siha_ddaProteins),0.01,mode))
  rownames(PROTREC_sihaddaprot3)<- PROTREC_sihaddaprot3[,1]
  PROTREC_sihaddaprot3[,2]<-as.numeric(as.character(PROTREC_sihaddaprot3[,2]))
  PROTREC_sihaddaprot3<- data.frame(PROTREC_sihaddaprot3)
  
  sihaddaprotrec811812val<- pairwise_recovery_protrec(PROTREC_sihaddaprot1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,1]!=0)],sihaddapep2[,1],complex_vector,protrecscoreset)
  sihaddaprotrec811813val<- pairwise_recovery_protrec(PROTREC_sihaddaprot1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,1]!=0)],sihaddapep3[,1],complex_vector,protrecscoreset)
  Siha_ddacrossrecovprotrec1<- rbind(Siha_ddacrossrecovprotrec1,sihaddaprotrec811812val[1:4])
  Siha_ddacrossrecovprotrec1<- rbind(Siha_ddacrossrecovprotrec1,sihaddaprotrec811813val[1:4])
  sihaddaprotrec812813val<- pairwise_recovery_protrec(PROTREC_sihaddaprot2,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,2]!=0)],sihaddapep3[,1],complex_vector,protrecscoreset)
  sihaddaprotrec812811val<- pairwise_recovery_protrec(PROTREC_sihaddaprot2,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,2]!=0)],sihaddapep1[,1],complex_vector,protrecscoreset)
  Siha_ddacrossrecovprotrec1<- rbind(Siha_ddacrossrecovprotrec1,sihaddaprotrec812811val[1:4])
  Siha_ddacrossrecovprotrec1<- rbind(Siha_ddacrossrecovprotrec1,sihaddaprotrec812813val[1:4])
  sihaddaprotrec813811val<- pairwise_recovery_protrec(PROTREC_sihaddaprot3,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,3]!=0)],sihaddapep1[,1],complex_vector,protrecscoreset)
  sihaddaprotrec813812val<- pairwise_recovery_protrec(PROTREC_sihaddaprot3,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,3]!=0)],sihaddapep2[,1],complex_vector,protrecscoreset)
  Siha_ddacrossrecovprotrec1<- rbind(Siha_ddacrossrecovprotrec1,sihaddaprotrec813811val[1:4])
  Siha_ddacrossrecovprotrec1<- rbind(Siha_ddacrossrecovprotrec1,sihaddaprotrec813812val[1:4])
  
  colnames(Siha_ddacrossrecovprotrec1)=c(rep("Protrec",4))
  rownames(Siha_ddacrossrecovprotrec1)=c("1--2","1--3","2--1","2--3","3--1","3--2")
  return(Siha_ddacrossrecovprotrec1)
}

pathh="./output_corum22"
for(protrecscoreset in c(0.4,0.5,0.6,0.7,0.8,0.95))
{
  for(mode in c("max","min","median","mean"))
  {
    resultt_h = get_result_sce_a_hela(hela_ddaprotrec,protrecscoreset,mode)
    write.csv(resultt_h,paste0(pathh,"/sce_a_",mode,"/result_hela_cross_score",protrecscoreset,"_5_",mode,".csv"))
    
    resultt_h1 = get_result_sce_a_hela(hela_ddaprotrec10,protrecscoreset,mode)
    write.csv(resultt_h1,paste0(pathh,"/sce_a_",mode,"/result_hela_cross_score",protrecscoreset,"_10_",mode,".csv"))
    
    resultt_h2 = get_result_sce_a_hela(hela_ddaprotrec20,protrecscoreset,mode)
    write.csv(resultt_h2,paste0(pathh,"/sce_a_",mode,"/result_hela_cross_score",protrecscoreset,"_20_",mode,".csv"))
    
    resultt_h3 = get_result_sce_a_hela(hela_ddaprotrec30,protrecscoreset,mode)
    write.csv(resultt_h3,paste0(pathh,"/sce_a_",mode,"/result_hela_cross_score",protrecscoreset,"_30_",mode,".csv"))
    
    resultt_h4 = get_result_sce_a_hela(hela_ddaprotrec40,protrecscoreset,mode)
    write.csv(resultt_h4,paste0(pathh,"/sce_a_",mode,"/result_hela_cross_score",protrecscoreset,"_40_",mode,".csv"))
    
    resultt_h5 = get_result_sce_a_hela(hela_ddaprotrec50,protrecscoreset,mode)
    write.csv(resultt_h5,paste0(pathh,"/sce_a_",mode,"/result_hela_cross_score",protrecscoreset,"_50_",mode,".csv"))
    
    resultt_h6 = get_result_sce_a_hela(hela_ddaprotrec200,protrecscoreset,mode)
    write.csv(resultt_h6,paste0(pathh,"/sce_a_",mode,"/result_hela_cross_score",protrecscoreset,"_200_",mode,".csv"))
    
    
    resultt_s = get_result_sce_a_siha(siha_ddaprotrec,protrecscoreset,mode)
    write.csv(resultt_s,paste0(pathh,"/sce_a_",mode,"/result_siha_cross_score",protrecscoreset,"_5_",mode,".csv"))
    
    resultt_s1 = get_result_sce_a_siha(siha_ddaprotrec10,protrecscoreset,mode)
    write.csv(resultt_s1,paste0(pathh,"/sce_a_",mode,"/result_siha_cross_score",protrecscoreset,"_10_",mode,".csv"))
    
    resultt_s2 = get_result_sce_a_siha(siha_ddaprotrec20,protrecscoreset,mode)
    write.csv(resultt_s2,paste0(pathh,"/sce_a_",mode,"/result_siha_cross_score",protrecscoreset,"_20_",mode,".csv"))
    
    resultt_s3 = get_result_sce_a_siha(siha_ddaprotrec30,protrecscoreset,mode)
    write.csv(resultt_s3,paste0(pathh,"/sce_a_",mode,"/result_siha_cross_score",protrecscoreset,"_30_",mode,".csv"))
    
    resultt_s4 = get_result_sce_a_siha(siha_ddaprotrec40,protrecscoreset,mode)
    write.csv(resultt_s4,paste0(pathh,"/sce_a_",mode,"/result_siha_cross_score",protrecscoreset,"_40_",mode,".csv"))
    
    resultt_s5 = get_result_sce_a_siha(siha_ddaprotrec50,protrecscoreset,mode)
    write.csv(resultt_s5,paste0(pathh,"/sce_a_",mode,"/result_siha_cross_score",protrecscoreset,"_50_",mode,".csv"))
    
    resultt_s6 = get_result_sce_a_siha(siha_ddaprotrec200,protrecscoreset,mode)
    write.csv(resultt_s6,paste0(pathh,"/sce_a_",mode,"/result_siha_cross_score",protrecscoreset,"_200_",mode,".csv"))
    
    print(protrecscoreset)
  }
}




#------sce b-------------------------
## sce b FCS--------------------------------
fcs <- function(data, complex_vector, sim_size, threshold)
{
  len=c()
  for(i in 1:length(complex_vector))
  {
    tp=length(complex_vector[[i]])
    if(!(tp %in% len) && tp>=5) len=append(len, tp)
  }
  testlis=list()
  for(o in 1:length(len))
  {
    test_mat=c()
    for(m in 1:sim_size)
    {
      set.seed(m+as.numeric(format(Sys.time(), "%S")))
      tmp=sample(unique(unlist(complex_vector)), len[o],replace=T)
      tmp=matrix(tmp,nrow=1)
      test_mat=rbind(test_mat,tmp)
    }
    testlis[[len[o]]]=test_mat
  }
  output_mat <- c()
  for (x in 1:ncol(data))
  {
    print(x)
    p_val <- c()
    for (i in 1:length(complex_vector))
    {
      size_complex <- length(complex_vector[[i]])
      
      if(size_complex < threshold)
      {
        p_val <- append(p_val, 1)
      }
      else
      {
        detected_prot <- rownames(data)[as.numeric(data[,x])!=0]
        obs_overlap <- length(intersect(detected_prot, complex_vector[[i]]))/size_complex
        #print(obs_overlap)
        test_intersects <- c()
        test_mt=testlis[[size_complex]]
        for (j in 1:nrow(test_mat))
        {
          test_intersects <- append(test_intersects, length(intersect(test_mt[j,], detected_prot))/ncol(test_mt))
        }
        
        p_val <- append(p_val,sum(as.numeric(test_intersects) >= as.numeric(obs_overlap))/length(test_intersects)) #this is the p-value for this complex current
      }
      #return(p_val)
    }
    output_mat <- rbind(output_mat, p_val)
  }
  colnames(output_mat) <- names(complex_vector)
  rownames(output_mat) <- colnames(data)
  return(output_mat)
}

fcs_prot_prob_assign <- function(cplx, p)
{
  complex_names_vec <- rep(names(cplx), lapply(cplx, length))
  data<- data.frame(table(complex_names_vec))
  data[,1]=as.numeric(as.character(data[,1]))
  data=data[order(data[,1]),]
  #save(data, file="D:/data.txt", ascii=TRUE)
  prot_vec <- unlist(cplx, use.names=F, recursive=F)
  #print(prot_vec)
  new_mat <- as.matrix(cbind(prot_vec, complex_names_vec, rep(p, data[,2])))
  #print(new_mat)
  #now we should have a list of proteins that we can now check individually
  prot_list <- unique(prot_vec) #the unique list of proteins
  output <- c()
  for (i in 1:length(prot_list))
  {
    output = rbind(output,cbind(prot_list[i], max(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]), 3]
    )
    )
    )
    )
    #output<- rbind(output,cbind(prot_list[i], mean(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]),3]))))
  }
  return(output)
}
## sce b HE--------------------------------
hgtest <- function(data, complex_vector, threshold)
{
  output_mat <- c()
  for (x in 1:ncol(data))
  {
    print(x)
    p_val <- c()
    for (i in 1:length(complex_vector))
    {
      size_complex <- as.numeric(length(complex_vector[[i]]))
      
      
      if(size_complex < threshold)
      {
        p_val <- append(p_val, 1)
      }
      else
      {
        detected_prot <- rownames(data)[as.numeric(data[,x])!=0] 
        obs_cplx_prot <- length(intersect(detected_prot, complex_vector[[i]])) #q
        total_prot <- as.numeric(length(union(detected_prot, unique(unlist(complex_vector, use.names=F, recursive=F))))) #N
        sample_size <- length(detected_prot) #k
        
        #phyper(q= 10-1, 10, 3458, 1468, lower.tail = FALSE)
        
        p_val <- append(p_val, phyper(q=obs_cplx_prot -1, m=size_complex, n=(total_prot- size_complex), k=sample_size, lower.tail=FALSE)) #this is the p-value for this complex current
      }
      #return(p_val)
    }
    output_mat <- rbind(output_mat, p_val)
  }
  
  colnames(output_mat) <- names(complex_vector)
  rownames(output_mat) <- colnames(data)
  return(output_mat)
}

hgtest_prot_prob_assign <- function(cplx, p)
{
  
  complex_names_vec <- rep(names(cplx), lapply(cplx, length))
  data<- data.frame(table(complex_names_vec))
  data[,1]=as.numeric(as.character(data[,1]))
  data=data[order(data[,1]),]
  #save(data, file="D:/data.txt", ascii=TRUE)
  prot_vec <- unlist(cplx, use.names=F, recursive=F)
  #print(prot_vec)
  new_mat <- as.matrix(cbind(prot_vec, complex_names_vec, rep(p, data[,2])))
  #now we should have a list of proteins that we can now check individually
  prot_list <- unique(prot_vec) #the unique list of proteins
  output <- c()
  for (i in 1:length(prot_list))
  {
    output<- rbind(output,cbind(prot_list[i], max(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]),3]))))
    #output<- rbind(output,cbind(prot_list[i], mean(unlist(new_mat[which(new_mat[,1]%in%prot_list[i]),3]))))
  }
  return(output)
}
## sce b GSEA-----------------------------
mat.ttest <- function(data)
{ output_mat<- c()
datat=as.data.frame(data)
for (i in 1:ncol(data))
{rep_pval<-c()
for (j in 1:nrow(data))
{ ttmp=as.matrix(datat[j,])
#print(ttmp)
pval <- rowTtest(as.matrix(ttmp),y=NULL, mu=datat[j,i])$p.value
#print(pval)
rep_pval<- append(rep_pval,pval)
}
output_mat<- cbind(output_mat,rep_pval) 
rep_pval=c()

}

rownames(output_mat)<- rownames(data)
colnames(output_mat)<- colnames(data)
return(output_mat) }

repgsea <- function(data, complex_vector)
{
  ttest_out <- mat.ttest(data)
  gsea_pvals <- c()
  for (k in 1:ncol(data))
  {
    ranks<- rank(ttest_out[,k])
    rep_gseapval<-c()
    for (x in 1:length(complex_vector))
    {if (length(ranks[which(names(ranks) %in% complex_vector[[x]])]) >= 1)
    {ks_pval <- ks.test(jitter(ranks[which(names(ranks) %in% complex_vector[[x]])]), jitter(ranks[which(!names(ranks) %in% complex_vector[[x]])]))$p.value
    rep_gseapval <- append(rep_gseapval, ks_pval)
    }
      else
      {rep_gseapval <- append(rep_gseapval, 1)
      }}
    gsea_pvals<- rbind(gsea_pvals,rep_gseapval) }
  # print(colnames(gsea_pvals))
  colnames(gsea_pvals)<-names(complex_vector)
  # print(colnames(gsea_pvals))
  rownames(gsea_pvals)<-colnames(data) 
  
  return(gsea_pvals) }

## sce b recovery function--------------------
pairwise_recovery <- function(predict_list, original_prot_list, check_prot_list, complex_vec)
{
  #extract proteins from complexes with significant p-values
  
  predict_list <- predict_list[as.numeric(predict_list) <= 0.05]
  
  add_proteins <- setdiff(unlist(complex_vec[names(predict_list)]), original_prot_list) #this is the set of recovered proteins associated with the significant complexes
  
  obs_intersect <- round(length(intersect(add_proteins, check_prot_list))/length( add_proteins), 3)
  if(length(add_proteins)==0) obs_intersect=0
  theoretical_vec <- c()
  
  for (i in 1:1000)
  {
    S_prime <- sample(unlist(complex_vec), length(add_proteins)) 
    add_proteins_prime <- setdiff(S_prime, original_prot_list)
    theoretical_intersect <- round(length(intersect(add_proteins_prime, check_prot_list))/length(add_proteins_prime), 3)
    if(!is.na(theoretical_intersect))
      theoretical_vec <- append(theoretical_vec, theoretical_intersect)
  }  
  
  p_value <- sum(theoretical_vec >= obs_intersect)/length(theoretical_vec)
  a=length(add_proteins)
  b=length(intersect(add_proteins, check_prot_list))
  #print(obs_intersect, p_value)
  if(obs_intersect==0)
  {
    p_value=0
    a=0
    b=0
    
  }
  
  return(c(obs_intersect, p_value, a, b, add_proteins,"a",intersect(add_proteins, check_prot_list)))
  #write.table(add_proteins, file="test.txt", col.names=F, row.names=F)
}

ntop_recovery_protrec <- function(prot_predict_list, original_prot_list, check_prot_list, complex_vec,aaa)
  # aaa is the number of predicted proteins with pval < 0.05 for the particular method
{
  
  prot_predict_list <- prot_predict_list[1:aaa,] # take up to aaa number of proteins in the predicted list
  
  add_proteins <- setdiff(prot_predict_list[,1], original_prot_list) 
  obs_intersect <- round(length(intersect(add_proteins, check_prot_list))/length(add_proteins), 3)
  
  theoretical_vec <- c()
  
  for (i in 1:1000)
  {
    S_prime <- sample(unlist(complex_vec), length(add_proteins)) 
    add_proteins_prime <- setdiff(S_prime, original_prot_list)
    theoretical_intersect <- round(length(intersect(add_proteins_prime, check_prot_list))/length(add_proteins_prime), 3)
    theoretical_vec <- append(theoretical_vec, theoretical_intersect)
  }  
  p_value <- sum(as.numeric(theoretical_vec >= obs_intersect))/length(theoretical_vec)
  
  return(c(obs_intersect, p_value, length(add_proteins), length(intersect(add_proteins, check_prot_list)), add_proteins))
  
}
## sce b RC------------------------------
#gsea doesnt need complex size info, so only run it once and reuse it
rc_ngsea_p <- repgsea(RC_N,complex_vector)
rc_cgsea_p <- repgsea(RC_C,complex_vector)

path="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_july/corum2022/RC/sce_b_"

#rcn loop
get_result_sce_b_rcn<-function(mode,complex_size,rc_nfcs, rc_nprotrec, rc_nhg){
  print("running...")
  print(paste0(path,mode,"/result_rcn_topn_",complex_size,"_",mode,".csv"))
  
  
  profcsmap=c()
  prohgmap=c()
  progseamap=c()
  for(i in 1:6)
  {
    PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_nprotrec[i,], rownames(RC_N),0.01,mode))
    rownames(PROTREC_prot)<- PROTREC_prot[,1]
    PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
    PROTREC_prot<- data.frame(PROTREC_prot)
    t1=data.frame(PROTREC_prot[order(PROTREC_prot[,2],decreasing=T),])
    rownames(t1)=rownames(PROTREC_prot)
    
    print("fcs started")
    t2=data.frame(fcs_prot_prob_assign(complex_vector,1-rc_nfcs[i,]))
    t2[,2]=1-as.numeric(t2[,2])
    t2=t2[as.numeric(t2[,2])<0.05,]
    l2=nrow(t2)
    
    print("he started")
    t3=data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_nhg[i,]))
    t3[,2]=1-as.numeric(t3[,2])
    t3=t3[as.numeric(t3[,2])<0.05,]
    l3=nrow(t3)
    
    print("gsea started")
    t4=data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_ngsea_p[i,]))
    t4[,2]=1-as.numeric(t4[,2])
    t4=t4[as.numeric(t4[,2])<0.05,]
    l4=nrow(t4)
    
    val<- pairwise_recovery(t(rc_nfcs)[,i],rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector,l2)
    profcsmap=rbind(profcsmap,cbind(rep(l2,4),val[1:4],valpro[1:4]))
    
    val<- pairwise_recovery(t(rc_nhg)[,i],rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector,l3)
    prohgmap=rbind(prohgmap,cbind(rep(l3,4),val[1:4],valpro[1:4]))
    
    val<- pairwise_recovery(t(rc_ngsea_p)[,i],rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector,l4)
    progseamap=rbind(progseamap,cbind(rep(l4,4),val[1:4],valpro[1:4]))
    print("1-6:")
    print(i)
  }
  for(i in 7:12)
  {
    PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_nprotrec[i,], rownames(RC_N),0.01,mode))
    rownames(PROTREC_prot)<- PROTREC_prot[,1]
    PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
    PROTREC_prot<- data.frame(PROTREC_prot)
    t1=data.frame(PROTREC_prot[order(PROTREC_prot[,2],decreasing=T),])
    rownames(t1)=rownames(PROTREC_prot)
    
    print("fcs started")
    t2=data.frame(fcs_prot_prob_assign(complex_vector,1-rc_nfcs[i,]))
    t2[,2]=1-as.numeric(t2[,2])
    t2=t2[as.numeric(t2[,2])<0.05,]
    l2=nrow(t2)
    
    print("he started")
    t3=data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_nhg[i,]))
    t3[,2]=1-as.numeric(t3[,2])
    t3=t3[as.numeric(t3[,2])<0.05,]
    l3=nrow(t3)
    
    print("gsea started")
    t4=data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_ngsea_p[i,]))
    t4[,2]=1-as.numeric(t4[,2])
    t4=t4[as.numeric(t4[,2])<0.05,]
    l4=nrow(t4)
    
    val<- pairwise_recovery(t(rc_nfcs)[,i],rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector,l2)
    profcsmap=rbind(profcsmap,cbind(rep(l2,4),val[1:4],valpro[1:4]))
    
    val<- pairwise_recovery(t(rc_nhg)[,i],rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector,l3)
    prohgmap=rbind(prohgmap,cbind(rep(l3,4),val[1:4],valpro[1:4]))
    
    val<- pairwise_recovery(t(rc_ngsea_p)[,i],rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_N)[(RC_N[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector,l4)
    progseamap=rbind(progseamap,cbind(rep(l4,4),val[1:4],valpro[1:4]))
    print("7-12:")
    print(i)
  }
  colnames(profcsmap)=c("number","fcs","protrec")
  colnames(prohgmap)=c("number","he","protrec")
  colnames(progseamap)=c("number","gsea","protrec")
  all1v1=c()
  all1v1=cbind(profcsmap,prohgmap,progseamap)
  write.csv(all1v1,paste0(path,mode,"/result_rcn_topn_",complex_size,"_",mode,".csv"))
  print(paste0("rcn_topn_",complex_size,"_",mode))
}

#rcc loop
get_result_sce_b_rcc<-function(mode,complex_size, rc_cfcs, rc_cprotrec, rc_chg){
  print("running...")
  print(paste0(path,mode,"/result_rcc_topn_",complex_size,"_",mode,".csv"))
  
  profcsmap=c()
  prohgmap=c()
  progseamap=c()
  for(i in 1:6)
  {
    PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[i,], rownames(RC_C),0.01,mode))
    rownames(PROTREC_prot)<- PROTREC_prot[,1]
    PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
    PROTREC_prot<- data.frame(PROTREC_prot)
    t1=data.frame(PROTREC_prot[order(PROTREC_prot[,2],decreasing=T),])
    rownames(t1)=rownames(PROTREC_prot)
    
    print("fcs started")
    t2=data.frame(fcs_prot_prob_assign(complex_vector,1-rc_cfcs[i,]))
    t2[,2]=1-as.numeric(t2[,2])
    t2=t2[as.numeric(t2[,2])<0.05,]
    l2=nrow(t2)
    
    print("he started")
    t3=data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_chg[i,]))
    t3[,2]=1-as.numeric(t3[,2])
    t3=t3[as.numeric(t3[,2])<0.05,]
    l3=nrow(t3)
    
    print("gsea started")
    t4=data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_cgsea_p[i,]))
    t4[,2]=1-as.numeric(t4[,2])
    t4=t4[as.numeric(t4[,2])<0.05,]
    l4=nrow(t4)
    
    val<- pairwise_recovery(t(rc_cfcs)[,i],rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector,l2)
    profcsmap=rbind(profcsmap,cbind(rep(l2,4),val[1:4],valpro[1:4]))
    
    val<- pairwise_recovery(t(rc_chg)[,i],rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector,l3)
    prohgmap=rbind(prohgmap,cbind(rep(l3,4),val[1:4],valpro[1:4]))
    
    val<- pairwise_recovery(t(rc_cgsea_p)[,i],rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i])!=0],complex_vector,l4)
    progseamap=rbind(progseamap,cbind(rep(l4,4),val[1:4],valpro[1:4]))
    print("1-6:")
    print(i)
  }
  for(i in 7:12)
  {
    PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[i,], rownames(RC_C),0.01,mode))
    rownames(PROTREC_prot)<- PROTREC_prot[,1]
    PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
    PROTREC_prot<- data.frame(PROTREC_prot)
    t1=data.frame(PROTREC_prot[order(PROTREC_prot[,2],decreasing=T),])
    rownames(t1)=rownames(PROTREC_prot)
    
    print("fcs started")
    t2=data.frame(fcs_prot_prob_assign(complex_vector,1-rc_cfcs[i,]))
    t2[,2]=1-as.numeric(t2[,2])
    t2=t2[as.numeric(t2[,2])<0.05,]
    l2=nrow(t2)
    
    print("he started")
    t3=data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_chg[i,]))
    t3[,2]=1-as.numeric(t3[,2])
    t3=t3[as.numeric(t3[,2])<0.05,]
    l3=nrow(t3)
    
    print("gsea started")
    t4=data.frame(hgtest_prot_prob_assign(complex_vector,1-rc_cgsea_p[i,]))
    t4[,2]=1-as.numeric(t4[,2])
    t4=t4[as.numeric(t4[,2])<0.05,]
    l4=nrow(t4)
    
    val<- pairwise_recovery(t(rc_cfcs)[,i],rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector,l2)
    profcsmap=rbind(profcsmap,cbind(rep(l2,4),val[1:4],valpro[1:4]))
    
    val<- pairwise_recovery(t(rc_chg)[,i],rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector,l3)
    prohgmap=rbind(prohgmap,cbind(rep(l3,4),val[1:4],valpro[1:4]))
    
    val<- pairwise_recovery(t(rc_cgsea_p)[,i],rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector)
    valpro=ntop_recovery_protrec(t1,rownames(RC_C)[(RC_C[,i])!=0],rownames(RC_peptides_uniq)[(RC_peptides_uniq[,i+12])!=0],complex_vector,l4)
    progseamap=rbind(progseamap,cbind(rep(l4,4),val[1:4],valpro[1:4]))
    print("7-12:")
    print(i)
  }
  colnames(profcsmap)=c("number","fcs","protrec")
  colnames(prohgmap)=c("number","he","protrec")
  colnames(progseamap)=c("number","gsea","protrec")
  all1v1=c()
  all1v1=cbind(profcsmap,prohgmap,progseamap)
  write.csv(all1v1,paste0(path,mode,"/result_rcc_topn_",complex_size,"_",mode,".csv"))
  print(paste0("rcc_topn_",complex_size,"_",mode))
}

for(sssize in c(5,10,20,30,40,50,200)){
  #RC_N
  print("fcs vector")
  rc_nfcs <- fcs(RC_N,complex_vector,1000,sssize)
  print("protrec vector")
  rc_nprotrec <- PROTREC_cplx_prob(RC_N, complex_vector, 0.01, sssize)
  print("he vector")
  rc_nhg<- hgtest(RC_N, complex_vector,sssize)
  
  print("finished generating RCN complex vectors")
  
  get_result_sce_b_rcn("max",sssize, rc_nfcs, rc_nprotrec, rc_nhg)
  get_result_sce_b_rcn("min",sssize, rc_nfcs, rc_nprotrec, rc_nhg)
  get_result_sce_b_rcn("median",sssize, rc_nfcs, rc_nprotrec, rc_nhg)
  get_result_sce_b_rcn("mean",sssize, rc_nfcs, rc_nprotrec, rc_nhg)
  
  #RC_C
  print("fcs vector")
  rc_cfcs <- fcs(RC_C,complex_vector,1000,sssize)
  print("protrec vector")
  rc_cprotrec <- PROTREC_cplx_prob(RC_C, complex_vector, 0.01, sssize)
  print("he vector")
  rc_chg<- hgtest(RC_C, complex_vector,sssize)
  
  print("finished generating RCC complex vectors")
  
  get_result_sce_b_rcc("max",sssize, rc_cfcs, rc_cprotrec, rc_chg)
  get_result_sce_b_rcc("min",sssize, rc_cfcs, rc_cprotrec, rc_chg)
  get_result_sce_b_rcc("median",sssize, rc_cfcs, rc_cprotrec, rc_chg)
  get_result_sce_b_rcc("mean",sssize, rc_cfcs, rc_cprotrec, rc_chg)
}




## sce b hela siha-----------------
#GSEA doesnt need complex size information
heladdagsea_p <- repgsea(Hela_ddaProteins,complex_vector)
sihaddagsea_p <- repgsea(Siha_ddaProteins,complex_vector)

get_result_sce_b_hela<-function(mode,complex_size,heladdaprotrec,heladdafcs,heladdahg){
  print("running...")
  print(paste0("hela",complex_size,"_",mode))
  
  profcsmap=c()
  prohgmap=c()
  progseamap=c()
  for(i in 1:3)
  {
    for(j in 1:2)
    {
      k=(i+j)%%3
      if(k==0)k=3
      
      PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-heladdaprotrec[i,], rownames(Hela_ddaProteins),0.01,mode))
      rownames(PROTREC_prot)<- PROTREC_prot[,1]
      PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
      PROTREC_prot<- data.frame(PROTREC_prot)
      t1=data.frame(PROTREC_prot[order(PROTREC_prot[,2],decreasing=T),])
      rownames(t1)=rownames(PROTREC_prot)
      
      print("fcs started")
      t2=data.frame(fcs_prot_prob_assign(complex_vector,1-heladdafcs[i,]))
      t2[,2]=1-as.numeric(t2[,2])
      t2=t2[as.numeric(t2[,2])<0.05,]
      l2=nrow(t2)
      
      print("he started")
      t3=data.frame(hgtest_prot_prob_assign(complex_vector,1-heladdahg[i,]))
      t3[,2]=1-as.numeric(t3[,2])
      t3=t3[as.numeric(t3[,2])<0.05,]
      l3=nrow(t3)
      
      print("gsea started")
      t4=data.frame(hgtest_prot_prob_assign(complex_vector,1-heladdagsea_p[i,]))
      t4[,2]=1-as.numeric(t4[,2])
      t4=t4[as.numeric(t4[,2])<0.05,]
      l4=nrow(t4)
      
      val<- pairwise_recovery(t(heladdafcs)[,i],rownames(Hela_ddaProteins)[(Hela_ddaProteins[,i])!=0],rownames(heladdapepuniq)[(heladdapepuniq[,k])!=0],complex_vector)
      valpro=ntop_recovery_protrec(t1,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,i])!=0],rownames(heladdapepuniq)[(heladdapepuniq[,k])!=0],complex_vector,l2)
      profcsmap=rbind(profcsmap,cbind(rep(l2,4),val[1:4],valpro[1:4]))
      
      val<- pairwise_recovery(t(heladdahg)[,i],rownames(Hela_ddaProteins)[(Hela_ddaProteins[,i])!=0],rownames(heladdapepuniq)[(heladdapepuniq[,k])!=0],complex_vector)
      valpro=ntop_recovery_protrec(t1,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,i])!=0],rownames(heladdapepuniq)[(heladdapepuniq[,k])!=0],complex_vector,l3)
      prohgmap=rbind(prohgmap,cbind(rep(l3,4),val[1:4],valpro[1:4]))
      
      val<- pairwise_recovery(t(heladdagsea_p)[,i],rownames(Hela_ddaProteins)[(Hela_ddaProteins[,i])!=0],rownames(heladdapepuniq)[(heladdapepuniq[,k])!=0],complex_vector)
      valpro=ntop_recovery_protrec(t1,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,i])!=0],rownames(heladdapepuniq)[(heladdapepuniq[,k])!=0],complex_vector,l4)
      progseamap=rbind(progseamap,cbind(rep(l4,4),val[1:4],valpro[1:4]))
    }
  }
  
  colnames(profcsmap)=c("number","fcs","protrec")
  colnames(prohgmap)=c("number","he","protrec")
  colnames(progseamap)=c("number","gsea","protrec")
  all1v1=c()
  all1v1=cbind(profcsmap,prohgmap,progseamap)
  write.csv(all1v1,paste0(path,mode,"/result_hela_topn_",complex_size,"_",mode,".csv"))
  print(paste0("hela",complex_size,"_",mode))
}

get_result_sce_b_siha<-function(mode,complex_size,sihaddaprotrec,sihaddafcs,sihaddahg){
  print("running...")
  print(paste0("siha",complex_size,"_",mode))
  
  profcsmap=c()
  prohgmap=c()
  progseamap=c()
  for(i in 1:3)
  {
    for(j in 1:2)
    {
      k=(i+j)%%3
      if(k==0)k=3
      
      PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-sihaddaprotrec[i,], rownames(Siha_ddaProteins),0.01,mode))
      rownames(PROTREC_prot)<- PROTREC_prot[,1]
      PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
      PROTREC_prot<- data.frame(PROTREC_prot)
      t1=data.frame(PROTREC_prot[order(PROTREC_prot[,2],decreasing=T),])
      rownames(t1)=rownames(PROTREC_prot)
      
      print("fcs started")
      t2=data.frame(fcs_prot_prob_assign(complex_vector,1-sihaddafcs[i,]))
      t2[,2]=1-as.numeric(t2[,2])
      t2=t2[as.numeric(t2[,2])<0.05,]
      l2=nrow(t2)
      
      print("he started")
      t3=data.frame(hgtest_prot_prob_assign(complex_vector,1-sihaddahg[i,]))
      t3[,2]=1-as.numeric(t3[,2])
      t3=t3[as.numeric(t3[,2])<0.05,]
      l3=nrow(t3)
      
      print("gsea started")
      t4=data.frame(hgtest_prot_prob_assign(complex_vector,1-sihaddagsea_p[i,]))
      t4[,2]=1-as.numeric(t4[,2])
      t4=t4[as.numeric(t4[,2])<0.05,]
      l4=nrow(t4)
      
      val<- pairwise_recovery(t(sihaddafcs)[,i],rownames(Siha_ddaProteins)[(Siha_ddaProteins[,i])!=0],rownames(sihaddapepuniq)[(sihaddapepuniq[,k])!=0],complex_vector)
      valpro=ntop_recovery_protrec(t1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,i])!=0],rownames(sihaddapepuniq)[(sihaddapepuniq[,k])!=0],complex_vector,l2)
      profcsmap=rbind(profcsmap,cbind(rep(l2,4),val[1:4],valpro[1:4]))
      
      val<- pairwise_recovery(t(sihaddahg)[,i],rownames(Siha_ddaProteins)[(Siha_ddaProteins[,i])!=0],rownames(sihaddapepuniq)[(sihaddapepuniq[,k])!=0],complex_vector)
      valpro=ntop_recovery_protrec(t1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,i])!=0],rownames(sihaddapepuniq)[(sihaddapepuniq[,k])!=0],complex_vector,l3)
      prohgmap=rbind(prohgmap,cbind(rep(l3,4),val[1:4],valpro[1:4]))
      
      val<- pairwise_recovery(t(sihaddagsea_p)[,i],rownames(Siha_ddaProteins)[(Siha_ddaProteins[,i])!=0],rownames(sihaddapepuniq)[(sihaddapepuniq[,k])!=0],complex_vector)
      valpro=ntop_recovery_protrec(t1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,i])!=0],rownames(sihaddapepuniq)[(sihaddapepuniq[,k])!=0],complex_vector,l4)
      progseamap=rbind(progseamap,cbind(rep(l4,4),val[1:4],valpro[1:4]))
    }
  }
  
  colnames(profcsmap)=c("number","fcs","protrec")
  colnames(prohgmap)=c("number","he","protrec")
  colnames(progseamap)=c("number","gsea","protrec")
  all1v1=c()
  all1v1=cbind(profcsmap,prohgmap,progseamap)
  write.csv(all1v1,paste0(path,mode,"/result_siha_topn_",complex_size,"_",mode,".csv"))
  print(paste0("siha",complex_size,"_",mode))
}

path="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_july/corum2022/helasiha/sce_b_"

for(sssize in c(5,10,20,30,40,50,200)){
  #hela
  print("protrec vector")
  heladdaprotrec <- PROTREC_cplx_prob(Hela_ddaProteins, complex_vector, 0.01, sssize)
  #FCS
  print("fcs vector")
  heladdafcs <- fcs(Hela_ddaProteins,complex_vector,1000,sssize)
  #HE
  print("he vector")
  heladdahg<- hgtest(Hela_ddaProteins, complex_vector,sssize)
  
  print("finished generating hela complex vectors")
  
  get_result_sce_b_hela("max",sssize,heladdaprotrec,heladdafcs,heladdahg)
  get_result_sce_b_hela("min",sssize,heladdaprotrec,heladdafcs,heladdahg)
  get_result_sce_b_hela("mean",sssize,heladdaprotrec,heladdafcs,heladdahg)
  get_result_sce_b_hela("median",sssize,heladdaprotrec,heladdafcs,heladdahg)
  
  #siha
  print("protrec vector")
  sihaddaprotrec <- PROTREC_cplx_prob(Siha_ddaProteins, complex_vector, 0.01, sssize)
  #FCS
  print("fcs vector")
  sihaddafcs <- fcs(Siha_ddaProteins,complex_vector,1000,sssize)
  #HE
  print("he vector")
  sihaddahg<- hgtest(Siha_ddaProteins, complex_vector,sssize)
  
  print("finished generating siha complex vectors")
  
  get_result_sce_b_siha("max",sssize,sihaddaprotrec,sihaddafcs,sihaddahg)
  get_result_sce_b_siha("min",sssize,sihaddaprotrec,sihaddafcs,sihaddahg)
  get_result_sce_b_siha("mean",sssize,sihaddaprotrec,sihaddafcs,sihaddahg)
  get_result_sce_b_siha("median",sssize,sihaddaprotrec,sihaddafcs,sihaddahg)
  
}



#------sce c-------------------------
## sce c two peptide rc----------
sce_c_get_plis_rc <- function(data_path, tissue_type){
  RC_peptides2= read.delim(data_path)
  if(tissue_type == 'rcn'){
    pp=cbind(RC_peptides2[,1:2],RC_peptides2[,5:10],RC_peptides2[,17:22])}
  else{
    pp=cbind(RC_peptides2[,1:2],RC_peptides2[,11:16],RC_peptides2[,23:28])
  }
  ppp=pp[pp[,2]==1,]
  pppp=c()
  for(i in 1:nrow(ppp))
  {
    l=length(which(is.na(ppp[i,3:14]) %in% "TRUE"))
    if(l<12)
    {
      pppp=rbind(pppp,ppp[i,1:2])
    }
  }
  lis=unique(pppp[,1])
  plis=c()
  for(i in 1:length(lis))
  {
    tmp=lis[i]
    g=which(pppp[,1] %in% tmp)
    if(length(g)>1) plis=append(plis,tmp)
  }
  return(plis)
}
#RC_N
plis_rcn <- sce_c_get_plis_rc("./peptides_rc.txt",'rcn')
#RC_C
plis_rcc <- sce_c_get_plis_rc("./peptides_rc.txt",'rcc')


## sce c two peptide hela siha------------------------------------------------------------------------------------------------------------
sce_c_get_ampeptides_helasiha <- function(ddapep){
  ampeptide=c()
  for (rowIndex in 1:nrow(ddapep))
  {
    if(as.numeric(ddapep[rowIndex,1]==1))
    {
      s = ddapep[rowIndex, 18]
      splitt = strsplit(s, ":", fixed = TRUE)
      splitm=unlist(splitt)
      if(length(splitm)>=1)
      {
        js=length(splitm)
        for(i in 1:length(splitm))
        {
          splits=strsplit(splitm[i],"|",fixed=TRUE)[1]
          splitS=splits[[1]]
          if (length(splitS) >= 1) 
          {
            subS = splitS[1]
            ampeptide=rbind(ampeptide,cbind(subS,ddapep[rowIndex,4],js))
          }
        }
      }
    }
  }
  ampeptide2=as.data.frame(ampeptide)
  ampeptide2=unique(ampeptide2)
  return(ampeptide2)
}

hela_ampeptide2 = sce_c_get_ampeptides_helasiha(Hela_ddapep)
siha_ampeptide2 = sce_c_get_ampeptides_helasiha(Siha_ddapep)

helapep <- read.csv("./DB_search_psm_hela.csv",header=TRUE)
sihapep <- read.csv("./DB_search_psm_siha.csv",header=TRUE)

sce_c_get_prolis_helasiha<-function(pep, ampeptide2){
  p=cbind(pep[,1],pep[,14])
  p=unique(p)
  plis=c()
  for(i in 1:nrow(p))
  {
    tmp=p[i,2]
    if(tmp=="")
    {
    }
    else
    {
      splitt = strsplit(tmp, ":", fixed = TRUE)
      splitm=unlist(splitt)
      if(length(splitm)==1)
      {
        splits=strsplit(splitm,"|",fixed=TRUE)[1]
        a=splits[[1]]
        plis=append(plis,a[1])
      }
    }
  }
  lis=table(plis)
  lis=data.frame(lis)
  lis=lis[lis[,2]>=2,]
  prolis=as.character(lis[,1])
  ###################
  am=cbind(ampeptide2[,1],as.numeric(ampeptide2[,3]))
  am=data.frame(am)
  am[,2]=as.numeric(am[,2])
  am=am[am[,2]==1,]
  g=table(am[,1])
  g=data.frame(g)
  g=g[g[,2]>1,]
  prolis2=as.character(g[,1])
  ###################
  prolis=union(prolis,prolis2)
  return(prolis)
}

hela_prolis = sce_c_get_prolis_helasiha(helapep, hela_ampeptide2)
siha_prolis = sce_c_get_prolis_helasiha(sihapep, siha_ampeptide2)


## sce c rc----------
get_result_sce_c_rcn<-function(rc_nprotrec,protrecscoreset,mode,size){
  RC_recov_protrec_NT1NT11 <-c()
  for (i in c(1:12))  #loop for NT1 -> NT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_nprotrec[i,], rownames(RC_N),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val1<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_N)[(RC_N[,i])!=0],plis_rcn,complex_vector,protrecscoreset)
  RC_recov_protrec_NT1NT11<- rbind(RC_recov_protrec_NT1NT11,val1[1:4])
  print(i)
  }
  print(paste0("./r_code_output/RC/twopeptide/result_rcn_twopept_",protrecscoreset,"_",size,"_",mode))
  return(RC_recov_protrec_NT1NT11)
}

get_result_sce_c_rcc<-function(rc_cprotrec,protrecscoreset,mode,size){
  RC_recov_protrec_CT1CT11 <-c()
  for (i in c(1:12))  #loop for CT1 -> CT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[i,], rownames(RC_C),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val1<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_C)[(RC_C[,i])!=0],plis_rcc,complex_vector,protrecscoreset)
  RC_recov_protrec_CT1CT11<- rbind(RC_recov_protrec_CT1CT11,val1[1:4])
  print(i)
  }
  print(paste0("./r_code_output/RC/twopeptide/result_rcc_twopept_",protrecscoreset,"_",size,"_",mode))
  return(RC_recov_protrec_CT1CT11)
}

path="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_july/corum2022/RC/sce_c_"

for(protrecscoreset in c(0.4,0.5,0.6,0.7,0.8,0.95))
{
  for(mode in c("max","min","median","mean"))
  {
    rcn_result5 = get_result_sce_c_rcn(rc_nprotrec,protrecscoreset,mode,5)
    write.csv(rcn_result5,paste0(path,mode,"/result_rcn_twopept_",protrecscoreset,"_5_",mode,".csv"))
    
    rcc_result5 = get_result_sce_c_rcc(rc_cprotrec,protrecscoreset,mode,5)
    write.csv(rcc_result5,paste0(path,mode,"/result_rcc_twopept_",protrecscoreset,"_5_",mode,".csv"))
    
    
    rcn_result10 = get_result_sce_c_rcn(rc_nprotrec10,protrecscoreset,mode,10)
    write.csv(rcn_result10,paste0(path,mode,"/result_rcn_twopept_",protrecscoreset,"_10_",mode,".csv"))
    
    rcc_result10 = get_result_sce_c_rcc(rc_cprotrec10,protrecscoreset,mode,10)
    write.csv(rcc_result10,paste0(path,mode,"/result_rcc_twopept_",protrecscoreset,"_10_",mode,".csv"))
    
    
    rcn_result20 = get_result_sce_c_rcn(rc_nprotrec20,protrecscoreset,mode,20)
    write.csv(rcn_result20,paste0(path,mode,"/result_rcn_twopept_",protrecscoreset,"_20_",mode,".csv"))
    
    rcc_result20 = get_result_sce_c_rcc(rc_cprotrec20,protrecscoreset,mode,20)
    write.csv(rcc_result20,paste0(path,mode,"/result_rcc_twopept_",protrecscoreset,"_20_",mode,".csv"))
    
    
    rcn_result30 = get_result_sce_c_rcn(rc_nprotrec30,protrecscoreset,mode,30)
    write.csv(rcn_result30,paste0(path,mode,"/result_rcn_twopept_",protrecscoreset,"_30_",mode,".csv"))
    
    rcc_result30 = get_result_sce_c_rcc(rc_cprotrec30,protrecscoreset,mode,30)
    write.csv(rcc_result30,paste0(path,mode,"/result_rcc_twopept_",protrecscoreset,"_30_",mode,".csv"))
    
    
    rcn_result40 = get_result_sce_c_rcn(rc_nprotrec40,protrecscoreset,mode,40)
    write.csv(rcn_result40,paste0(path,mode,"/result_rcn_twopept_",protrecscoreset,"_40_",mode,".csv"))
    
    rcc_result40 = get_result_sce_c_rcc(rc_cprotrec40,protrecscoreset,mode,40)
    write.csv(rcc_result40,paste0(path,mode,"/result_rcc_twopept_",protrecscoreset,"_40_",mode,".csv"))
    
    
    rcn_result50 = get_result_sce_c_rcn(rc_nprotrec50,protrecscoreset,mode,50)
    write.csv(rcn_result50,paste0(path,mode,"/result_rcn_twopept_",protrecscoreset,"_50_",mode,".csv"))
    
    rcc_result50 = get_result_sce_c_rcc(rc_cprotrec50,protrecscoreset,mode,50)
    write.csv(rcc_result50,paste0(path,mode,"/result_rcc_twopept_",protrecscoreset,"_50_",mode,".csv"))
    
    
    rcn_result200 = get_result_sce_c_rcn(rc_nprotrec200,protrecscoreset,mode,200)
    write.csv(rcn_result200,paste0(path,mode,"/result_rcn_twopept_",protrecscoreset,"_200_",mode,".csv"))
    
    rcc_result200 = get_result_sce_c_rcc(rc_cprotrec200,protrecscoreset,mode,200)
    write.csv(rcc_result200,paste0(path,mode,"/result_rcc_twopept_",protrecscoreset,"_200_",mode,".csv"))
  }
}


## sce c hela siha-------------------
get_result_sce_c_hela <- function(hela_ddaprotrec,protrecscoreset,mode,complex_size){
  prores=c()
  for(i in 1:3)
  {
    PROTREC_Heladdaprot1 <- data.frame(PROTREC_protprob(complex_vector, 1-hela_ddaprotrec[i,], rownames(Hela_ddaProteins),0.01,mode))
    rownames(PROTREC_Heladdaprot1)<- PROTREC_Heladdaprot1[,1] # cant set name here, why?
    PROTREC_Heladdaprot1[,2]<-as.numeric(as.character(PROTREC_Heladdaprot1[,2]))
    PROTREC_Heladdaprot1<- data.frame(PROTREC_Heladdaprot1)
    val<- pairwise_recovery_protrec(PROTREC_Heladdaprot1,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,i]!=0)],hela_prolis,complex_vector,protrecscoreset)
    prores=rbind(prores,val[1:4])
  }
  write.csv(prores,paste0(path,mode,"/result_hela_twopeptide",protrecscoreset,"_",complex_size,"_",mode,".csv"))
}

get_result_sce_c_siha <- function(siha_ddaprotrec,protrecscoreset,mode,complex_size){
  prores=c()
  for(i in 1:3)
  {
    PROTREC_Sihaddaprot1 <- data.frame(PROTREC_protprob(complex_vector, 1-siha_ddaprotrec[i,], rownames(Siha_ddaProteins),0.01,mode))
    rownames(PROTREC_Sihaddaprot1)<- PROTREC_Sihaddaprot1[,1]
    PROTREC_Sihaddaprot1[,2]<-as.numeric(as.character(PROTREC_Sihaddaprot1[,2]))
    PROTREC_Sihaddaprot1<- data.frame(PROTREC_Sihaddaprot1)
    val<- pairwise_recovery_protrec(PROTREC_Sihaddaprot1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,i]!=0)],siha_prolis,complex_vector,protrecscoreset)
    prores=rbind(prores,val[1:4])
  }
  write.csv(prores,paste0(path,mode,"/result_siha_twopeptide",protrecscoreset,"_",complex_size,"_",mode,".csv"))
}

path="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_july/corum2022/helasiha/sce_c_"

for(protrecscoreset in c(0.4,0.5,0.6,0.7,0.8,0.95))
{
  for(mode in c("max","min","median","mean"))
  {
    result = get_result_sce_c_hela(hela_ddaprotrec,protrecscoreset,mode,5)
    print(paste0("hela_cross_score",protrecscoreset,"_5_",mode))
    
    result1 = get_result_sce_c_hela(hela_ddaprotrec10,protrecscoreset,mode,10)
    print(paste0("hela_cross_score",protrecscoreset,"_10_",mode))
    
    result2 = get_result_sce_c_hela(hela_ddaprotrec20,protrecscoreset,mode,20)
    print(paste0("hela_cross_score",protrecscoreset,"_20_",mode))
    
    result3 = get_result_sce_c_hela(hela_ddaprotrec30,protrecscoreset,mode,30)
    print(paste0("hela_cross_score",protrecscoreset,"_30_",mode))
    
    result4 = get_result_sce_c_hela(hela_ddaprotrec40,protrecscoreset,mode,40)
    print(paste0("hela_cross_score",protrecscoreset,"_40_",mode))
    
    result5 = get_result_sce_c_hela(hela_ddaprotrec50,protrecscoreset,mode,50)
    print(paste0("hela_cross_score",protrecscoreset,"_50_",mode))
    
    result6 = get_result_sce_c_hela(hela_ddaprotrec200,protrecscoreset,mode,200)
    print(paste0("hela_cross_score",protrecscoreset,"_200_",mode))
    
    
    resultt = get_result_sce_c_siha(siha_ddaprotrec,protrecscoreset,mode,5)
    print(paste0("siha_cross_score",protrecscoreset,"_5_",mode))
    
    resultt1 = get_result_sce_c_siha(siha_ddaprotrec10,protrecscoreset,mode,10)
    print(paste0("siha_cross_score",protrecscoreset,"_10_",mode))
    
    resultt2 = get_result_sce_c_siha(siha_ddaprotrec20,protrecscoreset,mode,20)
    print(paste0("siha_cross_score",protrecscoreset,"_20_",mode))
    
    resultt3 = get_result_sce_c_siha(siha_ddaprotrec30,protrecscoreset,mode,30)
    print(paste0("siha_cross_score",protrecscoreset,"_30_",mode))
    
    resultt4 = get_result_sce_c_siha(siha_ddaprotrec40,protrecscoreset,mode,40)
    print(paste0("siha_cross_score",protrecscoreset,"_40_",mode))
    
    resultt5 = get_result_sce_c_siha(siha_ddaprotrec50,protrecscoreset,mode,50)
    print(paste0("siha_cross_score",protrecscoreset,"_50_",mode))
    
    resultt6 = get_result_sce_c_siha(siha_ddaprotrec200,protrecscoreset,mode,200)
    print(paste0("siha_cross_score",protrecscoreset,"_200_",mode))
    
    print(protrecscoreset)
  }
}

#------sce d-------------------------
## sce d hela siha-------
get_result_sce_d_hela<-function(hela_ddaprotrec, protrecscoreset, mode, complex_size){
  Hela_ddadiacrossrecovprotrec1 <-c()
  PROTREC_heladdaprot1 <- data.frame(PROTREC_protprob(complex_vector, 1-hela_ddaprotrec[1,], rownames(Hela_ddaProteins),0.01,mode))
  rownames(PROTREC_heladdaprot1)<- PROTREC_heladdaprot1[,1]
  PROTREC_heladdaprot1[,2]<-as.numeric(as.character(PROTREC_heladdaprot1[,2]))
  PROTREC_heladdaprot1<- data.frame(PROTREC_heladdaprot1)
  
  PROTREC_heladdaprot2 <- data.frame(PROTREC_protprob(complex_vector, 1-hela_ddaprotrec[2,], rownames(Hela_ddaProteins),0.01,mode))
  rownames(PROTREC_heladdaprot2)<- PROTREC_heladdaprot2[,1]
  PROTREC_heladdaprot2[,2]<-as.numeric(as.character(PROTREC_heladdaprot2[,2]))
  PROTREC_heladdaprot2<- data.frame(PROTREC_heladdaprot2)
  
  PROTREC_heladdaprot3 <- data.frame(PROTREC_protprob(complex_vector, 1-hela_ddaprotrec[3,], rownames(Hela_ddaProteins),0.01,mode))
  rownames(PROTREC_heladdaprot3)<- PROTREC_heladdaprot3[,1]
  PROTREC_heladdaprot3[,2]<-as.numeric(as.character(PROTREC_heladdaprot3[,2]))
  PROTREC_heladdaprot3<- data.frame(PROTREC_heladdaprot3)
  
  heladdaprotrec811diaval<- pairwise_recovery_protrec(PROTREC_heladdaprot1,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,1]!=0)],hela_dia_prot_list1[,1],complex_vector, protrecscoreset)
  Hela_ddadiacrossrecovprotrec1<- rbind(Hela_ddadiacrossrecovprotrec1,heladdaprotrec811diaval[1:4])
  heladdaprotrec811diaval<- pairwise_recovery_protrec(PROTREC_heladdaprot1,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,1]!=0)],hela_dia_prot_list2[,1],complex_vector, protrecscoreset)
  Hela_ddadiacrossrecovprotrec1<- rbind(Hela_ddadiacrossrecovprotrec1,heladdaprotrec811diaval[1:4])
  heladdaprotrec812diaval<- pairwise_recovery_protrec(PROTREC_heladdaprot2,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,2]!=0)],hela_dia_prot_list1[,1],complex_vector, protrecscoreset)
  Hela_ddadiacrossrecovprotrec1<- rbind(Hela_ddadiacrossrecovprotrec1,heladdaprotrec812diaval[1:4])
  heladdaprotrec812diaval<- pairwise_recovery_protrec(PROTREC_heladdaprot2,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,2]!=0)],hela_dia_prot_list2[,1],complex_vector, protrecscoreset)
  Hela_ddadiacrossrecovprotrec1<- rbind(Hela_ddadiacrossrecovprotrec1,heladdaprotrec812diaval[1:4])
  heladdaprotrec813diaval<- pairwise_recovery_protrec(PROTREC_heladdaprot3,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,3]!=0)],hela_dia_prot_list1[,1],complex_vector, protrecscoreset)
  Hela_ddadiacrossrecovprotrec1<- rbind(Hela_ddadiacrossrecovprotrec1,heladdaprotrec813diaval[1:4])
  heladdaprotrec813diaval<- pairwise_recovery_protrec(PROTREC_heladdaprot3,rownames(Hela_ddaProteins)[(Hela_ddaProteins[,3]!=0)],hela_dia_prot_list2[,1],complex_vector, protrecscoreset)
  Hela_ddadiacrossrecovprotrec1<- rbind(Hela_ddadiacrossrecovprotrec1,heladdaprotrec813diaval[1:4])
  
  colnames(Hela_ddadiacrossrecovprotrec1)=c(rep("Protrec",4))
  write.csv(Hela_ddadiacrossrecovprotrec1,paste0(path,mode,"/result_hela_ddadia",protrecscoreset,"_",complex_size,"_",mode,".csv"))
}

get_result_sce_d_siha<-function(siha_ddaprotrec, protrecscoreset, mode, complex_size){
  Siha_ddadiacrossrecovprotrec1 <-c()
  PROTREC_sihaddaprot1 <- data.frame(PROTREC_protprob(complex_vector, 1-siha_ddaprotrec[1,], rownames(Siha_ddaProteins),0.01,mode))
  rownames(PROTREC_sihaddaprot1)<- PROTREC_sihaddaprot1[,1]
  PROTREC_sihaddaprot1[,2]<-as.numeric(as.character(PROTREC_sihaddaprot1[,2]))
  PROTREC_sihaddaprot1<- data.frame(PROTREC_sihaddaprot1)
  
  PROTREC_sihaddaprot2 <- data.frame(PROTREC_protprob(complex_vector, 1-siha_ddaprotrec[2,], rownames(Siha_ddaProteins),0.01,mode))
  rownames(PROTREC_sihaddaprot2)<- PROTREC_sihaddaprot2[,1]
  PROTREC_sihaddaprot2[,2]<-as.numeric(as.character(PROTREC_sihaddaprot2[,2]))
  PROTREC_sihaddaprot2<- data.frame(PROTREC_sihaddaprot2)
  
  PROTREC_sihaddaprot3 <- data.frame(PROTREC_protprob(complex_vector, 1-siha_ddaprotrec[3,], rownames(Siha_ddaProteins),0.01,mode))
  rownames(PROTREC_sihaddaprot3)<- PROTREC_sihaddaprot3[,1]
  PROTREC_sihaddaprot3[,2]<-as.numeric(as.character(PROTREC_sihaddaprot3[,2]))
  PROTREC_sihaddaprot3<- data.frame(PROTREC_sihaddaprot3)
  
  sihaddaprotrec811diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,1]!=0)],siha_dia_prot_list1[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec811diaval[1:4])
  sihaddaprotrec811diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,1]!=0)],siha_dia_prot_list2[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec811diaval[1:4])
  sihaddaprotrec811diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot1,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,1]!=0)],siha_dia_prot_list3[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec811diaval[1:4])
  
  
  sihaddaprotrec812diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot2,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,2]!=0)],siha_dia_prot_list1[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec812diaval[1:4])
  sihaddaprotrec812diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot2,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,2]!=0)],siha_dia_prot_list2[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec812diaval[1:4])
  sihaddaprotrec812diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot2,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,2]!=0)],siha_dia_prot_list3[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec812diaval[1:4])
  
  sihaddaprotrec813diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot3,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,3]!=0)],siha_dia_prot_list1[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec813diaval[1:4])
  sihaddaprotrec813diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot3,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,3]!=0)],siha_dia_prot_list2[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec813diaval[1:4])
  sihaddaprotrec813diaval<- pairwise_recovery_protrec(PROTREC_sihaddaprot3,rownames(Siha_ddaProteins)[(Siha_ddaProteins[,3]!=0)],siha_dia_prot_list3[,1],complex_vector, protrecscoreset)
  Siha_ddadiacrossrecovprotrec1<- rbind(Siha_ddadiacrossrecovprotrec1,sihaddaprotrec813diaval[1:4])
  
  colnames(Siha_ddadiacrossrecovprotrec1)=c(rep("Protrec",4))
  write.csv(Siha_ddadiacrossrecovprotrec1,paste0(path,mode,"/result_siha_ddadia",protrecscoreset,"_",complex_size,"_",mode,".csv"))
}

path="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_july/corum2022/helasiha/sce_d_"

for(protrecscoreset in c(0.4,0.5,0.6,0.7,0.8,0.95))
{
  for(mode in c("max","min","median","mean"))
  {
    result = get_result_sce_d_hela(hela_ddaprotrec,protrecscoreset,mode,5)
    print(paste0("hela_cross_score",protrecscoreset,"_5_",mode))
    
    result1 = get_result_sce_d_hela(hela_ddaprotrec10,protrecscoreset,mode,10)
    print(paste0("hela_cross_score",protrecscoreset,"_10_",mode))
    
    result2 = get_result_sce_d_hela(hela_ddaprotrec20,protrecscoreset,mode,20)
    print(paste0("hela_cross_score",protrecscoreset,"_20_",mode))
    
    result3 = get_result_sce_d_hela(hela_ddaprotrec30,protrecscoreset,mode,30)
    print(paste0("hela_cross_score",protrecscoreset,"_30_",mode))
    
    result4 = get_result_sce_d_hela(hela_ddaprotrec40,protrecscoreset,mode,40)
    print(paste0("hela_cross_score",protrecscoreset,"_40_",mode))
    
    result5 = get_result_sce_d_hela(hela_ddaprotrec50,protrecscoreset,mode,50)
    print(paste0("hela_cross_score",protrecscoreset,"_50_",mode))
    
    result6 = get_result_sce_d_hela(hela_ddaprotrec200,protrecscoreset,mode,200)
    print(paste0("hela_cross_score",protrecscoreset,"_200_",mode))
    
    
    resultt = get_result_sce_d_siha(siha_ddaprotrec,protrecscoreset,mode,5)
    print(paste0("siha_cross_score",protrecscoreset,"_5_",mode))
    
    resultt1 = get_result_sce_d_siha(siha_ddaprotrec10,protrecscoreset,mode,10)
    print(paste0("siha_cross_score",protrecscoreset,"_10_",mode))
    
    resultt2 = get_result_sce_d_siha(siha_ddaprotrec20,protrecscoreset,mode,20)
    print(paste0("siha_cross_score",protrecscoreset,"_20_",mode))
    
    resultt3 = get_result_sce_d_siha(siha_ddaprotrec30,protrecscoreset,mode,30)
    print(paste0("siha_cross_score",protrecscoreset,"_30_",mode))
    
    resultt4 = get_result_sce_d_siha(siha_ddaprotrec40,protrecscoreset,mode,40)
    print(paste0("siha_cross_score",protrecscoreset,"_40_",mode))
    
    resultt5 = get_result_sce_d_siha(siha_ddaprotrec50,protrecscoreset,mode,50)
    print(paste0("siha_cross_score",protrecscoreset,"_50_",mode))
    
    resultt6 = get_result_sce_d_siha(siha_ddaprotrec200,protrecscoreset,mode,200)
    print(paste0("siha_cross_score",protrecscoreset,"_200_",mode))
    
    print(protrecscoreset)
  }
}

#------sce e-------------------------
## sce e rcc----------
get_result_sce_e_rcc <- function(rc_cprotrec,protrecscoreset,mode){
  RC_recov_protrec_CT1CT2 <- c()
  for (k in c(1:6))  #loop for CT1 -> CT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[k,], rownames(RC_C),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_C)[(RC_C[,k])!=0],RC_tpu,complex_vector,protrecscoreset)
  RC_recov_protrec_CT1CT2<- rbind(RC_recov_protrec_CT1CT2,val[1:4])
  print(k)
  }  
  
  RC_recov_protrec_CT2CT1 <- c()
  for (l in c(7:12))  #loop for CT2 -> CT2_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[l,], rownames(RC_C),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_C)[(RC_C[,l])!=0],RC_tpu,complex_vector,protrecscoreset)
  RC_recov_protrec_CT2CT1<- rbind(RC_recov_protrec_CT2CT1,val[1:4])   
  print(l)
  }  
  
  
  RC_recov_cross_protrec<- data.frame(RC_recov_protrec_CT1CT2,RC_recov_protrec_CT2CT1)  
  
  colnames(RC_recov_cross_protrec)<- c(rep("CT1->CT2",4),rep("CT2->CT1",4))
  
  return(RC_recov_cross_protrec)
}

pathh="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_dec/CODE_AND_DATA/output/"

for(protrecscoreset in c(0.4,0.5,0.6,0.7,0.8,0.95)) 
{
  for(mode in c("max","min","median","mean"))
  {
    result = get_result_sce_e_rcc(rc_cprotrec,protrecscoreset,mode)
    write.csv(result,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_5_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",5,",mode))
    
    result1 = get_result_sce_e_rcc(rc_cprotrec10,protrecscoreset,mode)
    write.csv(result1,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_10_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",10,",mode))
    
    result2 = get_result_sce_e_rcc(rc_cprotrec20,protrecscoreset,mode)
    write.csv(result2,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_20_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",20,",mode))
    
    result3 = get_result_sce_e_rcc(rc_cprotrec30,protrecscoreset,mode)
    write.csv(result3,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_30_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",30,",mode))
    
    result4 = get_result_sce_e_rcc(rc_cprotrec40,protrecscoreset,mode)
    write.csv(result4,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_40_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",40,",mode))
    
    result5 = get_result_sce_e_rcc(rc_cprotrec50,protrecscoreset,mode)
    write.csv(result5,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_50_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",50,",mode))
    
    result6 = get_result_sce_e_rcc(rc_cprotrec200,protrecscoreset,mode)
    write.csv(result6,paste0(pathh,"cross_recoveryrate/result_rc_cross_score",protrecscoreset,"_200_",mode,".csv"))
    print(paste0("RC,cross_recoveryrate,",protrecscoreset,",200,",mode))
  }
}

## sce e lc------
get_result_sce_e_lc <- function(lc_tprotrec,protrecscoreset,mode){
  
  LC_recov_protrec_CT1CT2 <- c()
  for (k in c(1:19))  #loop for CT1 -> CT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-lc_tprotrec[k,], rownames(LC_T),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(LC_T)[(LC_T[,k])!=0],Liver_tpu,complex_vector,protrecscoreset)
  LC_recov_protrec_CT1CT2<- rbind(LC_recov_protrec_CT1CT2,val[1:4])
  print(k)
  }  
  
  LC_recov_protrec_CT2CT1 <- c()
  for (l in c(20:38))  #loop for CT2 -> CT2_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-lc_tprotrec[l,], rownames(LC_T),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(LC_T)[(LC_T[,l])!=0],Liver_tpu,complex_vector,protrecscoreset)
  LC_recov_protrec_CT2CT1<- rbind(LC_recov_protrec_CT2CT1,val[1:4])   
  print(l)
  }  
  
  
  LC_recov_cross_protrec<- data.frame(LC_recov_protrec_CT1CT2,LC_recov_protrec_CT2CT1)  
  
  colnames(LC_recov_cross_protrec)<- c(rep("CT1->CT2",4),rep("CT2->CT1",4))
  
  return(LC_recov_cross_protrec)
}

pathh="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_dec/CODE_AND_DATA/output/cross_recoveryrate/corum2018/LC/"

for(protrecscoreset in c(0.4,0.5,0.6,0.7,0.8,0.95)) 
{
  for(mode in c("max","min","median","mean"))
  {
    result = get_result_sce_e_lc(lc_tprotrec,protrecscoreset,mode)
    write.csv(result,paste0(pathh,"result_lc_cross_score",protrecscoreset,"_5_",mode,".csv"))
    print(paste0("LC,cross_recoveryrate,",protrecscoreset,",5,",mode))
    
    result1 = get_result_sce_e_lc(lc_tprotrec10,protrecscoreset,mode)
    write.csv(result1,paste0(pathh,"result_lc_cross_score",protrecscoreset,"_10_",mode,".csv"))
    print(paste0("LC,cross_recoveryrate,",protrecscoreset,",10,",mode))
    
    result2 = get_result_sce_e_lc(lc_tprotrec20,protrecscoreset,mode)
    write.csv(result2,paste0(pathh,"result_lc_cross_score",protrecscoreset,"_20_",mode,".csv"))
    print(paste0("LC,cross_recoveryrate,",protrecscoreset,",20,",mode))
    
    result3 = get_result_sce_e_lc(lc_tprotrec30,protrecscoreset,mode)
    write.csv(result3,paste0(pathh,"result_lc_cross_score",protrecscoreset,"_30_",mode,".csv"))
    print(paste0("LC,cross_recoveryrate,",protrecscoreset,",30,",mode))
    
    result4 = get_result_sce_e_lc(lc_tprotrec40,protrecscoreset,mode)
    write.csv(result4,paste0(pathh,"result_lc_cross_score",protrecscoreset,"_40_",mode,".csv"))
    print(paste0("LC,cross_recoveryrate,",protrecscoreset,",40,",mode))
    
    result5 = get_result_sce_e_lc(lc_tprotrec50,protrecscoreset,mode)
    write.csv(result5,paste0(pathh,"result_lc_cross_score",protrecscoreset,"_50_",mode,".csv"))
    print(paste0("LC,cross_recoveryrate,",protrecscoreset,",50,",mode))
    
    result6 = get_result_sce_e_lc(lc_tprotrec200,protrecscoreset,mode)
    write.csv(result6,paste0(pathh,"result_lc_cross_score",protrecscoreset,"_200_",mode,".csv"))
    print(paste0("LC,cross_recoveryrate,",protrecscoreset,",200,",mode))
  }
}

## sce e hcc------
get_result_sce_e_hcc <- function(hcc_protrec,protrecscoreset,mode){
  
  HCC_recov_protrec_CT1CT2 <- c()
  for (k in c(1:15))  #loop for CT1 -> CT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-hcc_protrec[k,], rownames(HCC),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)
  
  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(HCC)[(HCC[,k])!=0],Liver_tpu,complex_vector,protrecscoreset)
  HCC_recov_protrec_CT1CT2<- rbind(HCC_recov_protrec_CT1CT2,val[1:4])
  print(k)
  }  
  
  HCC_recov_cross_protrec<- data.frame(HCC_recov_protrec_CT1CT2)  
  
  colnames(HCC_recov_cross_protrec)<- c(rep("CCLE veri",4))
  
  return(HCC_recov_cross_protrec)
}

pathh="W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_dec/CODE_AND_DATA/output/cross_recoveryrate/corum2018/HCC/"

for(protrecscoreset in c(0.4,0.5,0.6,0.7,0.8,0.95)) 
{
  for(mode in c("max","min","median","mean"))
  {
    result = get_result_sce_e_hcc(hcc_protrec,protrecscoreset,mode)
    write.csv(result,paste0(pathh,"result_hcc_cross_score",protrecscoreset,"_5_",mode,".csv"))
    print(paste0("HCC,cross_recoveryrate,",protrecscoreset,",5,",mode))
    
    result1 = get_result_sce_e_hcc(hcc_protrec10,protrecscoreset,mode)
    write.csv(result1,paste0(pathh,"result_hcc_cross_score",protrecscoreset,"_10_",mode,".csv"))
    print(paste0("HCC,cross_recoveryrate,",protrecscoreset,",10,",mode))
    
    result2 = get_result_sce_e_hcc(hcc_protrec20,protrecscoreset,mode)
    write.csv(result2,paste0(pathh,"result_hcc_cross_score",protrecscoreset,"_20_",mode,".csv"))
    print(paste0("HCC,cross_recoveryrate,",protrecscoreset,",20,",mode))
    
    result3 = get_result_sce_e_hcc(hcc_protrec30,protrecscoreset,mode)
    write.csv(result3,paste0(pathh,"result_hcc_cross_score",protrecscoreset,"_30_",mode,".csv"))
    print(paste0("HCC,cross_recoveryrate,",protrecscoreset,",30,",mode))
    
    result4 = get_result_sce_e_hcc(hcc_protrec40,protrecscoreset,mode)
    write.csv(result4,paste0(pathh,"result_hcc_cross_score",protrecscoreset,"_40_",mode,".csv"))
    print(paste0("HCC,cross_recoveryrate,",protrecscoreset,",40,",mode))
    
    result5 = get_result_sce_e_hcc(hcc_protrec50,protrecscoreset,mode)
    write.csv(result5,paste0(pathh,"result_hcc_cross_score",protrecscoreset,"_50_",mode,".csv"))
    print(paste0("HCC,cross_recoveryrate,",protrecscoreset,",50,",mode))
    
    result6 = get_result_sce_e_hcc(hcc_protrec200,protrecscoreset,mode)
    write.csv(result6,paste0(pathh,"result_hcc_cross_score",protrecscoreset,"_200_",mode,".csv"))
    print(paste0("HCC,cross_recoveryrate,",protrecscoreset,",200,",mode))
  }
}


#-----FCS filtering----------
rm(list=ls()) # remove all objects in the env/workspace

## methods--------
countgraph <- function(x,y,z, cplx, g) 
  #where x is the probability of proteins observed and in sample 1, y is the names in sample 1, 
  #z is the names in sample 2, cplx information and g is the set of identified proteins in sample 1
  # x is predicted prot' protrec values, y is names of predicted prot, z is check list, g is original list
{
  output <- c()
  names(x) <- y 
  levels <- seq(0, 1, by=0.1) #the distribution of probability levels
  for (i in 1:length(levels)) #from 0 to 1
  {
    total <- length(names(x)[x >= levels[i] & x < (levels[i] + 0.1)]) # total number of predicted prot within this range
    original <- length(intersect(names(x)[x >= levels[i] & x < (levels[i] + 0.1)], g)) # how many were in original screen
    validated <- length(intersect(names(x)[x >= levels[i] & x < (levels[i] + 0.1)], setdiff(z, g))) # validated prot
    
    total <- total - original - validated # unvalidated
    
    output <- rbind(output, cbind(levels[i], total, validated, original))
  }
  return(output)
}



# compared to the countgraph, this is more like a cumulative data
precision_recall_2 <- function(x,y,z, cplx, g) 
  #where x is the probability of pred proteins, y is the names of pred proteins, z is the names in cross batch check list
{
  output <- c()
  names(x) <- y 
  levels <- seq(0, 1, by=0.1) #the distribution of probability levels
  for (i in 1:length(levels)) #from 0 to 1
  {
    total <- length(names(x)[x >= (levels[i])])
    original <- length(intersect(names(x)[x >= (levels[i])], g))
    validated <- length(intersect(names(x)[x >= (levels[i])], setdiff(z, g)))
    
    total <- total - original - validated # unvalidated
    
    output <- rbind(output, cbind(levels[i], total, validated, original))
  }
  return(output)
}


#Take the FCS significant complexes, rank them by the tibtech approach
get_sig_complex <- function(replicate){
  FCS_sig_complex_rep <- as.vector(replicate[which(replicate <= 0.05)])
  return(FCS_sig_complex_rep)
}

get_mat_n <-function(replicate){
  fcs_normal_mat_rep <- fcs_normal_mat[which(rownames(fcs_normal_mat) %in% unique(unlist(complex_vector[names(replicate)]))),]
  #                                         names of proteins in fcs results                          names of complexes in replicate
  #                                                                        get the constituent proteins of these complexes
  #                                    get those proteins in both fcs results and those that make up complexes in the sample
  return(fcs_normal_mat_rep)
}

get_mat_c <-function(replicate){
  fcs_cancer_mat_rep <- fcs_cancer_mat[which(rownames(fcs_cancer_mat) %in% unique(unlist(complex_vector[names(replicate)]))),]
  #                                         names of proteins in fcs results                          names of complexes in replicate
  #                                                                        get the constituent proteins of these complexes
  #                                    get those proteins in both fcs results and those that make up complexes in the sample
  return(fcs_cancer_mat_rep)
}

get_protprob_normal_mat <- function(replicate){
  protprob_normal_mat_rep <- normal_mat[rownames(normal_mat)%in%unique(unlist(complex_vector[names(replicate)])),]
  return(protprob_normal_mat_rep)
}

get_protprob_cancer_mat <- function(replicate){
  protprob_cancer_mat_rep <- cancer_mat[rownames(cancer_mat)%in%unique(unlist(complex_vector[names(replicate)])),]
  return(protprob_cancer_mat_rep)
}

change_name <- function(replicate){
  rownames(replicate) <- replicate[,1]
  replicate <- replicate[,-1]
}

order_stuff <- function(replicate){
  return(replicate[order(replicate[,1], decreasing = T),])
}

## fcs filtering-----------------
get_average_values_n <- function(FCS_sig_complex_n){
  # bar plot
  fcs_n_total <- c()
  protrec_n_total <- c()
  
  fcs_n_validated <- c()
  protrec_n_validated <- c()
  
  fcs_n_original <- c()
  protrec_n_original <- c()
  
  # line plot
  ave_filtered_n_v1 <- c()
  ave_unfiltered_n_v1 <- c()
  
  ave_filtered_n_total <- c()
  ave_unfiltered_n_total <- c()
  
  ave_filtered_n_validated <- c()
  ave_unfiltered_n_validated <- c()
  
  ave_filtered_n_original <- c()
  ave_unfiltered_n_original <- c()
  
  # data
  for(i in c(1:6)){
    # for fsc_normal_mat:
    # c1t1 = 1, c2t1 = 2, ..., c1t2 = 7, c2t2 = 8, ..., c6t2 = 12
    # n1t1 = 1, n2t1 = 2, ..., n1t2 = 7, n2t2 = 8, ..., n6t2 = 12
    
    # rc_peptides_unique and rc:
    # c1t1 = 1+6, c2t1 = 2+6, ..., c6t1 = 6+6
    # c1t2 = 1+18, c2t2 =2+18, ..., c6t2 = 6+18
    # n1t1 = 1, n2t1 = 2, ..., n6t1 = 6
    # n1t2 = 1+12, n2t2 = 2+12, ..., n6t2 = 6+12
    
    # common data
    protprob_normal_mat <- get_protprob_normal_mat(get_sig_complex(FCS_sig_complex_n[i,]))
    fcs_normal_mat <- get_mat_n(get_sig_complex(FCS_sig_complex_n[i,]))
    
    protprob_normal_mat2 <- get_protprob_normal_mat(get_sig_complex(FCS_sig_complex_n[i+6,]))
    fcs_normal_mat2 <- get_mat_n(get_sig_complex(FCS_sig_complex_n[i+6,]))
    
    # bar plot data
    fcsT12 <- change_name(countgraph(fcs_normal_mat[,i], rownames(fcs_normal_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+12])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i]))]))
    protprobT12 <- change_name(countgraph(protprob_normal_mat[,i], rownames(protprob_normal_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+12])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i]))]))
    fcsT21 <- change_name(countgraph(fcs_normal_mat2[,i+6], rownames(fcs_normal_mat2),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+12])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+12]))]))
    protprobT21 <- change_name(countgraph(protprob_normal_mat2[,i+6], rownames(protprob_normal_mat2),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+12])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+12]))]))
    
    # line plot data
    protprob_T12_fcs_filtered <- order_stuff(precision_recall_2(protprob_normal_mat[,i], rownames(protprob_normal_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+12])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i]))]))
    protprob_T12_unfiltered <- order_stuff(precision_recall_2(normal_mat[,i], rownames(normal_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+12])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i]))]))
    protprob_T21_fcs_filtered <- order_stuff(precision_recall_2(protprob_normal_mat2[,i+6], rownames(protprob_normal_mat2),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+12])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+12]))]))
    protprob_T21_unfiltered <- order_stuff(precision_recall_2(normal_mat[,i+6], rownames(normal_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+12])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+12]))]))
    
    # bar plot results
    fcs_n_total <- cbind(fcs_n_total,fcsT12[,1],fcsT21[,1])
    protrec_n_total <- cbind(protrec_n_total,protprobT12[,1],protprobT21[,1])
    
    fcs_n_validated <- cbind(fcs_n_validated,fcsT12[,2],fcsT21[,2])
    protrec_n_validated <- cbind(protrec_n_validated,protprobT12[,2],protprobT21[,2])
    
    fcs_n_original <- cbind(fcs_n_original,fcsT12[,3],fcsT21[,3])
    protrec_n_original <- cbind(protrec_n_original,protprobT12[,3],protprobT21[,3])
    
    # line plot results
    ave_filtered_n_v1 <- cbind(ave_filtered_n_v1,protprob_T12_fcs_filtered[,1],protprob_T21_fcs_filtered[,1])
    ave_unfiltered_n_v1 <- cbind(ave_unfiltered_n_v1,protprob_T12_unfiltered[,1],protprob_T21_unfiltered[,1])
    
    ave_filtered_n_total <- cbind(ave_filtered_n_total,protprob_T12_fcs_filtered[,2],protprob_T21_fcs_filtered[,2])
    ave_unfiltered_n_total <- cbind(ave_unfiltered_n_total,protprob_T12_unfiltered[,2],protprob_T21_unfiltered[,2])
    
    ave_filtered_n_validated <- cbind(ave_filtered_n_validated,protprob_T12_fcs_filtered[,3],protprob_T21_fcs_filtered[,3])
    ave_unfiltered_n_validated <- cbind(ave_unfiltered_n_validated,protprob_T12_unfiltered[,3],protprob_T21_unfiltered[,3])
    
    ave_filtered_n_original <- cbind(ave_filtered_n_original,protprob_T12_fcs_filtered[,4],protprob_T21_fcs_filtered[,4])
    ave_unfiltered_n_original <- cbind(ave_unfiltered_n_original,protprob_T12_unfiltered[,4],protprob_T21_unfiltered[,4])
  }
  
  # bar plot average
  ave_fcs_n <- c()
  ave_fcs_n <- cbind(ave_fcs_n, rowMeans(fcs_n_total),rowMeans(fcs_n_validated),rowMeans(fcs_n_original))
  colnames(ave_fcs_n) <- c('total','validated','original')
  
  ave_protrec_n <- c()
  ave_protrec_n <- cbind(ave_protrec_n, rowMeans(protrec_n_total),rowMeans(protrec_n_validated),rowMeans(protrec_n_original))
  colnames(ave_protrec_n) <- c('total','validated','original')
  
  # line plot average
  ave_filtered_n <- c()
  ave_filtered_n <- cbind(ave_filtered_n, rowMeans(ave_filtered_n_v1),rowMeans(ave_filtered_n_total),rowMeans(ave_filtered_n_validated),rowMeans(ave_filtered_n_original))
  colnames(ave_filtered_n) <- c('v1','total','validated','original')
  ave_unfiltered_n <- c()
  ave_unfiltered_n <- cbind(ave_unfiltered_n, rowMeans(ave_unfiltered_n_v1),rowMeans(ave_unfiltered_n_total),rowMeans(ave_unfiltered_n_validated),rowMeans(ave_unfiltered_n_original))
  colnames(ave_unfiltered_n) <- c('v1','total','validated','original')
  
  
  # plot bar plot
  pdf("W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_dec/CODE_AND_DATA_organized/output/try/fcs_sig_complexes_protprob_methods_barplot_rcn.pdf",width=8,height=6)
  par(mfrow=c(1,2))
  barplot(t(ave_fcs_n), col=c("paleturquoise", "royalblue2", "dodgerblue4"), main="FCS (RCN)", ylab="Count", xlab="Prob(Protein p is present but unreported)", cex.lab=1.2, font=2)
  barplot(t(ave_protrec_n), col=c("paleturquoise", "royalblue2", "dodgerblue4"), main="PROTREC (RCN)", ylab="Count", xlab="Prob(Protein p is present but unreported)", cex.lab=1.2,font=2)
  
  # legend
  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(-0.45,0.9,legend = c("Original", "Validated", "Unvalidated"), col= c("dodgerblue4", "royalblue2","paleturquoise"), lwd = 10, xpd = TRUE, horiz = TRUE, cex = 0.8, seg.len=1, bty = 'n')
  dev.off()
  
  # plot line plot
  pdf("W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_dec/CODE_AND_DATA_organized/output/try/protprob_filteredvsunfiltered_rcn.pdf",width=9,height=8)
  par(mar = c(5, 5, 12, 5))
  plot(ave_filtered_n[,1], ave_filtered_n[,3]/(ave_filtered_n[,2] + ave_filtered_n[,3]), type="o", col="purple", ylim=c(0, 1), xlab="PROTREC score", ylab="% of valdiated proteins",cex.lab=1.7)
  lines(ave_unfiltered_n[,1], ave_unfiltered_n[,3]/(ave_unfiltered_n[,2] + ave_unfiltered_n[,3]), type="o", col="blue")
  par(new=T)
  
  plot(ave_unfiltered_n[,1], (ave_unfiltered_n[,2] + ave_unfiltered_n[,3]), axes=F, xlab=NA, ylab=NA, col="blue", pch=2)
  mtext('# of proteins with at least x PROTREC score', side=4, line=3,  outer=F, cex=1.7)
  axis(side=4)
  par(new=T)
  
  plot(ave_filtered_n[,1], (ave_filtered_n[,2] + ave_filtered_n[,3]), axes=F, xlab=NA, ylab=NA, col="purple", pch=2, ylim=range((ave_unfiltered_n[,2] + ave_unfiltered_n[,3]), (ave_filtered_n[,2] + ave_filtered_n[,3])))
  # use this title for supplementary figures but not main text figure
  title(main = "PROTREC performance w/ and w/o \nFCS filtering (RCN)", 
        cex.main = 2, line = 8)
  legend("topleft",xpd=TRUE,inset=c(0.005,-0.3),c("non-filtered # of validated proteins","fcs-filtered # of validated proteins"), lty=c(1,1), pch=1, col=c("blue", "purple"), cex=1)
  legend("topleft",xpd=TRUE,inset=c(0.005,-0.15),c("non-filtered # of proteins with at least x PROTREC score","fcs-filtered # of proteins with at least x PROTREC score"), lty=c(0,0), pch=2, col=c("blue", "purple"), cex=1)
  dev.off()
}
get_average_values_c <- function(FCS_sig_complex_c){
  # bar plot
  fcs_c_total <- c()
  protrec_c_total <- c()
  
  fcs_c_validated <- c()
  protrec_c_validated <- c()
  
  fcs_c_original <- c()
  protrec_c_original <- c()
  
  # line plot
  ave_filtered_c_v1 <- c()
  ave_unfiltered_c_v1 <- c()
  
  ave_filtered_c_total <- c()
  ave_unfiltered_c_total <- c()
  
  ave_filtered_c_validated <- c()
  ave_unfiltered_c_validated <- c()
  
  ave_filtered_c_original <- c()
  ave_unfiltered_c_original <- c()
  
  # data
  for(i in c(1:6)){
    # for fsc_cancer_mat:
    # c1t1 = 1, c2t1 = 2, ..., c1t2 = 7, c2t2 = 8, ..., c6t2 = 12
    # n1t1 = 1, n2t1 = 2, ..., n1t2 = 7, n2t2 = 8, ..., n6t2 = 12
    
    # rc_peptides_unique and rc:
    # c1t1 = 1+6, c2t1 = 2+6, ..., c6t1 = 6+6
    # c1t2 = 1+18, c2t2 =2+18, ..., c6t2 = 6+18
    # n1t1 = 1, n2t1 = 2, ..., n6t1 = 6
    # n1t2 = 1+12, n2t2 = 2+12, ..., n6t2 = 6+12
    
    # separate out the samples and filter out those with 1-prob <= 0.05, which means probability > 0.95
    # check which complexes are present in the sample, get the constituent proteins of these complexes
    # get the proteins in all samples from protrec results
    # get the intersection of the two sets of proteins
    protprob_cancer_mat <- get_protprob_cancer_mat(get_sig_complex(FCS_sig_complex_c[i,]))
    fcs_cancer_mat <- get_mat_c(get_sig_complex(FCS_sig_complex_c[i,]))
    
    protprob_cancer_mat2 <- get_protprob_cancer_mat(get_sig_complex(FCS_sig_complex_c[i+6,]))
    fcs_cancer_mat2 <- get_mat_c(get_sig_complex(FCS_sig_complex_c[i+6,]))
    
    # bar plot data
    fcsT12 <- change_name(countgraph(fcs_cancer_mat[,i], rownames(fcs_cancer_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+18])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+6])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+6]))]))
    protprobT12 <- change_name(countgraph(protprob_cancer_mat[,i], rownames(protprob_cancer_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+18])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+6])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+6]))]))
    fcsT21 <- change_name(countgraph(fcs_cancer_mat2[,i+6], rownames(fcs_cancer_mat2),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+18])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+6])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+18]))]))
    protprobT21 <- change_name(countgraph(protprob_cancer_mat2[,i+6], rownames(protprob_cancer_mat2),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+18])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+6])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+18]))]))
    
    # line plot data
    protprob_T12_fcs_filtered <- order_stuff(precision_recall_2(protprob_cancer_mat[,i], rownames(protprob_cancer_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+18])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+6])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+6]))]))
    protprob_T12_unfiltered <- order_stuff(precision_recall_2(cancer_mat[,i], rownames(cancer_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+18])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+6])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+6]))]))
    protprob_T21_fcs_filtered <- order_stuff(precision_recall_2(protprob_cancer_mat2[,i+6], rownames(protprob_cancer_mat2),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+18])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+6])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+18]))]))
    protprob_T21_unfiltered <- order_stuff(precision_recall_2(cancer_mat[,i+6], rownames(cancer_mat),union(rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+18])], rownames(RC_peptides_uniq)[!is.na(RC_peptides_uniq[,i+6])]),unique(as.vector(unlist(complex_vector))), rownames(RC)[which(!is.na(RC[,i+18]))]))
    
    # bar plot results
    fcs_c_total <- cbind(fcs_c_total,fcsT12[,1],fcsT21[,1])
    protrec_c_total <- cbind(protrec_c_total,protprobT12[,1],protprobT21[,1])
    
    fcs_c_validated <- cbind(fcs_c_validated,fcsT12[,2],fcsT21[,2])
    protrec_c_validated <- cbind(protrec_c_validated,protprobT12[,2],protprobT21[,2])
    
    fcs_c_original <- cbind(fcs_c_original,fcsT12[,3],fcsT21[,3])
    protrec_c_original <- cbind(protrec_c_original,protprobT12[,3],protprobT21[,3])
    
    # line plot results
    ave_filtered_c_v1 <- cbind(ave_filtered_c_v1,protprob_T12_fcs_filtered[,1],protprob_T21_fcs_filtered[,1])
    ave_unfiltered_c_v1 <- cbind(ave_unfiltered_c_v1,protprob_T12_unfiltered[,1],protprob_T21_unfiltered[,1])
    
    ave_filtered_c_total <- cbind(ave_filtered_c_total,protprob_T12_fcs_filtered[,2],protprob_T21_fcs_filtered[,2])
    ave_unfiltered_c_total <- cbind(ave_unfiltered_c_total,protprob_T12_unfiltered[,2],protprob_T21_unfiltered[,2])
    
    ave_filtered_c_validated <- cbind(ave_filtered_c_validated,protprob_T12_fcs_filtered[,3],protprob_T21_fcs_filtered[,3])
    ave_unfiltered_c_validated <- cbind(ave_unfiltered_c_validated,protprob_T12_unfiltered[,3],protprob_T21_unfiltered[,3])
    
    ave_filtered_c_original <- cbind(ave_filtered_c_original,protprob_T12_fcs_filtered[,4],protprob_T21_fcs_filtered[,4])
    ave_unfiltered_c_original <- cbind(ave_unfiltered_c_original,protprob_T12_unfiltered[,4],protprob_T21_unfiltered[,4])
  }
  
  # bar plot average
  ave_fcs_c <- c()
  ave_fcs_c <- cbind(ave_fcs_c, rowMeans(fcs_c_total),rowMeans(fcs_c_validated),rowMeans(fcs_c_original))
  colnames(ave_fcs_c) <- c('total','validated','original')
  
  ave_protrec_c <- c()
  ave_protrec_c <- cbind(ave_protrec_c, rowMeans(protrec_c_total),rowMeans(protrec_c_validated),rowMeans(protrec_c_original))
  colnames(ave_protrec_c) <- c('total','validated','original')
  
  # line plot average
  ave_filtered_c <- c()
  ave_filtered_c <- cbind(ave_filtered_c, rowMeans(ave_filtered_c_v1),rowMeans(ave_filtered_c_total),rowMeans(ave_filtered_c_validated),rowMeans(ave_filtered_c_original))
  colnames(ave_filtered_c) <- c('v1','total','validated','original')
  ave_unfiltered_c <- c()
  ave_unfiltered_c <- cbind(ave_unfiltered_c, rowMeans(ave_unfiltered_c_v1),rowMeans(ave_unfiltered_c_total),rowMeans(ave_unfiltered_c_validated),rowMeans(ave_unfiltered_c_original))
  colnames(ave_unfiltered_c) <- c('v1','total','validated','original')
  
  
  # plot bar plot
  pdf("W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_dec/CODE_AND_DATA_organized/output/try/fcs_sig_complexes_protprob_methods_barplot_rcc.pdf",width=8,height=6)
  par(mfrow=c(1,2))
  barplot(t(ave_fcs_c), col=c("paleturquoise", "royalblue2", "dodgerblue4"), main="FCS (RCC)", ylab="Count", xlab="Prob(Protein p is present but unreported)", cex.lab=1.2, font=2)
  barplot(t(ave_protrec_c), col=c("paleturquoise", "royalblue2", "dodgerblue4"), main="PROTREC (RCC)", ylab="Count", xlab="Prob(Protein p is present but unreported)", cex.lab=1.2,font=2)
  
  # legend
  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(-0.45,0.9,legend = c("Original", "Validated", "Unvalidated"), col= c("dodgerblue4", "royalblue2","paleturquoise"), lwd = 10, xpd = TRUE, horiz = TRUE, cex = 0.8, seg.len=1, bty = 'n')
  dev.off()
  
  # plot line plot
  pdf("W:/AAAANot_in_WD/SBS_proj/methodx/new_draft_dec/CODE_AND_DATA_organized/output/try/protprob_filteredvsunfiltered_rcc.pdf",width=9,height=8)
  par(mar = c(5, 5, 12, 5))
  plot(ave_filtered_c[,1], ave_filtered_c[,3]/(ave_filtered_c[,2] + ave_filtered_c[,3]), type="o", col="purple", ylim=c(0, 1), xlab="PROTREC score", ylab="% of valdiated proteins",cex.lab=1.7)
  lines(ave_unfiltered_c[,1], ave_unfiltered_c[,3]/(ave_unfiltered_c[,2] + ave_unfiltered_c[,3]), type="o", col="blue")
  par(new=T)
  
  plot(ave_unfiltered_c[,1], (ave_unfiltered_c[,2] + ave_unfiltered_c[,3]), axes=F, xlab=NA, ylab=NA, col="blue", pch=2)
  mtext('# of proteins with at least x PROTREC score', side=4, line=3,  outer=F, cex=1.7)
  axis(side=4)
  par(new=T)
  
  plot(ave_filtered_c[,1], (ave_filtered_c[,2] + ave_filtered_c[,3]), axes=F, xlab=NA, ylab=NA, col="purple", pch=2, ylim=range((ave_unfiltered_c[,2] + ave_unfiltered_c[,3]), (ave_filtered_c[,2] + ave_filtered_c[,3])))
  # use this title for supplementary figures but not main text figure
  title(main = "PROTREC performance w/ and w/o \nFCS filtering (RCC)", 
        cex.main = 2, line = 8)
  legend("topleft",xpd=TRUE,inset=c(0.005,-0.3),c("non-filtered # of validated proteins","fcs-filtered # of validated proteins"), lty=c(1,1), pch=1, col=c("blue", "purple"), cex=1)
  legend("topleft",xpd=TRUE,inset=c(0.005,-0.15),c("non-filtered # of proteins with at least x PROTREC score","fcs-filtered # of proteins with at least x PROTREC score"), lty=c(0,0), pch=2, col=c("blue", "purple"), cex=1)
  dev.off()
}

get_average_values_n(FCS_sig_complex_n)
get_average_values_c(FCS_sig_complex_c)