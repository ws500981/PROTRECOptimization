###DATA######
#' @title This is the RC protein data included in my package
#' @name RC
#' @docType data
#'
#'
#'
NULL

#' @title This is the RC cancer protein data included in my package
#' @name RC_C
#' @docType data
#'
#'
#'
NULL

#' @title This is the RC normal protein data included in my package
#' @name RC_N
#' @docType data
#'
#'
#'
NULL

#' @title This is the RC peptide data included in my package
#' @name RC_peptides_uniq
#' @docType data
#'
#'
#'
NULL

#' @title This is the complex_vector included in my package
#' @name complex_vector
#' @docType data
#' @references \url{http://mips.helmholtz-muenchen.de/corum/}
#'
#'
#'
NULL

#' @title This is the RC_ccle dataset included in my package
#' @name RC_ccle
#' @docType data
#' @references \url{https://pubmed.ncbi.nlm.nih.gov/31978347/}
#'
#'
#'
NULL

#' @title This is the Liver_ccle dataset included in my package
#' @name Liver_ccle
#' @docType data
#' @references \url{https://pubmed.ncbi.nlm.nih.gov/31978347/}
#'
#'
#'
NULL

#' @title This is the Hela_ddaProteins dataset included in my package
#' @name Hela_ddaProteins
#' @docType data
#'
#'
#'
NULL

#' @title This is the Siha_ddaProteins dataset included in my package
#' @name Siha_ddaProteins
#' @docType data
#'
#'
#'
NULL

#' @title This is the heladdapepuniq dataset included in my package
#' @name heladdapepuniq
#' @docType data
#'
#'
#'
NULL

#' @title This is the sihaddapepuniq dataset included in my package
#' @name sihaddapepuniq
#' @docType data
#'
#'
#'
NULL

#' @title This is the hela_dia_prot_list1 dataset included in my package
#' @name hela_dia_prot_list1
#' @docType data
#'
#'
#'
NULL

#' @title This is the hela_dia_prot_list2 dataset included in my package
#' @name hela_dia_prot_list2
#' @docType data
#'
#'
#'
NULL

#' @title This is the siha_dia_prot_list1 dataset included in my package
#' @name siha_dia_prot_list1
#' @docType data
#'
#'
#'
NULL

#' @title This is the siha_dia_prot_list2 dataset included in my package
#' @name siha_dia_prot_list2
#' @docType data
#'
#'
#'
NULL

#' @title This is the siha_dia_prot_list3 dataset included in my package
#' @name siha_dia_prot_list3
#' @docType data
#'
#'
#'
NULL

#' @title This is the LC_T dataset included in my package
#' @name LC_T
#' @docType data
#'
#'
#'
NULL

#' @title This is the HCC dataset included in my package
#' @name HCC
#' @docType data
#'
#'
#'
NULL

#' @title This is the fcs_normal_mat dataset included in my package
#' @name fcs_normal_mat
#' @docType data
#'
#'
#'
NULL

#' @title This is the fcs_cancer_mat dataset included in my package
#' @name fcs_cancer_mat
#' @docType data
#'
#'
#'
NULL

#' @title This is the protprob_normal_mat dataset included in my package
#' @name protprob_normal_mat
#' @docType data
#'
#'
#'
NULL

#' @title This is the protprob_cancer_mat dataset included in my package
#' @name protprob_cancer_mat
#' @docType data
#'
#'
#'
NULL

#' @title This is the FCS_sig_complex_n dataset included in my package
#' @name FCS_sig_complex_n
#' @docType data
#'
#'
#'
NULL

#' @title This is the FCS_sig_complex_c dataset included in my package
#' @name FCS_sig_complex_c
#' @docType data
#'
#'
#'
NULL


###PROTREC functions########
#' @title This assigns probabilities to complexes based on PROTREC
#' @name PROTREC_cplx_prob
#' @param data Proteomics expression data
#' @param complex protein complex
#' @param fdr false discovery rate of proteomics screen
#' @param threshold the minimum size of the complex required
#' @return output_mat A matrix of PROTREC (1-probabilities) assigned to complexes for the dataset
#' @export
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

#' @title This assigns probabilities to individual proteins based on PROTREC
#' @name PROTREC_protprob
#' @param cplx complex vector
#' @param p is the complex-based PROTREC probability, typically should be the 1 minus the result from certain sample of PROTREC_cplx_prob
#' @param prot is the list of observed proteins in the screen
#' @param fdr is the false discovery rate of the screen
#' @param mode is the mode of PROTREC score selection, default max
#' @return output The vector of PROTREC probabilities for the proteins of a sample
#' @export
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

###Recovery Rate########
#for checking recovery and providing a significance value for the recovery for PROTREC in scenarios a, c, d, and e
#' @title This works out the recovery significance for PROTREC
#' @name pairwise_recovery_protrec
#' @param prot_predict_list A vector of significant proteins for a given sample
#' @param original_prot_list The original set of proteins observed for that particular sample
#' @param check_prot_list A set of proteins observed in a second replicate
#' @param complex_vec A list of complex objects
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @return output A vector of five pieces of information: the observed proportion of overlap,
#' the significance of this overlap (p-value), the number of verified proteins, the total number of predicted missing proteins,
#' and the list of validated proteins separated by 'a'.
#' @export
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

#for checking recovery and providing a significance value for the recovery for FCS,HE and GSEA
#' @title This works out the recovery significance for FCS,HE and GSEA
#' @name pairwise_recovery
#' @param predict_list A vector of significant proteins for a given sample, default set 0.05 p-value as significant
#' @param original_prot_list The original set of proteins observed for that particular sample
#' @param check_prot_list A set of proteins observed in a second replicate
#' @param complex_vec A list of complex objects
#' @return output A vector of five pieces of information: the observed proportion of overlap,
#' the significance of this overlap (p-value), the number of verified proteins, the total number of predicted missing proteins,
#' and the list of validated proteins separated by 'a'.
#' @export
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


#for checking recovery and providing a significance value for the recovery for PROTREC in scenarios b
#' @title This works out the recovery significance for PROTREC
#' @name ntop_recovery_protrec
#' @param prot_predict_list A vector of significant proteins for a given sample
#' @param original_prot_list The original set of proteins observed for that particular sample
#' @param check_prot_list A set of proteins observed in a second replicate
#' @param complex_vec A list of complex objects
#' @param aaa The number of predicted proteins with pval < 0.05 for the particular method
#' @return output A vector of five pieces of information: the observed proportion of overlap,
#' the significance of this overlap (p-value), the number of verified proteins, the total number of predicted missing proteins,
#' and the list of validated proteins separated by 'a'.
#' @export
ntop_recovery_protrec <- function(prot_predict_list, original_prot_list, check_prot_list, complex_vec)
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

###FCS functions########
#' @title This assigns probabilities to complexes based on FCS
#' @name fcs
#' @param data Proteomics expression data
#' @param complex_vector protein complex
#' @param sim_size how many iterations for simulation, recommend 1000
#' @param threshold the minimum size of the complex required
#' @return output_mat A matrix of FCS probabilities assigned to complexes for the dataset
#' @export
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

#' @title This assigns probabilities to individual proteins based on FCS
#' @name fcs_prot_prob_assign
#' @param cplx complex vector
#' @param p is the complex-based FCS probability, typically should be the 1 minus the result from certain sample of fcs
#' @return The vector of FCS probabilities for the proteins of a sample
#' @export
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


###Hypergeometric Enrichment(HE) functions########
#' @title This assigns probabilities to complexes based on HE
#' @name hgtest
#' @param data Proteomics expression data
#' @param complex_vector protein complex
#' @param threshold the minimum size of the complex required
#' @return output_mat A matrix of HE probabilities assigned to complexes for the dataset
#' @export
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

#' @title This assigns probabilities to individual proteins based on HE
#' @name hgtest_prot_prob_assign
#' @param cplx complex vector
#' @param p is the complex-based HE probability, typically should be the 1 minus the result from certain sample of hgtest
#' @return output The vector of HE probabilities for the proteins of a sample
#' @export
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


###Gene Set Enrichment Analysis(GSEA) functions########
#' @title This performs t-test inside GSEA function
#' @name mat.ttest
#' @param data Proteomics expression data
#' @return output_mat A matrix containing all p-values after t-test
#' @export
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

#' @title This assigns probabilities to complexes based on GSEA
#' @name repgsea
#' @param data Proteomics expression data
#' @param complex_vector protein complex
#' @return gsea_pvals A matrix of GSEA probabilities assigned to complexes for the dataset
#' @export
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





###PROTREC hyperparameter tuning########
#' @title This performs hyperparameter tuning of PROTREC based on scenario a using the RC dataset
#' @name get_result_sce_a_rc
#' @param rc_nprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the RCN dataset
#' @param rc_cprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the RCC dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @return output A matrix with 16 columns and 6 rows: each four columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins. The first four columns
#' correspond to RCN dataset using sample 1 (validated using sample 2), the second four columns correspond to
#' RCN dataset using sample 2 (validated using sample 1). The next four columns correspond to RCC dataset
#' using sample 1 (validated using sample 2), the last four columns correspond to RCC dataset using
#' sample 2 (validated using sample 1). The 6 rows represent the 6 technical replicates of each sample.
#' @export
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


#' @title This performs hyperparameter tuning of PROTREC based on scenario a using the Hela dataset
#' @name get_result_sce_a_hela
#' @param hela_ddaprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the Hela dda dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @return output A matrix with 4 columns and 6 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 6 rows represent the 6 permutations of cross-validation using the samples.
#' @export
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

#' @title This performs hyperparameter tuning of PROTREC based on scenario a using the Siha dataset
#' @name get_result_sce_a_siha
#' @param hela_ddaprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the Siha dda dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @return output A matrix with 4 columns and 6 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 6 rows represent the 6 permutations of cross-validation using the samples.
#' @export
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

#' @title This performs hyperparameter tuning of PROTREC based on scenario b using the RCN dataset
#' @name get_result_sce_b_rcn
#' @param mode is the mode of PROTREC score selection, default max
#' @param complex_size is the minimum size of the complexes used
#' @param rc_nfcs A vector of significant proteins for a given sample using FCS, default set 0.05 p-value as significant
#' @param rc_nprotrec A vector of significant proteins for a given sample using PROTREC, default set 0.05 p-value as significant
#' @param rc_nhg A vector of significant proteins for a given sample using HE, default set 0.05 p-value as significant
#' @return output A matrix with 9 columns and 48 rows: each three columns contain information about
#' the number N of predicted proteins based on a given method, the results from the given method, and the results from PROTREC.
#' Each four rows contain information about the observed proportion of overlap (the PROTREC score),
#' the significance of this overlap (p-value), the total number of predicted missing proteins, and the number of
#' validated proteins. The twelve sets of 4 rows that make up the 48 rows correspond to the permutations of
#' cross-validation using the samples.
#' @export
get_result_sce_b_rcn<-function(mode,complex_size,rc_nfcs, rc_nprotrec, rc_nhg){
  print("running...")

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
  return(all1v1)
}

#' @title This performs hyperparameter tuning of PROTREC based on scenario b using the RCC dataset
#' @name get_result_sce_b_rcc
#' @param mode is the mode of PROTREC score selection, default max
#' @param complex_size is the minimum size of the complexes used
#' @param rc_cfcs A vector of significant proteins for a given sample using FCS, default set 0.05 p-value as significant
#' @param rc_cprotrec A vector of significant proteins for a given sample using PROTREC, default set 0.05 p-value as significant
#' @param rc_chg A vector of significant proteins for a given sample using HE, default set 0.05 p-value as significant
#' @return output A matrix with 9 columns and 48 rows: each three columns contain information about
#' the number N of predicted proteins based on a given method, the results from the given method, and the results from PROTREC.
#' Each four rows contain information about the observed proportion of overlap (the PROTREC score),
#' the significance of this overlap (p-value), the total number of predicted missing proteins, and the number of
#' validated proteins. The twelve sets of 4 rows that make up the 48 rows correspond to the permutations of
#' cross-validation using the samples.
#' @export
get_result_sce_b_rcc<-function(mode,complex_size, rc_cfcs, rc_cprotrec, rc_chg){
  print("running...")

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
  return(all1v1)
}

#' @title This performs hyperparameter tuning of PROTREC based on scenario b using the Hela dataset
#' @name get_result_sce_b_hela
#' @param mode is the mode of PROTREC score selection, default max
#' @param complex_size is the minimum size of the complexes used
#' @param heladdaprotrec A vector of significant proteins for a given sample using PROTREC, default set 0.05 p-value as significant
#' @param heladdafcs A vector of significant proteins for a given sample using FCS, default set 0.05 p-value as significant
#' @param heladdahg A vector of significant proteins for a given sample using HE, default set 0.05 p-value as significant
#' @return output A matrix with 9 columns and 24 rows: each three columns contain information about
#' the number N of predicted proteins based on a given method, the results from the given method, and the results from PROTREC.
#' Each four rows contain information about the observed proportion of overlap (the PROTREC score),
#' the significance of this overlap (p-value), the total number of predicted missing proteins, and the number of
#' validated proteins. The six sets of 4 rows that make up the 24 rows correspond to the permutations of
#' cross-validation using the samples.
#' @export
get_result_sce_b_hela<-function(mode,complex_size,heladdaprotrec,heladdafcs,heladdahg){
  print("running...")

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
  return(all1v1)
}

#' @title This performs hyperparameter tuning of PROTREC based on scenario b using the Siha dataset
#' @name get_result_sce_b_siha
#' @param mode is the mode of PROTREC score selection, default max
#' @param complex_size is the minimum size of the complexes used
#' @param sihaddaprotrec A vector of significant proteins for a given sample using PROTREC, default set 0.05 p-value as significant
#' @param sihaddafcs A vector of significant proteins for a given sample using FCS, default set 0.05 p-value as significant
#' @param sihaddahg A vector of significant proteins for a given sample using HE, default set 0.05 p-value as significant
#' @return output A matrix with 9 columns and 24 rows: each three columns contain information about
#' the number N of predicted proteins based on a given method, the results from the given method, and the results from PROTREC.
#' Each four rows contain information about the observed proportion of overlap (the PROTREC score),
#' the significance of this overlap (p-value), the total number of predicted missing proteins, and the number of
#' validated proteins. The six sets of 4 rows that make up the 24 rows correspond to the permutations of
#' cross-validation using the samples.
#' @export
get_result_sce_b_siha<-function(mode,complex_size,sihaddaprotrec,sihaddafcs,sihaddahg){
  print("running...")

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
  return(all1v1)
}



#' @title This performs hyperparameter tuning of PROTREC based on scenario c using the RCN dataset
#' @name get_result_sce_c_rcn
#' @param rc_nprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the RCN dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @param size is the minimum size of the complexes used
#' @return output A matrix with 4 columns and 12 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 12 rows represent the 12 permutations of cross-validation using the samples.
#' @export
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
  print(paste0("result_rcn_twopept_",protrecscoreset,"_",size,"_",mode))
  return(RC_recov_protrec_NT1NT11)
}

#' @title This performs hyperparameter tuning of PROTREC based on scenario c using the RCC dataset
#' @name get_result_sce_c_rcc
#' @param rc_cprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the RCC dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @param size is the minimum size of the complexes used
#' @return output A matrix with 4 columns and 12 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 12 rows represent the 12 permutations of cross-validation using the samples.
#' @export
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
  print(paste0("result_rcc_twopept_",protrecscoreset,"_",size,"_",mode))
  return(RC_recov_protrec_CT1CT11)
}


#' @title This performs hyperparameter tuning of PROTREC based on scenario c using the Hela dataset
#' @name get_result_sce_c_hela
#' @param hela_ddaprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the Hela dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @param complex_size is the minimum size of the complexes used
#' @return output A matrix with 4 columns and 3 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 3 rows represent the validation results using the samples.
#' @export
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
  return(prores)
}

#' @title This performs hyperparameter tuning of PROTREC based on scenario c using the Siha dataset
#' @name get_result_sce_c_siha
#' @param siha_ddaprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the Siha dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @param complex_size is the minimum size of the complexes used
#' @return output A matrix with 4 columns and 3 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 3 rows represent the validation results using the samples.
#' @export
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
  return(prores)
}


#' @title This performs hyperparameter tuning of PROTREC based on scenario d using the Hela dataset
#' @name get_result_sce_d_hela
#' @param hela_ddaprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the Siha dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @param complex_size is the minimum size of the complexes used
#' @return output A matrix with 4 columns and 6 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 6 rows represent the cross-validation results using the samples.
#' @export
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
  return(Hela_ddadiacrossrecovprotrec1)
}

#' @title This performs hyperparameter tuning of PROTREC based on scenario d using the Siha dataset
#' @name get_result_sce_d_siha
#' @param siha_ddaprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the Siha dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @param complex_size is the minimum size of the complexes used
#' @return output A matrix with 4 columns and 9 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 9 rows represent the cross-validation results using the samples.
#' @export
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
  return(Siha_ddadiacrossrecovprotrec1)
}


#' @title This performs hyperparameter tuning of PROTREC based on scenario e using the RCC dataset
#' @name get_result_sce_e_rcc
#' @param rc_cprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the RCC dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @return output A matrix with 8 columns and 6 rows: each four columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The first four columns correspond to sample 1 (validated using sample 2), the second four columns
#' correspond to sample 2 (validated using sample 1). The 6 rows represent the 6 technical replicates
#' of each sample.
#' @export
get_result_sce_e_rcc <- function(rc_cprotrec,protrecscoreset,mode){
  RC_recov_protrec_CT1CT2 <- c()
  for (k in c(1:6))  #loop for CT1 -> CT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[k,], rownames(RC_C),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)

  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_C)[(RC_C[,k])!=0],RC_ccle,complex_vector,protrecscoreset)
  RC_recov_protrec_CT1CT2<- rbind(RC_recov_protrec_CT1CT2,val[1:4])
  print(k)
  }

  RC_recov_protrec_CT2CT1 <- c()
  for (l in c(7:12))  #loop for CT2 -> CT2_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-rc_cprotrec[l,], rownames(RC_C),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)

  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(RC_C)[(RC_C[,l])!=0],RC_ccle,complex_vector,protrecscoreset)
  RC_recov_protrec_CT2CT1<- rbind(RC_recov_protrec_CT2CT1,val[1:4])
  print(l)
  }


  RC_recov_cross_protrec<- data.frame(RC_recov_protrec_CT1CT2,RC_recov_protrec_CT2CT1)

  colnames(RC_recov_cross_protrec)<- c(rep("CT1->CT2",4),rep("CT2->CT1",4))

  return(RC_recov_cross_protrec)
}


#' @title This performs hyperparameter tuning of PROTREC based on scenario e using the LC dataset
#' @name get_result_sce_e_lc
#' @param lc_tprotrec A matrix of PROTREC (1-probabilities) assigned to complexes for the LC dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @return output A matrix with 8 columns and 19 rows: each four columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The first four columns correspond to sample 1 (validated using sample 2), the second four columns
#' correspond to sample 2 (validated using sample 1). The 19 rows represent the 19 technical replicates.
#' @export
get_result_sce_e_lc <- function(lc_tprotrec,protrecscoreset,mode){

  LC_recov_protrec_CT1CT2 <- c()
  for (k in c(1:19))  #loop for CT1 -> CT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-lc_tprotrec[k,], rownames(LC_T),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)

  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(LC_T)[(LC_T[,k])!=0],Liver_ccle,complex_vector,protrecscoreset)
  LC_recov_protrec_CT1CT2<- rbind(LC_recov_protrec_CT1CT2,val[1:4])
  print(k)
  }

  LC_recov_protrec_CT2CT1 <- c()
  for (l in c(20:38))  #loop for CT2 -> CT2_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-lc_tprotrec[l,], rownames(LC_T),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)

  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(LC_T)[(LC_T[,l])!=0],Liver_ccle,complex_vector,protrecscoreset)
  LC_recov_protrec_CT2CT1<- rbind(LC_recov_protrec_CT2CT1,val[1:4])
  print(l)
  }


  LC_recov_cross_protrec<- data.frame(LC_recov_protrec_CT1CT2,LC_recov_protrec_CT2CT1)

  colnames(LC_recov_cross_protrec)<- c(rep("CT1->CT2",4),rep("CT2->CT1",4))

  return(LC_recov_cross_protrec)
}


#' @title This performs hyperparameter tuning of PROTREC based on scenario e using the HCC dataset
#' @name get_result_sce_e_hcc
#' @param hcc_protrec A matrix of PROTREC (1-probabilities) assigned to complexes for the HCC dataset
#' @param protrecscoreset The threshold of PROTREC score, default 0.95
#' @param mode is the mode of PROTREC score selection, default max
#' @return output A matrix with 4 columns and 15 rows: the columns contain information about
#' the observed proportion of overlap (the PROTREC score), the significance of this overlap (p-value),
#' the total number of predicted missing proteins, and the number of validated proteins.
#' The 15 rows represent the 15 samples.
#' @export
get_result_sce_e_hcc <- function(hcc_protrec,protrecscoreset,mode){

  HCC_recov_protrec_CT1CT2 <- c()
  for (k in c(1:15))  #loop for CT1 -> CT1_pep
  {PROTREC_prot <- data.frame(PROTREC_protprob(complex_vector, 1-hcc_protrec[k,], rownames(HCC),0.01,mode))
  rownames(PROTREC_prot)<- PROTREC_prot[,1]
  PROTREC_prot[,2]<-as.numeric(as.character(PROTREC_prot[,2]))
  PROTREC_prot<- data.frame(PROTREC_prot)

  val<- pairwise_recovery_protrec(PROTREC_prot,rownames(HCC)[(HCC[,k])!=0],Liver_ccle,complex_vector,protrecscoreset)
  HCC_recov_protrec_CT1CT2<- rbind(HCC_recov_protrec_CT1CT2,val[1:4])
  print(k)
  }

  HCC_recov_cross_protrec<- data.frame(HCC_recov_protrec_CT1CT2)

  colnames(HCC_recov_cross_protrec)<- c(rep("CCLE veri",4))

  return(HCC_recov_cross_protrec)
}

###FCS filtering for PROTREC########
#' @title This produces the distribution of PROTREC scores for the predicted, validated and original proteins
#' @name countgraph
#' @param x the probability of proteins observed in the sample
#' @param y the names of predicted proteins in the sample
#' @param z the names in the check list
#' @param cplx the complex information
#' @param g the set of proteins in the original list
#' @return output four columns contain information about
#' the PROTREC score range, the total number of predicted proteins, the number of validated proteins, and
#' the total number of proteins in the original list.
#' @export
countgraph <- function(x,y,z, cplx, g)
{
  output <- c()
  names(x) <- y
  levels <- seq(0, 1, by=0.1) #the distribution of probability levels
  for (i in 1:length(levels))
  {
    total <- length(names(x)[x >= levels[i] & x < (levels[i] + 0.1)]) # total number of predicted prot within this range
    original <- length(intersect(names(x)[x >= levels[i] & x < (levels[i] + 0.1)], g)) # how many were in original screen
    validated <- length(intersect(names(x)[x >= levels[i] & x < (levels[i] + 0.1)], setdiff(z, g))) # validated prot

    total <- total - original - validated # unvalidated

    output <- rbind(output, cbind(levels[i], total, validated, original))
  }
  return(output)
}

#' @title This produces the cumulative distribution of PROTREC scores for the predicted, validated and original proteins
#' @name precision_recall_2
#' @param x the probability of proteins observed in the sample
#' @param y the names of predicted proteins in the sample
#' @param z the names in the check list
#' @param cplx the complex information
#' @param g the set of proteins in the original list
#' @return output four columns contain information about
#' the PROTREC score range, the total number of predicted proteins, the number of validated proteins, and
#' the total number of proteins in the original list.
#' @export
precision_recall_2 <- function(x,y,z, cplx, g)
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

#' @title Takes the FCS significant complexes and ranks them by the tibtech approach
#' @name get_sig_complex
#' @param replicate the FCS significant complexes derived using a given replicate
#' @return output the filtered and ranked FCS significant complexes.
#' @export
get_sig_complex <- function(replicate){
  FCS_sig_complex_rep <- as.vector(replicate[which(replicate <= 0.05)])
  return(FCS_sig_complex_rep)
}

#' @title obtain the FCS significant complexes that overlap with the complexes present in fcs_normal_mat
#' @name get_mat_n
#' @param replicate the FCS significant complexes derived using a given replicate
#' @return output the those proteins in both FCS results and those that make up complexes in fcs_normal_mat
#' @export
get_mat_n <-function(replicate){
  fcs_normal_mat_rep <- fcs_normal_mat[which(rownames(fcs_normal_mat) %in% unique(unlist(complex_vector[names(replicate)]))),]
  return(fcs_normal_mat_rep)
}

#' @title obtain the FCS significant complexes that overlap with the complexes present in fcs_cancer_mat
#' @name get_mat_c
#' @param replicate the FCS significant complexes derived using a given replicate
#' @return output the those proteins in both FCS results and those that make up complexes in fcs_cancer_mat
#' @export
get_mat_c <-function(replicate){
  fcs_cancer_mat_rep <- fcs_cancer_mat[which(rownames(fcs_cancer_mat) %in% unique(unlist(complex_vector[names(replicate)]))),]
  return(fcs_cancer_mat_rep)
}

#' @title obtain the FCS significant complexes that overlap with the complexes present in normal_mat
#' @name get_protprob_normal_mat
#' @param replicate the FCS significant complexes derived using a given replicate
#' @return output the those proteins in both FCS results and those that make up complexes in normal_mat
#' @export
get_protprob_normal_mat <- function(replicate){
  protprob_normal_mat_rep <- normal_mat[rownames(normal_mat)%in%unique(unlist(complex_vector[names(replicate)])),]
  return(protprob_normal_mat_rep)
}

#' @title obtain the FCS significant complexes that overlap with the complexes present in cancer_mat
#' @name get_protprob_cancer_mat
#' @param replicate the FCS significant complexes derived using a given replicate
#' @return output the those proteins in both FCS results and those that make up complexes in cancer_mat
#' @export
get_protprob_cancer_mat <- function(replicate){
  protprob_cancer_mat_rep <- cancer_mat[rownames(cancer_mat)%in%unique(unlist(complex_vector[names(replicate)])),]
  return(protprob_cancer_mat_rep)
}

#' @title modify the distribution
#' @name change_name
#' @param replicate the distribution of PROTREC scores for the predicted, validated and original proteins
#' @return output the modified variable
#' @export
change_name <- function(replicate){
  rownames(replicate) <- replicate[,1]
  replicate <- replicate[,-1]
}

#' @title sort the cumulative distribution
#' @name order_stuff
#' @param replicate the cumulative distribution of PROTREC scores for the predicted, validated and original proteins
#' @return output the sorted variable
#' @export
order_stuff <- function(replicate){
  return(replicate[order(replicate[,1], decreasing = T),])
}

#' @title This performs fcs-filtering of complexes before applying the PROTREC method using the RCN dataset
#' @name get_average_values_n
#' @param FCS_sig_complex_n A vector containing proteins associated with significant complexes in FCS
#' @return output A bar plot and a line plot of the comparison between FCS-filtered PROTREC performance
#' and non-FCS-filtered PROTREC performance.
#' @export
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
  pdf("./output/section4/fcs_sig_complexes_protprob_methods_barplot_rcn.pdf",width=8,height=6)
  par(mfrow=c(1,2))
  barplot(t(ave_fcs_n), col=c("paleturquoise", "royalblue2", "dodgerblue4"), main="FCS (RCN)", ylab="Count", xlab="Prob(Protein p is present but unreported)", cex.lab=1.2, font=2)
  barplot(t(ave_protrec_n), col=c("paleturquoise", "royalblue2", "dodgerblue4"), main="PROTREC (RCN)", ylab="Count", xlab="Prob(Protein p is present but unreported)", cex.lab=1.2,font=2)

  # legend
  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(-0.45,0.9,legend = c("Original", "Validated", "Unvalidated"), col= c("dodgerblue4", "royalblue2","paleturquoise"), lwd = 10, xpd = TRUE, horiz = TRUE, cex = 0.8, seg.len=1, bty = 'n')
  dev.off()

  # plot line plot
  pdf("./output/section4/protprob_filteredvsunfiltered_rcn.pdf",width=9,height=8)
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

#' @title This performs fcs-filtering of complexes before applying the PROTREC method using the RCC dataset
#' @name get_average_values_c
#' @param FCS_sig_complex_c A vector containing proteins associated with significant complexes in FCS
#' @return output A bar plot and a line plot of the comparison between FCS-filtered PROTREC performance
#' and non-FCS-filtered PROTREC performance.
#' @export
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
  pdf("./output/section4/fcs_sig_complexes_protprob_methods_barplot_rcc.pdf",width=8,height=6)
  par(mfrow=c(1,2))
  barplot(t(ave_fcs_c), col=c("paleturquoise", "royalblue2", "dodgerblue4"), main="FCS (RCC)", ylab="Count", xlab="Prob(Protein p is present but unreported)", cex.lab=1.2, font=2)
  barplot(t(ave_protrec_c), col=c("paleturquoise", "royalblue2", "dodgerblue4"), main="PROTREC (RCC)", ylab="Count", xlab="Prob(Protein p is present but unreported)", cex.lab=1.2,font=2)

  # legend
  par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(-0.45,0.9,legend = c("Original", "Validated", "Unvalidated"), col= c("dodgerblue4", "royalblue2","paleturquoise"), lwd = 10, xpd = TRUE, horiz = TRUE, cex = 0.8, seg.len=1, bty = 'n')
  dev.off()

  # plot line plot
  pdf("./output/section4/protprob_filteredvsunfiltered_rcc.pdf",width=9,height=8)
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
