# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# hello <- function() {
#  print("Hello, world!")
# }


#main function giving estimates for missing mutations in each copy number state and for each peak,
#for missing mutations cumulated over all peaks, the signatures and estimates for the numbers of observations
#in different peaks depending on their mut ational context

missing_mutations <- function(vcf_path, copy_number_path, p, c=min(SNV$allele_counts)) {
    vcf <- readVcfAsVRanges(vcf_path, "hg19")
    copy_number = read.csv(copy_number_path)
    select = select_func(vcf,copy_number)
    relevant = relevant_copy_number_states(select)
    SNV = SNV_info(select,relevant)
    ML=max_lik(p,SNV,relevant)
    numbers = find_missing(SNV,relevant,ML,p,c)
    sigs = signatures(select)
    emp = empirical_distr(relevant,SNV,p)
    peak_0_given_bin = peak_0_given_bin(relevant,SNV,ML,numbers,emp,p)
    peak_0_given_af = peak_0_given_af(SNV,peak_0_given_bin)
    prob_of_peaks = prob_of_peaks(relevant,SNV,ML,peak_0_given_af,p)
    mutational_contexts_in_peaks = mutational_contexts_in_peaks(SNV, prob_of_peaks)
    return(list(missing=numbers,in_each_peak=round(apply(numbers[,(3:(ncol(numbers)-1))],2,sum)),signatures=sigs,
                mutational_contexts_in_peaks=mutational_contexts_in_peaks))
}


get_peak_centres <- function(Total_CN,Minor_CN,p) {
    return_vector = c()
    for (i in 1:(Total_CN-Minor_CN)){
        return_vector[i]=(i*p)/(p*Total_CN+2*(1-p))
    }
    return(return_vector)
}


select_func <- function(vcf,copy_number) {
    vcf <- keepSeqlevels(vcf, c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                                "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
    vcf <- vcf[(length(vcf)/2+1):length(vcf)]
    vcf <- vcf[apply(vcf@softFilterMatrix,1,prod)==1]
    x <- GRanges(copy_number)
    seqlevels(x) <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                      "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
    x <- keepSeqlevels(x, c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                            "chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
    mc <- mutationContext(vcf,BSgenome.Hsapiens.UCSC.hg19)
    SMT <- paste(subseq(elementMetadata(mc)$context,1,1), rep("[",length(mc)),
                 subseq(elementMetadata(mc)$alteration,1,1), rep(">",length(mc)),
                 subseq(elementMetadata(mc)$alteration,2,2), rep("]",length(mc)),
                 subseq(elementMetadata(mc)$context,3,3), sep="")
    return(list(vcf,x,SMT))
}

#we exclude all copy number states with less than 20 observations
relevant_copy_number_states <- function(select) {
    f <- function(row) {
        if(length(subsetByOverlaps(select[[1]],select[[2]][which(select[[2]]@elementMetadata$Total_CN==row[1]
                                                                 & select[[2]]@elementMetadata$Minor_CN==row[2])]))>19)
            tab=c(row[1],row[2])
    }
    rel=data.frame(matrix(unlist(apply(unique(select[[2]]@elementMetadata),1,f)),ncol=2,byrow=TRUE))
    colnames(rel)=c("Total_CN","Minor_CN")
    return(rel[order(rel$Total_CN,rel$Minor_CN),])
}

#gives data frame with SNVs, their respective allele counts, read counts, allele fractions, copy number states
#and mutational contexts
SNV_info <- function(select,relevant) {
    VAF = select[[1]]$VariantAlleleFrequency
    VAC = select[[1]]$VariantAlleleCount
    RC  = select[[1]]$ReadCount
    SNVs_in_copy_number_states <- function(row) {
        regions_with_desired_cns <- select[[2]][which(select[[2]]@elementMetadata$Total_CN==row[1] & select[[2]]@elementMetadata$Minor_CN==row[2])]
        indices <- as.data.frame(findOverlaps(select[[1]],regions_with_desired_cns))$queryHits   #gives indices of SNVs in desired copy number state
        h = cbind(as.data.frame(cbind(VAF,VAC,RC)[indices,]), mut = select[[3]][indices])
        delete_outliers = mean(h$RC)-3*sd(h$RC) <= h$RC & h$RC <= mean(h$RC)+3*sd(h$RC)          #delete outliers in number of reads
        number_of_reads = h$RC[delete_outliers]
        allele_fractions = h$VAF[delete_outliers]
        allele_counts = h$VAC[delete_outliers]
        mutational_context = h$mut[delete_outliers]
        h_result = cbind(as.data.frame(cbind(Total_CN=row[1],Minor_CN=row[2],allele_fractions,allele_counts,
                                             number_of_reads)),mutational_context)
    }
    return(as.data.frame(do.call("rbind",apply(relevant,1,SNVs_in_copy_number_states))))
}


max_lik <- function(p,SNV,relevant) {
    f <- function(row) {
        peak_centres = get_peak_centres(row[1],row[2],p)
        #the starting values for the maximum likelihood estimation shall be chosen in a senseful way:
        find_closest_peaks <- function(m) {max(which(abs(peak_centres-m)==min(abs(peak_centres-m))))}
        #closest is a vector showing, for each allele fraction in the copy number state given, which peak it is closest to
        closest=mapply(find_closest_peaks,SNV[which(SNV$Total_CN==row[1] & SNV$Minor_CN==row[2]),]$allele_fractions)
        #if there is more than one peak, find a good starting value, otherwise the weight shall be 1
        #we use a log transformation to make sure the resulting weights are within [0,1] and sum to 1
        if (row[1]-row[2]>1) {
            a=rep(0,length(peak_centres))
            for (i in 1:length(closest)) a[closest[i]]=a[closest[i]]+1
            b=(a/sum(a))+0.0001    #so we do not take log(0)
            c=(b/sum(b))
            start=log(c[1:(length(c)-1)]/c[length(c)])} else start=numeric(0)
        #only consider observations on the left hand side of peak 1 as we might have an extra "peak 0" violating the mixture of binomials
        SNV_cut=SNV[SNV$Total_CN==row[1] & SNV$Minor_CN==row[2] & SNV$allele_fractions>=peak_centres[1],]
        likelihood <- function(weights) {
            weights_mod=c(exp(weights)/(1+sum(exp(weights))),1/(1+sum(exp(weights))))
            g1=0
            for (j in 1:length(weights_mod))
                g1=g1+weights_mod[j]*dbinom(SNV_cut$allele_counts,SNV_cut$number_of_reads,peak_centres[j])
            g2=1
            for (k in 1:length(weights_mod))
                g2=g2-(weights_mod[k]*pbinom(round(peak_centres[1]*SNV_cut$number_of_reads),
                                             SNV_cut$number_of_reads,peak_centres[k]))
            g3=g1/g2
            return(-sum(log(g3)))
        }
        fit = optim(start,likelihood,method="BFGS")
        numbers = c(peak=c(exp(fit$par)/(1+sum(exp(fit$par))), 1/(1+sum(exp(fit$par))),
                           rep(0,max(relevant$Total_CN-relevant$Minor_CN)-length(fit$par)-1)),row[1],row[2])}
    numbers = transpose(apply(relevant,1,f))
    return(as.data.frame(numbers))
}




#find_missing where we estimate number of missing mutations in each peak separately
find_missing <- function(SNV,relevant,ML,p,c){
    h <- function(row) {
        peak_centres = get_peak_centres(row[1],row[2],p)
        SNV_cut = SNV[SNV$Total_CN==row[1] & SNV$Minor_CN==row[2] & SNV$allele_fractions>=peak_centres[1],]
        weights = ML[ML$Total_CN==row[1] & ML$Minor_CN==row[2],]
        prob_1 = 0
        #P(allele_fractions>peak_centre_1)
        for (i in 1:(row[1]-row[2]))
            prob_1 = prob_1+sum((as.numeric(weights[i])*(1-pbinom(round(peak_centres[1]*SNV_cut$number_of_reads),
                                                                  SNV_cut$number_of_reads,peak_centres[i])))*
                                    (1/length(SNV_cut$number_of_reads)))
        total_number_estimated = round(nrow(SNV_cut)/prob_1)     #estimated total number of observations we would have with the mixture of binomials
        prob_2 =rep(0,max(relevant$Total_CN-relevant$Minor_CN))
        #vector containing P(allele_counts<min(SNV$allele_counts)) for each peak
        for (i in 1:(row[1]-row[2]))
            prob_2[i] = sum(pbinom(c,SNV_cut$number_of_reads,peak_centres[i])*
                                (1/length(SNV_cut$number_of_reads)))
        numbers_in_peaks = rep(0,max(relevant$Total_CN-relevant$Minor_CN))
        #vector containing estimated total number of observations we would have in each peak with mixture of binomials; without missing mutations
        for (i in 1:(row[1]-row[2]))
            numbers_in_peaks[i] = total_number_estimated*as.numeric(weights[i])
        missing_numbers = c(row[1],row[2],peak=numbers_in_peaks*prob_2,total_number_estimated=total_number_estimated)
        #vector containing numbers of missing mutations in each peak and total number estimated to be actually there
    }
    return(as.data.frame(transpose(apply(relevant,1,h))))
}


signatures <- function(select) {
    SMT <- factor(select[[3]],levels=colnames(signatures.cosmic))
    #contains mutational contexts for all SNVs in vcf
    SMTdf <- as.data.frame(t(matrix(table(SMT))))
    colnames(SMTdf) <- colnames(signatures.cosmic)
    rownames(SMTdf) <- 1
    sigs <- whichSignatures(tumor.ref = SMTdf, signatures.ref =
                                signatures.cosmic,
                            sample.id = 1, contexts.needed = TRUE)
    return(sigs$weights)
}


empirical_distr <- function(relevant,SNV,p) {
    breaks = seq(0,1,by=1/100)
    g <- function(row) {
        peak_centres = get_peak_centres(row[1],row[2],p)
        SNV_need = SNV[SNV$Total_CN==row[1] & SNV$Minor_CN==row[2],]
        m = rep(0,100)
        for (j in 1:100)
            m[j] = length(which(SNV_need$allele_fractions>=breaks[j] &
                                    SNV_need$allele_fractions<breaks[j+1]))
        res = c(prob=m/sum(m),row[1],row[2])
    }
    return(as.data.frame(transpose(apply(relevant,1,g))))   #data frame with empirical probabilities, based on a histogram with 100 bins,
    #for each copy number state
}


#The output of this function is a matrix that, depending on the copy number state, gives P("peak_0"|bin), where we consider 100 bins
peak_0_given_bin <- function(relevant,SNV,ML,find,emp,p) {
    breaks=seq(0,1,by=1/100)
    g <- function(row) {
        peak_centres = get_peak_centres(row[1],row[2],p)
        weights = ML[ML$Total_CN==row[1] & ML$Minor_CN==row[2],]
        SNV_need = SNV[SNV$Total_CN==row[1] & SNV$Minor_CN==row[2],]
        SNV_cut = SNV[SNV$Total_CN==row[1] & SNV$Minor_CN==row[2] & SNV$allele_fractions>=peak_centres[1],]
        emp_need = emp[emp$Total_CN==row[1] & emp$Minor_CN==row[2],]
        b = rep(0,100)
        for (i in 1:100) {
            for (j in 1:(row[1]-row[2])) {
                b[i] = b[i]+sum(as.numeric(weights[j])*pbinom(round(breaks[i+1]*SNV_cut$number_of_reads),SNV_cut$number_of_reads,peak_centres[j])-
                                    as.numeric(weights[j])*pbinom(round(breaks[i]*SNV_cut$number_of_reads),SNV_cut$number_of_reads,peak_centres[j]))
            }
        }
        res = b/sum(b)  #probability distribution based on mixture of binomials, cumulated for each of the 100 bins (red line)
        total_number_estimated = find[find$Total_CN==row[1] & find$Minor_CN==row[2],]$total_number_estimated #estimated total number of SNVs with
        #mixture of binomials
        n = round(res*total_number_estimated)      #estimated number of observations in each bin under the red line (mixture of binomials)
        m = round(as.numeric(emp_need[1:100]*nrow(SNV_need)))   #estimated number of observations in each bin under the blue line (empirical)
        m_cut = c(m[1:floor(peak_centres[1]*100)],rep(0,100-length(m[1:floor(peak_centres[1]*100)])))
        n_cut = c(n[1:floor(peak_centres[1]*100)],rep(0,100-length(n[1:floor(peak_centres[1]*100)])))
        prob=rep(0,100)
        prob[m_cut==0 | n_cut>m_cut] = 0
        prob[m_cut!=0 & n_cut<m_cut] = (m_cut[m_cut!=0 & n_cut<m_cut]-n_cut[m_cut!=0 & n_cut<m_cut])/m_cut[m_cut!=0 & n_cut<m_cut]
        #P("peak_0"|bin) for each bin; depending on copy number state
        erg=c(bin=prob,row[1],row[2])
    }
    return(as.data.frame(transpose(apply(relevant,1,g))))
}


#gives P("peak_0"|allele_fraction) for all allele fractions
peak_0_given_af <- function(SNV,peak_0_given_bin) {
    breaks=seq(0,1,by=1/100)
    f <- function(row) {
        prob_peak0_need=peak_0_given_bin[peak_0_given_bin$Total_CN==row[1] & peak_0_given_bin$Minor_CN==row[2],]
        a=as.numeric(prob_peak0_need[which.max(breaks[breaks<=row[3]])])
    }
    return(apply(SNV[,1:5],1,f))
}


#gives matrix with rows P("peak_0"|allele_fraction),...,P(peak_(Total_CN-Minor_CN)|allele_fraction) for all allele fractions
prob_of_peaks <- function(relevant,SNV,ML,peak_0_given_af,p) {
    t = as.data.frame(table(SNV$number_of_reads))   #empirical probability distribution of number of reads
    colnames(t) <- c("read_count","frequency")
    breaks = seq(0,1,by=1/100)
    f <- function(row) {
        peak_centres = get_peak_centres(row[1],row[2],p)
        weights = ML[ML$Total_CN==row[1] & ML$Minor_CN==row[2],]
        denominator = 0
        for (i in 1:length(peak_centres))
            denominator = denominator+as.numeric(weights[i])*dbinom(row[4],row[5],peak_centres[i])*(t[which(t$read_count==row[5]),]$freq)/nrow(SNV)
        vector = c()
        for (i in 1:length(peak_centres))
            vector = c(vector,(as.numeric(weights[i])*dbinom(row[4],row[5],peak_centres[i])*((t[which(t$read_count==row[5]),]$freq)/nrow(SNV))/
                                   denominator)*as.numeric(1-row[6]))
        vector = c(row[6],vector,rep(0,max(relevant$Total_CN-relevant$Minor_CN)-length(vector)))
    }
    result = data.frame(matrix(transpose(unlist(apply(cbind(SNV[,1:5],peak_0_given_af),1,f))),
                               ncol=max(relevant$Total_CN-relevant$Minor_CN)+1))
    colnames(result) <- c(paste("peak",0:max(relevant$Total_CN-relevant$Minor_CN, sep="")))
    return(result)
}



mutational_contexts_in_peaks <- function(SNV, prob_of_peaks) {
    z = cbind(SNV,prob_of_peaks) #data frame with copy number state, allele fractions, allele counts,
    #number of reads, mutational context and probabilities
    a  = matrix(0,nrow=96,ncol=ncol(z)-6)
    for (i in 1:(ncol(z)-6)) {
        for (j in 1:96) {
            a[j,i] = sum(z[z$mutational_context==levels(z$mutational_context)[j],][,i+6])
        }
    }
    u = a/apply(a,1,sum)
    for (j in 1:96)
        u[j,] = round(u[j,]*nrow(z[z$mutational_context==levels(z$mutational_context)[j],]))
    return(u) #96x(Total_CN-Minor_CN+1) matrix containing estimated numbers of mutational contexts within
    #each peak
}



