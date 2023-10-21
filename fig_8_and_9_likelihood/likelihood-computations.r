
## Likelihood simulations for standard coalescent with three lineages (Section 5.3)
## Two scenarios:

## packages used: phyclust, xtable

library(phyclust)  # this is to make ms simulator by R. Hudson available within R




### functions needed to simulate data and compute likelihoods


sim.mutations.per.lineage<- function(nsim,nsam=3,max_sam=12,theta=1.0){
    # Coalescent 
    # six branches: a, b, c, ab, ac, bc (only one doubleton branch is nonzero)
    # works with nsim >= 1, nsam>=3
    # * if nsam = 3, then nsim independent loci  are simulated, each with sample site nsam = 3 
    # * if nsam >3, then mutational patterns computed for subsets of size 3 taken from the nsam genes
    # either for each subset, or for max_sam randomly chosen subsets, whatever is smaller
  
  
    # nsam>=3 genes needed per locus:
    if(nsam < 3)error("nsam needs to be at least 3")
  
    
 
    if(nsam==3)   {
      branch.mut=matrix(NA,nrow=nsim,ncol=6)
      for(i in 1:nsim) {
    
      ret_ms <- ms(nsam, opts=eval(paste0("-t ",theta)))  
      if(substring(ret_ms[3],11,11) =="0")branch.mut[i,] <- rep(0,6)  else 
      {
       haplotypes<-ret_ms[(length(ret_ms)-nsam+1):length(ret_ms)]
      
       a <- as.vector(as.numeric(unlist(strsplit(haplotypes[1], ""))))
       b <- as.vector(as.numeric(unlist(strsplit(haplotypes[2], ""))))
       c <- as.vector(as.numeric(unlist(strsplit(haplotypes[3], ""))))
       ab<- a+b>1
       ac<- a+c>1
       bc<- b+c>1
       seg.sites <-cbind(a,b,c,ab,ac,bc)               
       sums <- apply(seg.sites,2,sum)
        branch.mutations<-c(sums[1] -sum(sums[4:5]),
                        sums[2] -sum(sums[c(4,6)]),
                        sums[3] -sum(sums[c(5,6)]),sums[4:6])
        branch.mut[i,] <-branch.mutations}
    
    } } else
      {
    
    nr_subsets  <- choose(nsam,3)
    if(nr_subsets<=max_sam) {

      all_subsets <- combn(1:nsam, 3, simplify = FALSE) 
      }
    else
      {
      #message("as there are more than max_sam subsets we randomly sample subsets")
      nr_subsets  <- max_sam
      all_subsets <- lapply(rep(nsam,max_sam),sample,size=3) 
      }
   
    
    
    branch.mut=matrix(NA,nrow=nsim*nr_subsets,ncol=6)
    
    for(i in 1:nsim) {
      ret_ms    <- ms(nsam, opts=eval(paste0("-t ",theta))) 
      if(substring(ret_ms[3],11,11) =="0")branch.mut[i,] <- rep(0,6)  else 
      {
        haplotypes<-ret_ms[(length(ret_ms)-nsam+1):length(ret_ms)]
       
        for(j in 1:nr_subsets) {
  
          subset_haplotypes <-haplotypes[all_subsets[[j]]]
          a <- as.vector(as.numeric(unlist(strsplit(subset_haplotypes[1], ""))))
          b <- as.vector(as.numeric(unlist(strsplit(subset_haplotypes[2], ""))))
          c <- as.vector(as.numeric(unlist(strsplit(subset_haplotypes[3], ""))))
        
          not_segregating<- apply(cbind(a,b,c),1,sum)==3
          a[not_segregating]= b[not_segregating]= c[not_segregating]=0
          ab<- a+b>1
          ac<- a+c>1
          bc<- b+c>1
          seg.sites <-cbind(a,b,c,ab,ac,bc)               
          sums <- apply(seg.sites,2,sum)
          branch.mutations<-c(sums[1] -sum(sums[4:5]),
                          sums[2] -sum(sums[c(4,6)]),
                          sums[3] -sum(sums[c(5,6)]),sums[4:6])
      
          if(is.na(branch.mutations[1]))  branch.mut[(i-1)*nr_subsets,] <- rep(0,6)  else {
            branch.mut[(i-1)*nr_subsets + j,] <-branch.mutations}
          }}
    
        }
    }
    
    

   branch.mut
  }



loglik <- function(mutational_patterns,lambda_from=0,lambda_to=1,lambda_by=0.005){
  
  # compute likelihoods
  # for each mutational pattern (row of matrix mutational_patterns) and each lambda in a grid of parameter values
  
  lambda_grid <- seq(lambda_from,lambda_to,by=lambda_by)
  nsamp= dim(mutational_patterns)[1]
  logliks <- matrix(NA,nrow=nsamp,ncol=length(lambda_grid))
  for(j in 1:length(lambda_grid)){
    for(i in 1:nsamp){
      logliks[i,j]=log(likelihood2(sort_mut_counts(mutational_patterns[i,]),lambda_grid[j]))
    }}
    logliks
    
}


 


  
 likelihood1<- function(m1,m2,m3,m4,lambda){
   # likelihood computation for one mutational pattern
   #m1,m2,m3 are the frequencies for the branches a, b, c; m4 is doubleton branch
   # c is long branch, or if m4 is 0 all branches are taken as possible 
 
   j<-0:m3
   
   lambda^(m1+m2+m3+m4)/(factorial(m1)*factorial(m2)*factorial(m3)*
                           factorial(m4))*sum(choose(m3,j)*factorial(m4+j)*
                                                factorial(m1+m2+m3-j)/(1+2*lambda)^(m4+j+1)/(3+3*lambda)^(m1+m2+m3-j+1)) 
   
   
 }
 
 
 likelihood2<- function(ms,lambda){
   # likelihood computation for one mutational pattern
   #m1,m2,m3 are the frequencies for the branches a, b, c; m4 is doubleton branch
   # c is either long branch, or if m4 is 0 all branches are taken as possible 
   # long branches

   m1=ms[1];m2=ms[2];m3=ms[3];m4=ms[4];
   j<-0:m3
   if (m4>0)  likelihood1(m1,m2,m3,m4,lambda)
   else
   {likelihood1(m1,m2,m3,m4,lambda)+likelihood1(m1,m3,m2,m4,lambda)+likelihood1(m3,m2,m1,m4,lambda)}
 }
 
 
 
 
 sort_mut_counts <- function(mut.pat){
   # condense mutational counts provided by function sim.mutations.per.lineage
   # put long branch at position 3, short branches at positions 1:2,
   # and doubleton at position 4  
   doubleton <- mut.pat[4:6]>0
   if(is.na(sum(doubleton)))print(mut.pat)
   if(sum(doubleton)==0) ms= c(mut.pat[1:3],0)  else
   {
     
     ms=c(mut.pat[c(1:3)[-c(3:1)[doubleton] ]],mut.pat[c(3:1)[doubleton]],
          sum(mut.pat[4:6]))
   }
   ms
 }
 
 
 simulate.MLE <- function(nsimul=5,theta=2.1,nsim=20,nsam=3,maxsam=100,theta_by=0.005){
   # simulations to compute MLE and Watterson's estimator
   # nsimul simulation runs are carried out
   # theta = 2*lambda is the true standard scaled mutation parameter
   # each simulation run involves nsim independent loci, each with a sample of size nsam genes,
   # maxsam is the maximum sample size for estimating the composite likelihood when nsam > 3
   # theta_by is the grid width for a grid between theta-2 and theta+2 used to search for the MLEE
   # grid parameter value giving the largest likelihood is taken as MLE
   MLEs =rep(NA,nsimul)
   Watterson_theta= rep(NA,nsimul)
   lambda_from=theta/2-1;lambda_to=theta/2+2;lambda_by=theta_by/2
   lambda_grid=seq(from=lambda_from,to=lambda_to,by=theta_by/2)  
   logliks = matrix(NA,nrow=nsimul,ncol=length(lambda_grid))
   for(k in 1:nsimul){
     mutational_patterns<-sim.mutations.per.lineage(nsim=nsim,nsam=nsam,max_sam=maxsam,theta = theta)  
     logliks[k,] = apply(loglik(mutational_patterns,
                        lambda_from,lambda_to,lambda_by),2,sum)
     #MLE
     MLEs[k]=2*lambda_grid[which.max(logliks[k,])]
     # Watterson's theta
     Watterson_theta[k]=mean(apply(mutational_patterns,1,sum))*2/3
   }
   
   list(estimates=cbind(MLE=MLEs,Watterson=Watterson_theta), thetas=lambda_grid*2,
        likelihoods=logliks)
   
   
 }
 
 



# plot likelihood for Figure 7


ms=c(0,1,1,3)

theta_grid=seq(2,4,by=0.1)
log_lik=rep(NA,length(theta_grid))
for(i in 1:length(theta_grid)) log_lik[i]=log(likelihood2(ms,theta_grid[i]/2))

pdf("likelihood-plotB.pdf")
plot(theta_grid,log_lik,type="l",lwd=2,xlab=expression(theta),xlim=c(2,4),ylim=-c(7.6,7.2),ylab="log-likelihood")
watterson = sum(ms)*2/3
watterson
lines(c(watterson,watterson),c(-7.6,log_lik[findInterval(watterson,theta_grid)]),lty=2,lwd=2)
MLE = theta_grid[which.max(log_lik)]
MLE
lines(c(MLE,MLE),c(-7.6,log_lik[which.max(log_lik)]),lty=3,lwd=2)
text(x=c(MLE-0.193, watterson+0.01), y=c(-7.55, -7.55), pos=4, labels=c('MLE', 'Watterson'))
dev.off()

# For Figure 8 and Table 2


# 1 locus 20 lineages

res<-simulate.MLE(nsimul=100,theta=3,nsim=1,nsam=20,maxsam=1000,
                  theta_by=0.02)                                         # computation takes long!!

# 20 independent loci

res2<-simulate.MLE(nsimul=100,theta=3,nsim=20,nsam=3,maxsam=1000,
                   theta_by=0.02)

# n=3, 1 locus

res1<-simulate.MLE(nsimul=500,theta=3,nsim=1,nsam=3,maxsam=1000,
                   theta_by=0.02)




# result summaries

c(apply(res$estimates,2,mean), apply(res2$estimates,2,mean), apply(res1$estimates,2,mean))
c(apply(res$estimates,2,sd), apply(res2$estimates,2,sd),  apply(res1$estimates,2,sd))

MSE   = apply((res$estimates-3)^2,2,mean)  # 1 locus 20 lineages
MSE2  = apply((res2$estimates-3)^2,2,mean) # 20 independent loci
MSE1  = apply((res1$estimates-3)^2,2,mean) # n=3, 1 locus


# Table 2

library(xtable)

print(xtable(sqrt(rbind(MSE,MSE2,MSE1)),digits=3))


# plot likelihood functions

# Figure 8a

pdf("20-loci-likelihood-plot.pdf")
plot(res2$thetas,res2$likelihoods[1,],type="l",lwd=2,xlab=expression(theta),ylab="log-likelihood",xlim=c(1.5,5),ylim=c(-122,-115))
watterson = res2$estimates[1,2]

lines(c(watterson,watterson),c(-129.5,res2$likelihoods[1,findInterval(watterson,res2$thetas)]),lty=2,lwd=2)
MLE = res2$thetas[which.max(res2$likelihoods[1,])]
MLE
lines(c(MLE,MLE),c(-129.5,res2$likelihoods[1,which.max(res2$likelihoods[1,])]),lty=3,lwd=2)
text(x=c(MLE, watterson-0.6), y=c(-120, -120), pos=4, labels=c('MLE', 'Watterson'))
dev.off()

# Figure 8b

pdf("comp-likelihood-plot.pdf")

plot(res$thetas,res$likelihoods[1,],type="l",lwd=2,xlab=expression(theta),ylab="composite log-likelihood")
watterson = res$estimates[1,2]

lines(c(watterson,watterson),c(-8500,res$likelihoods[1,findInterval(watterson,res$thetas)]),lty=2,lwd=2)
MLE = res$thetas[which.max(res$likelihoods[1,])]
MLE
lines(c(MLE,MLE),c(-8500,res$likelihoods[1,which.max(res$likelihoods[1,])]),lty=3,lwd=2)
text(x=c(MLE-0.55, watterson+0.02), y=c(-8200, -8200), pos=4, labels=c('MLE', 'Watterson'))
dev.off()




