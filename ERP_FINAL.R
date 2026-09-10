###
# File description: SAN-10 ERP
#Question - Do multidimensional and non-linear models provide an improved representation of political ideology in roll-call data?
# Last modified: 01/09/2026
###

#Clear directory
rm( list = ls() )

#Add packages
library(rstudioapi)
library(ggplot2)
library(readr)
library(dplyr)
library(plotly)
library(rjags)
library(coda)

#Set working directory
setwd(dirname(getActiveDocumentContext()$path))

#Load data
vm<-read.table("votematrix-2024.dat", header=TRUE, sep="\t", quote="",
               fill=TRUE, check.names=FALSE)
write.csv(vm, file="votematrix-2024.csv", row.names=FALSE)
View(vm)

nrow(vm)
ncol(vm)

#Format matrix
vm<-subset(vm, select=-c(rowid))
vm<-vm[, names(vm) != ""]
vm<-vm%>%select(voteno, everything())
vm<-vm[order(vm$voteno), ]
vm[, -1][vm[, -1]==-9]<-NA
vm[, -1][vm[, -1]==2]<-1
vm[, -1][vm[, -1]==4]<-0
vm[, -1][vm[, -1]==5]<-0
vm[, -1][vm[, -1]==3]<-NA


#MP matrix
mp_info<-read.table("votematrix-2024.txt", header=FALSE, sep="\t", skip=20,
                    quote="", fill=TRUE, check.names=FALSE)
colnames(mp_info)<-c("mpid", "firstname", "surname", "party", "url")
mp_info<-mp_info[, c("mpid", "firstname", "surname", "party")]
View(mp_info)

#Merge and transpose matrices
mp_info$mpid_key<-paste0("mpid", mp_info$mpid)
metadata<-vm[, c("voteno", "date", "Bill")]
votes_only<-vm[, !names(vm) %in% c("voteno", "date", "Bill")]
vm_t<-as.data.frame(t(votes_only))
colnames(vm_t)<-vm$voteno
vm_t$mpid_key<-rownames(vm_t)
rownames(vm_t)<-NULL
vm_merged<-merge(mp_info, vm_t, by="mpid_key")
View(vm_merged)

vm_merged$surname[vm_merged$firstname=="Brendan"&vm_merged$surname=="O&"]<-"O'Hara"
vm_merged$surname[vm_merged$firstname=="Neil"&vm_merged$surname=="O&"]<-"O'Brien"
vm_merged$party[vm_merged$firstname=="Brendan"&vm_merged$surname=="O'Hara"]<-"SNP"
vm_merged$party[vm_merged$firstname=="Neil"&vm_merged$surname=="O'Brien"]<-"Con"
str(vm_merged)

#EDA
length(vm_merged)
dim(vm_merged)
sum(vm_merged[, -(1:5)]==0, na.rm=TRUE)
sum(vm_merged[, -(1:5)]==1, na.rm=TRUE)
sum(is.na(vm_merged[, -(1:5)]))
mean(is.na(vm_merged[, -(1:5)]))

nrow(mp_info)
table(vm_merged$party)
vm_merged[duplicated(paste(vm_merged$firstname, vm_merged$surname))|
            duplicated(paste(vm_merged$firstname, vm_merged$surname), fromLast=TRUE),
          c("mpid","firstname","surname","party")]

mp_votes<-rowSums(!is.na(vm_merged[, -(1:5)]))
pdf(file="votes_cast.pdf", width=8, height=6)
hist(mp_votes, main="Figure 1.1 - Number of Votes Cast per MP", xlab="Number of Votes Cast", xlim=c(0, 200))
dev.off()

aye_rate<-rowMeans(vm_merged[, -(1:5)]==1, na.rm=TRUE)
pdf(file="aye_votes.pdf", width=8, height=6)
hist(aye_rate, main="Figure 1.2 - Proportion of Aye Votes by MP", xlab="Aye Rate", ylim=c(0, 500))
dev.off()

vm_parties<-table(vm_merged$party)
vm_parties
sum(vm_parties)
as.data.frame(vm_parties)
party_colours<-c("Alliance"="gold", "Con"="royalblue",
                 "DUP"="darkred", "Green"="green", "Independent"="grey", "Lab"="red",
                 "LDem"="orange", "PC"="darkgreen", "Reform UK"="cyan", "SDLP"="olivedrab",
                 "SNP"="yellow", "Traditional Unionist Voice"="purple4", "UUP"="navy")
pdf(file="seats_distribution.pdf", width=8, height=6)
barplot(vm_parties, main="Figure 1.3 - Seats Distribution By Parties", xlab="Party", ylab="Number of MPs",
        col=party_colours[names(vm_parties)], las=2)
dev.off()

#Multidimensional Scaling
votes_mds<-vm_merged[, -(1:5)]
for(i in 1:ncol(votes_mds)){
  if(is.numeric(votes_mds[[i]])){
    votes_mds[[i]][is.na(votes_mds[[i]])]<-mean(votes_mds[[i]], na.rm=TRUE)
  }
}
distance_votes<-dist(votes_mds)
mds2<-cmdscale(distance_votes, k=2)
mds2_matrix<-as.matrix(mds2)
df_mds2<-data.frame(Name=paste(vm_merged$firstname, vm_merged$surname),
                    Party=vm_merged$party,
                    xdim=mds2_matrix[, 1], ydim=mds2_matrix[, 2])

pdf(file="mds2.pdf", width=8, height=6)
ggplot(df_mds2, aes(x=xdim, y=ydim, colour=Party))+
  geom_point(size=2)+
  scale_colour_manual(values=party_colours)+
  theme_minimal()+
  labs(title="Figure 1.4 - Multi-Dimensional Scaling (2D) for MP Voting Records")
dev.off()

mds3<-cmdscale(distance_votes, k=3)
mds3_matrix<-as.matrix(mds3)
df_mds3<-data.frame(Name=paste(vm_merged$firstname, vm_merged$surname),
                    Party=vm_merged$party,
                    xdim=mds3_matrix[, 1], ydim=mds3_matrix[, 2], zdim=mds3_matrix[, 3])
plot_ly(df_mds3, x=~xdim, y=~ydim, z=~zdim, color=~Party, colors=party_colours,
        text=~Name, type="scatter3d", mode="markers", marker=list(size=3),
        title="Figure 1.5 - Multi-Dimensional Scaling (3D) for MP Voting Records")




#Initial IRT model
#Preparation
jagsvotes<-as.matrix(vm_merged[, -(1:5)])
N<-nrow(jagsvotes)
J<-ncol(jagsvotes)
votescast<-rowSums(!is.na(jagsvotes))
orderedmps<-order(votescast, decreasing=TRUE)
fixedjags<-jagsvotes[orderedmps, , drop=FALSE]
jagsdata<-list(N=N, J=J,Y=fixedjags)
vm_fixed<-vm_merged[orderedmps, ]


#1PL:
opl_model<-"model{for(i in 1:N){
theta[i]~dnorm(0,1)}
for(j in 1:J){
b[j]~dnorm(0,1)}
for(i in 1:N){
for(j in 1:J){Y[i,j]~dbern(p[i,j])
logit(p[i,j])<-theta[i]-b[j]}}}"

opl_data<-list(Y=fixedjags, N=N, J=J)

opl_mcmc<-jags.model(textConnection(opl_model), data=opl_data, n.chains=4, n.adapt=2000)
update(opl_mcmc, 5000)
opl_posterior<-coda.samples(opl_mcmc, variable.names=c("theta", "b"), n.iter=10000)

summary(opl_posterior)
gelman.diag(opl_posterior)
effectiveSize(opl_posterior)
summary(effectiveSize(opl_posterior))

opl_summary<-summary(opl_posterior)$statistics
opl_theta<-opl_summary[grep("^theta", rownames(opl_summary)), "Mean"]
opl_b<-opl_summary[
  grep("^b", rownames(opl_summary)), "Mean"]

opl_votes_cast<-rowSums(!is.na(opl_data$Y))
opl_ideology<-data.frame(MP=paste(vm_fixed$firstname,
                                  vm_fixed$surname),
                         Party=vm_fixed$party, Votes=opl_votes_cast,
                         Ideology=opl_theta)
View(opl_ideology)

opl_ideology<-opl_ideology[order(opl_ideology$Ideology), ]
View(opl_ideology)

pdf(file="1pl_graph.pdf", width=8, height=6)
ggplot(opl_ideology, aes(x=Ideology, y=reorder(MP, Ideology), colour=Party))+
  geom_point(size=2)+
  scale_colour_manual(values=party_colours)+
  labs(x="Estimated Ideology", y="", title="One-Dimensional IRT (1PL) Ideology Estimates")+
  theme_minimal()+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
dev.off()

pdf(file="1pl_jitter_graph.pdf", width=8, height=6)
ggplot(opl_ideology, aes(Ideology, Party, colour=Party))+
  scale_colour_manual(values=party_colours)+
  geom_jitter(height=.2)
dev.off()


#2PL:
tpl_model<-"model{for(i in 1:N){
theta[i]~dnorm(0,1)}
for(j in 1:J){
b[j]~dnorm(0,1)
a[j]~dlnorm(0,4)}
for(i in 1:N){
for(j in 1:J){Y[i,j]~dbern(p[i,j])
logit(p[i,j])<-a[j]*theta[i]-b[j]}}}"

tpl_data<-list(Y=fixedjags, N=N, J=J)

tpl_mcmc<-jags.model(textConnection(tpl_model), data=tpl_data, n.chains=4, n.adapt=2000)
update(tpl_mcmc, 5000)
tpl_posterior<-coda.samples(tpl_mcmc, variable.names=c("theta", "a", "b"), n.iter=10000)

summary(tpl_posterior)
gelman.diag(tpl_posterior)
effectiveSize(tpl_posterior)
summary(effectiveSize(tpl_posterior))

tpl_summary<-summary(tpl_posterior)$statistics
tpl_theta<-tpl_summary[grep("^theta", rownames(tpl_summary)), "Mean"]
tpl_a<-tpl_summary[grep("^a", rownames(tpl_summary)), "Mean"]
tpl_b<-tpl_summary[grep("^b", rownames(tpl_summary)), "Mean"]

tpl_votes_cast<-rowSums(!is.na(tpl_data$Y))
tpl_ideology<-data.frame(MP=paste(vm_fixed$firstname, vm_fixed$surname),
                         Party=vm_fixed$party, Votes=tpl_votes_cast,
                         Ideology=tpl_theta)
View(tpl_ideology)

tpl_ideology<-tpl_ideology[order(tpl_ideology$Ideology), ]
View(tpl_ideology)

tpl_vote_parameters<-data.frame(Vote=as.numeric(colnames(fixedjags)),
                                Discrimination=tpl_a,
                                Difficulty=tpl_b)
tpl_vote_parameters<-merge(tpl_vote_parameters, metadata, by.x="Vote",
                           by.y="voteno", all.x=TRUE)
View(tpl_vote_parameters)

pdf(file="2pl_graph.pdf", width=8, height=6)
ggplot(tpl_ideology, aes(x=Ideology, y=reorder(MP, Ideology), colour=Party))+
  geom_point(size=2)+
  scale_colour_manual(values=party_colours)+
  labs(x="Estimated Ideology", y="", title="One-Dimensional IRT (2PL) Ideology Estimates")+
  theme_minimal()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

pdf(file="2pl.parameters.pdf", width=8, height=6)
ggplot(tpl_vote_parameters, aes(x=Difficulty, y=Discrimination))+
  geom_point()+
  theme_minimal()+
  labs(title="Vote Parameters", x="Difficulty (b)", y="Discrimination (a)")
dev.off()





#Multidimensional Models:
#2D IRT:
#Anchor Selection:
valid_votes<-apply(fixedjags, 2, function(x) length(unique(na.omit(x)))>1)
Y<-fixedjags[, valid_votes, drop=FALSE]
vote_n<-colSums(!is.na(Y))
aye<-colMeans(Y, na.rm=TRUE)
vote_quality<-vote_n*(1-2*abs(aye-0.5))
vote_cor<-cor(Y, use="pairwise.complete.obs")
vote_cor[is.na(vote_cor)]<-1
pairs<-combn(ncol(Y), 2)
scores<-apply(pairs, 2, function(pair) {
  vote_quality[pair[1]]*vote_quality[pair[2]]*(1-abs(vote_cor[pair[1], pair[2]]))
})
anchors_2d<-which(valid_votes)[pairs[, which.max(scores)]]
anchors_2d
#Anchors are 14 and 51
#Model
mirt2_model <- "model{for(i in 1:N){
theta1[i]~dnorm(0,1)
theta2[i]~dnorm(0,1)}
a1[14]<-1
a2[14]<-0
b[14]~dnorm(0,1)
a1[51]<-0
a2[51]<-1
b[51]~dnorm(0,1)
for(j in 1:13){b[j]~dnorm(0,1)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)}
for(j in 15:50){b[j]~dnorm(0,1)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)}
for(j in 52:J){b[j]~dnorm(0,1)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)}
for(i in 1:N){
for(j in 1:J){Y[i,j] ~ dbern(p[i,j])
logit(p[i,j])<-(a1[j]*theta1[i])+(a2[j]*theta2[i])-b[j]}}}"

mirt2_full_data<-list(Y=fixedjags, N=nrow(fixedjags), J=ncol(fixedjags))

mirt2_mcmc_full<-jags.model(textConnection(mirt2_model), data=mirt2_full_data, n.chains=3, n.adapt=2000)
update(mirt2_mcmc_full, 5000)
mirt2_posterior_full<-coda.samples(mirt2_mcmc_full, variable.names=c("theta1", "theta2", "a1", "a2", "b"), n.iter=5000)

summary(mirt2_posterior_full)
mirt2_unfixed<-names(which(apply(mirt2_posterior_full[[1]],2,var)>0&
                             apply(mirt2_posterior_full[[2]],2,var)>0&
                             apply(mirt2_posterior_full[[3]],2,var)>0))
gelman.diag(mirt2_posterior_full[,mirt2_unfixed])
effectiveSize(mirt2_posterior_full)
summary(effectiveSize(mirt2_posterior_full[,mirt2_unfixed]))

mirt2_summary_full<-summary(mirt2_posterior_full)$statistics
mirt2_theta1_full<-mirt2_summary_full[grep("^theta1", rownames(mirt2_summary_full)), "Mean"]
mirt2_theta2_full<-mirt2_summary_full[grep("^theta2", rownames(mirt2_summary_full)), "Mean"]
mirt2_a1_full<-mirt2_summary_full[grep("^a1", rownames(mirt2_summary_full)), "Mean"]
mirt2_a2_full<-mirt2_summary_full[grep("^a2", rownames(mirt2_summary_full)), "Mean"]
mirt2_b_full<-mirt2_summary_full[grep("^b", rownames(mirt2_summary_full)), "Mean"]

mirt2_votes_cast<-rowSums(!is.na(mirt2_full_data$Y))
mirt2_ideology_full<-data.frame(MP=paste(vm_fixed$firstname, 
                                         vm_fixed$surname),
                                Party=vm_fixed$party,
                                Votes=mirt2_votes_cast,
                                Theta1=mirt2_theta1_full,
                                Theta2=mirt2_theta2_full)
mirt2_ideology_full <- mirt2_ideology_full[order(mirt2_ideology_full$Theta1),]
View(mirt2_ideology_full)

mirt2_vote_parameters_full<-data.frame(Vote=1:ncol(fixedjags),
                                       Discrimination1=mirt2_a1_full,
                                       Discrimination2=mirt2_a2_full,
                                       Difficulty=mirt2_b_full)

mirt2_vote_parameters_full<-merge(mirt2_vote_parameters_full, metadata,
                                  by.x="Vote", by.y="voteno", all.x=TRUE)
View(mirt2_vote_parameters_full)

top5_dimension1<-mirt2_vote_parameters_full[
  order(-abs(mirt2_vote_parameters_full$Discrimination1)), ][1:5, ]
top5_dimension2<-mirt2_vote_parameters_full[
  order(-abs(mirt2_vote_parameters_full$Discrimination2)), ][1:5, ]
View(top5_dimension1)
View(top5_dimension2)


pdf(file="mirt2_graph.pdf", width=8, height=6)
ggplot(mirt2_ideology_full, aes(x=Theta1, y=Theta2, colour=Party))+
  geom_point()+
  scale_colour_manual(values=party_colours)+
  labs(title="2D-MIRT Estimate Graph")
dev.off()


#3D IRT:
#3D Anchor Selection:
anchor_threshold<-quantile(vote_n, 0.50)
candidates<-which(vote_n>=anchor_threshold)
triplets<-combn(candidates, 3)
scores<-apply(triplets, 2, function(x) {
  correlations<-abs(vote_cor[x, x][upper.tri(vote_cor[x, x])])
  mean(correlations)
})
best_triplet<-triplets[, which.min(scores)]
anchors_3d<-which(valid_votes)[best_triplet]
anchors_3d
round(vote_cor[best_triplet, best_triplet], 3)
#ANCHORS ARE 10, 51, 194


mirt3_model<-"model{for(i in 1:N){
theta1[i]~dnorm(0,1)
theta2[i]~dnorm(0,1)
theta3[i]~dnorm(0,1)}
a1[10]<-1
a2[10]<-0
a3[10]<-0
b[10]~dnorm(0,1)
a1[51]<-0
a2[51]<-1
a3[51]<-0
b[51]~dnorm(0,1)
a1[194]<-0
a2[194]<-0
a3[194]<-1
b[194]~dnorm(0,1)
for(j in 1:9){b[j]~dnorm(0,1)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)
a3[j]~dnorm(0,4)}
for(j in 11:50){b[j] ~ dnorm(0,1)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)
a3[j]~dnorm(0,4)}
for(j in 52:193){b[j]~dnorm(0,1)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)
a3[j]~dnorm(0,4)}
for(j in 195:J){b[j]~dnorm(0,1)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)
a3[j]~dnorm(0,4)}
for(i in 1:N){
for(j in 1:J){Y[i,j]~dbern(p[i,j])
logit(p[i,j])<-a1[j]*theta1[i]+a2[j]*theta2[i]+a3[j]*theta3[i]-b[j]}}}"

mirt3_full_data<-list(Y=fixedjags, N=nrow(fixedjags), J=ncol(fixedjags))

mirt3_mcmc_full<-jags.model(textConnection(mirt3_model), data=mirt3_full_data, n.chains=3, n.adapt=2000)
update(mirt3_mcmc_full, 5000)
mirt3_posterior_full<-coda.samples(mirt3_mcmc_full, variable.names=c("theta1","theta2","theta3","a1","a2","a3","b"), n.iter=5000)

summary(mirt3_posterior_full)
mirt3_unfixed<-names(which(apply(mirt3_posterior_full[[1]],2,var)>0&
                             apply(mirt3_posterior_full[[2]],2,var)>0&
                             apply(mirt3_posterior_full[[3]],2,var)>0))
gelman.diag(mirt3_posterior_full[, mirt3_unfixed])
effectiveSize(mirt3_posterior_full)
summary(effectiveSize(mirt3_posterior_full[, mirt3_unfixed]))

mirt3_summary_full<-summary(mirt3_posterior_full)$statistics
mirt3_theta1_full<-mirt3_summary_full[grep("^theta1",rownames(mirt3_summary_full)), "Mean"]
mirt3_theta2_full<-mirt3_summary_full[grep("^theta2",rownames(mirt3_summary_full)), "Mean"]
mirt3_theta3_full<-mirt3_summary_full[grep("^theta3",rownames(mirt3_summary_full)), "Mean"]
mirt3_a1_full<-mirt3_summary_full[grep("^a1",rownames(mirt3_summary_full)), "Mean"]
mirt3_a2_full<-mirt3_summary_full[grep("^a2",rownames(mirt3_summary_full)), "Mean"]
mirt3_a3_full<-mirt3_summary_full[grep("^a3",rownames(mirt3_summary_full)), "Mean"]
mirt3_b_full<-mirt3_summary_full[grep("^b",rownames(mirt3_summary_full)), "Mean"]

mirt3_votes_cast<-rowSums(!is.na(mirt3_full_data$Y))
mirt3_ideology_full<-data.frame(MP=paste(vm_fixed$firstname,
                                         vm_fixed$surname),
                                Party=vm_fixed$party,
                                Votes=mirt3_votes_cast, Theta1=mirt3_theta1_full,
                                Theta2=mirt3_theta2_full,
                                Theta3=mirt3_theta3_full)
mirt3_ideology_full<-mirt3_ideology_full[order(mirt3_ideology_full$Theta1),]
View(mirt3_ideology_full)

mirt3_vote_parameters_full<-data.frame(Vote=1:ncol(fixedjags), Discrimination1=mirt3_a1_full, Discrimination2=mirt3_a2_full, 
                                       Discrimination3=mirt3_a3_full, Difficulty=mirt3_b_full)
View(mirt3_vote_parameters_full)

mirt3_vote_parameters_full<-merge(mirt3_vote_parameters_full, metadata, by.x="Vote", by.y="voteno", all.x=TRUE)
mirt3_dimension1<-mirt3_vote_parameters_full[order(-abs(mirt3_vote_parameters_full$Discrimination1)),]
mirt3_dimension2<-mirt3_vote_parameters_full[order(-abs(mirt3_vote_parameters_full$Discrimination2)),]
mirt3_dimension3<-mirt3_vote_parameters_full[order(-abs(mirt3_vote_parameters_full$Discrimination3)),]
View(head(mirt3_dimension1, 20))
View(head(mirt3_dimension2, 20))
View(head(mirt3_dimension3, 20))

plot_ly(mirt3_ideology_full, x=~Theta1, y=~Theta2, z=~Theta3, color=~Party, 
        colors=party_colours, text=~MP, type="scatter3d",
        mode="markers", marker=list(size=2))

pdf(file="mirt3_graph.pdf", width=8, height=6)
ggplot(mirt3_vote_parameters_full, aes(x=Discrimination1, y=Discrimination2, colour=Discrimination3))+
  geom_point(size=2)+
  theme_minimal()+
  labs(title="3D MIRT Vote Discrimination Parameters", x="Dimension 1", y="Dimension 2")
dev.off()




#Non-Linear Models
#Anchors are 14 and 51 again (2D)
#Model:
circ_model<-"model{phi[1]<-3.1415926535
theta1[1]<-cos(phi[1])
theta2[1]<-sin(phi[1])
phi[211]<-1.5708
theta1[211]<-cos(phi[211])
theta2[211]<-sin(phi[211])
phi[221]<-0
theta1[221]<-cos(phi[221])
theta2[221]<-sin(phi[221])
phi[239]<--1.5708
theta1[239]<-cos(phi[239])
theta2[239]<-sin(phi[239])
for(i in 2:210){phi[i]~dunif(-3.1415926535, 3.1415926535)
theta1[i]<-cos(phi[i])
theta2[i]<-sin(phi[i])}
for(i in 212:220){phi[i]~dunif(-3.1415926535, 3.1415926535)
theta1[i]<-cos(phi[i])
theta2[i]<-sin(phi[i])}
for(i in 222:238){phi[i]~dunif(-3.1415926535, 3.1415926535)
theta1[i]<-cos(phi[i])
theta2[i]<-sin(phi[i])}
for(i in 240:N){phi[i]~dunif(-3.1415926535, 3.1415926535)
theta1[i]<-cos(phi[i])
theta2[i]<-sin(phi[i])}
a1[14]<-1
a2[14]<-0
b[14]~dnorm(0,0.2)
a1[51]<-0
a2[51]<-1
b[51]~dnorm(0,0.2)
for(j in 1:13){b[j]~dnorm(0,0.2)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)}
for(j in 15:50){b[j] ~ dnorm(0,0.2)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)}
for(j in 52:J){b[j]~dnorm(0,0.2)
a1[j]~dnorm(0,4)
a2[j]~dnorm(0,4)}
for(i in 1:N){
for(j in 1:J){Y[i,j]~dbern(p[i,j])
logit(p[i,j])<-a1[j]*theta1[i]+a2[j]*theta2[i]-b[j]}}}"

circ_data<-list(Y=fixedjags, N=nrow(fixedjags), J=ncol(fixedjags))

circ_mcmc<-jags.model(textConnection(circ_model), data=circ_data, n.chains=2, n.adapt=2000)
update(circ_mcmc, 5000)
circ_posterior<-coda.samples(circ_mcmc,
                             variable.names=c("phi","theta1","theta2","a1","a2","b"),
                             n.iter=5000)

summary(circ_posterior)
circ_unfixed<-setdiff(grep("^(a1|a2|b|phi|theta1|theta2)",
                           varnames(circ_posterior),
                           value=TRUE),
                      c("a1[14]","a2[14]", "a1[51]","a2[51]"))
circ_unfixed<-names(which(apply(circ_posterior[[1]], 2, var)>0&
                            apply(circ_posterior[[2]], 2, var)>0))
gelman.diag(circ_posterior[, circ_unfixed])

rhat<-sapply(circ_unfixed, function(x)
  gelman.diag(circ_posterior[,x],
              multivariate=FALSE,
              autoburnin=FALSE)$psrf[1,"Point est."])

max(rhat)
summary(rhat)
sort(rhat, decreasing=TRUE)[1:20]

effectiveSize(circ_posterior[,circ_unfixed])
summary(effectiveSize(circ_posterior[,circ_unfixed]))

circ_summary<-summary(circ_posterior)$statistics
circ_phi<-circ_summary[grep("^phi",rownames(circ_summary)),"Mean"]
circ_a1<-circ_summary[grep("^a1",rownames(circ_summary)),"Mean"]
circ_a2<-circ_summary[grep("^a2",rownames(circ_summary)),"Mean"]
circ_b<-circ_summary[grep("^b",rownames(circ_summary)),"Mean"]
circ_theta1<-cos(circ_phi)
circ_theta2<-sin(circ_phi)


circ_ideology<-data.frame(MP=paste(vm_fixed$firstname,vm_fixed$surname),
                          Party=vm_fixed$party,
                          Phi=circ_phi,
                          Theta1=circ_theta1,
                          Theta2=circ_theta2)

View(circ_ideology)

circ_votes_cast<-rowSums(!is.na(fixedjags))
data.frame(MP=circ_ideology$MP,
           Party=circ_ideology$Party,
           Votes=circ_votes_cast,
           Phi=circ_phi)


circ_vote_parameters<-data.frame(Vote=as.numeric(colnames(fixedjags)),
                                 Discrimination1=circ_a1,
                                 Discrimination2=circ_a2,
                                 Difficulty=circ_b)
circ_vote_parameters<-merge(circ_vote_parameters,
                            metadata,
                            by.x="Vote",
                            by.y="voteno",
                            all.x=TRUE)

circ_dimension1<-circ_vote_parameters[
  order(-abs(circ_vote_parameters$Discrimination1)),]
circ_dimension2<-circ_vote_parameters[
  order(-abs(circ_vote_parameters$Discrimination2)),]

View(head(circ_dimension1,20))
View(head(circ_dimension2,20))

circle<-data.frame(x=cos(seq(0,2*pi,length.out=250)),
                   y=sin(seq(0,2*pi,length.out=250)))

pdf(file="circular_graph.pdf",width=8,height=6)
ggplot()+
  geom_path(data=circle,aes(x,y),linewidth=0.4)+
  geom_point(data=circ_ideology, aes(Theta1,Theta2,colour=Party), size=2)+
  coord_equal()+
  scale_colour_manual(values=party_colours)+
  theme_minimal()+
  labs(title="Circular IRT Graph", x="Dimension 1", y="Dimension 2")
dev.off()


#Spherical:
#Anchors are 10, 51 and 194 (3D again)
#Model:
sphere_model <- "model{phi[1]<-3.1415926535
psi[1]<-0
theta1[1]<-cos(psi[1])*cos(phi[1])
theta2[1]<-cos(psi[1])*sin(phi[1])
theta3[1]<-sin(psi[1])
phi[221]<-0
psi[221]<-0
theta1[221]<-cos(psi[221])*cos(phi[221])
theta2[221]<-cos(psi[221])*sin(phi[221])
theta3[221]<-sin(psi[221])
phi[409]<-1.5708
psi[409]<-0
theta1[409]<-cos(psi[409])*cos(phi[409])
theta2[409]<-cos(psi[409])*sin(phi[409])
theta3[409]<-sin(psi[409])
for(i in 2:220){phi[i]~dunif(-3.1415926535, 3.1415926535)
psi[i]~dunif(-1.5708, 1.5708)
theta1[i]<-cos(psi[i])*cos(phi[i])
theta2[i]<-cos(psi[i])*sin(phi[i])
theta3[i]<-sin(psi[i])}
for(i in 222:408){phi[i]~dunif(-3.1415926535, 3.1415926535)
psi[i]~dunif(-1.5708, 1.5708)
theta1[i]<-cos(psi[i])*cos(phi[i])
theta2[i]<-cos(psi[i])*sin(phi[i])
theta3[i]<-sin(psi[i])}
for(i in 410:N){phi[i]~dunif(-3.1415926535, 3.1415926535)
psi[i]~dunif(-1.5708, 1.5708)
theta1[i]<-cos(psi[i])*cos(phi[i])
theta2[i]<-cos(psi[i])*sin(phi[i])
theta3[i]<-sin(psi[i])}
a1[10]<-1
a2[10]<-0
a3[10]<-0
b[10]~dnorm(0, 0.2)
a1[51]<-0
a2[51]<-1
a3[51]<-0
b[51]~dnorm(0, 0.2)
a1[194]<-0
a2[194]<-0
a3[194]<-1
b[194]~dnorm(0, 0.2)
for(j in 1:9){b[j]~dnorm(0, 0.2)
a1[j]~dnorm(0, 4)
a2[j]~dnorm(0, 4)
a3[j]~dnorm(0, 4)}
for(j in 11:50){b[j]~dnorm(0, 0.2)
a1[j]~dnorm(0, 4)
a2[j]~dnorm(0, 4)
a3[j]~dnorm(0, 4)}
for(j in 52:193){b[j]~dnorm(0, 0.2)
a1[j]~dnorm(0, 4)
a2[j]~dnorm(0, 4)
a3[j]~dnorm(0, 4)}
for(j in 195:J){b[j]~dnorm(0, 0.2)
a1[j]~dnorm(0, 4)
a2[j]~dnorm(0, 4)
a3[j]~dnorm(0, 4)}
for(i in 1:N){
for(j in 1:J){Y[i,j]~dbern(p[i,j])
logit(p[i,j])<-a1[j]*theta1[i]+a2[j]*theta2[i]+a3[j]*theta3[i]-b[j]}}}"

sphere_data<-list(Y=fixedjags, N=nrow(fixedjags), J=ncol(fixedjags))
sphere_mcmc<-jags.model(textConnection(sphere_model), data=sphere_data, n.chains=4, n.adapt=2000)
update(sphere_mcmc, 5000)
sphere_posterior<-coda.samples(sphere_mcmc, variable.names=c("phi", "psi", "a1", "a2", "a3", "b"), n.iter=5000)

summary(sphere_posterior)
effectiveSize(sphere_posterior)
sphere_unfixed<-names(which(apply(sphere_posterior[[1]], 2, var)>0&
                              apply(sphere_posterior[[2]], 2, var)>0&
                              apply(sphere_posterior[[3]], 2, var)>0&
                              apply(sphere_posterior[[4]], 2, var)>0))
gelman.diag(sphere_posterior[, sphere_unfixed])
sphere_summary<-summary(sphere_posterior)$statistics
sphere_phi<-sphere_summary[grep("^phi", rownames(sphere_summary)), "Mean"]
sphere_psi<-sphere_summary[grep("^psi", rownames(sphere_summary)), "Mean"]
sphere_a1<-sphere_summary[grep("^a1",  rownames(sphere_summary)), "Mean"]
sphere_a2<-sphere_summary[grep("^a2",  rownames(sphere_summary)), "Mean"]
sphere_a3<-sphere_summary[grep("^a3",  rownames(sphere_summary)), "Mean"]
sphere_b<-sphere_summary[grep("^b",   rownames(sphere_summary)), "Mean"]

sphere_theta1<-cos(sphere_psi)*cos(sphere_phi)
sphere_theta2<-cos(sphere_psi)*sin(sphere_phi)
sphere_theta3<-sin(sphere_psi)

sphere_ideology<-data.frame(MP=paste(vm_fixed$firstname, vm_fixed$surname),
                            Party=vm_fixed$party,
                            Phi=sphere_phi,
                            Psi=sphere_psi,
                            Theta1=sphere_theta1,
                            Theta2=sphere_theta2,
                            Theta3=sphere_theta3)
View(sphere_ideology)

sphere_votes_cast<-rowSums(!is.na(fixedjags))
data.frame(MP=sphere_ideology$MP,Party=sphere_ideology$Party,
           Votes=sphere_votes_cast,
           Phi=sphere_phi,
           Psi=sphere_psi)

sphere_vote_parameters<-data.frame(Vote=as.numeric(colnames(fixedjags)),
                                   Discrimination1=sphere_a1,
                                   Discrimination2=sphere_a2,
                                   Discrimination3=sphere_a3,
                                   Difficulty=sphere_b)

sphere_vote_parameters<-merge(sphere_vote_parameters, metadata,
                              by.x="Vote", by.y="voteno", all.x=TRUE)
sphere_dimension1<-sphere_vote_parameters[order(-abs(sphere_vote_parameters$Discrimination1)), ]
sphere_dimension2<-sphere_vote_parameters[order(-abs(sphere_vote_parameters$Discrimination2)), ]
sphere_dimension3<-sphere_vote_parameters[order(-abs(sphere_vote_parameters$Discrimination3)), ]

View(head(sphere_dimension1, 20))
View(head(sphere_dimension2, 20))
View(head(sphere_dimension3, 20))

plot_ly(sphere_ideology, x=~Theta1, y=~Theta2, z=~Theta3,
        color=~Party, colors=party_colours, text=~MP, type="scatter3d",
        mode="markers", marker=list(size=3))

#END OF CODE