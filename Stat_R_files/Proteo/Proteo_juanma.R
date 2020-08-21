library("nlme") #pour la fonction lme
options(width=160)
library(plyr)
library(reshape2)

library(ggplot2)
library(Hmisc)
library(IDPmisc)


# Sans pepitides en commun  -----------------------------------------------
setwd("/media/juanma/JUANMA/Stage M2/git_repo/Proteo/")

# lecture du fichier de sortie quanti de masschroq
tab.q=read.table("peptides_q1_fractiona1.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#v?rification de la structure et du nb de peptides
str(tab.q) 
length(unique(tab.q$peptide)) 
# cr?ation de la colonne peptide-charge
tab.q$peptiz=paste(tab.q$peptide,tab.q$z,sep='-') 

tab.prot=read.table("fractiona1_proteins.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# v?rification de la structure (nb de lignes) et du nb de peptides
str(tab.prot) #23025 obs 3var
length(unique(tab.prot$peptide)) #20041 length(unique(tab.prot$peptide))=length(unique(tab.q$peptide)) on peut continuer

# homog?n?isation de tab.q et tab.prot
tab.prot=tab.prot[tab.prot$peptide %in% unique(tab.q$peptide),]		
length(unique(tab.prot$peptide)) #20041, comme dans tab.q	

msrunfile=unique(tab.q$msrun)
msrunfile=msrunfile[order(msrunfile)]

Stade=rev(c("bulk","bulk","bulk","bulk",rep(c("0 DPA","13 DPA","26 DPA","39 DPA","55 DPA","76 DPA","118 DPA","179 DPA","222 DPA"),each=3)))
Rep=rev(c(1,2,3,4,rep(c(1,2,3),9)))

stade_rep=  rev(c("bulk_1","bulk_2","bulk_3","bulk_4","0_1","0_2","0_3","13_1","13_2","13_3","26_1","26_2","26_3","39_1","39_2","39_3",
                  "55_1","55_2","55_3","76_1","76_2","76_3","118_1","118_2","118_3","179_1","179_2","179_3","222_1","222_2","222_3"))

indice=match(tab.q$msrun, msrunfile)
tab.q$Stade=Stade[indice]
tab.q$Rep=Rep[indice]
tab.q$stade_rep=stade_rep[indice]


################## calcul du biais=ref-obj
sample=unique(tab.q$msrun) #29
ref.sample="msruna1"
nbrun=length(sample) #29


################## Lissage du biais
mat.obj=tapply(tab.q$logarea,list(tab.q$peptiz,tab.q$msrun),FUN=sum)
ref=mat.obj[,ref.sample]
mat.ref=as.data.frame(matrix(rep(ref,nbrun),ncol=nbrun)) 
dim(mat.ref) #23383obs    29var

tab.biais=mat.obj-mat.ref
colnames(tab.biais) <- colnames(mat.obj)
tab.biais=stack(tab.biais)
biais.peptiz <- rep(rownames(mat.obj),nbrun)
tab.biais=cbind.data.frame(biais.peptiz, tab.biais)
colnames(tab.biais)=c("peptiz","biais","msrun")
dim(tab.biais)#678107 x 3
length(unique(tab.biais$peptiz)) #23383

indice=match(paste(tab.biais$msrun, tab.biais$peptiz), paste(tab.q$msrun, tab.q$peptiz))
tab.biais$rt=tab.q$rt[indice]

# on ordonne les lignes de tab.biais par msrun et par rt

tab.biais=tab.biais[order(tab.biais$msrun, tab.biais$rt),]
tab.biais=tab.biais[!is.na(tab.biais$rt),]#v?rification :mm nb de lignes que tab.q
tab.biais$correc=NA
compteur=unique(tab.biais$msrun)

#pdf("./figures/Normalisation_biais.pdf")
for (i in 1:length(compteur)) {
  #calcul de la correction
  select.msrun= tab.biais$msrun==compteur[i]	
  select.na <- is.na(tab.biais$biais)
  select.inf <- is.infinite(tab.biais$biais)
  
  cour <- tab.biais[select.msrun,]
  cour.nomiss <- tab.biais[select.msrun & !select.na & !select.inf,]
  
  # fonction de lissage de la courbe biais en fonction du rt
  yyy=smooth.spline(cour.nomiss$rt,cour.nomiss$biais,spar=0.5)
  # on ajoute le biais liss? (qui est le facteur de correction pour la normalisation) ? tab.biais
  indice=match(tab.biais$rt[select.msrun & !select.na], yyy$x)
  tab.biais$correc[select.msrun & !select.na] = yyy$y[indice]
  
  # repr?sentation graphique du biais et de la courbe de lissage				
  rt.min=min(cour.nomiss$rt)
  rt.max=max(cour.nomiss$rt)
  biais.min=min(cour.nomiss$biais)
  biais.max=max(cour.nomiss$biais)
  plot(cour.nomiss$rt,cour.nomiss$biais,main=paste("msrun = ", compteur[i], sep=""),xlab="RT",ylab="Biais (log)",xlim=c(rt.min,rt.max),ylim=c(biais.min,biais.max))
  lines(yyy$x,yyy$y,col="red",type="l")
  
  
  #calcul de la correction pour les peptides pr?sents dans le msrun ? normaliser et absents du msrun de r?f?rence: on prend la valeur de correction du peptide dont le temps de r?tention est juste avant
  smo=tab.biais$correc[select.msrun]
  j=1
  while(is.na(smo[j])) {smo[j]=0; j=j+1}
  while (j < length(smo)) {
    while(!is.na(smo)[j] & j<length(smo)) {last=smo[j]; j=j+1}
    while(is.na(smo[j]) & j<length(smo)) {smo[j]=last; j=j+1}
  }  
  smo[length(smo)]=smo[length(smo)-1] 
  tab.biais$correc[select.msrun] =smo
  
}
dev.off()

################## reunir l'info logQ et correc et tab.q
indice=match(paste(tab.q$msrun, tab.q$peptiz), paste(tab.biais$msrun, tab.biais$peptiz))
tab.q$correc=tab.biais$correc[indice]
head(tab.q)
dim(tab.q) # 549754 on v?rifie qu'on a bien le nb de ligne attendu (comme le tab.biais)

# on corrige les intensit?s de peptiz avec le facteur de correction	
tab.q$logQnorm=tab.q$logarea+tab.q$correc	


# Retire les donn?es des bulks qui ne sont plus n?cessaires ? l'analyse
tab.q=tab.q[-which(tab.q$msrun=="msrunb31"),] #retire bulk1
tab.q=tab.q[-which(tab.q$msrun=="msrunb30"),] #retire bulk2
tab.q=tab.q[-which(tab.q$msrun=="msrunb29"),] #retire bulk3
tab.q=tab.q[-which(tab.q$msrun=="msrunb28"),] #retire bulk3

dim(tab.q)#495691     22var

############# Suppression des peptides communs ? plusieurs proteines

#on compte le nb de prot auxquelles appartient chaque peptide
count = as.data.frame(table(tab.prot$peptide))
colnames(count)=c("peptide","Freq_prot")

# on ne conserve que les peptides qui n'appartiennent qu'? une seule prot?ine
good_peptides=count$peptide[count$Freq_prot==1] #17824 elements (attention, good_peptides est construit avec tab.prot qui n'a pas ?t? pr?alablement homog?n?is? avec tabok)
tab.q=tab.q[tab.q$peptide %in% good_peptides,] # 427255obs 22var

# on ajoute l'info prot?ine dans tab.q (sans risque de dupliquer des lignes puisque maintenant il n'y a plus de peptides communs)
tab.q = merge(tab.q,tab.prot,"peptide", all.x=TRUE, all.y=FALSE)
dim(tab.q) #431047 24var
length(unique(tab.q$peptide)) #17821
length(unique(tab.q$protein)) #2726


##################  Enlever les proteines avec moins de 2 peptides et pas des peptidez 
temp=unique(tab.q[,c("protein", "peptide")])
count=table(temp$protein)
bad.prot=names(count[count<2]) 
tab.q=tab.q[!tab.q$protein %in% bad.prot,] 
dim(tab.q) #426610 24var
length(unique(tab.q$peptide)) #17612
length(unique(tab.q$protein)) #2517

tab.q$Rep <- gsub("4","3",tab.q$Rep) #remplace les 4 pas un 3 pour faciliter la mod?lisation 
tab.q$stade_rep <- gsub("_4","_3",tab.q$stade_rep)

################## vérification des contrastes
options(contrasts=c("contr.sum", "contr.poly")) #on modifie les contrastes par défaut

options("contrasts") #on vérifie qu'on a bien les bons contrastes

################## estimer les abondances de prot
proteines=unique(tab.q$protein)

tab.peptiz=NULL
tab.rep=NULL
tab.effect.msrun=NULL
tab.residue=NULL
tab.intercept=NULL
tab.stade=NULL

for (i in 1:length(proteines)){
  tryCatch({ #tryCatch sert à attrapper les messages d'erreur sans que les calculs de la boucle de n'interrompent
    myprot=proteines[i]
    sub= tab.q[tab.q$protein==myprot,]
    sub$Stade=as.factor(sub$Stade)
    sub$Rep=as.factor(sub$Rep)
    npep=length(unique(sub$peptiz))
    nStade=length(unique(sub$Stade))
    nrep=length(unique(sub$Rep))
    
    model = lme(logQnorm ~ Stade+Rep+peptiz, random=~1|msrun, data=sub ) 
    
    fixe=fixef(model)
    N=length(fixe)
    
    estimation=unique(sub[,c("Stade", "protein")])
    estimation=estimation[order(estimation$Stade),]
    intercept=rep(fixe[1], nrow(estimation))
    estimation$intercept=intercept
    tab.intercept=rbind.data.frame(tab.intercept, estimation)
    
    tab.stade=rbind.data.frame(tab.stade, cbind.data.frame(Stade.effect=c(fixe[2:nStade], -sum(fixe[2:nStade])), Stade=estimation$Stade, protein=rep(myprot, nStade)))
    
    tab.rep=rbind.data.frame(tab.rep, cbind.data.frame(rep.effect=c(fixe[(nStade+1):(nStade+2)], -sum(fixe[(nStade+1):(nStade+2)]) ), Rep=c("1", "2", "3"), protein=rep(myprot, nrep)))
    
    tab.peptiz=rbind.data.frame(tab.peptiz, cbind.data.frame(peptidez.effect=c(fixe[(nStade+nrep):(length(fixe))], -sum(fixe[(nStade+nrep):(length(fixe))])), peptiz=levels(as.factor(as.character(sub$peptiz))), protein=myprot))
    
    effect.msrun=model$coefficients$random$msrun
    tab.effect.msrun=rbind.data.frame(tab.effect.msrun, cbind.data.frame(msrun=rownames(effect.msrun), effet.msrun=effect.msrun[,1], protein=rep(myprot, nrow(effect.msrun))))
    
    tab.residue=rbind.data.frame(tab.residue, cbind.data.frame(msrun=sub$msrun, protein=sub$protein, peptiz=sub$peptiz, residue=resid(model)))
  },error=function(e){cat("Error for protein", myprot, ":", conditionMessage(e), "\n")}) #fin de tryCatch
}

# Error for protein c140.a2.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein c150.a6.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein c583.a2.a1 : nlminb problem, convergence error code = 1
# message = singular convergence (7) 
# Error for protein d1616.a2.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d1776.a2.a1 : the leading minor of order 1 is not positive definite 
# Error for protein d2114.a2.a1 : contrasts can be applied only to factors with 2 or more levels 
# Error for protein d2426.a2.a1 : the leading minor of order 1 is not positive definite 
# Error for protein d2571.a1.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3376.a1.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3407.a1.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3417.a1.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3420.a1.a1 : Singularity in backsolve at level 0, block 1 

tab.q=merge(tab.q, tab.intercept, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.stade, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.rep, by=c("Rep", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.peptiz, by=c("peptiz", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.effect.msrun, by=c("msrun", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.residue, by=c("msrun", "protein", "peptiz"), all.x=TRUE)
#tab.q$ab.estime=apply(tab.q[,c("intercept", "Stade.effect", "rep.effect")], 1, sum) ligne ajout?e qui cr?e des probl?mes

tab.sans.pep.commun=tab.q
# save(tab.sans.pep.commun, file="./data_kiwi_Juan.RData")


# Sans pep communs, sans pep avec sdRT trop grand -------------------------

setwd("/media/juanma/JUANMA/Stage M2/git_repo/Proteo/")

# lecture du fichier de sortie quanti de masschroq
tab.q=read.table("peptides_q1_fractiona1.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#v?rification de la structure et du nb de peptides
str(tab.q) 
length(unique(tab.q$peptide)) 
# cr?ation de la colonne peptide-charge
tab.q$peptiz=paste(tab.q$peptide,tab.q$z,sep='-') 

tab.prot=read.table("fractiona1_proteins.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# v?rification de la structure (nb de lignes) et du nb de peptides
str(tab.prot) #23025 obs 3var
length(unique(tab.prot$peptide)) #20041 length(unique(tab.prot$peptide))=length(unique(tab.q$peptide)) on peut continuer

# homog?n?isation de tab.q et tab.prot
tab.prot=tab.prot[tab.prot$peptide %in% unique(tab.q$peptide),]		
length(unique(tab.prot$peptide)) #20041, comme dans tab.q	

msrunfile=unique(tab.q$msrun)
msrunfile=msrunfile[order(msrunfile)]

Stade=rev(c("bulk","bulk","bulk","bulk",rep(c("0 DPA","13 DPA","26 DPA","39 DPA","55 DPA","76 DPA","118 DPA","179 DPA","222 DPA"),each=3)))
Rep=rev(c(1,2,3,4,rep(c(1,2,3),9)))
stade_rep=  rev(c("bulk_1","bulk_2","bulk_3","bulk_4","0_1","0_2","0_3","13_1","13_2","13_3","26_1","26_2","26_3","39_1","39_2","39_3",
                  "55_1","55_2","55_3","76_1","76_2","76_3","118_1","118_2","118_3","179_1","179_2","179_3","222_1","222_2","222_3"))

indice=match(tab.q$msrun, msrunfile)
tab.q$Stade=Stade[indice]
tab.q$Rep=Rep[indice]
tab.q$stade_rep=stade_rep[indice]

# calcul de l'?cart-type du temps de r?tention 
total_rt <-tapply(tab.q$rt,list(tab.q$peptide),FUN=sd, na.rm=TRUE)
total_rt=total_rt[-which(is.na(total_rt))] # quand il n'y a qu'1 valeur par peptide, il ne peut pas calculer d'?cart-type et renvoie NA
#repr?sentation de la distribution des valeurs
hist(total_rt,nclass=500,freq=T,xlab="Deviation standart en sec", xlim=c(0,100))
# suppression des peptides pr?sentant trop de variabilit?
bad_peptides = names(total_rt[total_rt>20]) #213
tab.q=tab.q[!tab.q$peptide %in% bad_peptides,]
dim(tab.q) #549754     20
length(unique(tab.q$peptide)) #19825

################## calcul du biais=ref-obj
sample=unique(tab.q$msrun) #29
ref.sample="msruna1"
nbrun=length(sample) #29


################## Lissage du biais
mat.obj=tapply(tab.q$logarea,list(tab.q$peptiz,tab.q$msrun),FUN=sum)
ref=mat.obj[,ref.sample]
mat.ref=as.data.frame(matrix(rep(ref,nbrun),ncol=nbrun)) 
dim(mat.ref) #23383obs    29var

tab.biais=mat.obj-mat.ref
colnames(tab.biais) <- colnames(mat.obj)
tab.biais=stack(tab.biais)
biais.peptiz <- rep(rownames(mat.obj),nbrun)
tab.biais=cbind.data.frame(biais.peptiz, tab.biais)
colnames(tab.biais)=c("peptiz","biais","msrun")
dim(tab.biais)#678107 x 3
length(unique(tab.biais$peptiz)) #23383

indice=match(paste(tab.biais$msrun, tab.biais$peptiz), paste(tab.q$msrun, tab.q$peptiz))
tab.biais$rt=tab.q$rt[indice]

# on ordonne les lignes de tab.biais par msrun et par rt

tab.biais=tab.biais[order(tab.biais$msrun, tab.biais$rt),]
tab.biais=tab.biais[!is.na(tab.biais$rt),]#v?rification :mm nb de lignes que tab.q
tab.biais$correc=NA
compteur=unique(tab.biais$msrun)

#pdf("./figures/Normalisation_biais.pdf")
for (i in 1:length(compteur)) {
  #calcul de la correction
  select.msrun= tab.biais$msrun==compteur[i]	
  select.na <- is.na(tab.biais$biais)
  select.inf <- is.infinite(tab.biais$biais)
  
  cour <- tab.biais[select.msrun,]
  cour.nomiss <- tab.biais[select.msrun & !select.na & !select.inf,]
  
  # fonction de lissage de la courbe biais en fonction du rt
  yyy=smooth.spline(cour.nomiss$rt,cour.nomiss$biais,spar=0.5)
  # on ajoute le biais liss? (qui est le facteur de correction pour la normalisation) ? tab.biais
  indice=match(tab.biais$rt[select.msrun & !select.na], yyy$x)
  tab.biais$correc[select.msrun & !select.na] = yyy$y[indice]
  
  # repr?sentation graphique du biais et de la courbe de lissage				
  rt.min=min(cour.nomiss$rt)
  rt.max=max(cour.nomiss$rt)
  biais.min=min(cour.nomiss$biais)
  biais.max=max(cour.nomiss$biais)
  plot(cour.nomiss$rt,cour.nomiss$biais,main=paste("msrun = ", compteur[i], sep=""),xlab="RT",ylab="Biais (log)",xlim=c(rt.min,rt.max),ylim=c(biais.min,biais.max))
  lines(yyy$x,yyy$y,col="red",type="l")
  
  
  #calcul de la correction pour les peptides pr?sents dans le msrun ? normaliser et absents du msrun de r?f?rence: on prend la valeur de correction du peptide dont le temps de r?tention est juste avant
  smo=tab.biais$correc[select.msrun]
  j=1
  while(is.na(smo[j])) {smo[j]=0; j=j+1}
  while (j < length(smo)) {
    while(!is.na(smo)[j] & j<length(smo)) {last=smo[j]; j=j+1}
    while(is.na(smo[j]) & j<length(smo)) {smo[j]=last; j=j+1}
  }  
  smo[length(smo)]=smo[length(smo)-1] 
  tab.biais$correc[select.msrun] =smo
  
}
dev.off()

################## reunir l'info logQ et correc et tab.q
indice=match(paste(tab.q$msrun, tab.q$peptiz), paste(tab.biais$msrun, tab.biais$peptiz))
tab.q$correc=tab.biais$correc[indice]
head(tab.q)
dim(tab.q) # 549754 on v?rifie qu'on a bien le nb de ligne attendu (comme le tab.biais)

# on corrige les intensit?s de peptiz avec le facteur de correction	
tab.q$logQnorm=tab.q$logarea+tab.q$correc	

tab.q=tab.q[-which(tab.q$msrun=="msrunb31"),] #retire bulk1
tab.q=tab.q[-which(tab.q$msrun=="msrunb30"),] #retire bulk2
tab.q=tab.q[-which(tab.q$msrun=="msrunb29"),] #retire bulk3
tab.q=tab.q[-which(tab.q$msrun=="msrunb28"),] #retire bulk3

#on compte le nb de prot auxquelles appartient chaque peptide
count = as.data.frame(table(tab.prot$peptide))
colnames(count)=c("peptide","Freq_prot")

# on ne conserve que les peptides qui n'appartiennent qu'? une seule prot?ine
good_peptides=count$peptide[count$Freq_prot==1] #17824 elements (attention, good_peptides est construit avec tab.prot qui n'a pas ?t? pr?alablement homog?n?is? avec tabok)
tab.q=tab.q[tab.q$peptide %in% good_peptides,] # 427255obs 22var

# on ajoute l'info prot?ine dans tab.q (sans risque de dupliquer des lignes puisque maintenant il n'y a plus de peptides communs)
tab.q = merge(tab.q,tab.prot,"peptide", all.x=TRUE, all.y=FALSE)
dim(tab.q) #427255obs 24var
length(unique(tab.q$peptide)) #17662
length(unique(tab.q$protein)) #2723


##################  Enlever les proteines avec moins de 2 peptides et pas des peptidez 
temp=unique(tab.q[,c("protein", "peptide")])
count=table(temp$protein)
bad.prot=names(count[count<2]) #213 prot
tab.q=tab.q[!tab.q$protein %in% bad.prot,] 
dim(tab.q) #422744obs 24var
length(unique(tab.q$peptide)) #17449
length(unique(tab.q$protein)) #2510

tab.q$Rep <- gsub("4","3",tab.q$Rep) #remplace les 4 pas un 3 pour faciliter la mod?lisation 
tab.q$stade_rep <- gsub("_4","_3",tab.q$stade_rep)


################## vérification des contrastes
options(contrasts=c("contr.sum", "contr.poly")) #on modifie les contrastes par défaut

options("contrasts") #on vérifie qu'on a bien les bons contrastes

################## estimer les abondances de prot
proteines=unique(tab.q$protein)

tab.peptiz=NULL
tab.rep=NULL
tab.effect.msrun=NULL
tab.residue=NULL
tab.intercept=NULL
tab.stade=NULL

for (i in 1:length(proteines)){
  tryCatch({ #tryCatch sert à attrapper les messages d'erreur sans que les calculs de la boucle de n'interrompent
    myprot=proteines[i]
    sub= tab.q[tab.q$protein==myprot,]
    sub$Stade=as.factor(sub$Stade)
    sub$Rep=as.factor(sub$Rep)
    npep=length(unique(sub$peptiz))
    nStade=length(unique(sub$Stade))
    nrep=length(unique(sub$Rep))
    
    model = lme(logQnorm ~ Stade+Rep+peptiz, random=~1|msrun, data=sub ) 
    
    fixe=fixef(model)
    N=length(fixe)
    
    estimation=unique(sub[,c("Stade", "protein")])
    estimation=estimation[order(estimation$Stade),]
    intercept=rep(fixe[1], nrow(estimation))
    estimation$intercept=intercept
    tab.intercept=rbind.data.frame(tab.intercept, estimation)
    
    tab.stade=rbind.data.frame(tab.stade, cbind.data.frame(Stade.effect=c(fixe[2:nStade], -sum(fixe[2:nStade])), Stade=estimation$Stade, protein=rep(myprot, nStade)))
    
    tab.rep=rbind.data.frame(tab.rep, cbind.data.frame(rep.effect=c(fixe[(nStade+1):(nStade+2)], -sum(fixe[(nStade+1):(nStade+2)]) ), Rep=c("1", "2", "3"), protein=rep(myprot, nrep)))
    
    tab.peptiz=rbind.data.frame(tab.peptiz, cbind.data.frame(peptidez.effect=c(fixe[(nStade+nrep):(length(fixe))], -sum(fixe[(nStade+nrep):(length(fixe))])), peptiz=levels(as.factor(as.character(sub$peptiz))), protein=myprot))
    
    effect.msrun=model$coefficients$random$msrun
    tab.effect.msrun=rbind.data.frame(tab.effect.msrun, cbind.data.frame(msrun=rownames(effect.msrun), effet.msrun=effect.msrun[,1], protein=rep(myprot, nrow(effect.msrun))))
    
    tab.residue=rbind.data.frame(tab.residue, cbind.data.frame(msrun=sub$msrun, protein=sub$protein, peptiz=sub$peptiz, residue=resid(model)))
  },error=function(e){cat("Error for protein", myprot, ":", conditionMessage(e), "\n")}) #fin de tryCatch
}

tab.q=merge(tab.q, tab.intercept, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.stade, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.rep, by=c("Rep", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.peptiz, by=c("peptiz", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.effect.msrun, by=c("msrun", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.residue, by=c("msrun", "protein", "peptiz"), all.x=TRUE)
#tab.q$ab.estime=apply(tab.q[,c("intercept", "Stade.effect", "rep.effect")], 1, sum) ligne ajout?e qui cr?e des probl?mes

tab.sans.sdRT=tab.q
# load("./effet_filtre/data_isma_propre.RData")
# save(tab.sans.pep.commun,tab.sans.sdRT, file="./data_kiwi_Juan.RData")


# Filtre sur peptides repetees  -------------------------------------------
setwd("/media/juanma/JUANMA/Stage M2/git_repo/Proteo/")

# lecture du fichier de sortie quanti de masschroq
tab.q=read.table("peptides_q1_fractiona1.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#v?rification de la structure et du nb de peptides
str(tab.q) 
length(unique(tab.q$peptide)) 
# cr?ation de la colonne peptide-charge
tab.q$peptiz=paste(tab.q$peptide,tab.q$z,sep='-') 

tab.prot=read.table("fractiona1_proteins.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# v?rification de la structure (nb de lignes) et du nb de peptides
str(tab.prot) #23025 obs 3var
length(unique(tab.prot$peptide)) #20041 length(unique(tab.prot$peptide))=length(unique(tab.q$peptide)) on peut continuer

# homog?n?isation de tab.q et tab.prot
tab.prot=tab.prot[tab.prot$peptide %in% unique(tab.q$peptide),]		
length(unique(tab.prot$peptide)) #20041, comme dans tab.q	

msrunfile=unique(tab.q$msrun)
msrunfile=msrunfile[order(msrunfile)]

Stade=rev(c("bulk","bulk","bulk","bulk",rep(c("0 DPA","13 DPA","26 DPA","39 DPA","55 DPA","76 DPA","118 DPA","179 DPA","222 DPA"),each=3)))
Rep=rev(c(1,2,3,4,rep(c(1,2,3),9)))
stade_rep=  rev(c("bulk_1","bulk_2","bulk_3","bulk_4","0_1","0_2","0_3","13_1","13_2","13_3","26_1","26_2","26_3","39_1","39_2","39_3",
                  "55_1","55_2","55_3","76_1","76_2","76_3","118_1","118_2","118_3","179_1","179_2","179_3","222_1","222_2","222_3"))

indice=match(tab.q$msrun, msrunfile)
tab.q$Stade=Stade[indice]
tab.q$Rep=Rep[indice]
tab.q$stade_rep=stade_rep[indice]

# calcul de l'?cart-type du temps de r?tention 
total_rt <-tapply(tab.q$rt,list(tab.q$peptide),FUN=sd, na.rm=TRUE)
total_rt=total_rt[-which(is.na(total_rt))] # quand il n'y a qu'1 valeur par peptide, il ne peut pas calculer d'?cart-type et renvoie NA
#repr?sentation de la distribution des valeurs
hist(total_rt,nclass=500,freq=T,xlab="Deviation standart en sec", xlim=c(0,100))
# suppression des peptides pr?sentant trop de variabilit?
bad_peptides = names(total_rt[total_rt>20]) #213
tab.q=tab.q[!tab.q$peptide %in% bad_peptides,]
dim(tab.q) #549754     20
length(unique(tab.q$peptide)) #19825


################## calcul du biais=ref-obj
sample=unique(tab.q$msrun) #29
ref.sample="msruna1"
nbrun=length(sample) #29


################## Lissage du biais
mat.obj=tapply(tab.q$logarea,list(tab.q$peptiz,tab.q$msrun),FUN=sum)
ref=mat.obj[,ref.sample]
mat.ref=as.data.frame(matrix(rep(ref,nbrun),ncol=nbrun)) 
dim(mat.ref) #23383obs    29var

tab.biais=mat.obj-mat.ref
colnames(tab.biais) <- colnames(mat.obj)
tab.biais=stack(tab.biais)
biais.peptiz <- rep(rownames(mat.obj),nbrun)
tab.biais=cbind.data.frame(biais.peptiz, tab.biais)
colnames(tab.biais)=c("peptiz","biais","msrun")
dim(tab.biais)#678107 x 3
length(unique(tab.biais$peptiz)) #23383

indice=match(paste(tab.biais$msrun, tab.biais$peptiz), paste(tab.q$msrun, tab.q$peptiz))
tab.biais$rt=tab.q$rt[indice]

# on ordonne les lignes de tab.biais par msrun et par rt

tab.biais=tab.biais[order(tab.biais$msrun, tab.biais$rt),]
tab.biais=tab.biais[!is.na(tab.biais$rt),]#v?rification :mm nb de lignes que tab.q
tab.biais$correc=NA
compteur=unique(tab.biais$msrun)

#pdf("./figures/Normalisation_biais.pdf")
for (i in 1:length(compteur)) {
  #calcul de la correction
  select.msrun= tab.biais$msrun==compteur[i]	
  select.na <- is.na(tab.biais$biais)
  select.inf <- is.infinite(tab.biais$biais)
  
  cour <- tab.biais[select.msrun,]
  cour.nomiss <- tab.biais[select.msrun & !select.na & !select.inf,]
  
  # fonction de lissage de la courbe biais en fonction du rt
  yyy=smooth.spline(cour.nomiss$rt,cour.nomiss$biais,spar=0.5)
  # on ajoute le biais liss? (qui est le facteur de correction pour la normalisation) ? tab.biais
  indice=match(tab.biais$rt[select.msrun & !select.na], yyy$x)
  tab.biais$correc[select.msrun & !select.na] = yyy$y[indice]
  
  # repr?sentation graphique du biais et de la courbe de lissage				
  rt.min=min(cour.nomiss$rt)
  rt.max=max(cour.nomiss$rt)
  biais.min=min(cour.nomiss$biais)
  biais.max=max(cour.nomiss$biais)
  plot(cour.nomiss$rt,cour.nomiss$biais,main=paste("msrun = ", compteur[i], sep=""),xlab="RT",ylab="Biais (log)",xlim=c(rt.min,rt.max),ylim=c(biais.min,biais.max))
  lines(yyy$x,yyy$y,col="red",type="l")
  
  
  #calcul de la correction pour les peptides pr?sents dans le msrun ? normaliser et absents du msrun de r?f?rence: on prend la valeur de correction du peptide dont le temps de r?tention est juste avant
  smo=tab.biais$correc[select.msrun]
  j=1
  while(is.na(smo[j])) {smo[j]=0; j=j+1}
  while (j < length(smo)) {
    while(!is.na(smo)[j] & j<length(smo)) {last=smo[j]; j=j+1}
    while(is.na(smo[j]) & j<length(smo)) {smo[j]=last; j=j+1}
  }  
  smo[length(smo)]=smo[length(smo)-1] 
  tab.biais$correc[select.msrun] =smo
  
}
dev.off()

################## reunir l'info logQ et correc et tab.q
indice=match(paste(tab.q$msrun, tab.q$peptiz), paste(tab.biais$msrun, tab.biais$peptiz))
tab.q$correc=tab.biais$correc[indice]
head(tab.q)
dim(tab.q) # 549754 on v?rifie qu'on a bien le nb de ligne attendu (comme le tab.biais)

# on corrige les intensit?s de peptiz avec le facteur de correction	
tab.q$logQnorm=tab.q$logarea+tab.q$correc		



#-------------------- Suppression des peptides communs ? plusieurs prot?ines, des peptides peu r?p?tables 

#Retire les donn?es des bulks qui ne sont plus n?cessaires ? l'analyse
tab.q=tab.q[-which(tab.q$msrun=="msrunb31"),] #retire bulk1
tab.q=tab.q[-which(tab.q$msrun=="msrunb30"),] #retire bulk2
tab.q=tab.q[-which(tab.q$msrun=="msrunb29"),] #retire bulk3
tab.q=tab.q[-which(tab.q$msrun=="msrunb28"),] #retire bulk3


############# Suppression des peptides communs ? plusieurs proteines

#on compte le nb de prot auxquelles appartient chaque peptide
count = as.data.frame(table(tab.prot$peptide))
colnames(count)=c("peptide","Freq_prot")

# on ne conserve que les peptides qui n'appartiennent qu'? une seule prot?ine
good_peptides=count$peptide[count$Freq_prot==1] #17824 elements (attention, good_peptides est construit avec tab.prot qui n'a pas ?t? pr?alablement homog?n?is? avec tabok)
tab.q=tab.q[tab.q$peptide %in% good_peptides,] # 427255obs 22var

# on ajoute l'info prot?ine dans tab.q (sans risque de dupliquer des lignes puisque maintenant il n'y a plus de peptides communs)
tab.q = merge(tab.q,tab.prot,"peptide", all.x=TRUE, all.y=FALSE)
dim(tab.q) #427811obs 24var
length(unique(tab.q$peptide)) #17662
length(unique(tab.q$protein)) #2723

############# Suppression des peptides non répétables: un peptide doit être présent dans au moins 2 rep sur 3 et pour au moins 6 points de la gamme (pour pouvoir calculer des corrélations correctes) 

# pour chaque quantité d'UPS, on compte le nb de rep dans lesquels un peptiz est présent
tempFRIM=unique(tab.q[,c("peptiz", "msrun", "Rep", "Stade")])
countFRIM=table(tempFRIM$peptiz, tempFRIM$Stade)  #20640 x 9

# on sélectionne les peptiz qui ne sont présent dans moins de 2 rep dans au moins 1 échantillon
countFRIM2=countFRIM>1 # selection des pep présent dans au moins 2 rep partour
countFRIM2=apply(countFRIM2,1, sum)
bad.peptiz=names(countFRIM2[countFRIM2<9]) #213
countFRIM=countFRIM[!rownames(countFRIM) %in% bad.peptiz,] #10304

# on sélectionne les peptiz qui sont présent dans 3 rep dans au moins 28 échantillons => on tolère 15,15% de données manquantes
count2=countFRIM>2 # selection des pep présent dans au moins 2 rep partour
count2=apply(count2,1, sum)
bad.peptiz=names(count2[count2<0]) #213
countFRIM=countFRIM[!rownames(countFRIM) %in% bad.peptiz,] 	#10304

# on supprime les bad.peptiz de tabok
tab.q=tab.q[tab.q$peptiz %in% rownames(countFRIM),]
dim(tab.q) # 257066    24
length(unique(tab.q$peptide)) #9106
length(unique(tab.q$protein)) #2131
length(unique(tab.q$peptiz)) #10304

##################  Enlever les proteines avec moins de 2 peptides et pas des peptidez 
temp=unique(tab.q[,c("protein", "peptiz")])
count=table(temp$protein)
bad.prot=names(count[count<2]) #213 prot
tab.q=tab.q[!tab.q$protein %in% bad.prot,] 
dim(tab.q) #243817obs 24var
length(unique(tab.q$peptide)) #8562
length(unique(tab.q$protein)) #1587

tab.q$Rep <- gsub("4","3",tab.q$Rep) #remplace les 4 pas un 3 pour faciliter la mod?lisation 
tab.q$stade_rep <- gsub("_4","_3",tab.q$stade_rep)


################## vérification des contrastes
options(contrasts=c("contr.sum", "contr.poly")) #on modifie les contrastes par défaut

options("contrasts") #on vérifie qu'on a bien les bons contrastes

################## estimer les abondances de prot
proteines=unique(tab.q$protein)

tab.peptiz=NULL
tab.rep=NULL
tab.effect.msrun=NULL
tab.residue=NULL
tab.intercept=NULL
tab.stade=NULL

for (i in 1:length(proteines)){
  tryCatch({ #tryCatch sert à attrapper les messages d'erreur sans que les calculs de la boucle de n'interrompent
    myprot=proteines[i]
    sub= tab.q[tab.q$protein==myprot,]
    sub$Stade=as.factor(sub$Stade)
    sub$Rep=as.factor(sub$Rep)
    npep=length(unique(sub$peptiz))
    nStade=length(unique(sub$Stade))
    nrep=length(unique(sub$Rep))
    
    model = lme(logQnorm ~ Stade+Rep+peptiz, random=~1|msrun, data=sub ) 
    
    fixe=fixef(model)
    N=length(fixe)
    
    estimation=unique(sub[,c("Stade", "protein")])
    estimation=estimation[order(estimation$Stade),]
    intercept=rep(fixe[1], nrow(estimation))
    estimation$intercept=intercept
    tab.intercept=rbind.data.frame(tab.intercept, estimation)
    
    tab.stade=rbind.data.frame(tab.stade, cbind.data.frame(Stade.effect=c(fixe[2:nStade], -sum(fixe[2:nStade])), Stade=estimation$Stade, protein=rep(myprot, nStade)))
    
    tab.rep=rbind.data.frame(tab.rep, cbind.data.frame(rep.effect=c(fixe[(nStade+1):(nStade+2)], -sum(fixe[(nStade+1):(nStade+2)]) ), Rep=c("1", "2", "3"), protein=rep(myprot, nrep)))
    
    tab.peptiz=rbind.data.frame(tab.peptiz, cbind.data.frame(peptidez.effect=c(fixe[(nStade+nrep):(length(fixe))], -sum(fixe[(nStade+nrep):(length(fixe))])), peptiz=levels(as.factor(as.character(sub$peptiz))), protein=myprot))
    
    effect.msrun=model$coefficients$random$msrun
    tab.effect.msrun=rbind.data.frame(tab.effect.msrun, cbind.data.frame(msrun=rownames(effect.msrun), effet.msrun=effect.msrun[,1], protein=rep(myprot, nrow(effect.msrun))))
    
    tab.residue=rbind.data.frame(tab.residue, cbind.data.frame(msrun=sub$msrun, protein=sub$protein, peptiz=sub$peptiz, residue=resid(model)))
  },error=function(e){cat("Error for protein", myprot, ":", conditionMessage(e), "\n")}) #fin de tryCatch
}

# Error for protein c289.a2.a1 : NA/NaN/Inf in foreign function call (arg 1) 

tab.q=merge(tab.q, tab.intercept, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.stade, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.rep, by=c("Rep", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.peptiz, by=c("peptiz", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.effect.msrun, by=c("msrun", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.residue, by=c("msrun", "protein", "peptiz"), all.x=TRUE)
#tab.q$ab.estime=apply(tab.q[,c("intercept", "Stade.effect", "rep.effect")], 1, sum) ligne ajout?e qui cr?e des probl?mes



tab.sans.pepRep=tab.q
# Filtre sur peptides mal correlees ---------------------------------------

setwd("/media/juanma/JUANMA/Stage M2/git_repo/Proteo/")

# lecture du fichier de sortie quanti de masschroq
tab.q=read.table("peptides_q1_fractiona1.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#v?rification de la structure et du nb de peptides
str(tab.q) 
length(unique(tab.q$peptide)) 
# cr?ation de la colonne peptide-charge
tab.q$peptiz=paste(tab.q$peptide,tab.q$z,sep='-') 

tab.prot=read.table("fractiona1_proteins.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# v?rification de la structure (nb de lignes) et du nb de peptides
str(tab.prot) #23025 obs 3var
length(unique(tab.prot$peptide)) #20041 length(unique(tab.prot$peptide))=length(unique(tab.q$peptide)) on peut continuer

# homog?n?isation de tab.q et tab.prot
tab.prot=tab.prot[tab.prot$peptide %in% unique(tab.q$peptide),]		
length(unique(tab.prot$peptide)) #20041, comme dans tab.q	

msrunfile=unique(tab.q$msrun)
msrunfile=msrunfile[order(msrunfile)]

Stade=rev(c("bulk","bulk","bulk","bulk",rep(c("0 DPA","13 DPA","26 DPA","39 DPA","55 DPA","76 DPA","118 DPA","179 DPA","222 DPA"),each=3)))
Rep=rev(c(1,2,3,4,rep(c(1,2,3),9)))
stade_rep=  rev(c("bulk_1","bulk_2","bulk_3","bulk_4","0_1","0_2","0_3","13_1","13_2","13_3","26_1","26_2","26_3","39_1","39_2","39_3",
                  "55_1","55_2","55_3","76_1","76_2","76_3","118_1","118_2","118_3","179_1","179_2","179_3","222_1","222_2","222_3"))

indice=match(tab.q$msrun, msrunfile)
tab.q$Stade=Stade[indice]
tab.q$Rep=Rep[indice]
tab.q$stade_rep=stade_rep[indice]

# calcul de l'?cart-type du temps de r?tention 
total_rt <-tapply(tab.q$rt,list(tab.q$peptide),FUN=sd, na.rm=TRUE)
total_rt=total_rt[-which(is.na(total_rt))] # quand il n'y a qu'1 valeur par peptide, il ne peut pas calculer d'?cart-type et renvoie NA
#repr?sentation de la distribution des valeurs
hist(total_rt,nclass=500,freq=T,xlab="Deviation standart en sec", xlim=c(0,100))
# suppression des peptides pr?sentant trop de variabilit?
bad_peptides = names(total_rt[total_rt>20]) #213
tab.q=tab.q[!tab.q$peptide %in% bad_peptides,]
dim(tab.q) #549754     20
length(unique(tab.q$peptide)) #19825


################## calcul du biais=ref-obj
sample=unique(tab.q$msrun) #29
ref.sample="msruna1"
nbrun=length(sample) #29


################## Lissage du biais
mat.obj=tapply(tab.q$logarea,list(tab.q$peptiz,tab.q$msrun),FUN=sum)
ref=mat.obj[,ref.sample]
mat.ref=as.data.frame(matrix(rep(ref,nbrun),ncol=nbrun)) 
dim(mat.ref) #23383obs    29var

tab.biais=mat.obj-mat.ref
colnames(tab.biais) <- colnames(mat.obj)
tab.biais=stack(tab.biais)
biais.peptiz <- rep(rownames(mat.obj),nbrun)
tab.biais=cbind.data.frame(biais.peptiz, tab.biais)
colnames(tab.biais)=c("peptiz","biais","msrun")
dim(tab.biais)#678107 x 3
length(unique(tab.biais$peptiz)) #23383

indice=match(paste(tab.biais$msrun, tab.biais$peptiz), paste(tab.q$msrun, tab.q$peptiz))
tab.biais$rt=tab.q$rt[indice]

# on ordonne les lignes de tab.biais par msrun et par rt

tab.biais=tab.biais[order(tab.biais$msrun, tab.biais$rt),]
tab.biais=tab.biais[!is.na(tab.biais$rt),]#v?rification :mm nb de lignes que tab.q
tab.biais$correc=NA
compteur=unique(tab.biais$msrun)

#pdf("./figures/Normalisation_biais.pdf")
for (i in 1:length(compteur)) {
  #calcul de la correction
  select.msrun= tab.biais$msrun==compteur[i]	
  select.na <- is.na(tab.biais$biais)
  select.inf <- is.infinite(tab.biais$biais)
  
  cour <- tab.biais[select.msrun,]
  cour.nomiss <- tab.biais[select.msrun & !select.na & !select.inf,]
  
  # fonction de lissage de la courbe biais en fonction du rt
  yyy=smooth.spline(cour.nomiss$rt,cour.nomiss$biais,spar=0.5)
  # on ajoute le biais liss? (qui est le facteur de correction pour la normalisation) ? tab.biais
  indice=match(tab.biais$rt[select.msrun & !select.na], yyy$x)
  tab.biais$correc[select.msrun & !select.na] = yyy$y[indice]
  
  # repr?sentation graphique du biais et de la courbe de lissage				
  rt.min=min(cour.nomiss$rt)
  rt.max=max(cour.nomiss$rt)
  biais.min=min(cour.nomiss$biais)
  biais.max=max(cour.nomiss$biais)
  plot(cour.nomiss$rt,cour.nomiss$biais,main=paste("msrun = ", compteur[i], sep=""),xlab="RT",ylab="Biais (log)",xlim=c(rt.min,rt.max),ylim=c(biais.min,biais.max))
  lines(yyy$x,yyy$y,col="red",type="l")
  
  
  #calcul de la correction pour les peptides pr?sents dans le msrun ? normaliser et absents du msrun de r?f?rence: on prend la valeur de correction du peptide dont le temps de r?tention est juste avant
  smo=tab.biais$correc[select.msrun]
  j=1
  while(is.na(smo[j])) {smo[j]=0; j=j+1}
  while (j < length(smo)) {
    while(!is.na(smo)[j] & j<length(smo)) {last=smo[j]; j=j+1}
    while(is.na(smo[j]) & j<length(smo)) {smo[j]=last; j=j+1}
  }  
  smo[length(smo)]=smo[length(smo)-1] 
  tab.biais$correc[select.msrun] =smo
  
}
dev.off()

################## reunir l'info logQ et correc et tab.q
indice=match(paste(tab.q$msrun, tab.q$peptiz), paste(tab.biais$msrun, tab.biais$peptiz))
tab.q$correc=tab.biais$correc[indice]
head(tab.q)
dim(tab.q) # 549754 on v?rifie qu'on a bien le nb de ligne attendu (comme le tab.biais)

# on corrige les intensit?s de peptiz avec le facteur de correction	
tab.q$logQnorm=tab.q$logarea+tab.q$correc		



#-------------------- Suppression des peptides communs ? plusieurs prot?ines, des peptides peu r?p?tables et mal corr?l?s

#Retire les donn?es des bulks qui ne sont plus n?cessaires ? l'analyse
tab.q=tab.q[-which(tab.q$msrun=="msrunb31"),] #retire bulk1
tab.q=tab.q[-which(tab.q$msrun=="msrunb30"),] #retire bulk2
tab.q=tab.q[-which(tab.q$msrun=="msrunb29"),] #retire bulk3
tab.q=tab.q[-which(tab.q$msrun=="msrunb28"),] #retire bulk3

dim(tab.q)#490070     22var

############# Suppression des peptides communs ? plusieurs proteines

#on compte le nb de prot auxquelles appartient chaque peptide
count = as.data.frame(table(tab.prot$peptide))
colnames(count)=c("peptide","Freq_prot")

# on ne conserve que les peptides qui n'appartiennent qu'? une seule prot?ine
good_peptides=count$peptide[count$Freq_prot==1] #17824 elements (attention, good_peptides est construit avec tab.prot qui n'a pas ?t? pr?alablement homog?n?is? avec tabok)
tab.q=tab.q[tab.q$peptide %in% good_peptides,] # 427255obs 22var

# on ajoute l'info prot?ine dans tab.q (sans risque de dupliquer des lignes puisque maintenant il n'y a plus de peptides communs)
tab.q = merge(tab.q,tab.prot,"peptide", all.x=TRUE, all.y=FALSE)
dim(tab.q) #427811obs 24var
length(unique(tab.q$peptide)) #17662
length(unique(tab.q$protein)) #2723

############# Suppression des peptides non répétables: un peptide doit être présent dans au moins 2 rep sur 3 et pour au moins 6 points de la gamme (pour pouvoir calculer des corrélations correctes) 

# pour chaque quantité d'UPS, on compte le nb de rep dans lesquels un peptiz est présent
temp=unique(tab.q[,c("peptiz", "msrun", "Rep", "Stade")])
count=table(temp$peptiz, temp$Stade)  #20640 x 9

# on sélectionne les peptiz qui ne sont présent dans moins de 2 rep dans au moins 1 échantillon
count2=count>1 # selection des pep présent dans au moins 2 rep partour
count2=apply(count2,1, sum)
bad.peptiz=names(count2[count2<9]) #213
count=count[!rownames(count) %in% bad.peptiz,] #10304

# on sélectionne les peptiz qui sont présent dans 3 rep dans au moins 28 échantillons => on tolère 15,15% de données manquantes
count2=count>2 # selection des pep présent dans au moins 2 rep partour
count2=apply(count2,1, sum)
bad.peptiz=names(count2[count2<0]) #213
count=count[!rownames(count) %in% bad.peptiz,] 	#10304

# on supprime les bad.peptiz de tabok
tab.q=tab.q[tab.q$peptiz %in% rownames(count),]
dim(tab.q) # 257066    24
length(unique(tab.q$peptide)) #9106
length(unique(tab.q$protein)) #2131
length(unique(tab.q$peptiz)) #10304




##################  Enlever les proteines avec moins de 2 peptides et pas des peptidez 
temp=unique(tab.q[,c("protein", "peptiz")])
count=table(temp$protein)
bad.prot=names(count[count<2]) #213 prot
tab.q=tab.q[!tab.q$protein %in% bad.prot,] 
dim(tab.q) #243817obs 24var
length(unique(tab.q$peptide)) #8562
length(unique(tab.q$protein)) #1587


################## Enlever les peptides qui ne sont pas bien correles. Attention, ici on élimine bcp de peptides, même avec des critères peu stringeant car comme il n'y a pas de variation de qté dans les protéines de levure, on s'attend à ce que la corrélation soit de 0.Donc à faire uniquement pour les protéines UPS (sinon on enlève quasi toutes les prot de levure)!!!
tab.q$Stade <- gsub("7.7","07.7",tab.q$Stade)

#calculer les correlations pour toutes les protéines
proteines = unique(tab.q$protein) #1587 proteines
pepoktot=NULL
rCutoff=0.8
pvalCutoff=0.01
for (i in 1:length(proteines)){
  pepok=NULL
  sub=tab.q[tab.q$protein==proteines[i],]
  toto=tapply(sub$logQnorm,list(sub$Stade,sub$peptiz),FUN=mean,na.rm=TRUE)
  
  cor.cal=rcorr(as.matrix(NaRV.omit(toto)), type="pearson")
  diag(cor.cal$r)=NA
  tab.r=stack(as.data.frame(cor.cal$r))
  colnames(tab.r)=c("r", "pep1")
  tab.r$pep1=as.character(tab.r$pep1)
  tab.r$pep2=rep(rownames(cor.cal$r), ncol(cor.cal$r))
  tab.p=stack(as.data.frame(cor.cal$P))
  tab.r$pval=tab.p[,1]
  tab.n=stack(as.data.frame(cor.cal$n))
  tab.r$pval=tab.p[,1]
  tab.r=na.omit(tab.r)
  
  # looking for ref.pep= peptides with the highest nb of coefficients of correlation superior or equal to the mean of the positives coefficients of correlation	
  rmean=mean(tab.r$r[tab.r$r>0])
  if (!is.na(rmean)){	
    best.cor=tab.r[tab.r$r>=rmean,]
    best.cor=na.omit(best.cor)
    count=table(best.cor$pep1)
    ref.pep.temp=names(count[count==max(count)])
    rmean.pep=tapply(best.cor$r, list(best.cor$pep1), mean)
    rmean.pep=rmean.pep[names(rmean.pep) %in% ref.pep.temp]
    ref.pep=names(rmean.pep[rmean.pep==max(rmean.pep)])[1]
    tab.ref.pep=tab.r[tab.r$pep1==ref.pep,]
    if (length(tab.ref.pep$pep2[tab.ref.pep$r>rCutoff & tab.ref.pep$pval<pvalCutoff])>0){
      pepok=c(ref.pep, tab.ref.pep$pep2[tab.ref.pep$r>rCutoff & tab.ref.pep$pval<pvalCutoff])
    }else{
      pepok=NULL
    }
  }
  pepoktot=c(pepoktot, pepok)
  
}

tab.q=tab.q[tab.q$peptiz %in% pepoktot,]
dim(tab.q) # 169072
length(unique(tab.q$peptide)) # 5955
length(unique(tab.q$peptiz)) #6721
length(unique(tab.q$protein)) #1254


tab.q$Rep <- gsub("4","3",tab.q$Rep) #remplace les 4 pas un 3 pour faciliter la mod?lisation 
tab.q$stade_rep <- gsub("_4","_3",tab.q$stade_rep)

#-------------------- Calcul de l'abondances des protéines par la méthode modélisation


################## vérification des contrastes
options(contrasts=c("contr.sum", "contr.poly")) #on modifie les contrastes par défaut

options("contrasts") #on vérifie qu'on a bien les bons contrastes

################## estimer les abondances de prot
proteines=unique(tab.q$protein)

tab.peptiz=NULL
tab.rep=NULL
tab.effect.msrun=NULL
tab.residue=NULL
tab.intercept=NULL
tab.stade=NULL

for (i in 1:length(proteines)){
  tryCatch({ #tryCatch sert à attrapper les messages d'erreur sans que les calculs de la boucle de n'interrompent
    myprot=proteines[i]
    sub= tab.q[tab.q$protein==myprot,]
    sub$Stade=as.factor(sub$Stade)
    sub$Rep=as.factor(sub$Rep)
    npep=length(unique(sub$peptiz))
    nStade=length(unique(sub$Stade))
    nrep=length(unique(sub$Rep))
    
    model = lme(logQnorm ~ Stade+Rep+peptiz, random=~1|msrun, data=sub ) 
    
    fixe=fixef(model)
    N=length(fixe)
    
    estimation=unique(sub[,c("Stade", "protein")])
    estimation=estimation[order(estimation$Stade),]
    intercept=rep(fixe[1], nrow(estimation))
    estimation$intercept=intercept
    tab.intercept=rbind.data.frame(tab.intercept, estimation)
    
    tab.stade=rbind.data.frame(tab.stade, cbind.data.frame(Stade.effect=c(fixe[2:nStade], -sum(fixe[2:nStade])), Stade=estimation$Stade, protein=rep(myprot, nStade)))
    
    tab.rep=rbind.data.frame(tab.rep, cbind.data.frame(rep.effect=c(fixe[(nStade+1):(nStade+2)], -sum(fixe[(nStade+1):(nStade+2)]) ), Rep=c("1", "2", "3"), protein=rep(myprot, nrep)))
    
    tab.peptiz=rbind.data.frame(tab.peptiz, cbind.data.frame(peptidez.effect=c(fixe[(nStade+nrep):(length(fixe))], -sum(fixe[(nStade+nrep):(length(fixe))])), peptiz=levels(as.factor(as.character(sub$peptiz))), protein=myprot))
    
    effect.msrun=model$coefficients$random$msrun
    tab.effect.msrun=rbind.data.frame(tab.effect.msrun, cbind.data.frame(msrun=rownames(effect.msrun), effet.msrun=effect.msrun[,1], protein=rep(myprot, nrow(effect.msrun))))
    
    tab.residue=rbind.data.frame(tab.residue, cbind.data.frame(msrun=sub$msrun, protein=sub$protein, peptiz=sub$peptiz, residue=resid(model)))
  },error=function(e){cat("Error for protein", myprot, ":", conditionMessage(e), "\n")}) #fin de tryCatch
}

# Error for protein c289.a2.a1 : NA/NaN/Inf in foreign function call (arg 1) 

tab.q=merge(tab.q, tab.intercept, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.stade, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.rep, by=c("Rep", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.peptiz, by=c("peptiz", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.effect.msrun, by=c("msrun", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.residue, by=c("msrun", "protein", "peptiz"), all.x=TRUE)
#tab.q$ab.estime=apply(tab.q[,c("intercept", "Stade.effect", "rep.effect")], 1, sum) ligne ajout?e qui cr?e des probl?mes

tab.sans.pepCor=tab.q
# save(tab.sans.pep.commun,tab.sans.sdRT,tab.sans.pepRep,tab.sans.pepCor, file="./data_Proteo_Kiwi_Juan.RData")
# Aucun filtre ------------------------------------------------------------

setwd("/media/juanma/JUANMA/Stage M2/git_repo/Proteo/")

# lecture du fichier de sortie quanti de masschroq
tab.q=read.table("peptides_q1_fractiona1.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#v?rification de la structure et du nb de peptides
str(tab.q) 
length(unique(tab.q$peptide)) 
# cr?ation de la colonne peptide-charge
tab.q$peptiz=paste(tab.q$peptide,tab.q$z,sep='-') 

tab.prot=read.table("fractiona1_proteins.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
# v?rification de la structure (nb de lignes) et du nb de peptides
str(tab.prot) #23025 obs 3var
length(unique(tab.prot$peptide)) #20041 length(unique(tab.prot$peptide))=length(unique(tab.q$peptide)) on peut continuer

# homog?n?isation de tab.q et tab.prot
tab.prot=tab.prot[tab.prot$peptide %in% unique(tab.q$peptide),]		
length(unique(tab.prot$peptide)) #20041, comme dans tab.q	

msrunfile=unique(tab.q$msrun)
msrunfile=msrunfile[order(msrunfile)]

Stade=rev(c("bulk","bulk","bulk","bulk",rep(c("0 DPA","13 DPA","26 DPA","39 DPA","55 DPA","76 DPA","118 DPA","179 DPA","222 DPA"),each=3)))
Rep=rev(c(1,2,3,4,rep(c(1,2,3),9)))
stade_rep=  rev(c("bulk_1","bulk_2","bulk_3","bulk_4","0_1","0_2","0_3","13_1","13_2","13_3","26_1","26_2","26_3","39_1","39_2","39_3",
                  "55_1","55_2","55_3","76_1","76_2","76_3","118_1","118_2","118_3","179_1","179_2","179_3","222_1","222_2","222_3"))

indice=match(tab.q$msrun, msrunfile)
tab.q$Stade=Stade[indice]
tab.q$Rep=Rep[indice]
tab.q$stade_rep=stade_rep[indice]

################## calcul du biais=ref-obj
sample=unique(tab.q$msrun) #29
ref.sample="msruna1"
nbrun=length(sample) #29


################## Lissage du biais
mat.obj=tapply(tab.q$logarea,list(tab.q$peptiz,tab.q$msrun),FUN=sum)
ref=mat.obj[,ref.sample]
mat.ref=as.data.frame(matrix(rep(ref,nbrun),ncol=nbrun)) 
dim(mat.ref) #23383obs    29var

tab.biais=mat.obj-mat.ref
colnames(tab.biais) <- colnames(mat.obj)
tab.biais=stack(tab.biais)
biais.peptiz <- rep(rownames(mat.obj),nbrun)
tab.biais=cbind.data.frame(biais.peptiz, tab.biais)
colnames(tab.biais)=c("peptiz","biais","msrun")
dim(tab.biais)#678107 x 3
length(unique(tab.biais$peptiz)) #23383

indice=match(paste(tab.biais$msrun, tab.biais$peptiz), paste(tab.q$msrun, tab.q$peptiz))
tab.biais$rt=tab.q$rt[indice]

# on ordonne les lignes de tab.biais par msrun et par rt

tab.biais=tab.biais[order(tab.biais$msrun, tab.biais$rt),]
tab.biais=tab.biais[!is.na(tab.biais$rt),]#v?rification :mm nb de lignes que tab.q
tab.biais$correc=NA
compteur=unique(tab.biais$msrun)

#pdf("./figures/Normalisation_biais.pdf")
for (i in 1:length(compteur)) {
  #calcul de la correction
  select.msrun= tab.biais$msrun==compteur[i]	
  select.na <- is.na(tab.biais$biais)
  select.inf <- is.infinite(tab.biais$biais)
  
  cour <- tab.biais[select.msrun,]
  cour.nomiss <- tab.biais[select.msrun & !select.na & !select.inf,]
  
  # fonction de lissage de la courbe biais en fonction du rt
  yyy=smooth.spline(cour.nomiss$rt,cour.nomiss$biais,spar=0.5)
  # on ajoute le biais liss? (qui est le facteur de correction pour la normalisation) ? tab.biais
  indice=match(tab.biais$rt[select.msrun & !select.na], yyy$x)
  tab.biais$correc[select.msrun & !select.na] = yyy$y[indice]
  
  # repr?sentation graphique du biais et de la courbe de lissage				
  rt.min=min(cour.nomiss$rt)
  rt.max=max(cour.nomiss$rt)
  biais.min=min(cour.nomiss$biais)
  biais.max=max(cour.nomiss$biais)
  plot(cour.nomiss$rt,cour.nomiss$biais,main=paste("msrun = ", compteur[i], sep=""),xlab="RT",ylab="Biais (log)",xlim=c(rt.min,rt.max),ylim=c(biais.min,biais.max))
  lines(yyy$x,yyy$y,col="red",type="l")
  
  
  #calcul de la correction pour les peptides pr?sents dans le msrun ? normaliser et absents du msrun de r?f?rence: on prend la valeur de correction du peptide dont le temps de r?tention est juste avant
  smo=tab.biais$correc[select.msrun]
  j=1
  while(is.na(smo[j])) {smo[j]=0; j=j+1}
  while (j < length(smo)) {
    while(!is.na(smo)[j] & j<length(smo)) {last=smo[j]; j=j+1}
    while(is.na(smo[j]) & j<length(smo)) {smo[j]=last; j=j+1}
  }  
  smo[length(smo)]=smo[length(smo)-1] 
  tab.biais$correc[select.msrun] =smo
  
}
dev.off()

################## reunir l'info logQ et correc et tab.q
indice=match(paste(tab.q$msrun, tab.q$peptiz), paste(tab.biais$msrun, tab.biais$peptiz))
tab.q$correc=tab.biais$correc[indice]
head(tab.q)
dim(tab.q) # 549754 on v?rifie qu'on a bien le nb de ligne attendu (comme le tab.biais)

# on corrige les intensit?s de peptiz avec le facteur de correction	
tab.q$logQnorm=tab.q$logarea+tab.q$correc		



#-------------------- Suppression des peptides communs ? plusieurs prot?ines, des peptides peu r?p?tables et mal corr?l?s

#Retire les donn?es des bulks qui ne sont plus n?cessaires ? l'analyse
tab.q=tab.q[-which(tab.q$msrun=="msrunb31"),] #retire bulk1
tab.q=tab.q[-which(tab.q$msrun=="msrunb30"),] #retire bulk2
tab.q=tab.q[-which(tab.q$msrun=="msrunb29"),] #retire bulk3
tab.q=tab.q[-which(tab.q$msrun=="msrunb28"),] #retire bulk3

dim(tab.q)#490070     22var

# on ajoute l'info prot?ine dans tab.q (sans risque de dupliquer des lignes puisque maintenant il n'y a plus de peptides communs)
tab.q = merge(tab.q,tab.prot,"peptide", all.x=TRUE, all.y=FALSE)
dim(tab.q) #427811obs 24var
length(unique(tab.q$peptide)) #17662
length(unique(tab.q$protein)) #2723

##################  Enlever les proteines avec moins de 2 peptides et pas des peptidez 
temp=unique(tab.q[,c("protein", "peptide")])
count=table(temp$protein)
bad.prot=names(count[count<2]) 
tab.q=tab.q[!tab.q$protein %in% bad.prot,] 
dim(tab.q) #584250 24var
length(unique(tab.q$peptide)) #20034
length(unique(tab.q$protein)) #2727


#-------------------- Calcul de l'abondances des protéines par la méthode modélisation

tab.q$Rep <- gsub("4","3",tab.q$Rep) #remplace les 4 pas un 3 pour faciliter la mod?lisation 
tab.q$stade_rep <- gsub("_4","_3",tab.q$stade_rep)

################## vérification des contrastes
options(contrasts=c("contr.sum", "contr.poly")) #on modifie les contrastes par défaut

options("contrasts") #on vérifie qu'on a bien les bons contrastes

################## estimer les abondances de prot
proteines=unique(tab.q$protein)

tab.peptiz=NULL
tab.rep=NULL
tab.effect.msrun=NULL
tab.residue=NULL
tab.intercept=NULL
tab.stade=NULL

for (i in 1:length(proteines)){
  tryCatch({ #tryCatch sert à attrapper les messages d'erreur sans que les calculs de la boucle de n'interrompent
    myprot=proteines[i]
    sub= tab.q[tab.q$protein==myprot,]
    sub$Stade=as.factor(sub$Stade)
    sub$Rep=as.factor(sub$Rep)
    npep=length(unique(sub$peptiz))
    nStade=length(unique(sub$Stade))
    nrep=length(unique(sub$Rep))
    
    model = lme(logQnorm ~ Stade+Rep+peptiz, random=~1|msrun, data=sub ) 
    
    fixe=fixef(model)
    N=length(fixe)
    
    estimation=unique(sub[,c("Stade", "protein")])
    estimation=estimation[order(estimation$Stade),]
    intercept=rep(fixe[1], nrow(estimation))
    estimation$intercept=intercept
    tab.intercept=rbind.data.frame(tab.intercept, estimation)
    
    tab.stade=rbind.data.frame(tab.stade, cbind.data.frame(Stade.effect=c(fixe[2:nStade], -sum(fixe[2:nStade])), Stade=estimation$Stade, protein=rep(myprot, nStade)))
    
    tab.rep=rbind.data.frame(tab.rep, cbind.data.frame(rep.effect=c(fixe[(nStade+1):(nStade+2)], -sum(fixe[(nStade+1):(nStade+2)]) ), Rep=c("1", "2", "3"), protein=rep(myprot, nrep)))
    
    tab.peptiz=rbind.data.frame(tab.peptiz, cbind.data.frame(peptidez.effect=c(fixe[(nStade+nrep):(length(fixe))], -sum(fixe[(nStade+nrep):(length(fixe))])), peptiz=levels(as.factor(as.character(sub$peptiz))), protein=myprot))
    
    effect.msrun=model$coefficients$random$msrun
    tab.effect.msrun=rbind.data.frame(tab.effect.msrun, cbind.data.frame(msrun=rownames(effect.msrun), effet.msrun=effect.msrun[,1], protein=rep(myprot, nrow(effect.msrun))))
    
    tab.residue=rbind.data.frame(tab.residue, cbind.data.frame(msrun=sub$msrun, protein=sub$protein, peptiz=sub$peptiz, residue=resid(model)))
  },error=function(e){cat("Error for protein", myprot, ":", conditionMessage(e), "\n")}) #fin de tryCatch
}

# Error for protein d2114.a2.a1 : contrasts can be applied only to factors with 2 or more levels 
# Error for protein d2571.a1.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3316.a2.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3349.a1.a1 : the leading minor of order 1 is not positive definite 
# Error for protein d3376.a1.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3407.a1.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3417.a1.a1 : Singularity in backsolve at level 0, block 1 
# Error for protein d3420.a1.a1 : Singularity in backsolve at level 0, block 1 

tab.q=merge(tab.q, tab.intercept, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.stade, by=c("Stade", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.rep, by=c("Rep", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.peptiz, by=c("peptiz", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.effect.msrun, by=c("msrun", "protein"), all.x=TRUE)
tab.q=merge(tab.q, tab.residue, by=c("msrun", "protein", "peptiz"), all.x=TRUE)
#tab.q$ab.estime=apply(tab.q[,c("intercept", "Stade.effect", "rep.effect")], 1, sum) ligne ajout?e qui cr?e des probl?mes

tab.sans.filtre=tab.q
load("./effet_filtre/data_isma_propre.RData")

tab.sans.filtre=tab.sans.filtre[,!colnames(tab.sans.filtre) %in% c("maxintensity", "rtbegin", "rtend", "isotope", "mods", "ninumber","nirank", "niratio", "msrunfile", "group", "mz", "sequence", "z")]
tab.sans.pep.commun=tab.sans.pep.commun[,!colnames(tab.sans.pep.commun) %in% c("maxintensity", "rtbegin", "rtend", "isotope", "mods", "ninumber","nirank", "niratio", "msrunfile", "group", "mz", "sequence", "z", "protein_description")]
tab.sans.sdRT=tab.sans.sdRT[,!colnames(tab.sans.sdRT) %in% c("maxintensity", "rtbegin", "rtend", "isotope", "mods", "ninumber","nirank", "niratio", "msrunfile", "group", "mz", "sequence", "z", "protein_description")]
tab.sans.pepRep=tab.sans.pepRep[,!colnames(tab.sans.pepRep) %in% c("maxintensity", "rtbegin", "rtend", "isotope", "mods", "ninumber","nirank", "niratio", "msrunfile", "group", "mz", "sequence", "z", "protein_description")]
tab.sans.pepCor=tab.sans.pepCor[,!colnames(tab.sans.pepCor) %in% c("maxintensity", "rtbegin", "rtend", "isotope", "mods", "ninumber","nirank", "niratio", "msrunfile", "group", "mz", "sequence", "z", "protein_description")]

save(tab.sans.filtre,tab.sans.pep.commun,tab.sans.sdRT,tab.sans.pepRep,tab.sans.pepCor,file="./data_kiwi_juan.RData")


length(unique(tab.sans.filtre$peptiz))#prot=2727 / peptide=20034 / peptiz= 23660
length(unique(tab.sans.pep.commun$peptiz))#prot=2517 / peptide=17612 / peptiz= 20610
length(unique(tab.sans.sdRT$peptiz))#prot=2510 / peptide=17449 / peptiz= 20400
length(unique(tab.sans.pepRep$peptiz))#prot=1587 / peptide=8562 / peptiz= 9760
length(unique(tab.sans.pepCor$peptiz))#prot=1254 / peptide=5955 / peptiz= 6721



# Ibaq calculs ------------------------------------------------------------
require("readxl")
require("readODS")
setwd("/media/juanma/JUANMA/Stage M2/Analyse_stats/Proteo/")
# setwd("D:/Stage M2/Analyse_stats/Proteo/")

rm(list=ls())
load("./data_kiwi_Juan.RData")
prot_data=read.csv("20191217_Baldet_72_kiwi_nsaf_good_sample_proteins.csv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
prot_data$MW<-as.numeric(sub(",",".",prot_data$MW,fixed=T))
prot_data2<-read_ods("20191217_Baldet_72_kiwi_clean_sample.ods",sheet = "proteins") ## Nouvelles données corrigées?? 
colnames(prot_data2)<-colnames(prot_data)
prot_data2$MW<-as.numeric(sub(",",".",prot_data2$MW,fixed=T))
# prot_data<-prot_data[!grepl("UPS",prot_data$accession),] ## Enlever proteines UPS
quant_data<-read_excel("Protein_total.xlsx")
info_data<-read_ods("metadata.ods")
info_data<-info_data[-which(info_data$stade=="bulk"),]
info_data["Protein (mg/gFW)"]<-quant_data$`Protein (mg/gFW)`
# conv_factor <- read.table("./data_brutes/facteur_multiplicateur.txt",header=T,sep="\t",dec=".")
conv_factor<-info_data
conv_factor$stade <- gsub("7.7","07.7",conv_factor$stade)
prot_data<-prot_data2 ## Choisir entre nouvelle ou ancienne table MW 


# A faire successivement sans oublierde changer le nom ? la fin
# tabok=tab.sans.filtre 
# rm(list=ls())
# 
# tabok=tab.sans.pep.commun
# rm(list=ls())
# 
# tabok=tab.sans.sdRT 
# rm(list=ls())
# 
# tabok=tab.sans.pepRep 
# rm(list=ls())
# 
# tabok=tab.sans.pepCor 

etape2<-function(tabok,prot_data,conv_factor){
  #### calcul de l'iBAQ
  # attention, la somme des intensités se fait sur les abondances non transformées par le log. Il faut donc que tu prennes 10^(tabok$logQnorm). Voilà ce qui est écrit dans le supplementary material de l'article de schwanhausser:
  #MaxQuant computes protein intensities as the sum of all identified peptide intensities (maximum detector peak intensities of the peptide elution profile, including all peaks in the isotope cluster). Protein intensities were divided by the number of theoretically observable peptides (calculated by in silico protein digestion with a PERL script, all fully tryptic peptides between 6 and 30 amino acids were counted while missed cleavages were neglected). The resulting “iBAQ” intensities were log-transformed and plotted against known log-transformed absolute molar amounts of the spiked-in standard proteins (UPS2 standard). Linear regression was used to fit iBAQ intensities to absolute standard protein amounts.
  
  ibaq=tapply(10^(tabok$logQnorm), list(tabok$protein, tabok$msrun), sum)
  length(rownames(ibaq))
  prot=rep(rownames(ibaq), ncol(ibaq))
  ibaq=stack(as.data.frame(ibaq))
  colnames(ibaq)=c("sum", "msrun")
  ibaq$protein=prot
  ibaq=merge(ibaq, unique(tabok[, c("msrun", "Rep", "Stade")]), "msrun")
  ibaq=merge(ibaq, unique(prot_data[,c("Protein.ID", "Theoretical.number.of.tryptic.peptides")]), by.x="protein", by.y="Protein.ID")
  ibaq$ibaq=ibaq$sum/ibaq$Theoretical.number.of.tryptic.peptides
  
  
  ibaq_fmolgFW = ibaq
  sumRatioIbaq = aggregate(ibaq_fmolgFW$ibaq, by=list(msrun=ibaq_fmolgFW$msrun), FUN=sum,na.rm=T)
  ibaq_fmolgFW <- merge(ibaq_fmolgFW,sumRatioIbaq,by="msrun")
  colnames(ibaq_fmolgFW)[ncol(ibaq_fmolgFW)] <- "sumRatio"
  ibaq_fmolgFW$Ratio = ibaq_fmolgFW$ibaq/ibaq_fmolgFW$sumRatio
  ibaq_fmolgFW <- merge(ibaq_fmolgFW,prot_data[,c("Protein.ID","MW")],by.x="protein",by.y="Protein.ID")
  ibaq_fmolgFW$MW <- 10^3*(ibaq_fmolgFW$MW)
  ibaq_fmolgFW$Ratio2 <- ibaq_fmolgFW$Ratio/ibaq_fmolgFW$MW
  ibaq_fmolgFW <- merge(ibaq_fmolgFW,conv_factor[,c("msrun","Protein (mg/gFW)")],by="msrun")
  ibaq_fmolgFW$conc <- ibaq_fmolgFW$Ratio2*ibaq_fmolgFW$`Protein (mg/gFW)`
  ibaq_fmolgFW$ibaq_fmolgFW <-  ibaq_fmolgFW$conc*10^15
  ibaq_fmolgFW <- ibaq_fmolgFW[,c("protein", "msrun","ibaq_fmolgFW")]
  ibaq_fmolgFW <- merge(ibaq_fmolgFW,conv_factor[,c("msrun","stade")],by="msrun")
  ibaq_fmolgFW <- merge(ibaq_fmolgFW,prot_data[,c("Protein.ID","accession")],by.x="protein",by.y="Protein.ID")
  
  #### calcul du Top3 
  
  #selection des 3 peptides les plus intenses
  med=as.data.frame(tapply(tabok$logQnorm, list(tabok$peptiz), median))
  colnames(med)="median"
  med$peptiz=rownames(med)
  med=merge(med, unique(tabok[,c("peptiz", "protein")]), "peptiz")
  
  prot=unique(med$protein)
  pep.int=NULL
  for (i in 1:length(prot)){
    sub=med[med$protein==prot[i],]
    if (nrow(sub)>2){
      sub=sub[order(sub$median, decreasing=TRUE),]
      pep.int=c(pep.int, sub$peptiz[1:3])
    }
  }
  
  top3=tabok[tabok$peptiz %in% pep.int,]
  dim(top3) #142546
  length(unique(top3$protein)) #2162
  
  TOP3 <-tapply(10^(top3$logQnorm), list(top3$protein, top3$msrun), mean, na.rm=TRUE)
  prot=rep(rownames(TOP3), ncol(TOP3))
  TOP3=stack(as.data.frame(TOP3))
  colnames(TOP3)=c("top3", "msrun")
  TOP3$protein=prot
  tab.index=merge(TOP3, ibaq, c("protein", "msrun"), all.x=TRUE, all.y=TRUE)
  
  TOP3_fmolgFW = TOP3
  sumRatioTOP3 = aggregate(TOP3_fmolgFW$top3, by=list(msrun=TOP3_fmolgFW$msrun), FUN=sum,na.rm=T)
  TOP3_fmolgFW <- merge(TOP3_fmolgFW,sumRatioTOP3,by="msrun")
  colnames(TOP3_fmolgFW)[ncol(TOP3_fmolgFW)] <- "sumRatio"
  TOP3_fmolgFW$Ratio = TOP3_fmolgFW$top3/TOP3_fmolgFW$sumRatio
  TOP3_fmolgFW <- merge(TOP3_fmolgFW,prot_data[,c("Protein.ID","MW")],by.x="protein",by.y="Protein.ID")
  TOP3_fmolgFW$MW <- 10^3*(TOP3_fmolgFW$MW)
  TOP3_fmolgFW$Ratio2 <- TOP3_fmolgFW$Ratio/TOP3_fmolgFW$MW
  TOP3_fmolgFW <- merge(TOP3_fmolgFW,conv_factor[,c("msrun","Protein (mg/gFW)")],by="msrun")
  TOP3_fmolgFW$conc <- TOP3_fmolgFW$Ratio2*TOP3_fmolgFW$`Protein (mg/gFW)`
  TOP3_fmolgFW$TOP3_fmolgFW <-  TOP3_fmolgFW$conc*10^15
  TOP3_fmolgFW <- TOP3_fmolgFW[,c("protein", "msrun","TOP3_fmolgFW")]
  TOP3_fmolgFW <- merge(TOP3_fmolgFW,conv_factor[,c("msrun","stade")],by="msrun")
  TOP3_fmolgFW <- merge(TOP3_fmolgFW,prot_data[,c("Protein.ID","accession")],by.x="protein",by.y="Protein.ID")
  
  index.conc <- merge(TOP3_fmolgFW, ibaq_fmolgFW, c("protein", "msrun","stade","accession"), all.x=TRUE, all.y=TRUE)
  
  
  #### calcul de average sur non log
  ave <-tapply(10^(tabok$logQnorm), list(tabok$protein, tabok$msrun), mean, na.rm=TRUE)
  prot=rep(rownames(ave), ncol(ave))
  ave=stack(as.data.frame(ave))
  colnames(ave)=c("average", "msrun")
  ave$protein=prot
  tab.index=merge(tab.index, ave, c("protein", "msrun"), all.x=TRUE, all.y=TRUE)
  
  ave_fmolgFW = ave
  sumRatioave = aggregate(ave_fmolgFW$average, by=list(msrun=ave_fmolgFW$msrun), FUN=sum,na.rm=T)
  ave_fmolgFW <- merge(ave_fmolgFW,sumRatioave,by="msrun")
  colnames(ave_fmolgFW)[ncol(ave_fmolgFW)] <- "sumRatio"
  ave_fmolgFW$Ratio = ave_fmolgFW$average/ave_fmolgFW$sumRatio
  ave_fmolgFW <- merge(ave_fmolgFW,prot_data[,c("Protein.ID","MW")],by.x="protein",by.y="Protein.ID")
  ave_fmolgFW$MW <- 10^3*(ave_fmolgFW$MW)
  ave_fmolgFW$Ratio2 <- ave_fmolgFW$Ratio/ave_fmolgFW$MW
  ave_fmolgFW <- merge(ave_fmolgFW,conv_factor[,c("msrun","Protein (mg/gFW)")],by="msrun")
  ave_fmolgFW$conc <- ave_fmolgFW$Ratio2*ave_fmolgFW$`Protein (mg/gFW)`
  ave_fmolgFW$ave_fmolgFW <-  ave_fmolgFW$conc*10^15
  ave_fmolgFW <- ave_fmolgFW[,c("protein", "msrun","ave_fmolgFW")]
  ave_fmolgFW <- merge(ave_fmolgFW,conv_factor[,c("msrun","stade")],by="msrun")
  ave_fmolgFW <- merge(ave_fmolgFW,prot_data[,c("Protein.ID","accession")],by.x="protein",by.y="Protein.ID")
  
  index.conc <- merge(index.conc, ave_fmolgFW, c("protein", "msrun","stade","accession"), all.x=TRUE, all.y=TRUE)
  
  #### calcul de average sur log
  avelog <-tapply(tabok$logQnorm, list(tabok$protein, tabok$msrun), mean, na.rm=TRUE)
  prot=rep(rownames(avelog), ncol(avelog))
  avelog=stack(as.data.frame(avelog))
  colnames(avelog)=c("averageLog", "msrun")
  avelog$protein=prot
  tab.index=merge(tab.index, avelog, c("protein", "msrun"), all.x=TRUE, all.y=TRUE)
  
  avelog_fmolgFW = avelog
  avelog_fmolgFW$averageLog[!is.finite(avelog_fmolgFW$averageLog)] <- NA
  sumRatioavelog = aggregate(10^(avelog_fmolgFW$averageLog), by=list(msrun=avelog_fmolgFW$msrun), FUN=sum,na.rm=T)
  avelog_fmolgFW <- merge(avelog_fmolgFW,sumRatioavelog,by="msrun")
  colnames(avelog_fmolgFW)[ncol(avelog_fmolgFW)] <- "sumRatio"    
  avelog_fmolgFW$Ratio = 10^(avelog_fmolgFW$averageLog)/avelog_fmolgFW$sumRatio
  avelog_fmolgFW <- merge(avelog_fmolgFW,prot_data[,c("Protein.ID","MW")],by.x="protein",by.y="Protein.ID")
  avelog_fmolgFW$MW <- 10^3*(avelog_fmolgFW$MW)
  avelog_fmolgFW$Ratio2 <- avelog_fmolgFW$Ratio/avelog_fmolgFW$MW
  avelog_fmolgFW <- merge(avelog_fmolgFW,conv_factor[,c("msrun","Protein (mg/gFW)")],by="msrun")
  avelog_fmolgFW$conc <- avelog_fmolgFW$Ratio2*avelog_fmolgFW$`Protein (mg/gFW)`
  avelog_fmolgFW$avelog_fmolgFW <-  avelog_fmolgFW$conc*10^15
  avelog_fmolgFW <- avelog_fmolgFW[,c("protein", "msrun","avelog_fmolgFW")]
  avelog_fmolgFW <- merge(avelog_fmolgFW,conv_factor[,c("msrun","stade")],by="msrun")
  avelog_fmolgFW <- merge(avelog_fmolgFW,prot_data[,c("Protein.ID","accession")],by.x="protein",by.y="Protein.ID")
  
  index.conc <- merge(index.conc, avelog_fmolgFW, c("protein", "msrun","stade","accession"), all.x=TRUE, all.y=TRUE)
  
  
  #### calcul de modelisation sans rep
  tabok$ab.estime.res=tabok$logQnorm-tabok$peptidez.effect-tabok$effet.msrun-tabok$rep.effect
  model=tapply(tabok$ab.estime.res, list(tabok$protein, tabok$msrun), mean, na.rm=TRUE)
  prot=rep(rownames(model), ncol(model))
  model=stack(as.data.frame(model))
  colnames(model)=c("model_sans_rep", "msrun")
  model$protein=prot
  tab.index=merge(tab.index, model, c("protein", "msrun"), all.x=TRUE, all.y=TRUE)
  
  modelssRep_fmolgFW = model
  sumRatioModelssRep = aggregate(10^(modelssRep_fmolgFW$model_sans_rep), by=list(msrun=modelssRep_fmolgFW$msrun), FUN=sum,na.rm=T)
  modelssRep_fmolgFW <- merge(modelssRep_fmolgFW,sumRatioModelssRep,by="msrun")
  colnames(modelssRep_fmolgFW)[ncol(modelssRep_fmolgFW)] <- "sumRatio"    
  modelssRep_fmolgFW$Ratio = 10^(modelssRep_fmolgFW$model_sans_rep)/modelssRep_fmolgFW$sumRatio
  modelssRep_fmolgFW <- merge(modelssRep_fmolgFW,prot_data[,c("Protein.ID","MW")],by.x="protein",by.y="Protein.ID")
  modelssRep_fmolgFW$MW <- 10^3*(modelssRep_fmolgFW$MW)
  modelssRep_fmolgFW$Ratio2 <- modelssRep_fmolgFW$Ratio/modelssRep_fmolgFW$MW
  modelssRep_fmolgFW <- merge(modelssRep_fmolgFW,conv_factor[,c("msrun","Protein (mg/gFW)")],by="msrun")
  modelssRep_fmolgFW$conc <- modelssRep_fmolgFW$Ratio2*modelssRep_fmolgFW$`Protein (mg/gFW)`
  modelssRep_fmolgFW$modelssRep_fmolgFW <-  modelssRep_fmolgFW$conc*10^15
  modelssRep_fmolgFW <- modelssRep_fmolgFW[,c("protein", "msrun","modelssRep_fmolgFW")]
  modelssRep_fmolgFW <- merge(modelssRep_fmolgFW,conv_factor[,c("msrun","stade")],by="msrun")
  modelssRep_fmolgFW <- merge(modelssRep_fmolgFW,prot_data[,c("Protein.ID","accession")],by.x="protein",by.y="Protein.ID")
  
  index.conc <- merge(index.conc, modelssRep_fmolgFW, c("protein", "msrun","stade","accession"), all.x=TRUE, all.y=TRUE)
  
  
  #### calcul de modelisation avec rep
  tabok$ab.estime.res=tabok$logQnorm-tabok$peptidez.effect-tabok$effet.msrun
  model=tapply(tabok$ab.estime.res, list(tabok$protein, tabok$msrun), mean, na.rm=TRUE)
  prot=rep(rownames(model), ncol(model))
  model=stack(as.data.frame(model))
  colnames(model)=c("model_avec_rep", "msrun")
  model$protein=prot
  tab.index=merge(tab.index, model, c("protein", "msrun"), all.x=TRUE, all.y=TRUE)
  
  modelAcRep_fmolgFW = model
  sumRatioModelssRep = aggregate(10^(modelAcRep_fmolgFW$model_avec_rep), by=list(msrun=modelAcRep_fmolgFW$msrun), FUN=sum,na.rm=T)
  modelAcRep_fmolgFW <- merge(modelAcRep_fmolgFW,sumRatioModelssRep,by="msrun")
  colnames(modelAcRep_fmolgFW)[ncol(modelAcRep_fmolgFW)] <- "sumRatio"    
  modelAcRep_fmolgFW$Ratio = 10^(modelAcRep_fmolgFW$model_avec_rep)/modelAcRep_fmolgFW$sumRatio
  modelAcRep_fmolgFW <- merge(modelAcRep_fmolgFW,prot_data[,c("Protein.ID","MW")],by.x="protein",by.y="Protein.ID")
  modelAcRep_fmolgFW$MW <- 10^3*(modelAcRep_fmolgFW$MW)
  modelAcRep_fmolgFW$Ratio2 <- modelAcRep_fmolgFW$Ratio/modelAcRep_fmolgFW$MW
  modelAcRep_fmolgFW <- merge(modelAcRep_fmolgFW,conv_factor[,c("msrun","Protein (mg/gFW)")],by="msrun")
  modelAcRep_fmolgFW$conc <- modelAcRep_fmolgFW$Ratio2*modelAcRep_fmolgFW$`Protein (mg/gFW)`
  modelAcRep_fmolgFW$modelAcRep_fmolgFW <-  modelAcRep_fmolgFW$conc*10^15
  modelAcRep_fmolgFW <- modelAcRep_fmolgFW[,c("protein", "msrun","modelAcRep_fmolgFW")]
  modelAcRep_fmolgFW <- merge(modelAcRep_fmolgFW,conv_factor[,c("msrun","stade")],by="msrun")
  modelAcRep_fmolgFW <- merge(modelAcRep_fmolgFW,prot_data[,c("Protein.ID","accession")],by.x="protein",by.y="Protein.ID")
  
  index.conc <- merge(index.conc, modelAcRep_fmolgFW, c("protein", "msrun","stade","accession"), all.x=TRUE, all.y=TRUE)
  index.conc <- melt(index.conc,id=c("protein","msrun","stade","accession"))
  colnames(index.conc)[3] <- "index"
  return(list("tab"=tab.index,"conc"=index.conc))
}
noFil<-etape2(tab.sans.filtre,prot_data,conv_factor)
index.sans.filtre=noFil$tab
conc.sans.filtre=noFil$conc

noPep<-etape2(tab.sans.pep.commun,prot_data,conv_factor)
index.sans.pep.commun=noPep$tab
conc.sans.pep.commun=noPep$conc

# ups_id<-prot_data[grepl("UPS",prot_data$accession),"Protein.ID"]
# tab.sans.sdRT<-tab.sans.sdRT[-which(tab.sans.sdRT$protein %in% ups_id),]
filsdRT<-etape2(tab.sans.sdRT,prot_data,conv_factor)
index.sans.sdRT=filsdRT$tab
conc.sans.sdRT=filsdRT$conc

nopepRep<-etape2(tab.sans.pepRep,prot_data,conv_factor)
index.sans.pepRep=nopepRep$tab
conc.sans.pepRep=nopepRep$conc

nopepCor<-etape2(tab.sans.pepCor,prot_data,conv_factor)
index.sans.pepCor=nopepCor$tab
conc.sans.pepCor=nopepCor$conc

save(index.sans.filtre,index.sans.pep.commun,index.sans.sdRT,index.sans.pepRep ,index.sans.pepCor, file="./kiwi_juan_propre_correc_sansUPS_nouvellle_tableMW.RData")

save(conc.sans.filtre,conc.sans.pep.commun,conc.sans.sdRT,conc.sans.pepRep ,conc.sans.pepCor, file="./kiwi_juan_Conc(fmolgFW)_index_correc_sansUPS_nouvelle_tableMW.RData")

##############              Precision = CV(%) entre les replicas


load("./kiwi_juan_Conc(fmolgFW)_index.RData")
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

mylist=list(conc.sans.filtre,conc.sans.pep.commun,conc.sans.sdRT,conc.sans.pepRep ,conc.sans.pepCor)
type_quanti=c("sans.filtre", "sans.pep.commun", "sans.sdRT", "sans.pepRep", "sans.pepCor")
index= c("ibaq", "top3", "average", "averageLog", "model_sans_rep", "model_avec_rep")
conv_factor <- read.table("./data_brutes/facteur_multiplicateur.txt",header=T,sep="\t",dec=".")

stade=c("07.7","15","21.7","28","34.3","41.3","48.5","50.3","53.0")
index= c("ibaq_fmolgFW", "TOP3_fmolgFW", "ave_fmolgFW", "avelog_fmolgFW", "modelssRep_fmolgFW", "modelAcRep_fmolgFW")
index <- index[-5]
index_title <- c("iBAQ","TOP3","Average","Average-Log", "Model")

CVtot=NULL
for (elem in 1:5){
  tab.index=mylist[[elem]]
  colnames(tab.index) <- c("protein","msrun","stade","index","value")
  tab.index$amin=paste0(tab.index$index,"_", tab.index$stade)
  
  tab.sd=tapply(tab.index$value, list(tab.index$protein, tab.index$amin), sd)
  tab.moy=tapply(tab.index$value, list(tab.index$protein, tab.index$amin), mean)
  CV=100*(tab.sd/tab.moy)
  prot=rep(rownames(CV), ncol(CV))
  CV=stack(as.data.frame(CV))
  colnames(CV)=c("CV", "amin")
  CV$protein=prot
  CV=merge(CV, unique(tab.index[,c("protein", "amin", "stade", "index")]), c("protein", "amin"))
  CV$quanti=type_quanti[elem]
  CVtot=rbind.data.frame(CVtot, CV)
}

svg("./CV entre Rep_all prot_all index.svg",width=17, height = 24)
CVtot$quanti <- factor(CVtot$quanti, levels= c("sans.filtre","sans.pep.commun","sans.sdRT","sans.pepRep","sans.pepCor"))
CVindex <- CVtot[CVtot$index==index[i],]
CVindex$quanti <- as.factor(CVindex$quanti)
plot4 <- ggplot(CVindex,aes(x= as.factor(quanti),y=CV))+geom_boxplot()+
  ylab("CV (%) between replicates")+xlab("\n")+ ggtitle(index_title[i]) + theme(legend.position="none")+
  theme(axis.text=element_text(size=25),axis.title.y = element_text(size=25),
        plot.title=element_text(size=27,hjust = 0.5))+
  scale_x_discrete(labels=c("sans.filtre" = "Data 0", "sans.pep.commun" = "Data 1", "sans.sdRT" = "Data 2", "sans.pepRep" = "Data 3","sans.pepCor" = "Data 4"))
multiplot(plot1,plot3,plot5,plot2,plot4,cols=2)#l'ordre desplots est fait pour avoir ibaq,top3,average,averageLog et Model
dev.off()


svg("./CV entre Rep_all prot_all index par stade2.svg",width=28, height = 35)
CVindex <- CVtot[!CVtot$index=="modelssRep_fmolgFW",]
CVindex$index <- factor(CVindex$index, levels= c("ibaq_fmolgFW","TOP3_fmolgFW","ave_fmolgFW","avelog_fmolgFW","modelAcRep_fmolgFW"))
CVindex$quanti <- factor(CVindex$quanti, levels= c("sans.filtre","sans.pep.commun","sans.sdRT","sans.pepRep","sans.pepCor"))
levels(CVindex$quanti) <- c("Data 0", "Data 1", "Data 2","Data 3","Data 4")

plot<- ggplot(aes(index,CV,fill=index),data=CVindex)+geom_boxplot()+facet_grid(facets=stade~quanti)+ 
  theme(axis.text=element_text(size=40),
        axis.title.y=element_text(size=50),
        strip.text.x = element_text(size = 50),
        strip.text.y = element_text(size = 50),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ theme(panel.spacing = unit(1, "lines"))+
  ylab("CV (%) between replicates \n")
print(plot)
dev.off()

#scale_x_discrete(labels=c("ibaq_fmolgFW"="iBAQ","TOP3_fmolgFW"="TOP3","ave_fmolgFW"="Average","avelog_fmolgFW"="Average-Log","modelAcRep_fmolgFW"="Model"))+



##############              Precision = CV(%) entre les proteines ribosomales avec les abondances


load("./effet_filtre/data_isma_propre.RData")
load("./effet_filtre/index_isma_propre.RData")

protRib <- read.table("./RPL-RPS_FRIM.txt",header=T, sep="\t")
protRib_RPL <- protRib[protRib$categorie %in% "RPL",]
protRib_RPL <- unique(protRib_RPL$solyc)
prot_data=read.table("F:/These/Proteomic Gif Sur Yvette/iBAQ_FRIM_24012017/data_brutes/protein_dataFRIM.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
temp <- prot_data[prot_data$Solyc %in% protRib_RPL,]
protRib_RPL <- unique(temp$Protein.ID)#69
prot.rib.com <- intersect(index.sans.pepCor$protein,protRib_RPL)#38
index.sans.pepCor$Stade <- gsub("07.7","7.7",index.sans.pepCor$Stade)

mylist=list(index.sans.filtre, index.sans.pep.commun, index.sans.sdRT,index.sans.pepRep,index.sans.pepCor)
type_quanti=c("sans.filtre", "sans.pep.commun", "sans.sdRT", "sans.pepRep", "sans.pepCor")
index= c("ibaq", "top3", "average", "averageLog", "model_sans_rep", "model_avec_rep")
stade=c("7.7","15","21.7","28","34.3","41.3","48.5","50.3","53")
ratioRiballcom <- NULL

for (elem in 1:5){
  tab.index=mylist[[elem]]
  tab.indexribcom=tab.index[tab.index$protein %in% protRib_RPL,]
  tab.indexribcom$averageLog <- 10^tab.indexribcom$averageLog
  tab.indexribcom$model_sans_rep <- 10^tab.indexribcom$model_sans_rep
  tab.indexribcom$model_avec_rep <- 10^tab.indexribcom$model_avec_rep
  
  
  tab.indexribcom <- melt(tab.indexribcom[,c(1,3,6,8:12)],id=c("Stade","protein"))
  tab.indexribcom$Stade <- as.factor( tab.indexribcom$Stade)
  colnames( tab.indexribcom) <- c("Stade","protein","index","value")
  # ref <- tab.index$index == c("averageLog","model_sans_rep","model_avec_rep")
  # tab.index$value[ref] <- 10^(tab.index$value[ref])
  
  for (j in 1:9){
    
    subcom <-tab.indexribcom[tab.indexribcom$Stade == stade[j],]
    subcom <- ddply(subcom,.(index),transform,mean=mean(value,na.rm=T),ecty=sd(value,na.rm=T))
    subcom$CV <- subcom$ecty*100/subcom$mean
    subcom$quanti <- type_quanti[elem]
    ratioRiballcom <- rbind(ratioRiballcom,subcom)
  }}

dim(ratioRiballcom)#29016obs 7var

index= c("ibaq", "top3", "average", "averageLog", "model_sans_rep", "model_avec_rep")
index <- index[-5]
index_title <- c("iBAQ","Top3","Average","AverageLog", "Model")

svg(paste("./effet_filtre/Rib_RibprotCom_CV.svg"),width=18,height=22)
sub=ratioRiballcom[ratioRiballcom$index==index[i] ,]
sub$quanti <- factor(sub$quanti, levels= c("sans.filtre","sans.pep.commun","sans.sdRT","sans.pepRep","sans.pepCor"))
plot3 <- ggplot(sub,aes(x= as.factor(quanti),y=CV))+geom_boxplot()+theme(legend.title =  element_blank())+
  ylab("CV (%) between ribosomal protein \n")+xlab("\n")+ ggtitle(index_title[i]) +
  theme(axis.text=element_text(size=17, colour = "black"),
        axis.title.x=element_text(size=20, colour = "black"),
        axis.title.y=element_text(size=20, colour = "black"),
        plot.title=element_text(size=20, colour = "black"))+
  scale_x_discrete(labels=c("sans.filtre" = "", "sans.pep.commun" = "", "sans.sdRT" = "", "sans.pepRep" = "","sans.pepCor" = ""))
multiplot(plot1,plot3,plot5,plot2,plot4,cols=2)
dev.off()

##############              Precision = CV(%) entre les proteines ribosomales

load("./effet_filtre/Conc(fmolgFW)_index.RData")
protRib <- read.table("./RPL-RPS_FRIM.txt",header=T, sep="\t")
protRib_RPL <- protRib[protRib$categorie %in% "RPL",]
protRib_RPL <- unique(protRib_RPL$solyc)
prot_data=read.table("F:/These/Proteomic Gif Sur Yvette/iBAQ_FRIM_24012017/data_brutes/protein_dataFRIM.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
temp <- prot_data[prot_data$Solyc %in% protRib_RPL,]
protRib_RPL <- unique(temp$Protein.ID)#69
prot.rib.com <- intersect(conc.sans.pepCor$protein,protRib_RPL)#38
stade=c("07.7","15","21.7","28","34.3","41.3","48.5","50.3","53")

#merge pour avoir les description des proteines
prot_data=read.table("F:/These/Proteomic Gif Sur Yvette/iBAQ_FRIM_24012017/data_brutes/protein_dataFRIM.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
mylist=list(conc.sans.filtre,conc.sans.pep.commun,conc.sans.sdRT,conc.sans.pepRep ,conc.sans.pepCor)


stade=c("07.7","15","21.7","28","34.3","41.3","48.5","50.3","53")
index= c("ibaq_fmolgFW", "TOP3_fmolgFW", "ave_fmolgFW", "avelog_fmolgFW", "modelssRep_fmolgFW", "modelAcRep_fmolgFW")
index <- index[-5]
index_title <- c("iBAQ","Top3","Average","AverageLog", "Model")
type_quanti=c("sans.filtre", "sans.pep.commun", "sans.sdRT", "sans.pepRep", "sans.pepCor")

CVRib=NULL
for (elem in 1:5){
  tab.index=mylist[[elem]]
  colnames(tab.index) <- c("protein","msrun","Stade","index","value")
  tab.indexribcom=tab.index[tab.index$protein %in% protRib_RPL,]
  
  for (j in 1:9){
    
    subcom <-tab.indexribcom[tab.indexribcom$Stade == stade[j],]
    subcom <- ddply(subcom,.(index),transform,mean=mean(value,na.rm=T),ecty=sd(value,na.rm=T))
    subcom$CV <- subcom$ecty*100/subcom$mean
    subcom$quanti <- type_quanti[elem]
    CVRib <- rbind(CVRib,subcom)
  }}

dim(CVRib)#39468obs 7var

CVRib <- CVRib[!CVRib$index =="modelssRep_fmolgFW",]
CVRib$quanti <- factor(CVRib$quanti, levels= c("sans.filtre","sans.pep.commun","sans.sdRT","sans.pepRep","sans.pepCor"))


svg(paste("./effet_filtre/Rib_RibprotCom_CV_conc.svg"),width=18,height=22)
sub=CVRib[CVRib$index==index[i] ,]
sub$quanti <- factor(sub$quanti, levels= c("sans.filtre","sans.pep.commun","sans.sdRT","sans.pepRep","sans.pepCor"))
plot5 <- ggplot(sub,aes(x= as.factor(quanti),y=CV))+geom_boxplot()+theme(legend.title =  element_blank())+
  ylab("CV (%) between ribosomal protein \n")+xlab("\n")+ ggtitle(index_title[i]) +
  theme(axis.text=element_text(size=17, colour = "black"),
        axis.title.x=element_text(size=20, colour = "black"),
        axis.title.y=element_text(size=20, colour = "black"),
        plot.title=element_text(size=20, colour = "black"))+
  scale_x_discrete(labels=c("sans.filtre" = "", "sans.pep.commun" = "", "sans.sdRT" = "", "sans.pepRep" = "","sans.pepCor" = ""))
multiplot(plot1,plot3,plot5,plot2,plot4,cols=2)
dev.off()

test <- conc.sans.pepCor[conc.sans.pepCor$variable %in% "modelAcRep_fmolgFW",]
test_df <- dcast(test,protein~msrun, value.var = "value")
test_df$na <- apply(test_df, 1, function(x) sum(is.na(x)))
test_df <- test_df[!test_df$na == 26,]


#TOP3_fmolgFW       ibaq_fmolgFW       ave_fmolgFW        avelog_fmolgFW      modelAcRep_fmolgFW    modelssRep_fmolgFW

## iBAQ     top3     average      averageLog  Model
#
# sans filtre           2727      2455    2727          2727        2707
# filtre shared pep     2517      2162    2517          2517        2502
# filtre RT             2510      2155    2510          2510        2494
# filtre occurance      1587      1201    1587          1587        1586 
# fltre outlier         1254      911     1254          1254        1253




##############              plot peptides frim par filter (A finir)
load("./effet_filtre/data_isma_propre.RData")
load("./effet_filtre/index_isma_propre.RData")

mylist=list(index.sans.filtre, index.sans.pep.commun, index.sans.sdRT,index.sans.pepRep,index.sans.pepCor)
mylist2=list(tab.sans.filtre, tab.sans.pep.commun, tab.sans.sdRT,tab.sans.pepRep,tab.sans.pepCor)
palette(c("green", "darkorange","magenta","cyan", "darkblue", "darkorange", "red", "purple", "darkgreen", "black", sample(colors(), 50)))

prot <- unique(tab.sans.filtre$protein)
type_quanti=c("sans.filtre", "sans.pep.commun", "sans.sdRT", "sans.pepRep", "sans.pepCor")
stade_title<- c("07.7","15","21.7","28","34.3","41.3","48.5","50.3","53")
Data_title <- c("Data 1","Data 2","Data 3","Data 4","Data 5")
index= c("ibaq", "top3", "average", "averageLog", "model_sans_rep", "model_avec_rep")
index <- index[-5]

svg(paste("./plot_intensity_peptides_FRIM.svg"),width=28,height=39.2)
par(mfrow=c(6,5))
for (i in 1:12){
  
  par(mgp=c(3,0.7,0),mar=c(5.1,4.9,4.1,4.9)) 
  tabok=mylist2[[1]]
  sub=tabok[tabok$protein==prot[i],]
  sub$peprep=paste0(sub$peptiz, sub$Rep)
  a=unique(sub[,c("peprep", "peptide")])
  b=cbind.data.frame(peptide=unique(sub$peptide), color=palette()[1:length(unique(sub$peptide))])
  a=merge(a, b, "peptide")
  toto=tapply(sub$logQnorm, list(sub$peprep, sub$Stade), mean, na.rm=TRUE)
  temp=cbind.data.frame(peprep=rownames(toto), bidon=1)
  temp=merge(temp, a, "peprep")
  temp=temp[order(temp$peprep),]
  matplot(t(toto), type="l", lty=1, main=paste0(prot[i]), col=temp$color,ylab="", ylim=c(5,10),xaxt='n',yaxt='n',cex.lab=5,cex.main=2.5)
  #legend("topleft", legend=c("ibaq", "top3", "average", "averageLog", "model"), pch=c(17, 8, 18,16,15))
  axis(1, at=c(1:9), labels=F) 
  
  #mtext("Concentration of UPS (nM)", side=1, line=3.5, cex=1.2)
  mtext("Peptide ion intensity (log10 transformed)", side=2, line=2.5, cex=1.4)
  axis(2, cex.axis=1.8) 
  #axis(side=4, at= c(5,6,7,8,9,10),labels=format(10^(c(5,6,7,8,9,10)),scientific=T), cex.axis=1.8,las=1)
  axis(side=4, at= c(5,6,7,8,9,10),labels=F) 
  text(11.9,c(5,6,7,8,9,10),labels=format(10^(c(5,6,7,8,9,10)),scientific=T), cex=1.9,las=1,srt=-90,xpd=T)
  #text(12.5,7.5,"Protein abundance", cex=1.8,srt=-90,xpd=T)
  
  tab.index=mylist[[1]]
  
  ##?a bloque parici 
  sub=tab.index[tab.index$protein==prot[i],]
  sub$Stade <- as.factor(sub$Stade)
  toto_temp <- sub[,c("top3","Stade")]
  toto= aggregate(toto_temp, list(toto_temp$Stade), mean, na.rm=TRUE)
  par(new=TRUE)
  plot(log10(toto$top3), ylim=c(5,10), type="b", pch=8, axes=FALSE, xlab="", ylab="",lwd=2) 
  
  
  sub=tab.index[tab.index$protein==prot[i],]
  toto=tapply(sub$value[sub$index=="averageLog"], list(sub$Stade[sub$index=="averageLog"]), mean)
  par(new=TRUE)
  plot(log10(toto), ylim=c(5,10), type="b", pch=16, axes=FALSE, xlab="", ylab="",lwd=2) #rond
  
  sub=tab.index[tab.index$protein==prot[i],]
  toto=tapply(sub$value[sub$index=="ibaq"], list(sub$Stade[sub$index=="ibaq"]), mean)
  par(new=TRUE)
  plot(log10(toto), ylim=c(5,10),  type="b", pch=17, axes=FALSE, xlab="", ylab="",lwd=2) #triangle
  
  sub=tab.index[tab.index$protein==prot[i],]
  toto=tapply(sub$value[sub$index=="model_avec_rep"], list(sub$Stade[sub$index=="model_avec_rep"]), mean)
  par(new=TRUE)
  plot(log10(toto), ylim=c(5,10),  type="b", pch=15, axes=FALSE, xlab="", ylab="",lwd=2) # carré
  
  sub=tab.index[tab.index$protein==prot[i],]
  toto=tapply(sub$value[sub$index=="average"], list(sub$Stade[sub$index=="average"]), mean)
  par(new=TRUE)
  plot(log10(toto), ylim=c(5,10),type="b", pch=18, axes=FALSE, xlab="", ylab="",lwd=2) #losange
  
  for (elem in 2:5){
    tabok=mylist2[[elem]]
    
    sub=tabok[tabok$protein==prot[i],]
    if (nrow(sub)>0){
      sub$peprep=paste0(sub$peptiz, sub$Rep)
      toto=tapply(sub$logQnorm, list(sub$peprep, sub$Stade), mean, na.rm=TRUE)
      temp=cbind.data.frame(peprep=rownames(toto), bidon=1)
      temp=merge(temp, a, "peprep")
      temp=temp[order(temp$peprep),]
      matplot(t(toto), type="l", lty=1, main=paste0(descriptionUPS[i]), col=temp$color, ylab="",ylim=c(5,10),xaxt='n',yaxt='n',cex.lab=2,cex.main=2.5)
      axis(1, at=c(1:11), labels=F) 
      #mtext("Concentration of UPS1 (fmol.?l-1,logscale)", side=1, line=3.5, cex=1.2)
      #mtext("Intensity (Log10-transformed)", side=2, line=2.5, cex=1.2)
      axis(2, cex.axis=1.8) 
      axis(side=4, at= c(5,6,7,8,9,10),labels=F) 
      text(11.9,c(5,6,7,8,9,10),labels=format(10^(c(5,6,7,8,9,10)),scientific=T), cex=1.9,las=1,srt=-90,xpd=T)
      if (elem==5){    text(12.7,7.5,"Protein abundance", cex=2.5,srt=-90,xpd=T)
      }else{NULL}
      tab.index=mylist[[elem]]
      
      sub=tab.index[tab.index$protein==prot[i],]
      toto=tapply(sub$value[sub$index=="top3"], list(sub$Stade[sub$index=="top3"]), mean, na.rm=TRUE)
      par(new=TRUE)
      plot(log10(toto), ylim=c(5,10), type="b", pch=8, axes=FALSE, xlab="", ylab="",lwd=2) #rond
      
      sub=tab.index[tab.index$protein==prot[i],]
      toto=tapply(sub$value[sub$index=="averageLog"], list(sub$Stade[sub$index=="averageLog"]), mean)
      par(new=TRUE)
      plot(log10(toto), ylim=c(5,10), type="b", pch=16, axes=FALSE, xlab="", ylab="",lwd=2) #rond
      
      sub=tab.index[tab.index$protein==prot[i],]
      toto=tapply(sub$value[sub$index=="ibaq"], list(sub$Stade[sub$index=="ibaq"]), mean)
      par(new=TRUE)
      plot(log10(toto), ylim=c(5,10), type="b", pch=17, axes=FALSE, xlab="", ylab="",lwd=2) #triangle
      
      sub=tab.index[tab.index$protein==prot[i],]
      toto=tapply(10^(sub$value[sub$index=="model_avec_rep"]), list(sub$Stade[sub$index=="model_avec_rep"]), mean)
      par(new=TRUE)
      plot(log10(toto), ylim=c(5,10), type="b", pch=15, axes=FALSE, xlab="", ylab="",lwd=2) # carré
      
      sub=tab.index[tab.index$protein==prot[i],]
      toto=tapply(sub$value[sub$index=="average"], list(sub$Stade[sub$index=="average"]), mean)
      par(new=TRUE)
      plot(log10(toto), ylim=c(5,10), type="b", pch=18, axes=FALSE, xlab="", ylab="",lwd=2) #losange
      
      
      
    }else{
      plot.new()
    }
  }
}
dev.off()


# CV entre methodes -------------------------------------------------------

library(goeveg)
require("ggplot2")
cv_nan<-function(x) return(cv(x,na.rm = T)*100)
cv_meth<-aggregate(conc.sans.sdRT$value,by=list(conc.sans.sdRT$variable,conc.sans.sdRT$index),cv_nan)
colnames(cv_meth)<-c("Methode","Stade","CV")
cv_meth<-cv_meth[sort(as.numeric(cv_meth$Stade),index.return=T)$ix,]
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
ggplot(cv_meth, aes(fill=Methode, y=CV, x=factor(Stade,level=unique(sort(as.numeric(cv_meth$Stade)))))) +
  geom_bar(position="dodge", stat="identity")+theme+scale_fill_manual(values=c("red","blue","green","orange","grey","black"))+ylab("CV (%)")+xlab("Stade")
