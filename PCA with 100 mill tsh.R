#_______Volátiles con áreas post filtro de 100 millones para pulque joven_____________________#
#This analysis is carried out after filtering out those detected compounds with less tha 100,000,000 Arbitrary Units in order to
#create a Principal Component Analysis that explained a greater portion of the original variability in the data set. Original detected compounds 
#can be shared, contact information in my main profile page.

vol100=read.csv("IX. Areas post filtrado con umbral de 100 mill.csv")
#vol100

estacion=row.names(vol100)
estacion
#Which latin number associates with which sampling region is specified on
#C:/Users/VT/Documents/UAQ/Proyecto de investigación/Tesis/Construcción del proyecto/Construcción tesis/Resultados/Cromatografía/Segundo acercamiento de volátiles/Heatmap & PCA with 100 mill tsh
#in "Relacion de muestreos con numeros latinos asignados en R.xlsx" In this case,
#1=4.1, 2=4.2, 3=1.1, 4=1.2, 5=2.1, 6=2.2, 7=3.1, 8=3.2, 9=5.1 and 10=5.2. Further details
#available at the cited article or GitHub information.

#We can easily visualize the name of each X sampling region.

names(vol100)
#to see the head of the data set columns. We are analyzing 53 features or variables.

apply(vol100, 2, mean)
#apply allows us to apply a certain function (mean in this case); 1 will be used
#if we are trying to compute the mean of the rows and 2 for computing the mean of the columns.
##WATCH OUT: the means are vastly different=WE NEED MEAN ZERO AND VARIANCE 1 TO 
#PROCEED WITH PCA

apply(vol100, 2, var)
#Let´s check out the variance of each feature/variable. 
#They vary greatly too.

#Now is important to clarify that the units are not different among the
#features/factors/variables studied in this case.

#However, first, we need to scale; if not, most of the principal components
#would be driven by the features with the greatest mean and variance.

#Let´s standardize by adjusting mean zero and standard deviation 1 with the following command:

pr.out=prcomp(vol100, scale=TRUE)
#Now the variables have mean zero by the prcomp comand and standard deviation (S.D.) 
#1 with scale=TRUE.

names(pr.out)
#center and scale refer to the means and S.D. of the variables used for scaling prior
#to PCA .

#The rotation matrix provides the PCA loadings, each column contains the corresponding 
#principal component loading vector.

#The X refers to the principal component score vectors.

pr.out$center

pr.out$scale

pr.out$rotation

#We got multiple PCs (10), which is expected as there will be min(n-1,p) principal components for
#observations (n) and features (p).

#In this case, it is partially possible to define which compounds mostly drive the first
#two PCs. We describe the compounds with the greatest influence on each component
#and in parenthesis its loading vector and its sign.

#PC1: propanol 2-methyl (-0.208), butanol 3-methyl(-0,22), octanoic acid ethyl ester (-0.23),
#ispentyl hexanoate (-0.21), n-caprylic acid siobutyl ester (-0.23), octanoic acid 3-methyl butyl ester (-0,20),
#phenylethyl alcohol (-0.22), octanoic acid (-0.20), hexadecanoic acid ethyl ester (-0.21), 
#octanoic acid 2-phenylethyl ester (-0.21), cyclopentadecanone 2-hydroxy (-0.21)

#PC2:  octadecandienoic acid (0.23), ehtyl 9-hexadecenoate (0.21), heptadecanoic acid
#ethyl ester (0.21), octadecenoic acid ethyl ester (0.20), pentadecanoic acid ethyl ester (0.21),
#tetradecanoic acid ethyl ester (-0.21), benzenepropanoic acid ethyl ester (0.20), acetic acid 2-phenyl ester (-0.19),
#cyclodecane (0.23), decanoic acid propyl ester (0.23), hexanoic acid ethyl ester (0.21).



#Let's check now if only two PCs suffice to explain at least 50% of the 
#original variability of the data. *Spoiler alert*: they don't, but they get pretty 
#close...

#Let´s plot the first two Principal Components (PC):

biplot(pr.out, scale=0)
#scale=0 ensures the arrows are scaled to represent the loadings.
#We can change the sign of the PCs without any relevant consequence, but let's not
#do it for now...

#pr.out$rotation=-pr.out$rotation
#pr.out$x=-pr.out$x
#biplot(pr.out, scale=0)

#The prcomp function also gives us the S.D. of each P.C., being those:

pr.out$sdev

#The variance explained by each P.C. is obtained by squaring:

pr.var=pr.out$sdev^2
pr.var

#Now, the proportion of variance explained by each P.C. is obtained after dividing
#the variance explained by each P.C. by the total variance explained by the X number of
#P.C.s:

pve=pr.var/sum(pr.var)
pve

#Then, the 1st P.C. explains 29.49% of the variance in the data, the 2nd the 16.05% and so forth.

#Finally, we can plot the PVE explained by each P.C., as well as the cumulative PVE:

par(mfrow=c(1,2))
plot(pve, xlab="Principal Component", ylab="Proportion of Variance Explained", ylim=c(0,1), type="b")
plot(cumsum(pve), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", ylim=c(0,1), type="b")
#cumsum computes the cumulative sum of the elements of a numeric vector. E.g. a=c(1,2,8,-3)
#cumsum(a)=1 3 11 8


#Returning to our analysis, above we have our scree plot to figure out the optimal number of
#P.C.s we'll require.In this case, two may suffice (they both explain 45.54 % of the variance in data).

##########CLUSTERING BY K-MEANS ############

library(ggplot2)

clvol100mill=pr.out$x[,1:2]
#We select the first two score vectors of each of the two main Principal
#Components!! Then, we can proceed with k-means clustering!!!

clvol100mill

set.seed(19)
wcss <- vector()
for(i in 1:9){
  wcss[i] <- sum(kmeans(clvol100mill, i)$withinss)
}

wcss

ggplot() + geom_point(aes(x = 1:9, y = wcss), color = 'blue') + 
  geom_line(aes(x = 1:9, y = wcss), color = 'blue') + 
  ggtitle("Elbow method for volatiles K-means") + 
  xlab('Number of Centroids k') + 
  ylab('WCSS')

#We´d need to use 6 or 7 centroids...

#It means that only a certain amount of samples group together, meaning
#most of the sampled pulque's aroma is unique. Curious enough, the number of centroids
#doesn't match the number of them found analyzing physical and microbiological parameters!!
#(5 centroids with WCSS=3.51 were found over that analysis) :O..

set.seed(20)
km.out=kmeans(clvol100mill, 6, iter.max=1000, nstart=50)
#Here, we specify to work with only 6 clusters. It is suggested to use nstar=20 or 50, and 
##iter.max refers to the maximum iterations to apply to the algorithm and 
#nstart to the amount of centroid groups that are being employed internally or
#multiple initial cluster assignments (20 or 50).
#The last consideration will allow us to find a global optimum, rather than a local one,
#measured in the reduction of within cluster sumed squares, as the nstart value increases.

km.out$cluster
#Here we check the assignments of the n observations. We see how the data is sepa-
#rated into distinctive groups (as expected).

km.out
#We can see the sumary of the clustering run...

km.out$tot.withinss

#Now, we can plot our results...

plot(clvol100mill, col=(km.out$cluster +1), main="K.Means Clustering with K=6", xlab="Principal Component 1 (29.49%)", ylab="Principal Component 2 (16.05%)", pch=20, cex=2)

#Let´s try to incorporate labels to the clustering plot:

text(clvol100mill, row.names(vol100), cex=0.6, pos=4, col="red")


#Let´s try 7 clusters...

set.seed(21)
km.out=kmeans(clvol100mill, 7, iter.max=1000, nstart=50)
#Here, we specify to work with only 10 clusters. It is suggested to use nstar=20 or 50, and 
##iter.max refers tot the maximum iterations to apply to the algorithm and 
#nstart to the amount of centroid groups that are being employed internally or
#multiple initial cluster assignments (20 or 50).
#The last consideration will allow us to find a global optimum, rather than a local one,
#measured in the reduction of within cluster sumed squares as the nstart value increases.

km.out$cluster
#Here we check the assignments of the n observations. We see how the data is sepa-
#rated into distinctive groups (as expected).

km.out
#We can see the sumary of the clustering run...

km.out$tot.withinss

#Now, we can plot our results...

plot(clvol100mill, col=(km.out$cluster +1), main="K.Means Clustering with K=7", xlab="Principal Component 1(29.49%)", ylab="Principal Component 2(16.05%)", pch=20, cex=2)

text(clvol100mill, row.names(vol100), cex=0.6, pos=4, col="red")

#Based on the Elbow method, the difference on the wcss
#between 6 or 7 clusters is 1.98. Therefore, 
#6 clusters could already describe/group the data properly, considering also the 
#sampling findings.

#_#_#_#_#_____________________#######_____________________________#_#_#_#_#
