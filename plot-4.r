tab1=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_2\\nrg_tm2.txt")
tab2=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_2\\nrg_ds2.txt")

t=tab1[,1]#/10**13
e=tab1[,2]#/10**20

plot(t,e,pch=16,cex=0.7,col="green",xlab="t", ylab="e")

r2=tab2[,1]
e=tab2[,2]#/10**20

plot(r2,e,pch=16,cex=0.7,col="green",xlab="r2", ylab="e")
