tab1=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_2\\jacobi1.txt")

t=tab1[,1]

cj=tab1[,2]

plot(t,cj,pch=16,cex=0.7,col="red",xlab="t", ylab="cj")
