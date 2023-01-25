tab1=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_2\\dis_tm1.txt")

t=tab1[,1]

x=tab1[,2]
y=tab1[,3]


#plot delle distanze-x 
plot(t,x,pch=16,cex=0.7,col="red",xlab="t", ylab="x (AU)")
      #,xlim=c(0,3),ylim=c(5,7.5))

#plot delle distanze-y
plot(t,y,pch=16,cex=0.7,col="red",xlab="t", ylab="y (AU)")
      #,xlim=c(0,3),ylim=c(9.5,12)