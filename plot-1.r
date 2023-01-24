tab1=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_2\\3corpi.txt")

x=tab1[,2]
y=tab1[,4]

vx=tab1[,3]
vy=tab1[,5]

#calcolo il delta t
t=tab1[,1]

#plot della posizione
plot(x,y,pch=16,cex=0.7,col="blue",xlab="x", ylab="y")

