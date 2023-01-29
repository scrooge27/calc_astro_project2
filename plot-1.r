tab1=read.table("C:\\Users\\simop\\Desktop\\UNI\\II anno\\calc-astro\\project_2\\pos_tm2.txt")

t=tab1[,1]

x=tab1[,2]
y=tab1[,4]

vx=tab1[,3]
vy=tab1[,5]

#plot della traiettoria
plot(x,y,pch=16,cex=0.7,col="blue",xlab="x (AU)", ylab="y (AU)")

#plot della posizione-x (non richiesto)
#plot(t,x,pch=16,cex=0.7,col="blue",xlab="t (T)", ylab="x (AU)")

#plot della posizione-y (non richiesto)
#plot(t,y,pch=16,cex=0.7,col="blue",xlab="t (T)", ylab="y (AU)")
