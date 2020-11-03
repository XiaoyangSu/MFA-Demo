library(dplyr)
library(ggplot2)
library(minqa)


# Figure 2 Solving process of MFA
# This code shows how the optimum solution of fluxes were found by minimizing the sum of squared residue (SSR) value
# function calculating the labeling pattern when combining x and y
CP <- function(x,y) {
  nx <- length(x)
  ny <- length(y)
  n <- nx+ny-1
  Tx <- matrix(0,ncol=ny,nrow=n)
  Product <- rep(0,n)
  for (i in 1:ny) {
    Tx[i:(i+nx-1),i] <- x 
  }
  Product <- Tx %*% y
  return(t(Product))
}

# set the labeling pattern of glc[1-3] and glc[4-6] (use 1,2-13C2 glucose as tracer)
Glc123<- c(0,0,1,0)
Glc456<- c(1,0,0,0)

Glc<-rbind(Glc123,Glc456)
#function calculating the labeling pattern (MID) of FBP given the value of f3 and f5
MID <- function(f) {

  EMUMatrix_11 <- c(-(100+f[1]),0,f[1],0)
  EMUMatrix_12 <- c(0,-(100+f[1]),0,f[1])
  EMUMatrix_13 <- c((100+f[1]),0,-(100+f[1]+f[2]),f[2])
  EMUMatrix_14 <- c(0,(100+f[1]),(100+f[2]),-(200+f[1]+f[2]))
  
  EMUMatrix_1 <- matrix(c(EMUMatrix_11,EMUMatrix_12,EMUMatrix_13,EMUMatrix_14),nrow=4,ncol=4,byrow=T)
  
  EMUMatrix_1B <- matrix(c(-100,0,0,0,0,-100,0,0),nrow=4,ncol=2) %*% Glc
  
  EMU_Solve_1 <- solve(EMUMatrix_1) %*% EMUMatrix_1B
  
  FBP <- (CP(EMU_Solve_1[3,],EMU_Solve_1[4,])*f[1]+c(0,0,1,rep(0,4))*100)/(100+f[1])
  return(FBP)
}

#function calculating the squared residue (SSR) value given the value of f3 and f5
Res <- function(f) {
  sum((MID(f)-FBP_Label)^2)/0.003^2
}

# When f3=50 and f5=150, predict the labeling pattern of FBP (M+0,M+1...M+6)
TestFlux <- c(50,150)
FBP_Label <- MID(TestFlux) 

###############Optimization Path######
#Set the starting point as f3=150 and f5=200, use bobyqa algorithm to resolve f3 and f5
fluxes<-bobyqa(c(150,200),lower=c(0,0),upper=c(500,500),fn=Res)
#Set the starting point as f3=170 and f5=240, use bobyqa algorithm to resolve f3 and f5
fluxes.170.240<-bobyqa(c(170,240),lower=c(0,0),upper=c(500,500),fn=Res,control=list(iprint=3))
######################################
#Generate a map of SSR value, with all possible combinations of f3 (range in 0~200)and f5 (range in 0~900)
FluxSearch <- expand.grid(F3 = seq(0, 200, by = 0.5), F5 = seq(0, 900, by=1.5))
FluxSearch <- as.matrix(FluxSearch)
FluxSearch <- cbind(FluxSearch,matrix(0,nrow=nrow(FluxSearch),ncol=8))

for (i in 1:nrow(FluxSearch)) {
  Results <- MID(FluxSearch[i,1:2])
  FluxSearch[i,3:9] <- Results
}

for (i in 1:nrow(FluxSearch)) {
  FluxSearch[i,10] <- sum((FluxSearch[i,3:9]-FBP_Label)^2)/0.003^2
}

colnames(FluxSearch)[10] <- c("Res")

Optimization.Path1 <- matrix(0,nrow=20,ncol=2)
Optimization.Path1[,1] <- c(seq(170,50,by=-10),rep(50,7)) 
Optimization.Path1[,2] <- c(240,360,350,340,320,310,290,280,260,230,210,180,rep(150,8))
Optimization.Path1 <- data.frame(Optimization.Path1)
colnames(Optimization.Path1) <- c("Path1.f3","Path1.f5")

# plot the map with optimization path
FluxSearch %>% data.frame() %>% tbl_df() %>% 
  #filter(Res<=qchisq(0.95,df=1)) %>%
  ggplot(aes(x=F3,y=F5)) + 
  geom_tile(aes(fill=log10(1+Res))) + 
  geom_path(aes(x=Path1.f3,y=Path1.f5),data=Optimization.Path1) +
  scale_fill_gradientn(colors=rainbow(7),name="") + 
  theme_classic(base_size=25) +
  scale_x_continuous(~f[3],limits = c(0,200)) +
  scale_y_continuous(~f[5],limits=c(0,900)) +
  theme(axis.line =element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = 'black', size = 1),
        legend.position = "none",legend.key.size = unit(8,units = "mm"))
ggsave(plot,filename=paste("FBP_Large_",i,".png",sep=""),width = 7,height=7,dpi=1200)


#Figure 3. Evaluation on the tracer selection for MFA 
#This section generate the region of confidence interval
#Plot the SSR map with a threshold (95% CI)

plot <- FluxSearch %>% data.frame() %>% tbl_df() %>% 
  filter(Res<=qchisq(0.95,df=1)*100/9) %>% #this is the threshold of 95% CI
  ggplot(aes(x=F3,y=F5)) + 
  geom_tile(aes(fill=Res)) + 
  scale_fill_gradientn(colors=rainbow(7),name="") + 
  theme_classic(base_size=25) +
  scale_x_continuous(~f[3],limits = c(0,200)) +
  scale_y_continuous(~f[5],limits=c(0,900)) +
  theme(axis.line =element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = 'black', size = 1),
        legend.position = c(0.8,0.3),legend.key.size = unit(8,units = "mm"))

ggsave(plot,filename=paste("50-Glc",".png",sep=""),width = 7,height = 6,dpi=600)

# To plot the graph from different tracer, change the labeling pattern in Line 23 and rerun the entire code.
