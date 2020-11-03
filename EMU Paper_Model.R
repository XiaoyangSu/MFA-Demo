library(dplyr)
library(ggplot2)
library(numDeriv)

# Plotting flux confidence region for the network from Antoniewicz et al. (2007). 

# Describing the labeling pattern of the tracer.
# In this example, we use 2-13C1-A as the tracer.
# We need the labeling patterns of A2, A3, A23 and A123
Tracer <- list(
  rbind(c(0,1),c(1,0)),
  c(0,1,0),
  c(0,1,0,0)
)

# Function calculating the labeling pattern when combining x and y
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

# Function that generates the labeling of all EMUs
# MID function takes the free fluxes, i.e. c(f2,f4) as the input

MID <- function(f) {
  
  # f1 is fixed as 100
  # f3,f5 and f6 are calculated from f1, f2 and f4
  V <- rep(0,6)
  V[1] <- 100
  V[2] <- f[1]
  V[4] <- f[2]
  
  V[3] <- V[2]-V[1]+2*V[4]
  V[5] <- V[4]
  V[6] <- V[1]-V[4]
  
  # If the combination of f2 and f4 makes other flux negative, the zero labeling is returned
  if(min(V)<=0) {
    EMU_Solve_1 <- matrix(0,ncol=2,nrow=5)
    EMU_Solve_2 <- matrix(0,ncol=3,nrow=2)
    EMU_Solve_3 <- matrix(0,ncol=4,nrow=3)
    
    colnames(EMU_Solve_1) <- c("M+0","M+1")
    colnames(EMU_Solve_2) <- c("M+0","M+1","M+2")
    colnames(EMU_Solve_3) <- c("M+0","M+1","M+2","M+3")
    
    rownames(EMU_Solve_1) <- c("C1","B2","D2","B3","D3")
    rownames(EMU_Solve_2) <- c("D23","B23")
    rownames(EMU_Solve_3) <- c("F123","D123","B123")
    return(list(EMU_Solve_1,EMU_Solve_2,EMU_Solve_3))
  }
  
  # The EMU matrices are constructed same as the paper Antoniewicz et al. (2007)
  EMUMatrix_11<- matrix(c(-V[4],0,0,0,V[5]),nrow=5,ncol=1)
  EMUMatrix_12<- matrix(c(V[4],-(V[1]+V[3]),V[2],0,0),nrow=5,ncol=1)
  EMUMatrix_13<- matrix(c(0,V[3],-(V[2]+V[5]),0,0),nrow=5,ncol=1)
  EMUMatrix_14<- matrix(c(0,0,V[5],-(V[1]+V[3]),V[2]),nrow=5,ncol=1)
  EMUMatrix_15<- matrix(c(0,0,0,V[3],-(V[2]+V[5])),nrow=5,ncol=1)
  
  EMUMatrix_1 <- cbind(EMUMatrix_11,EMUMatrix_12,EMUMatrix_13,EMUMatrix_14,EMUMatrix_15)
  EMUMatrix_1B <- matrix(c(0,-V[1],0,0,0,0,0,0,-V[1],0),nrow=5,ncol=2) %*% Tracer[[1]]
  
  EMU_Solve_1 <- solve(EMUMatrix_1) %*% EMUMatrix_1B
  
  EMU_Input_2 <- rbind(CP(EMU_Solve_1[4,],EMU_Solve_1[1,]),Tracer[[2]])
  
  EMU_Solve_2 <- solve(matrix(c(-(V[5]+V[2]),V[3],V[2],-(V[1]+V[3])),nrow=2,ncol=2)) %*% 
    matrix(c(-V[5],0,0,-V[1]),nrow=2,ncol=2) %*% EMU_Input_2
  
  EMU_Input_3 <- rbind(CP(EMU_Solve_2[2,],EMU_Solve_1[1,]),Tracer[[3]])
  
  EMU_Solve_3 <- solve(matrix(c(-(V[5]+V[2]-V[3]),0,0,(V[5]+V[2]-V[3]),-(V[5]+V[2]),V[3],0,V[2],-(V[1]+V[3])),nrow=3,ncol=3)) %*% 
    matrix(c(0,-V[5],0,0,0,-V[1]),nrow=3,ncol=2) %*% EMU_Input_3
  
  # Label the metabolites and the fractions
  colnames(EMU_Solve_1) <- c("M+0","M+1")
  colnames(EMU_Solve_2) <- c("M+0","M+1","M+2")
  colnames(EMU_Solve_3) <- c("M+0","M+1","M+2","M+3")
  
  rownames(EMU_Solve_1) <- c("C1","B2","D2","B3","D3")
  rownames(EMU_Solve_2) <- c("D23","B23")
  rownames(EMU_Solve_3) <- c("F123","D123","B123")
  
  # Return the results of labeling
  return(list(EMU_Solve_1,EMU_Solve_2,EMU_Solve_3))
}

# The test flux combination is the one that is used in the paper 
Testflux <- c(110,20)

# Generate the labeling patterns
F_Label <- MID(Testflux)[[3]][1,]
B_Label <- MID(Testflux)[[3]][3,]

# Generate the flux combinations
FluxSearch <- expand.grid(F2 = seq(50, 500, by = 1), F4 = seq(10, 40, by=0.1))

FluxSearch <- as.matrix(FluxSearch)

FluxSearch <- cbind(FluxSearch,matrix(0,nrow=nrow(FluxSearch),ncol=11))

# Generate the labeling patterns for the flux combinations
for (i in 1:nrow(FluxSearch)) {
  Results <- MID(FluxSearch[i,1:2])[[3]]
  FluxSearch[i,3:6] <- Results[2,]
  FluxSearch[i,7:10] <- Results[3,]
}

for (i in 1:nrow(FluxSearch)) {
  FluxSearch[i,11] <- sum((FluxSearch[i,3:10]-c(F_Label,B_Label))^2)/0.003^2
}

for (i in 1:nrow(FluxSearch)) {
  FluxSearch[i,12] <- sum((FluxSearch[i,3:6]-c(F_Label))^2)/0.003^2
}

for (i in 1:nrow(FluxSearch)) {
  FluxSearch[i,13] <- sum((FluxSearch[i,7:10]-c(B_Label))^2)/0.003^2
}

colnames(FluxSearch)[11:13] <- c("BF_Res","F_Res","B_Res")

# Plot the confidence regions
Plot_BF <- FluxSearch %>% data.frame() %>% tbl_df() %>% filter(BF_Res<=qchisq(0.99,df=1)) %>%
  ggplot(aes(x=F2,y=F4)) + 
  geom_tile(aes(fill=BF_Res)) + 
  scale_fill_gradientn(colors=rainbow(7),name="") + 
  theme_classic(base_size=25) +
  scale_x_continuous(~f[2],limits = c(50,500)) +
  scale_y_continuous(~f[4],limits=c(10,40)) +
  theme(axis.line =element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = 'black', size = 1),
        legend.position = c(0.8, 0.4),legend.key.size = unit(8,units = "mm"))


Plot_F <- FluxSearch %>% data.frame() %>% tbl_df() %>% filter(F_Res<=qchisq(0.95,df=1)) %>%
  ggplot(aes(x=F2,y=F4)) + 
  geom_tile(aes(fill=F_Res)) + 
  scale_fill_gradientn(colors=rainbow(7),name="") + 
  theme_classic(base_size=25) +
  scale_x_continuous(~f[2],limits = c(50,500)) +
  scale_y_continuous(~f[4],limits=c(10,40)) +
  theme(axis.line =element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = 'black', size = 1),
        legend.position = c(0.8, 0.3),legend.key.size = unit(8,units = "mm"))

Plot_B <- FluxSearch %>% data.frame() %>% tbl_df() %>% filter(B_Res<=qchisq(0.95,df=1)) %>%
  ggplot(aes(x=F2,y=F4)) + 
  geom_tile(aes(fill=B_Res)) + 
  scale_fill_gradientn(colors=rainbow(7),name="") + 
  theme_classic(base_size=25) +
  scale_x_continuous(~f[2],limits = c(50,500)) +
  scale_y_continuous(~f[4],limits=c(10,40)) +
  theme(axis.line =element_line(colour = 'black', size = 2),
        axis.ticks = element_line(colour = 'black', size = 1),
        legend.position = c(0.8, 0.4),legend.key.size = unit(8,units = "mm"))

ggsave(Plot_BF,filename=paste("Plot_BF",".png",sep=""),width = 7,height = 6,dpi=600)
ggsave(Plot_F,filename=paste("Plot_F",".png",sep=""),width = 7,height = 6,dpi=600)
ggsave(Plot_B,filename=paste("Plot_B",".png",sep=""),width = 7,height = 6,dpi=600)


