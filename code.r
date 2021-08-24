  #For sample size 150
  T1=NULL
  T2=NULL
  T3=NULL
  set.seed(10)
  for(k in 1:10000)
  { 
    a=0.9   
    b=0.8#....value of beta#   
    G=2   
    lambda_A=2.5   
    lambda_B=2.5   
    lambda_AB=1.5   
    lambda_BA=1.5 
    X=matrix(nrow=150,ncol=1)  
    for(i in 1:150)   
    {     
      X[i,]=runif(1,min=0,max=1)   
    }   
    X  
    u=matrix(nrow=150,ncol=1)   
    for(i in 1:150)   
    {     
      if(X[i,]<=(1/3))     
      {      
        u[i,]=0     
      }     
      else if((X[i,]>(1/3))&&(X[i,]<=(2/3)))     
      {       
        u[i,]=1     
      }    
      else     
      {      
        u[i,]=2     
      }   
    }  
    u   
    X1=matrix(nrow=75,ncol=1)#....matrix for storing values of X_A#
    for(i in 1:75)   
    {    
      X1[i,]=rpois(1,lambda_A*(a^(G-u[i,])))   
    }   
    X1   
    X2=matrix(nrow=75,ncol=1)#....matrix for storing values of X_B#
    for(j in 1:75)   
    {    
      X2[j,]=rpois(1,lambda_B*(a^(G-u[j+75,])))   
    }   
    X2   
    X3=matrix(nrow=75,ncol=1)#...matrix for storing values of X_AB#
    for(i in 1:75)   
    {     
      X3[i,]=rpois(1,(((lambda_AB+(0.8)*(lambda_AB-lambda_A))*a^(G-u[i,])+(0.8)*X1[i,])/(1+0.8)))   
    }   
    X3   
    X4=matrix(nrow=75,ncol=1)#...matrix for storing values of X_BA#
    for(j in 1:75)   
    {     
      X4[j,]=rpois(1,(((lambda_BA+(0.8)*(lambda_BA-lambda_B))*a^(G-u[j+75,])+(0.8)*X2[j,])/(1+0.8)))   
    }   
    X4   
    S1=sum(X1)   
    S2=sum(X2)   
    S3=sum(X3)   
    S4=sum(X4)   
    N=150   
    N1=75   
    N2=75   
    N3=75   
    N4=75   
    S1   
    S2   
    S3   
    S4   
    c=(1/3)*(1+a+(a^2))#....value of common probability Pie#   
    c   
    lambdahat_A=(S1+0.5)/(c*(N1+1))   
    lambdahat_B=(S2+0.5)/(c*(N2+1))   
    lambdahat_AB=(S3+0.5)/(c*(N3+1))   
    lambdahat_BA=(S4+0.5)/(c*(N4+1))   
    lambdahat_A   
    lambdahat_B   
    lambdahat_AB   
    lambdahat_BA   
    Z=(1/2)*(log(lambdahat_B)-log(lambdahat_A)+log(lambdahat_AB)-log(lambdahat_BA))   
    lambda1hat=(lambdahat_A+lambdahat_B)/2
    lambda2hat=(lambdahat_AB+lambdahat_BA)/2
    d=((1/3)*(a^4+a^2+1))-((1/9)*(a^2+a+1)^2)#....value of pie*#
    d   
    sigma22hat=((1/c)^2*((lambda2hat*c)+((b^2*lambda1hat*c+(1+b)^2*lambda2hat^2*d)/(1+b)^2)))
    sigma11hat=((lambda1hat*(c+lambda1hat*d))/c^2)
    sigma12hat=(((b/(1+b))*c+lambda2hat*d)*(lambda1hat/c^2))
    sigma0hat=(sigma11hat/(2*lambda1hat^2))+(sigma22hat/(2*lambda2hat^2))+(sigma12hat/(lambda1hat*lambda2hat))
    T1[k]=(sqrt(N)*Z)/(sqrt(2)*sqrt(sigma0hat))
    D1=(sqrt(N)*(lambdahat_B-lambdahat_A))/(2*sqrt(sigma11hat))
    D2=(sqrt(N)*(lambdahat_AB-lambdahat_BA))/(2*sqrt(sigma22hat))
    r=-(sigma12hat/sqrt(sigma11hat*sigma22hat))#....value of estimated correlation coefficient#   
    if((D1>0)&&(D2>0))   
    {
      T2[k]=sqrt(D1^2+(D2^2)-(2*r*D1*D2))/sqrt(1-(r^2))   
    }   
    else if((D1<=0)&&(D2>=D1))   
    {     
      T2[k]=(D2-(r*D1))/sqrt(1-(r^2))   
    }   
    else   
    {     
      T2[k]=(D1-(r*D2))/sqrt(1-(r^2))   
    }   
    T3[k]=max(D1,D2) 
  } 
  #Finding the upper 5% point when the test statistic is T2
  c21=sort(T2)[9500]
  #Finding the upper 5% point when the test statistic is T3
  c31=sort(T3)[9500]
  #Checking empirical power for the same sample size by taking an alternative 
  T1=NULL 
  T2=NULL 
  T3=NULL 
  set.seed(10) 
  for(k in 1:10000) 
  {   
    a=0.9   
    b=0.8#....value of beta#   
    G=2   
    lambda_A=2.5   
    lambda_B=3   
    lambda_AB=2   
    lambda_BA=1.5   
    X=matrix(nrow=150,ncol=1)   
    for(i in 1:150)   
    {    
      X[i,]=runif(1,min=0,max=1)   
    }   
    X   
    u=matrix(nrow=150,ncol=1)   
    for(i in 1:150)   
    {     
      if(X[i,]<=(1/3))     
      {       
        u[i,]=0     
      }    
      else if((X[i,]>(1/3))&&(X[i,]<=(2/3)))     
      {       
        u[i,]=1     
      }    
      else     
      {      
        u[i,]=2     
      }   
    }   
    u   
    X1=matrix(nrow=75,ncol=1)#....matrix for storing values of X_A#
    for(i in 1:75)   
    {     
      X1[i,]=rpois(1,lambda_A*(a^(G-u[i,])))   
    }   
    X1   
    X2=matrix(nrow=75,ncol=1)#....matrix for storing values of X_B#
    for(j in 1:75)   
    {    
      X2[j,]=rpois(1,lambda_B*(a^(G-u[j+75,])))   
    }   
    X2   
    X3=matrix(nrow=75,ncol=1)#...matrix for storing values of X_AB#
    for(i in 1:75)   
    {    
      X3[i,]=rpois(1,(((lambda_AB+(0.8)*(lambda_AB-lambda_A))*a^(G-u[i,])+(0.8)*X1[i,])/(1+0.8)))   
    }   
    X3  
    X4=matrix(nrow=75,ncol=1)#...matrix for storing values of X_BA#
    for(j in 1:75)   
    {     
      X4[j,]=rpois(1,(((lambda_BA+(0.8)*(lambda_BA-lambda_B))*a^(G-u[j+75,])+(0.8)*X2[j,])/(1+0.8)))   
    }   
    X4   
    S1=sum(X1)   
    S2=sum(X2)   
    S3=sum(X3)   
    S4=sum(X4)   
    N=150   
    N1=75   
    N2=75   
    N3=75   
    N4=75   
    S1   
    S2   
    S3   
    S4   
    c=(1/3)*(1+a+(a^2))#....value of common probability Pie#   
    c   
    lambdahat_A=(S1+0.5)/(c*(N1+1))   
    lambdahat_B=(S2+0.5)/(c*(N2+1))   
    lambdahat_AB=(S3+0.5)/(c*(N3+1))  
    lambdahat_BA=(S4+0.5)/(c*(N4+1))   
    lambdahat_A   
    lambdahat_B   
    lambdahat_AB   
    lambdahat_BA   
    Z=(1/2)*(log(lambdahat_B)-log(lambdahat_A)+log(lambdahat_AB)-log(lambdahat_BA))   
    lambda1hat=(lambdahat_A+lambdahat_B)/2
    lambda2hat=(lambdahat_AB+lambdahat_BA)/2
    d=((1/3)*(a^4+a^2+1))-((1/9)*(a^2+a+1)^2)#....value of pie*#
    d   
    sigma22hat=((1/c)^2*((lambda2hat*c)+((b^2*lambda1hat*c+(1+b)^2*lambda2hat^2*d)/(1+b)^2)))
    sigma11hat=((lambda1hat*(c+lambda1hat*d))/c^2)
    sigma12hat=(((b/(1+b))*c+lambda2hat*d)*(lambda1hat/c^2))
    sigma0hat=(sigma11hat/(2*lambda1hat^2))+(sigma22hat/(2*lambda2hat^2))+(sigma12hat/(lambda1hat*lambda2hat))
    T1[k]=(sqrt(N)*Z)/(sqrt(2)*sqrt(sigma0hat))
    D1=(sqrt(N)*(lambdahat_B-lambdahat_A))/(2*sqrt(sigma11hat))
    D2=(sqrt(N)*(lambdahat_AB-lambdahat_BA))/(2*sqrt(sigma22hat))
    r=-(sigma12hat/sqrt(sigma11hat*sigma22hat))#....value of estimated correlation coefficient#   
    if((D1>0)&&(D2>0))   
    {     
      T2[k]=sqrt(D1^2+(D2^2)-(2*r*D1*D2))/sqrt(1-(r^2))   
    }   
    else if((D1<=0)&&(D2>=D1))   
    {     
      T2[k]=(D2-(r*D1))/sqrt(1-(r^2))   
    }   
    else   
    {    
      T2[k]=(D1-(r*D2))/sqrt(1-(r^2))   
    }   
    T3[k]=max(D1,D2) 
  }
  s=0 
  for(i in 1:10000) 
  {   
    if(T1[i]>1.645)   
    {     
      s=s+1   
    } 
  } 
  p11=(s/10000)#...Empirical Power#
  c=0 
  for(i in 1:10000) 
  {   
    if(T2[i]>c21)   
    {     
      c=c+1   
    } 
  } 
  p21=(c/10000)#...Empirical Power#
  r=0 
  for(i in 1:10000) 
  {   
    if(T3[i]>c31)   
    {     
      r=r+1   
    } 
  } 
  p31=(r/10000)#...Empirical Power#
  #For Sample Size 200
  T1=NULL
  T2=NULL 
  T3=NULL 
  set.seed(10) 
  for(k in 1:10000) 
  { 
    a=0.9 
    b=0.8#....value of beta#
    G=2 
    lambda_A=2.5 
    lambda_B=2.5 
    lambda_AB=1.5 
    lambda_BA=1.5 
    X=matrix(nrow=200,ncol=1)
    for(i in 1:200) 
    {
      X[i,]=runif(1,min=0,max=1) 
    }
    X
    u=matrix(nrow=200,ncol=1) 
    for(i in 1:200)
    {
      if(X[i,]<=(1/3)) 
      {
        u[i,]=0
      }
      else if((X[i,]>(1/3))&&(X[i,]<=(2/3)))
      {
        u[i,]=1 
      }
      else 
      { 
        u[i,]=2 
      }
    }
    u
    X1=matrix(nrow=100,ncol=1)#....matrix for storing values of X_A#
    for(i in 1:100)
    {
      X1[i,]=rpois(1,lambda_A*(a^(G-u[i,]))) 
    }
    X1
    X2=matrix(nrow=100,ncol=1)#....matrix for storing values of X_B# 
    for(j in 1:100) 
    {
      X2[j,]=rpois(1,lambda_B*(a^(G-u[j+100,])))
    }
    X2
    X3=matrix(nrow=100,ncol=1)#...matrix for storing values of X_AB# 
    for(i in 1:100)
    {
      X3[i,]=rpois(1,(((lambda_AB+(0.8)*(lambda_AB-lambda_A))*a^(G-u[i,])+(0.8)*X1[i,])/(1+0.8))) 
    }
    X3
    X4=matrix(nrow=100,ncol=1)#...matrix for storing values of X_BA# 
    for(j in 1:100)
    {
      X4[j,]=rpois(1,(((lambda_BA+(0.8)*(lambda_BA-lambda_B))*a^(G-u[j+100,])+(0.8)*X2[j,])/(1+0.8))) 
    }
    X4
    S1=sum(X1) 
    S2=sum(X2)
    S3=sum(X3)
    S4=sum(X4)
    N=200;N1=100;N2=100
    N3=100;N4=100
    S1
    S2
    S3
    S4 
    c=(1/3)*(1+a+(a^2))#....value of common probability Pie#
    c
    lambdahat_A=(S1+0.5)/(c*(N1+1))
    lambdahat_B=(S2+0.5)/(c*(N2+1))
    lambdahat_AB=(S3+0.5)/(c*(N3+1))
    lambdahat_BA=(S4+0.5)/(c*(N4+1))
    lambdahat_A 
    lambdahat_B 
    lambdahat_AB
    lambdahat_BA 
    Z=(1/2)*(log(lambdahat_B)-log(lambdahat_A)+log(lambdahat_AB)-log(lambdahat_BA))
    lambda1hat=(lambdahat_A+lambdahat_B)/2
    lambda2hat=(lambdahat_AB+lambdahat_BA)/2
    d=((1/3)*(a^4+a^2+1))-((1/9)*(a^2+a+1)^2)#....value of pie*#
    d
    sigma22hat=((1/c)^2*((lambda2hat*c)+((b^2*lambda1hat*c+(1+b)^2*lambda2hat^2*d)/(1+b)^2)))
    sigma11hat=((lambda1hat*(c+lambda1hat*d))/c^2)
    sigma12hat=(((b/(1+b))*c+lambda2hat*d)*(lambda1hat/c^2))
    sigma0hat=(sigma11hat/(2*lambda1hat^2))+(sigma22hat/(2*lambda2hat^2))+(sigma12hat/(lambda1hat*lambda2hat))
    T1[k]=(sqrt(N)*Z)/(sqrt(2)*sqrt(sigma0hat))
    D1=(sqrt(N)*(lambdahat_B-lambdahat_A))/(2*sqrt(sigma11hat))
    D1
    D2=(sqrt(N)*(lambdahat_AB-lambdahat_BA))/(2*sqrt(sigma22hat))
    D2
    r=-(sigma12hat/sqrt(sigma11hat*sigma22hat))#....value of estimated correlation coefficient#
    if((D1>0)&&(D2>0))
    {
      T2[k]=sqrt(D1^2+(D2^2)-(2*r*D1*D2))/sqrt(1-(r^2))
    }
    else if((D1<=0)&&(D2>=D1))
    {
      T2[k]=(D2-(r*D1))/sqrt(1-(r^2)) 
    }
    else
    {
      T2[k]=(D1-(r*D2))/sqrt(1-(r^2))
    }
    T2[k]
    T3[k]=max(D1,D2)
    T3[k]
  }
  #Finding the upper 5% point when the test statistic is T2
  c22=sort(T2)[9500]
  #Finding the upper 5% point when the test statistic is T3
  c32=sort(T3)[9500]
  #Checking empirical power for the same sample size by taking an alternative
  T1=NULL
  T2=NULL
  T3=NULL
  set.seed(10) 
  for(k in 1:10000)
  { 
    a=0.9   
    b=0.8#....value of beta#   
    G=2   
    lambda_A=2.5   
    lambda_B=3  
    lambda_AB=2  
    lambda_BA=1.5   
    X=matrix(nrow=200,ncol=1)   
    for(i in 1:200)   
    {     
      X[i,]=runif(1,min=0,max=1)   
    }   
    X   
    u=matrix(nrow=200,ncol=1)   
    for(i in 1:200)   
    {     
      if(X[i,]<=(1/3))     
      {      
        u[i,]=0     
      }    
      else if((X[i,]>(1/3))&&(X[i,]<=(2/3)))     
      {       
        u[i,]=1     
      }     
      else     
      {       
        u[i,]=2     
      }   
    }   
    u   
    X1=matrix(nrow=100,ncol=1)#....matrix for storing values of X_A#   
    for(i in 1:100)   
    {    
      X1[i,]=rpois(1,lambda_A*(a^(G-u[i,])))   
    }   
    X1  
    X2=matrix(nrow=100,ncol=1)#....matrix for storing values of X_B#  
    for(j in 1:100)  
    {     
      X2[j,]=rpois(1,lambda_B*(a^(G-u[j+100,])))   
    }
    X2   
    X3=matrix(nrow=100,ncol=1)#...matrix for storing values of X_AB#   
    for(i in 1:100)   
    {    
      X3[i,]=rpois(1,(((lambda_AB+(0.8)*(lambda_AB-lambda_A))*a^(G-u[i,])+(0.8)*X1[i,])/(1+0.8)))   
    }   
    X3   
    X4=matrix(nrow=100,ncol=1)#...matrix for storing values of X_BA#   
    for(j in 1:100)   
    {     
      X4[j,]=rpois(1,(((lambda_BA+(0.8)*(lambda_BA-lambda_B))*a^(G-u[j+100,])+(0.8)*X2[j,])/(1+0.8)))   
    }   
    X4   
    S1=sum(X1)   
    S2=sum(X2)   
    S3=sum(X3)   
    S4=sum(X4)  
    N=200;N1=100;N2=100   
    N3=100;N4=100   
    S1   
    S2   
    S3   
    S4   
    c=(1/3)*(1+a+(a^2))#....value of common probability Pie#   
    c   
    lambdahat_A=(S1+0.5)/(c*(N1+1))   
    lambdahat_B=(S2+0.5)/(c*(N2+1))   
    lambdahat_AB=(S3+0.5)/(c*(N3+1))   
    lambdahat_BA=(S4+0.5)/(c*(N4+1))   
    lambdahat_A  
    lambdahat_B   
    lambdahat_AB   
    lambdahat_BA   
    Z=(1/2)*(log(lambdahat_B)-log(lambdahat_A)+log(lambdahat_AB)-log(lambdahat_BA))   
    lambda1hat=(lambdahat_A+lambdahat_B)/2
    lambda2hat=(lambdahat_AB+lambdahat_BA)/2
    d=((1/3)*(a^4+a^2+1))-((1/9)*(a^2+a+1)^2)#....value of pie*#
    d   
    sigma22hat=((1/c)^2*((lambda2hat*c)+((b^2*lambda1hat*c+(1+b)^2*lambda2hat^2*d)/(1+b)^2)))
    sigma11hat=((lambda1hat*(c+lambda1hat*d))/c^2)
    sigma12hat=(((b/(1+b))*c+lambda2hat*d)*(lambda1hat/c^2))
    sigma0hat=(sigma11hat/(2*lambda1hat^2))+(sigma22hat/(2*lambda2hat^2))+(sigma12hat/(lambda1hat*lambda2hat))
    T1[k]=(sqrt(N)*Z)/(sqrt(2)*sqrt(sigma0hat))
    D1=(sqrt(N)*(lambdahat_B-lambdahat_A))/(2*sqrt(sigma11hat))
    D1   
    D2=(sqrt(N)*(lambdahat_AB-lambdahat_BA))/(2*sqrt(sigma22hat))
    D2   
    r=-(sigma12hat/sqrt(sigma11hat*sigma22hat))#....value of estimated correlation coefficient#   
    if((D1>0)&&(D2>0))   
    {     
      T2[k]=sqrt(D1^2+(D2^2)-(2*r*D1*D2))/sqrt(1-(r^2))   
    }   
    else if((D1<=0)&&(D2>=D1))   
    {    
      T2[k]=(D2-(r*D1))/sqrt(1-(r^2))   
    }  
    else   
    {    
      T2[k]=(D1-(r*D2))/sqrt(1-(r^2))   
    }   
    T2[k]   
    T3[k]=max(D1,D2)   
    T3[k] 
  }
  s=0 
  for(i in 1:10000) 
  {   
    if(T1[i]>1.645)   
    {     
      s=s+1   
    } 
  }
  p12=(s/10000)#...Empirical Power#
  c=0 
  for(i in 1:10000)
  {   
    if(T2[i]>c22)   
    {     
      c=c+1   
    } 
  }
  p22=(c/10000)#...Empirical Power#
  r=0 
  for(i in 1:10000)
  {   
    if(T3[i]>c32)   
    {     
      r=r+1  
    } 
  } 
  p32=(r/10000)#...Empirical Power#
  #For Sample Size 300 
  T1=NULL 
  T2=NULL 
  T3=NULL 
  set.seed(10) 
  for(k in 1:10000) 
  {   
    a=0.9   
    b=0.8#....value of beta#   
    G=2   
    lambda_A=2.5  
    lambda_B=2.5  
    lambda_AB=1.5   
    lambda_BA=1.5   
    X=matrix(nrow=300,ncol=1)   
    for(i in 1:300)   
    {     
      X[i,]=runif(1,min=0,max=1)   
    }   
    X   
    u=matrix(nrow=300,ncol=1)   
    for(i in 1:300)   
    {     
      if(X[i,]<=(1/3))     
      {       
        u[i,]=0     
      }    
      else if((X[i,]>(1/3))&&(X[i,]<=(2/3)))     
      {       
        u[i,]=1     
      }    
      else     
      {       
        u[i,]=2     
      }   
    }   
    u   
    X1=matrix(nrow=150,ncol=1)#....matrix for storing values of X_A#   
    for(i in 1:150)   
    {     
      X1[i,]=rpois(1,lambda_A*(a^(G-u[i,])))   
    }   
    X1   
    X2=matrix(nrow=150,ncol=1)#....matrix for storing values of X_B#   
    for(j in 1:150)   
    {     
      X2[j,]=rpois(1,lambda_B*(a^(G-u[j+150,])))   
    }   
    X2   
    X3=matrix(nrow=150,ncol=1)#...matrix for storing values of X_AB#   
    for(i in 1:150)   
    {     
      X3[i,]=rpois(1,(((lambda_AB+(0.8)*(lambda_AB-lambda_A))*a^(G-u[i,])+(0.8)*X1[i,])/(1+0.8)))   
    }   
    X3   
    X4=matrix(nrow=150,ncol=1)#...matrix for storing values of X_BA#   
    for(j in 1:150)   
    {     
      X4[j,]=rpois(1,(((lambda_BA+(0.8)*(lambda_BA-lambda_B))*a^(G-u[j+150,])+(0.8)*X2[j,])/(1+0.8)))   
    }   
    X4   
    S1=sum(X1)   
    S2=sum(X2)   
    S3=sum(X3)   
    S4=sum(X4)   
    N=300   
    N1=150   
    N2=150   
    N3=150   
    N4=150   
    S1   
    S2   
    S3   
    S4   
    c=(1/3)*(1+a+(a^2))#....value of common probability Pie#   
    c   
    lambdahat_A=(S1+0.5)/(c*(N1+1))   
    lambdahat_B=(S2+0.5)/(c*(N2+1))   
    lambdahat_AB=(S3+0.5)/(c*(N3+1))   
    lambdahat_BA=(S4+0.5)/(c*(N4+1))   
    lambdahat_A   
    lambdahat_B   
    lambdahat_AB   
    lambdahat_BA   
    Z=(1/2)*(log(lambdahat_B)-log(lambdahat_A)+log(lambdahat_AB)-log(lambdahat_BA))   
    lambda1hat=(lambdahat_A+lambdahat_B)/2
    lambda2hat=(lambdahat_AB+lambdahat_BA)/2
    d=((1/3)*(a^4+a^2+1))-((1/9)*(a^2+a+1)^2)#....value of pie*#
    d   
    sigma22hat=((1/c)^2*((lambda2hat*c)+((b^2*lambda1hat*c+(1+b)^2*lambda2hat^2*d)/(1+b)^2)))
    sigma11hat=((lambda1hat*(c+lambda1hat*d))/c^2)
    sigma12hat=(((b/(1+b))*c+lambda2hat*d)*(lambda1hat/c^2))
    sigma0hat=(sigma11hat/(2*lambda1hat^2))+(sigma22hat/(2*lambda2hat^2))+(sigma12hat/(lambda1hat*lambda2hat))
    T1[k]=(sqrt(N)*Z)/(sqrt(2)*sqrt(sigma0hat))
    D1=(sqrt(N)*(lambdahat_B-lambdahat_A))/(2*sqrt(sigma11hat))
    D2=(sqrt(N)*(lambdahat_AB-lambdahat_BA))/(2*sqrt(sigma22hat))
    r=-(sigma12hat/sqrt(sigma11hat*sigma22hat))#....value of estimated correlation coefficient#   
    if((D1>0)&&(D2>0))   
    {     
      T2[k]=sqrt(D1^2+(D2^2)-(2*r*D1*D2))/sqrt(1-(r^2))   
    }   
    else if((D1<=0)&&(D2>=D1))   
    {     
      T2[k]=(D2-(r*D1))/sqrt(1-(r^2))   
    }   
    else   
    {     
      T2[k]=(D1-(r*D2))/sqrt(1-(r^2))   
    }   
    T3[k]=max(D1,D2) 
  }
  #Finding the upper 5% point when the test statistic is T2
  c23=sort(T2)[9500]
  #Finding the upper 5% point when the test statistic is T3
  c33=sort(T3)[9500]
  #Checking empirical power for the same sample size by taking an alternative 
  T1=NULL 
  T2=NULL 
  T3=NULL 
  set.seed(10) 
  for(k in 1:10000) 
  {   
    a=0.9   
    b=0.8#....value of beta#   
    G=2   
    lambda_A=2.5  
    lambda_B=3  
    lambda_AB=2  
    lambda_BA=1.5  
    X=matrix(nrow=300,ncol=1)  
    for(i in 1:300)   
    {    
      X[i,]=runif(1,min=0,max=1)   
    }   
    X   
    u=matrix(nrow=300,ncol=1)   
    for(i in 1:300)   
    {     
      if(X[i,]<=(1/3))    
      {
        u[i,]=0     
      }    
      else if((X[i,]>(1/3))&&(X[i,]<=(2/3)))     
      {       
        u[i,]=1     
      }     
      else     
      {       
        u[i,]=2     
      }   
    }   
    u   
    X1=matrix(nrow=150,ncol=1)#....matrix for storing values of X_A#   
    for(i in 1:150)   
    {     
      X1[i,]=rpois(1,lambda_A*(a^(G-u[i,])))   
    }   
    X1   
    X2=matrix(nrow=150,ncol=1)#....matrix for storing values of X_B#   
    for(j in 1:150)   
    {     
      X2[j,]=rpois(1,lambda_B*(a^(G-u[j+150,])))   
    }   
    X2   
    X3=matrix(nrow=150,ncol=1)#...matrix for storing values of X_AB#   
    for(i in 1:150)   
    {     
      X3[i,]=rpois(1,(((lambda_AB+(0.8)*(lambda_AB-lambda_A))*a^(G-u[i,])+(0.8)*X1[i,])/(1+0.8)))   
    }   
    X3   
    X4=matrix(nrow=150,ncol=1)#...matrix for storing values of X_BA#   
    for(j in 1:150)   
    {    
      X4[j,]=rpois(1,(((lambda_BA+(0.8)*(lambda_BA-lambda_B))*a^(G-u[j+150,])+(0.8)*X2[j,])/(1+0.8)))   
    }   
    X4  
    S1=sum(X1)   
    S2=sum(X2)   
    S3=sum(X3)   
    S4=sum(X4)   
    N=300   
    N1=150   
    N2=150   
    N3=150   
    N4=150   
    S1   
    S2   
    S3   
    S4   
    c=(1/3)*(1+a+(a^2))#....value of common probability Pie#   
    c   
    lambdahat_A=(S1+0.5)/(c*(N1+1))   
    lambdahat_B=(S2+0.5)/(c*(N2+1))   
    lambdahat_AB=(S3+0.5)/(c*(N3+1))   
    lambdahat_BA=(S4+0.5)/(c*(N4+1))   
    lambdahat_A   
    lambdahat_B   
    lambdahat_AB   
    lambdahat_BA   
    Z=(1/2)*(log(lambdahat_B)-log(lambdahat_A)+log(lambdahat_AB)-log(lambdahat_BA))   
    lambda1hat=(lambdahat_A+lambdahat_B)/2
    lambda2hat=(lambdahat_AB+lambdahat_BA)/2
    d=((1/3)*(a^4+a^2+1))-((1/9)*(a^2+a+1)^2)#....value of pie*#
    d   
    sigma22hat=((1/c)^2*((lambda2hat*c)+((b^2*lambda1hat*c+(1+b)^2*lambda2hat^2*d)/(1+b)^2)))
    sigma11hat=((lambda1hat*(c+lambda1hat*d))/c^2)
    sigma12hat=(((b/(1+b))*c+lambda2hat*d)*(lambda1hat/c^2))
    sigma0hat=(sigma11hat/(2*lambda1hat^2))+(sigma22hat/(2*lambda2hat^2))+(sigma12hat/(lambda1hat*lambda2hat))
    T1[k]=(sqrt(N)*Z)/(sqrt(2)*sqrt(sigma0hat))
    D1=(sqrt(N)*(lambdahat_B-lambdahat_A))/(2*sqrt(sigma11hat))
    D2=(sqrt(N)*(lambdahat_AB-lambdahat_BA))/(2*sqrt(sigma22hat))
    r=-(sigma12hat/sqrt(sigma11hat*sigma22hat))#....value of estimated correlation coefficient#   
    if((D1>0)&&(D2>0))   
    {     
      T2[k]=sqrt(D1^2+(D2^2)-(2*r*D1*D2))/sqrt(1-(r^2))   
    }   
    else if((D1<=0)&&(D2>=D1))   
    {    
      T2[k]=(D2-(r*D1))/sqrt(1-(r^2))   
    }   
    else   
    {     
      T2[k]=(D1-(r*D2))/sqrt(1-(r^2))   
    }   
    T3[k]=max(D1,D2) 
  } 
  s=0 
  for(i in 1:10000) 
  {  
    if(T1[i]>1.645)   
    {     
      s=s+1   
    } 
  }
  p13=(s/10000)#...Empirical Power#
  c=0
  for(i in 1:10000)
  {  
    if(T2[i]>c23)   
    {     
      c=c+1   
    }
  }
  p23=(c/10000)#...Empirical Power#
  r=0 
  for(i in 1:10000) 
  {  
    if(T3[i]>c33)   
    {    
      r=r+1  
    } 
  } 
  p33=(r/10000)#...Empirical Power#
  #Plotting of power Curves
  power=matrix(c(p11,p12,p13,p21,p22,p23,p31,p32,p33),nrow=3,byrow=T)
  p1=power[1,] 
  p2=power[2,] 
  p3=power[3,] 
  power 
  n=c(150,200,300)
  #Plotting on the same graph 
  plot(n,p1,type='o',main="Power Curves",xlab='Sample Size',ylab='Power',col='blue',pch='o',lty=1,ylim=c(0,1.4))
  points(n,p2,col="red",pch="*")
  lines(n,p2,col="red",lty=2)
  points(n,p3,col="brown",pch="+")
  lines(n,p3,col="brown",lty=3)
  legend("topright",c("T1","T2","T3"),col=c("blue","red","brown"),bty="n",lty=280:300)
  