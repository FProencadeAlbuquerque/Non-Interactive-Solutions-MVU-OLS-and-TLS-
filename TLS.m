Pcarga='90';

arquivo=strcat('Dados230kVCarga',Pcarga);

load(arquivo)

Rex=5;

Xex=48.8;

Bex=3.371*10^(-4);
 
m=1;

for j=1:m
    
    Ve=Vemod(1000)*exp(i*Vefase(1000)*pi/180);
        
    Ie=Iemod(1000)*exp(i*Iefase(1000)*pi/180);
        
    Vs=Vsmod(1000)*exp(i*Vsfase(1000)*pi/180);
        
    Is=Ismod(1000)*exp(i*Isfase(1000)*pi/180);

    sigma=0.002;

    n=100000;

    b1=(real(Ie+Is))/(-imag(Ve+Vs));
    
    rng(42+j)
    
    w1=sigma*abs(real(Ve-Vs))*normrnd(0,1,[1,n]);
    
    rng(42+j+1)
   
    w2=sigma*abs(imag(Ve-Vs))*normrnd(0,1,[1,n]);
    
    rng(42+j+2)
    
    w3=sigma*abs(real(Ve-Vs))*normrnd(0,1,[1,n]);
    
    rng(42+j+3)
   
    w4=sigma*abs(imag(Ve-Vs))*normrnd(0,1,[1,n]);
    
    rng(42+j+4)
    
    w5=sigma*normrnd(0,1,[1,n]);
    
    rng(42+j+5)
    
    w6=sigma*normrnd(0,1,[1,n]);
    
    rng(42+j+6)
   
    w7=sigma*normrnd(0,1,[1,n]);
    
    rng(42+j+7)
    
    w8=sigma*normrnd(0,1,[1,n]);

    DeltaH11=abs(real(Ie)+b1*imag(Ve))*w5;
    DeltaH12=abs(-imag(Ie)+b1*real(Ve))*w6;
    DeltaH21=abs(imag(Ie)-b1*real(Ve))*w6;
    DeltaH22=abs(real(Ie)+b1*imag(Ve))*w5;
    DeltaH31=abs(-b1*imag(Vs)-real(Is))*w7;
    DeltaH32=abs(-b1*real(Vs)+imag(Is))*w8;
    DeltaH41=abs(b1*real(Vs)-imag(Is))*w8;
    DeltaH42=abs(-b1*imag(Vs)-real(Is))*w7;
        
    B1=2*b1;
    
    rng(42+j+8)
    
    r1=normrnd(0,sigma*B1,[1,n]);

    B=zeros(1,n);
    R=zeros(1,n);
    Z=zeros(1,n);

    for k=1:n
    
        B(k)=((2*real(Ie+Is))/(-imag(Ve+Vs)))+r1(k);
    
        Ytil=[w1(k);w2(k);w1(k);w2(k)];
    
        Y=[real(Ve-Vs); imag(Ve-Vs); real(Ve-Vs);imag(Ve-Vs)]+Ytil;
    
        Hbar=[(imag(Ve)*b1+real(Ie)) (-imag(Ie)+real(Ve)*b1); (-real(Ve)*b1+imag(Ie)) (imag(Ve)*b1+real(Ie)); (-real(Is)-imag(Vs)*b1) (-real(Vs)*b1+imag(Is)); (-imag(Is)+real(Vs)*b1) (-real(Is)-imag(Vs)*b1)];
    
        Htil=[DeltaH11(k) DeltaH12(k); DeltaH21(k) DeltaH22(k); DeltaH31(k) DeltaH32(k); DeltaH41(k) DeltaH42(k)];

        H=Hbar+Htil;
    
        L=[H Y];
        [U,S,V]=svd(L);
        VXY=V(1:2,3);
        VYY=V(3,3);
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%% Total least Square %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
         X=-VXY*VYY^(-1);
         R(k)=X(1,1);
         Z(k)=X(2,1);
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        Bcal(j)=mean(B);
        DesvB(j)=abs((Bex-Bcal(j))/(Bex))*100;
        stdB(j)=std(B);
        Xcal(j)=mean(Z);
        DesvX(j)=abs((Xex-Xcal(j))/Xex)*100;
        stdX(j)=std(Z);
        Rcal(j)=mean(R);
        DesvR(j)=abs((Rex-Rcal(j))/Rex)*100;
        stdR(j)=std(R);
        
 end
    
 fprintf ('Desvio Relativo de R= %.8f Devio padrão de R= %.8f \n', mean(DesvR),mean(stdR))
 fprintf ('Desvio Relativo de X= %.8f Devio padrão de X= %.8f \n', mean(DesvX),mean(stdX))
 fprintf ('Desvio Relativo de B= %.8f Devio padrão de B= %.8f *10^(-5) \n', mean(DesvB),mean(stdB)*10^5)   
    


