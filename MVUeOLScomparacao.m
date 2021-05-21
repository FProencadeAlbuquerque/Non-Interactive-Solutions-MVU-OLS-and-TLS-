%%%%%%%%%%%% Tratamento dos sinais de tempo contínuo obtidos no simulink %%%%%%%%%%%%%%%%%

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

    sigma=0.12;

    modr2=sigma*abs(real(Ve-Vs));
    modr3=sigma*abs(imag(Ve-Vs));
    modr4=sigma*abs(real(Ve-Vs));
    modr5=sigma*abs(imag(Ve-Vs));

    n=100000;

    B=zeros(1,n);
    R1=zeros(1,n);
    Z1=zeros(1,n);

        for k=1:n
        
        %Montando a matriz do problema%
    
        b1=(real(Ie+Is))/(-imag(Ve+Vs));
    
        H=[(imag(Ve)*b1+real(Ie)) (-imag(Ie)+real(Ve)*b1); (-real(Ve)*b1+imag(Ie)) (imag(Ve)*b1+real(Ie)); (-real(Is)-imag(Vs)*b1) (-real(Vs)*b1+imag(Is)); (-imag(Is)+real(Vs)*b1) (-real(Is)-imag(Vs)*b1)];
    
        C=[modr2^2 0 0 0;0 modr3^2 0 0;0 0 modr4^2 0;0 0 0 modr5^2];
        
        CL=[modr2 0 0 0;0 modr3 0 0;0 0 modr4 0;0 0 0 modr5];
        
        wb=Bex*sigma;
        
        rng(42+k)
        
        r1=normrnd(0,wb,[1,n]);
        
        rng(42+k+1)
    
        e1=CL*randn(4,1);
        
        b=((2*real(Ie+Is))/(-imag(Ve+Vs)))+r1(k);
    
        Y=[real(Ve-Vs);imag(Ve-Vs);real(Ve-Vs);imag(Ve-Vs)]+e1;
        
        %%%%%%%%%%% Solucao MVU %%%%%%%%%%%%%%%%%%%%%%%
        
        X1=(transpose(H)*C^(-1)*H)^(-1)*transpose(H)*C^(-1)*Y;
        
        %%%%% Método dos Mínimos quadrados ordinários %%%%
        
        X2=(transpose(H)*H)^(-1)*transpose(H)*Y;
        
        R1(k)=X1(1,1);
        Z1(k)=X1(2,1);
        B(k)=b;
        R2(k)=X2(1,1);
        Z2(k)=X2(2,1);
        
    end
    
    %%%%%%%%%% valores para o MVU
        
    Bcal(j)=mean(B);
    DesvB(j)=abs((Bex-Bcal(j))/(Bex))*100;
    stdB(j)=std(B);
    Xcal(j)=mean(Z1);
    DesvX(j)=abs((Xex-Xcal(j))/Xex)*100;
    stdX(j)=std(Z1);
    Rcal(j)=mean(R1);
    DesvR(j)=abs((Rex-Rcal(j))/Rex)*100;
    stdR(j)=std(R1);
    
    %%%%%%%% valores para o OLS
    
    Xcal2(j)=mean(Z2);
    DesvX2(j)=abs((Xex-Xcal2(j))/Xex)*100;
    stdX2(j)=std(Z2);
    Rcal2(j)=mean(R2);
    DesvR2(j)=abs((Rex-Rcal2(j))/Rex)*100;
    stdR2(j)=std(R2);
        
        
end

%%
fprintf ('Valores Obtidos para o caso MVU \n')
fprintf ('Desvio Relativo de R= %.8f Devio padrão de R= %.8f \n', mean(DesvR),mean(stdR))
fprintf ('Desvio Relativo de X= %.8f Devio padrão de X= %.8f \n', mean(DesvX),mean(stdX))
fprintf ('Desvio Relativo de B= %.4f Devio padrão de B= %.8f *10^(-5) \n', mean(DesvB),mean(stdB)*10^5)
fprintf ('Valores Obtidos para o caso OLS \n')
fprintf ('Desvio Relativo de R= %.8f Devio padrão de R= %.8f \n', mean(DesvR2),mean(stdR2))
fprintf ('Desvio Relativo de X= %.8f Devio padrão de X= %.8f \n', mean(DesvX2),mean(stdX2))    
