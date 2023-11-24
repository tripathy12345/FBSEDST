function [a3, alfa, fre, y]=FBSE_S_transform(x,fs)
f=x;
N=length(f);
nb=(1:N);
% f=f';
MM=N;
if exist('alfa') == 0
    x=2;
    alfa=zeros(1,MM);
    for i=1:MM
        ex=1;
        while abs(ex)>.00001
            ex=-besselj(0,x)/besselj(1,x);
            x=x-ex;
        end
        alfa(i)=x;
        fprintf('Root # %g  = %8.5f ex = %9.6f \n',i,x,ex)
        x=x+pi;
    end
end
a=N;
for m1=1:MM
D(m1,:)=besselj(0,alfa(m1)/a*nb);
end
for m1=1:MM
    a3(m1)=(2/(a^2*(besselj(1,alfa(m1))).^2))*sum(nb.*f.*D(m1,:));
end
% plot(fre,a3)
invfk=[1./alfa]';
W=2*pi*0.5*repmat(alfa,N,1).*repmat(invfk,1,N);
G=exp((-W.^2)/2); %Gaussian in FBSE domain
% Compute Toeplitz matrix with the shifted fbse domain
HW=toeplitz(a3(1:N)',a3);
HW=[HW(1:N,:)];
zz=HW.*G;
% % Compute FBSE-Stockwell Transform
y=abs(zz*D); %%%%FBSE-DST time-frequency matrix
fre=(alfa*fs)/(2*pi*N);
end