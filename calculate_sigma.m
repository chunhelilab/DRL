function [xx]=calculate_sigma(x,par,signal,kk,d)

%par is parameters
%d is diffusion coefficient
%the xx is the vector of sigma which is Row Major Order


%%lead in the parameters
LIF=signal(1);
CH=signal(2);
PD=signal(3);
g0=0.01;
a=par(1);
b=par(2);
S=par(3);
n=par(4);
k=par(5);

syms stat3 gbx2 klf4 klf2 sall4 tcf3 tfcp2l1 oct4 mekerk esrrb nanog sox2; 
 Ajac=jacobian([b*Hr(PD,S,n)+g0-k*mekerk;
    a*Ha(sox2,S,n)+a*Ha(klf2,S,n)+b*Hr(esrrb,S,n)+g0-k*oct4;
    a*Ha(nanog,S,n)+a*Ha(sall4,S,n)+g0-k*sox2;
    a*Ha(oct4,S,n)+a*Ha(klf2,S,n)+b*Hr(mekerk,S,n)+g0-k*nanog;
    a*Ha(tfcp2l1,S,n)+a*Ha(nanog,S,n)+b*Hr(tcf3,S,n)+g0-k*esrrb;
    a*Ha(klf4,S,n)+a*Ha(sall4,S,n)+g0-k*klf2;
    a*Ha(stat3,S,n)+a*Ha(klf4,S,n)+a*Ha(esrrb,S,n)+b*Hr(tcf3,S,n)+b*Hr(oct4,S,n)+g0-k*tfcp2l1;
    a*Ha(stat3,S,n)+a*Ha(gbx2,S,n)+g0-k*klf4;
    a*Ha(stat3,S,n)+g0-k*gbx2;
    b*Hr(CH,S,n)+b*Hr(mekerk,S,n)+g0-k*tcf3;
    a*Ha(tfcp2l1,S,n)+g0-k*sall4;
    a*Ha(LIF,S,n)+g0-k*stat3],[mekerk,oct4,sox2,nanog,esrrb,klf2,tfcp2l1,klf4,gbx2,tcf3,sall4,stat3]);
Ajac=subs(Ajac,{'mekerk','oct4','sox2','nanog','esrrb','klf2','tfcp2l1','klf4','gbx2','tcf3','sall4','stat3'},{x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12)});
Ajac=double(Ajac);

% A*sigma+sigma*A'+2D


P=zeros(kk^2,kk^2);  %coefficient matrix

%%the initial of coeffiicient matrix
for i=0:(kk-1)
    P(i*kk+1:i*kk+kk,i*kk+1:i*kk+kk)=P(i*kk+1:i*kk+kk,i*kk+1:i*kk+kk)+Ajac;
end

for m=0:kk-1
    for i=1:kk
        for j=1:kk
            P(m*kk+i,(j-1)*kk+i)=P(m*kk+i,(j-1)*kk+i)+Ajac(m+1,j);
        end
    end
end

B=zeros(kk^2,1);
for i=1:kk
    B((i-1)*kk+i)=-2*d;
end


xx=P\B;

end
function H=Ha(X,S,n)
H=X.^n./(S.^n+X.^n);
end
function H=Hr(X,S,n)
H=S.^n./(X.^n+S.^n);
end