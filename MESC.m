function f=MESC(t,x,p,signal)
LIF=signal(1);
CH=signal(2);
PD=signal(3);
g0=0.01;
a(1:18)=p(1);
b(1:8)=p(2);
S=p(3);
n=p(4);
k=p(5);

% mekerk=x(1);
% oct4=x(2);
% sox2=x(3);
% nanog=x(4);
% esrrb=x(5);
% klf2=x(6);
% tfcp2l1=x(7);
% klf4=x(8);
% gbx2=x(9);
% tcf3=x(10);
% sall4=x(11);
% stat3=x(12);

f(1,:)=b(1)*Hr(PD,S,n)+g0-k*x(1,:);
f(2,:)=a(1)*Ha(x(3,:),S,n)+a(2)*Ha(x(6,:),S,n)+b(2)*Hr(x(5,:),S,n)+g0-k*x(2,:);
f(3,:)=a(3)*Ha(x(4,:),S,n)+a(4)*Ha(x(11,:),S,n)+g0-k*x(3,:);
f(4,:)=a(5)*Ha(x(2,:),S,n)+a(6)*Ha(x(6,:),S,n)+b(3)*Hr(x(1,:),S,n)+g0-k*x(4,:);
f(5,:)=a(7)*Ha(x(7,:),S,n)+a(8)*Ha(x(4,:),S,n)+b(4)*Hr(x(10,:),S,n)+g0-k*x(5,:);
f(6,:)=a(9)*Ha(x(8,:),S,n)+a(10)*Ha(x(11,:),S,n)+g0-k*x(6,:);
f(7,:)=a(11)*Ha(x(12,:),S,n)+a(12)*Ha(x(8,:),S,n)+a(13)*Ha(x(5,:),S,n)+b(5)*Hr(x(10,:),S,n)+b(6)*Hr(x(2,:),S,n)+g0-k*x(7,:);
f(8,:)=a(14)*Ha(x(12,:),S,n)+a(15)*Ha(x(9,:),S,n)+g0-k*x(8,:);
f(9,:)=a(16)*Ha(x(12,:),S,n)+g0-k*x(9,:);
f(10,:)=b(7)*Hr(CH,S,n)+b(8)*Hr(x(1,:),S,n)+g0-k*x(10,:);
f(11,:)=a(17)*Ha(x(7,:),S,n)+g0-k*x(11,:);
f(12,:)=a(18)*Ha(LIF,S,n)+g0-k*x(12,:);
end

function H=Ha(X,S,n)
H=X.^n./(S.^n+X.^n);
end
function H=Hr(X,S,n)
H=S.^n./(X.^n+S.^n);
end