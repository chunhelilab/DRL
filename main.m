%This is an implementation of the DRL(dimension reduction of landscape)
%method on MESC model.
clear
cycle_index=3000;  %% The number of random initial conditions to the ODEs to be solved
par=[3.8 0.4 2.9 4 1];   %% The parameters of the ODE
signal=[3,3,3];  %% Signal
d=0.04;  %% The diffusion coefficient 
N=12; %% The dimension of the system
tic() %% Time
%% Solve the ODEs, calculate the paths and actions;
[xx,sigma,n,ycell,action]=Solver(cycle_index,par,signal,d);

index=size(n,1);  %% The number of the stable states
alpha=zeros(index,1);  %% The weight of the stable states
sigma0=cell(index,1);  %% The covariance of the Gaussian density function
mu=zeros(index,N);  %% The mean value of the Gaussian density function

for i=1:index
   %The mean value of each stable state
   mu(i,:)=xx(n(i,1),:); 
   %The covariance of each stable state
   sigma0{i}=reshape(sigma(n(i,1),:),N,N)';  
   %The weight of each stable state
   alpha(i)=n(i,2)/sum(n(:,2)); 
end

%% DRL
%Calculate the mean value 
Mu=0;
for i=1:index
    Mu=Mu+alpha(i)*mu(i,:);
end
%Calculate the covariance
Sigma=-Mu'*Mu; 
for i=1:index
    Sigma=Sigma+alpha(i)*(sigma0{i}+mu(i,:)'*mu(i,:));
end

%Calculate the eigenvalues and eigenvectors of the covariance
[V,D] = eigs(Sigma,2);

if sign(V(:,1)'*[1,1,1,1,1,1,1,1,1,1,1,1]')<0
    V(:,1)=-V(:,1);
end
if sign(V(:,2)'*[1,1,1,1,1,1,1,1,1,1,1,1]')<0
    V(:,2)=-V(:,2);
end

%%%Calculate the covariance and mean value after dimension reduction
sigma0_pca=cell(index,1);
mu_pca=zeros(index,2);
for i=1:index
   mu_pca(i,:)=V'*mu(i,:)';
   sigma0_pca{i}=V'*sigma0{i}*V;
end

%% plot the landscape
y_max=[24,12]; %% Range of the landscape
y_min=[0,0];
step=(y_max-y_min)/100; %% Length of the step
[a1,a2]=meshgrid(y_min(1):step(1):y_max(1),y_min(2):step(2):y_max(2)); %% Grid
[s1,s2]=size(a1);
P=zeros(s1,s2);
z=zeros(s1,s2);
for kk=1:index
    sig=sigma0_pca{kk};
    x_wen=mu_pca(kk,:);
    for i=1:s1
        for j=1:s2
            z(i,j)=multivariate_normal_distribution([a1(i,j);a2(i,j)],x_wen',sig,2);  %% Normal distribution
        end
    end

    P=P+z*alpha(kk);
end
P=P/sum(sum(P));
surf(a1,a2,-log(max(P,10^-100)));   %% Plot landscape
shading interp
xlabel('PC1')
ylabel('PC2')
zlabel('U')
axis([0 22 0 11 0 240])

for i=1:size(n,1)
    A(i)=floor((mu_pca(i,1)-y_min(1))/step(1))+1;
    B(i)=floor((mu_pca(i,2)-y_min(2))/step(2))+1;
end
 hold on

%Plot the grid
for i=1:floor(size(a1,1)/4)
    plot3(a1(4*i-1,:),a2(4*i-1,:),-log(max(P(4*i-1,:),10^-100)),'Color',[0.4 0.4 0.4],'LineWidth',0.01);
end
for i=1:floor(size(a1,2)/4)
    plot3(a1(:,4*i-1),a2(:,4*i-1),-log(max(P(:,4*i-1),10^-100)),'Color',[0.4 0.4 0.4],'LineWidth',0.01);
end

%% Plot the paths
%Calculate the paths after dimension reduction
y12=V'*ycell{1,2};
y21=V'*ycell{2,1};
view([-25,75])
hold on
z3path=griddata(a1,a2,-log(max(P,10^-100)),y12(1,:),y12(2,:));
plot3(y12(1,:),y12(2,:),z3path+3,'w','LineWidth',2);

z3path=griddata(a1,a2,-log(max(P,10^-100)),y21(1,:),y21(2,:));
plot3(y21(1,:),y21(2,:),z3path+3,'Color',[0.85,0.43,0.83],'LineWidth',2);

view([-25 75])
set(gcf,'outerposition', [100 100 800 650]);
toc()