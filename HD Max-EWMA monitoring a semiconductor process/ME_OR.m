tic;
clc;
clear;
load('Phase_I_data.mat');
load('Phase_II_data.mat');

lambda=0.2;
UCL_OR=2.798;
n=5;
p=248;

t=20;
m=1442;

N=m+n;
mu_W=n*(N+1)/2;
sigma_W=m*n*(N+1)/12;

if mod(N,2)==0;
    mu_AB=n*N/4;
    sigma_AB=m*n*(N^2-4)/(48*(N-1));
else
    mod(N,2)==1;
    mu_AB=n*(N^2-1)/(4*N);
    sigma_AB=m*n*(N+1)*(N^2+3)/(48*N^2);
end


OR=zeros(1,p);
OR_0=zeros(m,p);
OR_1=zeros(n,p);
x_m=Phase_I_data(1:m,:);
d_0_OR= zeros(size(x_m,1),1);


% d_0_OR=pdist2(x_m,OR,'minkowski',1);
% d_0_OR=pdist2(x_m,OR);
  d_0_OR=pdist2(x_m,OR,'chebychev');


Y_n = zeros(n,p,t);
Y_n(:,:,1)=Phase_II_data(1:5,:);
for tk=2:t
     Y_n(:,:,tk)=Phase_II_data((tk-1)*5+1:tk*5,:);
end
 d_1_OR = zeros(n,1,t);
 N1= zeros(m+n,1,t);

T1k_OR = zeros(1,t);

T2k_OR = zeros(1,t);

LK_OR = zeros(1,t);
Z=zeros(1,t+1);
Z1(1)=1; 
Z2(1)=1; 
Z(1)=1; 
for tk=1:t  

  % d_1_OR(:,1,tk) =pdist2(Y_n(:,:,tk),OR,'minkowski',1);
  % d_1_OR(:,1,tk) =pdist2(Y_n(:,:,tk),OR,'minkowski',2);
   d_1_OR(:,1,tk) =pdist2(Y_n(:,:,tk),OR,'chebychev');
  
N1(:,:,tk) = sort([d_0_OR;d_1_OR(:,:,tk)]);
 
[~,N1_order(:,:,tk)]=ismember(d_1_OR(:,:,tk),N1(:,:,tk));
    T1k_OR(tk) = (sum(N1_order(:,:,tk))-mu_W)^2/sigma_W;
     N1_ORDER=sum(abs( N1_order(:,:,tk)-(N+1)/2));
    T2k_OR(tk)=(N1_ORDER-mu_AB)^2/sigma_AB;

     Z1(tk+1)=lambda*T1k_OR(tk)+(1-lambda)*Z1(tk);
     Z2(tk+1)=lambda*T2k_OR(tk)+(1-lambda)*Z2(tk);
     Z(tk+1)=max(Z1(tk+1),Z2(tk+1));
end 
 a=UCL_OR*ones(21,1);
 figure(1)
 box on 
 hold on
 i=[0:1:20];
 plot(i,a,'r-', 'LineWidth', 2);
 plot(i,Z,'bo-', 'LineWidth', 1);
 grid on;
 
% xlabel('$t$','interpreter','latex');
% ylabel('$M^{1}_t$','interpreter','latex');
% text(20.5,UCL_OR,'$UCL=2.798$','interpreter','latex')
% set(gca,'xtick',0:5:20,'ytick',0:1:4);
% axis([0 26 0 4]);
% hold off

% xlabel('$t$','interpreter','latex');
% ylabel('$M^{2}_t$','interpreter','latex');
% text(20.5,UCL_OR,'$UCL=2.798$','interpreter','latex')
% set(gca,'xtick',0:5:20,'ytick',0:1:4);
% axis([0 26 0 4]);
% hold off

 xlabel('$t$','interpreter','latex');
ylabel('$M^{+\infty}_t$','interpreter','latex');
text(20.5,UCL_OR,'$UCL=2.798$','interpreter','latex')
set(gca,'xtick',0:5:20,'ytick',0:1:4);
axis([0 26 0 4]);
hold off


toc;




