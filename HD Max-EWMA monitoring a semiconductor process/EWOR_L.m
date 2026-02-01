tic;
clc;
clear;
load('Phase_I_data.mat');
load('Phase_II_data.mat');


lambda=0.2;
UCL_OR=4.125;
n=5;
p=248;
t=20;
m=1442;

OR=zeros(1,p);
x_m=Phase_I_data(1:m,:);
d_0_OR= zeros(size(x_m,1),1);
d_0_OR=pdist2(x_m,OR);

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
Z(1)=2; 
for tk=1:t  
d_1_OR(:,1,tk) =pdist2(Y_n(:,:,tk),OR);
  
N1(:,:,tk) = sort([d_0_OR;d_1_OR(:,:,tk)]);
 
T1k_OR(tk) = Wilcoxon(N1(:,:,tk), d_1_OR(:,:,tk));
  
   
 T2k_OR(tk) = AB(N1(:,:,tk), d_1_OR(:,:,tk));
  
 LK_OR(tk)= T1k_OR(tk)+ T2k_OR(tk);

  Z(tk+1)=lambda*LK_OR(tk)+(1-lambda)*Z(tk);
end 
 a=UCL_OR*ones(21,1);
 figure(1)
 box on 
 hold on
 i=[0:1:20];
 plot(i,a,'r-', 'LineWidth', 2);
 plot(i,Z,'bo-', 'LineWidth', 1);
  grid on;
xlabel('$t$','interpreter','latex');
ylabel('$Z_t$','interpreter','latex');
text(20.5,UCL_OR,'$UCL=4.125$','interpreter','latex')
set(gca,'xtick',0:5:20,'ytick',0:2:8);
axis([0 26 0 8]);
hold off

toc;




