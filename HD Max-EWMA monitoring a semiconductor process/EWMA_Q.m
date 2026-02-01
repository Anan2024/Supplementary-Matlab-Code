tic;
clear;
warning off;

lambda=0.2;
LCL=0.39;

m=1442;
n=5;
p=248;

load('Phase_I_data.mat');
load('Phase_II_data.mat'); 

t=20;

N=m+n;

x_m=Phase_I_data(1:m,:);
Y_n = zeros(n,p,t);
Y_n(:,:,1)=Phase_II_data(1:5,:);
for tk=2:t
     Y_n(:,:,tk)=Phase_II_data((tk-1)*5+1:tk*5,:);
end


A=cov(x_m);
B=mean(x_m);
d_0_OR=pdist2(x_m,B,'mahalanobis').^2;
RS=1./(1+d_0_OR);

d_1_OR = zeros(n,1,t);
Z=zeros(1,t+1);
Z(1)=0.5;
for tk=1:t   
   C=zeros(1,n);
     d_1_OR(:,1,tk) =pdist2(Y_n(:,:,tk),B,'mahalanobis',A).^2;
     RT=1./(1+d_1_OR(:,1,tk));
    for zz=1:n 
     C(zz)=sum(RS<=RT(zz))/m;
    end
     Q=sum(C)/n;
     Z(tk+1)=lambda*Q+(1-lambda)*Z(tk);
     

end 

a=LCL*ones(21,1);
 figure(1)
 box on 
 hold on
 i=[0:1:20];
 plot(i,a,'r-', 'LineWidth', 2);
 plot(i,Z,'bo-', 'LineWidth', 1);
 grid on;
 
xlabel('$t$','interpreter','latex');
ylabel('$Q_t$','interpreter','latex');
text(20.5,LCL,'$LCL=0.39$','interpreter','latex')
set(gca,'xtick',0:5:20,'ytick',0:0.2:0.6);
axis([0 26 0 0.6]);
hold off


toc;




