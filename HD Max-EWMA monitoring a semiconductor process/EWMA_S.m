tic;
clc;
clear;

rho=0.5;
p = 248; 
lambda = 0.2; 

m = 1442; 
n =5;
t=20;
UCL=68.31;
load('Phase_I_data.mat');
load('Phase_II_data.mat');
x_m=Phase_I_data(1:m,:);

[theta0, A0] = computeAEMMedian(x_m);

Y_n = zeros(n,p,t);
Y_n(:,:,1)=Phase_II_data(1:5,:);
for tk=2:t
     Y_n(:,:,tk)=Phase_II_data((tk-1)*5+1:tk*5,:);
end


Z = zeros(p, t);
V=zeros(n,p);
Q = zeros(1, t);
for tk=1:t
    for tkk=1:n
        V(tkk,:)=(A0 * ((Y_n(tkk,:,tk))'-theta0)./ norm(A0 * ((Y_n(tkk,:,tk))'-theta0)));
    end
   Z(:,tk+1) =(1-lambda).*Z(:,tk)+lambda.* mean(V)';
   Q(tk+1)=(2-lambda)/(lambda)*p*Z(:,tk+1)'*Z(:,tk+1);

end

 a=UCL*ones(21,1);
 figure(1)
 box on 
 hold on
 i=[0:1:20];
 plot(i,a,'r-', 'LineWidth', 2);
 plot(i,Q,'bo-', 'LineWidth', 1);
  grid on;
 
xlabel('$t$','interpreter','latex');
ylabel('$Z_t$','interpreter','latex');
text(20.5,UCL,'$UCL=68.31$','interpreter','latex')
set(gca,'xtick',0:5:20,'ytick',0:20:100);
axis([0 26 0 100]);
hold off

toc;

