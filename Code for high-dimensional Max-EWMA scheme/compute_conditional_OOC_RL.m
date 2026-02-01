function RL = compute_conditional_OOC_RL(H, lambda, p, m, n, T_0, d,delta_1,delta_2)

OR=zeros(1,p);
RL=T_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mean and std of Wilcoxon statistic
N=m+n;
mu_W=n*(N+1)/2;
sigma_W=m*n*(N+1)/12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%mean and std of AB statistic
if mod(N,2)==0;
    mu_AB=n*N/4;
    sigma_AB=m*n*(N^2-4)/(48*(N-1));
else
    mod(N,2)==1;
    mu_AB=n*(N^2-1)/(4*N);
    sigma_AB=m*n*(N+1)*(N^2+3)/(48*N^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho=0.5;
mu_0=zeros(1,p);
sigma_0=ones(p,p)*rho;
sigma_0(logical(eye(size(sigma_0))))=1;
x_m = mvnrnd(mu_0,sigma_0,m);

d_0_OR= zeros(size(x_m,1),1);
d_0_OR=pdist2(x_m,OR,'minkowski',d);

Y_n = zeros(n,p,T_0);
mu_1=zeros(1,p)+delta_1;
sigma_1=ones(p,p)*rho;
sigma_1(logical(eye(size(sigma_1))))=1;
sigma_1=sigma_1.*delta_2^2;

for tk=1:T_0
    Y_n(:,:,tk)=mvnrnd(mu_1,sigma_1,n);
end
d_1_OR = zeros(n,1,T_0);
N1= zeros(m+n,1,T_0);
T1k_OR = zeros(1,T_0);
T2k_OR = zeros(1,T_0);
LK_OR = zeros(1,T_0);
Z1=zeros(1,T_0+1);
Z2=zeros(1,T_0+1);
Z=zeros(1,T_0+1);
Z1(1)=1;
Z2(1)=1;
Z(1)=1;
for tk=1:T_0     
    d_1_OR(:,1,tk) =pdist2(Y_n(:,:,tk),OR,'minkowski',d);
    N1(:,:,tk) = sort([d_0_OR;d_1_OR(:,:,tk)]);
    [~,N1_order(:,:,tk)]=ismember(d_1_OR(:,:,tk),N1(:,:,tk));
    T1k_OR(tk) = (sum(N1_order(:,:,tk))-mu_W)^2/sigma_W;
    N1_ORDER=sum(abs( N1_order(:,:,tk)-(N+1)/2));
    T2k_OR(tk)=(N1_ORDER-mu_AB)^2/sigma_AB;
    Z1(tk+1)=lambda*T1k_OR(tk)+(1-lambda)*Z1(tk);
    Z2(tk+1)=lambda*T2k_OR(tk)+(1-lambda)*Z2(tk);
    Z(tk+1)=max(Z1(tk+1),Z2(tk+1));
    if  Z(tk+1)>H
        RL=tk;
        break
    end
end
end
