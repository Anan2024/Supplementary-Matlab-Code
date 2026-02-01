function MRL = calculate_OOC_MRL(H, lambda, p, r, m, n, T_0, d,delta_1,delta_2)

RL_values=ones(1,r)*T_0;
parfor ii=1:r
    RL_values(ii)=compute_conditional_OOC_RL(H, lambda, p, m, n, T_0, d,delta_1,delta_2);
end
MRL=median(RL_values);
end
