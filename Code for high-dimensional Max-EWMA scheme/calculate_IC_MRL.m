function MRL = calculate_IC_MRL(H, lambda, p, r, m, n, T_0, d)

for ii=1:r
    RL_values(ii)=compute_conditional_IC_RL(H, lambda, p, m, n, T_0, d);
end
MRL=median(RL_values);
end
