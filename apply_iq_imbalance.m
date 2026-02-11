function y = apply_iq_imbalance(x, g, phi)
    % 经典 IQ 不平衡模型：y = alpha*x + beta*conj(x)
    a = (1+g);
    b = (1-g);
    alpha = 0.5*(a*exp(-1j*phi/2) + b*exp( 1j*phi/2));
    beta  = 0.5*(a*exp(-1j*phi/2) - b*exp( 1j*phi/2));
    y = alpha*x + beta*conj(x);
end