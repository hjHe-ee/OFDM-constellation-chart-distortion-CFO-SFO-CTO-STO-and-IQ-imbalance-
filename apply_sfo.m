function y = apply_sfo(x, rho)
    % Fs_rx = Fs*(1+rho)，接收采样更快 -> 取样时刻更密，等效离散索引 tq = n/(1+rho)
    N  = numel(x);
    t  = (0:N-1).';
    tq = t/(1+rho);

    % 若 rho<0（更慢），tq 可能超界，需要外推；这里先给 0 外推
    y = interp1(t, x, tq, 'pchip', 0);
end