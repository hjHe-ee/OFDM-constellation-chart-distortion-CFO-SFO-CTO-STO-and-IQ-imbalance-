function s = qammod_fallback(idx, M)
    % 简易方形QAM映射（非Gray），并归一化到单位平均功率
    m = sqrt(M);
    if abs(m-round(m)) > 0
        error('Fallback supports square QAM only (M=4,16,64,...)');
    end
    m = round(m);
    I = mod(idx, m);
    Q = floor(idx / m);

    % 映射到 [-m+1, ..., m-1] 的奇数电平
    re = 2*I - (m-1);
    im = 2*Q - (m-1);
    s  = re + 1j*im;

    % 单位平均功率归一化
    s = s ./ sqrt(mean(abs(s).^2));
end