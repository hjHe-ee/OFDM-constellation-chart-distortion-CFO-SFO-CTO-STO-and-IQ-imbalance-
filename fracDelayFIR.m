function [y, gd] = fracDelayFIR(x, D, L)
%FRACDELAYFIR  toolbox-free 分数延时 FIR（windowed-sinc）
% y[n] ≈ x[n - (gd + D)]，其中 gd=(L-1)/2 是滤波器固有群时延
%
% 输入:
%   x : 输入序列（列向量）
%   D : 分数延时（可以为负，负值=提前）
%   L : 滤波器长度（必须为奇数）
%
% 输出:
%   y  : 滤波后序列（已在尾部补零，便于取样不越界）
%   gd : 群时延 = (L-1)/2

    if mod(L,2) ~= 1
        error('Filter length L must be odd.');
    end
    x = x(:);
    gd = (L-1)/2;

    n = (0:L-1).';
    m = n - gd;  % m = -gd : +gd

    % windowed-sinc 系数：h[m] = sinc(m - D) * w[m]
    h = sinc(m - D);

    % 自己实现 Hamming 窗（避免依赖工具箱）
    w = 0.54 - 0.46*cos(2*pi*n/(L-1));
    h = h .* w;

    % DC增益归一化（避免幅度偏差）
    h = h / sum(h);

    % 尾部补零以防取样越界（尤其 idx0+Nfft-1+gd）
    xpad = [x; zeros(L,1)];
    y = filter(h, 1, xpad);
end