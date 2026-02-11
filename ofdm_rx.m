function rxSyms = ofdm_rx(rxWave, Nfft, Ncp, usedSc, mode, ctoInt, stoFrac)
    % mode:
    % 'none' : 正常去CP
    % 'cto'  : FFT窗起点提前 ctoInt 个采样（整数）
    % 'sto'  : FFT窗起点提前 stoFrac 个采样（小数），用插值取样
    L = Nfft + Ncp;
    Nsym = floor(numel(rxWave)/L);
    rxWave = rxWave(1:Nsym*L);
    rxMat  = reshape(rxWave, L, Nsym);

    Ndata   = numel(usedSc);
    rxSyms  = zeros(Ndata*Nsym, 1);
    usedIdx = usedSc + Nfft/2 + 1;

    for s = 1:Nsym
        r = rxMat(:,s);

        switch lower(mode)
            case 'cto'
                % FFT窗起点：理想是 Ncp+1；提前 ctoInt -> Ncp+1-ctoInt
                start = (Ncp + 1) - ctoInt;
                if start < 1, error('CTO too large: start<1. Keep ctoInt<=Ncp.'); end
                idx = start + (0:Nfft-1);
                rfft_in = r(idx);

            case 'sto'
                % 小数提前：toolbox-free 分数延时 FIR（windowed-sinc）
                % 目标：rfft_in[n] = r( start0 + n ), start0 = Ncp - stoFrac (0-based)
                start0 = (Ncp) - stoFrac;   % 0-based
                if start0 < 0
                    error('STO too large: keep stoFrac<=Ncp.');
                end

                startInt = floor(start0);           % 整数部分
                frac     = start0 - startInt;       % 小数部分 in [0,1)
                D        = -frac;                   % 需要的“等效延时”（负值=提前）

                % 分数延时滤波器长度（奇数，越长越接近理想sinc，运算更大）
                Lfd  = 41;                          % 推荐 21~81 之间
                [y, gd] = fracDelayFIR(r, D, Lfd);   % y ≈ r[n - (gd + D)]，gd=(Lfd-1)/2

                % 关键点：当 D=-frac 时，取 y(startInt + gd + n) 就对齐到 r(startInt+frac+n)
                idx0 = startInt + gd + 1;           % 1-based
                rfft_in = y(idx0 + (0:Nfft-1));


            otherwise
                rfft_in = r(Ncp+1 : Ncp+Nfft);
        end

        Y = fft(rfft_in, Nfft);
        Yshift = fftshift(Y);

        ys = Yshift(usedIdx);
        rxSyms((s-1)*Ndata+1:s*Ndata) = ys;
    end
end