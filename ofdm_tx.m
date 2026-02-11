function [txWave, qamSyms] = ofdm_tx(Nfft, Ncp, M, Nsym, usedSc)
    % QAM symbols (unit average power). If Communications Toolbox is missing,
    % fall back to simple square QAM mapping (Gray not guaranteed).
    Ndata = numel(usedSc);
    dataIdx = randi([0 M-1], Ndata*Nsym, 1);

    if exist('qammod','file') == 2
        qamSyms = qammod(dataIdx, M, 'UnitAveragePower', true);
    else
        qamSyms = qammod_fallback(dataIdx, M);
    end

    qamMat = reshape(qamSyms, Ndata, Nsym);

    L = Nfft + Ncp;
    txSymTime = zeros(L, Nsym);

    usedIdx = usedSc + Nfft/2 + 1;   % fftshift 后的索引（DC 在 Nfft/2+1）
    for s = 1:Nsym
        Xshift = zeros(Nfft,1);
        Xshift(usedIdx) = qamMat(:,s);

        X = ifftshift(Xshift);       % 转回 IFFT 的自然频域排列
        xt = ifft(X, Nfft);

        txSymTime(:,s) = [xt(end-Ncp+1:end); xt]; % 加 CP
    end

    txWave = txSymTime(:);
end