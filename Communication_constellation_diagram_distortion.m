%% OFDM 星座图畸变演示（AWGN 信道）
% 生成并绘制：
% 1) 仅 AWGN（理想）
% 2) CFO 载波频偏
% 3) SFO 采样频偏
% 4) CTO 粗定时偏移（整数采样，FFT窗“提前”落在CP内）
% 5) STO 细定时偏移（小数采样，FFT窗“提前”落在CP内）
% 6) IQ 不平衡
%
% 说明：
% - CTO/STO 这里按 OFDM 常见定义：FFT 窗在 CP 内“提前” → 不产生ISI，但会造成子载波线性相位斜坡与（细偏移时）一定ICI。
% - CFO 这里用“归一化到子载波间隔”的频偏 epsilon：f_off = epsilon * Δf，其中 Δf = Fs/Nfft
% - SFO 这里用采样率相对偏差 rho：Fs_rx = Fs*(1+rho)

clear; close all; clc;

%% ---------------- 参数 ----------------
Nfft  = 64;
Ncp   = 16;
M     = 16;          % 16-QAM
Nsym  = 100;         % OFDM 符号数（想让SFO更明显可增大，但别让累计漂移超过CP太多）
SNRdB = 20;          % AWGN SNR（按“每个采样点”的平均功率定义）

% 使用子载波（类似 802.11a 的 52 个，去掉 DC）
usedSc = [-26:-1 1:26];         % 以 DC 为 0 的子载波编号
Ndata  = numel(usedSc);

% 各类非理想参数（可自行调大/调小以观察更明显的畸变）
epsCFO  = 0.08;      % CFO：归一化频偏（单位：子载波间隔）
rhoSFO  = 1e-3;      % SFO：采样频偏（1e-3=1000 ppm）；太大可能引入ISI/越过CP
ctoInt  = 5;         % CTO：整数采样的“提前”偏移（samples），必须 <= Ncp
stoFrac = 0.3;      % STO：小数采样的“提前”偏移（samples），建议 < 1 且 <= Ncp

% IQ 不平衡参数
gIQ   = 0.15;                 % 幅度不平衡（15%）
phiIQ = deg2rad(15);           % 相位不平衡（15°）

%% ---------------- 发射端：生成OFDM波形 ----------------
[txWave, txQamSyms] = ofdm_tx(Nfft, Ncp, M, Nsym, usedSc);

% 加 AWGN
awgn_add = @(x,snrdb) x + sqrt(mean(abs(x).^2) / (10^(snrdb/10)) / 2) .* ...
                     (randn(size(x)) + 1j*randn(size(x)));

%% ---------------- Case 1：理想（仅AWGN） ----------------
rx1 = awgn_add(txWave, SNRdB);
rxSyms1 = ofdm_rx(rx1, Nfft, Ncp, usedSc, 'none', 0, 0);

%% ---------------- Case 2：CFO ----------------
n = (0:numel(txWave)-1).';
rx2 = txWave .* exp(1j*2*pi*epsCFO*n/Nfft);   % epsilon * Δf 的离散实现
rx2 = awgn_add(rx2, SNRdB);
rxSyms2 = ofdm_rx(rx2, Nfft, Ncp, usedSc, 'none', 0, 0);

%% ---------------- Case 3：SFO ----------------
% 模拟：接收机采样更快 Fs_rx=Fs*(1+rho) -> 相当于对原信号按 n/(1+rho) 取样
rx3 = apply_sfo(txWave, rhoSFO);
rx3 = awgn_add(rx3, SNRdB);
rxSyms3 = ofdm_rx(rx3, Nfft, Ncp, usedSc, 'none', 0, 0);

%% ---------------- Case 4：CTO（整数定时偏移，提前落在CP内） ----------------
rx4 = awgn_add(txWave, SNRdB);
rxSyms4 = ofdm_rx(rx4, Nfft, Ncp, usedSc, 'cto', ctoInt, 0);

%% ---------------- Case 5：STO（小数定时偏移，提前落在CP内） ----------------
rx5 = awgn_add(txWave, SNRdB);
rxSyms5 = ofdm_rx(rx5, Nfft, Ncp, usedSc, 'sto', 0, stoFrac);

%% ---------------- Case 6：IQ 不平衡 ----------------
rx6 = apply_iq_imbalance(txWave, gIQ, phiIQ);
rx6 = awgn_add(rx6, SNRdB);
rxSyms6 = ofdm_rx(rx6, Nfft, Ncp, usedSc, 'none', 0, 0);

%% ---------------- 绘图 ----------------
figure('Color','w');
titles = { ...
  sprintf('Ideal (AWGN only), SNR=%g dB', SNRdB), ...
  sprintf('CFO: \\epsilon=%.3f (in subcarrier spacing)', epsCFO), ...
  sprintf('SFO: \\rho=%.0f ppm', rhoSFO*1e6), ...
  sprintf('CTO: %d samples early (within CP)', ctoInt), ...
  sprintf('STO: %.2f samples early (within CP)', stoFrac), ...
  sprintf('IQ imbalance: g=%.2f, \\phi=%.1f°', gIQ, rad2deg(phiIQ))};

rxAll = {rxSyms1, rxSyms2, rxSyms3, rxSyms4, rxSyms5, rxSyms6};

for k = 1:numel(rxAll)
    subplot(2,3,k);
    plot(real(rxAll{k}), imag(rxAll{k}), '.', 'MarkerSize', 6);
    grid on; axis square;
    xlabel('I'); ylabel('Q'); title(titles{k});
    xlim([-2 2]); ylim([-2 2]);   % 16QAM 单位平均功率下大致范围
end












