% Figures and tables of the report. Not part of the given code.
% Same steps and functions as receiverSC.m, for every capture and setting of the report.
% Run from Scripts_Module6. Needs the four captures in ../captures/wireless/.

clear(); close("all"); addpath("functions/");
p.Rs = 25e6; p.Fs = 20e9; p.fc = 2.4e9; p.sps = p.Fs / p.Rs; p.nPilot = 60; p.nData = 300; % As in receiverSC.m.
p.band = [p.fc - 0.875*p.Rs, p.fc + 0.875*p.Rs];                                       % Passband (Step 2).
p.bits = eval(regexp(fileread("receiverSC.m"), "\[[01, ]+\]\.'", "match", "once"));   % Same bits as receiverSC.m.
caps = {"6b_16psk_gray", 16, "psk"; "6b_32psk_gray", 32, "psk"; "18b_16qam_gray", 16, "qam"; "18b_32qam_gray", 32, "qam"};
Ls = [1 2 3 5 10 50];

%Q3%%%%%
y = load("../captures/wireless/6b_16psk_gray.mat").Trace_1_VOLT;
[P0, f] = frequencyDomain(y, p.Fs, "single");
P1 = frequencyDomain(bandpass(y, p.band, p.Fs), p.Fs, "single");                  % Default steepness (0.85).
P2 = frequencyDomain(bandpass(y, p.band, p.Fs, Steepness=0.95), p.Fs, "single");  % Our steepness.
k = abs(f - p.fc) > 100e6 & abs(f - p.fc) < 400e6;                                % Next to the signal band.
fprintf("Noise floor drop next to the band: %.1f dB (0.85), %.1f dB (0.95)\n", 10 * log10(median(P0(k)) ./ [median(P1(k)), median(P2(k))]));
fig = figure();
plot(f / 1e9, 10 * log10(P0)); hold("on"); plot(f / 1e9, 10 * log10(P1)); plot(f / 1e9, 10 * log10(P2));
xlim([f(1) f(end)] / 1e9); xlabel("Frequency (GHz)"); ylabel("Power Spectrum (dB)");
legend("Before filtering", "After filtering, steepness 0.85 (default)", "After filtering, steepness 0.95");
exportgraphics(fig, "../figures/A/psd_filter.pdf");
%Q3%%%%%

%% Receiver up to the equalizer, once per capture, without and with the filter
for c = 1:4
    for useFilter = [false true]
        rx(c, useFilter + 1) = front(p, "../captures/wireless/" + caps{c, 1} + ".mat", caps{c, 2}, caps{c, 3}, useFilter);
    end
end

%% Results: {capture, filter, equalizer, L}. Equalizer NaN: no equalizer, only scaled to unit power.
res = [];
%Q4%%%%% %Q6%%%%%
for c = 1:4
    for useFilter = [false true]
        res = [res, back(rx(c, useFilter + 1), c, useFilter, NaN, NaN)];
    end
end
%Q4%%%%% %Q6%%%%%
%Q12%%%%% %Q14%%%%% %Q15%%%%%
for c = 1:4
    for L = Ls
        res = [res, back(rx(c, 2), c, true, "zf", L), back(rx(c, 2), c, true, "mmse", L)];
    end
end
%Q12%%%%% %Q14%%%%% %Q15%%%%%

psk = [1 2]; qam = [3 4]; % Rows of caps.
filt = {false, NaN, NaN, "No filter"; true, NaN, NaN, "Filter"};
eqL  = @(eq) {true, NaN, NaN, "No equalizer"; true, eq, 1, upper(eq) + ", L = 1"; true, eq, 3, upper(eq) + ", L = 3"; true, eq, 10, upper(eq) + ", L = 10"; true, eq, 50, upper(eq) + ", L = 50"};

%Q4%%%%%
plotGrid(res, caps, psk, filt, "../figures/A/rx_psk_filter.pdf");
filterTable(res, caps, psk, "../tables/filter_psk.tex");
%Q4%%%%%

%Q6%%%%%
plotGrid(res, caps, qam, filt, "../figures/A/rx_qam_filter.pdf");
filterTable(res, caps, qam, "../tables/filter_qam.tex");
%Q6%%%%%

%Q12%%%%%
plotGrid(res, caps, psk, eqL("zf"), "../figures/B/rx_psk_zf.pdf");
eqTable(res, caps, psk, Ls, "zf", "../tables/zf_psk.tex");
%Q12%%%%%

%Q14%%%%%
plotGrid(res, caps, psk, eqL("mmse"), "../figures/B/rx_psk_mmse.pdf");
eqTable(res, caps, psk, Ls, "mmse", "../tables/mmse_psk.tex");
%Q14%%%%%

%Q15%%%%%
plotGrid(res, caps, qam, {true, NaN, NaN, "No equalizer"; true, "zf", 10, "ZF, L = 10"; true, "mmse", 10, "MMSE, L = 10"}, "../figures/B/rx_qam_eq.pdf");
eqTable(res, caps, qam, Ls, ["zf" "mmse"], "../tables/eq_qam.tex");
%Q15%%%%%

%% Receiver steps (as in receiverSC.m)
function r = front(p, file, M, modScheme, useFilter) % From the capture to the received pilot and data symbols.
    yDSO = load(file).Trace_1_VOLT;
    symbolsHeader = symbolMapping([1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1].', 2, "psk"); % 18-bit MF sequence, as in receiverSC.m.
    xHeader = iqModulation(pulseShaping(symbolsHeader, p.sps), p.fc, p.Fs);
    r.M = M; r.modScheme = modScheme;
    r.sPilot = symbolMapping(p.bits(1 : p.nPilot * log2(M)), M, modScheme);
    r.sData  = symbolMapping(p.bits(1 : p.nData * log2(M)), M, modScheme);
    r.bits   = p.bits(1 : p.nData * log2(M));
    yFiltered = yDSO;
    if useFilter, yFiltered = bandpass(yDSO, p.band, p.Fs, Steepness=0.95); end
    nHeader = numel(xHeader); nPilot = (p.nPilot + 10) * p.sps; nData = (p.nData + 10) * p.sps;
    yFrame = frameSynchronization(yFiltered, xHeader, nHeader + nPilot + nData);
    c = xcorr(yFiltered, xHeader); [~, i] = max(c); r.start = i - numel(yFiltered) + 1;          % Frame start (sample), Step 5,
    if r.start + numel(yFrame) - 1 > numel(yFiltered)                                               % found as in frameSynchronization.m.
        [~, i] = max(c(1 : i - numel(yFrame) + 1)); r.start = i - numel(yFiltered) + 1;
    end
    assert(isequal(yFrame, yFiltered(r.start + (0 : numel(yFrame) - 1))));      % Same as frameSynchronization.m.
    rFrame = iqDemodulation(yFrame, p.fc, p.Fs, p.Rs);
    sHeader = correlationReceiver(rFrame(1 : nHeader), p.sps, numel(symbolsHeader));
    rFrame = rFrame * exp(-1j * angle(sum(sHeader .* conj(symbolsHeader)))); % Static phase correction.
    r.pilot = correlationReceiver(rFrame(nHeader + (1 : nPilot)), p.sps, p.nPilot);
    r.data  = correlationReceiver(rFrame(nHeader + nPilot + (1 : nData)), p.sps, p.nData);
end

function out = back(r, c, useFilter, eq, L) % Equalizer, demapping, BER, SER and EVM.
    if isequaln(eq, NaN)
        s = r.data / rms(r.data) * rms(r.sData);
    else
        s = equalization(r.data, r.pilot, r.sPilot, eq, L);
    end
    bits = symbolDemapping(s, r.M, r.modScheme);
    err = bit2int(r.bits, log2(r.M)) ~= bit2int(bits, log2(r.M));
    out = struct(c=c, useFilter=useFilter, eq=eq, L=L, sym=s, err=err, ber=mean(r.bits ~= bits), ser=mean(err), evm=rms(abs(s - r.sData)) * 100, start=r.start);
end

%% Figures and tables
function r = pick(res, c, useFilter, eq, L)
    r = res(arrayfun(@(x) isequaln({x.c, x.useFilter, x.eq, x.L}, {c, useFilter, eq, L}), res));
end

function s = fmtBer(x), s = replace(sprintf("%.4f", x), "0.0000", "0"); end

function filterTable(res, caps, rows, file) % BER, EVM and frame start without and with the filter.
    out = "";
    for c = rows
        n = pick(res, c, false, NaN, NaN); f = pick(res, c, true, NaN, NaN);
        out = out + sprintf("%d-%s & %s & %s & %.1f & %.1f & %d & %d \\\\\n", caps{c, 2}, upper(caps{c, 3}), fmtBer(n.ber), fmtBer(f.ber), n.evm, f.evm, n.start, f.start);
    end
    writelines(out + "\hline", file);
end

function eqTable(res, caps, rows, Ls, eqs, file) % BER and EVM without equalizer and for each L.
    out = "";
    for c = rows
        n = pick(res, c, true, NaN, NaN);
        out = out + sprintf("%d & --", caps{c, 2}) + join(repmat(string(sprintf(" & %s & %.1f", fmtBer(n.ber), n.evm)), 1, numel(eqs)), "") + " \\" + newline;
        for L = Ls
            out = out + sprintf("%d & %d", caps{c, 2}, L);
            for eq = eqs
                r = pick(res, c, true, eq, L); out = out + sprintf(" & %s & %.1f", fmtBer(r.ber), r.evm);
            end
            out = out + " \\" + newline;
        end
        out = out + "\hline" + newline;
    end
    writelines(strip(out, "right", newline), file);
end

function plotGrid(res, caps, rows, panels, file) % Constellations: one row per capture, one column per setting.
    fig = figure(Units="centimeters", Position=[2 2 4.2 * size(panels, 1) 4.6 * numel(rows)]);
    tiledlayout(numel(rows), size(panels, 1), TileSpacing="compact", Padding="compact");
    for c = rows
        [M, m] = caps{c, 2:3};
        for j = 1:size(panels, 1)
            r = pick(res, c, panels{j, 1:3});
            nexttile;
            plot(r.sym(~r.err), ".", MarkerSize=5); hold("on");
            plot(r.sym(r.err), "r.", MarkerSize=5);
            plot(symbolMapping(int2bit(0 : M - 1, log2(M)), M, m), "kx", LineWidth=1, MarkerSize=5);
            axis("equal"); xlim([-1.6 1.6]); ylim([-1.6 1.6]); set(gca, XTick=[], YTick=[], Box="on");
            title(sprintf("%d-%s, %s", M, upper(m), panels{j, 4}), FontSize=7, FontWeight="normal");
            subtitle(sprintf("SER %.1f%%, EVM %.1f%%", r.ser * 100, r.evm), FontSize=6);
        end
    end
    exportgraphics(fig, file, ContentType="vector");
end
