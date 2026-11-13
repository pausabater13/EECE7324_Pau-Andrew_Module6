function rDataEqualized = equalization(rData, r, s, equalizerType, L)
    %EQUALIZATION Applies a FIR filter to equalize the channel impulse response.
    %
    % Input:
    %   rData             [samples or symbols] Received data.
    %   r                 [samples or symbols] Received pilot.
    %   s                 [samples or symbols] Original pilot.
    %   equalizerType     Type of equalizer ("zf" or "mmse").
    %   L                 [samples or symbols] Length of the FIR equalizer filter.
    %
    % Output:
    %   rDataEqualized    [samples or symbols] Equalized data.

    % Least Squares (LS) Channel Impulse Response (CIR) Estimation.
    S = convmtx(s, L);
    S = S(1 : numel(s), :);
    h = NaN;
    H = fft(h);

    if strcmp(equalizerType, "zf")
        E = NaN; % Zero Forcing (ZF) Channel Equalizer.
    elseif strcmp(equalizerType, "mmse")
        sigma2 = var(r - S * h);    % Noise variance estimation.
        Ps     = var(s);            % Power of the transmitted signal estimation.

        E = NaN; % Minimum Mean Square Error (MMSE) Equalizer.
    else
        error("Choose a valid equalizer type.");
    end

    e = ifft(E);
    rDataEqualized = filter(e, 1, rData); % Applies the equalizer filter.
end
