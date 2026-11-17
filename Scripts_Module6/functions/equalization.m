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
    %Q9%%%%%
    h = (S' * S) \ (S' * r); % LS estimate h = (S^H S)^-1 S^H r.
    %Q9%%%%%
    H = fft(h);

    if strcmp(equalizerType, "zf")
        %Q11%%%%%
        E = 1 ./ H; % Zero Forcing (ZF) Channel Equalizer.
        %Q11%%%%%
    elseif strcmp(equalizerType, "mmse")
        sigma2 = var(r - S * h);    % Noise variance estimation.
        Ps     = var(s);            % Power of the transmitted signal estimation.

        %Q13%%%%%
        E = conj(H) ./ (abs(H).^2 + sigma2 / Ps); % Minimum Mean Square Error (MMSE) Equalizer.
        %Q13%%%%%
    else
        error("Choose a valid equalizer type.");
    end

    e = ifft(E);
    rDataEqualized = filter(e, 1, rData); % Applies the equalizer filter.
end
