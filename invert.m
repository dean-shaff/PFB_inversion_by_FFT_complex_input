function invert(Nchan,OS_Nu,OS_De,Nin,fname_in,fname_compare,compare_offset)
% Combines mutiple sub-channel pass-band chunks (from an oversampled PFB
% that have had their transition bands discarded) into a single contiguous
% block, then inverse FFTs.  The output is compared with the original in
% both the time and frequency domains.
% Ian Morrison
% 21-4-16

% load and concatenate the chunks
chan = 1;
load(strcat(fname_in,int2str(chan),'.mat'));
FFFF = F1((length(F1)/2)+1:length(F1)); % upper half is first part of FFFF
for chan = 2 : Nchan,
    load(strcat(fname_in,int2str(chan),'.mat'));
    FFFF = [FFFF; F1];
end;
chan = 1;
load(strcat(fname_in,int2str(chan),'.mat'));
FFFF = [FFFF; F1(1:(length(F1)/2))]; % lower half is last part of FFFF

len = length(FFFF);

save('N_channels','FFFF');

figure;%set(gcf,'Visible', 'off');
subplot(211); plot((1:len),abs(FFFF)); box on; grid on; title('FFFF Mag');
subplot(212); plot((1:len),angle(FFFF)); box on; grid on; title('FFFF Phase'); xlabel('time');

% back transform
z1 = (ifft((FFFF), len))./(OS_Nu/OS_De);  % re-scale by OS factor
length(z1)
len
figure;%set(gcf,'Visible', 'off');
subplot(211); plot((1:len),real(z1(1:len))); box on; grid on; title('z1 Real');
subplot(212); plot((1:len),imag(z1(1:len))); box on; grid on; title('z1 Imag'); xlabel('time');

figure;%set(gcf,'Visible', 'off');
subplot(211); plot((1:len),10.0*log10(abs(real(z1(1:len)))+1e-12)); box on; grid on; title('z1 Real - Log scale');
axis([1 len -100 10]);
subplot(212); plot((1:len),10.0*log10(abs(imag(z1(1:len)))+1e-12)); box on; grid on; title('z1 Imag - Log scale'); xlabel('time');
axis([1 len -100 10]);


% Compare with original
% =====================

fid_in = fopen(fname_compare);

% Read original input to PFB
Vstream = single(fread(fid_in, 2*Nin, 'single'));

if feof(fid_in)
    error('Error - hit end of input file!');
end;

% Parse real and imag components
Vstream = reshape(Vstream, 2, []);
Vdat = complex(Vstream(1,:), Vstream(2,:));

length(Vstream)

figure;%set(gcf,'Visible', 'off');
subplot(211); plot((1:Nin),real(Vdat(1,1:Nin))); box on; grid on;
title('Vin Real');
subplot(212); plot((1:Nin),imag(Vdat(1,1:Nin))); box on; grid on;
title('Vin Imag'); xlabel('time');

% Time domain comparison with original input - good for integer sample
% delays when those delays have been sync'd out
if (1)
    centre_Vdat = Nin/4;
    centre_z1 = centre_Vdat + compare_offset;

    plot_range = 25;
    figure; %set(gcf,'Visible', 'off');
    t_plot = (-plot_range+1:plot_range);
    z1_plot = z1(centre_z1-plot_range+1:centre_z1+plot_range);
    Vdat_plot = Vdat(1,centre_Vdat-plot_range+1:centre_Vdat+plot_range);
    subplot(311); plot(t_plot,...
      real(z1_plot),...
      t_plot,...
      real(Vdat_plot));
    box on; grid on; title('z1 vs Vdat Real'); legend({'z1 Real', 'Vdat Real'});
    subplot(312); plot(t_plot,...
      imag(z1_plot),...
      t_plot,...
      imag(Vdat_plot));
    box on; grid on; title('z1 vs Vdat Imag'); legend({'z1 Imag', 'Vdat Imag'});
    xlabel('time');

    xcorr_z1_Vdat = xcorr(z1_plot, Vdat_plot);
    subplot(313); plot((1:length(xcorr_z1_Vdat)), real(xcorr_z1_Vdat), (1:length(xcorr_z1_Vdat)), imag(xcorr_z1_Vdat));
    box on; grid on; title('Cross correlation z1 vs Vdat');
    % calculate total RMS error over Nsamp samples
    Nsamp = 200;
    Realerr = 0;
    Imagerr = 0;

    for j = centre_Vdat-(Nsamp/2):centre_Vdat+(Nsamp/2-1),
        Realerr = Realerr + (real(Vdat(1,j)) - real(z1(j+compare_offset)))^2;
        Imagerr = Imagerr + (imag(Vdat(1,j)) - imag(z1(j+compare_offset)))^2;
    end;

    RMSerr_real = (Realerr/Nsamp)^0.5
    RMSerr_imag = (Imagerr/Nsamp)^0.5

end;


% Frequency domain comparison with original input - good for fractional
% sample delays to avoid the need for interpolation of the input time series
if (1)
    Vdat_shift = Vdat(1,1-compare_offset:len-compare_offset);
    VDAT = fft(Vdat_shift);
%     figure;
%     subplot(211); plot((1:len),abs(VDAT)); box on; grid on; title('VDAT Mag');
%     subplot(212); plot((1:len),angle(VDAT)); box on; grid on; title('VDAT Phase'); xlabel('time');

    % cross-power spectrum of FFFF and a similar length VDAT
    CP = FFFF.*transpose(conj(VDAT));
    figure;%set(gcf,'Visible', 'off');
    subplot(211); plot((1:len),abs(CP),'Linewidth',2);
    box on; grid on; title('Cross-Power Mag'); set(gca, 'YScale', 'log')
    subplot(212); plot((1:len),angle(CP),'Linewidth',2);
    box on; grid on; title('Cross-Power Phase'); xlabel('time');

    % sum the total cross-power magnitude
    total_cp = sum(abs(CP))

  end;

return
end
