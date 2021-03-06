function PFB_channelizer_complex(Nchan,OS_Nu,OS_De,Nin,Nblocks,fname_in,fname_out)
%
% Takes as input a data file generated by "gen_test_vector.m".
% Assumes one polarization, complex valued single-precision inputs.
% Passes the input through the PFB and stores the output of all channels
% to file.
% The following PFB parameters are selectable:
%     - critically sampled or oversampled
%     - number of PFB channels
%     - oversampling factor
%     - prototype filter response (designed separately and its coefficients
%       provided in a file)
%
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% I. Morrison      31-Jul-2015  Original version
%
% I. Morrison      25-Oct-2015  Complex input version
%
% ----------------------------------------------------------------------

% Define globals common also to CS_PFB() / OS_PFB() sub-functions
global L; global Nu; global M; global L_M; global fname_pfb;

% Number of channels in filter-bank
L = Nchan;

% OS factor numerator
Nu = OS_Nu;

% PFB type
pfb_type = 1; % 0 for critically sampled, 1 for oversampled

if pfb_type == 0,
    OS = 1;
    M = L;
else
    OS = OS_Nu/OS_De;
    M = (L*OS_De)/OS_Nu;
end

L_M = L - M; % Overlap

% PFB prototype filter coefficients file name
fname_pfb = 'Prototype_FIR.mat';

%=======================================

% Open input file
fid_in = fopen(fname_in);

% Open files for writing
for i = 1 : L
    fid_out(i) = fopen(strcat(fname_out,int2str(i),'.dump'), 'w');
end;
fid_out_all = fopen(strcat(fname_out, 'all', '.dump'),'w');
fwrite(fid_out_all, char('0'*zeros(4096, 1)), 'char');

% Initialise output
y2 = zeros(L,Nin/M);
dat_out_per_chan_len = Nin/M;
y_out = zeros(L,2*dat_out_per_chan_len*Nblocks);

%===============
% Main loop
% Read input blocks and filter

for ii = 1 : Nblocks

    % Print loop number
    fprintf('Loop # %i of %i\n', ii, Nblocks);

    % Read stream of voltages into a single column
    Vstream = single(fread(fid_in, 2*Nin, 'single'));

    if feof(fid_in)
        error('Error - hit end of input file!');
    end;

    % Parse real and imag components
    Vstream = reshape(Vstream, 2, []);
    Vdat = complex(Vstream(1,:), Vstream(2,:));

    % Evaluate the channel outputs
    for n = 1 : Nin/M
        if pfb_type == 0,
            y2(:,n) = CS_PFB(Vdat(1,(n-1)*L+1:n*L));
        else
            y2(:,n) = OS_PFB(Vdat(1,(n-1)*M+1:1:n*M));
        end;
    end;
    % y2(:, 14)
    %Write each output channel's samples to its own file
    s = (ii - 1) * dat_out_per_chan_len * 2 + 1;
    e = ii * dat_out_per_chan_len * 2;
    for i = 1 : L
        % Interleave real/imag, selecting the particular channel
        yy = y2(i,:);
        z = [real(transpose(yy)), imag(transpose(yy))];
        dat = reshape(transpose(z),2*Nin/M,1);
        y_out(i,s:e) = dat;
        fwrite(fid_out(i), dat, 'single');
    end;

end;

fclose(fid_in);
for i = 1 : L
    fclose(fid_out(i));
end;
fwrite(fid_out_all, reshape(y_out, L*2*dat_out_per_chan_len*Nblocks, 1),'single');
fclose(fid_out_all);

return
end



% CS-PFB
% Critically sampled Polyphase Filter-Bank Channelizer function, based on
% code by Thushara Kanchana Gunaratne, RO/RCO, NSI-NRC, Canada, 2015-03-05
function y = CS_PFB(x)

global L; global fname_pfb;

%Declaration and Initialization of Input Mask
%As Persistence Variables
persistent n h xM;
if isempty(n)

    %Loading the Prototype Filter as an initiation task
    %This Will NOT repeat in subsequent runs
    FiltCoefStruct = load(fname_pfb);
    h = FiltCoefStruct.h;

    %Initiate the Input Mask that is multiplied with the Filter mask
    xM = complex(zeros(1,length(h)));
    %Initiate the Output mask
    yP = complex(zeros(L,1));

    %Control Index - Initiation
    n = 0;

end; %End if

%Multiplying the Indexed Input Mask and Filter Mask elements and
%accumulating
for k = 1 : L
    yP(k,1) = sum(xM(k:L:end).*h(k:L:end));
end; % For k

%The Linear Shift of Input through the FIFO
%Shift the Current Samples by M to the Right
xM(1,L+1:end) = xM(1,1:end-L);
%Assign the New Input Samples for the first M samples
xM(1,1:L) = fliplr(x);%Note the Flip (Left-Right) place the Newest sample
                      % to the front

% FFT stage - real case
% % %Note the Input Signal is Real-Valued. Hence, only half of the output
% % %Channels are Independent. The Packing Method is used here. However,
% % %any Optimized Real IFFT Evaluation Algorithm Can be used in its place
% % %Evaluating the Cross-Stream (i.e. column wise) IDFT using Packing
% % %Method
% % %The Complex-Valued Sequence of Half Size
% y2C = yP(1:2:end) + 1j*yP(2:2:end);
% %The Complex IDFT of LC=L/2 Points
% IFY2C = L*L/2*ifft(y2C);
% %
% y(1:L/2) = (0.5*((IFY2C+conj(circshift(flipud(IFY2C),[+1,0])))...
%             - 1j*exp(2j*pi*(0:1:L/2-1).'/L).*...
%               (IFY2C-conj(circshift(flipud(IFY2C),[+1,0])))));
% % [0,+1]
% y(L/2+1) = 0.5*((IFY2C(1)+conj(IFY2C(1)) + 1j*(IFY2C(1)-conj(IFY2C(1)))));
%
% y(L/2+2:L) = conj(fliplr(y(2:L/2)));

% FFT stage - complex case
y = L*L/2*ifft(yP);
% y = ifft(yP);
% y = fft(yP);

%Changing the Control Index
n = n+1;

end %Function CS_PFB



% OS-PFB
% Oversampled Polyphase Filter-Bank Channelizer function, based on code by
% Thushara Kanchana Gunaratne, RO/RCO, NSI-NRC, Canada, 2015-03-05
% with mods by Gianni to deal with arbitrary OS factors
function y = OS_PFB(x)
global L; global Nu; global M; global L_M; global fname_pfb;

% Choose way circular shifts are done
PFB_code = 1;  % 0 for Thushara's, 1 for Gianni's

%Declaration and Initialization of Input Mask
%As Persistance Variables
persistent n h xM;
if isempty(n)

    %Loading the Prototype Filter as an initiation task
    %This Will NOT repeat in subsequent runs
    FiltCoefStruct = load(fname_pfb);
    h = FiltCoefStruct.h;

    %Initiate the Input Mask that is multiplied with the Filter mask
    xM = zeros(1,length(h));
    %Initiate the Output mask
    yP = zeros(L,1);

    %Control Index - Initiation
    n = 0;

end; %End if

%Multiplying the Indexed Input Mask and Filter Mask elements and
%accumulating
for k = 1 : L
    yP(k,1) = sum(xM(k:L:end).*h(k:L:end));
end; % For k

%The Linear Shift of Input through the FIFO
%Shift the Current Samples by M to the Right
xM(1,M+1:end) = xM(1,1:end-M);
%Assign the New Input Samples for the first M samples
xM(1,1:M) = fliplr(x);%Note the Flip (Left-Right) place the Newest sample
                      % to the front

%Performing the Circular Shift to Compensate the Shift in Band Center
%Frequencies
if(PFB_code == 0)   % Thushara's code
    if n == 0
        y1S = yP;
    else
        y1S = [yP((Nu-n)*L_M+1:end); yP(1:(Nu-n)*L_M)];
    end;
end;
if(PFB_code == 1)   % Gianni's code
    y1S=circshift(yP,[n 0]);
end;

% %Evaluating the Cross-Stream (i.e. column wise) IDFT
% yfft = L*L*(ifft(yP));%
%
% %Modulating the Channels (i.e. FFT Outputs) to compensate the shift in the
% %center frequency
% %y = yfft.*exp(2j*pi*(1-M/L)*n*(0:1:L-1).');
% y = yfft.*exp(-2j*pi*M/L*n*(0:1:L-1).');

% FFT stage - real case
% % %Note the Input Signal is Real-Valued. Hence, only half of the output
% % %Channels are Independent. The Packing Method is used here. However,
% % %any Optimized Real IFFT Evaluation Algorithm Can be used in its place
% % %Evaluating the Cross-Stream (i.e. column wise) IDFT using Packing
% % %Method
% % %The Complex-Valued Sequence of Half Size
% y2C = y1S(1:2:end) + 1j*y1S(2:2:end);
% %The Complex IDFT of LC=L/2 Points
%
% IFY2C = L*L/2*ifft(y2C);
% %
% y(1:L/2) = (0.5*((IFY2C+conj(circshift(flipud(IFY2C),[+1,0])))...
%             - 1j*exp(2j*pi*(0:1:L/2-1).'/L).*...
%              (IFY2C-conj(circshift(flipud(IFY2C),[+1,0])))));
% % [0,+1]
% y(L/2+1) = 0.5*((IFY2C(1)+conj(IFY2C(1)) + 1j*(IFY2C(1)-conj(IFY2C(1)))));
%
% y(L/2+2:L) = conj(fliplr(y(2:L/2)));


% FFT stage - complex case
y = L*L*ifft(y1S);
% y = L*ifft(y1S);
% y = fft(y1S);


%Changing the Control Index
if(PFB_code == 0)   % Thushara's code
    n = n+1;
    n = mod(n,Nu);
end;
if(PFB_code == 1)   % Gianni's code
    n=mod(n+L_M,L);
end;

end %Function OS_PFB
