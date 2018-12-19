fprintf('\nTest of OS-PFB Inversion via FFT\n');

%% GLOBAL PARAMETERS

% Number of PFB output channels - power of 2, min OS_Nu, max 256
N = 64;

% PFB oversampling factor
OS_Nu = 32;  % numerator - should be a sub-multiple of N
OS_De = 27;  % denominator

% Width of PFB channel passband in MHz = spacing of PFB output channels
fine_chan_passband = 0.8;

% Length of forward FFT to process fine channels
ffft_length = 2^10;

% Length of test vector blocks (spacing of impusles)
block_length = 2*N*ffft_length;


%% GENERATE TEST VECTOR (input to PFB)

test_vector_filename = 'test_vec.dump';

Wave_type = 0;  % 0 for pulsar, 1 for impulse
impulse_offset = block_length/4;  % location of impulse within each block
impulse_width = 1;  % number of samples width of impusle
f_sample_out = N*fine_chan_passband;  % sample rate in MHz
period = 0.001;  % simulated pulsar period in seconds
noise = 1.0;  % sets SNR of simulated pulsar signal

fprintf('\nGenerating test vector...\n');
gen_test_vector_complex(Wave_type,impulse_offset,impulse_width,block_length,1,f_sample_out,period,noise,test_vector_filename);


%% DESIGN PFB PROTOTYPE FILTER
% function design_PFB(Nchan,OS_Nu,OS_De,Ntaps-1,ffft_len,display)  % subtract 1 from num taps because design functions adds 1
% taps_per_chan = 8;
% Ntaps = N*taps_per_chan - 1;  % must be odd
Ntaps = 385;

display = 0;    % 1 to display filter design plot, 0 otherwise

fprintf('\nDesigning PFB prototype filter...\n');
if (display == 1)
    fprintf('\nPress any key to continue...\n');
end;
design_PFB(N,OS_Nu,OS_De,Ntaps-1,ffft_length,display);  % subtract 1 from num taps because design functions adds 1


%% PFB Channelize - one block
% function PFBchannelizer(Nchan,OS_Nu,OS_De,Nin,Nblocks,fname_in,fname_out)
% minimum Nin is (block_length/OS_factor) - can be longer
fprintf('\nChannelizing...\n');
%PFB_channelizer_complex(N,OS_Nu,OS_De,OS_De*block_length/OS_Nu,1,test_vector_filename,'fine_channel_');
PFB_channelizer_CSIRO(N,OS_Nu,OS_De,OS_De*block_length/OS_Nu,1,test_vector_filename,'fine_channel_');


%% PROCESS EACH FINE CHANNEL
input_offset = 128;  % number of samples to drop at the start of the PFB output data, to ensure impulse within window
equalise_ripple = 1;  % 1 to equalise PFB ripple, 0 to not
fprintf('\nProcessing each channel...\n');
for chan = 1:N,
    fprintf('channel %d\n', chan);
    % function fine_chan_proc(chan,Nin,OS_Nu,OS_De,input_offset,fname_in,fname_out,equalise_ripple)
    fine_chan_proc(chan,ffft_length,OS_Nu,OS_De,input_offset,strcat('fine_channel_',int2str(chan),'.dump'),strcat('chunk_',int2str(chan),'.mat'),equalise_ripple);
end;


%% Combine chunks, back-transform and compare to original
% function invert(Nchan,Nin,fname_in,fname_compare,compare_offset)
fprintf('\nCombining channels and back transforming...\n');
compare_offset = -(Ntaps-1)/2 - (OS_De*N/OS_Nu)*input_offset;
invert(N,OS_Nu,OS_De,block_length,'chunk_',test_vector_filename,compare_offset);

fprintf('\nDone! Press any key to close plots and exit...\n\n');
pause;
close all;
clear all;
