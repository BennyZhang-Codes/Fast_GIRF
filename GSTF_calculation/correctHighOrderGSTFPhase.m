function correctHighOrderGSTFPhase(H, phaseAtZero, corr_delay, referenceChannel)
% Correct High-Order GSTF phase using the physical self-term as reference.
%
% Arguments:
%   H:                GIRF_matrix object
%   phaseAtZero:      reference phase from H_FFT
%   corr_delay:       timing correction used by Fast_GIRF [s]
%   referenceChannel: spatial self-term channel
%                     x input -> 2, y input -> 3, z input -> 4
%
% This function reproduces GIRF_matrix.correct_GSTF_phase() while allowing
% the reference channel to be selected for the full spatial basis. The
% original GIRF_matrix.m is therefore kept unchanged for backwards
% compatibility with the two-slice Fast_GIRF workflow.

if referenceChannel < 1 || referenceChannel > size(H.gstf,2)
    error('referenceChannel is outside the available GSTF channels.');
end

linear_array = (1:1:size(H.gstf,1))*pi;
corrected_phase = angle(H.gstf.*exp(1i*(H.f_axis*2*pi*corr_delay+linear_array)).');
corrected_phase_atZero = corrected_phase(floor(size(corrected_phase,1)/2)+3,referenceChannel);

if abs(corrected_phase_atZero - phaseAtZero) > 1
    linear_array = (0:1:size(H.gstf,1)-1)*pi;
end

H.gstf = H.gstf.*exp(1i*(H.f_axis*2*pi*corr_delay+linear_array)).';

end
