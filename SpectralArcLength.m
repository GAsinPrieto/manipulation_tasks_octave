%--------------------------------------------------------------------------
% EUROBENCH PROJECT H2020
%
% Benchmarking Scenario: Standing During Manipulation Skills
%
% Performance Indicators: Range of Motion, Movement Smoothness and Postural Balance 
%
%--------------------------------------------------------------------------
%% Spectral Arc Length. SPARC. Code AUTHOR: Sivakumar Balasubramanian. DATE: July 02, 2014.

% The function assumes that the movement speed profile is already filtered
% and segmented.

function S = SpectralArcLength (speed, Ts, parameters)

% speed: Is the speed profile of the movement. This is a 1xN row vector. N
% is the total number of points in the speed profile. 
%
% Ts: Sampling time in seconds. (Default value = 0.01 sec)
% NOTE: If your data was not sampled at 100 Hz, you must enter the
% appropriate sampling time.
%
% parameters: This contains the parameters to be used spectral arc length
% computation. This is a 1x2 column vector. 
%   - parameter(1): The amplitude threshold to be used to choose the
%   cut-off frequency. Default value = 0.05
%   - parameter (2): Maximum cut-off frequency for the spectral arc length
%   calculation. Default value = 20 Hz
% To represent the maximum frequency component of a movement. This was
% chosen to cover both normal and abnormal motor behaviour. 
% You can use a value lower than 20 Hz if you are aware of the maximum
% frequency component in the movement of interest. 
%   - parameter (3): Zero padding index. It controls the resolution of the
%   movement spectrum. Default value = 4
%
% Output: {S}
% S: Smoothness of the given movement. 


% Check input arguments.

switch nargin
    case 0
        disp('Error! Input at least the movement speed profile for the smoothness calculation.');
        fprintf('\n');
        help('SpectralArcLength');
        return;
    case 1
        % Default sampling time.
        Ts = 1/100; % 10 ms
        parameters = [0.05, 20, 4];
    case 2
        % Default parameters are used for the spectral arc length calculation.
        parameters = [0.05, 20, 4];
end;

% Check if the input argument are of the appropriate dimensions.

% Speed profile
sz = size(speed); %size returns a row vector

%if (sz(1)== 1) && (sz(2) > 1)
    %disp('Error! Speed must be a row vector.');
    %fprintf('\n');
    %help('SpectralArcLength');
    %return;
%end

% Sampling time
if ~isscalar(Ts)
    disp('Error! Ts must be a row scalar.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end;

% Parameters
%if length(parameters)~= 3
    %disp('Error! Parameter should have three elements.');
    %fprintf('\n');
    %help('SpectralArcLength');
    %return;
%end;

% Calculate the spectrum of the speed profile
N = length(speed);
Nfft = 2^nextpow2(ceil(log2(N))+ parameters(3));
speedSpectrum = abs(fft(speed, Nfft));
    
% Normalize spectrum with respect to the DC component
freq = 0:(1/Ts)*(1/Nfft):(1/Ts)*((Nfft-1)/Nfft);
speedSpectrum = speedSpectrum / max(speedSpectrum);

% Choose the spectrum that is always above the amplitude threshold and
% within the cut-off frequency
inxFc = find((freq(1:end) <= parameters(2)) & (speedSpectrum (1:end) >= parameters(1)), 1, 'last');

% Calculate the Spectral Arc Length
% 1. Select the spectrum of interest
speedSpectrum = speedSpectrum(1:inxFc);
% 2. Calculate the incremental arc lengths
dArcLengths = sqrt((1/(inxFc-1))^2 + (diff(speedSpectrum)).^2);
% 3. Compute movement smoothness S
S = - sum(dArcLengths);
return;


    
    