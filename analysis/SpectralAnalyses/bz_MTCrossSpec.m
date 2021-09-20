function [Sc,Cmat,Ctot,Cvec,Cent,f] = bz_MTCrossSpec(lfp,varargin)
% MTCrossSpec - Compute Multi-taper cross-spectral matrix
%
%   USAGE
%
%   [Sc,Cmat,Ctot,Cvec,Cent,f] = bz_MTCrossSpec(lfp,varargin)
%
%   INPUTS
%   lfp - buzcode lfp data struct
%   <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   sampling rate (in Hz) (default = from timestamps if
%                   available, otherwise 1250Hz)
%     'range'       frequency range (in Hz) (default = all)
%     'window'      duration (in s) of the time window (default = 5)
%     'overlap'     overlap between successive windows (default = window/2)
%     'step'        step between successive windows (default = window/2)
%     'tapers'      relative resolution and order of the tapers [NW K]
%                   (default = [3 5])
%     'pad'         FFT padding (see help for <a href="matlab:help cohgramc">cohgramc</a>) (default = 0)
%     'show'        plot results (default = 'off')
%     'cutoffs'     cutoff values for color plot (default = [0 1])
%    =========================================================================
%
%   OUTPUT
%   Sc - Cross Spectral matrix frequency x channels x channels
%   Cmat - Coherence matrix frequency x channels x channels
%   Ctot - Total coherence : SV(1)^2/sum(SV^2)
%   Cvec - Leading eigenvector
%   Cent - A different measure of total coherence : GM/AM of SV^2s
%   f - frequencies
%
%
%
%  DEPENDENCIES
%
%    This function requires the <a href="http://www.chronux.org">chronux</a> toolbox.
%
% Created by Pablo Abad 2021
%

% Make sure Chronux is installed
CheckChronux('CrossSpecMatc')

% Default Params
p = inputParser();
addParameter(p,'frequency',1250,@isnumeric);
addParameter(p,'range',[],@isnumeric);
addParameter(p,'window',5,@isnumeric);
addParameter(p,'overlap',2.5,@isnumeric);
addParameter(p,'step',2.5,@isnumeric);
addParameter(p,'tapers',[3 5], @isnumeric);
addParameter(p,'pad',0,@isnumeric);
addParameter(p,'show','on',@isstr);
addParameter(p,'cutoffs',[0 1], @isnumeric);

parse(p,varargin{:});

frequency = p.Results.frequency;
range = p.Results.range;
window = p.Results.window;
overlap = p.Results.overlap;
step = p.Results.step;
tapers = p.Results.tapers;
pad = p.Results.pad;
show = p.Results.show;
cutoffs = p.Results.cutoffs;

% Define params for chornux functions
params.Fs = frequency;
if ~isempty(range)
    params.fpass = range;
end
params.tapers = tapers;
params.pad = pad;

[Sc,Cmat,Ctot,Cvec,Cent,f] = CrossSpecMatc(double(lfp.data),[window window-overlap], params);

if strcmpi(show,'on')
    figure,
    
    
    
    
end

end

