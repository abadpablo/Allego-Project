function [skaggs]= bz_SkaggsIndex (map,varargin)

% Skaggs for place cells, modifiyed from UCL
% Adapted from SkaggsIndex by Pablo Abad 23/07 to resemble buzCode
% 
% USAGE
%   skaggs = bz_SkaggsIndex(map,<options>);
%
% INPUT
%   map         map obtained using <a href="matlab:help Map">Map</a>
%   options     optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'basepath'    By default pwd
%     'minTime'     minimum time in a bin to be included for analysis(default = 0)
%     'unSmooth'    including or not the analysis of the unsmoothed firing
%                   map (default true)
%     'duration'    duration of the recording (default obtained by
%                   sessionInfo)
%     'verbose'     display processing information (default = 'off')
%    =========================================================================
% OUTPUT
% Skaggs values 
%
% Pablo Abad 07/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'unSmooth',true,@islogical);
addParameter(p,'duration',[],@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
minTime = p.Results.minTime;
unSmooth = p.Results.unSmooth;
duration = p.Results.duration;

%% Skaggs Modulated by total time
nanmask=(map.z>0)|(map.time>minTime);
duration = sum(map.time(nanmask)); % duration only for the speed filtered map
meanRate = sum(map.z(nanmask).*map.time(nanmask))/duration;

if meanRate<=0.0
    fprintf(' mean rate is %f \n', meanRate);
    skaggs.bits_per_sec = NaN;
    skaggs.bits_per_spike = NaN;
    return;
end

p_x = map.time ./ duration;
p_r = map.z ./ meanRate;
% avoid finding 0.log(0) terms in sum
dummy = p_x .* map.z;
idx = dummy > 0;
skaggs.bitsPerSec = sum(dummy(idx).*log2(p_r(idx)));
skaggs.bitsPerSpike = skaggs.bitsPerSec/meanRate;

%% Unsmoothed
if unSmooth
    nanmask_Uns=(map.zUnSmooth>0)|(map.time>minTime);
    duration_Uns = sum(map.time(nanmask)); % duration only for the speed filtered map
    meanRate_Uns = sum(map.zUnSmooth(nanmask_Uns).*map.time(nanmask_Uns))/duration_Uns;

    if meanRate_Uns<=0.0
        fprintf(' mean rate is %f \n', meanRate);
        skaggs.bits_per_sec_Uns = NaN;
        skaggs.bits_per_spike_Uns = NaN;
        return;
    end

    p_x_Uns = map.time ./ duration_Uns;
    p_r_Uns = map.zUnSmooth ./ meanRate_Uns;
    % avoid finding 0.log(0) terms in sum
    dummy_Uns = p_x_Uns .* map.zUnSmooth;
    idx_Uns = dummy_Uns > 0;
    skaggs.bitsPerSec_Uns = sum(dummy_Uns(idx_Uns).*log2(p_r_Uns(idx_Uns)));
    skaggs.bitsPerSpike_Uns = skaggs.bitsPerSec_Uns/meanRate_Uns;
end


