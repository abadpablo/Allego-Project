function [periodic] = bz_Power2D(firingMaps,id,varargin)
%
%   Power Spectrum in the 2D firing map
%   Function based on Jorge Brotons Power2D
%   This function calculates the power spectrum of a firing matrix or in
%   general 2D data matrix
%   
%   USAGE
%   periodic = bz_Power2D(firingMaps,id,varargin)
%
%   INPUTS
%   firingMaps
%   id: id of the neuron
%
%
%
%
%
% Defaults Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'showFig',true,@islogical);

parse(p,varargin{:})
basepath = p.Results.basepath;
showFig = p.Results.showFig;




end

