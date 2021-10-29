function [ax1,ax2] = bz_plotStability(spikes,id,varargin)
%   Function adapted from plot_stabilityV2 from Jorge Brotons to OFS and
%   tint data
%   
%   Adapted from Pablo Abad to meet buzcode requirements
%   plot_stability - Plot amplitude and firing rate of a set of spikes over time
%
%   USAGE
%   [ax1,ax2] = bz_plotStability(spiketimes)
%
%   Description:  
%   Creates two axes in the same position to simultaneously plot two 
%   aspects of the stability of a given set of spikes.  In the upper half
%   and using the right y-axis is a scatter plot of the peak-to-peak amplitude 
%   of waveforms taken from the full duration of the recording.  For efficiency
%   a maximum number of randomly drawn data points is used as set by
%   spikes.params.display.max_scatter.  In the lower half and using the left
%   y-axis is a histogram of the firing rate as a function of time.  The 
%   unwrapped time is used (see unwrap_time.m).  The size of the bins used
%   is set by spikes.params.display.stability_bin_size.
%   
%
%   
%
%
%
%
%
% Default params 
p = inputParser;

addParameter(p,'binsize',1,@isnumeric);
addParameter(p,'numwaves',600,@isnumeric);
addParameter(p,'timestamps',[],@isnumeric);

parse(p,varargin{:})
binsize = p.Results.binsize;
numwaves = p.Results.numwaves;
timestamps = p.Results.timestamps;

cla
ax1 = gca;

% We re ploting two axes on top of each other
% The first axes will show the firing rate as a function of time
TotTime = timestamps(2)-timestamps(1);

% bin data
tlims = [0 TotTime];
tlims = timestamps;
num_bins = round(diff(tlims) / binsize);
edges = linspace(tlims(1),tlims(2),num_bins+1);
spiketimes = spikes.times{id};
[status,interval,index] = InIntervals(spiketimes,timestamps);
spiketimes = spiketimes(status);
n = histc(spiketimes,edges);
n(end) = [];
vals = n/mean(diff(edges));

% show histogram
hold on
% bar(edges(1:end-1) + mean(edges(1:2)),vals,1.0); shading flat;
bar(edges(1:end-1), vals, 1.0)
hold off
set(gca,'XLim',tlims,'YLim',[0 2*max(get(gca,'YLim'))])
yticks = get(gca,'YTick'); set(gca,'YTick',yticks(yticks <=max(yticks)/2))
xlabel('Time(sec)')
ylabel('Firing rate (Hz)')


%% Second pass is scatter plot of amplitude versus time
% Get amplitude over time
amp = spikes.amplitudes{id};
amp = amp(status == 1);

% prettify axes
ax2 = axes('Parent',get(ax1,'Parent'),'Unit',get(ax1,'Unit'),'Position',get(ax1,'Position'));
hold on
% l = scatter(spiketimes(ind),amp(ind));
l = scatter(spiketimes,amp);
hold off
l1 = line( get(ax2,'XLim'), [ 0 0] );
set(l1,'Color',[ 0 0 0])
hold off
set(l,'Marker','.','MarkerEdgeColor',[.3 .5 .3])
set(ax2,'Xlim',tlims,'YLim',max( get(ax2,'YLim') ) *[-1 1],'XTickLabel',{},'YAxisLocation','right');
yticks = get(ax2,'YTick'); set(ax2,'YTick',yticks(yticks>=0))
ylabel('Amplitude')
set(ax2,'Color','none')

% linke properties of the 2 axes together
linkaxes([ax1 ax2],'x');
linkprop([ax1 ax2],'Position');
linkprop([ax1 ax2],'Unit');

end

