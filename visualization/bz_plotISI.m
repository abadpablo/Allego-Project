function [] = bz_plotISI(spikes,id,varargin)

% Default Params
p = inputParser;

addParameter(p,'timestamps',[],@isnumeric)
addParameter(p,'binsize',1,@isnmeric);
addParameter(p,'isi_maxlag',[],@isnumeric);
addParameter(p,'shadow',1,@isnumeric);
addParameter(p,'refractory_period',1,@isnumeric)
addParameter(p,'trialDuration',[],@isnumeric)
addParameter(p,'show_isi',1,@isnumeric)

parse(p,varargin{:})

timestamps = p.Results.timestamps;
bin_size = p.Results.binsize;
isi_maxlag = p.Results.isi_maxlag;
shadow = p.Results.shadow;
refractory_period = p.Results.refractory_period;
trialDuration = p.Results.trialDuration;
show_isi = p.Results.show_isi;
trialDuration = timestamps(2)-timestamps(1);

spiketimes = spikes.times{id};
[status,interval,index] = InIntervals(spiketimes,timestamps);
spiketimes = spiketimes(status);


data.spiketimes = spiketimes;
data.show_isi = show_isi;
data.isi_maxlag = isi_maxlag;
data.autocorr_maxlag = data.isi_maxlag;
data.shadow = shadow;
data.refractory_period = refractory_period;
data.corr_bin_size = bin_size;
data.isi_bin_size = bin_size;
if isi_maxlag >=1
    meanFr = length(spiketimes)/trialDuration;
    isis = diff(spiketimes);
    maxlag = isi_maxlag;
    bins = bin_size;
    [n,x] = histc(isis*1000,linspace(0,1000*maxlag,bins));
    
    if isi_maxlag >0.5
        intervalNumber = length(spiketimes)-1;
        [MaxIsh PosMax] = max(n);
        rp = refractory_period;
        bins = round(1000*maxlag/bin_size);

        % make plot 
        isis = diff(spiketimes);
        isis = isis(isis<=maxlag);
        [n,x] = histc(isis*1000,linspace(0,1000*maxlag,bins));

        ymax = max(n)+1;
        values = n(1:10);
        Values1 = sum(values1);
        ProbBurst = Values1*100/length(spiketimes-1);
        ll = length(spiketimes);
        dif_max = 0.015;
        dif_spk = diff(spiketimes);
        no_burst = find(dif_spk > dif_max);
        count = 0;
        ind=1;
        if  (length(no_burst)+1~=length(spiketimes))
            aa=dif_spk;
            aa(no_burst)=0;
            pp=find(aa);
            if (no_burst(ind)-pp(1) >= 2)
                burst_spk=spiketimes(pp(1):no_burst(ind));
                cont=cont+1;
                reg_burst(cont).nspk=length(burst_spk);
                reg_burst(cont).isi=mean(diff(burst_spk));
                reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                reg_burst(cont).burst=burst_spk;
                ind=ind+1; 
            end
            while (ind < length(no_burst))
              if (no_burst(ind+1)-no_burst(ind) > 1) %% este >1 no tiene que ver con el numero de spikes
                 burst_spk=spiketimes(no_burst(ind)+1:no_burst(ind+1));
                 if (length(burst_spk)>=2)
                     cont=cont+1;
                     reg_burst(cont).nspk=length(burst_spk);
                     reg_burst(cont).isi=mean(diff(burst_spk));
                     reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                     reg_burst(cont).burst=burst_spk;
                 end
                 ind=ind+1;
              else
                 ind=ind+1;
              end
            end
            if (ind == length(no_burst))
                if (no_burst(ind) < length(spiketimes)-1)
                 burst_spk=spiketimes(no_burst(ind)+1:length(spiketimes));
                 cont=cont+1;
                 reg_burst(cont).nspk=length(burst_spk);
                 reg_burst(cont).isi=mean(diff(burst_spk));
                 reg_burst(cont).duration=burst_spk(end)-burst_spk(1);
                 reg_burst(cont).burst=burst_spk;
                end
            end
            %%%%% Parte que calcula el IBI:
            if (cont>2)
                for jj=1:cont-1
                    t1=reg_burst(jj).burst(end);
                    t2=reg_burst(jj+1).burst(1);
                    dif_ibi(jj)=t2-t1;
                end
            else
                dif_ibi=[];
            end
            %%%%%%% el promedio del numero de spikes por burst:
              kk=[reg_burst.nspk];
            %   sprintf('numero promedio de spikes: %2.5g', mean(kk))
              n_spike= mean(kk);
            %%%% el promedio de ISI en cada burst es:
              kk1=[reg_burst.isi];
            %   sprintf('el promedio de ISI por el registro: %2.5g ms', mean(kk1)*1000)
              prom_isi=mean(kk1)*1000;
            %%%% la duracion:
              kk2=[reg_burst.duration];
            %   sprintf('la duracion promedio: %2.5g ms', mean(kk2)*1000)
              prom_duracion=mean(kk2)*1000;
            %%%% el IBI promedio:
            %   sprintf('el IBI promedio: %2.5g s', mean(dif_ibi))
              if isempty(dif_ibi)
                  prom_ibi=0;
              else
                  prom_ibi=mean(dif_ibi);
              end
        % numero de bursts
          num_burst=cont;
        else
            %   sprintf('%s ','neurona no es burster')
              num_burst=cont;
              n_spike= 0;
              prom_isi=0;
              prom_duracion=0;
              prom_ibi=0;
        end
    numbursrNorm=(num_burst*100)/length(spiketimes);
    % XCorr
    [cross,lags] = pxcorr(spiketimes,spiketimes, round(1000/data.corr_bin_size), maxlag);
    autoco=cross; 
    MeanAutoc=mean(autoco(1001:end));
    trough1=mean(autoco(1051:1071));
    peak1=mean(autoco(1101:1141));
    TMI=abs(trough1-peak1)/(trough1+peak1);

    end
end
        
%% estimate poisson contamination for y-axis
[expected,lb,ub,RPV] = ss_rpv_contaminationV2(spiketimes,spiketimes,refractory_period,trialDuration,shadow);
    
if isempty(lb)
    data.ystr = [num2str(RPV) ' RPVs (' num2str(round(expected*100)) '%)' ];
else    
    data.ystr = [num2str(RPV) ' RPVs (' num2str(round(lb*100)) '-' num2str(round(expected*100)) '-' num2str(round(ub*100)) '%)' ];
end
  
set(gca,'UserData', data,'Tag','isi' )
% write updating method
update_isi( [], [], show_isi, gca);   
    

% callback that updates display on axes
function update_isi( hObject, event, displaymode, ax)

        % get display mode
        data = get( ax,'UserData');
        data.show_isi = displaymode;
        set(ax,'UserData',data)
        set(gcf,'CurrentAxes',ax )
        
        % update display
        make_isi;
    
        % set context menu - allow switch between ISI and autocorrelation
        % as well as a global switch for all plot_isi instances on the current figure
        cmenu = uicontextmenu;
        item(1) = uimenu(cmenu, 'Label', 'Show ISI', 'Callback', {@update_isi, 1,ax} );
        item(2) = uimenu(cmenu, 'Label', 'Show autocorr', 'Callback',{@update_isi, 0,ax}  );
        item(3) = uimenu(cmenu, 'Label', 'Use this style on all ISIs in figure', 'Callback',{@impose_all,displaymode,ax},'Separator','on'  );
        set(item(2-displaymode), 'Checked', 'on');    
        set(ax,'UIContextMenu', cmenu )

end

% callback to impose display mode on all plot_isi axes in this figure
function impose_all(hObject,event,displaymode,ax)
        [o,h] = gcbo;
        my_axes = findobj(h,'Tag','isi');
        my_axes = setdiff(my_axes,ax);
        for j = 1:length(my_axes), update_isi([],[],displaymode,my_axes(j)); end
end    
    
% plots the ISI or autocorrelation
function make_isi     

    data = get( gca,'UserData');
    spiketimes = data.spiketimes;
    shadow = data.shadow;
    rp     = data.refractory_period;
    cla reset
    
    % ISI case
    if data.show_isi
        maxlag = data.isi_maxlag;
        bins = round(1000* maxlag/data.isi_bin_size );

        % make plot
        isis = diff(spiketimes);   
        isis = isis(isis <= maxlag); 
        [n,x] =hist(isis*1000,linspace(0,1000*maxlag,bins));
        ymax   = max(n)+1;
%       values1=n(1:10);
% Values1=sum(values1);
% ProbBurst=Values1*100/length(spiketimes-1);
%%        
        
        % make patches to represent shadow and refractory period

        patch([0 shadow shadow 0 ], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none');
        patch([shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        hold on,    b2 = bar(x,n,1.0); hold off
        set(b2,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])

        % update axes
        set(gca,'YLim',[0 ymax],'XLim',[0 1000*maxlag])
        xlabel('Interspike interval (msec)')
       ylabel({'No. of spikes',data.ystr})
 
    else
        maxlag = data.autocorr_maxlag;
        
        % calculate autocorrelation
        if length(spiketimes) > 1
          [cross,lags] = pxcorr(spiketimes,spiketimes, round(1000/data.corr_bin_size), maxlag);
        else
            cross = 0;  lags = 0;
        end
        cross(find(lags==0)) = 0;
        
        % place patches to represent shadow and refractory period
        ymax = max(cross) + 1;
        patch(shadow*[-1 1 1 -1], [0 0 ymax ymax], [.5 .5 .5],'EdgeColor', 'none');
        patch([shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        patch(-[shadow [rp rp] shadow ], [0 0 ymax ymax], [1 0 0],'EdgeColor','none');
        
        % plot autocorrelation histogram
        hold on, bb = bar(lags*1000,cross,1.0); hold off;  
        set(bb,'FaceColor',[0 0 0 ],'EdgeColor',[0 0 0 ])

        % set axes
        set(gca, 'XLim', maxlag*1000*[-1 1]);
        set(gca,'YLim',[0 ymax])
        xlabel( 'Time lag (msec)')
        ylabel({'Autocorrelation (Hz)',data.ystr})
        
    end
    set(gca,'Tag','isi','UserData',data)
     
end

end

