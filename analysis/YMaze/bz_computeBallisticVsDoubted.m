function [ballistic, ballistic_perc, doubted, doubted_perc] = bz_computeBallisticVsDoubted(tran)
global filename;
global filepath;
global files2analyse;
global pathFigures;
% Let's compute when the animal performs a ballistic trial (when it goes
% direct from one ROI to another one and when the animal performs a doubted
% trial (when it goes out of the ROI and then it enters again, either it is
% the center of the Maze or if it is one arm

doubted_trials = {'1010','2020','3030','4040'}
win = 4;
counter_doubted = 0;
counter_ballistic = 0;

for i=1:2:length(tran)-4
    
    if i==1
        epoch = tran(1:4)
    else
        epoch = tran(i:i+win-1)
    end
    
    
    cj{i} = contains(doubted_trials,epoch);
    
end

for i=1:length(cj)
    
    if ~isempty(cj{i})
    
        if find(cj{i} == 1)

            counter_doubted = counter_doubted + 1;
        else
            counter_ballistic = counter_ballistic + 1;

        end
    end
    
end


ballistic = counter_ballistic;
doubted = counter_doubted;

total = ballistic + doubted;

ballistic_perc = (ballistic/total)*100;
doubted_perc = (doubted/total)*100;

figure,
c = categorical({'Ballistic %', 'Doubted %'})
bar(c, [ballistic_perc doubted_perc])
title(['Ballistic ', num2str(ballistic_perc), '%', '   ' ,'Doubted ' ,num2str(doubted_perc),'%'])

saveas(gcf,'Ballistic vs Doubted.png');

end
