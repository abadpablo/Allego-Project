%(c) Jorge R. Brotons-Mas 
% Notch filter for 50 and 100 Hz
% Uses a notch filter for 50 hz or 100 Hz depending on the number of inputs
% in the functuon 
% CA1=Notched(CA1,Args.Fs,50,100); % NOTCH FILTER FOR 50 an 100 Hz Arguments (Signal, Frequency 1 , Frequency 2)

function [dataN]=Notched(data,Args,Fr1,Fr2,Fr3);
if nargin<3
    disp('error not enough arguments, check for FS etc...')
elseif nargin==3
%     CA1 =filter(b,a,CA1);
    [Hd b a] = NotchIrFir50 (Args.Fs);
    dataN=filter(b,a,data);
elseif nargin==4
    [Hd b a] = NotchIrFir50(Args.Fs) ;
    dataN=filter(b,a,data);

    [Hd b a] = NotchIrFir100(Args.Fs);
    dataN=filter(b,a,dataN);

elseif nargin == 5
    [Hd b a] = NotchIrFir50(Args.Fs) ;
    dataN=filter(b,a,data);

    [Hd b a] = NotchIrFir100(Args.Fs);
    dataN=filter(b,a,dataN);
    
    [Hd b a] = NotchIrFir150(Args.Fs);
    dataN=filter(b,a,dataN);
    
end
