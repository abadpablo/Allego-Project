function timeInROIs = bz_timeSpentROIs(stem,left,right,center,varargin)

p = inputParser();

addParameter(p,'Fs',30,@isnumeric);
parse(p,varargin{:});
Fs = p.Results.Fs;


stem_len = length(stem);
left_len = length(left);
right_len = length(right);
center_len = length(center);

time_stem = stem_len/Fs;
time_left = left_len/Fs;
time_right = right_len/Fs;
time_center = center_len/Fs;

time = [time_stem time_left time_right time_center];
c = categorical({'StemArm', 'LeftArm', 'RightArm', 'Center'});

bar(c,time), ylabel('Seconds (s) '),title('Time spent in ROIs');

timeInROIs = {time_stem,time_left, time_right, time_center};
end

