function [data, time] = ecoupling_load_subj(hdr_file)
%LOAD_SUBJ_DATA   Load data for all experiment phases.
%
%  [data, time] = ecoupling_load_subj(hdr_file)

subj_dir = fileparts(hdr_file);

load(hdr_file);
periods = hdr.par.periods;
data = struct();
time = [];
for i = 1:length(periods)
    data.(periods{i}) = [];
    mat_files = hdr.data.(periods{i}).mat;
    for j = 1:length(mat_files)
        [parent, filename, ext] = fileparts(mat_files{j});
        [orig_dir, sub_dir] = fileparts(parent);
        filepath = fullfile(subj_dir, sub_dir, [filename ext]);
        
        if ~exist(filepath, 'file')
            continue
        end
        s = load(filepath);
        data.(periods{i}) = [data.(periods{i}); s.data];
        period = periods(i);
        run = j;
        start = {datestr(s.time.start_time, 30)};
        finish = {datestr(s.time.finish_time, 30)};
        duration = s.time.duration;
        r = table(period, run, duration, start, finish);
        time = [time; r];
    end
end
