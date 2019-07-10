function ecoupling_convert(ses_dir, bids_dir, sub, ses)
%ECOUPLING_CONVERT   Convert eCoupling data to BIDS format.
%
%  ecoupling_convert(ses_dir, bids_dir, sub, ses)
%
%  INPUTS
%  ses_dir - char
%      Path to directory with behavioral data for a session of eCoupling.
%      Must contain file named header.mat.
%
%  bids_dir - char
%      Path to the main directory where BIDS file structure should be
%      created.
%
%  sub - char
%      Subject code.
%
%  ses - char
%      Session code.

% load data from all runs
hdr_file = fullfile(ses_dir, 'header.mat');
[data, time] = ecoupling_load_subj(hdr_file);

% define fields with the sync square timing
% on: white to black (TTL 0)
% off: black to white (TTL 1)
sync = struct();
sync.know = {'fam_onset' 'fix_onset'};
sync.prac_disp = {'stim_onset' 'fix_onset'};
sync.prac_study = {'item1_onset' 'fix_onset'};
sync.prac_test = {'test_onset' 'fix_onset'};
sync.prac_match = {'cue_onset' 'fix_onset'};
sync.disp = {'stim_onset' 'fix_onset'};
sync.study = {'item1_onset' 'fix_onset'};
sync.test = {'test_onset' 'fix_onset'};
sync.match = {'cue_onset' 'fix_onset'};

rt = struct();
rt.know = {'fam_rt' 'fam'};
rt.prac_disp = {'rt' 'resp'};
rt.prac_study = {};
rt.prac_test = {'rt' 'resp'};
rt.prac_match = {'rt' 'resp'};
rt.disp = {'rt' 'resp'};
rt.study = {};
rt.test = {'rt' 'resp'};
rt.match = {'rt' 'resp'};

% fields to exclude if they already exist
exclude = {'fam_rt' 'fam' 'rt' 'resp' ...
           'onset' 'duration' 'trial_type' 'response_time' 'response'};

f = fieldnames(data);
for i = 1:length(f)
    d = data.(f{i});
    if isempty(d)
        continue
    end
    
    fp = sync.(f{i});
    fr = rt.(f{i});
    
    c_task = regexp(f{i}, '_', 'split');
    task = c_task{end};
    task_full = f{i};
    
    if ismember('run', d.Properties.VariableNames)
        run = d.run;
    else
        run = ones(size(d, 1), 1);
    end

    % extract required fields
    urun = unique(run);
    onset = [];
    duration = [];
    trial_type = repmat({task}, [size(d, 1) 1]);
    response_time = [];
    response = [];
    offset = [];
    for j = 1:length(urun)
        dr = d(run==urun(j),:);
        start = dr.(fp{1});
        finish = dr.(fp{2});
        
        onset = [onset; start];
        duration = [duration; finish - start];
        offset = [offset; finish];
        if ~isempty(fr)
            response_time = [response_time; dr.(fr{1})];
            response = [response; dr.(fr{2})];
        end
    end

    % exclude fields that are cells or are included as standard fields
    f2 = d.Properties.VariableNames;
    include = true(1, length(f2));
    for j = 1:length(f2)
        if ~isnumeric(d.(f2{j}))
            include(j) = false;
        end
        if ismember(f2{j}, exclude)
            include(j) = false;
        end
    end
    
    % combine new fields and included old fields
    if ~isempty(response)
        t_new = table(onset, duration, trial_type, response, response_time);
    else
        t_new = table(onset, duration, trial_type);
    end
    t_full = [t_new d(:,f2(include))];

    full_onset = onset;
    full_offset = offset;
    for j = 1:length(urun)
        t_run = t_full(run==urun(j),:);
        
        % construct standard filepath
        task_dir = fullfile(bids_dir, ['sub-' sub], ['ses-' ses], 'func');
        if ~exist(task_dir, 'dir')
            mkdir(task_dir)
        end
        filebase = sprintf('sub-%s_ses-%s_task-%s_run-%02d', ...
                           sub, ses, task_full, urun(j));
        
        % write out events
        filepath = fullfile(task_dir, [filebase '_events.tsv']);
        writetable(t_run, filepath, 'Delimiter', '\t', 'FileType', 'text');

        % get sync information in order
        run_onset = full_onset(run==urun(j));
        run_offset = full_offset(run==urun(j));
        onset = [run_onset; run_offset];
        signal = [zeros(size(run_onset)); ones(size(run_offset))];
        
        % sort in chronological order
        [~, ind] = sort(onset);
        onset = onset(ind);
        signal = signal(ind);

        % write out send signals
        t_send = table(onset, signal);
        filepath = fullfile(task_dir, [filebase '_send.tsv']);
        writetable(t_send, filepath, 'Delimiter', '\t', 'FileType', 'text');
    end
end
