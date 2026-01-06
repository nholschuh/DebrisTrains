% script run_array
%
% Script for running array (usually just used for debugging).
%
% Authors: Nick Holschuh and John Paden
%
% See also: run_master.m, master.m, run_array.m, array.m, load_sar_data.m,
% array_proc.m, array_task.m, array_combine_task.m

%% User Setup
% =====================================================================



param_ssheet_name = 'accum_param_2023_Antarctica_BaslerMKB.xlsx';
params = read_param_xls(ct_filename_param(param_ssheet_name),'');

cluster_type = 'slurm';

% Process each of the segments
ctrl_chain = {};

range_limit = [4,6,8,10,12,14,16,22];


for rl_ind = 1:length(range_limit)
    rl = range_limit(rl_ind);


    param_override = [];

    %%%%%%%%%%%%%%%%%%% Here we set the appropriate out path
    param_override.array.in_path = 'sar';
    param_override.array.tomo_en = 1;
    param_override.array.Nsrc = 3;
    param_override.array.Nsv = 128;
    param_override.array.method = 'music';
    param_override.array.out_path = ['music3D_ndh_aw_',num2str(rl)];
    param_override.array.imgs = {[2*ones([rl 1]),(1:rl).']};
    param_override.array.img_comb = [];


    %dbstop if error;
    param_override.cluster.type = cluster_type;
    % param_override.cluster.rerun_only = true;
    % param_override.cluster.desired_time_per_job  = 10*60;
    % param_override.cluster.cpu_time_mult  = 4;
    % param_override.cluster.mem_mult  = 2;
    % param_override.cluster.max_jobs_active       = 1;
    % param_override.cluster.qsub_submit_arguments = '-q debug -m n -l nodes=1:ppn=1:dcwan:dc2,pmem=%dmb,walltime=%d:00';
    param_override.cluster.max_time_per_job = 172800*5;
    param_override.cluster.cpu_time_mult  = 3;
    param_override.cluster.mem_mult  = 3;
    param_override.cluster.mem_mult_mode = 'auto';
    param_override.cluster.max_cpu_time_mode = 'truncate';
    param_override.cluster.max_mem_mode = 'truncate';
    param_override.cluster.stop_on_error = 0;

    %% Automated Section
    % =====================================================================

    % Input checking
    global gRadar;
    if exist('param_override','var')
        param_override = merge_structs(gRadar,param_override);
    else
        param_override = gRadar;
    end


    for param_idx = 1:length(params)
        param = params(param_idx);
        if param.cmd.array
            ctrl_chain{end+1} = array(param,param_override);
        end
    end
end
cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
