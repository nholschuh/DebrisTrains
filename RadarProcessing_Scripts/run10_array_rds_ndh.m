% script run_array
%
% Script for running array (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_array.m, array.m, load_sar_data.m,
% array_proc.m, array_task.m, array_combine_task.m

%% User Setup
% =====================================================================
param_override = [];

if exist('param_ssheet_name') == 0
    param_ssheet_name = 'rds_param_2011_Greenland_P3.xlsx';
end
if exist('cluster_type') == 0
    cluster_type = 'slurm';
    all_ping_pong = 0;
    rds_settings_flag = 0;
end


%%%%%%%%%%%%%%%%%%% This loop handles multiple segids selected
params = read_param_xls(ct_filename_param(param_ssheet_name),'');
if exist('seg_ids') ~= 0
    params = ct_set_params(params,'cmd.array',0);
    seg_opts = {params.day_seg};
    for ndh_loop = 1:length(seg_ids)
        [segidflag, activate_ind] = strcmp_ndh(seg_opts,seg_ids{ndh_loop});
        if segidflag == 1
            params(activate_ind).cmd.array = 1;
            params(activate_ind).cmd.frms = eval(frame_inds{ndh_loop});


            %%%%%%%%%%%%%%%%%%%%%% This handles just left or just right
            if all_ping_pong == 1
                for img_ind = 1:length(params(activate_ind).array.imgs)
                    half = length(params(activate_ind).array.imgs{img_ind}(:,1))/2;
                    params(activate_ind).array.imgs{img_ind} = params(activate_ind).array.imgs{img_ind}(1:half,:);
                end
            elseif all_ping_pong == 2
                for img_ind = 1:length(params(activate_ind).array.imgs)
                    half = length(params(activate_ind).array.imgs{img_ind}(:,1))/2;
                    params(activate_ind).array.imgs{img_ind} = params(activate_ind).array.imgs{img_ind}(half+1:end,:);
                end
            end
        end
    end
else
  params = ct_set_params(params,'cmd.array',0);
end

%%%%%%%%%%%%%%%%%%% Here we set the appropriate out path
param_override.array.in_path = 'sar';

% param_override.array.out_path = 'standard_ndh';
% param_override.array.out_path = 'standard_tukey_ndh';
% param_override.array.out_path = 'mvdr_ndh';
param_override.array.out_path = 'music3D_ndh';
% param_override.array.out_path = 'standard_phase_ndh';




if param_override.array.out_path(1) == 's';

    if length(param_override.array.out_path) > 15
        param_override.array.tomo_en = 0;
        params = ct_set_params(params,'array.method','standardphase');
        params = ct_set_params(params,'array.bin_rng',0);
        params = ct_set_params(params,'array.line_rng',0);
        params = ct_set_params(params,'array.dbin',1);
        params = ct_set_params(params,'array.dline',1);
        params = ct_set_params(params,'array.ft_over_sample',1);

    else
        param_override.array.tomo_en = 0;
        param_override.array.method = 'standard';
        param_override.array.out_path = 'standard_ndh';
        if noise_flag == 2
            param_override.array.in_path = 'sar_ndh_burstremoved';
            param_override.array.out_path = 'standard_ndh_burstremoved';
        end
    end


elseif param_override.array.out_path(1) == 'm';
    param_override.array.tomo_en = 1;
    param_override.array.Nsrc = 3;
    param_override.array.Nsv = 64;
    param_override.array.method = 'music';
    param_override.array.out_path = 'music3D_ndh_v2';
    if exist('noise_flag') == 1
        if noise_flag == 2
            param_override.array.in_path = 'sar_ndh_burstremoved';
            param_override.array.out_path = 'music3D_ndh_burstremoved';
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%% This handles just left or just right
if all_ping_pong == 1
    param_override.array.out_path = [param_override.array.out_path,'_ping'];
elseif all_ping_pong == 2
    param_override.array.out_path = [param_override.array.out_path,'_pong'];
end



if exist('narrowband_flag')
    if narrowband_flag == 1
        nb_adjust
        if length(params.array.in_path) > 0
            params = ct_set_params(params,'array.in_path',['nb_',params.sar.out_path]);
        else
            params = ct_set_params(params,'array.in_path',['nb_sar']);
        end
        param_override.array.out_path = ['nb_',param_override.array.out_path];
    end
end


cmd_method = 'array';
if rds_settings_flag == 1
    rds_settings_ndh;
end

%dbstop if error;
param_override.cluster.type = cluster_type;
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 10*60;
% param_override.cluster.cpu_time_mult  = 4;
% param_override.cluster.mem_mult  = 2;
% param_override.cluster.max_jobs_active       = 1;
% param_override.cluster.qsub_submit_arguments = '-q debug -m n -l nodes=1:ppn=1:dcwan:dc2,pmem=%dmb,walltime=%d:00';
param_override.cluster.max_time_per_job = 172800*5;
param_override.cluster.cpu_time_mult  = 10;
param_override.cluster.mem_mult  = 10;
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

% Process each of the segments
ctrl_chain = {};
for param_idx = 1:length(params)
    param = params(param_idx);
    if param.cmd.array
        ctrl_chain{end+1} = array(param,param_override);
    end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
