% script run_sar
%
% Script for running sar (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_sar.m, sar.m, sar_task.m,
%   sar_coord_task.m

%% User Setup
% =====================================================================
param_override = [];

if exist('param_ssheet_name') == 0
  param_ssheet_name = 'rds_param_2011_Greenland_P3.xlsx'
if exist('cluster_type') == 0
  cluster_type = 'debug';
end

%%%%%%%%%%%%%%%%%%% This loop handles multiple segids selected
params = read_param_xls(ct_filename_param(param_ssheet_name),'');
if exist('seg_ids') ~= 0
  params = ct_set_params(params,'cmd.sar',0);
  seg_opts = {params.day_seg};
  for ndh_loop = 1:length(seg_ids)
      [segidflag, activate_ind] = strcmp_ndh(seg_opts,seg_ids{ndh_loop},1);
      if segidflag == 1
         params(activate_ind).cmd.sar = 1;
         params(activate_ind).cmd.frms = eval(frame_inds{ndh_loop});
      end
  end
else
  params = ct_set_params(params,'cmd.sar',1);
end

%%%%%%%%%%%%%%%%%%% Here we set the appropriate out path
params = ct_set_params(params,'sar.out_path',['sar_ndh']);

if exist('narrowband_flag')
  if narrowband_flag == 1
    nb_adjust
    if length(params.sar.out_path) > 0
      params = ct_set_params(params,'sar.out_path',['nb_',params.sar.out_path]);
    else
        params = ct_set_params(params,'sar.out_path',['nb_sar']);
    end
  end
end


cmd_method = 'sar';

if rds_settings_flag == 1 
  rds_settings_ndh;
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if strcmpi(cmd_method,'generic')
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
  elseif ~param.cmd.(cmd_method)
    continue;
  end
  

  %% radar.wfs
  for wf = 1:length(params(param_idx).radar.wfs)
    %params(param_idx).radar.wfs(wf).deconv.en = true;
    if noise_flag == 1
        params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
        params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold';
        %params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold_tukey';
    elseif noise_flag == 2
        params(param_idx).radar.wfs(wf).burst.en = true;
        params(param_idx).radar.wfs(wf).bad_value = NaN;
        params(param_idx).radar.wfs(wf).burst.fn = 'analysis_burst';
        params = ct_set_params(params,'sar.out_path',['sar_ndh_burstremoved']);
    elseif noise_flag == 3
        params(param_idx).radar.wfs(wf).burst.en = true;
        params(param_idx).radar.wfs(wf).bad_value = NaN;
        params(param_idx).radar.wfs(wf).burst.fn = 'analysis_burst_2';
        pparams = ct_set_params(params,'sar.out_path',['sar_ndh_burstremoved_2']);
    else
        params(param_idx).radar.wfs(wf).coh_noise_method = [];
        params(param_idx).radar.wfs(wf).coh_noise_arg.fn = [];
        params(param_idx).radar.wfs(wf).burst.en = false;
    end
  end
end

%dbstop if error;
param_override.cluster.type = cluster_type;
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 120*60;
param_override.cluster.max_time_per_job = 172800*1;
param_override.cluster.cpu_time_mult  = 100;
param_override.cluster.mem_mult  = 20;
param_override.cluster.mem_mult_mode = 'auto';
param_override.cluster.max_cpu_time_mode = 'truncate';
param_override.cluster.max_mem_mode = 'truncate';
param_override.cluster.stop_on_error = 0;
% param_override.cluster.mem_mult  = 2;
% param_override.cluster.max_jobs_active       = 1;
% param_override.cluster.qsub_submit_arguments = '-q debug -m n -l nodes=1:ppn=1:dcwan:dc2,pmem=%dmb,walltime=%d:00';

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
  if param.cmd.sar
    ctrl_chain{end+1} = sar(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
