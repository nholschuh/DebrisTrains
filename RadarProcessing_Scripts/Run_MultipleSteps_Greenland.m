%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run_MultipleSteps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rds_settings_flag = 0;


%%cluster_type = 'debug';
cluster_type = 'slurm';
noise_flag = 0;
narrowband_flag = 0;


%%%%%%%%%%% PetermannLines
seg_id_list= {'20110429_01','20110429_02','20110507_01','20110507_02'};
frame_id_list = {'[]','[]','[]','[]'};  

process_steps = [0 1];
%%%%%%%%%%%%%%%% 0 1 2 3 4 5 6 7 8 9  10
%%%%%%%%%%%%%%%% 1 2 3 4 5 6 7 8 9 10 11
% 0/1 - Sar Processing
% 1/2 - Array Processing



for kkkk = 1:length(seg_id_list)
    seg_ids = seg_id_list{kkkk};
    if isstr(frame_id_list{kkkk})
        frame_ids = eval(frame_id_list{kkkk});
    else
        frame_ids = frame_id_list{kkkk};
    end
    param_ssheet_name = ['rds_param_',seg_id_list{kkkk}(1:4),'_Greenland_P3.xls'];
    
    chain_num_start = max(cluster_get_chain_list);
    if isempty(chain_num_start)
        chain_num_start = 0;
    end
    chain_links = [];
    counter = 0;
    
    if rds_settings_flag == 1
        copyfile(source_file,'rds_settings_ndh.m')
    end
    
    %%
    %clearvars -except stop* noise_flag process_steps counter chain* cluster* source_file
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% 9 -- sar process
    if process_steps(10) == 1
        counter = counter+1;
        run09_sar_rds_ndh
        [chain_links(counter).chain,chain_links(counter).chain_fn_links] = cluster_load_chain(chain_num_start+counter);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% 10 -- array process
    if process_steps(11) == 1
        counter = counter+1;
        run10_array_rds_ndh
        [chain_links(counter).chain,chain_links(counter).chain_fn_links] = cluster_load_chain(chain_num_start+counter);
    end
    
    %% Intermediate run step!!!0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ctrl_chain = cell(0);
    for i = 1:length(chain_links)
        if length(ctrl_chain) == 0
            ctrl_chain{1} = [chain_links(i).chain{1,:}];
        else
            ctrl_chain{1} = [ctrl_chain{1} chain_links(i).chain{1,:}];
        end
    end
    
    ctrl_chain = cluster_run(ctrl_chain);
    
    counter = 0;
    chain_links = [];
    ctrl_chain = {};
    chain_num_start = max(cluster_get_chain_list);
    if isempty(chain_num_start)
        chain_num_start = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
end
