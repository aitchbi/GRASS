function n = hb_get_hcp_task_length(task)
%
%
%
% March 2020
% Hamid Behjat


switch task
    case {'tfMRI_EMOTION_LR','tfMRI_EMOTION_RL','EMOTION'}
        n = 176;
    case {'tfMRI_GAMBLING_LR','tfMRI_GAMBLING_RL','GAMBLING'}
        n = 253;
    case {'tfMRI_LANGUAGE_LR','tfMRI_LANGUAGE_RL','LANGUAGE'}
        n = 316;
    case {'tfMRI_MOTOR_LR','tfMRI_MOTOR_RL','MOTOR'}
        n = 284;
    case {'tfMRI_RELATIONAL_LR','tfMRI_RELATIONAL_RL','RELATIONAL'}
        n = 232;
    case {'tfMRI_SOCIAL_LR','tfMRI_SOCIAL_RL','SOCIAL'}
        n = 274;
    case {'tfMRI_WM_LR','tfMRI_WM_RL','WM'}
        n = 405;
    case {'rfMRI_REST1_LR','rfMRI_REST2_LR',...
            'rfMRI_REST1_RL','rfMRI_REST2_RL','REST1','REST2'}
        n = 1200;
end

end