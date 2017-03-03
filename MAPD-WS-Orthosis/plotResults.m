% Supplementary code for the study 
% "Optimal control based stiffness identification of an ankle-foot orthosis
% using a predictive walking model", by M. Sreenivasa, M. Millard, M. Felis
% , K. Mombaur & S.I. Wolf
% Contact: M. Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>,
% Heidelberg University, Germany
%
% This program reads the results from the MAPD-WS-Orthosis OCP and plots them. 
% Use the boolean variables below to control which results you wish to plot

clear; 
clf; 
clc;

bPlot_q = 1;
bPlot_qvel = 0;
bPlot_muscleTorques = 0;
bPlot_muscleActivationsExcitations = 0;
bPlot_footContactInfo = 0;
bPlot_orthosisTorques = 0;
bPlot_stages = 1;
bCalcAndPrintFittingErrors = 1;

res_path = 'results';
lnWdt = 2;
aug_dat  = dlmread([res_path,'/pathWalker2d_augmented.txt'],',');
stage_change_idx = find(diff(aug_dat(:,end)) > 0);
stage_change_time = [0;aug_dat(stage_change_idx,1);aug_dat(end,1)];

if bPlot_q
    
    end_cutoff = 1;
    
    csv_dat = dlmread([res_path,'/pathWalker2d.csv'],',',13,0);
    csv_orig = dlmread([res_path,'/expTraj.csv'],',');
    q = csv_dat(1:end-end_cutoff,2:end-1);
    q_orig = csv_orig(1:end-end_cutoff,2:end-1);
    timeStamp = csv_dat(1:end-end_cutoff,1);
    timeStamp_orig = csv_orig(1:end-end_cutoff,1);
        
    if bCalcAndPrintFittingErrors
        % Compute q vector differences
        for i = 1:length(timeStamp_orig)
            curr_time = timeStamp_orig(i);
            diff_time = abs(timeStamp-repmat(curr_time, length(timeStamp), 1));
            [tmp(i),best_idx] = min(diff_time);
            err(i,:) = abs(q_orig(i,:) - q(best_idx,:));
        end
        Abs_Max_Errors = [max(err(:,1:2)) max(err(:,3:end))*180/pi];
        Abs_Mean_Errors = [mean(err(:,1:2)) mean(err(:,3:end))*180/pi];
        RMS_Errors = [rms(err(:,1:2)) rms(err(:,3:end))*180/pi];
        pelvis_rms = rms([err(:,3)])*180/pi;
        hip_rms = rms([err(:,4);err(:,7)])*180/pi;
        knee_rms = rms([err(:,5);err(:,8)])*180/pi;
        ankle_rms = rms([err(:,6);err(:,9)])*180/pi;
        torso_rms = rms([err(:,10)])*180/pi;
        disp(sprintf ('Fitting RMS (pelvis, hips, knees, ankles, torso) =\n %.2f, %.2f, %.2f, %.2f, %.2f',...
            pelvis_rms, hip_rms, knee_rms, ankle_rms, torso_rms));
    end
    
    dim_q = 10;
    labels_q = {'Pos X','Pos Z','q Pelvis','q R Hip','q R Knee','q R Ankle','q L Hip','q L Knee','q L Ankle','q Torso'};
    for i = 1:dim_q
        subplot(4,3,i); hold on;
        plot(timeStamp, q(:,i), '-k');
        plot(timeStamp_orig, q_orig(:,i), '--b');
        if bPlot_stages
            plot([timeStamp(stage_change_idx) timeStamp(stage_change_idx)], [-2 2], '--k');
        end
        xlim([0 timeStamp(end)]);
        ylim([min(q_orig(:,i))-0.3 max(q_orig(:,i))+0.3]);
        ylabel(labels_q(i));
    end
end

if bPlot_muscleActivationsExcitations
    muscleActivations = aug_dat(:,12:25);
    muscleExcitations = aug_dat(:,26:39);
    timeStamp = aug_dat(:,1);
    dim_muscleActivations = 14;
    labels_muscleActivations = {'Ext R Hip','Flex R Hip','Ext R Knee','Flex R Knee','Ext R Ankle','Flex R Ankle','Ext L Hip','Flex L Hip','Ext L Knee','Flex L Knee','Ext L Ankle','Flex L Ankle','Ext Torso','Flex Torso'};
    for i = 1:dim_muscleActivations
        subplot(4,4,i); hold on;
        plot(timeStamp, muscleActivations(:,i), '-r');
        plot(timeStamp, muscleExcitations(:,i), '-b');
        plot([0 timeStamp(end)], [0 0], '--k');
        plot([0 timeStamp(end)], [1 1], '--k');
        title(labels_muscleActivations(i));
        if bPlot_stages
            plot([timeStamp(stage_change_idx) timeStamp(stage_change_idx)], [0 1], '--k');
        end
        xlim([0 timeStamp(end)])
        ylim([-0.1 1.1]);
    end
    
    timeStamp(end)
    sum(sum(muscleActivations))/timeStamp(end)
    
end

if bPlot_orthosisTorques   
    end_cutoff = 5;
    csv_dat = dlmread([res_path,'/pathWalker2d.csv'],',',13,0);
    csv_orig = dlmread([res_path,'/expTraj.csv'],',');
    q = csv_dat(1:end-end_cutoff,2:end-1);
    q_orig = csv_orig(1:end-end_cutoff,2:end-1);
    timeStamp = csv_dat(1:end-end_cutoff,1);
    timeStamp_orig = csv_orig(1:end-end_cutoff,1);
    
    timeStamp_aug = aug_dat(:,1);
    muscleTorques_total = aug_dat(:,40:53);
    muscleTorques_passive = aug_dat(:,54:67);
    muscleTorques_active = muscleTorques_total - muscleTorques_passive;
    orthosis_torques =  aug_dat(:,92:93);    
   
    subplot(3,1,1); hold on;
    plot(timeStamp, q(:,6), '-b');
    plot(timeStamp, q(:,9), '-r');
    plot(timeStamp_orig, q_orig(:,6), '--b');
    plot(timeStamp_orig, q_orig(:,9), '--r');
    plot([timeStamp(1) timeStamp(end)], [0 0], '--k')
    if bPlot_stages
        plot([timeStamp(stage_change_idx) timeStamp(stage_change_idx)], [-20 20], '--k');
    end
    ylim([-0.5 0.25]);
    ylabel('< Flex - Ankle Angle - Ext >');
    legend('R','L');
    
    subplot(3,1,2); hold on;
    plot(timeStamp_aug, muscleTorques_total(:,5), '-b');
    plot(timeStamp_aug, muscleTorques_total(:,11), '-r');
    plot([timeStamp_aug(1) timeStamp_aug(end)], [0 0], '--k')
    if bPlot_stages
        plot([timeStamp(stage_change_idx) timeStamp(stage_change_idx)], [-500 500], '--k');
    end
    ylim([-2 50]);
    ylabel('< Flex - Ankle Muscle Torques - Ext >');
    legend('R','L');
   
    subplot(3,1,3); hold on;
    plot(timeStamp_aug, orthosis_torques(:,1), '-b');
    plot(timeStamp_aug, orthosis_torques(:,2), '-r');   
    plot([timeStamp_aug(1) timeStamp_aug(end)], [0 0], '--k');
    if bPlot_stages
        plot([timeStamp(stage_change_idx) timeStamp(stage_change_idx)], [-500 500], '--k');
    end
    ylim([-75 75]);
    ylabel('Orthosis Torques');
    legend('R','L');
end

if bPlot_qvel
    qvel = aug_dat(:,2:11);
    timeStamp = aug_dat(:,1);
    dim_qvel = 10;
    labels_qvel = {'Vel X','Vel Z','qvel Pelvis','qvel R Hip','qvel R Knee','qvel R Ankle','qvel L Hip','qvel L Knee','qvel L Ankle','qvel Torso'};
    for i = 1:dim_qvel
        subplot(4,3,i); hold on;
        plot(timeStamp, qvel(:,i), '-k');
        if bPlot_stages
            plot([timeStamp(stage_change_idx) timeStamp(stage_change_idx)], [-15 15], '--k');
        end
        ylim([min(qvel(:,i))-0.2 max(qvel(:,i))+0.2])
        xlim([0 timeStamp(end)])
        ylabel(labels_qvel(i));
    end
end

if bPlot_muscleTorques
    muscleTorques_total = aug_dat(:,40:53);
    muscleTorques_passive = aug_dat(:,54:67);
    muscleTorques_active = muscleTorques_total - muscleTorques_passive;
    timeStamp = aug_dat(:,1);
    dim_muscleTorques = 14;
    labels_muscleTorques = {'Ext R Hip','Flex R Hip','Ext R Knee','Flex R Knee','Ext R Ankle','Flex R Ankle','Ext L Hip','Flex L Hip','Ext L Knee','Flex L Knee','Ext L Ankle','Flex L Ankle','Ext Torso','Flex Torso'};
    for i = 1:dim_muscleTorques
        subplot(4,4,i); hold on;
        plot(timeStamp, muscleTorques_total(:,i), '-k');
        plot(timeStamp, muscleTorques_active(:,i), '-r');
        if bPlot_stages
            plot([timeStamp(stage_change_idx) timeStamp(stage_change_idx)], [-50 50], '--k');
        end
        ylim([min(muscleTorques_total(:,i))-2 max(muscleTorques_total(:,i))+2])
        ylabel(labels_muscleTorques(i));
        xlim([0 timeStamp(end)])
    end
end

if bPlot_footContactInfo
    
    stage_change_name = { 'L Toe Off', 'R Heel Off', 'L Heel On', 'L Toe On', 'Right Toe Off', 'Left Heel Off', 'R Heel On', 'R Toe On', 'L Toe Off'};
    
    rHeelPos = aug_dat(:,68:70);
    rHalxPos = aug_dat(:,71:73);
    lHeelPos = aug_dat(:,74:76);
    lHalxPos = aug_dat(:,77:79);
    rHeelForce = aug_dat(:,80:82);
    rHalxForce = aug_dat(:,83:85);
    lHeelForce = aug_dat(:,86:88);
    lHalxForce = aug_dat(:,89:91);
    timeStamp = aug_dat(:,1);
    
    patchColor = [0.6 0.8 0.6; 0.6 0.6 0.8; 0.8 0.6 0.6];
    
       subplot(3,1,1); hold on; title ('Foot contact positions');
    plot(timeStamp, lHeelPos(:,3), '-r', 'linewidth', lnWdt);
    plot(timeStamp, lHalxPos(:,3), '--r', 'linewidth', lnWdt);
    plot(timeStamp, rHeelPos(:,3), '-b', 'linewidth', lnWdt);
    plot(timeStamp, rHalxPos(:,3), '--b', 'linewidth', lnWdt);
    legend('lHeel','lHalx','rHeel','rHalx');    
    plot([0 timeStamp(end)], [0 0], '--k');
    xlim([0 timeStamp(end)])
    if bPlot_stages
        plot([timeStamp(stage_change_idx) timeStamp(stage_change_idx)], [-50 50], '--k');
    end
    ylim([-0.05 0.2]);
    
    subplot(3,1,2); hold on;  title ('Anterior-Posterior forces');
    plot(timeStamp, lHeelForce(:,1)+lHalxForce(:,1), '-r', 'linewidth', lnWdt);
    plot(timeStamp, rHeelForce(:,1)+rHalxForce(:,1), '-b', 'linewidth', lnWdt);
    plot([0 timeStamp(end)], [0 0], '--k');
    xlim([0 timeStamp(end)]);
    if bPlot_stages
        plot([stage_change_time stage_change_time], [-100 100], '--k');
    end
    ylim([-100 100]);
    
    subplot(3,1,3); hold on; title ('Vertical forces');
    plot(timeStamp, lHeelForce(:,3), '-r');
    plot(timeStamp, lHalxForce(:,3), '--r');
    plot(timeStamp, rHeelForce(:,3), '-b');
    plot(timeStamp, rHalxForce(:,3), '--b');
    plot(timeStamp, lHeelForce(:,3)+lHalxForce(:,3), '-r', 'linewidth', lnWdt);
    plot(timeStamp, rHeelForce(:,3)+rHalxForce(:,3), '-b', 'linewidth', lnWdt);
    legend('lHeel','lHalx','rHeel','rHalx', 'L GRF', 'R GRF');
    plot([0 timeStamp(end)], [0 0], '--k');
    xlim([0 timeStamp(end)]);
    if bPlot_stages
        plot([stage_change_time stage_change_time], [-100 420], '--k');
    end
    ylim([-75 450]);
end