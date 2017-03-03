% Read and plot IK and ID results %
clear; clc; clf;

bPlot_trials = 0;
bPlot_avg = 1;

modelType = {'barefeet','orthosis'};
modelColor = {'r','g'};
modelLine = {'-','-','-'};
patchColor = [0.8 0.6 0.6; 0.6 0.8 0.6];

rootDir = '../';
patientWeight = 24.7;
[b,a] = butter(3,0.2,'low');
fntSz = 12; lnWidth = 2;

% Barefeet C3Ds
usefiles(1).ik = [3053171 3053173:3053184];
usefiles(1).id = [3053176 3053178 3053179 3053183 3053184];
% Orthosis C3Ds
usefiles(2).ik = [3053138:3053148 3053150];
usefiles(2).id = [3053141 3053142 3053144 3053146];

for model = 1:length(modelType)
    for trial = 1:length(usefiles(model).ik)
        c3d_fileName = [rootDir,'Data/c3d_',char(modelType(model)),'/',int2str(usefiles(model).ik(trial)),'.c3d'];
        [l_step,r_step] = fnc_getSteppingEvents (c3d_fileName, false);
        ik_fileName = [rootDir,'InvKinDyn/results/',char(modelType(model)),'_',int2str(usefiles(model).ik(trial)),'.csv'];
        [reps_ik(:,:,trial),reps_pos(:,:,trial),reps_t(:,:,trial)] = fnc_plot_ik_2d (ik_fileName, l_step, r_step, bPlot_trials, char(modelColor(model)));
    end
    for trial = 1:length(usefiles(model).id)
        c3d_fileName = [rootDir,'Data/c3d_',char(modelType(model)),'/',int2str(usefiles(model).id(trial)),'.c3d'];
        [l_step,r_step] = fnc_getSteppingEvents(c3d_fileName, true);
        id_fileName = [rootDir,'InvKinDyn/results/',char(modelType(model)),'_',int2str(usefiles(model).id(trial)),'.torque'];
        reps_id(:,:,trial) = fnc_plot_id_2d (id_fileName, l_step, r_step, patientWeight, bPlot_trials, char(modelColor(model)));
        display([int2str(usefiles(model).id(trial)),' ',char(modelType(model))]);
    end
    
    modelAvg_ik(model,:,:) = filtfilt(b,a,squeeze(nanmean(reps_ik,3)));
    modelStd_ik(model,:,:) = filtfilt(b,a,squeeze(nanstd(reps_ik,0,3)));
    modelAvg_id(model,:,:) = filtfilt(b,a,squeeze(nanmean(reps_id,3)));
    modelStd_id(model,:,:) = filtfilt(b,a,squeeze(nanstd(reps_id,0,3)));
    modelAvg_pos(model,:,:) = filtfilt(b,a,squeeze(nanmean(reps_pos,3)));
    modelStd_pos(model,:,:) = filtfilt(b,a,squeeze(nanstd(reps_pos,0,3)));
    modelAvg_t(model,:,:) = squeeze(nanmean(reps_t,3));
    
    timeStamps = [1:100];
    
    if bPlot_avg
        subplot(2,6,1); hold on;
        patch([timeStamps fliplr(timeStamps)],[-modelAvg_ik(model,:,1)+modelStd_ik(model,:,1) fliplr(-modelAvg_ik(model,:,1)-modelStd_ik(model,:,1))], patchColor(model,:),'linestyle','none');
        ylim([-20 60]);
        title('Left');
        ylabel('<-- Ext.    (\circ)    Flex.   -->');
        
        subplot(2,6,2); hold on;
        patch([timeStamps fliplr(timeStamps)],[-modelAvg_ik(model,:,4)+modelStd_ik(model,:,4) fliplr(-modelAvg_ik(model,:,4)-modelStd_ik(model,:,4))], patchColor(model,:),'linestyle','none');
        ylim([-20 60]);
        title('Right');
        
        subplot(2,6,3); hold on;
        patch([timeStamps fliplr(timeStamps)],[modelAvg_ik(model,:,2)+modelStd_ik(model,:,2) fliplr(modelAvg_ik(model,:,2)-modelStd_ik(model,:,2))], patchColor(model,:),'linestyle','none');
        ylim([-5 80]);
        title('Left');
        
        subplot(2,6,4); hold on;
        patch([timeStamps fliplr(timeStamps)],[modelAvg_ik(model,:,5)+modelStd_ik(model,:,5) fliplr(modelAvg_ik(model,:,5)-modelStd_ik(model,:,5))], patchColor(model,:),'linestyle','none');
        ylim([-5 80]);
        title('Right');
        
        subplot(2,6,5); hold on;
        patch([timeStamps fliplr(timeStamps)],[-modelAvg_ik(model,:,3)+modelStd_ik(model,:,3) fliplr(-modelAvg_ik(model,:,3)-modelStd_ik(model,:,3))], patchColor(model,:),'linestyle','none');
        ylim([-15 30]);
        title('Left');
        
        subplot(2,6,6); hold on;
        patch([timeStamps fliplr(timeStamps)],[-modelAvg_ik(model,:,6)+modelStd_ik(model,:,6) fliplr(-modelAvg_ik(model,:,6)-modelStd_ik(model,:,6))], patchColor(model,:),'linestyle','none');
        ylim([-15 30]);
        title('Right');
        
        subplot(2,6,7); hold on;
        patch([timeStamps fliplr(timeStamps)],[modelAvg_id(model,:,1)+modelStd_id(model,:,1) fliplr(modelAvg_id(model,:,1)-modelStd_id(model,:,1))], patchColor(model,:),'linestyle','none');
        ylabel('<-- Flex.  (Nm/Kg)   Ext. -->'); xlabel('% Step');
        ylim([-1.25 1.25]);
        title('Left');
        
        subplot(2,6,8); hold on;
        patch([timeStamps fliplr(timeStamps)],[modelAvg_id(model,:,4)+modelStd_id(model,:,4) fliplr(modelAvg_id(model,:,4)-modelStd_id(model,:,4))], patchColor(model,:),'linestyle','none');
        ylim([-1.25 1.25]);
        title('Right');
        
        subplot(2,6,9); hold on;
        patch([timeStamps fliplr(timeStamps)],[modelAvg_id(model,:,2)+modelStd_id(model,:,2) fliplr(modelAvg_id(model,:,2)-modelStd_id(model,:,2))], patchColor(model,:),'linestyle','none');
        ylim([-1.25 1.25]);
        title('Left');
        
        subplot(2,6,10); hold on;
        patch([timeStamps fliplr(timeStamps)],[modelAvg_id(model,:,5)+modelStd_id(model,:,5) fliplr(modelAvg_id(model,:,5)-modelStd_id(model,:,5))], patchColor(model,:),'linestyle','none');
        ylim([-1.25 1.25]);
        title('Right');
        
        subplot(2,6,11); hold on;
        patch([timeStamps fliplr(timeStamps)],[modelAvg_id(model,:,3)+modelStd_id(model,:,3) fliplr(modelAvg_id(model,:,3)-modelStd_id(model,:,3))], patchColor(model,:),'linestyle','none');
        ylim([-1.25 1.25]);
        title('Left');
        
        subplot(2,6,12); hold on;
        patch([timeStamps fliplr(timeStamps)],[modelAvg_id(model,:,6)+modelStd_id(model,:,6) fliplr(modelAvg_id(model,:,6)-modelStd_id(model,:,6))], patchColor(model,:),'linestyle','none');
        ylim([-1.25 1.25]);
        title('Right');
    end
end

% Plot averages on top for better visualization
if bPlot_avg
    for model = 1:length(modelType)
        subplot(2,6,1); hold on;
        plot([1:100],-modelAvg_ik(model,:,1),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,2); hold on;
        plot([1:100],-modelAvg_ik(model,:,4),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,3); hold on;
        plot([1:100],modelAvg_ik(model,:,2),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,4); hold on;
        plot([1:100],modelAvg_ik(model,:,5),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,5); hold on;
        plot([1:100],-modelAvg_ik(model,:,3),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,6); hold on;
        plot([1:100],-modelAvg_ik(model,:,6),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        % ID
        subplot(2,6,7); hold on;
        plot([1:100],modelAvg_id(model,:,1),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,8); hold on;
        plot([1:100],modelAvg_id(model,:,4),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,9); hold on;
        plot([1:100],modelAvg_id(model,:,2),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,10); hold on;
        plot([1:100],modelAvg_id(model,:,5),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,11); hold on;
        plot([1:100],modelAvg_id(model,:,3),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        subplot(2,6,12); hold on;
        plot([1:100],modelAvg_id(model,:,6),'k','linewidth',lnWidth,'linestyle',char(modelLine(model)),'color',char(modelColor(model)));
        
        legend('Barefeet 2D','Neuroswing 2D','Location','South');
    end
    for i=1:12
        subplot(2,6,i); hold on; set(gca,'fontsize',fntSz,'xtick',[0 60 100]);
        plot([0 100],[0 0],'--k');
        plot([60 60],[-100 100],'--k');
        xlim([0 100]);
        axis square;
    end
end