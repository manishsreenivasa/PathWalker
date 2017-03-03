function [id_out_res] = fnc_plot_id_2d (file_id, l_step, r_step, patientWeight, bPlot_trials, modelColor)

id_in = dlmread(file_id);
id_in = id_in/patientWeight;
l_stepDuration = l_step(2)-l_step(1);
r_stepDuration = r_step(2)-r_step(1);

if (~isnan(r_stepDuration))
    
    id_out_r = id_in(r_step(1):r_step(2),5:7);
    id_out_r_res = resample(id_out_r,100,length(id_out_r),0);
    
    if bPlot_trials
        subplot(2,6,8); hold on;
        plot([1:100],id_out_r_res(:,1),'k','linewidth',1,'linestyle','--', 'color', modelColor);
        
        subplot(2,6,10); hold on;
        plot([1:100],id_out_r_res(:,2),'k','linewidth',1,'linestyle','--', 'color', modelColor);
        
        subplot(2,6,12); hold on;
        plot([1:100],id_out_r_res(:,3),'k','linewidth',1,'linestyle','--', 'color', modelColor);
    end
else
    id_out_r_res = NaN(100,3);
end
if (~isnan(l_stepDuration))
    id_out_l = id_in(l_step(1):l_step(2),8:10);
    id_out_l_res = resample(id_out_l,100,length(id_out_l),0);
    
    if bPlot_trials
        subplot(2,6,7); hold on;
        plot([1:100],id_out_l_res(:,1),'k','linewidth',1,'linestyle','--', 'color', modelColor);
        
        subplot(2,6,9); hold on;
        plot([1:100],id_out_l_res(:,2),'k','linewidth',1,'linestyle','--', 'color', modelColor);
        
        subplot(2,6,11); hold on;
        plot([1:100],id_out_l_res(:,3),'k','linewidth',1,'linestyle','--', 'color', modelColor);
    end
else
    id_out_l_res = NaN(100,3);
end
id_out_res = [id_out_l_res id_out_r_res];