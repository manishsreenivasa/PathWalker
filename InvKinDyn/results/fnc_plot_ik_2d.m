function [ik_out_res,pelvis_pos_xz_res,t_out] = fnc_plot_ik_2d (file_ik, l_step, r_step, bPlot_trials, modelColor)

ik_in   = dlmread(file_ik);
ik_in(:,4:end)   = ik_in(:,4:end)*180/pi;
l_stepDuration = l_step(2)-l_step(1);
r_stepDuration = r_step(2)-r_step(1);

if (~isnan(r_stepDuration))
    ik_out_r = ik_in(r_step(1):r_step(2),5:7);
    t_r = ik_in(r_step(2),1)-ik_in(r_step(1),1);
    ik_out_r_res = resample(ik_out_r,100,length(ik_out_r),0);
    
    pelvis_pos_xz = ik_in(r_step(1):r_step(2),2:4);
    pelvis_pos_xz_res = resample(pelvis_pos_xz,100,length(pelvis_pos_xz),0);
    pelvis_pos_xz_res(:,1) = pelvis_pos_xz_res(:,1)-pelvis_pos_xz_res(1,1);
    
    if bPlot_trials
        subplot(2,6,2); hold on;
        plot([1:100],-ik_out_r_res(:,1),'k','linewidth',1,'linestyle','--', 'color', modelColor);
        
        subplot(2,6,4); hold on;
        plot([1:100],ik_out_r_res(:,2),'k','linewidth',1,'linestyle','--', 'color', modelColor);
        
        subplot(2,6,6); hold on;
        plot([1:100],-ik_out_r_res(:,3),'k','linewidth',1,'linestyle','--', 'color', modelColor);
    end
else
    ik_out_r_res = NaN(100,3);
    pelvis_pos_xz_res = NaN(100,3);
    t_r = NaN;
end
if (~isnan(l_stepDuration))
    
    ik_out_l = ik_in(l_step(1):l_step(2),8:10);
    t_l = ik_in(l_step(2),1)-ik_in(l_step(1),1);
    ik_out_l_res = resample(ik_out_l,100,length(ik_out_l),0);
    
    pelvis_pos_xz = ik_in(l_step(1):l_step(2),2:4);
    pelvis_pos_xz_res = resample(pelvis_pos_xz,100,length(pelvis_pos_xz),0);
    pelvis_pos_xz_res(:,1) = pelvis_pos_xz_res(:,1)-pelvis_pos_xz_res(1,1);
    
    if bPlot_trials
        subplot(2,6,1); hold on;
        plot([1:100],-ik_out_l_res(:,1),'k','linewidth',1,'linestyle','--', 'color', modelColor);
        
        subplot(2,6,3); hold on;
        plot([1:100],ik_out_l_res(:,2),'k','linewidth',1,'linestyle','--', 'color', modelColor);
        
        subplot(2,6,5); hold on;
        plot([1:100],-ik_out_l_res(:,3),'k','linewidth',1,'linestyle','--', 'color', modelColor);
    end
else
     ik_out_l_res = NaN(100,3);
     pelvis_pos_xz_res = NaN(100,3);
     t_l = NaN;
end
ik_out_res = [ik_out_l_res ik_out_r_res];
t_out = [t_l t_r];