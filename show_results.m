fig = figure();
fig.Position = [300,700,1300,1200];
%%% CPG output %%%
pos = [0.08,0.73,0.88,0.25];
str_pos = [0.507,0.7,0.01,0.01];
ax1 = subplot('Position',pos);
annotation('textbox',str_pos,'String',"(a)",'EdgeColor','none','FontSize',20)
plots = plot_cpg_output(t_sim,x,init_time,trans_time,T_conv1,T_conv2);
ylim(cla.output_lim)
labels = ["$\theta_1$ (LF)","$\theta_2$ (LM)","$\theta_3$ (LH)","$\theta_4$ (RF)","$\theta_5$ (RM)","$\theta_6$ (RH)"];
hold on
plot(t_sim,sigma,'k--') % Threshold
lgd = legend(plots,labels);
lgd.Orientation = "horizontal";
lgd.NumColumns = 3;
lgd.FontSize = 16;
box on

%%% Gaits %%%
pos = [0.08,0.37,0.88,0.3];
str_pos = [0.507,0.34,0.01,0.01];
if min(gaits == [1,2,3])
    gait_pos = [0.07,0.025,0.08,0.025,0.09];
elseif min(gaits == [3,2,1])
    gait_pos = [0.06,0.025,0.08,0.025,0.09];
end

ax2 = subplot('Position',pos);
annotation('textbox',str_pos,'String',"(b)",'EdgeColor','none','FontSize',20)

color_wave = [0.5,0,1];
color_tetrapod = 0.9*[0,1,1];
color_tripod = [0,1,0.5];
colors = [color_wave;color_tetrapod;color_tripod];
color_trans = 0.7*[1,1,1];

% Init
index = 1:length(t_init_sim);
plot_gaits(t_sim,x_out,sigma,index,colors(gaits(1),:));
hold on

t_end = init_time;
y_pos = 1.2*min(X_0(1,:));
[x_begin,y_begin] = utils.axesPos2figPos([0,y_pos],ax1);
[x_end,y_end] = utils.axesPos2figPos([t_end,y_pos],ax1);
annotation("doublearrow",[x_begin,x_end],[y_begin,y_end])
annotation("textbox",[x_begin+gait_pos(1),y_begin-0.01,0.01,0.01],'String',gait_name(gaits(1)),"EdgeColor",'none')

% Trans1
index_trans1 = round(T_conv1/dt);
index = length(t_init_sim) + (1:index_trans1);
plot_gaits(t_sim,x_out,sigma,index,color_trans);
hold on

t_begin = init_time;
t_end = init_time + T_conv1;
y_pos = 1.2*max(X_0(1,:));
[x_begin,y_begin] = utils.axesPos2figPos([t_begin,y_pos],ax1);
[x_end,y_end] = utils.axesPos2figPos([t_end,y_pos],ax1);
annotation("doublearrow",[x_begin,x_end],[y_begin,y_end])
annotation("textbox",[x_begin+gait_pos(2),y_begin+0.015,0.01,0.01],'String',"transition1","EdgeColor",'none')

% Gait2
index = length(t_init_sim) + (index_trans1+1:length(t_trans_sim1));
plot_gaits(t_sim,x_out,sigma,index,colors(gaits(2),:));
hold on

t_begin = init_time + T_conv1;
t_end = init_time + trans_time;
y_pos = 1.2*min(X_0(1,:));
[x_begin,y_begin] = utils.axesPos2figPos([t_begin,y_pos],ax1);
[x_end,y_end] = utils.axesPos2figPos([t_end,y_pos],ax1);
annotation("doublearrow",[x_begin,x_end],[y_begin,y_end])
annotation("textbox",[x_begin+gait_pos(3),y_begin-0.01,0.01,0.01],'String',gait_name(gaits(2)),"EdgeColor",'none')

% Trans2
index_trans2 = round(T_conv2/dt);
index = length([t_init_sim,t_trans_sim1]) + (1:index_trans2);
plot_gaits(t_sim,x_out,sigma,index,color_trans);
hold on

t_begin = init_time + trans_time;
t_end = init_time + trans_time + T_conv2;
y_pos = 1.2*max(X_0(1,:));
[x_begin,y_begin] = utils.axesPos2figPos([t_begin,y_pos],ax1);
[x_end,y_end] = utils.axesPos2figPos([t_end,y_pos],ax1);
annotation("doublearrow",[x_begin,x_end],[y_begin,y_end])
annotation("textbox",[x_begin+gait_pos(4),y_begin+0.015,0.01,0.01],'String',"transition2","EdgeColor",'none')

% Gait3
index = length([t_init_sim,t_trans_sim1]) + (index_trans2+1:length(t_trans_sim2)-2);
plot_gaits(t_sim,x_out,sigma,index,colors(gaits(3),:));
hold on

t_begin = init_time + trans_time + T_conv2;
t_end = init_time + 2*trans_time;
y_pos = 1.2*min(X_0(1,:));
[x_begin,y_begin] = utils.axesPos2figPos([t_begin,y_pos],ax1);
[x_end,y_end] = utils.axesPos2figPos([t_end,y_pos],ax1);
annotation("doublearrow",[x_begin,x_end],[y_begin,y_end])
annotation("textbox",[x_begin+gait_pos(5),y_begin-0.01,0.01,0.01],'String',gait_name(gaits(3)),"EdgeColor",'none')

% Contact
filename = "data/" + cla.name + "_" + savename + "_height_foot_tip.mat";
if exist(filename,"file")
    load(filename); 
    contact = height_foot_tip([4,5,6,1,2,3],:) > 0.005;
    plot_contact(contact,dt)
end

xlim([t_sim(1),t_sim(end)])
xlines(init_time,trans_time,T_conv1,T_conv2)
xtick_settings(init_time,trans_time,T_conv1,T_conv2)

yline([0,1,2,3,4,5],'LineWidth',1)
yline(0.8+[0,1,2,3,4,5],'LineWidth',1)

ylim([0,5.8])
yticks(0.4+[0,1,2,3,4,5])
yticklabels(["RH","RM","RF","LH","LM","LF"])
ax2.TickLength = [0,0];
box on

%%% alpha %%%
pos = [0.08,0.23,0.88,0.08];
str_pos = [0.507,0.2,0.01,0.01];
ax3 = subplot('Position',pos);
annotation('textbox',str_pos,'String',"(c)",'EdgeColor','none','FontSize',20)
plot(t_sim,mod(alpha_direct,2*pi),'r')
hold on
plot(t_sim,mod(alpha_sim,2*pi),'k--')
yticks([1,2,3,4]*pi/3)
yticklabels(["$\frac{1}{3}\pi$","$\frac{2}{3}\pi$","$\pi$","$\frac{4}{3}\pi$"])
ylim([pi/3,4*pi/3])
xlines(init_time,trans_time,T_conv1,T_conv2)
xtick_settings(init_time,trans_time,T_conv1,T_conv2)
ylabel("$\alpha$")
ax3.FontSize = 20;
grid on

%%% beta %%%
pos = [0.08,0.07,0.88,0.1];
str_pos = [0.507,0.02,0.01,0.01];
ax4 = subplot('Position',pos);
annotation('textbox',str_pos,'String',"(d)",'EdgeColor','none','FontSize',20)
color = [0,95,255]/255;
plot(t_sim,mod(beta_direct,2*pi),'Color',color)
hold on
plot(t_sim,mod(beta_sim,2*pi),'k--')
yticks([0,1,2,3,4]*pi/3)
yticklabels(["$0$","$\frac{1}{3}\pi$","$\frac{2}{3}\pi$","$\pi$","$\frac{4}{3}\pi$"])
ylim([0,4*pi/3])
xlines(init_time,trans_time,T_conv1,T_conv2)
xtick_settings(init_time,trans_time,T_conv1,T_conv2)
ylabel("$\beta$")
xlabel("Time (s)")
ax4.FontSize = 20;
grid on

clearvars ax1 ax2 ax3 ax4 lgd plots 
clearvars color color_wave color_tetrapod color_tripod color_trans
clearvars gait_pos pos str_pos x_begin x_end y_begin y_end y_pos

%% plot functions
function p = plot_cpg_output(t_sim,x,init_time,trans_time,T_conv1,T_conv2)
    xlines(init_time,trans_time,T_conv1,T_conv2)
    hold on
    for i = 1:6
        p(i) = plot(t_sim,x(2*i-1,:));
        hold on
    end
    xtick_settings(init_time,trans_time,T_conv1,T_conv2)

    ylabel("CPG output")
end

function c = color_order(color)
    mat = [color;1,1,1];
    c = repmat(mat,[6,1]);
end

function plot_gaits(t,x_out,sigma,index,color)
    swing_stance = x_out(:,index) > sigma(index);
    mat = [0.8*swing_stance(6,:);ones(size(t(index)))-0.8*swing_stance(6,:);...
        0.8*swing_stance(5,:);ones(size(t(index)))-0.8*swing_stance(5,:);...
        0.8*swing_stance(4,:);ones(size(t(index)))-0.8*swing_stance(4,:);...
        0.8*swing_stance(3,:);ones(size(t(index)))-0.8*swing_stance(3,:);...
        0.8*swing_stance(2,:);ones(size(t(index)))-0.8*swing_stance(2,:);...
        0.8*swing_stance(1,:);0.8*ones(size(t(index)))-0.8*swing_stance(1,:)...
        ];
    
    colors = color_order(color);
    
    ax_init = area(t(index),mat');
    for i = 1:12
        ax_init(i).FaceColor = colors(i,:);
    end
end

function plot_contact(contact,dt)
    orange = [255,192,0]/255;
    l = size(contact,2);
    
    for k = 1:6 
        i_TD = -1;
        i_LO = -1;
        if contact(k,1) == 1 % Swing
            i_TD = find(contact(k,:) == 0,1,"first"); % Index of t_TD
            l_swing = i_TD;
            pgon = polyshape([i_LO,i_LO,i_LO+l_swing,i_LO+l_swing]*dt,6-k+[0 +0.8 +0.8 0]);
            plot(pgon,'EdgeColor',orange,'FaceColor',"none",'LineWidth',3,'LineStyle','-')
            hold on
        elseif contact(k,1) == 0 % Stance
            i_LO = find(contact(k,:) == 1,1,"first");
        end

        for i = 1:l
            % Swing
            if i == i_LO % Liftoff
                i_TD = find(contact(k,i_LO:end) == 0,1,"first") + i_LO - 1; % Index of t_TD
                if i_TD
                    l_swing = i_TD - i_LO + 1;
                    pgon = polyshape([i_LO,i_LO,i_LO+l_swing,i_LO+l_swing]*dt,6-k+[0 +0.8 +0.8 0]);
                    plot(pgon,'EdgeColor',orange,'FaceColor',"none",'LineWidth',3,'LineStyle','-')
                    hold on
                else
                    pgon = polyshape([i_LO,i_LO,l-1,l-1]*dt,6-k+[0 +0.8 +0.8 0]);
                    plot(pgon,'EdgeColor',orange,'FaceColor',"none",'LineWidth',3,'LineStyle','-')
                    hold on
                end
            elseif i == i_TD % Touchdown
                i_LO = find(contact(k,i_TD:end) == 1,1,"first") + i_TD - 1; % Index of t_LO
            end
        end
    end
end

function xlines(init_time,trans_time,T_conv1,T_conv2)
    xline(init_time,'LineWidth',1)
    xline(init_time+T_conv1,'LineWidth',1)
    xline(init_time+trans_time,'LineWidth',1)
    xline(init_time+trans_time+T_conv2,'LineWidth',1)
end

function xtick_settings(init_time,trans_time,T_conv1,T_conv2)
    xticks([0,init_time,init_time+T_conv1,init_time+trans_time,...
        init_time+trans_time+T_conv2,init_time+2*trans_time])
    xtickformat('%.1f')
end