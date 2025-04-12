classdef utils
    methods
    end
    methods (Static)
        function plot_lc(p,cla,p_rep)
            if exist("p_rep","var")
                plot(p_rep(1,:),p_rep(2,:),'r')
                hold on
                plot(p(1,:),p(2,:),'k--')
            else
                plot(p(1,:),p(2,:),'r')
            end
            xlabel("$x$")
            ylabel("$y$")
            xlim(cla.x_lim)
            ylim(cla.y_lim)
            xticks(cla.x_tick)
            yticks(cla.y_tick)
            axis square
        end

        function plot_lc_sol(theta,X_0)
            color_plot = ["r","b","g"];
            component = ["x","y","z"];
            for i = 1:size(X_0,1)
                plot(theta,X_0(i,:),color_plot(i),'DisplayName',"$X_{0," + component(i) + "}(\theta)$")
                hold on
            end
            hold off
            legend('Location','north')
            xlim([0,2*pi])
            range = (max(X_0,[],'all') - min(X_0,[],'all')) / 1.5;
            ylim([min(X_0,[],'all') - range,max(X_0,[],'all') + range])
            xticks([0,pi/2,pi,3*pi/2,2*pi])
            xticklabels(["$0$","$\frac{1}{2}\pi$","$\pi$","$\frac{3}{2}\pi$","$2\pi$"])
            xlabel("$\theta$")
            ylabel("$X_{0}(\theta)$")
        end
        
        function plot_psf(theta,Z)
            color_plot = ["r","b"];
            component = ["x","y"];
            for i = 1:size(Z,1)
                plot(theta,Z(i,:),color_plot(i),'DisplayName',"$Z_{" + component(i) + "}(\theta)$")
                hold on
            end
            hold off
            legend('Location','southwest')
            xlim([0,2*pi])
            range = (max(Z,[],'all') - min(Z,[],'all')) / 3;
            ylim([min(Z,[],'all') - range,max(Z,[],'all') + range])
            xticks([0,pi/2,pi,3*pi/2,2*pi])
            xticklabels(["$0$","$\frac{1}{2}\pi$","$\pi$","$\frac{3}{2}\pi$","$2\pi$"])
            xlabel("$\theta$")
            ylabel("$Z(\theta)$")
        end

        function fig = show_lc_psf(t,X_0,Z,cla)
            fig = figure();
            fig.Position(3:4) = [1000,350];

            pos = [0.05,0.25,0.2,0.7];
            subplot('Position',pos)
            str_pos = [0.13,0,0.1,0.1];
            utils.plot_lc(X_0,cla)
            annotation('textbox',str_pos,'String',"(a)",'EdgeColor','none','FitBoxToText','on','FontName','Times New Roman')

            pos = [0.35,0.25,0.25,0.7];
            subplot('Position',pos)
            str_pos = [0.456,0,0.1,0.1];
            utils.plot_lc_sol(t,X_0)
            annotation('textbox',str_pos,'String',"(b)",'EdgeColor','none','FitBoxToText','on','FontName','Times New Roman')

            pos = [0.7,0.25,0.25,0.7];
            subplot('Position',pos)
            str_pos = [0.806,0,0.1,0.1];
            utils.plot_psf(t,Z)
            annotation('textbox',str_pos,'String',"(c)",'EdgeColor','none','FitBoxToText','on','FontName','Times New Roman')

            fontsize(fig,24,"points")
        end
        
        function plot_pcf_general(phi,Gamma,phi_star,phi_bar,y_label,plot_color)
            if ~exist("plot_color",'var')
                plot_color = 'r';
            end

            if phi_star == 0 || phi_star == pi
                ticks = [-pi,0,pi];
                ticklabels = ["$-\pi$","$0$","$\pi$"];
            elseif phi_star == 2*pi/3
                ticks = [-pi,-phi_star,0,phi_star,pi];
                ticklabels = ["$-\pi$","$-\frac{2}{3}\pi$","$0$","$\frac{2}{3}\pi$","$\pi$"];
            end

            yline(0,'k','LineWidth',1)
            hold on
            xline(phi_star,'k--','LineWidth',1)
            plot(phi,Gamma,plot_color)
            hold on
            scatter(phi_star,0,1000,'k.')
            hold on
            plot(phi_bar,0,'ok')
            
            xticks(ticks)
            xticklabels(ticklabels)
            xlim([ticks(1),ticks(end)])
            ylim([min(Gamma)*1.3,max(Gamma)*1.3])
            xlabel("$\varphi$")
            ylabel(y_label)
            xtickangle(0)
            box on
        end

        function fig = show_pcf2(phi,Gamma1,Gamma2,phi_star,phi_bar)
            fig = figure();
            fig.Position(3:4) = [700,350];
            
            pos = [0.1,0.25,0.35,0.7];
            str_pos = [0.25,0.01,0.1,0.1];
            subplot('Position',pos)
            utils.plot_pcf_general(phi,Gamma1,phi_star(1),phi_bar(1),"$\Gamma_{\mathrm{odd}}(\varphi)$");
            annotation('textbox',str_pos,'String',"(a)",'EdgeColor','none','FitBoxToText','on','FontName','Times New Roman')
            box on

            pos = [0.55,0.25,0.35,0.7];
            str_pos = [0.7,0.01,0.1,0.1];
            subplot('Position',pos)
            utils.plot_pcf_general(phi,Gamma2,phi_star(2),phi_bar(2),"$\Gamma_{\mathrm{even}}(\varphi)$");
            annotation('textbox',str_pos,'String',"(b)",'EdgeColor','none','FitBoxToText','on','FontName','Times New Roman')
            box on
        end

        function [xFig,yFig] = axesPos2figPos(data,handle_axes)
            x = data(1);
            y = data(2);
            handle_axes.Units = 'normalize';
            axesPos = handle_axes.Position;
            %  axesPos(1): x position of axes in figure
            %  axesPos(2): y position of axes in figure
            %  axesPos(3): Width of axes in figure scale
            %  axesPos(4): Height of axes in figure scale
            widthData = handle_axes.XLim(2)-handle_axes.XLim(1);
            heightData = handle_axes.YLim(2)-handle_axes.YLim(1);
            xmin = handle_axes.XLim(1);
            ymin = handle_axes.YLim(1);
        
            xFig = (x-xmin)/widthData*axesPos(3) + axesPos(1);
            yFig = (y-ymin)/heightData*axesPos(4) + axesPos(2);
        end

        function save_fig(save,fig,save_name,cla)
            if save == 1
                if exist("cla",'var')
                    exportgraphics(fig,"figs/" + cla.name + "_" + save_name + ".pdf")
                else
                    exportgraphics(fig,"figs/" + save_name + ".pdf")
                end
            end
        end
    end
end