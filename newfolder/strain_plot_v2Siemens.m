% function []=strain_plot_v2Siemens(data,stdev,units,suffix,limits,dt)
% 
% 
% color=[0.9,0,0;0,0.7,0;0,0,0.95];
% t = [0:dt:dt*(size(data,1)-1)];
% t_int = [0:dt/10:dt*(size(data,1)-1)];
% 
%         figure
%         hold on
%         intData = interp1(1:size(data,1),data(:,1),1:0.1:size(data,1),'spline');
%         p1=plot(t_int,intData,'-','Color',color(:,3));
%         p12=plot(t,data(:,1),'.','Color',color(:,3));
% 
%         ylabel(units,'Interpreter','latex','FontSize',16)
%         filename=[suffix,'.eps'];
% 
%         ylim(limits)
% 
% 
% 
%         errorbar(t,data(:,1),stdev(:,1),'.','Color',color(:,3));
% 
%         xlabel('$time \quad [\mathrm{s}]$','Interpreter','latex','FontSize',16)
%         %title(ID,'Interpreter','latex');
% 
%         xlim([0 dt*(size(data,1)+1)])
%         set(gcf,'color','w');
%         set(gca,'TickLabelInterpreter','latex','FontSize',14);
% 
% 
%         export_fig(filename,'-eps');
%         hold off
%         close
% 
% 
% end
function [] = strain_plot_v2Siemens(data, stdev, units, suffix, limits, dt)

    color = [0.9, 0, 0; 0, 0.7, 0; 0, 0, 0.95];
    t = 0:dt:dt*(size(data,1)-1);
    t_int = 0:dt/10:dt*(size(data,1)-1);

    figure
    hold on

    % Interpolated curve
    intData = interp1(1:size(data,1), data(:,1), 1:0.1:size(data,1), 'spline');
    plot(t_int, intData, '-', 'Color', color(:,3));
    plot(t, data(:,1), '.', 'Color', color(:,3));

    % Error bars
    errorbar(t, data(:,1), stdev(:,1), '.', 'Color', color(:,3));

    % Labels and limits
    ylabel(units, 'Interpreter', 'latex', 'FontSize', 16)
    xlabel('$time \quad [\mathrm{s}]$', 'Interpreter', 'latex', 'FontSize', 16)
    ylim(limits)
    xlim([0 dt*(size(data,1)+1)])
    set(gcf, 'Color', 'w');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);

    % Save in all three formats
    filename_eps = [suffix, '.eps'];
    filename_pdf = [suffix, '.pdf'];
    filename_jpg = [suffix, '.jpg'];

    export_fig(filename_eps, '-eps');
    export_fig(filename_pdf, '-pdf');
    export_fig(filename_jpg, '-jpg', '-r300');

    hold off
    close

end
