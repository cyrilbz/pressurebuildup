clear all 
close all

%% Program to compare data between cases
%% post_process must have been done first on each case !

root_dr = pwd ; 

dir_names = {'\no_sugar_fluxes', ...
    '\no_ray_sugar_fluxes', ...
    '\with_all_sugar_fluxes'} ; 

 my_legend = {'No sugar fluxes','With VAC-vessel sugar fluxes','With all sugar fluxes'} ;

%Set width and height of output figure
width_cm = 15;
height_cm = 12;
% Set the screen DPI (replace this with your actual screen DPI)
screen_dpi = 96;
% Convert centimeters to inches
width_inches = width_cm / 2.54;
height_inches = height_cm / 2.54;
% Convert inches to pixels using the screen DPI
width_pixels = width_inches * screen_dpi;
height_pixels = height_inches * screen_dpi;
position = [100, 100, width_pixels,height_pixels]; % Set width and height as desired

%%%% loop over the different directories to load all the data you need
%%%% using 'name.mat' 
%%%% native variables available : t,H,rf,rv,uroot,ufv,upv,vp,pp,nsv,nslc,
%%%% ngv,Vbark_cell,pbark,nsb 
%%%% processed variables :
%%%% th,pwv,Cs_v,Cs_p,Cs_bark,Temp_save,p_mean,D,Dx,Tm_bark,Qbx,ru

for i=1:length(dir_names)
mydr = strcat(root_dr,dir_names{i}) ;
cd(mydr)
th(i) = load('th.mat') ;
ru(i) = load('ru.mat') ;
p_mean(i) = load('p_mean.mat') ;
Cs_bark(i) = load('Cs_bark.mat') ;
Cs_p(i) = load('Cs_p.mat') ;
Cs_v(i) = load('Cs_v.mat') ;
pwv(i) = load('pwv.mat') ;
pp(i) = load('pp.mat') ;
pbark(i) = load('pbark.mat') ;
D(i) = load('D.mat') ;
Dx(i) = load('Dx.mat') ;
Qbx(i) = load('Qbx.mat') ;
Tm_bark(i) = load('Tm_bark.mat') ;
Vbark_cell(i) = load('Vbark_cell.mat') ;
end

cd(root_dr)

%%% define some values you want to export
dd=zeros(length(dir_names),1) ;
pmean=zeros(length(dir_names),1) ;
p_end=zeros(length(dir_names),1) ;
Tmb=zeros(length(dir_names),1) ;

 c=hsv(length(dir_names)) ; % color map
 c= [c ; c ] ; % color map
 %my_styles= ["-","-","-","-","-","-","-"] ;
  my_styles= ["-",":","-","-","-","-","-"] ;
  
 %my_styles=[my_styles; my_styles] ;
 
%%% plot for loop and additional processing
for i=1:length(dir_names)
    %%% instant at which you want to export some values
    ind = find(th(i).th==24) ;

    %% diameter comparison
    set (gcf,'color','white')
    buff = D(i).D ;
    buff = (buff - buff(1))*10^6 ; % diameter variations in µm
    
    figure(1) 
    hold on
    f1 = plot(th(i).th/24,buff,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    xlabel('time [days]','interpreter','latex','FontSize',18)
    ylabel('Stem diameter variations [$\mu$ m]','interpreter','latex','FontSize',18)
    %xlim([0 14])
    box on
    hold off
    dd(i) = buff(ind) ; % store diameter shrinkage
    
    %% mean pressure comparison
    buff = p_mean(i).p_save ;
    figure(2)
    set (gcf,'color','white')
    hold on
    ylim([0 300])
    colormap('cool');
%     myc = linspace(min(buff(1:1*8640:end)), max(buff(1:1*8640:end)), length(buff(1:1*8640:end))); % Color values based on x-coordinates
%     scatter(th(i).th(1:1*8640:end)/24, buff(1:1*8640:end), [], myc, 'filled'); % Use scatter to assign gradient colors
    %f2 = plot(th(i).th(1:8640:end)/24,buff(1:8640:end),'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    f2 = plot(th(i).th,buff,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    xlabel('time [h]','interpreter','latex','FontSize',18)
    ylabel('$\bar p_w^v$ [kPa]','interpreter','latex','FontSize',18)
    box on
    hold off
    pmean(i)=buff(ind) ;
    p_end(i)=buff(end) ;
    
    %% radial pressure gradient comparison
    buff = pwv(i).pwv ;
    set (gcf,'color','white')
    figure(3)
    hold on
    semilogy(ru(i).ru*1000,buff(end,:)/1000,'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    set(gca, 'YScale', 'log')
    set (gcf,'color','white')
    xlabel('r [mm]','interpreter','latex','FontSize',18)
    ylabel('$p_w^v$ (t=72h) [kPa]','interpreter','latex','FontSize',18)
    box on
    hold off
    
    %% Vessel mean sugar content
    buff = Cs_v(i).Cs_v ;
    Cs_vmean(i,:)=sum(buff,2)/length(ru(i).ru) ;
    Cs_v_end(i) = Cs_vmean(i,end) ;
    set (gcf,'color','white')
    figure(4)
    hold on
    semilogy(th(i).th,Cs_vmean(i,:),'Color',c(i,:),'Linewidth',1.5,'Linestyle',my_styles(i)) ;
    %set(gca, 'YScale', 'log')
    set (gcf,'color','white')
    xlabel('t [h]','interpreter','latex','FontSize',18)
    ylabel('$\bar C_s^{v}$  [mol/m3]','interpreter','latex','FontSize',18)
    box on
    hold off
end
% save('Tmb_24h.mat','Tmb')
% save('pmean_24h.mat','pmean')
% save('dd_24h.mat','dd')

% save('Ds_Dsray7-new','char_list')
% save('p_end_Dsray7-new.mat','p_end')
% save('Cs_v_end_Dsray7-new.mat','Cs_v_end')

% save('Cs_list','char_list')
% save('p_end-Ds7.mat','p_end')
% save('Cs_v_end_Dsray7-new.mat','Cs_v_end')


%% radial pressure gradient comparison
%     buff = pwv(1).pwv ;
%     set (gcf,'color','white')
%     figure(11)
%     my_color=cool(length(th(1).th));
%    for i=1:2*8640:length(th(1).th)
%         plot(ru(1).ru*1000,buff(i,:)/1000,'Color',my_color(i,:),'Linewidth',1.5) ; 
%         hold on
%    end
%    xlim([3.25 6])
% % Add a colorbar
% colorbarPosition = [0.23, 0.75, 0.23, 0.03]; % [x, y, width, height]
% % Add a colorbar inside the axes
% colormap(my_color)
% %cbar = colorbar;
% cbar = colorbar('Location', 'southoutside');
% cbar.Position = colorbarPosition;
% desiredTicks = linspace(min(th(1).th), max(th(1).th)/24, 2); % Example tic
% cbarTicks = linspace(0, 1, length(desiredTicks)); % Normalize the ticks to the colormap range
% cbar.TickLabels = arrayfun(@(x) sprintf('%.0f', x), desiredTicks, 'UniformOutput', false);
% cbar.Ticks = cbarTicks;
%     ax = gca;
%     ax.FontSize = 18;       
%     ax.TickLabelInterpreter = 'latex' ;
% set(cbar,'TickLabelInterpreter','latex')
% set(cbar,'FontSize',14)
% text(0.23,0.85,'t [days]','Units','normalized','FontSize',14,'interpreter','latex')
% 
%     set(gca, 'YScale', 'log')
%     set (gcf,'color','white')
%     xlabel('r [mm]','interpreter','latex','FontSize',18)
%     ylabel('$p_w^v$ [kPa]','interpreter','latex','FontSize',18)
%     box on
%     hold off


%% add legends
for i=2
    figure(i)
    lgnd(i) = legend(my_legend,'interpreter','latex','location','north');
    ax = gca;
    ax.FontSize = 18;       
    ax.TickLabelInterpreter = 'latex' ;
    %lgnd(i).Title.String = {'$P_s^{ray}$ [m/s]'} ;
    legend boxoff
    %set(ax,'Color','parula')
end
for i=1:4
    figure(i)
    
    %lgnd(i) = legend(my_legend,'interpreter','latex','location','north');
    ax = gca;
    ax.FontSize = 18;       
    ax.TickLabelInterpreter = 'latex' ;
%     lgnd(i).Title.String = {'$P_s^{ray}$', '[m/s]'} ;
%     legend boxoff
%     set(ax,'Color','parula')
end

% %% additional plots
% Text = max(5 -th(1).th,-10) ; 
% Text2 = max(5 -8*th(1).th,-10) ; 
% figure(5)
% hold on
% plot(th(1).th,Text,'k:','DisplayName','$T_{ext} (-1K/h)$')
% plot(th(1).th,Text2,'k-.','DisplayName','$T_{ext} (-8K/h)$')
% % set(lgnd(5),'color','none')
% text(0.15,0.95,'c)','Units','normalized','FontSize',14,'interpreter','latex')
% hold off
% %print -depsc fig6-Tm_dTCs.eps
% 
% figure(3)
% %set(lgnd(3),'Position',[0.2, 0.4, 0.3, 0.3])
% text(0.025,0.95,'c)','Units','normalized','FontSize',14,'interpreter','latex')
% fig=figure(3);
% set(fig, 'Position', position);
% set(gcf,'PaperPositionMode','auto')
% 
% figure(1)
% % %hold on
%  text(0.1,0.95,'c)','Units','normalized','FontSize',14,'interpreter','latex')
%  ax= gca;
% % lgnd(1).NumColumns = 3 ;
%  
% fig=figure(1);
% set(fig, 'Position', position);
% set(gcf,'PaperPositionMode','auto')
% %print -depsc fig1-diam.eps
%  
%   figure(2)
% % %hold on
%  text(0.025,0.95,'a)','Units','normalized','FontSize',14,'interpreter','latex')
%  ax= gca;
%  fig=figure(2);
% set(fig, 'Position', position);
% set(gcf,'PaperPositionMode','auto')
% 
% figure(11)
% text(0.025,0.95,'b)','Units','normalized','FontSize',14,'interpreter','latex')
% fig=figure(11);
% set(fig, 'Position', position);
% set(gcf,'PaperPositionMode','auto')
