clc, clear, close all
% Sequence of data extraction:
% 1. Experiment
% 2. DIC (VIC3D)
% 3. Preprocessing: filter data and rotate point cloud to global frame
% 4. Postprocessing I: generate maps of curvature and rotation
% 5. Postprocessing II: extract fold location and fold angle
% 6. Comparison: classify experimental data / generate statistics



%% 4. Postprocessing of data from DIC

% post1_extract_curvature([30,1]);

%% 5. EXTRACT FOLD EVOLUTION AND STATISTICS

strip = 3;
tN = 5;
rN = 1;

tmax = 400;
% plotOn = false;
plotOn = true;

update = false;
% update = true;

loadSettings = true;
% loadSettings = false;


if loadSettings
   if exist(sprintf('%s\\Test %d\\exp_test%d_run%d.mat',pwd,strip,tN,rN))
       disp('Loading saved parameters...')
        load(sprintf('%s\\Test %d\\exp_test%d_run%d.mat',pwd,strip,tN,rN))    
   else
       disp('Warning: settings not found! Applying default values...')
   end
end



%<<<<<<<<<<<<<<<<<<<<<<<<<<< Parameters >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
if ~exist('settings','var')
    
%     thmode = 'max';
    thmode = 'median';
    
    dlam = 0.1; % Trust region for fold tracking
    kmin = 5e-3; % Minimum fold curvature
    minInitialHeight = 0.005;
    dlam_maxdist = 0.1;
    excludeRange = 0.05;
    max_disc = 20;
    
    
    bypass_td_calculation = false;
%     bypass_td_calculation = true;
    i1 = 295;
    i2 = i1;
    
    
    
else
    thmode = settings.thmode;
    dlam = settings.dlam;
    kmin = settings.kmin;
    minInitialHeight = settings.minInitialHeight;
    dlam_maxdist = settings.dlam_maxdist;
    excludeRange = settings.excludeRange;
    max_disc = settings.max_disc;
    bypass_td_calculation = settings.bypass_td_calculation;
    i1 = settings.te_ind(1);
    i2 = settings.te_ind(2);
end


battens = [-1:0.5:1];
close all

    load(sprintf('%s\\Test %d\\post_test%d_run%d.mat',pwd,strip,tN,rN))
    tds = tdic;
    load(sprintf('%s\\Test %d\\test%d_run%d_coords.mat',pwd,strip,tN,rN))

    
    window = 100;
    th1 = smoothdata(th1,2,'lowess',window);
    th2 = smoothdata(th2,2,'lowess',window);
    
%     lambda2 = lambda2 / s2(end);
    %==========================================================================
    % 5.1 Statistics
    %==========================================================================
    Np1 = size(lambda1,2);
    Np2 = size(lambda2,2);
    Nsteps = size(k1,1);
    
    %--------------------------------------------------------------------------
    % 5.1.2 Initial deployment time
    %--------------------------------------------------------------------------
    % <<<<<<<<<<<<<<<<<<<<< Parameters >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    NframeAve = 20; % # of frames used to compute displacement RMS error
    
    % Compute RMS error between nominally identical frames
    disp1 = sqrt(mean((Z1dic - Z1dic(1,:)).^2,2,'omitnan'));
    disp2 = sqrt(mean((Z2dic - Z2dic(1,:)).^2,2,'omitnan'));
    % disp1 = sqrt(mean((th1 - th1(1,:)).^2,2,'omitnan'));
    % disp2 = sqrt(mean((th2 - th2(1,:)).^2,2,'omitnan'));
    
    
    
    % Define thereshold as mean + 3 sigma of RMS error
    threshold1 = mean(disp1(1:NframeAve)) + 3*std(disp1(1:NframeAve));
    threshold2 = mean(disp2(1:NframeAve)) + 3*std(disp2(1:NframeAve));
    
    [threshold1,threshold2 ] = deal(0.5);
    
    % Find initiation of deployment
    ind1 = find(disp1>=threshold1,1,'first');
    ind2 = find(disp2>=threshold2,1,'first');
    
    
    t0 = tdic(min([ind1,ind2]));
    t = tds - t0;
    t = 1e3*t;
    
    %==========================================================================
    % 5.2 Fold angle
    %==========================================================================
    % <<<<<<<<<<<<<<<<<<< Parameter >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Lave = 0.5; % Length to use for angle calculation
    
    
    % 5.2.1 Compute fold angle
    if strncmp(thmode,'median',6)
        th1v_raw = [median(th1(:,lambda1(1,:)<-1 + Lave),2,'omitnan'),...
            median(th1(:,lambda1(1,:)>1-Lave),2,'omitnan')];
        
        th2v_raw = [median(th2(:,lambda2(1,:)<-1 + Lave),2,'omitnan'),...
            median(th2(:,lambda2(1,:)>1-Lave),2,'omitnan')];
    else
        
        
        th1v_raw = [min(th1(:,1:floor(Np1/2)),[],2),...
            max(th1(:,floor(Np1/2)+1:Np1),[],2)];
        
        th2v_raw = [min(th2(:,1:floor(Np2/2)),[],2),...
            max(th2(:,floor(Np2/2)+1:Np2),[],2)];
    end
    
    th1v = smoothdata(th1v_raw,1,'movmedian',10);
    th2v = smoothdata(th2v_raw,'movmedian',10);
    
    
    %--------------------------------------------------------------------------
    % 5.2.2 End of deployment
    %--------------------------------------------------------------------------
    
    deth1 = abs(th1v(:,1)-th1v(:,2));
    deth2 = abs(th2v(:,1)-th2v(:,2));
    [~,ind1] = min(deth1);
    [~,ind2] = min(deth2);
    threshold1 = mean(deth1(ind1:end),'omitnan') + 3*std(deth1(ind1:end),'omitnan');
    threshold2 = mean(deth2(ind2:end),'omitnan') + 3*std(deth2(ind2:end),'omitnan');
    
    figure
    plot(deth1,'LineWidth',2);
    hold on
    plot(deth2,'LineWidth',2);
    grid on
    xlabel('Frame #')
    ylabel('\Delta\theta [deg]')
    
    
    set(gca,'FontSize',14)
    set(gcf,'Units','centimeters')
    set(gcf,'Position',[38,12,12,8])
    
    
    if ~bypass_td_calculation
        i1 = find(deth1<=threshold1,1,'first');
        i2 = find(deth2<=threshold2,1,'first');
    end

%     
% %     
    % CHECKPOINT: verify that complete deployment occurs
    if ~isempty(i1) && ~isempty(i2)
        tf = t(max([i1,i2]));
    else
        tf = t(end);
        i2 = length(t);
        i1 = length(t);
    end
    
    
    %     tf = t(i2);
    
    if plotOn
        % USER INTERFACE: Plot angle map
        fig3 = figure;
        surf(t,lambda1',th1','EdgeColor','none','FaceColor','interp')
        hold on
        plot3([0,0],[-1,1],[60,60],'-k')
        plot3(tf*[1,1],[-1,1],[60,60],'-k')
        xlim([0,tf])
        xlim([0,tmax])
        xlabel('Time [ms]')
        ylabel('\xi')
        h = colorbar('Xtick',[-90:15:90]);
        set(get(h,'Title'),'String', '\theta [deg]')
        colormap('jet')
        view([0,90])
        caxis([-50,50])
        set(gcf,'Units','centimeters')
        set(gcf,'Position',[2,2,15,10])
        set(gca,'FontSize',14)
%         title('Longeron 1')
        ylim([-1,1])
        yticks([-1:0.25:1])
        
        for JJ = 1: length(battens)
            plot3([0,tmax],battens(JJ)*[1,1],[60,60],'-k','LineWidth',1)
        end
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fig4 = figure;
        surf(t,lambda2',th2','EdgeColor','none','FaceColor','interp')
        hold on
        plot3([0,0],[-1,1],[60,60],'-k')
        plot3(tf*[1,1],[-1,1],[60,60],'-k')
        xlim([0,tf])
        xlim([0,tmax])
        xlabel('Time [ms]')
        ylabel('\xi')
        h = colorbar('Xtick',[-90:15:90]);
        set(get(h,'Title'),'String', '\theta [deg]')
        colormap('jet')
        view([0,90])
        caxis([-50,50])
        set(gcf,'Units','centimeters')
        set(gcf,'Position',[2,12,15,10])
        set(gca,'FontSize',14)
%         title('Longeron 2')
        ylim([-1,1])
        yticks([-1:0.25:1])
        for JJ = 1: length(battens)
            plot3([0,tmax],battens(JJ)*[1,1],[60,60],'-k','LineWidth',1)
        end
    end
    
    %==========================================================================
    % 5.3 Fold location
    %==========================================================================
%     minDistance = dlam_maxdist/mean(lambda1(1,2:end)-lambda1(1,1:end-1));
       
    
    fold1 = track_fold(k1,lambda1, dlam, kmin, minInitialHeight, dlam_maxdist, excludeRange,max_disc);
    fold2 = track_fold(k2,lambda2, dlam, kmin, minInitialHeight, dlam_maxdist, excludeRange,max_disc);
    Nfolds1 = size(fold1,2);
    Nfolds2 = size(fold2,2);

    
    % USER INTERFACE: Plot curvature map
    if plotOn
        fig1 = figure;
        surf(t,lambda1', k1','EdgeColor','none','FaceColor','interp')
        hold on
        plot3([0,0],[-1,1],[1,1],'k','LineWidth',2)
        plot3(tf*[1,1],[-1,1],[1,1],'-k')
        
        xlabel('Time [ms]')
        ylabel('\xi')
        h = colorbar;
        set(get(h,'Title'),'String', '\kappa [mm^{-1}]')
        colormap('jet')
        view([0,90])
        caxis([0,0.015])
        
        hold on
        for JJ = 1: length(battens)
            plot3([0,tmax],battens(JJ)*[1,1],[1,1],'-w','LineWidth',1)
        end
%         plot3(t, fold1,ones(Nsteps,1),'-r','LineWidth',2)
        set(gcf,'Units','centimeters')
        set(gcf,'Position',[14,2,15,10])
        set(gca,'FontSize',14)
        xlim([0,tf])
        xlim([0,tmax])
%         title('Longeron 1')
        ylim([-1,1])
        yticks([-1:0.25:1])

        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        fig2 = figure;
        surf(t,lambda2', k2','EdgeColor','none','FaceColor','interp')
        hold on
        plot3([0,0],[-1,1],[1,1],'k','LineWidth',2)
        plot3(tf*[1,1],[-1,1],[1,1],'-k')
        
        xlabel('Time [ms]')
        ylabel('\xi')
        h = colorbar;
        set(get(h,'Title'),'String', '\kappa [mm^{-1}]')
        colormap('jet')
        view([0,90])
        caxis([0,0.02])
        
        hold on
        for JJ = 1: length(battens)
            plot3([0,tmax],battens(JJ)*[1,1],[1,1],'-w','LineWidth',1)
        end
%         plot3(t, fold2,ones(Nsteps,1),'-r','LineWidth',2)
        set(gcf,'Units','centimeters')
        set(gcf,'Position',[14,12,15,10])
        set(gca,'FontSize',14)
        xlim([0,tf])
        xlim([0,tmax])
        ylim([-1,1])
        
%         title('Longeron 2')
        yticks([-1:0.25:1])

        %==========================================================================
        % 5.4 Plot fold locations and fold angle
        %==========================================================================
        % Plot theta vs time
        figure
        hold on
        plot(t,th1v,'-b','LineWidth',2)
        hold on
        plot(t, th2v,'-r','LineWidth',2)
        text(tf-50,-45,'t_d','FontSize',14,'Rotation',90,'FontWeight','bold')
%         text(te-50,-250,'t_e','FontSize',14,'Rotation',90,'FontWeight','bold')

        xlabel('Time [ms]')
        ylabel('\theta [deg]')
        grid on
        
%         plot(t(i1)*[1,1],[-50,50],'--k','LineWidth',2)
        plot(tf*[1,1],[-50,50],'--r','LineWidth',2)
        xticks([-600:200:2000])
        yticks([-90:15:90])
        xlim([0,tf])
        xlim([0,tmax])
        plot([0,0],[-60,60],'-k')
        plot(tf*[1,1],[-60,60],'-k')
        set(gcf,'Units','centimeters')
        set(gcf,'Position',[26,12,12,8])
        set(gca,'FontSize',14)
        
        % Plot fold location vs time
        figure
        hold on
        h = plot(t,fold1,'b','LineWidth',2);
        hold on
        g =plot(t, fold2,'r','LineWidth',2);
        xlabel('Time [ms]')
        ylabel('\xi')
        grid on
        xlim([0,tf])
        xlim([0,tmax])
        yticks([-1:0.25:1])
        plot(t(i1)*[1,1],[-1,1],'--k','LineWidth',2)
        plot(t(i2)*[1,1],[-1,1],'--r','LineWidth',2)
        plot([0,0],[-1,1],'-k')
        plot(tf*[1,1],[-1,1],'-k')
        set(gcf,'Units','centimeters')
        set(gcf,'Position',[26,2,12,8])
        set(gca,'FontSize',14)
        for JJ = 1: length(battens)
            plot([0,tmax],battens(JJ)*[1,1],'-k','LineWidth',1)
        end
        
    end
    %==========================================================================
    % 5.5 Load marker data
    %==========================================================================
    load(sprintf('%s\\Test %d\\test%d_run%d_coords.mat',pwd,strip,tN, rN))
    m1 = markers(:,1:3);
    m2 = markers(:,4:6);
    tmarker = tmarker - t0;
    tmarker = 1e3*tmarker(1:size(m1,1));
    m1f = m1(tmarker>0,:);
    m2f = m2(tmarker>0,:);
    m1f(m1f(:,3)>200,:)=nan;
    m2f(m2f(:,3)>200,:)=nan;
    m1(m1(:,3)>200,:)=nan;
    m2(m2(:,3)>200,:)=nan; 
    
    
    ind1 = find(m1(:,3)==max(m1f(:,3)),1);
    ind2 = find(m2(:,3)==max(m2f(:,3)),1);
    te = tmarker(max([ind1,ind2]));
    
    figure
    plot(tmarker, m1(:,3),'LineWidth',2)
    hold on
    plot(tmarker, m2(:,3),'LineWidth',2)
    plot(tf*[1,1],[-700,100],'--r','LineWidth',2);
    plot(te*[1,1],[-700,100],'--k','LineWidth',2);
    
    text(tf-50,-250,'t_d','FontSize',14,'Rotation',90,'FontWeight','bold')
    text(te-50,-250,'t_e','FontSize',14,'Rotation',90,'FontWeight','bold')
    xlabel('Time [ms]')
    ylabel('H [mm]')
    ylim([-700,100])
    set(gcf,'Units','centimeters')
        set(gcf,'Position',[38,2,12,8])
        set(gca,'FontSize',14)
    grid on
    xlim([0,1200])
    xlim([0,tmax])
    
    
    
    
    %==========================================================================
    % 5.6 Compute statistics
    %==========================================================================
    [initialFold1,foldPropagation1,...
        initialFold2,foldPropagation2] = deal(zeros(1,4));
    
    initialFold1(1:Nfolds1) = fold1(1,:);
    initialFold2(1:Nfolds2) = fold2(1,:);
    
    foldPropagation1(1:Nfolds1) = max(abs(fold1(t>0,:)-fold1(t==0,:)))
    foldPropagation2(1:Nfolds2) = max(abs(fold2(t>0,:)-fold2(t==0,:)))
    
    
    [Long1,Long2] = deal(table);
    Long1.initialAngle = th1v(1,:)
    Long1.initialFold = initialFold1;
    Long1.foldPropagation = foldPropagation1;
    
    Long2.initialAngle = th2v(1,:)
    Long2.initialFold = initialFold2;
    Long2.foldPropagation = foldPropagation2;
    
    deploymentTime = tf
    completeDeployment = max([min(abs(th1v)),min(abs(th2v))])<5;
    
    if update
        T = table_update(strip,tN,rN,Long1,Long2,deploymentTime,completeDeployment);
        
        % Save settings
        settings = struct;
        settings.thmode = thmode;
        settings.dlam = dlam;
        settings.kmin = kmin;
        settings.minInitialHeight = minInitialHeight;
        settings.dlam_maxdist = dlam_maxdist;
        settings.excludeRange = excludeRange;
        settings.max_disc = max_disc;
        settings.bypass_td_calculation = bypass_td_calculation;
        if bypass_td_calculation
            settings.te_ind = [i1,i2];
        else
            settings.te_ind = [0,0];
        end
        
        save(sprintf('%s\\Test %d\\exp_test%d_run%d.mat',pwd,strip,tN,rN),...
            'th1v','th2v','t','fold1','fold2','m1','m2','tmarker','tf','te','settings')
        
    end
    % end
