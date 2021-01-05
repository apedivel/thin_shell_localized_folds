function post1_extract_curvatures(strip,cases,realTimePlot)

% Strip length
L1 = 1125;
L2 = L1;

dt = 1;

% Check input
if nargin<3
    realTimePlot = true;
end

if nargin <2
%     cases = [1:43,1:43,1:43;
%         ones(1,43),2*ones(1,43),3*ones(1,43)]';
    tstart = 1;
    tend = 8;
    Nexp = tend-tstart+1;
    cases = [tstart:tend,tstart:tend,tstart:tend;
        ones(1,Nexp),2*ones(1,Nexp),3*ones(1,Nexp)]';
    cases(cases(:,1)==4 & cases(:,2)==3,:)=[];
    cases(end+1,:) = [2,4];
    cases = [5,1];
end

if nargin <1
    error('Specify the strip!')
end

% Number of experiment data to extract
N = size(cases,1);


% Iterate over each experiment
for CC = 1: N
    tN = cases(CC,1);
    rN = cases(CC,2);
    try
        %==========================================================================
        % 4.0 Load analysis
        %==========================================================================
        load(sprintf('%s\\Test %d\\test%d_run%d_coords.mat',pwd,strip, tN,rN))
        
        %==========================================================================
        % 4.1 Compute curvilinear abscissa for each point in the cloud
        % <<<<<<<<<<<<<<<<<<<<< User parameters >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        window = 50; % Number of points to pick for the moving window
        % movingPlot = false;
        
        % Initialization
        Nsteps = length(tdic);
        Np1 = size(X1dic,2);
        Np2 = size(X2dic,2);
        
        %--------------------------------------------------------------------------
        % 4.1.1 Compute midpoints
        %--------------------------------------------------------------------------
        tic
        X1m_noisy = movmedian([X1dic(1,:);Z1dic(1,:)],window,2, 'omitnan');
        X1m = smoothdata(X1m_noisy,2,'lowess',window);
        
        %""""""""""""""""""""""""""
        X2m_noisy = movmedian([X2dic(1,:);Z2dic(1,:)],window,2, 'omitnan');
        X2m = smoothdata(X2m_noisy,2,'lowess',window);
        toc
        
        %--------------------------------------------------------------------------
        % 4.1.2 Compute curvilinear abscissa
        %--------------------------------------------------------------------------
        dx = X1m(1,2:end)-X1m(1,1:end-1);
        dz = X1m(2,2:end)-X1m(2,1:end-1);
        
        s1m = zeros(Np1,1);
        
        for II = 2:Np1
            s1m(II) = s1m(II-1) + sqrt(dx(II-1)^2 + dz(II-1)^2) ;
        end
        
        % Find center of the strip
        Xc01 = markers(4,1);
        ind0 =  find(abs(X1m(1,:)-Xc01)==min(abs(X1m(1,:)-Xc01)));
        s1m = s1m - s1m(ind0(1));
        
        
        %""""""""""""""""""""""""""
        dx = X2m(1,2:end)-X2m(1,1:end-1);
        dz = X2m(2,2:end)-X2m(2,1:end-1);
        
        s2m = zeros(Np2,1);
        
        for II = 2:Np2
            s2m(II) = s2m(II-1) + sqrt(dx(II-1)^2 + dz(II-1)^2) ;
        end
        Xc02 = markers(1,1);
        ind0 =  find(abs(X2m(1,:)-Xc02)==min(abs(X2m(1,:)-Xc02)));
        s2m = s2m - s2m(ind0(1));
        
        
        %--------------------------------------------------------------------------
        % 4.1.3 Compute curvilinear abscissa for each point in the cloud
        %--------------------------------------------------------------------------
        s1 = zeros(Np1,1);
        t1 = zeros(2,Np1);
        f1 = false(Np1,1);
        % window = 10;
        tic
        
        for II = 1: Np1
            
            Xp = [X1dic(1,II); Z1dic(1,II)];
            if ~isnan(Xp)
                f1(1:end)  = 0;
                %------------------------------------------------------------------
                % 4.1.3.1 Find closest point Q on the midline
                dX = Xp - X1m;
                dist = sum(dX.^2,1);
                closest = find(dist==min(dist));
                Xq = X1m(:,closest);
                
                %------------------------------------------------------------------
                % 4.1.3.2 Find local orientation (using PCA)
                if closest <= window/2
                    f1(1 : closest + window/2) = 1;
                elseif closest > Np1 - window/2
                    f1(closest- window/2: Np1) =1 ;
                else
                    f1(closest - window/2 : closest + window/2) = 1;
                end
                R = pca(X1m(:,f1)');
                t1(:,II) = R(:,1);
                
                % CHECKPOINT: Check orientation of the tangent vector
                if t1(1,II)<0
                    t1(:,II) = -t1(:,II);
                end
                
                %------------------------------------------------------------------
                % 4.1.3.3 Project P on tangent vector
                s1(II) = s1m(closest) + dot(Xp-Xq,t1(:,II));
            end
        end
        
        %--------------------------------------------------------------------------
        % 4.1.3.4. Sort data by curvilinear abscissa
        %--------------------------------------------------------------------------
        [s1,ind1] = sort(s1);
        X1dic = X1dic(:,ind1);
        Y1dic = Y1dic(:,ind1);
        Z1dic = Z1dic(:,ind1);
        
        
        %""""""""""""""""""""""""""
        s2 = zeros(Np2,1);
        t2 = zeros(2,Np2);
        f2 = false(Np2,1);
        % window = 10;
        tic
        for II = 1: Np2
            
            Xp = [X2dic(1,II); Z2dic(1,II)];
            if ~isnan(Xp)
                f2(1:end)  = 0;
                %------------------------------------------------------------------
                % 4.1.3.1 Find closest point Q on the midline
                dX = Xp - X2m;
                dist = sum(dX.^2,1);
                closest = find(dist==min(dist));
                Xq = X2m(:,closest);
                
                %------------------------------------------------------------------
                % 4.1.3.2 Find local orientation (using PCA)
                if closest <= window/2
                    f2(1 : closest + window/2) = 1;
                elseif closest > Np2 - window/2
                    f2(closest- window/2: Np2) =1 ;
                else
                    f2(closest - window/2 : closest + window/2) = 1;
                end
                R = pca(X2m(:,f2)');
                t2(:,II) = R(:,1);
                
                % CHECKPOINT: Check orientation of the tangent vector
                if t2(1,II)<0
                    t2(:,II) = -t2(:,II);
                end
                
                %------------------------------------------------------------------
                % 4.1.3.3 Project P on tangent vector
                s2(II) = s2m(closest) + dot(Xp-Xq,t2(:,II));
            end
        end
        
        %--------------------------------------------------------------------------
        % 4.1.3.4. Sort data by curvilinear abscissa
        %--------------------------------------------------------------------------
        [s2,ind2] = sort(s2);
        X2dic = X2dic(:,ind2);
        Y2dic = Y2dic(:,ind2);
        Z2dic = Z2dic(:,ind2);
        toc
        

        
        %==========================================================================
        % 4.2 Compute curvatures and orientation at each time step
        %==========================================================================
        % <<<<<<<<<<<<<<<<<<<<< Parameters >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        window = 50;
        ds = 1;
        time = zeros(Nsteps,1);
        max_dist = 20;
        
        %--------------------------------------------------------------------------
        % Define new curvilinear abscissa with uniform spacing
        %--------------------------------------------------------------------------
        s1m = [s1m(1):ds:s1m(end)];
        s2m = [s2m(1):ds:s2m(end)];
        Np1 = length(s1m);
        Np2 = length(s2m);

  
        %-------------------------------------
        [t1,t10] = deal(zeros(2,Np1));
        [f1,f10] = deal(false(Np1,1));
        [k1,k1_raw] = deal(zeros(Nsteps,Np1-2));
        [th1,th1_raw] = deal(zeros(Nsteps,Np1));
        
        [t2,t20] = deal(zeros(2,Np2));
        [f2,f20] = deal(false(Np2,1));
        [k2,k2_raw] = deal(zeros(Nsteps,Np2-2));
        [th2,th2_raw] = deal(zeros(Nsteps,Np2));
        
        X1m_noisy = zeros([2,Np1]);
        X2m_noisy = zeros([2,Np2]);
        
        closest1 = abs(s1m-s1)<=max_dist;
        closest2 = abs(s2m-s2)<=max_dist;

        for JJ  =1 :dt: Nsteps
            t1 = t10;
            t2 = t20;
            %----------------------------------------------------------------------
            % 4.2.1 Compute coordinates of the midline
            %-------------------------------------------------------------------
            for KK = 1: Np1
                X1m_noisy(:,KK) = median([X1dic(JJ,closest1(:,KK)); Z1dic(JJ,closest1(:,KK))],2,'omitnan');
            end
            X1m = smoothdata(X1m_noisy,2,'lowess',window);

            for KK = 1: Np2
                X2m_noisy(:,KK) = median([X2dic(JJ,closest2(:,KK)); Z2dic(JJ,closest2(:,KK))],2,'omitnan');
            end

            X2m = smoothdata(X2m_noisy,2,'lowess',window);
            
            tic
            for II = 1: Np1
                f1 = f10;
                %------------------------------------------------------------------
                % 4.2.2 Get points on the midline in the vicinity of II-th point
                %------------------------------------------------------------------
                if II <= window/2
                    f1(1 : II + window/2) = 1;
                elseif II > Np1 - window/2
                    f1(II- window/2: Np1) =1 ;
                else
                    f1(II - window/2 : II + window/2) = 1;
                end
                
                %-------------------------------------------------------------------
                % 4.2.3 Local tangent
                %-------------------------------------------------------------------
                R = pca(X1m(:,f1)');
                
                % CHECKPOINT: Enough points are being used for the tangent
                % estimation
                if size(R) == [2 2]
                    t1(:,II) = R(:,1);
                    
                    % CHECKPOINT: Check orientation of the tangent vector
                    if t1(1,II)<0
                        t1(:,II) = -t1(:,II);
                    end
                end
            end
                    %--------------------------------------------------------------
                    % 4.2.4 Local orientation
                    %--------------------------------------------------------------
                    th1_raw(JJ,:) = atan2d(t1(2,:),t1(1,:));
                    th1(JJ,:) = th1_raw(JJ,:);
%                         th1(JJ,:) = smoothdata(th1_raw(JJ,:),'lowess',window);
                    
                    %--------------------------------------------------------------
                    % 4.2.5 Longitudinal curvature
                    %--------------------------------------------------------------
%                     s1_int = [s1m(JJ,1):1:s1m(JJ,end)];
%                     t1_int(1,:) = interp1(s1m(JJ,:),t1(1,:),s1_int);
%                     t1_int(2,:) = interp1(s1m(JJ,:),t1(2,:),s1_int);
                    
%                     k1_raw(JJ,:) = vecnorm(t1(:,3:end)-t1_int(:,1:end-2),1)./(s1_int(JJ,3:end)-s1_int(JJ,1:end-2));
                    k1_raw(JJ,:) = vecnorm(t1(:,3:end)-t1(:,1:end-2),1)./(s1m(3:end)-s1m(1:end-2));
                    k1(JJ,:) = smoothdata(k1_raw(JJ,:),'lowess',window);
%                     k1(JJ,:) = k1_raw(JJ,:);
%                     k1(JJ,:) = smoothdata(k1_raw(JJ,:),'lowess',window);
%                 else
%                     th1(JJ,:) = nan(1,Np1);
%                     k1(JJ,:) = nan(1,Np1-2);
%                     
%                 end
%             end
            
            
            for II = 1: Np2
                f2 = f20;
                %------------------------------------------------------------------
                % 4.2.2 Get points on the midline in the vicinity of II-th point
                %------------------------------------------------------------------
                if II <= window/2
                    f2(1 : II + window/2) = 1;
                elseif II > Np2 - window/2
                    f2(II- window/2: Np2) =1 ;
                else
                    f2(II - window/2 : II + window/2) = 1;
                end
                
                %-------------------------------------------------------------------
                % 4.2.3 Local tangent
                %-------------------------------------------------------------------
                R = pca(X2m(:,f2)');
                
                % CHECKPOINT: Enough points are being used for the tangent
                % estimation
                if size(R) == [2 2]
                    t2(:,II) = R(:,1);
                    
                    % CHECKPOINT: Check orientation of the tangent vector
                    if t2(1,II)<0
                        t2(:,II) = -t2(:,II);
                    end
                end
            end
                    %--------------------------------------------------------------
                    % 4.2.4 Local orientation
                    %--------------------------------------------------------------
                    th2_raw(JJ,:) = atan2d(t2(2,:),t2(1,:));
                    th2(JJ,:) = th2_raw(JJ,:);
%                         th2(JJ,:) = smoothdata(th1_raw(JJ,:),'lowess',window);
                    
                    %--------------------------------------------------------------
                    % 4.2.5 Longitudinal curvature
                    %--------------------------------------------------------------
                    k2_raw(JJ,:) = vecnorm(t2(:,3:end)-t2(:,1:end-2),1)./(s2m(3:end)-s2m(1:end-2));
                    k2(JJ,:) = k2_raw(JJ,:);
                    %                         k2(JJ,:) = smoothdata(k2_raw(JJ,:),'lowess',window);
%                 else
%                     th2(JJ,:) = nan(1,Np2);
%                     k2(JJ,:) = nan(1,Np2-2);
%                     
%                 end
%             end
            
            % USER INTERFACE: Plot geometry, interpolation, curvature and
            % orientation
            if realTimePlot
                
                subplot(3,2,1)
                plot(X1dic(JJ,:), Z1dic(JJ,:),'.k')
                hold on
                plot(X1m(1,:), X1m(2,:),'r')
                hold off
                axis([-600 600 -600 100])
                
                
                subplot(3,2,3)
                plot(s1m,th1(JJ,:))
                hold on
                plot(s1m,th1_raw(JJ,:))
                axis([-600 600 -60 60])
                hold off
                
                subplot(3,2,5)
                plot(s1m(2:end-1), k1(JJ,:))
                xlim([-600,600])
                
                
                subplot(3,2,2)
                plot(X2dic(JJ,:), Z2dic(JJ,:),'.k')
                hold on
                plot(X2m(1,:), X2m(2,:),'r')
                hold off
                axis([-600 600 -600 100])
                
                
                subplot(3,2,4)
                plot(s2m,th2(JJ,:))
                hold on
                plot(s2m,th2_raw(JJ,:))
                axis([-600 600 -60 60])
                hold off
                
                subplot(3,2,6)
                plot(s2m(2:end-1), k2(JJ,:))
                xlim([-600,600])
                
                drawnow
            end
            
            % USER INTERFACE: print elapsed time
            time(JJ) = toc;
            disp(sprintf('Test %d run %d - Step %d of %d completed... t = %.1f s',tN,rN,JJ,Nsteps,sum(time)))
            
            
        end
        
        % USER INTERFACE: print total time
        disp(sprintf('Curvature extraction completed! Total elapsed time %d min %.0f s',...
            floor(sum(time)/60), mod(sum(time),60)))
        
        
        % Remove extra points in the orientation matrix
        th1 = th1(:,2:end-1);
        th2 = th2(:,2:end-1);
        
        % 5.1.1 Non dimensional length
%         L1 = max(s1)-min(s1);
        lambda1 = 2*s1m(:,2:end-1)/L1;
        
        
%         L2 = max(s2)-min(s2);
        lambda2 = 2*s2m(:,2:end-1)/L2;
        k1 = smoothdata(k1,2,'lowess',window);
        k2 = smoothdata(k2,2,'lowess',window);
        
        
        % Remove empty rows when downsampling
        tdic = tdic(1:dt:end);
        lambda1 = lambda1(1:dt:end,:);
        lambda2 = lambda2(1:dt:end,:);
        k1 = k1(1:dt:end,:);
        k2 = k2(1:dt:end,:);
        th1 = th1(1:dt:end,:);
        th2 = th2(1:dt:end,:);
        
        save(sprintf('%s\\Test %d\\post_test%d_run%d.mat',pwd,strip,tN,rN),'lambda1','lambda2','s1','s2','k1','k2','th1','th2','tdic')
    catch
        
        warning(['Test %d run %d could not be extracted! Manual check required.'])
        
    end
    
end
end
