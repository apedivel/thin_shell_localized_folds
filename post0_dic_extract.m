clc, clear, close all


%*************************************************************************
% Parameters
%*************************************************************************
strip = 3;


% Input path (location of mat files from VIC 3D)
path = sprintf('G:\\HS camera\\Experiment %d\\',strip);
% path = sprintf('C:\\Users\\anto_\\OneDrive\\Desktop\\HS\\',strip);

% Output path (location of processed mat files)
output_path = sprintf('%s\\Test %d\\',pwd,strip);

% overwrite = false;
overwrite = true;

% Select range of tests
tstart = 6;
tend = 6;

% % Get camera frequency from readme file
readme = load([output_path,'\readme.dat']);
f = readme(:,5); % Frequency

% f = 500*ones(4,1);
dt = 1./f;



for tN = tstart:tend
    for rN = 1:1
        folder = sprintf('test%d_run%d\\',tN,rN);
        %--------------------------------------------------------------------------
        % Load marker coordinates
        %--------------------------------------------------------------------------
        cs_data = csvread([path,folder,'coords.csv'],2,0);
        
        % Extract coordinate system markers
        cs_raw = reshape(cs_data(1,1:9),3,[]);
        
        % Origin of the coordinate system
        or = cs_raw(:,1);
        
        % Translate coordinate system points
        cs = cs_raw - or;
        
        % Define X axis (along direction connecting the pulleys)
        X = cs(:,2)/norm(cs(:,2));
        
        % Define Y axis (orthonormal to X)
        Y = cs(:,3) - dot(cs(:,3),X)*X;
        Y = Y/norm(Y);
        
        % Define Z axis
        Z = cross(X,Y);
        
        % Rotation matrix to change coordinates
        R = inv([X,Y,Z]);
        
        % Markers
        marker_raw = cs_data(:,10:end);
        markers = zeros(size(marker_raw));
        N = size(marker_raw,1);
        Np = size(marker_raw,2)/3;
        for II = 1: Np
           markers(:,3*(II-1)+1:3*II) = (R*(marker_raw(:,3*(II-1)+1:3*II)'- or ))';
        end
        tmarker = [0:dt(tN):Np-1]';
        
        
        % Find all files in the folder
        files = dir([path,folder]);
        
        % Initialize files
        count = 0;
        [X1dic,Y1dic,Z1dic,X2dic,Y2dic,Z2dic] = deal([]);
        for II = 1: length(files)
            % Check if file contains point cloud data
            file = files(II).name;
            ext = strsplit(file,'.');
            fn = ext{1};
            ext = ext{2};
            
            if ~isempty(ext) &&  strcmp(ext,'mat') && strcmp(fn(1:3),'dep')
                count = count + 1;
                load([path,folder,file])
                
                % Longeron 1 coordinates (in global frame)
                x1 = X - or(1) + U;
                y1 = Y - or(2) + V;
                z1 = Z - or(3)+ W;
                
                % Longeron 2 coordinates (in global frame)
                x2 = X_0 - or(1) + U_0;
                y2 = Y_0 - or(2) + V_0;
                z2 = Z_0 - or(3) + W_0;
                
                % Remove void points
                filter1 = Z~=max(max(Z));
                filter2 = Z_0~=max(max(Z_0));
                
                % Reshape data into a vector
                x1 = reshape(x1(filter1),1,[]);
                y1 = reshape(y1(filter1),1,[]);
                z1 = reshape(z1(filter1),1,[]);
                
                x2 = reshape(x2(filter2),1,[]);
                y2 = reshape(y2(filter2),1,[]);
                z2 = reshape(z2(filter2),1,[]);
                
                
                % Rotate points in new coordinate system
                c1 = R*[x1;y1;z1];
                x1rot = c1(1,:);
                y1rot = c1(2,:);
                z1rot = c1(3,:);
                
                c2 =R*[x2;y2;z2];
                x2rot = c2(1,:);
                y2rot = c2(2,:);

                z2rot = c2(3,:);
                
                X1dic(count,:) = x1rot;
                Y1dic(count,:) = y1rot;
                Z1dic(count,:) = z1rot;
                
                X2dic(count,:) = x2rot;
                Y2dic(count,:) = y2rot;
                Z2dic(count,:) = z2rot;
                
                
                % Remove outliers
                if count>1
                    lost_points1 = U==Uold1 | V==Vold1 | W==Wold1;
                    lost_points1 = reshape(lost_points1(filter1),1,[]);
                    lost_points2 = U_0==Uold2 | V_0==Vold2 | W_0==Wold2;
                    lost_points2 = reshape(lost_points2(filter2),1,[]);
                    
                    X1dic(count,lost_points1) = nan;
                    Y1dic(count,lost_points1) = nan;
                    Z1dic(count,lost_points1) = nan;
                    X2dic(count,lost_points2) = nan;
                    Y2dic(count,lost_points2) = nan;
                    Z2dic(count,lost_points2) = nan;

                end
                
                

                
                Uold1 = U;
                Vold1 = V;
                Wold1 = W;
                
                Uold2 = U_0;
                Vold2 = V_0;
                Wold2 = W_0;

            end
            

            
        end
        
        tdic = [0:dt(tN):(count-1)*dt(tN)];
        
        % Final filter
%         filter1 = Z1dic(1,:)<-500 | Y1dic(1,:)>450;
%         filter2 = Z2dic(1,:)<-500;
        filter1 = Y1dic(1,:)>450 | X1dic(1,:)<-600;
        filter2 = Y2dic(1,:)>450| X2dic(1,:)<-600;

        X1dic(:,filter1) = [];
        Y1dic(:,filter1) = [];
        Z1dic(:,filter1) = [];
        X2dic(:,filter2) = [];
        Y2dic(:,filter2) = [];
        Z2dic(:,filter2) = [];
        

        
        filename = folder(1:end-1);
        if ~exist([output_path,filename,'_coords.mat'],'file') || overwrite
            
            save([output_path,filename,'_coords.mat'],'tdic','X1dic','Y1dic','Z1dic',...
                'X2dic','Y2dic','Z2dic','tmarker','markers','or','R')
            disp(sprintf('%s exported!',filename))
        else
            disp(sprintf('%s already exists...',filename))
        end

        
        % Plot time sequence
  
        N = min([size(X1dic,1),size(markers,1)]);

        
        for II =1:N
            plot3(X1dic(II,:),Y1dic(II,:),Z1dic(II,:),'.')
            hold on
            plot3(X2dic(II,:),Y2dic(II,:),Z2dic(II,:),'.')
            mark = reshape(markers(II,:),3,[]);
            plot3(mark(1,:),mark(2,:),mark(3,:),'ok','LineWidth',3)
            hold off
            grid on
            axis equal
            xlim([-600,600])
            ylim([0,600])
            zlim([-800,100])
            pause(0.001)
        end
        hold off
    end
end