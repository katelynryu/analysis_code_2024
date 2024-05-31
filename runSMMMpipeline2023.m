%% PIPELINE
% Metadata
fps=80;

pixel_size = 1/1.97; % mm
noi = [1,2,5,6,7,8,10,12]; % index body parts of interest
threshold_quantile=0.98;
% colors for paws
colors_paws = [205,32,41;61,80,157;166,75,156;142,205,212]./255; %  RF,LF,RH,LH

%TSC1Sbase = 'PR-OFT-TSC1-social-dbase-01022023.mat';
% TSC1Sbase = 'SMALL-new-dbase-012024.mat';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing of SLEAP predicitions 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PreStep1: CREAT DATABASE as a matlab structure

%d = uigetdir([],'Pick a root directory');
d = 'Z:\behdata\2021-09-TSC1-social\SocialTsc1\analysisFiles\SMALL_last2';
%path_vid = 'Y:\behdata\2021-09-TSC1-social\';
path_vid = 'Z:\behdata\2023-09-OFTsocialapproach\';
files = dir(fullfile(d,'**\*.*'));
files = arrayfun(@(x) fullfile(x.folder,x.name),files,'UniformOutput',false);
files = files(endsWith(files,'.h5'));
%Name of file to hold the database
%TSC1Sbase = 'PR-OFT-TSC1-social-dbase-01022023.mat';
TSC1Sbase = 'SMALL-new-dbaselast2.mat';

%Draw each session in directory into database
for i= 1:length(files)
    %Only selecting the text files for analysis
    if exist(TSC1Sbase,'file')
        %load(TSC1Sbase)
        newIdx = length(dbaseA)+1;
        %Check if file already loaded and plot state
        currentFnames = {dbaseA.fname};
    else
        dbaseA = struct;
        newIdx = 1;

        currentFnames = [];
    end

    fname = files{i};

    %Get animal and identifiers
    [fpath,file,~] = fileparts(fname);
    identifiers = strsplit(file,'-');
    %fileID = strcat(identifiers{2}(1:4),'-',identifiers{3}(1:2));
    fileID = strcat(identifiers{2}(1:4));
    identifiers = strsplit(fpath,'\');
    condition = identifiers{end};


    sleapout_nam = fname;% [d,'\',file,'.h5'];
    tracks = h5read(sleapout_nam, '/tracks');
    node_names = h5read(sleapout_nam, '/node_names');
    occupancy = h5read(sleapout_nam, '/track_occupancy');
    point_scores = h5read(sleapout_nam,'/point_scores');
    tracking_scores = h5read(sleapout_nam,'/tracking_scores');

    meta = h5info(sleapout_nam);

    fprintf(1, 'File: %s \n',sleapout_nam);
    fprintf(1, 'Number of tracks: %d \n',size(tracks,4))
    fprintf(1, 'Number of body parts tracked: %d \n',size(node_names,1))
    fprintf(1, 'With scores for predictions: %d \n',size(meta.Datasets,1)>4)

    for ii=1:size(tracks,4)
        fprintf(1, 'Track #%d - frames with predictions %d / %d \n' ,ii ,sum(occupancy(ii,:)),size(tracks,1));
        for n=1:length(noi)
            fprintf(1, '#%s - frames with NaNs %d / %d \n' ,node_names(noi(n),1) ,sum(isnan(tracks(:,noi(n),1,ii))),size(tracks,1));
        end
    end
     time_frames = h5read([path_vid,'OFTsocialgroup-',fileID,'.h5'], '/pg0_time');% '-00.h5' for WTWT
     time_frames = time_frames';


    if any(strcmp(fname,currentFnames))
        return;
    else
        %Pack in new dbase element and save
        dbaseA(newIdx).fileID = fileID;
        dbaseA(newIdx).fname = fname;

        dbaseA(newIdx).condition = condition;
        dbaseA(newIdx).SLEAPtracks = tracks;
        dbaseA(newIdx).SLEAPscores = point_scores;
        dbaseA(newIdx).time = time_frames;
        dbaseA(newIdx).meta.behvideo = [path_vid,'OFTsocialgroup-',fileID,'.h5']; %'-00.h5' for WTWT
        dbaseA(newIdx).meta.node_names = node_names;
        dbaseA(newIdx).meta.tracking_scores = tracking_scores;
        dbaseA(newIdx).meta.occupancy = occupancy;

        %save(TSC1Sbase,'dbaseA');
    end
    clear fileID file fname fpath condition time_frames trials
end
%% VERIFY DATABASE and quality of SLEAP predictions
load(TSC1Sbase)

%% PreStep2: Creating filtered tracks and finding frames with good predicitons for traning

for i=1:length(dbase)
    dbase(i).tracks = zeros (size(dbase(i).SLEAPtracks));
    for m = 1:size(dbase(i).SLEAPtracks,4) %mouse
        for p=1:size(dbase(i).SLEAPtracks,2)%noi %body part

            for c = 1:size(dbase(i).SLEAPtracks,3) %coordinate
                nans = isnan(dbase(i).SLEAPtracks(:,p,c,m)); %finding nans
                x_raw = dbase(i).SLEAPtracks(:,p,c,m); % get one position over time
                cutoff=quantile(abs(diff(x_raw)),threshold_quantile); %  quantile threshold for between-frame position displacement
                [x_new]=filter_position(x_raw,dbase(i).time,cutoff); % The outliers are replaced by linear interpolations
                dbase(i).tracks(:,p,c,m) = x_new;

            end
        end
    end
end
clear p s i c m

%% PreStep3: Eliminating ghost mice

for i=1:length(dbase)
    for m = 1:size(dbase(i).SLEAPtracks,4) %mouse
        for p=1:size(dbase(i).SLEAPtracks,2)%noi %body part
            for c = 1:size(dbase(i).SLEAPtracks,3) %coordinate
                % eliminating ghost mice


                %                 dbase(i).tracks(:,p,c,m) = temp_yline;
                %                 % y axis
                %                 temp_yline = dbase(i).tracks(:,p,c,m);
                %                 temp_yline(temp_yline>1100 )=NaN;
                %                 temp_yline(temp_yline<180 )=NaN;
                %                 dbase(i).tracks(:,p,c,m) = temp_yline;

                %x and yaxis
                temp_xline = dbase(i).tracks(:,p,c,m);
                if c==1 % x axis
                    temp_xline(temp_xline>970 )=NaN;
                    temp_xline(temp_xline<70 )=NaN;

                    %making sure other coordinate is also nan
                    otherc = dbase(i).tracks(:,p,2,m);
                    otherc(isnan(temp_xline))=NaN;
                    dbase(i).tracks(:,p,2,m)=otherc;

                elseif c==2 % y axis
                    temp_xline(temp_xline>1100 )=NaN;
                    temp_xline(temp_xline<180 )=NaN;

                    %making sure other coordinate is also nan
                    otherc = dbase(i).tracks(:,p,1,m);
                    otherc(isnan(temp_xline))=NaN;
                    dbase(i).tracks(:,p,1,m)=otherc;
                end
                dbase(i).tracks(:,p,c,m)=temp_xline;
                %                 %convert that WHOLE mouse to nan
                %                %x and yaxis
                %                 temp_xline = dbase(i).tracks(:,p,c,m);
                %                 if c==1 % x axis
                %                     temp_xline(temp_xline>970 )=-999;
                %                     temp_xline(temp_xline<70 )=-999;
                %                 elseif c==2 % y axis
                %                     temp_xline(temp_xline>1100 )=-999;
                %                     temp_xline(temp_xline<180 )=-999;
                %                 end
                %                 if sum(temp_xline==-999)>0 %if there is ghost mouse
                %                     for noi=1:size(dbase(i).SLEAPtracks,2)
                %                         for coor = 1:2
                %                         temp=dbase(i).tracks(:,noi,coor,m);
                %                         temp(temp_xline==-999)=NaN;
                %                         dbase(i).tracks(:,noi,coor,m)=temp;
                %                         end
                %                     end
                %                 else
                %                     dbase(i).tracks(:,p,c,m)=temp_xline;
                %                 end
            end
        end
    end
end
clear p s i c m

%% PreStep4: Fill small gaps in filtered predictions with linear interpolation
for i=1:length(dbase)
    for m = 1:size(dbase(i).SLEAPtracks,4) %mouse
        for p=1:size(dbase(i).SLEAPtracks,2)%noi %body part
            for c = 1:size(dbase(i).SLEAPtracks,3) %coordinate
                nans = isnan(dbase(i).SLEAPtracks(:,p,c,m)); %finding nans
                dbase(i).meta.frNoSLEEAP(p,m) =  sum(nans)/length(nans); % fraction of nans per position
                stretch =bwconncomp(nans); %finding nan stretches
                length_stretch = zeros(1,stretch.NumObjects);
                for s=1:length(stretch.PixelIdxList)
                    %length_stretch(s,1) = length(stretch.PixelIdxList{1,s});
                    length_stretch(1,s) = length(stretch.PixelIdxList{1,s});
                end

                %LESS THAN 3 - FILL MISSING
                %finding frame index
                idx_temp = find(length_stretch<3);
                frames_idx = cell2mat(stretch.PixelIdxList(1,idx_temp)');
                fm = fillmissing(dbase(i).SLEAPtracks(:,p,c,m),'linear');
                dbase(i).tracks(frames_idx, p, c, m) = fm(frames_idx,1);
                clear frames_idx idx_temp fm

                %3 frames TO 1 sec - reconstruction with fillgaps local
%                 idx_temp  = find(length_stretch>3 & length_stretch<fps);
%                 frames_idx = cell2mat(stretch.PixelIdxList(1,idx_temp)');
%                 fm1 = fillgaps(dbase(i).SLEAPtracks(:,p,c,m),20*fps);
%                 dbase(i).tracks(frames_idx, p, c, m) = fm1(frames_idx,1);

                clear frames_idx idx_temp fm1 x_new x_raw
            end
        end
        clear frames_idx nans
        dbase(i).badframes(:,m) = any(isnan(squeeze(dbase(i).tracks(:, noi,2,m))),2);
        %define box, make nan for outside
        %dbase(i).badframes(:,m) = any(isnan(squeeze(dbase(i).tracks(:, noi(1:end-1),2,m))),2);% without tail tip
        dbase(i).BadFract(1,m)=sum(dbase(i).badframes(:,m))/length(dbase(i).badframes);

    end
end

clear p s i c m


% % %% finding touch frames (ghost)
% % noi=[1 2 10 13]; % nose, neck, tail base, mid body
% % for i=33:length(dbase)
% %     dbase(i).touch = zeros(size(dbase(i).SLEAPtracks(:,1)));
% %     for m = 1:size(dbase(i).SLEAPtracks,4) %mouse
% %         % for p=1:size(dbase(i).SLEAPtracks,2)%noi %body part
% %         % for c = 1:size(dbase(i).SLEAPtracks,3) %coordinate
% %         for f = 1:size(dbase(i).SLEAPtracks,1)
% %             %if one track is nan - touch state
% %             if sum(isnan(dbase(i).tracks(f,:,1,1)))>5
% %                 dbase(i).touch(f,1)=1;
% %             elseif sum(isnan(dbase(i).tracks(f,:,1,2)))>5
% %                 dbase(i).touch(f,1)=1;
% %             end
% %             %if they are closer than 50 pixels
% %             %             for n1=noi
% %             %                 for n2=noi
% %             %                     if returnDist(squeeze(dbase(i).tracks(f, n1,:,1)), ...
% %             %                             squeeze(dbase(i).tracks(f, n2,:,2))) < 50
% %             %                         dbase(i).touch(f,1)=1;
% %             %                     end
% %             %                 end
% %             %             end
% %         end
% %     end
% % end

%% PreStep5: finding touch frames (overlap)

nodes=[1 2 3 4 5 6 7 8 9 10 13]; %all body parts except tail tip and midtail
for i=1:length(dbase) %video
    dbase(i).touch = zeros(size(dbase(i).SLEAPtracks(:,1)));
    for f = 1:size(dbase(i).SLEAPtracks,1) %frames
        %finding coordinates
        %mouse 1
        x1=rmmissing(dbase(i).tracks(f,nodes,1,1)'); %rmmissing takes out nans
        y1=rmmissing(dbase(i).tracks(f,nodes,2,1)');

        %mouse 2
        x2=rmmissing(dbase(i).tracks(f,nodes,1,2)');
        y2=rmmissing(dbase(i).tracks(f,nodes,2,2)');

        if length(x1)>2 && length(x2)>2
            %mouse 1
            k1=boundary(x1,y1, 0.3);
            p1=polyshape(x1(k1), y1(k1));

            %mouse 2
            k2=boundary(x2,y2, 0.3);
            p2=polyshape(x2(k2), y2(k2));

            %polyvec=[p1 p2];
            TF = overlaps(p1, p2); %determines of shapes overlap
            if TF ==1
                dbase(i).touch(f,1)=1;
            end
        else
        end

    end
end
% percent of touch frames
for i=1:length(dbase)
    dbase(i).touchFract = (sum(dbase(i).touch/length(dbase(i).touch)))*100;
end

%% PreStep6: exlude bad and touch frames 
%combining bad and touch
for i=1:length(dbase)
    tempbad = or( logical(dbase(i).badframes(:,1)),  logical(dbase(i).badframes(:,2)));
    dbase(i).excludedFrames = or(tempbad, logical(dbase(i).touch));
end

%% CHECK - plotting distance histograms
%CHECK - distance distribution without bad and touch frames
conditions = nominal({dbase.condition});
unicon = unique(conditions);
cond_dist = cell(length(unicon),2);
%plotting distance matrices
for j=1:length(unicon)
    temp_idx = find(conditions==unicon(j));
    alldist = [];
    for i=temp_idx
        tempdist=dbase(i).AAdist/(1.97*10);
        fdist = tempdist(dbase(i).excludedFrames==0); %filtered dist
        alldist = [alldist; fdist];

    end
    subplot(length(unicon),1,j)
    histogram(alldist,'Normalization','pdf')
    xlim([0 60])
end

%% CHECK - finding min and maxs
minmax = zeros(32,4);
for i=1:length(dbase)
    minmax(i,1) = min(min(min(dbase(i).tracks(:,:,1,:))));
    minmax(i,2) = max(max(max(dbase(i).tracks(:,:,1,:))));
    minmax(i,3) = min(min(min(dbase(i).tracks(:,:,2,:))));
    minmax(i,4) = max(max(max(dbase(i).tracks(:,:,2,:))));
end
%% CHECK - plotting all tracks to find outliers
noi=10;
for i=33:length(dbase)
    plot(dbase(i).tracks(:,noi,1,1), dbase(i).tracks(:,noi,2,1))
    hold on;
end
% hold on; xline(970)
% hold on; xline(60)
% hold on; yline(170)
% hold on; yline(1100)
%% inspect filtered tracks
for i=33%:length(dbase)
    figure
    subplot(2,1,1)
    plot(dbase(i).time-dbase(i).time(1,1),dbase(i).SLEAPtracks(:, noi,2,1),"LineWidth",1)
    ylim([0 1200])
    prepfig
    set(gca,'PlotBoxAspectRatio',[1 0.35 1])
    subplot(2,1,2)
    plot(dbase(i).time-dbase(i).time(1,1),dbase(i).SLEAPtracks(:, noi,2,2),"LineWidth",1)
    ylim([0 1200])
    xlabel('Time, sec')
    prepfig
    set(gca,'PlotBoxAspectRatio',[1 0.35 1])
    sgtitle(dbase(i).condition)
    %saveas(gcf,['Tracks_filtered_',dbase(i).fileID,'.pdf'])
    close(gcf)

    figure
    timex = dbase(i).time-dbase(i).time(1,1);
    plot3(dbase(i).SLEAPtracks(1:9400, 10,1,1),dbase(i).SLEAPtracks(1:9400, 10,2,1),timex(1:9400),'b',"LineWidth",1)
    hold on
    plot3(dbase(i).SLEAPtracks(1:9400, 10,1,2),dbase(i).SLEAPtracks(1:9400, 10,2,2),timex(1:9400),'r',"LineWidth",1)
    prepfig
    set(gca,'PlotBoxAspectRatio',[1 1 3])
    xlabel('Position,px')
    zlabel('Time, sec')
    title(dbase(i).condition)
    %saveas(gcf,['Centroids_time_',dbase(i).fileID,'.pdf'])
    %close(gcf)
end

%% Basic analysis1: centroid analysis
for t=1:length(dbase)
    dbase(t).AAdist = zeros(size(dbase(t).time));
    for i=1:length(dbase(t).time)
        dbase(t).AAdist(i,1) = returnDist(squeeze(dbase(t).tracks(i, 10,:,1)),squeeze(dbase(t).tracks(i, 10,:,2)));
    end
end

% creat mean pdf of distances

conditions = nominal({dbase.condition});
unicon = unique(conditions);
cond_dist = cell(length(unicon),2);
for i=1:length(unicon)
    temp_idx = find(conditions==unicon(i));
    cond_dist{i,1} = string(unicon(i));
    for j=1:length(temp_idx)
        cond_dist{i,2}(1,j) = median(dbase(temp_idx(1,j)).AAdist,'omitnan'); % distance between pairs
        % how many pixels away
        cond_dist{i,3}(1,j) = sum(dbase(temp_idx(1,j)).AAdist>550)/(length(dbase(temp_idx(1,j)).time)-sum(dbase(temp_idx(1,j)).badframes(:,1)));
        h=histogram(dbase(temp_idx(1,j)).AAdist(~isnan(dbase(temp_idx(1,j)).AAdist),1),100, 'Normalization','pdf');
        temp_h= h.Values;
        cond_dist{i,4}(:,j) = temp_h;
        clear temp_h h

    end
end

%creating density plot
plot(mean(cond_dist{1, 4}'),'r')
hold on
plot(mean(cond_dist{3, 4}'),'b')
prepfig;
plot(mean(cond_dist{2, 4}'),'g')
ylabel('Cumulative density')
xlabel('Distance')

% creat mean pdf of distances
figure
for i=1:6
    subplot(3,2,i)
    bar(mean(cell2mat(cond_dist(i, 4)),2), 'k','EdgeColor','k')
    title(cond_dist(i,1))
    prepfig;
    set(gca,'PlotBoxAspectRatio',[1 0.5 1])
end

% plot franction in non-social
X = [cell2mat(cond_dist(1,3))'; cell2mat(cond_dist(2,3))';cell2mat(cond_dist(3,3))'];
grp = [ones(size(cond_dist{1,2}))'; 2*ones(size(cond_dist{2,2}))';  3*ones(size(cond_dist{3,2}))'];

[~,~,stats] = anovan(X, {grp});
multcompare(stats);
figure;
boxplot(X,grp,'Widths',0.3)
ylabel('Fraction of time with no social interactions')
yt = get(gca, 'YTick');
axis([xlim    min(yt)-min(yt)*0.2  max(yt)*1.2])
xt = get(gca, 'XTick');
prepfig;
xticklabels(cond_dist(:,1))
hold on
c= zeros(length(grp),3);
c(grp==1,1:3) =  repmat([0 0.4470 0.7410],sum(grp==1),1);
c(grp==2,1:3) =  repmat([0.8500 0.3250 0.0980],sum(grp==2),1);
c(grp==3,1:3) =  repmat([0.9290 0.6940 0.1250],sum(grp==3),1);

swarmchart(grp,X,25,c,'filled','MarkerFaceAlpha',0.5)
plot(xt([1 3]), [1 1]*max(yt)*1.1, '-k',  mean(xt([1 3])), max(yt)*1.15, '*k')
set(gca,'PlotBoxAspectRatio',[1 2 1])



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate training dataset
% by smartsampling from the single mouse UMAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trainset Step1: reembed to WT single mouse space

%%
load('Z:\mkislin\DTSC\cohort_optoUR\trainingser_sLEAPWT.mat') 
load('Z:\mkislin\DTSC\cohort_optoUR\vecsVals_sLEAP-WT.mat')
%load('Y:\mkislin\DTSC\cohort_optoUR\template_WT.mat')
load( 'Y:\behdata\2021-09-TSC1-social\SocialTsc1\mat_files\template_WT.mat')
%load( 'Z:\behdata\2021-09-TSC1-social\SocialTsc1\temp-SMALL-new-upto31.mat')
%%
fps=80;

pixel_size = 1/1.97; % mm
noi = [1,2,5,6,7,8,10,12]; % index body parts of interest
threshold_quantile=0.98;
% colors for paws
colors_paws = [205,32,41;61,80,157;166,75,156;142,205,212]./255; %  RF,LF,RH,LH


% here we select body parts for futher analysis
noi = [1,2,5,6,7,8,10,12];
xIdx = noi;
yIdx = noi;

%  prepare matrix to calculate distances from body part positions
[X, Y] = meshgrid(xIdx,yIdx);
X = X(:); Y = Y(:);
IDX = find(X~=Y);

numProjections = 10; 
numModes = 10; pcaModes = numModes;
minF = .25; maxF = 20;

parameters = setRunParameters([]);
parameters.trainingSetSize = 400;
parameters.pcaModes = pcaModes;
parameters.samplingFreq = fps;
parameters.minF = minF;
parameters.maxF = maxF;
parameters.numModes = parameters.pcaModes;
numPerDataSet = parameters.trainingSetSize;
numPoints = 5000;
%%
marker=0;
for i = 103%:%length(dbase)
    marker = marker+1;
    if marker==4
        marker=0;z
        run_umap
    else
        for m = 1%:size(dbase(i).mtracks,4)
            tic

            fprintf(1,['Processing dist for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
            p1 = squeeze(dbase(i).mtracks(:,:,:,m));
            p1 = permute(p1, [2 3 1]);
            p1Dist = zeros(length(IDX),size(p1,3));
            for ii = 1:size(p1Dist,1)
                p1Dist(ii,:) = returnDist(squeeze(p1(X(ii),:,:)),squeeze(p1(Y(ii),:,:)));
            end
            p1Dsmooth = zeros(size(p1Dist));
            for ii = 1:size(p1Dist,1)
                p1Dsmooth(ii,:) = medfilt1(smooth(p1Dist(ii,:),'moving',5,'omitnan'),5, 'omitnan');
            end
            p1Dist = p1Dsmooth;

            fprintf(1,['Processing projections for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
            p2Dist = bsxfun(@minus,p1Dist,muv');
            projections_temp = p2Dist'*vecs(:,1:numProjections);
            for ii=1:pcaModes
                projections_temp(:,ii) = fillgaps(projections_temp(:,ii),20*fps);
            end

            [data,~] = findWavelets(projections_temp,numModes,parameters);
            nnData = log(data);
            nnData(nnData<-3) = -3;
            fprintf(1,['Processing umap transform for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
            % reembeding
            temp_umap = run_umap(nnData, 'template_file', 'Y:\behdata\2021-09-TSC1-social\SocialTsc1\mat_files\template_WT.mat');
            %temp_umap = run_umap(nnData, 'template_file', 'Y:\mkislin\DTSC\cohort_optoUR\template_WT.mat');

            umapClustWT = zeros(length(temp_umap),1);
            for ii =1:length(temp_umap)
                [~,idx_test] = pdist2(embedding,temp_umap(ii,:),'cityblock' ,'Smallest',1);
                umapClustWT(ii,1) = int16(clusterIdentifiers(idx_test));
                clear idx_test
            end
            dbase(i).umapClustWT(:,1:2,m) = temp_umap;
            dbase(i).umapClustWT(:,3,m) = umapClustWT;
            title(dbase(i).condition)
            sfname = strcat('Umap_',dbase(i).fileID,'_m',num2str(m),'.pdf');
            saveas(gcf,sfname)
            close(gcf)
            toc
            clear umapClustWT nnData data projections_temp p1 p1Dist p1Dsmooth p2Dist temp_umap
            % clear RAM
            % user = memory;
            % sprintf('%.0f MB',user.MemUsedMATLAB/(1024^2))
            % delete(gcp('nocreate'))
            % clear JAVA
            % user = memory;
            % sprintf('%.0f MB',user.MemUsedMATLAB/(1024^2))
        end
    end
end
%
save(Bugbasev2,'dbase');
%% Trainset Step2: sampling on the pre-existing data
% % Generate 20 equally-spaced points in 2D with pre-existing data
% tic
% %y = rand(10, 2);
% x = maximin(2000, 2, 'data', embedding(1:100:end,:),'cycles',1, 'iterations', 100);
% toc
% figure
% hold on
% plot(x(:, 1), x(:, 2), 'x')
% %plot(y(:, 1), y(:, 2), 'o')
% %plot(embedding(1:100:end,1), embedding(1:100:end,2), 'o')
%% Trainset Step2: Simple and faster sampling to get frames indexes in TrainingIdxRare
% Equal-size spectral clustering
yy = linspace(min(embedding(:,2)),max(embedding(:,2)),100);
xx = linspace(min(embedding(:,1)),max(embedding(:,1)),100);
[XX, YY] = meshgrid(xx,yy);

xx = embedding(1:100:end,1);
yy = embedding(1:100:end,2);
k = boundary(xx,yy);
TFin = inpolygon(XX(:), YY(:),xx(k), yy(k));
points = [XX(:), YY(:)];
%plot(points(TFin,1), points(TFin,2),'x')

temp_idx = find(TFin==1);
kk=0;
TrainingIdxRare = cell(length(dbase),2);
tic
for i = 1:length(dbase)
    for m = 1:size(dbase(i).umapClustWT,3)
        
        temp_idxRare=[];

        for ii=1:length(temp_idx)
            [~,idx_test] = pdist2(squeeze(dbase(i).umapClustWT(:,1:2,m)),points(temp_idx(ii,1),:),'cityblock' ,'Smallest',1);
            if dbase(i).excludedFrames(idx_test)==0
                temp_idxRare = [temp_idxRare; idx_test];
            else
            end
            clear idx_test
        end
        TrainingIdxRare{i,m}=temp_idxRare;
        kk=kk+length(temp_idxRare);
        
    end
end
toc

%% Trainset Step3: calculate distances and do online PCA to find modes

% here we select body parts for futher analysis
noi = [1,2,5,6,7,8,10,12];
xIdx = noi;
yIdx = noi;

%  prepare matrix to calculate distances from body part positions
[X, Y] = meshgrid(xIdx,yIdx);
X = X(:); Y = Y(:);
IDX = find(X~=Y);

batchSize = 6000;
firstBatch = true;
currentImage = 0;
TrainingDistances=[];
TrainingIdx=cell(length(dbase),2);

for i = 1:length(dbase)
    for m = 1:size(dbase(i).umapClustWT,3)

        rand_frames = TrainingIdxRare{i,m}; % frames based on the UmapWT
        % calculate distance between all body parts
        fprintf(1,['Processing dist for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);

        p1 = squeeze(dbase(i).tracks(:,:,:,m));
        p1 = permute(p1, [2 3 1]);
        p1Dist = zeros(length(IDX),size(p1,3));

        for ii = 1:size(p1Dist,1)
            p1Dist(ii,:) = returnDist(squeeze(p1(X(ii),:,:)),squeeze(p1(Y(ii),:,:)));
        end

        p1Dsmooth = zeros(size(p1Dist));
        for ii = 1:size(p1Dist,1)
            p1Dsmooth(ii,:) = medfilt1(smooth(p1Dist(ii,:),'moving',5,'omitnan'),5, 'omitnan');
        end

        p1Dist = p1Dsmooth(:,rand_frames); % only frames with distance between all body parts

        % for training
        % online PCA
        if firstBatch
            firstBatch = false;
            if size(p1Dist,2)<batchSize
                cBatchSize = size(p1Dist,2);
                XX = p1Dist';
            else
                cBatchSize = batchSize;
                XX = p1Dist(:,randperm(size(p1Dist,2),cBatchSize))';
            end

            currentImage = batchSize;
            muv = sum(XX);
            C = cov(XX).*batchSize + (muv'*muv)./cBatchSize;
        else

            if size(p1Dist,2)<batchSize
                cBatchSize = size(p1Dist,2);
                XX = p1Dist';
            else
                cBatchSize = batchSize;
                XX = p1Dist(:,randperm(size(p1Dist,2),cBatchSize))';
            end

            tempMu = sum(XX);
            muv = muv+tempMu;
            C = C + cov(XX).*cBatchSize + (tempMu'*tempMu)./cBatchSize;
            currentImage = currentImage + cBatchSize;
            fprintf(1,['Current Image = ' num2str(currentImage) '\n']);
        end

        TrainingDistances = [TrainingDistances p1Dist];
        TrainingIdx{i,m}=rand_frames;

        clear p1 p1Dist p1Dsmooth rand_frames

    end
end

save('tempC_stats_sLEAPSocialTsc1_Feb23.mat','C','muv','currentImage','TrainingDistances','TrainingIdx');

L = currentImage;
muv = muv./L;
C = C./L - muv'*muv;

[vecs,vals] = eig(C);
vals = flipud(diag(vals));
vecs = fliplr(vecs);
figure;
plot(1-(vals./sum(vals)))
prepfig;
xlabel('PCA component')
ylabel('var explained')

save('vecsVals_sLEAPSocialTsc1_Feb23.mat','vecs','vals','muv');

%% Trainset Step4: PCA projections (for each movie) and wavelet decomposition to generate training dataset
load('vecsVals_sLEAPSocialTsc1_Feb23.mat');
dataTrain = [];

numProjections = 10; % Note that in for social =12
numModes = 10; pcaModes = numModes;
minF = .25; maxF = 20;

parameters = setRunParameters([]);
parameters.trainingSetSize = 400;
parameters.pcaModes = pcaModes;
parameters.samplingFreq = fps;
parameters.minF = minF;
parameters.maxF = maxF;
parameters.numModes = parameters.pcaModes;
numPerDataSet = parameters.trainingSetSize;
numPoints = 5000;

for i = 1:length(dbase)
    for m =1:size(dbase(i).tracks,4)
        fprintf(1,['Processing dist for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
        p1 = squeeze(dbase(i).tracks(:,:,:,m));
        p1 = permute(p1, [2 3 1]);
        p1Dist = zeros(length(IDX),size(p1,3));
        for ii = 1:size(p1Dist,1)
            p1Dist(ii,:) = returnDist(squeeze(p1(X(ii),:,:)),squeeze(p1(Y(ii),:,:)));
        end
        p1Dsmooth = zeros(size(p1Dist));
        for ii = 1:size(p1Dist,1)
            p1Dsmooth(ii,:) = medfilt1(smooth(p1Dist(ii,:),'moving',5,'omitnan'),5, 'omitnan');
        end
        p1Dist = p1Dsmooth;

        fprintf(1,['Processing projections for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
        p2Dist = bsxfun(@minus,p1Dist,muv'); % PCA modes
        projections_temp = p2Dist'*vecs(:,1:numProjections);
        for ii=1:pcaModes
            projections_temp(:,ii) = fillgaps(projections_temp(:,ii),20*fps);
        end

        [data,~] = findWavelets(projections_temp,numModes,parameters);
        amps = sum(data,2);

        signalIdx = TrainingIdx{i,m};
        signalAmps = amps(signalIdx);
        nnData = log(data(signalIdx,:));
        nnData(nnData<-3) = -3;

        dataTrain = [dataTrain; nnData];
    end
end

%% Trainset Step4a: run uMap on the social training dataset
[uniqueA i j] = unique(dataTrain,'row','first');
indexToDupes = find(not(ismember(1:numel(dataTrain),i)));

[reduction, umapTsc1So, clusterIdentifiers, extras] = run_umap(uniqueA,...
    'min_dist' , 0.0001,'spread', 10, 'n_neighbors', 10,'n_components' ,2, 'metric', 'cityblock',...
    'cluster_output','numeric','cluster_detail','high',...
    'qf_tree' ,'true', 'save_template_file', 'template_Tsc1So_Feb23_n-neig10.mat');
%%
saveas(gcf,'Umap_Tsc1So_Feb23_n-neig10.pdf')
save('trainingser_sLEAP_SocialTsc1_8Feb23_n-neig10.mat', 'uniqueA','clusterIdentifiers','umapTsc1So','extras', 'reduction')
%% Trainset Step4b: clustering with DBScan            
%idx = dbscan(umapTsc1So.embedding,1,5);
%figure;
%gscatter(umapTsc1So.embedding(:,1),umapTsc1So.embedding(:,2),idx);

% D = pdist2(umapTsc1So.embedding(1:100:end,:),umapTsc1So.embedding(1:100:end,:));
% [idx, corepts] = dbscan(D,10,5,'Distance','precomputed');
% numGroups = length(unique(idx));
% figure;
% gscatter(umapTsc1So.embedding(1:100:end,1),umapTsc1So.embedding(1:100:end,2),idx,hsv(numGroups));

%% Instect Social Behavior embedding
load('Z:\mkislin\t-SNEmouse_opto-mini project\MouseMotionMapper\utilities\colormaps.mat')

figure;
boundaries =cell(13,1);
for clust=min(clusterIdentifiers)+1:max(clusterIdentifiers)
    temp_x =  umapTsc1So.embedding(clusterIdentifiers==clust,1);
    temp_y =  umapTsc1So.embedding(clusterIdentifiers==clust,2);
    kk = boundary(temp_x ,temp_y);
    boundaries{clust+1,1}(:,1) = temp_x(kk);
    boundaries{clust+1,1}(:,2)= temp_y(kk);
    hold on;
    plot(temp_x(kk),temp_y(kk));
    text(mean(temp_x(kk)),mean(temp_y(kk)),num2str(clust)) % add cluster IDs
end
xlim([-40 40]); ylim([-40 40]); prepfig;

figure; gscatter(umapTsc1So.embedding(:,1),umapTsc1So.embedding(:,2),clusterIdentifiers);
xlim([-40 40]); ylim([-40 40]); prepfig;
axis equal off xy

figure;
maxVal = max(max(abs(umapTsc1So.embedding)));
maxVal = round(maxVal * 1.2);
sigma = maxVal / 40;
numPoints = 501;
rangeVals = [-maxVal maxVal];

[xx,density] = findPointDensity(umapTsc1So.embedding,sigma,numPoints,rangeVals);

maxDensity = max(density(:));
imagesc(xx,xx,density)
xlim([-40 40]); ylim([-40 40]); prepfig;
caxis([0 maxDensity*0.5])
colormap(cmapLL);

figure;
%imagesc(xx,xx,density)
%colormap gray
contour(xx,xx,density,50)
hold on
scatter(umapTsc1So.embedding(1:50:end,1), umapTsc1So.embedding(1:50:end,2),3,[0 0.4470 0.7410],'filled',...
    'MarkerFaceAlpha',0.2, 'MarkerEdgeColor','none');
xlim([-40 40]); ylim([-40 40]); prepfig;
axis equal off xy

%% Check clustering with clustergram
a=zeros(max(clusterIdentifiers)+1,size(uniqueA,2));
for i=min(clusterIdentifiers):max(clusterIdentifiers)
a(i+1,:) = mean(uniqueA(clusterIdentifiers(1,:)==i,:));
end
% Clustergram
Y = pdist(a,'cityblock');
test=squareform(Y);
Z=linkage(Y,'average');
% colored clustergram
%plot(a','DisplayName','a')
cgo = clustergram(a,'Colormap',redbluecmap, 'ColumnLabelsRotate',0,'Cluster',1);
set(cgo, 'Linkage','complete','Dendrogram', 15);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reembedding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reembedding Step1: Load variables 

TSC1Sbase = 'PR-OFT-TSC1-social-dbase-01022023.mat';
load('vecsVals_sLEAPSocialTsc1_Feb23.mat');
%load(TSC1Sbase);
load('trainingser_sLEAP_SocialTsc1_8Feb23_n-neig10.mat', 'clusterIdentifiers')
load('trainingser_sLEAP_SocialTsc1_8Feb23_n-neig10.mat', 'umapTsc1So')

% Metadata
fps=80;
pixel_size = 1/1.97; % mm

numProjections = 12;
numModes = 12; pcaModes = numModes;
minF = .25; maxF = 20;

parameters = setRunParameters([]);
parameters.trainingSetSize = 400;
parameters.pcaModes = pcaModes;
parameters.samplingFreq = fps;
parameters.minF = minF;
parameters.maxF = maxF;
parameters.numModes = parameters.pcaModes;
numPerDataSet = parameters.trainingSetSize;
numPoints = 5000;

noi = [1,2,5,6,7,8,10,12];
xIdx = noi;
yIdx = noi;

%  prepare matrix to calculate distances from body part positions
[X Y] = meshgrid(xIdx,yIdx);
X = X(:); Y = Y(:);
IDX = find(X~=Y);

%% Reembedding Step2: reembeed in 2d umap

%run_umap; % make sure umap has all mex-files and run on the example data
embedding = umapTsc1So.embedding;

for i = 55:59%length(dbase)
    temp = h5read(TSC1Sbase, '/dbase/tracks', [1,i],[1,1]);
    temp_id = h5read(TSC1Sbase, '/dbase/fileID', [1,i],[1,1]);
    temp_condition = h5read(TSC1Sbase, '/dbase/condition', [1,i],[1,1]);
    dbase(i).condition=temp_condition;
    dbase(i).tracks = temp{1,1};
    dbase(i).fileID = temp_id{1,1};
    for m = 1:size(dbase(i).tracks,4)
        tic

        fprintf(1,['Processing dist for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
        p1 = squeeze(dbase(i).tracks(:,:,:,m));
        p1 = permute(p1, [2 3 1]);
        p1Dist = zeros(length(IDX),size(p1,3));
        for ii = 1:size(p1Dist,1)
            p1Dist(ii,:) = returnDist(squeeze(p1(X(ii),:,:)),squeeze(p1(Y(ii),:,:)));
        end
        p1Dsmooth = zeros(size(p1Dist));
        for ii = 1:size(p1Dist,1)
            p1Dsmooth(ii,:) = medfilt1(smooth(p1Dist(ii,:),'moving',5,'omitnan'),5, 'omitnan');
        end
        p1Dist = p1Dsmooth;

        fprintf(1,['Processing projections for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
        p2Dist = bsxfun(@minus,p1Dist,muv');
        projections_temp = p2Dist'*vecs(:,1:numProjections);
        for ii=1:pcaModes
            projections_temp(:,ii) = fillgaps(projections_temp(:,ii),20*fps);
        end

        [data,~] = findWavelets(projections_temp,numModes,parameters);
        nnData = log(data);
        nnData(nnData<-3) = -3;

        fprintf(1,['Processing umap transform for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
        % check the template
        temp_umap = run_umap(nnData, 'template_file', 'Y:\behdata\2021-09-TSC1-social\SocialTsc1\template_Tsc1So_Feb23_n-neig10.mat');
        for ii =1:length(temp_umap)
            [~,idx_test] = pdist2(embedding,temp_umap(ii,:),'cityblock' ,'Smallest',1);
            umapClust(ii,1) = int16(clusterIdentifiers(idx_test));
            clear idx_test
        end
        dbase(i).umapClust(:,1:2,m) = temp_umap;
        dbase(i).umapClust(:,3,m) = umapClust;
        title(dbase(i).condition)
        saveas(gcf,['Umap_',dbase(i).fileID,'_m',num2str(m),'.pdf'])
        close(gcf)
        toc
        clear umapClust nnData data projections_temp p1 p1Dist p1Dsmooth p2Dist temp_umap
        user = memory;
        sprintf('%.0f MB',user.MemUsedMATLAB/(1024^2))
        delete(gcp('nocreate'))
        clear JAVA
        user = memory;
        sprintf('%.0f MB',user.MemUsedMATLAB/(1024^2))
    end
    temp_umapClust = dbase(i).umapClust;
    save(['Umap_',num2str(i),'.mat'],  'temp_umapClust')
    clear temp_umapClust temp_id temp_condition 
end
%
%save(TSC1Sbase,'dbase');
%% load all umap mat to the dbase

for i=53:length(dbase)
    load(['Umap_',num2str(i),'.mat'],'temp_umapClust');
    dbase(i).umapClust = temp_umapClust;
    clear temp_umapClust
end

%% inspect density plots per mouse and average per condition
%maxVal = max(max(abs(UMAPembedding)));
% maxVal = 15; %round(maxVal * 1.2);
% sigma = maxVal / 40;
% numPoints = 1000;
% rangeVals = [-15 10];
load('Z:\mkislin\t-SNEmouse_opto-mini project\MouseMotionMapper\utilities\colormaps.mat')
fig_path = 'Y:\behdata\2021-09-TSC1-social\SocialTsc1\figures-Feb23\UmapDensityPerMouse\';

for i=53:length(dbase)
    for m=1:2
        UMAPembedding = squeeze(dbase(i).umapClust(:,1:2,m));
        maxVal = max(max(abs(UMAPembedding)));
        maxVal = round(maxVal * 1.2);
        sigma = maxVal / 40;
        numPoints = 501;
        rangeVals = [-maxVal maxVal];
        [xx,density] = findPointDensity(UMAPembedding,sigma,numPoints,rangeVals);
        figure;
        maxDensity = max(density(:));
        imagesc(xx,xx,density)
        xlim([-40 40]); ylim([-40 40]); prepfig;
        caxis([0 maxDensity*0.5])
        colormap(cmapLL);
        set(gca,'YDir','normal');
        title(dbase(i).condition)
        saveas(gcf,[fig_path,'Umap_',dbase(i).fileID,'_m',num2str(m),'.pdf'])
        close(gcf)
        
        dbase(i).UMapDensity(:,:,m) = density;
    end
end

%% inspect speed in umap and tSNE space with ethograms

% calculate density map for umap
maxVal = max(max(abs(embedding)));
maxVal = round(maxVal * 1.2);
sigma = maxVal / 40;
numPoints = 501;
rangeVals = [-maxVal maxVal];
[xx,density] = findPointDensity(embedding,sigma,numPoints,rangeVals);
maxDensity = max(density(:));
% calculate density map for tSNE
% maxVal = max(max(abs(ydata)));
% maxVal = round(maxVal * 1.2);
% sigma = maxVal / 40;
% numPoints = 501;
% rangeVals = [-maxVal maxVal];
% [xxSNE,densitySNE] = findPointDensity(ydata,sigma,numPoints,rangeVals);
% maxDensitySNE = max(densitySNE(:));

jpl=[1 5 6 7 8 10 12]; % index of joints to plot
% colors for body part labels
color1 = [221,28,119]./255;
color2 = [142,205,212]./255;
colors_behaviors = [1,1,1; 254,218,117; 250,126,30; 214,41,118; 102,194,165;150,47,191;252,141,98;141,160,203;231,138,195;166,216,84;179,179,179;120,120,120;80,80,80]./255;
fig_path = 'Y:\behdata\2021-09-TSC1-social\Predictions_RT\';
minframes = 40;

for f=1:52%length(dbase) %video
    clust_all = squeeze(dbase(f).umapClust(:,3,:));% for pclusts
    clust_select = largeBWConnComp((((clust_all(:,1)== 27) + (clust_all(:,2)==27))==2), minframes); % find epochs
    for s=1:length(clust_select.PixelIdxList)
        length_stretch = length(clust_select.PixelIdxList{1,s});
        if length_stretch>10

            % EDIT, TO SELECT FOR EACH INSTANCE GREATER THAN 40 FRAMES
            temp_str = min(clust_select.PixelIdxList{1,s}); % 1; %start frame
            temp_stp = max(clust_select.PixelIdxList{1,s}); % 4000; %end frame

            beh_vid = h5read(dbase(f).meta.behvideo,'/pg0', [1 1 1 temp_str], [Inf, Inf, Inf, temp_stp+1-temp_str],[1,1,1,1]);
            temp_tracks1 = squeeze(dbase(f).tracks(temp_str:temp_stp,:,:,1));
            temp_tt1 = squeeze(temp_tracks1(:,10,:));
            temp_tracks2 = squeeze(dbase(f).tracks(temp_str:temp_stp,:,:,2));
            temp_tt2 = squeeze(temp_tracks2(:,10,:));
            temp_time = dbase(f).time(temp_str:temp_stp,1);

            temp_umap1 = dbase(f).umapClust(temp_str:temp_stp,:,1);
            temp_umap2 = dbase(f).umapClust(temp_str:temp_stp,:,2);

%             temp_tSNE1 = squeeze(dbase(f).eV(temp_str:temp_stp,:,1));
%             temp_tSNE2 = squeeze(dbase(f).eV(temp_str:temp_stp,:,2));
% 
%             temp_pClusts1 = squeeze(dbase(f).pClusts(temp_str:temp_stp,1));
%             temp_pClusts2 = squeeze(dbase(f).pClusts(temp_str:temp_stp,2));

            avi_filename =['file-', dbase(f).fileID,'-',dbase(f).condition,'-frames-', num2str(temp_str),'-', num2str(temp_stp),'.avi'];
            avi_file = VideoWriter(avi_filename,'MPEG-4'); %'MPEG-4'
            avi_file.FrameRate= fps/8;
            avi_file.open();
            figure('Position',[100 50 1200 700])
            thm = tiledlayout(2,3,'TileSpacing','Compact','Padding','compact');

            for m=10:size(beh_vid,4)

                ax = nexttile;
                imshow(beh_vid(:,:,1,m))
                colormap(ax, "gray")
                hold on
                sz = linspace(1,5,10);
                for j=1:size(jpl,2)
                    %plot(tracks(m,j,1),tracks(m,j,2),'marker','s', 'MarkerSize', 2, 'MarkerEdgeColor', map(i,:))
                    scatter(temp_tracks1(m-9:m,jpl(j),1),temp_tracks1(m-9:m,jpl(j),2),sz,'filled','MarkerEdgeColor',color1,'MarkerFaceColor',color1)
                    hold on
                    scatter(temp_tracks2(m-9:m,jpl(j),1),temp_tracks2(m-9:m,jpl(j),2),sz,'filled','MarkerEdgeColor',color2,'MarkerFaceColor',color2)
                end

                axis equal off tight;
                %xlim([temp_tt1(m,1)-200, temp_tt1(m,1)+200])
                %ylim([temp_tt1(m,2)-200, temp_tt1(m,2)+200])
                title([dbase(f).condition, '-mouse-1']);

                % umap
                ax1 = nexttile;
                imagesc(xx,xx,density)
                caxis([0 maxDensity*0.5])
                colormap(ax1, flipud(pink))
                hold on
                scatter(temp_umap1(m-9:m,1),temp_umap1(m-9:m,2),sz,'filled','MarkerEdgeColor',color1,'MarkerFaceColor',color1)
                scatter(temp_umap2(m-9:m,1),temp_umap2(m-9:m,2),sz,'filled','MarkerEdgeColor',color2,'MarkerFaceColor',color2)
                axis equal off xy
                title('UMAP');

%                 ax2 =nexttile;
%                 sz = linspace(2,20,10);
%                 imagesc(xxSNE,xxSNE,densitySNE)
% 
%                 colormap(ax2, cmap1)
%                 caxis([0 maxDensitySNE*1.2])
%                 hold on
%                 scatter(temp_tSNE1(m-9:m,1),temp_tSNE1(m-9:m,2),sz,'filled','MarkerEdgeColor',color1,'MarkerFaceColor',color1)
%                 scatter(temp_tSNE2(m-9:m,1),temp_tSNE2(m-9:m,2),sz,'filled','MarkerEdgeColor',color2,'MarkerFaceColor',color2)
%                 axis equal off xy
%                 title('tSNE');

                ax0 = nexttile;
                imshow(beh_vid(:,:,1,m))
                colormap(ax0, "gray")
                hold on
                sz = linspace(1,5,10);
                for j=1:size(jpl,2)
                    %plot(tracks(m,j,1),tracks(m,j,2),'marker','s', 'MarkerSize', 2, 'MarkerEdgeColor', map(i,:))
                    scatter(temp_tracks1(m-9:m,jpl(j),1),temp_tracks1(m-9:m,jpl(j),2),sz,'filled','MarkerEdgeColor',color1,'MarkerFaceColor',color1)
                    hold on
                    scatter(temp_tracks2(m-9:m,jpl(j),1),temp_tracks2(m-9:m,jpl(j),2),sz,'filled','MarkerEdgeColor',color2,'MarkerFaceColor',color2)
                end

                axis equal off tight;
                %xlim([temp_tt2(m,1)-200, temp_tt2(m,1)+200])
                %ylim([temp_tt2(m,2)-200, temp_tt2(m,2)+200])
                title([dbase(f).condition, '-mouse-2']);

                ax3 = nexttile;

                imagesc([temp_umap1(:,3) temp_umap2(:,3)]')
                colormap(ax3, colors_behaviors)
                set(gca, 'PlotBoxAspectRatio', [1,0.5,1])
                xlim([m-2*fps+9, m+2*fps-9])
                hold on
                xline(ax3, [m-9, m-9], 'r--', 'LineWidth',2);
                box off
                set(gca,'TickDir','out')
                yticks([1,2])
                yticklabels({temp_umap1(m-9,3),temp_umap2(m-9,3)})
                xlabel('Frames')
                title('Behavior class');

%                 ax4 =nexttile;
%                 imagesc([temp_pClusts1 temp_pClusts2]')
%                 colormap(ax4, "lines")
% 
%                 set(gca, 'PlotBoxAspectRatio', [1,0.5,1])
%                 xlim([m-2*fps-9, m+2*fps-9])
%                 hold on
%                 xline(ax4, [m-9, m-9], 'r--', 'LineWidth',2);
%                 box off
%                 set(gca,'TickDir','out')
%                 yticks([1,2])
%                 yticklabels({temp_pClusts1(m-9,1),temp_pClusts2(m-9,1)})
%                 xlabel('Frames')
%                 title('Behavior clusters');

                %temp = getframe(gcf);
                %temp = imresize(temp.cdata, [700, 1200]);
                %avi_file.writeVideo(temp);
                avi_file.writeVideo(getframe(gcf));
                delete(nexttile(1));delete(nexttile(2));delete(nexttile(3));
                delete(nexttile(4));delete(nexttile(5));delete(nexttile(6));
                %clf(gcf,'reset')
                %clf('reset')

                %close all
            end
            clear m j temp temp_tt1 temp_tracks1 temp_tt2 temp_tracks2 temp_stp temp_str f1...
                temp_time temp_pClusts1 temp_pClusts2 temp_tSNE1 temp_tSNE2 temp_umap1 temp_umap2
            close all
            avi_file.close();
        end
    end
end


%% resetting
delete(nexttile(1))
delete(nexttile(2))
delete(nexttile(3))
delete(nexttile(4))
delete(nexttile(5))
delete(nexttile(6))


%% wavelet fingerprints
temp_amp=[];
avg_amp = cell(length(dbase),2);
for c = 8%:8
    for i=5%:length(dbase)
        for m =1:2
            for n=1:13
                temp_velocity=getVelocity(squeeze(dbase(i).tracks(:,n,:,m)));
                idx_cl=(dbase(i).umapClust(:,3,m)==c);
                idx_cl(dbase(i).badframes(:,m),1)=0;
                stretch =bwconncomp(idx_cl); %finding nan stretches
                for s=1:length(stretch.PixelIdxList)
                    %length_stretch(s,1) = length(stretch.PixelIdxList{1,s});
                    length_stretch = length(stretch.PixelIdxList{1,s});
                    %finding frame numbers
                    frames = stretch.PixelIdxList{1,s};
                    %transpose
                    frames=frames.';
                    if length_stretch>80
                        [temp_velocity_amp,~] = findWavelets(temp_velocity(frames,1),numModes,parameters);
                        temp_amp=[temp_amp temp_velocity_amp'];
                    else
                    end


                end
                avg_amp{i,m}(n,:)=nanmean(temp_amp,2);
            end


        end
    end
    figure;
    imagesc(avg_amp{i, m}(noi,:))
    caxis([0 0.2])
    prepfig;
    set(gca, 'PlotBoxAspectRatio', [1,0.5,1])
    drawnow()
end

%% Step7: Generate bradies

%% Step8: from clusters to classes



%% Save ethograms as database

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD CODE BELOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step1: learn a behavior space from animal pose

%% Step1a: calculate distances and do online PCA to find modes

% here we select body parts for futher analysis
noi = [1,2,5,6,7,8,10,12];
xIdx = noi;
yIdx = noi;

%  prepare matrix to calculate distances from body part positions
[X Y] = meshgrid(xIdx,yIdx);
X = X(:); Y = Y(:);
IDX = find(X~=Y);

% generate random sample of 50,000 frames

TrainingDistances=[];
TrainingIdx=cell(length(dbase),2);

batchSize = 2000;
firstBatch = true;
currentImage = 0;

for i = 1:length(dbase)
    good = find(dbase(i).excludedFrames(:) == 0);
    for m =1:size(dbase(i).tracks,4)

        %firstBatch = true;
        %good = find(dbase(i).excludedFrames(:) == 0); %array with all good frames
        rand_frames = randsample(good, batchSize-1); %1500 random frames per vid
        % calculate distance between all body parts
        fprintf(1,['Processing dist for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);

        p1 = squeeze(dbase(i).tracks(:,:,:,m));
        p1 = permute(p1, [2 3 1]);
        p1Dist = zeros(length(IDX),size(p1,3));

        for ii = 1:size(p1Dist,1)
            p1Dist(ii,:) = returnDist(squeeze(p1(X(ii),:,:)),squeeze(p1(Y(ii),:,:)));
        end

        p1Dsmooth = zeros(size(p1Dist));
        for ii = 1:size(p1Dist,1)
            p1Dsmooth(ii,:) = medfilt1(smooth(p1Dist(ii,:),'moving',5,'omitnan'),5, 'omitnan');
        end

        % find missing data
        [idx_nan_r,idx_nan_c]  = find(isnan(p1Dsmooth));
        rand_frames(ismember(rand_frames,unique(idx_nan_c)))=[];

        p1Dist = p1Dsmooth(:,rand_frames); % only frames with distance between all body parts

        % for training
        % online PCA
        if firstBatch
            firstBatch = false;
            if size(p1Dist,2)<batchSize
                cBatchSize = size(p1Dist,2);
                XX = p1Dist';
            else
                cBatchSize = batchSize;
                XX = p1Dist(:,randperm(size(p1Dist,2),cBatchSize))';
            end

            currentImage = batchSize;
            muv = sum(XX);
            C = cov(XX).*batchSize + (muv'*muv)./cBatchSize;
        else

            if size(p1Dist,2)<batchSize
                cBatchSize = size(p1Dist,2);
                XX = p1Dist';
            else
                cBatchSize = batchSize;
                XX = p1Dist(:,randperm(size(p1Dist,2),cBatchSize))';
            end

            tempMu = sum(XX);
            muv = muv+tempMu;
            C = C + cov(XX).*cBatchSize + (tempMu'*tempMu)./cBatchSize;
            currentImage = currentImage + cBatchSize;
            fprintf(1,['Current Image = ' num2str(currentImage) '\n']);
        end

        TrainingDistances = [TrainingDistances p1Dist ];
        TrainingIdx{i,m}=rand_frames;

        clear p1 p1Dist p1Dsmooth rand_frames

    end

end

save('tempC_stats_sLEAPSocialTsc1_Jan23withopto.mat','C','muv','currentImage','TrainingDistances','TrainingIdx');

L = currentImage;
muv = muv./L;
C = C./L - muv'*muv;

[vecs,vals] = eig(C);
vals = flipud(diag(vals));
vecs = fliplr(vecs);
figure;
plot(1-(vals./sum(vals)))
prepfig;
xlabel('PCA component')
ylabel('var explained')

save('vecsVals_sLEAPSocialTsc1_Jan23withopto.mat','vecs','vals','muv');

%% Step1b: PCA projections (for each movie) and wavelet decomposition to generate training dataset
load('vecsVals_sLEAPSocialTsc1_Jan23withopto.mat');
dataTrain = [];

numProjections = 10;
numModes = 10; pcaModes = numModes;
minF = .25; maxF = 20;

parameters = setRunParameters([]);
parameters.trainingSetSize = 400;
parameters.pcaModes = pcaModes;
parameters.samplingFreq = fps;
parameters.minF = minF;
parameters.maxF = maxF;
parameters.numModes = parameters.pcaModes;
numPerDataSet = parameters.trainingSetSize;
numPoints = 5000;

for i = 1:length(dbase)
    for m =1:size(dbase(i).tracks,4)
        fprintf(1,['Processing dist for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
        p1 = squeeze(dbase(i).tracks(:,:,:,m));
        p1 = permute(p1, [2 3 1]);
        p1Dist = zeros(length(IDX),size(p1,3));
        for ii = 1:size(p1Dist,1)
            p1Dist(ii,:) = returnDist(squeeze(p1(X(ii),:,:)),squeeze(p1(Y(ii),:,:)));
        end
        p1Dsmooth = zeros(size(p1Dist));
        for ii = 1:size(p1Dist,1)
            p1Dsmooth(ii,:) = medfilt1(smooth(p1Dist(ii,:),'moving',5,'omitnan'),5, 'omitnan');
        end
        p1Dist = p1Dsmooth;

        fprintf(1,['Processing projections for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
        p2Dist = bsxfun(@minus,p1Dist,muv'); % PCA modes
        projections_temp = p2Dist'*vecs(:,1:numProjections);
        for ii=1:pcaModes
            projections_temp(:,ii) = fillgaps(projections_temp(:,ii),20*fps);
        end

        [data,~] = findWavelets(projections_temp,numModes,parameters);
        amps = sum(data,2);

        signalIdx = TrainingIdx{i,m};
        signalAmps = amps(signalIdx);
        nnData = log(data(signalIdx,:));
        nnData(nnData<-3) = -3;

        dataTrain = [dataTrain; nnData];
    end
end

%% Step1c: clustering

% % %% k-means
% % 
% % tic
% % C100 = kmeans(dataTrain,100,'Replicates',20,'Distance','cityblock','MaxIter',100);
% % toc
% % %% tSNE
% % tic
% % ydata = tsne(dataTrain);
% % toc
% % 
% % maxVal = max(max(abs(ydata)));
% % maxVal = round(maxVal * 1.2);
% % sigma = maxVal / 40;
% % numPoints = 501;
% % rangeVals = [-maxVal maxVal];
% % 
% % [xx,density] = findPointDensity(ydata,sigma,numPoints,rangeVals);
% % 
% % maxDensity = max(density(:));
% % 
% % figure;
% % imagesc(xx,xx,density)
% % axis equal off xy
% % 
% % load('Z:\mkislin\t-SNEmouse_opto-mini project\MouseMotionMapper\utilities\colormaps.mat')
% % colormap(cmap1)
% % caxis([0 maxDensity * 1.2])
% % prepfig

%% find noncommon wavelets outside the initial umap clusters

% % dataTrainRare=[];
% % for i = 1:length(dbase)
% %     for m = 1:size(dbase(i).umapClustWT,3)
% %         tic
% % 
% %         fprintf(1,['Processing dist for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
% %         p1 = squeeze(dbase(i).tracks(:,:,:,m));
% %         p1 = permute(p1, [2 3 1]);
% %         p1Dist = zeros(length(IDX),size(p1,3));
% %         for ii = 1:size(p1Dist,1)
% %             p1Dist(ii,:) = returnDist(squeeze(p1(X(ii),:,:)),squeeze(p1(Y(ii),:,:)));
% %         end
% %         p1Dsmooth = zeros(size(p1Dist));
% %         for ii = 1:size(p1Dist,1)
% %             p1Dsmooth(ii,:) = medfilt1(smooth(p1Dist(ii,:),'moving',5,'omitnan'),5, 'omitnan');
% %         end
% %         p1Dist = p1Dsmooth;
% % 
% %         fprintf(1,['Processing projections for mouse ' num2str(m) ' file #' num2str(dbase(i).fileID) '\n']);
% %         p2Dist = bsxfun(@minus,p1Dist,muv');
% %         projections_temp = p2Dist'*vecs(:,1:numProjections);
% %         for ii=1:pcaModes
% %             projections_temp(:,ii) = fillgaps(projections_temp(:,ii),20*fps);
% %         end
% % 
% %         [data,~] = findWavelets(projections_temp,numModes,parameters);
% %         nnData = log(data);
% %         nnData(nnData<-3) = -3;
% %           
% %         idx_rare = TrainingIdxRare{i,m};
% % 
% %         dataTrainRare = [dataTrainRare; nnData(idx_rare,:)];
% %         TrainingIdxRare{i,m}=idx_rare;
% %         close(gcf)
% %         toc
% %         clear temp_umap in2 in4 in5 idx_rare umapClust nnData data projections_temp p1 p1Dist p1Dsmooth p2Dist temp_umap...
% %             zValues zCosts zGuesses inConvHull meanMax exitFlags outputStatistics eV
% % 
% %     end
% % end

%% final uMap
[reduction, umapTsc1So, clusterIdentifiers, extras] = run_umap(dataTrainCombo,...
    'min_dist' , 0.1, 'n_neighbors', 15,'n_components' ,2, 'metric', 'cityblock',...
    'cluster_output','numeric', 'cluster_detail','medium', 'save_template_file', 'template_Tsc1So_Jan23Combo.mat');

% 'min_dist' , 0.1, 'n_neighbors', 15,'n_components' ,2, 'metric', 'cityblock',...
%    'cluster_output','numeric', 'cluster_detail','medium', 'save_template_file', 'template_Tsc1So_Jan23Combo.mat');
%% save training results
save('trainingser_sLEAP_SocialTsc1_Feb21.mat','dataTrain','ydata','C100','umapTsc1So.embedding','clusterIdentifiers');


%% create local mex for pdist2
% codegen -config cfg findNearestCentroid -args {dataTrain,nnData(ii,:)}
