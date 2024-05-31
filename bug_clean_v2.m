%% CREAT DATABASE as a matlab structure

%d = uigetdir([],'Pick a root directory');
d = 'Y:\behdata\2023-06-SocDev\ELE-analysis';
path_vid = 'Y:\behdata\2023-06-SocDev\C57';
files = dir(fullfile(d,'**\*.*'));
files = arrayfun(@(x) fullfile(x.folder,x.name),files,'UniformOutput',false);
files = files(endsWith(files,'.h5'));
%Name of file to hold the database
Bugbasev2 = 'OFT-social-wbug-dbasev2.mat';

%read csv file 
parameters = readmatrix([d,'\Video List v2.csv'], OutputType="string");

%% 

%Draw each session in directory into database 
for i=104%:size(parameters,1)
    fileID = parameters(i, 1);
    camera = parameters(i, 2);
    week = parameters(i, 3);
    sex = parameters(i, 4);
    condition = parameters(i, 5);
    mouseID1 = parameters(i, 6);
    mouseID2 = parameters(i, 7);
    
    %Only selecting the text files for analysis
    if exist(Bugbasev2,'file')
        load(Bugbasev2)
        newIdx = length(dbase)+1;
        %Check if file already loaded and plot state
        currentFIDs = {dbase.fileID};
    else
        dbase = struct;
        newIdx = 1;
        currentFIDs = [];
    end

    %mice sleap files 
    sleapmice_nam = strcat(d,'\predictions-mice-all\OFTsocialgroup-0', fileID, '-90-mice.slp.h5');
    mtracks = h5read(sleapmice_nam, '/tracks');
    mnode_names = h5read(sleapmice_nam, '/node_names');
    moccupancy = h5read(sleapmice_nam, '/track_occupancy');
    mpoint_scores = h5read(sleapmice_nam,'/point_scores');
    mtracking_scores = h5read(sleapmice_nam,'/tracking_scores');

    metam = h5info(sleapmice_nam);

    fprintf(1, 'File: %s \n',sleapmice_nam);
    fprintf(1, 'Number of tracks: %d \n',size(mtracks,4))
    fprintf(1, 'Number of body parts tracked: %d \n',size(mnode_names,1))
    fprintf(1, 'With scores for predictions: %d \n',size(metam.Datasets,1)>4)

    for ii=1:size(mtracks,4)
        fprintf(1, 'Track #%d - frames with predictions %d / %d \n' ,ii ,sum(moccupancy(ii,:)),size(mtracks,1));
        %for n=1:length(noi)
            %fprintf(1, '#%s - frames with NaNs %d / %d \n' ,mnode_names(noi(n),1) ,sum(isnan(mtracks(:,noi(n),1,ii))),size(mtracks,1));
        %end
    end
    time_frames = h5read(strcat(path_vid, sex,'ale\week', week, '\OFTsocialgroup-0',fileID,'-90.h5'), '/pg0_time');% '-00.h5' for WTWT
    time_frames = time_frames';

  
    %bug sleap file
    sleapbug_nam = strcat(d,'\predictions-bug-all\OFTsocialgroup-0', fileID, '-90-bug.h5');
    btracks = h5read(sleapbug_nam, '/tracks');
    bnode_names = h5read(sleapbug_nam, '/node_names');
    boccupancy = h5read(sleapbug_nam, '/track_occupancy');
    bpoint_scores = h5read(sleapbug_nam,'/point_scores');
    btracking_scores = h5read(sleapbug_nam,'/tracking_scores');

    metab = h5info(sleapbug_nam);

    fprintf(1, 'File: %s \n',sleapbug_nam);
    fprintf(1, 'Number of tracks: %d \n',size(btracks,4))
    fprintf(1, 'Number of body parts tracked: %d \n',size(bnode_names,1))
    fprintf(1, 'With scores for predictions: %d \n',size(metab.Datasets,1)>4)

    for ii=1:size(btracks,4)
        fprintf(1, 'Track #%d - frames with predictions %d / %d \n' ,ii ,sum(boccupancy(ii,:)),size(btracks,1));
        %for n=1:length(noi)
            %fprintf(1, '#%s - frames with NaNs %d / %d \n' ,bnode_names(noi(n),1) ,sum(isnan(btracks(:,noi(n),1,ii))),size(mtracks,1));
        %end
    end

    if any(strcmp(fileID,currentFIDs))
        return;
    else
        %Pack in new dbase element and save
        dbase(newIdx).fileID = fileID;
        %dbase(newIdx).fname = fname;

        dbase(newIdx).condition = condition;
        dbase(newIdx).week = week;
        dbase(newIdx).camera = camera;
        dbase(newIdx).sex = sex;
        dbase(newIdx).mouseID1 = mouseID1;
        dbase(newIdx).mouseID2 = mouseID2;
        dbase(newIdx).SLEAPtracksm = mtracks;
        dbase(newIdx).SLEAPscoresm = mpoint_scores;
        dbase(newIdx).time = time_frames;
        dbase(newIdx).metam.behvideo = strcat(path_vid, sex,'ale\week', week, '\OFTsocialgroup-0',fileID,'-90.h5'); %[path_vid,'OFTsocialgroup-0',fileID,'-90.h5']; %'-00.h5' for WTWT
        dbase(newIdx).metam.node_names = mnode_names;
        dbase(newIdx).metam.tracking_scores = mtracking_scores;
        dbase(newIdx).metam.occupancy = moccupancy;

        dbase(newIdx).SLEAPtracksb = btracks;
        dbase(newIdx).SLEAPscoresb = bpoint_scores;
        dbase(newIdx).metab.node_names = bnode_names;
        dbase(newIdx).metab.tracking_scores = btracking_scores;
        dbase(newIdx).metab.occupancy = boccupancy;


        save(Bugbasev2,'dbase');
    end
    clear fileID file fpath condition time_frames trials

end

%% VERIFY DATABASE and quality of SLEAP predictions
load('OFT-social-wbug-dbasev2.mat')

%% constants
fps=80;

pixel_size = 1/1.97; % mm
noi = [1,2,5,6,7,8,10,12]; % index body parts of interest
threshold_quantile=0.9;
% colors for paws
colors_paws = [205,32,41;61,80,157;166,75,156;142,205,212]./255; %  RF,LF,RH,LH

%% bug
% creating filtered tracks and finding frames with good predicitons for traning

for i=103:length(dbase)
    dbase(i).tracks = zeros (size(dbase(i).SLEAPtracksb));
    %for m = 1:size(dbase(i).SLEAPtracksb,4) %mouse
        for p=1:size(dbase(i).SLEAPtracksb,2)%noi %body part

            for c = 1:size(dbase(i).SLEAPtracksb,3) %coordinates
                x_raw = dbase(i).SLEAPtracksb(:,p,c); % get one position over time
                sleap_scores = dbase(i).SLEAPscoresb;
                x_raw(sleap_scores(:, c) < 0.4) = NaN;
                
                cutoff=quantile(abs(diff(x_raw)),threshold_quantile); %  quantile threshold for between-frame position displacement
                [x_new]=repeat_filter_pos(x_raw,dbase(i).time,threshold_quantile, 500); % The outliers are replaced by linear interpolations
                dbase(i).tracks(:,p,c) = x_new;
                nans = isnan(x_new); %finding nans
                dbase(i).metab.frNoSLEEAP(p) =  sum(nans)/length(nans); % fraction of nans per position
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
                fm = fillmissing(dbase(i).SLEAPtracksb(:,p,c),'linear');
                dbase(i).tracks(frames_idx, p, c) = fm(frames_idx,1);
                clear frames_idx idx_temp fm

                idx_temp  = find(length_stretch>3 & length_stretch<5*fps);
                frames_idx = cell2mat(stretch.PixelIdxList(1,idx_temp)');
                fm1 = fillgaps(dbase(i).SLEAPtracksb(:,p, c), 20*fps);
                dbase(i).tracks(frames_idx, p,c) = fm1(frames_idx,1);

                clear frames_idx idx_temp fm1 x_new x_raw length_stretch

                x_raw = dbase(i).tracks(:, p,c); % get one position over time

                cutoff=quantile(abs(diff(x_raw)),threshold_quantile); %  quantile threshold for between-frame position displacement
                if cutoff > 200
                    cutoff=15.0;
                end
                [x_new]=repeat_filter_pos(x_raw,dbase(i).time,threshold_quantile, 500); % The outliers are replaced by linear interpolations
                x_new(x_new>1000 | x_new<200, 1)=NaN;
                dbase(i).tracks(:,p,c) = x_new;
                nans = isnan(x_new); %finding nans
                dbase(i).metab.frNoSLEEAP(p) =  sum(nans)/length(nans); % fraction of nans per position
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
                fm = fillmissing(dbase(i).tracks(:,p,c),'linear');
                dbase(i).tracks(frames_idx, p, c) = fm(frames_idx,1);
                clear frames_idx idx_temp fm.

                %3 frames TO 1 sec - reconstruction with fillgaps local
                idx_temp  = find(length_stretch>3 & length_stretch<fps);
                frames_idx = cell2mat(stretch.PixelIdxList(1,idx_temp)');
                fm1 = fillgaps(dbase(i).SLEAPtracksb(:,p,c),3*fps);
                dbase(i).tracks(frames_idx, p, c) = fm1(frames_idx,1);

                clear frames_idx idx_temp fm1 x_new x_raw
            end
                
        end
        clear frames_idx nans

 end
%end

clear p s i c m

%% mouse
% creating filtered tracks and finding frames with good predicitons for traning

for i=103:length(dbase)
    dbase(i).mtracks = zeros (size(dbase(i).SLEAPtracksm));
    for m = 1:size(dbase(i).SLEAPtracksm,4) %mouse
        for p=1:size(dbase(i).SLEAPtracksm,2)%noi %body part

            for c = 1:size(dbase(i).SLEAPtracksm,3) %coordinates
                x_raw = dbase(i).SLEAPtracksm(:,p,c,m); % get one position over time
                sleap_scores = dbase(i).SLEAPscoresm;
                x_raw(sleap_scores(:, c) < 0.4) = NaN;
                nans = isnan(x_raw); %finding nans
                cutoff=quantile(abs(diff(x_raw)),threshold_quantile); %  quantile threshold for between-frame position displacement
                % [x_new]=filter_position(x_raw,dbase(i).time,cutoff); % The outliers are replaced by linear interpolations
                [x_new]=repeat_filter_pos(x_raw,dbase(i).time, threshold_quantile, 100);
                dbase(i).mtracks(:,p,c,m) = x_new;
                dbase(i).metam.frNoSLEEAP(p,m) =  sum(nans)/length(nans); % fraction of nans per position
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
                fm = fillmissing(dbase(i).SLEAPtracksm(:,p,c,m),'linear');
                dbase(i).mtracks(frames_idx, p, c, m) = fm(frames_idx,1);
                clear frames_idx idx_temp fm

                % %3 frames TO 1 sec - reconstruction with fillgaps local
                % idx_temp  = find(length_stretch>3 & length_stretch<fps);
                % frames_idx = cell2mat(stretch.PixelIdxList(1,idx_temp)');
                % fm1 = fillgaps(dbase(i).SLEAPtracksm(:,p,c,m),3*fps);
                % dbase(i).mtracks(frames_idx, p, c, m) = fm1(frames_idx,1);
                % 
                clear frames_idx idx_temp fm1 x_new x_raw
            end
        end
        clear frames_idx nans

 
        dbase(i).badframesm(:,m) = any(isnan(squeeze(dbase(i).mtracks(:, noi,2,m))),2);
        %define box, make nan for outside 
        %dbase(i).badframes(:,m) = any(isnan(squeeze(dbase(i).tracks(:, noi(1:end-1),2,m))),2);% without tail tip
        dbase(i).BadFract(1,m)=sum(dbase(i).badframesm(:,m))/length(dbase(i).badframesm);

    end
end

clear p s i c m