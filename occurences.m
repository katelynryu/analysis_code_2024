
d = 'Y:\behdata\2023-06-SocDev\bug\code';
csvFiles = dir(fullfile(d, '*.csv'));

fullT = table();

combinedT = table(zeros(13, 1), ...
    [0; 1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12], zeros(13, 1), ...
    zeros(13, 1), string.empty(13, 0), string.empty(13, 0), ...
    'VariableNames', {'Mouse', 'State', 'Count', 'Fraction', 'Group', 'FileID'});

for i = 1:numel(csvFiles)
    filename = fullfile(d, csvFiles(i).name);
    fileID = csvFiles(i).name(1:3);
    tempT = readtable(filename);

    lastCol = tempT{:, end};
    valueCounts = tabulate(lastCol); 
    tempFraction = zeros(13, 1);

    for j=1:13
        idx = find(valueCounts(:, 1) == combinedT.State(j));

        if ~isempty(idx)
            combinedT.Count(j) = valueCounts(idx, 2);
            tempFraction(j) = valueCounts(idx, 3) / 100;
        end
        
        combinedT.Mouse(j) = i;

        for f = 1:length(dbase)
            if(strcmp(dbase(f).fileID, fileID))
                idx2 = f;
            end
        end

        combinedT.Group(j) = dbase(idx2).condition;
        combinedT.FileID(j) = fileID;
    end

    if  (any(tempFraction == 0))
        tempidx = find(tempFraction == 0);
        tempFraction(tempidx) = realmin;
        [~, maxidx] = max(tempFraction);
        tempFraction(maxidx) = tempFraction(maxidx) - length(tempidx) * realmin; 
    end
    
    combinedT.Fraction = tempFraction;

    fullT = [fullT; combinedT];
end


%%
writetable(fullT, "umap-combine-wcontrol.csv");


%% for control data

d = 'Y:\behdata\2023-06-SocDev\bug\code\F-WT-control';
csvFiles = dir(fullfile(d, '*.csv'));

combinedT = table(zeros(13, 1), ...
    [0; 1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12], zeros(13, 1), ...
    zeros(13, 1), string.empty(13, 0), string.empty(13, 0), ...
    'VariableNames', {'Mouse', 'State', 'Count', 'Fraction', 'Group', 'FileID'});

for i = 1:numel(csvFiles)
    filename = fullfile(d, csvFiles(i).name);
    fileID = csvFiles(i).name(1:6);
    tempT = readtable(filename);

    lastCol = tempT{:, end};
    valueCounts = tabulate(lastCol); 
    tempFraction = zeros(13, 1);

    for j=1:13
        idx = find(valueCounts(:, 1) == combinedT.State(j));

        if ~isempty(idx)
            combinedT.Count(j) = valueCounts(idx, 2);
            tempFraction(j) = valueCounts(idx, 3) / 100;
        end
        
        combinedT.Mouse(j) = i + 204;
        combinedT.Group(j) = "No Bug";
        combinedT.FileID(j) = fileID;
    end

    if  (any(tempFraction == 0))
        tempidx = find(tempFraction == 0);
        tempFraction(tempidx) = realmin;
        [~, maxidx] = max(tempFraction);
        tempFraction(maxidx) = tempFraction(maxidx) - length(tempidx) * realmin; 
    end
    
    combinedT.Fraction = tempFraction;

    fullT = [fullT; combinedT];
end
