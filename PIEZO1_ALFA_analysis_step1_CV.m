%% Script for the analysis of MINFLUX signals from DNA-Paint labelled ALFA-tagged PIEZO1
    % Dependencies: DBSCAN.m; LoadDataFromMFXFile.m; LoadImage2Array.m;
    % SetALFAanalysisParameters.m, PlotRawData.m; PIEZO1Superparticle;
    % PlotTrimerAnalysisResult.m ExtractClusterRAW
    
    clear all;
    close all;
TrimWithOneNeighbor = 0
Dist4 = 0;
Dist4List = [];
%% choose which data to analyse: 0 = soma, 1 = neurite etc.
    DataSource = 0;
    if DataSource == 0
        result.RefTable = 'FileListP1Halo_ALFA_CTL.xlsx';
    
        elseif DataSource == 1
        result.RefTable = 'FileListWedgeDel_ALFA_CTL.xlsx';

    end

%%  %%%%%%%%%%%%%%%%%%%%%%%% –CHOOSE OPTIONS – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotRaw             = false; % set to 'true' for raw data summary plot
    filter              = true; % set to 'true' if data should be filtered 
    plotFiltered        = true;
    FindClusters        = false;
    plotTraceAVGs       = true;
    MakeConfocalOverlay = false;
    aggregateTraces     = false;
    TrimTrace           = true;
    
    % create empty result arrays    
    result.Trimers = [];
    result.Angles = [];
    result.InterBlades =[];
    result.Counts = [];
    result.NumTraces = [];
    result.InterBladesWithNN = [];
    result.TrimerCount=0;
    result.TrimerCoord = [];
    result.LocPerTrace = [];
    numTracesperHour = [];
    result.NearestNeighborInCluster = [];
    result.NearestNeighborOutsideCluster = [];
    AllStdTrimer = [];
    result.AllOpenTOP = [];
    result.AllOpenBOTTOM = [];
%%  %%%%%%%%%%%%%%%%%%%%%%%% – Load data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if matches(result.RefTable,'FileListNeuritesFINAL.xlsx')
        result.color = [0.6 0.4 0.8];
    else
        result.color = [0 0.5 0];
    end
    % read reftable    
    data = readtable(result.RefTable);

%% loop through all files or select single file    
for FileIdx = 4%:size(data,1)
% for FileIdx =  1 % 1-9 for neurite; 1-10 soma (for Trimers: neurite #10 example, soma #2; for Clusters: soma #6)
    myfile = string(data{FileIdx,1});
    UpperZ = [data{FileIdx,2}];
    LowerZ = [data{FileIdx,3}];
    include = [data{FileIdx,4}];
        if include==0
            continue
        end
    % [traces_RAW]=LoadDataFromMFXFile_V3(myfile);
    [traces_RAW]=LoadDataFromMFXFile_V3(myfile);
%%  %%%%%%%%%%%%%%%%%%%%%%%% – set analysis parameters – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [p]=SetALFAanalysisParameters(); % set analysis paramters
    p.z_Upper_threshold = 1e9*UpperZ;  % Z-filt threshold
    p.z_Lower_threshold = 1e9*LowerZ;  % Z-filt threshold

%%  %%%%%%%%%%%%%%%%%%%%%%%% – create empty results arrays for iteration – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Trimers = [];
    TrimWithNN = [];
    InterBlade = [];
    Angles =[];
    TrimerCoord = [];
        tempALL = [];

%%  %%%%%%%%%%%%%%%%%%%%%%%% – plot raw data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotRaw
        traces_RAW(:,3) = 0.7*traces_RAW(:,3);
        PlotRawData(traces_RAW,myfile);
    end

%%  %%%%%%%%%%%%%%%%%%%%%%%% – filtering – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotRaw
    traces_FILT = traces_RAW;
    traces_FILT(:,1:3)=1e9*traces_FILT(:,1:3);
    else
    traces_FILT = traces_RAW;
    traces_FILT(:,1:3)=1e9*traces_FILT(:,1:3);
    traces_FILT(:,3) = traces_FILT(:,3)*0.7; % apply scaling factor of 0.7 to z-coordinates (ref. Gwosch et al. Nat Methods 2020, Ostersehlt LM Net Methods 2022)
    % [traced_FILT_drift] = mydriftcorrection(traces_FILT);
    % traces_FILT(:,1:3) = traced_FILT_drift;
    end
    %% aggregate traces
    if aggregateTraces
    p.loc_per_trace_threshold = 0; traces_Aggregate = []; j=1;
    [a, ~, ~] = unique(traces_FILT(:,4));
    agg = 3;
    for i=1:size(a,1)
    temp = traces_FILT(traces_FILT(:,4)==a(i,1),1:8);
        for k=1:agg:size(temp,1)
            if size(temp,1)<(agg+1)
                continue
            elseif k+(agg-1)>size(temp,1)
                traces_Aggregate(j,1:8) = mean(temp(k:size(temp,1),:),1);
            elseif k+(agg-1)<=size(temp,1)
                traces_Aggregate(j,1:8) = mean(temp(k:k+(agg-1),:));
            end
            j=j+1;
        end
    end
    traces_FILT = traces_Aggregate;
    clear a agg temp;
    end

    if TrimTrace
        j=1;
        [a, ~, ~] = unique(traces_FILT(:,4));
        for i=1:size(a,1)
        temp = traces_FILT(traces_FILT(:,4)==a(i,1),1:8);
        if size(temp,1)>1
        temp(1:2,:)=[];
        tempALL = cat(1,tempALL,temp);
        end
        end
        traces_FILT = tempALL;
        clear j temp tempALL;
    end


    % if TrimTrace
    %     j=1;
    %     [a, ~, ~] = unique(traces_FILT(:,4));
    %     for i=1:size(a,1)
    %     temp = traces_FILT(traces_FILT(:,4)==a(i,1),1:8);
    %     if size(temp,1)>2
    %     tempmean = [mean(temp(:,1),1) mean(temp(:,2),1) mean(temp(:,3),1)];
    %     dist1 = temp(1,1:3)-tempmean;
    %     dist2 = temp(2,1:3)-tempmean;
    %     dist3 = temp(3,1:3)-tempmean;
    %     dist1abs = sqrt(dist1(:,1)^2+dist1(:,2)^2+dist1(:,3)^2);
    %     dist2abs = sqrt(dist2(:,1)^2+dist2(:,2)^2+dist2(:,3)^2);
    %     dist3abs = sqrt(dist3(:,1)^2+dist3(:,2)^2+dist3(:,3)^2);
    %     if dist1abs > 10 && dist2abs > 10 && dist3abs > 10
    %         temp(1:3,:)=[];
    %     elseif dist1abs > 10 && dist2abs > 10
    %         temp(1:2,:)=[];
    %     elseif dist1abs > 10
    %         temp(1,:)=[];
    %     end
    %     tempALL = cat(1,tempALL,temp);
    %     end
    %     end
    %     traces_FILT = tempALL;
    %     clear j temp tempALL;
    % end




    if filter
        % Z-filtering
        traces_FILT(traces_FILT(:, 3) < p.z_Lower_threshold, :)= []; % filter by Z upper limit 160,:) = []; 
        traces_FILT(traces_FILT(:, 3) > p.z_Upper_threshold, :)= []; % filter by Z lower limit
        
        % XY-filtering
        %% remove traces originating from cell body
        
            if matches(result.RefTable,'FileListNeuritesFINAL.xlsx')
                if FileIdx == 5
                    mymessage = 'gefunden'
                    myothermessage = myfile
                    traces_FILT(traces_FILT(:, 1) > -3500,:)= [];  % optional X-right
                end
            end
        

        % efo, cfr filtering

        traces_FILT = traces_FILT(traces_FILT(:,5) <= p.efo_threshold, :);                          % filter by efo
        traces_FILT = traces_FILT(traces_FILT(:,6) <= p.cfr_threshold, :);                          % filter by cfr
            
        % time filtering 
        traces_FILT = traces_FILT(traces_FILT(:,8) <= p.time_threshold, :);                         % filter by recording time
        
        %% filter by standard deviation
        [~, ~, id_tid] = unique(traces_FILT(:,4));
        StDevALL=[accumarray(id_tid,traces_FILT(:,1),[],@std) accumarray(id_tid,traces_FILT(:,2),[],@std) accumarray(id_tid,traces_FILT(:,3),[],@std)];
       
        % create arrays for graph coloring purposes only 
        STDforeachiteration = zeros(size(id_tid,1),3);
        for Tidx = 1:size(id_tid,1)
            STDforeachiteration(Tidx,:)=StDevALL(id_tid(Tidx,1),:);
        end
        % StdMIN = max(STDforeachiteration,[],2);
        % % StdMIN = min(STDforeachiteration,[],2);
        StdMIN = mean(STDforeachiteration,2);
       
        % new filtering code
        traces_FILT = traces_FILT(StdMIN(:,1)<p.stdev_trace_threshold,:);




        %% filter by localisations per trace 
        [uv_tid, ~, id_tid] = unique(traces_FILT(:,4));                                                     % Unique elements and locations in third column  
        n_tid = histcounts(id_tid,size(uv_tid,1)); % How many of each?
        % create arrays for graph coloring purposes only 
        LocPerTrace=n_tid';
        Numtraceforeachiteration = id_tid;
        for Tidx = 1:size(id_tid,1)
            Numtraceforeachiteration(Tidx)=LocPerTrace(id_tid(Tidx,1));
        end
        

        % old filtering code
        traces_FILT = traces_FILT(ismember(traces_FILT(:,4), uv_tid(n_tid > p.loc_per_trace_threshold)),:); % Keep ones with more than threshold.
        STDforeachiteration = STDforeachiteration(Numtraceforeachiteration(:,1)>p.loc_per_trace_threshold,:);
        Numtraceforeachiteration = Numtraceforeachiteration(Numtraceforeachiteration(:,1)>p.loc_per_trace_threshold,:);
        
        numTracesperHour(FileIdx,1)=size(LocPerTrace,1);
        % LocHist = histcounts(LocPerTrace,"BinWidth",1)';
        result.LocPerTrace = cat(1,result.LocPerTrace,LocPerTrace);
        n_tid=n_tid(n_tid(:)>p.loc_per_trace_threshold)';
                 


    else % the following operations must be performed when no filtering is applied in order for the script to work
        p.loc_per_trace_threshold=0;
        [uv_tid, ~, id_tid] = unique(traces_FILT(:,4));
        StDevALL=[accumarray(id_tid,traces_FILT(:,1),[],@std) accumarray(id_tid,traces_FILT(:,2),[],@std) accumarray(id_tid,traces_FILT(:,3),[],@std)];
        n_tid = histcounts(id_tid,"BinWidth",1);                                                            % How many of each?
        traces_FILT = traces_FILT(ismember(traces_FILT(:,4), uv_tid(n_tid > p.loc_per_trace_threshold)),:); % Keep ones with more than threshold.
        n_tid=n_tid(n_tid(:)>p.loc_per_trace_threshold)';
        myNumLoc = size(n_tid);
        result.NumTraces = cat(1,result.NumTraces, myNumLoc);
        result.Counts = cat(1,result.Counts, n_tid);
    end

    % clean up
    [~, ~, Std_tid] = unique(traces_FILT(:,4));
    StDevALL=[accumarray(Std_tid,traces_FILT(:,1),[],@std) accumarray(Std_tid,traces_FILT(:,2),[],@std) accumarray(Std_tid,traces_FILT(:,3),[],@std)];
       
    clear StdMIN LocPerTrace



    %% plot filered data - optional
    if plotFiltered
        figure; % 3D
            scatter3(traces_FILT(:,1), traces_FILT(:,2), traces_FILT(:,3),5,traces_FILT(:,3), 'filled');  % color options: ClusterIDs, myXYZ(:,5)
            colormap hot; 
             axis equal; 
            title(strcat(myfile,' filtered data'));
            % view(90,0);
            set(gcf,'renderer','Painters');
        figure; % 2D
            scatter(traces_FILT(:,1), traces_FILT(:,2), 5,traces_FILT(:,3), 'filled');  % color options: ClusterIDs, myXYZ(:,5)
            colormap parula; axis equal; title(strcat(myfile,' filtered data'));
            % view(90,0);
            set(gcf,'renderer','Painters');
    end

%%  %%%%%%%%%%%%%%%%%%%%%%%% – calculate center of mass for each trace – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, ~, id_tid] = unique(traces_FILT(:,4)); 
    CenterX = [accumarray(id_tid,traces_FILT(:,1),[],@mean)];
    CenterY = [accumarray(id_tid,traces_FILT(:,2),[],@mean)];
    CenterZ = [accumarray(id_tid,traces_FILT(:,3),[],@mean)];
    OrigTID = [accumarray(id_tid,traces_FILT(:,4),[],@mean)];
    traces_AVG = [CenterX CenterY CenterZ];

    % clean up
    clear CenterX CenterY CenterZ;
    
%% %%%%%%%%%%%%%%%%%%%%%%%% - load confocal image and plot overlay with MFX data    
    if MakeConfocalOverlay 
        Xshift = [data{FileIdx,5}];
        Yshift = [data{FileIdx,6}];
        PixelSize = [data{FileIdx,7}];
        myImage = string(data{FileIdx,8});
        myXYZ = 1e-9*traces_AVG;
        [ImageArray, ImageArrayFULL] = LoadImage2Array(myImage,myXYZ,Xshift,Yshift,PixelSize);
        figure;
        tiledlayout(1,2);
        nexttile
            scatter(ImageArray(:,1),ImageArray(:,2),200,ImageArray(:,3),'square', 'filled')
            axis equal; axis tight; colormap gray; hold on;
            scatter(myXYZ(:,1), myXYZ(:,2),50,'red','filled');
        nexttile
        ClusterIDsRAW = dbscan(traces_FILT(:,1:3),100,100);
        scatter3(traces_FILT(:,1), traces_FILT(:,2), traces_FILT(:,3),10,ClusterIDsRAW(:,1), 'filled');  % color options: ClusterIDs, myXYZ(:,5)
        colormap parula; axis equal; title(strcat(myfile,' filtered data'));
    
        % plot full size image - optional
        scatter(ImageArrayFULL(:,1),ImageArrayFULL(:,2),1000,ImageArrayFULL(:,3),'square', 'filled')
        axis equal; colormap gray; hold on;
        scatter(myXYZ(:,1), myXYZ(:,2),50,'red','filled');
    end

%%  %%%%%%%%%%%%%%%%%%%%%%%% – run DBSCAN to find signals from the same fluorophore – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ClusterIDs = dbscan(traces_AVG(:,1:3),p.MaxDistProtamers,p.MinNumProtamers);
    
    % add column with cluster ID to main coordinates table
    traces_AVG_noMerge = [traces_AVG ClusterIDs n_tid OrigTID];
    
    % split clust_XYZ into 2 submatrices - one with single signals and one with multiple signals per channel 
    singleloc = traces_AVG_noMerge(traces_AVG_noMerge(:,4) == -1,:);
    multiloc = traces_AVG_noMerge(traces_AVG_noMerge(:,4) > -1,:);
    
    % calculate averages of clusters (i.e. signals from one channel)
    [~, ~, id] = unique(multiloc(:,4));
    MyMeanX = [accumarray(id,multiloc(:,1),[],@mean)];
    MyMeanY = [accumarray(id,multiloc(:,2),[],@mean)];
    MyMeanZ = [accumarray(id,multiloc(:,3),[],@mean)];
    OrigTIDmulti = [accumarray(id,multiloc(:,6),[],@min)];
    ClID = [accumarray(id,multiloc(:,4),[],@mean)];
    numLocs = [accumarray(id,multiloc(:,5),[],@sum)];
    multiloc = [MyMeanX, MyMeanY, MyMeanZ, ClID, numLocs,OrigTIDmulti];
    
    % merge submatrices into single matrix
    traces_AVG = cat(1,singleloc(:,1:6),multiloc(:,1:6));

    % clean up
    traces_AVG_noMerge = cat(2,traces_AVG_noMerge,StDevALL);
    clear MyMeanX MyMeanY MyMeanZ singleloc multiloc ClusterIDs ClID StDevALL numLocs

%%  %%%%%%%%%%%%%%%%%%%%%%%% – find trimers – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create distance matrix
    DistMatrix=squareform(pdist(traces_AVG(:,1:3)));
    [H, Hidx] = sort(DistMatrix,1);
    NearestNeighbor = H(2,:)';

    % find trimers
    for i=1:length(H)
        % check if third nearest neighbors are further away than p.MinTrimerDist
        if(H(4,i)>p.MinTrimerDist && H(4,Hidx(2,i))>p.MinTrimerDist && H(4,Hidx(3,i))>p.MinTrimerDist)
            % check if all three triangle sides are shorter than p.MaxInterBladeDist
            if(H(2,i)<p.MaxInterBladeDist && H(3,i)<p.MaxInterBladeDist && DistMatrix(Hidx(2,i),Hidx(3,i))<p.MaxInterBladeDist)    
                 Trimers = cat(1,Trimers, sort(Hidx(1:3,i))');
                 % Get NN for each trimer: 4th, 5th, 6th and average NN dist for 4 to 6
                 TrimWithNN = cat(1, TrimWithNN, [H(2,i) H(3,i) DistMatrix(Hidx(2,i),Hidx(3,i)) mean(H(4,i)) mean(H(5,i)) mean(H(6,i)) mean(H(4:6,i))]);
                 InterBlade = cat(1, InterBlade, [H(2,i) H(3,i) DistMatrix(Hidx(2,i),Hidx(3,i))]);
                 TrimerCoord = cat(1,TrimerCoord, [traces_AVG(Hidx(1,i),1:3) traces_AVG(Hidx(2,i),1:3) traces_AVG(Hidx(3,i),1:3)]);
            end
        else
            if(H(2,i)<p.MaxInterBladeDist && H(3,i)<p.MaxInterBladeDist && DistMatrix(Hidx(2,i),Hidx(3,i))<p.MaxInterBladeDist)
            myCount = 0;
            if(H(4,i)<p.MinTrimerDist)
                myCount=myCount+1;
                Dist4 = H(4,i);
            end
            if(H(4,Hidx(2,i))<p.MinTrimerDist)
                myCount=myCount+1;
                Dist4 = H(4,Hidx(2,i));
            end
            if(H(4,Hidx(3,i))<p.MinTrimerDist)
                myCount=myCount+1;
                Dist4 = H(4,Hidx(3,i));
            end  

            if(myCount==1)
                TrimWithOneNeighbor = TrimWithOneNeighbor+1;
                Dist4List = cat(1, Dist4List,Dist4);
            end
            end
        end
    end

    % evaluate trimers
    if size(Trimers,1)>0
        % find and remove duplicates
        [~,ia,~] = unique(Trimers(:,1:3),'rows'); % find duplicates
        Trimers = Trimers(ia,:); % remove duplicates
        TrimWithNN = TrimWithNN(ia,:);

        InterBlade = InterBlade(ia,:); % remove duplicates
        TrimerCoord = TrimerCoord(ia,:);

        % calculate angles
        for j = 1:size(Trimers,1)
            PointA = traces_AVG(Trimers(j,1),1:3);
            PointB = traces_AVG(Trimers(j,2),1:3);
            PointC = traces_AVG(Trimers(j,3),1:3);
            [TrimerAngles] = CalcTrimerAngle1(PointA, PointB, PointC);           
            Angles = cat(1,Angles, TrimerAngles);
        end
        MaxAngleTrimer = max(Angles,[],2);

        % filter by angle
            maxAngle = 120;
            Trimers=Trimers(MaxAngleTrimer(:,1)<maxAngle,:);        
            InterBlade=InterBlade(MaxAngleTrimer(:,1)<maxAngle,:);
            TrimerCoord=TrimerCoord(MaxAngleTrimer(:,1)<maxAngle,:);
            Angles=Angles(MaxAngleTrimer(:,1)<maxAngle,:);
            TrimWithNN=TrimWithNN(MaxAngleTrimer(:,1)<maxAngle,:); 

        
            % get STD of trimers
                if ~exist('IndivTrimersRAW_avg', 'var')
                IndivTrimersRAW_avg = struct();
                IndivTrimersRAW_avg.Trimer.A = [];
                IndivTrimersRAW_avg.Trimer.B = [];
                 IndivTrimersRAW_avg.Trimer.C = [];
                end
            for n = 1:size(Trimers,1)
                idA = traces_AVG(Trimers(n,1),6);
                idB = traces_AVG(Trimers(n,2),6);
                idC = traces_AVG(Trimers(n,3),6);
                SizeA = size(traces_FILT(traces_FILT(:,4)==idA,1),1);
                SizeB = size(traces_FILT(traces_FILT(:,4)==idB,1),1);
                SizeC = size(traces_FILT(traces_FILT(:,4)==idC,1),1);
                StdA = [std(traces_FILT(traces_FILT(:,4)==idA,1)) std(traces_FILT(traces_FILT(:,4)==idA,2)) std(traces_FILT(traces_FILT(:,4)==idA,3)) SizeA];
                StdB = [std(traces_FILT(traces_FILT(:,4)==idB,1)) std(traces_FILT(traces_FILT(:,4)==idB,2)) std(traces_FILT(traces_FILT(:,4)==idB,3)) SizeB];
                StdC = [std(traces_FILT(traces_FILT(:,4)==idC,1)) std(traces_FILT(traces_FILT(:,4)==idC,2)) std(traces_FILT(traces_FILT(:,4)==idC,3)) SizeC];
                StdTrimer = [StdA; StdB; StdC];
                AllStdTrimer = cat(1,AllStdTrimer, StdTrimer);
                % extract raw coordinates
                CoordA = traces_FILT(traces_FILT(:,4)==idA,1:3);
                CoordB = traces_FILT(traces_FILT(:,4)==idB,1:3);
                CoordC = traces_FILT(traces_FILT(:,4)==idC,1:3);
                IndivTrimersRAW.(strcat('Trimer',string(FileIdx),'sub',string(n))).A = CoordA;
                IndivTrimersRAW.(strcat('Trimer',string(FileIdx),'sub',string(n))).B = CoordB;
                IndivTrimersRAW.(strcat('Trimer',string(FileIdx),'sub',string(n))).C = CoordC;
                % extract avg raw coordinates
                CoordA_avg = mean(CoordA); 
                CoordB_avg = mean(CoordB);
                CoordC_avg = mean(CoordC);
                % IndivTrimersRAW_avg.(strcat('Trimer',string(FileIdx),'sub',string(n))).A = CoordA_avg;
                % IndivTrimersRAW_avg.(strcat('Trimer',string(FileIdx),'sub',string(n))).B = CoordB_avg;
                % IndivTrimersRAW_avg.(strcat('Trimer',string(FileIdx),'sub',string(n))).C = CoordC_avg;
                IndivTrimersRAW_avg.Trimer.A = [IndivTrimersRAW_avg.Trimer.A; CoordA_avg];
                IndivTrimersRAW_avg.Trimer.B = [IndivTrimersRAW_avg.Trimer.B; CoordB_avg];
                IndivTrimersRAW_avg.Trimer.C = [IndivTrimersRAW_avg.Trimer.C; CoordC_avg];

            end
            % clean up
            clear idA idB idC StdA StdB StdC StdTrimer
        
        % write results
            result.Trimers = cat(1,result.Trimers, Trimers);
            result.InterBlades = cat(1, result.InterBlades, InterBlade);
            result.TrimerCoord = cat(1,result.TrimerCoord,TrimerCoord);
            result.Angles = cat(1, result.Angles, Angles);
            result.InterBladesWithNN = cat(1, result.InterBladesWithNN,TrimWithNN);
    end

    % clean up
    clear TrimWithNN PointA PointB PointC InterBlade TrimerCoord Angles MaxAngleTrimer maxAngle TrimerAngles

%%  %%%%%%%%%%%%%%%%%%%%%%%% – plot processed data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotTraceAVGs
        figure;
        scatter3(traces_AVG(:,1), traces_AVG(:,2), traces_AVG(:,3),20, traces_AVG(:,3), 'filled');
        colormap parula; axis equal; title(myfile); hold on;
    
        % add trimers to plot
        for k=1:size(Trimers,1)
            result.TrimerCount=result.TrimerCount+1;
            scatter3(traces_AVG(Trimers(k,1),1), traces_AVG(Trimers(k,1),2), traces_AVG(Trimers(k,1),3),50, result.color, 'filled')
            scatter3(traces_AVG(Trimers(k,2),1), traces_AVG(Trimers(k,2),2), traces_AVG(Trimers(k,2),3),50, result.color, 'filled')
            scatter3(traces_AVG(Trimers(k,3),1), traces_AVG(Trimers(k,3),2), traces_AVG(Trimers(k,3),3),50, result.color, 'filled')
            text(traces_AVG(Trimers(k,1),1), traces_AVG(Trimers(k,1),2), traces_AVG(Trimers(k,1),3),sprintfc(' %d',result.TrimerCount))
        end
    end

%%  %%%%%%%%%%%%%%%%%%%%%%%% – identify and plot clusters - OPTIONAL  – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if FindClusters
        ClusterIDs2 = dbscan(traces_AVG(:,1:3),p.MaxDistWithinCluster,p.MinNumChannelsPerCluster);        

        ClusterSize = grouptransform(ClusterIDs2,ClusterIDs2,@numel);
        traces_AVG_Cluster = traces_AVG(ClusterIDs2>-1,:);
        traces_AVG_noCluster = traces_AVG(ClusterIDs2< 0,:);
        NNClusColor = ClusterIDs2(ClusterIDs2>-1,:);
        result.NearestNeighborInCluster = cat(1,result.NearestNeighborInCluster,NearestNeighbor(ClusterIDs2>-1,:));
        result.NearestNeighborOutsideCluster = cat(1,result.NearestNeighborOutsideCluster,NearestNeighbor(ClusterIDs2< 0,:));
        % filter by nearest neighbor < 40 nm
        result.NearestNeighborInCluster = result.NearestNeighborInCluster(result.NearestNeighborInCluster(:,1)<40);
        result.NearestNeighborOutsideCluster = result.NearestNeighborOutsideCluster(result.NearestNeighborOutsideCluster(:,1)<40);
        
        % extract data of individual clusters to structure
        for k =1:max(ClusterIDs2)
            % calculate nearest neighbor distances in cluster
            temp = traces_AVG(ClusterIDs2==k,:);
            DistMatrix=squareform(pdist(temp(:,1:3)));
            [H, Hidx] = sort(DistMatrix,1);
            NNClus = H(2,:)';
            IndivClusters.(strcat('Cluster_',string(FileIdx),'sub',string(k))) = cat(2,traces_AVG(ClusterIDs2==k,:),NNClus(:,1));
            %% extract cluster raw locs
                clusterMEANs = traces_AVG(ClusterIDs2==k,:);
                [ClusterRAW] = ExtractClusterRAW(clusterMEANs, traces_FILT);
                IndivClustersRAW.(strcat('Cluster_',string(FileIdx),'sub',string(k))) = ClusterRAW;
        end

        % create cluster label
        ClusterLabelX = [accumarray(NNClusColor,traces_AVG_Cluster(:,1),[],@mean)]+100;
        ClusterLabelY = [accumarray(NNClusColor,traces_AVG_Cluster(:,2),[],@mean)];
        ClusterLabelZ = [accumarray(NNClusColor,traces_AVG_Cluster(:,3),[],@mean)];
        ClusterColorNew = [accumarray(NNClusColor,NNClusColor(:,1),[],@mean)];

        figure;
        % plot clusters
        scatter3(traces_AVG_Cluster(:,1), traces_AVG_Cluster(:,2), traces_AVG_Cluster(:,3),20, NNClusColor(:,1), 'filled'); %'red', 'filled');
        colormap jet; hold on;

        % add cluster labels
        text(ClusterLabelX(:,1), ClusterLabelY(:,1), ClusterLabelZ(:,1),sprintfc(' %d',ClusterColorNew(:,1)),'FontSize',20)
        
        % plot non-clusters
        scatter3(traces_AVG_noCluster(:,1), traces_AVG_noCluster(:,2), traces_AVG_noCluster(:,3),20, [0.6 0.6 0.6], 'filled'); 
        axis equal;

        % clean up
        clear NNClus temp ClusterIDs2 ClusterLabelX ClusterLabelY ClusterLabelZ clusterMEANs ClusterRAW ClusterSize ClusterColorNew NearestNeighbor NNClusColor
    end

end % end of main loop

%% - post processing - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine filtered data into single array
traces_FILT = cat(2,traces_FILT,STDforeachiteration, Numtraceforeachiteration);

% create array with distances of nearest neighbors in trimers
result.NearestNeighborInTrimer = sort(result.InterBlades,2);
result.NearestNeighborInTrimer = cat(1,result.NearestNeighborInTrimer(:,1),result.NearestNeighborInTrimer(:,2),result.NearestNeighborInTrimer(:,1));

% plot trimer analysis results
[result] = PlotTrimerAnalysisResult(result);

% find pit-shaped clusters
if FindClusters
[result.NearestNeighborPitCluster, result.ClusterOpeningRadius, result.ClusterDepth,  result.AllOpenTOP, result.AllOpenBOTTOM] = FindPitShapedClusters(IndivClusters, IndivClusters);
end

%% – final clean up – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear uv_tid uv n_tid
% clear plotRaw  filter...
%     plotFiltered FindClusters...
%     FitSurface plotTraceAVGs MakeConfocalOverlay
% clear i j k K ia id id_tid H Hidx FileIdx
% clear myfile UpperZ LowerZ include result.TrimerCount Tidx DistMatrix Trimers data
% clear Std_tid STDforeachiteration DataSource Numtraceforeachiteration NearestNeighbor aggregateTraces TrimTrace SizeA SizeB SizeC 

