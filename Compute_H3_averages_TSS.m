AllFiles = dir('Occupancy*.mat');
noFiles = numel(AllFiles);
for f = 1:noFiles
    load(AllFiles(f).name, 'Occ')
    
    % Normalize occupancy to average 1 for every chromosome (neglecting the
    % reagions with very high occupancy, e.g. center of chr12)
    noChr = numel(Occ);
    for chr = 1:noChr
        tmp = Occ{chr};
        Filter = tmp > 10 * mean(tmp); % detect regions where Occ > 10 * mean(Occ)
        Filter = imdilate(Filter, strel('line', 101, 0));
        tmp(Filter) = nan;
        Occ{chr} = Occ{chr} / nanmean(tmp);
    end
    
    beforeTSS = 1000;
    afterTSS = 1000;

    DatasetLabel = AllFiles(f).name(11:end-4);
    
    chrLen = [230218,813184,316620,1531933,...
        576874,270161,1090940,562643,...
        439888,745751,666816,1078177,...
        924431,784333,1091291,948066];
    
    load('Yeast_annotations.mat', 'ORF', 'Chr', 'Watson', 'TSS');
    
    % Remove genes for which we don't have annotations
    idx = ~isnan(TSS);
    ORF = ORF(idx);
    Chr = Chr(idx);
    Watson = Watson(idx);
    TSS = TSS(idx);

    % Align all TSSs
    noGenes = numel(ORF);
    AlignedOcc = nan(noGenes, 1 + beforeTSS + afterTSS);
    for g = 1:noGenes
        if Watson(g)
            leftEdge = max([TSS(g) - beforeTSS, 1]);
            rightEdge = min([TSS(g) + afterTSS, chrLen(Chr(g))]);
            AlignedOcc(g, beforeTSS + 1 - (TSS(g) - leftEdge)...
                : beforeTSS + 1 + (rightEdge - TSS(g))) = ...
                Occ{1,Chr(g)}(leftEdge : rightEdge);
        else % Watson(g) == false
            leftEdge = max([TSS(g) - afterTSS, 1]);
            rightEdge = min([TSS(g) + beforeTSS, chrLen(Chr(g))]);
            AlignedOcc(g, beforeTSS + 1 - (rightEdge - TSS(g))...
                : beforeTSS + 1 + (TSS(g) - leftEdge)) = ...
                fliplr(Occ{1,Chr(g)}(leftEdge : rightEdge));
        end
    end

    Position = [-beforeTSS:afterTSS]';
    AvgOcc = nanmean(AlignedOcc)';
    
    XLSfilename = [DatasetLabel, '_avg_occ_TSS.xlsx'];
    xlswrite(XLSfilename, {'Position'}, 'Sheet1', 'A1')
    xlswrite(XLSfilename, {'Occupancy'}, 'Sheet1', 'B1')
    xlswrite(XLSfilename, Position, 'Sheet1', 'A2')
    xlswrite(XLSfilename, AvgOcc, 'Sheet1', 'B2')
    
    save([DatasetLabel, '_avg_occ_TSS.mat'], 'Position', 'AvgOcc')
    
    fprintf('File %s processed.\n', AllFiles(f).name)
end