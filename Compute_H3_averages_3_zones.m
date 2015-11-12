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
    
    beforeNDRcenter = 1000;
    afterNDRcenter = 1000;

    DatasetLabel = AllFiles(f).name(11:end-4);
    
    chrLen = [230218,813184,316620,1531933,...
        576874,270161,1090940,562643,...
        439888,745751,666816,1078177,...
        924431,784333,1091291,948066];
    
    load('Yeast_annotations.mat', 'ORF', 'Chr', 'Watson', 'NDRcenter', 'NDRwidth');
    
    % Remove genes for which we don't have annotations
    idx = (~isnan(NDRwidth) & (~isnan(NDRcenter)));
    ORF = ORF(idx);
    Chr = Chr(idx);
    Watson = Watson(idx);
    NDRcenter = NDRcenter(idx);
    NDRwidth = NDRwidth(idx);

    % Align all NDR centers
    noGenes = numel(ORF);
    AlignedOcc = nan(noGenes, 1 + beforeNDRcenter + afterNDRcenter);
    for g = 1:noGenes
        if Watson(g)
            leftEdge = max([NDRcenter(g) - beforeNDRcenter, 1]);
            rightEdge = min([NDRcenter(g) + afterNDRcenter, chrLen(Chr(g))]);
            AlignedOcc(g, beforeNDRcenter + 1 - (NDRcenter(g) - leftEdge)...
                : beforeNDRcenter + 1 + (rightEdge - NDRcenter(g))) = ...
                Occ{1,Chr(g)}(leftEdge : rightEdge);
        else % Watson(g) == false
            leftEdge = max([NDRcenter(g) - afterNDRcenter, 1]);
            rightEdge = min([NDRcenter(g) + beforeNDRcenter, chrLen(Chr(g))]);
            AlignedOcc(g, beforeNDRcenter + 1 - (rightEdge - NDRcenter(g))...
                : beforeNDRcenter + 1 + (NDRcenter(g) - leftEdge)) = ...
                fliplr(Occ{1,Chr(g)}(leftEdge : rightEdge));
        end
    end

    % Sort the genes according to NDRwidths
    [NDRwidth, idx] = sort(NDRwidth);
    SortedORF = ORF(idx);
    AlignedOcc = AlignedOcc(idx, :);
    
    % Compute the average occupancy at the 3 regions
    AvgOccMinus1 = nan(noGenes, 1);
    AvgOccNDR = nan(noGenes, 1);
    AvgOccPlus1 = nan(noGenes, 1);
    AvgOcc3zones = nan(noGenes, 1);
    for g = 1:noGenes
        AvgOccMinus1(g) = mean(AlignedOcc(g, 1+beforeNDRcenter-round(NDRwidth(g)/2)-147:1+beforeNDRcenter-round(NDRwidth(g)/2)-1));
        AvgOccNDR(g) = mean(AlignedOcc(g, 1+beforeNDRcenter-round(NDRwidth(g)/2):1+beforeNDRcenter+round(NDRwidth(g)/2)));
        AvgOccPlus1(g) = mean(AlignedOcc(g, 1+beforeNDRcenter+round(NDRwidth(g)/2)+1:1+beforeNDRcenter+round(NDRwidth(g)/2)+147));
        AvgOcc3zones(g) = mean(AlignedOcc(g, 1+beforeNDRcenter-round(NDRwidth(g)/2)-147:1+beforeNDRcenter+round(NDRwidth(g)/2)+147));
    end
    
    XLSfilename = [DatasetLabel, '_avg_occ_statistics.xlsx'];
    xlswrite(XLSfilename, {'ORF'}, 'Sheet1', 'A1')
    xlswrite(XLSfilename, {'-1'}, 'Sheet1', 'B1')
    xlswrite(XLSfilename, {'NDR'}, 'Sheet1', 'C1')
    xlswrite(XLSfilename, {'+1'}, 'Sheet1', 'D1')
    xlswrite(XLSfilename, {'-1,NDR,+1'}, 'Sheet1', 'E1')
    xlswrite(XLSfilename, SortedORF, 'Sheet1', 'A2')
    xlswrite(XLSfilename, AvgOccMinus1, 'Sheet1', 'B2')
    xlswrite(XLSfilename, AvgOccNDR, 'Sheet1', 'C2')
    xlswrite(XLSfilename, AvgOccPlus1, 'Sheet1', 'D2')
    xlswrite(XLSfilename, AvgOcc3zones, 'Sheet1', 'E2')
    
    save([DatasetLabel, '_avg_occ_statistics.mat'], 'SortedORF', 'AvgOccMinus1', 'AvgOccNDR', 'AvgOccPlus1', 'AvgOcc3zones')
    
    fprintf('File %s processed.\n', AllFiles(f).name)
end