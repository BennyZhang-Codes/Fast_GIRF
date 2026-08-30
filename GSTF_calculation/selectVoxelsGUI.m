function [ validVoxels, selection ] = selectVoxelsGUI(positions, magnitudeData, fieldData, numPE, numSlices, calcChannels, defaultSaveFile)
% Interactive voxel selection for spatially resolved high-order GIRF data.
%
% Arguments:
%   positions:       voxel positions [numVoxels, 3], in meters
%   magnitudeData:   magnitude time courses [numVoxels, numROP]
%   fieldData:       phase-derivative data [numVoxels, numROP, test gradients, ADCs]
%   numPE:           number of phase-encoding steps per in-plane direction
%   numSlices:       number of excited slices
%   calcChannels:    number of spherical-harmonic channels (4, 9, or 16)
%   defaultSaveFile: suggested file name for saving the voxel selection
%
% Returns:
%   validVoxels: logical mask of selected voxels
%   selection:   structure containing the selected mask and QC metrics
%
% The automatic selection follows the legacy Fast_GIRF magnitude criterion.
% The GUI allows manual addition/removal of voxels and reports the rank and
% normalized condition number of the spherical-harmonic probing matrix.

if nargin < 7
    defaultSaveFile = 'voxel_selection.mat';
end

numVoxels = size(positions,1);
if numVoxels ~= numPE*numPE*numSlices
    error('Voxel dimensions do not match numPE x numPE x numSlices.');
end
if size(magnitudeData,1) ~= numVoxels
    error('magnitudeData does not match the number of voxel positions.');
end
if size(fieldData,1) ~= numVoxels
    error('fieldData does not match the number of voxel positions.');
end

%% Automatic voxel selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx1 = min(10,size(magnitudeData,2));
idx2 = min(50,size(magnitudeData,2));
signalMetric = mean(magnitudeData(:,idx1:idx2),2);
maxMetric = max(signalMetric);

if maxMetric > 0
    autoMask = signalMetric >= 0.6*maxMetric;
else
    autoMask = false(numVoxels,1);
end
autoMask = autoMask & isfinite(signalMetric);
selected = autoMask;
currentVoxel = find(selected,1,'first');
if isempty(currentVoxel)
    currentVoxel = 1;
end

% Build the complete probing matrix once. Sub-selection is applied in the
% GUI callbacks, avoiding repeated calls to createProbingMatrix().
probingMatrixAll = createProbingMatrix(positions,calcChannels);

%% Create GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('Name','Fast High-Order GIRF - Voxel Selection', ...
    'NumberTitle','off', 'Color','w', 'Units','normalized', ...
    'Position',[0.08 0.08 0.84 0.80], 'CloseRequestFcn',@closeGUI);

sliceAxes = gobjects(numSlices,1);
sliceImages = gobjects(numSlices,1);

if numSlices <= 2
    nRows = 1;
else
    nRows = 2;
end
nCols = ceil(numSlices/nRows);

leftX = 0.04;
leftW = 0.43;
leftY = 0.18;
leftH = 0.76;
gapX = 0.025;
gapY = 0.07;
axW = (leftW-gapX*(nCols-1))/nCols;
axH = (leftH-gapY*(nRows-1))/nRows;

for sl = 1:numSlices
    row = floor((sl-1)/nCols);
    col = mod(sl-1,nCols);
    pos = [leftX + col*(axW+gapX), ...
           leftY + (nRows-1-row)*(axH+gapY), axW, axH];
    sliceAxes(sl) = axes('Parent',fig,'Position',pos);
    map = getSelectionMap(sl);
    sliceImages(sl) = imagesc(sliceAxes(sl),map.');
    set(sliceImages(sl),'ButtonDownFcn',@(src,event)toggleVoxel(src,event,sl));
    axis(sliceAxes(sl),'image');
    set(sliceAxes(sl),'XTick',1:numPE,'YTick',1:numPE,'YDir','normal');
    xlabel(sliceAxes(sl),'PE 1');
    ylabel(sliceAxes(sl),'PE 2');
    title(sliceAxes(sl),sprintf('Slice %d',sl));
    caxis(sliceAxes(sl),[0 2]);
end

% 0 = rejected, 1 = auto-selected, 2 = manually added
colormap(fig,[0.82 0.82 0.82; 0.42 0.72 0.46; 0.35 0.55 0.86]);

axMagnitude = axes('Parent',fig,'Position',[0.54 0.57 0.42 0.34]);
axField = axes('Parent',fig,'Position',[0.54 0.20 0.42 0.28]);

statsText = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.54 0.10 0.42 0.055], 'BackgroundColor','w', ...
    'HorizontalAlignment','left','FontSize',10,'FontWeight','bold');

uicontrol(fig,'Style','pushbutton','String','Auto Select','Units','normalized', ...
    'Position',[0.05 0.065 0.10 0.055],'Callback',@autoSelect);
uicontrol(fig,'Style','pushbutton','String','Reset','Units','normalized', ...
    'Position',[0.16 0.065 0.08 0.055],'Callback',@resetSelection);
uicontrol(fig,'Style','pushbutton','String','Load','Units','normalized', ...
    'Position',[0.25 0.065 0.08 0.055],'Callback',@loadSelection);
uicontrol(fig,'Style','pushbutton','String','Save','Units','normalized', ...
    'Position',[0.34 0.065 0.08 0.055],'Callback',@saveSelection);
uicontrol(fig,'Style','pushbutton','String','Confirm','Units','normalized', ...
    'Position',[0.83 0.035 0.13 0.07],'FontWeight','bold','Callback',@confirmSelection);

updateGUI();
uiwait(fig);

if isvalid(fig)
    delete(fig);
end

validVoxels = logical(selected(:));
selection = buildSelectionStructure();

%% GUI callbacks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function toggleVoxel(~,~,sliceIndex)
        cp = get(sliceAxes(sliceIndex),'CurrentPoint');
        pe1 = round(cp(1,1));
        pe2 = round(cp(1,2));
        if pe1 < 1 || pe1 > numPE || pe2 < 1 || pe2 > numPE
            return;
        end
        localIndex = pe1 + (pe2-1)*numPE;
        voxelIndex = (sliceIndex-1)*numPE*numPE + localIndex;
        selected(voxelIndex) = ~selected(voxelIndex);
        currentVoxel = voxelIndex;
        updateGUI();
    end

    function autoSelect(~,~)
        selected = autoMask;
        idx = find(selected,1,'first');
        if ~isempty(idx)
            currentVoxel = idx;
        end
        updateGUI();
    end

    function resetSelection(~,~)
        selected = false(numVoxels,1);
        updateGUI();
    end

    function loadSelection(~,~)
        [fileName,pathName] = uigetfile('*.mat','Load voxel selection',defaultSaveFile);
        if isequal(fileName,0)
            return;
        end
        data = load(fullfile(pathName,fileName));
        if isfield(data,'selection') && isfield(data.selection,'validVoxels')
            loadedMask = logical(data.selection.validVoxels(:));
        elseif isfield(data,'validVoxels')
            loadedMask = logical(data.validVoxels(:));
        else
            errordlg('Selected file does not contain a valid voxel mask.','Load selection');
            return;
        end
        if numel(loadedMask) ~= numVoxels
            errordlg('Saved voxel mask does not match the current dataset.','Load selection');
            return;
        end
        selected = loadedMask;
        idx = find(selected,1,'first');
        if ~isempty(idx)
            currentVoxel = idx;
        end
        updateGUI();
    end

    function saveSelection(~,~)
        selection = buildSelectionStructure(); %#ok<NASGU>
        validVoxels = selection.validVoxels; %#ok<NASGU>
        [fileName,pathName] = uiputfile('*.mat','Save voxel selection',defaultSaveFile);
        if isequal(fileName,0)
            return;
        end
        save(fullfile(pathName,fileName),'selection','validVoxels');
    end

    function confirmSelection(~,~)
        [rankA,condA] = calculateMatrixQC(selected);
        if sum(selected) < calcChannels || rankA < calcChannels
            warndlg(sprintf(['The selected voxels do not support the requested field expansion.\n' ...
                'Selected voxels: %d\nRank: %d / %d\nNormalized condition number: %.3g'], ...
                sum(selected),rankA,calcChannels,condA),'Voxel selection');
            return;
        end
        uiresume(fig);
    end

    function closeGUI(~,~)
        % Closing the GUI keeps the current selection. If the selection is
        % rank deficient, fall back to the automatic selection.
        [rankA,~] = calculateMatrixQC(selected);
        if rankA < calcChannels
            selected = autoMask;
        end
        uiresume(fig);
    end

%% GUI update helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function updateGUI()
        for sliceIndex = 1:numSlices
            set(sliceImages(sliceIndex),'CData',getSelectionMap(sliceIndex).');
        end

        cla(axMagnitude);
        plot(axMagnitude,magnitudeData(currentVoxel,:),'LineWidth',1);
        grid(axMagnitude,'on');
        xlabel(axMagnitude,'ADC sample');
        ylabel(axMagnitude,'Magnitude (a.u.)');
        title(axMagnitude,sprintf('Voxel %d - FID magnitude',currentVoxel));

        cla(axField);
        fieldTrace = reshape(fieldData(currentVoxel,:,:,:),size(fieldData,2),[]);
        if ~isempty(fieldTrace)
            plot(axField,fieldTrace(:,1),'LineWidth',1);
        end
        grid(axField,'on');
        xlabel(axField,'ADC sample');
        ylabel(axField,'d\phi/dt (rad/s)');
        title(axField,'Polarity-difference field trace (first test gradient)');

        [rankA,condA] = calculateMatrixQC(selected);
        if isfinite(condA)
            condText = sprintf('%.3g',condA);
        else
            condText = 'Inf';
        end
        set(statsText,'String',sprintf(['Selected voxels: %d    SH channels: %d    ' ...
            'Rank: %d/%d    Normalized cond(A): %s'], ...
            sum(selected),calcChannels,rankA,calcChannels,condText));
        drawnow;
    end

    function map = getSelectionMap(sliceIndex)
        idx = (sliceIndex-1)*numPE*numPE + (1:numPE*numPE);
        map = zeros(numPE,numPE);
        selectedSlice = selected(idx);
        autoSlice = autoMask(idx);
        values = zeros(numPE*numPE,1);
        values(selectedSlice & autoSlice) = 1;
        values(selectedSlice & ~autoSlice) = 2;
        map(:) = values;
    end

    function [rankA,condA] = calculateMatrixQC(mask)
        A = probingMatrixAll(mask,:);
        if isempty(A)
            rankA = 0;
            condA = Inf;
            return;
        end
        rankA = rank(A);
        columnNorm = sqrt(sum(abs(A).^2,1));
        columnNorm(columnNorm == 0) = 1;
        An = A ./ columnNorm;
        if size(An,1) < size(An,2) || rankA < size(An,2)
            condA = Inf;
        else
            condA = cond(An);
        end
    end

    function selectionOut = buildSelectionStructure()
        [rankA,condA] = calculateMatrixQC(selected);
        selectionOut = struct();
        selectionOut.validVoxels = logical(selected(:));
        selectionOut.positions = positions;
        selectionOut.signalMetric = signalMetric;
        selectionOut.autoMask = logical(autoMask(:));
        selectionOut.calcChannels = calcChannels;
        selectionOut.rankA = rankA;
        selectionOut.condA_normalized = condA;
        selectionOut.numPE = numPE;
        selectionOut.numSlices = numSlices;
    end
end
