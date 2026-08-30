function plotHighOrderGSTFs(H, ax_name, fLimit_kHz)
% Plot grouped High-Order GSTF magnitude and phase.
%
% Arguments:
%   H:          GIRF/GSTF object containing gstf and f_axis
%   ax_name:    driven gradient axis: 'x', 'y', or 'z'
%   fLimit_kHz: optional symmetric frequency limit in kHz (default 5)
%
% The Fast_GIRF spatial basis is ordered as:
%   1: B0
%   2-4: X, Y, Z
%   5-9: XY, YZ, Z2, XZ, X2-Y2
%   10-16: Y(3X2-Y2), XYZ, YZ2, Z3, XZ2, Z(X2-Y2), X(X2-3Y2)
%
% With a first-order input gradient, GSTF units are:
%   B0 term:       m
%   first order:   dimensionless
%   second order:  1/m
%   third order:   1/m^2

if nargin < 3 || isempty(fLimit_kHz)
    fLimit_kHz = 5;
end

if ~ismember(ax_name,{'x','y','z'})
    error('ax_name must be ''x'', ''y'', or ''z''.');
end

numChannels = size(H.gstf,2);
if ~ismember(numChannels,[4,9,16])
    error('High-Order GSTF must contain 4, 9, or 16 spatial channels.');
end

termNames = {'B0','X','Y','Z', ...
    'XY','YZ','Z2','XZ','X2-Y2', ...
    'Y(3X2-Y2)','XYZ','Y(4Z2-X2-Y2)','Z(2Z2-3X2-3Y2)', ...
    'X(4Z2-X2-Y2)','Z(X2-Y2)','X(X2-3Y2)'};

orderChannels = {1, 2:min(4,numChannels), 5:min(9,numChannels), 10:min(16,numChannels)};
orderNames = {'B0','First order','Second order','Third order'};
orderUnits = {'m','dimensionless','1/m','1/m^2'};

for order = 1:4
    channels = orderChannels{order};
    channels = channels(channels <= numChannels);
    if isempty(channels)
        continue;
    end

    figure('Name',['High-Order GSTF - ',orderNames{order}], ...
        'Units','normalized','Position',[0.12 0.12 0.60 0.68]);

    subplot(2,1,1);
    hold on;
    grid on;
    for channel = channels
        plot(H.f_axis/1000,abs(H.gstf(:,channel)), ...
            'DisplayName',termNames{channel});
    end
    xlim([-fLimit_kHz fLimit_kHz]);
    xlabel('Frequency (kHz)');
    ylabel(['Magnitude (',orderUnits{order},')']);
    title([upper(ax_name),' input - ',orderNames{order},' GSTF magnitude']);
    legend('Location','best');

    subplot(2,1,2);
    hold on;
    grid on;
    for channel = channels
        plot(H.f_axis/1000,angle(H.gstf(:,channel)), ...
            'DisplayName',termNames{channel});
    end
    xlim([-fLimit_kHz fLimit_kHz]);
    xlabel('Frequency (kHz)');
    ylabel('Phase (rad)');
    title([upper(ax_name),' input - ',orderNames{order},' GSTF phase']);
    legend('Location','best');
end

end
