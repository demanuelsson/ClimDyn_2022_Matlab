

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change to correct minus symbol
ax=gca;
set(ax,'ticklabelinterpreter','none')  %or 'tex' but not 'latex'
yticklabels(ax, strrep(yticklabels(ax),'-','–'));
xticklabels(ax, strrep(xticklabels(ax),'-','–'));
%cb = colorbar();
h.TickLabels = strrep(h.TickLabels, '-', '–');
%%%%%%%%%%%%%%%%%%%%%%%%%%%  


% −