function [ ] = graph2pdf_formatting(varargin)
%set different values of a graph


%% parse input
p = inputParser;

addParameter(p,'title', '');
addParameter(p,'xlabel', '');
addParameter(p,'ylabel', '');
addParameter(p,'legend', '');
addParameter(p,'fontSize', 16);
addParameter(p,'baseline',0);
    
parse(p, varargin{:});
%% determine objects to modifiy
current_figure = gcf;
all_axes = findall(current_figure, 'type', 'axes');

%add a zero line
if p.Results.baseline
    for i = 1:length(all_axes)
        refline(all_axes(i),0)
    end
end
    
%select a larger space of possible lane markings (color,shape)
set(0,'defaultAxesLineStyleOrder','-|--|:')

% title on top of graph
if ~isempty( p.Results.title )
    for i = 1:length(all_axes)
        title(all_axes(i),p.Results.title,'interpreter','tex')
    end
end

%label on x and y axis
if ~isempty( p.Results.xlabel )
    for i = 1:length(all_axes)
        xlabel(all_axes(i),p.Results.xlabel)
    end
end
if ~isempty( p.Results.ylabel )
    for i = 1:length(all_axes)
        ylabel(all_axes(i),p.Results.ylabel,'interpreter','tex')
    end
end

%legend entries for each lane
if ~isempty( p.Results.legend )
    for i = 1:length(all_axes)
        legend(all_axes(i),p.Results.legend )
    end
end

%set font size
for i = 1:length(all_axes)
    set(all_axes(i), 'FontSize', p.Results.fontSize)
end

%set xtick number [ ] to labels [ ]
%set(gca, 'XTick',1:numParameterType{param1}, 'XTickLabel',parameterTypes{param1})
%set(gca, 'XTick', []);
%set(gca, 'YTick', []);


