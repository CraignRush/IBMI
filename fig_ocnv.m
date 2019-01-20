%% Matlab figure Conversion

%% Debugging: Open Fotoresist Plot
%Photoresist_eval;
i= 1;
fig = get(groot,'CurrentFigure');
while ~isempty(fig)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% P A R A M E T E R S %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ENABLE_XTIC_LABEL_ROTATION = false;
    ENABLE_YTIC_LABEL_ROTATION = false;
    
    ENABLE_FIG_BACKGROUND = true;
    
    ENABLE_FIG_SAVING = true;
    
    FIGURENAME = num2str(i);
    i = i +1;
    FORMAT = 'epsc'; % svg, pdf
    
    
    %% Get handles
    fig = gcf;
    ax = gca;
    
    if ~isgraphics(fig) || ~isgraphics(ax)
        error('Please select the figure!');
    end
    
    
    %% Set bg-transparent
    if ENABLE_FIG_BACKGROUND
        fig.Color = [1,1,1];
    else
        fig.Color = 'none';
    end

    
    
    %% Save figure in folder

        %folder = uigetdir('Figures\','Please specify saving location');
        folder = 'Figures_HW2\';
        if strcmp(FORMAT,'epsc')
            
            path = strcat(folder,'\', FIGURENAME, '.', 'eps');
        end
        if strcmp(FORMAT,'pdf')
            print('-painters',fig,path,['-d' FORMAT],'-bestfit');
        else
            print('-painters',fig,path,['-d' FORMAT]);
        end
    close(fig);
    fig  = get(groot,'CurrentFigure');
end