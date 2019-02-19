function analyzer_GUI(input)
global fig_save plt
clearvars('-except', 'input', 'fig_save', 'plt');   close all;

try; properties.file = input.file; catch; error('input file not specified'); end
try; properties.plotFigure =input.plotFigure; catch; properties.plotFigure = false; end
try; properties.saveFigure = input.saveFigure; catch;properties.saveFigure = false; end

%global variable to control plotting of the fits during
%the characterization process
plt = properties.plotFigure;
if properties.saveFigure; fig_save = 'png'; end

 
fol = [pwd '\'];
try;
  
    msgbox({'Parameter characterization in progress...'});
    characterizeDocument([fol properties.file],@parameterCharacterization2);
    close all;
    msgbox({'Parameter characterization succesful!'});
catch;
end
end