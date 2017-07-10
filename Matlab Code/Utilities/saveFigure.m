function [] = saveFigure(folderName, figureName)
%Llamar a la función al final de la figura para tener el mismo formato en
%todas. Copiar el codigo comentado abajo para guardar la figura
%automaticamente en un archivo al que llamaremos desde al Lyx. Así, si se
%cambia algo en el matlab, el lyx lo substituira automaticamente por la
%nueva version.
%Añadidas modificaciones para combinar el formato con un guardado
%automático, se pasa el nombre de la figura como argumento.

     %Formating
     grafWidth   = 16;
     grafAR      = 0.6;
%      grid on
     set(gcf,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
     %Saving
     folderName = strrep(folderName,'.','');
     folderPath = [pwd filesep folderName];
     figureName = strrep(figureName,'.','');
     format_Grafico = [folderPath filesep figureName];
     exists=exist(folderPath,'dir');
     if isequal(exists,0)
         mkdir(folderPath);
         saveas(gcf,format_Grafico,'epsc'); 
     elseif isequal(exists,7)
         saveas(gcf,format_Grafico,'epsc');
     else
         disp('Imposible to save the Figure, check what can be failing.')
     end
     
end




    