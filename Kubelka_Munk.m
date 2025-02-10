clc;
clear;
close all;
%%
data_date = "24.12.10";

[bandgap_table] = load_data(data_date);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bandgap_table] = load_data(data_date)
    bandgap_table = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% material = Sb2S3  -> material_inx = 1 %%%%
    %%%% material = Sb2Se3 -> material_inx = 3 %%%%
    %%%% material = Bi2S3  -> material_inx = 4 %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Specify the root folder
    alloy.rootFolder = "./"+data_date;

    % Get all items in the directory
    items = dir(alloy.rootFolder);
    
    % Filter to get directories only
    isFolder = [items.isdir]; % Logical array indicating directories
    folderNames = {items(isFolder).name}; % Get names of directories
    
    % Remove '.' and '..' (current and parent directory)
    folderNames = folderNames(~ismember(folderNames, {'.', '..'}));
    folderNames = string(folderNames);


    for i= 1:length(folderNames)
        alloy.index = i;
        % Specify the folder where the .txt files are located
        folderPath = alloy.rootFolder + "/" + folderNames(i);
        % Get a list of all .txt files in the folder
        filePattern = fullfile(folderPath, "*.txt"); % Pattern to match .txt files
        txtFiles = dir(filePattern);
        
        for j=1:length(txtFiles)
            
            alloy.name = extractAfter(folderNames(i),3) + ' - ' +...
                         extractBetween(txtFiles(j).name, '- ',' -');

            fileID = fopen(folderPath + "/" + string(txtFiles(j).name), 'r');
            for k = 1:18
                fgetl(fileID);
            end
            data = [];
    
            for k = 1:1201
                line = fgetl(fileID); % Read a line from the file
            
                data = [data; str2num(line)]; % Append to the data matrix
            end
            
            wavelength = data(:,1);
            reflectance = data(:,2);
            
            % plot the diffuse reflectance
            figure
            plot(wavelength, reflectance)
            title("Diffuse Reflectance - " + alloy.name)
            xlabel("Wavelenght (nm)")
            ylabel("Diffuse Reflectance (%)")
            xlim([620 1240])
            saveas(gca, alloy.rootFolder + "/DR img/" + ...
                   num2str(alloy.index) + ". " + alloy.name+".jpg")

            % call the function to calculate and plot K-M
            % [bandgap] = all_functions(wavelength, reflectance, alloy);
            % bandgap_table = [bandgap_table; alloy.name, bandgap];
            

        end

    end
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%% Run all of the functions %%%%%%%%%%%%%%%%%%%%%%%%%
function [bandgap] = all_functions(wavelength, reflectance, alloy)
    energy_eV = 1240./wavelength;
    reflectance_normalized = reflectance/100;

    F = (1 - reflectance_normalized) .^ 2 ./ (2 * reflectance_normalized);
    % to the power of 2, in the case of direct bandgap
    % to the power of 0.5, in the case of indirect bandgap
    A = (energy_eV .* F) .^ 2;
    
    % Find bandgap with finding the inflection point of sigmoid function after
    % smoothening
    [bandgap, sharpest_slope] = fit_findIP(A, energy_eV, alloy);
    
    % Plot the Kubelka-Munk function
    plot_KM(A, energy_eV, bandgap, sharpest_slope, alloy)
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%% Find Bandgap by fitting  %%%%%%%%%%%%%%%%%%%%%%%%%

function [bandgap, sharpest_slope] = fit_findIP(A, energy_eV, alloy)

    switch alloy.index
        case {6, 5, 7}
            MinEnergyFit = 1.2;
            MaxEnergyFit = 1.6;
        case {4, 3, 2}
            MinEnergyFit = 1.2;
            MaxEnergyFit = 1.8;
        case 1
            MinEnergyFit = 1.4;
            MaxEnergyFit = 1.8;
    end
    
    energy_range = (energy_eV > MinEnergyFit)& (energy_eV < MaxEnergyFit);

    energy_short = energy_eV(energy_range);
    A_short = A(energy_range);

    [f, ~] = fit(energy_short, A_short, 'gompertz');

    figure
    plot(energy_eV, A, energy_eV, f(energy_eV))
    xlim([1 2])
    legend("data", "fitted")
    saveas(gcf, alloy.rootFolder + "/K-M img/" + ...
           num2str(alloy.index) + ". " + alloy.name + " - Sigmoid.jpg");

    sharpest_grad = max( gradient( f(energy_short)));
    
    indx = find(gradient(f(energy_short))==sharpest_grad);
    
    which_d = 2;

    if which_d == 1 
        w = 20; % window gap
        sharpest_slope = (A_short(indx + w) - A_short(indx))/...
                            (energy_short(indx + w) - energy_short(indx));

        y1 = A_short(indx);
    else
        w = 1; % window gap
        sharpest_slope = (f(energy_short(indx + w)) - f(energy_short(indx)) )/...
                            (energy_short(indx + w) - energy_short(indx));
    
        y1 = f(energy_short(indx));
    end
    
    x1 = energy_short(indx);

    
    bandgap = - y1 / sharpest_slope + x1;

end

%%
%%%%%%%%%%%%%%%%%%%%% Plot the Kubelka-Munk function  %%%%%%%%%%%%%%%%%%%%%
function [] = plot_KM(A, energy_eV, bandgap, sharpest_slope, alloy)
    
    switch alloy.index
        case {3, 1}
            shift = 0.1;
        case {6, 5, 4, 2, 7}
            shift = 0.2;
    end
    
    figure
    plot(energy_eV, A, 'LineWidth',1.5)
    xlim([1 2])
    title("Kubelka-Munk Function - "+ alloy.name)
    xlabel("energy (eV)")
    ylabel("{(F(R_{\infty})h{\nu})}^2")

    hold on
    lin_x = linspace(bandgap, bandgap+shift, 20);
    lin_y = sharpest_slope * (lin_x - bandgap);

    plot(lin_x, lin_y, 'LineStyle','--', 'LineWidth', 2, 'Color', 'red' )

    legend(alloy.name, "Fitted Line")

    hold off

    saveas(gcf, alloy.rootFolder + "/K-M img/" + ...
           num2str(alloy.index) + ". " + alloy.name + " - Line.jpg");

end
%%