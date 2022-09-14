% ----------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>
% Copyright (C) 2016, 2017    Nathan Schmidt, Gevorg Grigoryan
% ----------------------------------------------------------------------------
%
% FCOILSCAN Local Crick parameter analysis of an input structure, and
% accommodation index calculation.
%
% [C, P, AIP] = fcoilscan(FILE, TYPE, NC, WIN, PH, ODIR, XL)
%
% This function calculates the local Crick parameters for an input structure by
% fitting the Crick coiled coil equations at every residue position using a
% window size of WIN residues. Fitting is done by using the FCRICK function. If used
% for research, please cite:
%
% [1] G. Grigoryan, W. F. DeGrado, "Probing Designability via a Generalized Model of 
% Helical Bundle Geometry", J. Mol. Biol., 405(4): 1079-1100, 2011.
% [2] N. W. Schmidt, G. Grigoryan, W. F. DeGrado, "The accommodation index measures
% the perturbation associated with insertions and deletions in coiled-coils: Application
% to understand signaling in histidine kinases", Protein Science, 2017.
%
% ====== Inputs:
% FILE - input coiled-coil coordinates. Either a PDB file with multiple chains OR a flat
%        text file with X, Y, Z coordinates of CA atoms. In the latter case, must have
%        one line per atom, with each line having only three floating point values,
%        corresponding to X, Y, and Z-coordinates. Which input type is used is specified
%        by the following parameter.
% TYPE - input type, either 0 or 1. The value of 1 indicates that FILE is a PDB file, while
%        0 means that FILE is an XYZ file.
% NC   - number of chains in the input structure. Whether the input is specified as a PDB
%        or an XYZ file, all CA atoms will be combined into a single array, in the order
%        read, and will be split into NC equal-length chains. This means that the number of
%        CA atoms in the input must be an integer multiple of NC. If it is not, an error
%        will be thrown.
% WIN  - the length of the local window, in residues, over which local Crick parameter
%        fits will be performed (7 is recommended). Thus, a total of NC - WIN + 1 fits will
%        be performed.
% PH   - the ideal minorhelical phase change (in degrees) between neighboring CA positions.
%        For a 7/2 canonical coiled coil this value should be typically set to 360/3.5.
% ODIR - output directory where all results will be placed. This directory will be created
%        if it does not already exist.
% XL   - optional: if given and evaluates to true, will attempt to write results in an Excel
%        file using xlswrite (may not always work, depending on the operating system and
%        Matlab vs. Octave). The default behavior is not to write this file.
%
% ====== Outputs:
% C    - a cell array which lists chains, residue numbers and provides the AI profile fit
%        parameters and their errors.
% P    - a cell array of the local-window Crick parameters. The rows correspond to a residue
%        position (the starting position of the local window) and the columns correspond to
%        fit Crick parameter and fitting error (in Angstroms). 
% AIP  - lists minorhelical expected phase, observed phase of the fits, the phase difference
%        between expected and observed phases, and finally the AI at every residue position,
%        for each chain.
%
% ====== ODIR directory
% fcoilscan also creates a folder containing an Excel spreadsheet of fit data (if xlOut was
% set, see above), and two pdfs with plots of local Crick parameters and AI/phase profiles.
% The Excel file has three spreadsheets: 1) CoilFits, 2) Parameters, and 3) AI profiles as,
% which contain the data described above in C, P, and AIP, respectively. Also included in the
% ODIR directory are plots of the local Crick parameters: R0, R1, w0, w1, alpha, rise/residue,
% and plots of the AI profile and phase difference plot for each chain in the coiled coil.

function [coilfits, parameters, AIphases] = fcoilscan(file, coorType, Nchains, win, Ph_ideal, outDir, xlOut)
    if (~exist('xlOut', 'var')), xlOut = 0; end

    % In Octave, load the octave packages optim for leasqr and io for xlswrite.
    progV = ver;
    if (any(strcmpi('octave', {progV.Name})))
      pkg load optim
      if (xlOut), pkg load io; end
      graphics_toolkit gnuplot % in case the system does not support the default OpenGL
    end


    % cread CA coordinates from either a PDB file or a flat matrix of
    % pre-parsed coordinates
    M = readCA(file, coorType);

    % Make directory for the structure data
    if (~exist(outDir, 'dir')), mkdir(outDir); end

    % Check to make sure M is reasonable - must be divisible by Nchains and
    % contain no rows of zeros
    if mod(size(M, 1), Nchains) ~= 0
        error('Total number of bundle coords is not divisible by number of chains!')
    end

    % Create a folder for the PDB file. Put text files in the folder where
    % each file contains the CA coordinates corresponding to win residue
    % (usually 7) pieces of the coiled-coil for a rolling window along the
    % coil. So if coiled-coil is N residues I will have N-win+1 files -> for a
    % 11 residue coiled-coil split into 7 residue pieces I will have 5 files
    % containing residues 1-7, 2-8, 3-9, 4-10, & 5-11. Then fit each of the
    % pieces with fcrick.

    mkdir(outDir);
    N = size(M, 1)/Nchains;
    output = {};
    for a = 1:(N - (win-1))
        disp(sprintf('\nfitting window %d/%d :REPORT', a, (N - (win-1))));
        Calphas = [];
        for m = 1:Nchains
            Calphas = [Calphas; M(N*(m - 1) + a:N*(m - 1) + a + win - 1, :)];
        end
        [err XYZ pret] = fcrick(Calphas, Nchains, 'GENERAL-HLXPH', 2, 0, [], [], [], []);
        outpt = squeeze(struct2cell(pret));
        output(a,:) = [outpt(1,:) num2cell(err)];
    end
    disp(sprintf('\nextracting accommodation index data...:REPORT'));

    position = [1:N-win+1]';
    output = [['Starting Position' outpt(2,:) 'rmsd error']; ...
              [num2cell(position) output]];

    %Calculate the phase differences and accommodation indices
    Phase_ex = bsxfun(@plus, ((2:length(output(:,1))) - 2)'*Ph_ideal, ...
       cell2mat(output(2, find(strncmp(output(1,:),'ph1',3))))*180/pi);
    Phase_obs = [cell2mat(output(2:length(output(:,1)), ...
              find(strncmp(output(1,:),'ph1',3))))]*180/pi;
    % Calculate AI, cumulative phases must be calculated to get full AI, i.e.
    % in case the phase difference is such that it differs by an amount
    % more than ~180 degrees
    len_ph = length(Phase_obs(:,1))
    CumPhase_obs = cumsum([Phase_ex(1,:);
    mod(atan2(cosd(Phase_obs(1:len_ph-1,:)), sind(Phase_obs(1:len_ph-1,:)))*180/pi ...
    -  atan2(cosd(Phase_obs(2:len_ph,:)), sind(Phase_obs(2:len_ph,:)))*180/pi, 360)]);

    delPhase = Phase_ex - CumPhase_obs;
    AI = delPhase/Ph_ideal;
    % Next write the phase values to lie between [0 360] degrees.
    Phase_ex = mod(Phase_ex, 360);
    Phase_obs = mod(Phase_obs, 360);

    PhaseLabels = strcat({'Phase Expected ', 'Phase Observed ', ...
                     'Phase Difference ', 'AI '}, num2str(1));
    PhaseValues = [Phase_ex(:,1) Phase_obs(:,1) delPhase(:,1) AI(:,1)];

    for i = 2:Nchains
        PhaseLabels = [PhaseLabels; strcat({'Phase Expected ', 'Phase Observed ', ...
            'Phase Difference ', 'AI '}, num2str(i))];
        PhaseValues = [PhaseValues Phase_ex(:,i) Phase_obs(:,i) delPhase(:,i) AI(:,i)];
    end
    Phases = [reshape(PhaseLabels', [1, 4*Nchains]); num2cell(PhaseValues)];

    % Fit accommodation index curve, and calculate standard deviation of the
    % parameters
    options = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'MaxFunEvals', 100000, 'MaxIter', 1000, 'Jacobian', 'on');
    LOW = [-0.6 -0.1 1 0.01];
    UPP = [1.20 0.1 N (N + 1)/2];
    Std = cell(size(AI, 2), length(LOW));
    for i = 1:size(AI, 2)
        [xfit(i,:), ~, RES, ~, ~, ~, J] = lsqcurvefit(@erfnorm, [0 0 (N + 1)/2 10], position, AI(:,i), [LOW], [UPP], options);
        rows = size(AI, 1);
        cols = length(xfit(i, :));
        % A posteriori variance factor
        COVP = full(J'*J); % parameter covariance matrix
        if (rank(COVP) < size(COVP, 1))
            Std(i, :) = repmat({'underdefined'}, 1, length(cols));
        else
          Sigma_o = sqrt(sum(RES.^2)/max(rows-cols, eps)); % if there are more parameters than windows, precision will be infinity
          Precision_of_solved_parameters = Sigma_o * sqrt(diag(inv(COVP)));
          Std(i, :) = num2cell(Precision_of_solved_parameters);
        end
    end
    disp(sprintf('\npreparing output...:REPORT'));

    % Generate excel spreadsheet file with best fit Crick parameters, the
    % calculated AI profile and phase difference data, and the AI profile
    % fit parameters and errors.
    header = {'Chain' 'First Position' 'Second Position' ...
            'Insertion Index' 'Background' 'Midpoint Position' ...
            'Gaussian Accommodation, sigma (Residues)' ...
            'Insertion Index Error' 'Background Error' ...
            'Midpoint Position Error' 'Gaussian Accommodation, sigma, Error'};
    xlFile = [outDir '/CrickParam.xlsx'];
    if (xlOut), xlswrite(xlFile, header, 'CoilFits', strcat('A', num2str(1), ':', 'L', num2str(1))); end 
    for i = 1:Nchains
        fitdata(i, :) = {i 1 N xfit(i,1) xfit(i,2) xfit(i,3) xfit(i,4) ...
            Std{i, 1} Std{i, 2} Std{i, 3} Std{i, 4}};
        if (xlOut), xlswrite(xlFile, fitdata(i, :), 'CoilFits', strcat('A', num2str(i + 1), ':', 'L', num2str(i + 1))); end
    end
    if (xlOut)
        xlswrite(xlFile, output, 'Parameters', 'A1');
        xlswrite(xlFile, Phases, 'AIprofiles', 'A1');
    end

    % Finally generate Figures. Make a pdf consisting of a [3 by 2] array of
    % Crick parameter plots as a function of position along the coil.
    % Make a pdf consisting of a [2 by 1] array of the calculated AI and
    % associated fits for each chain in the coil on top, and their associated
    % phase changes on bottom.
    figure('visible', 'off');
    subplot(3,2,1)
    plot(position, cell2mat(output(2:length(output(:,1)), ...
                    find(strncmp(output(1,:),'R0',2)))),'b-o', 'linewidth', 2);
    xlabel('Position', 'fontweight', 'bold');
    ylabel('R_0 (A)', 'fontweight', 'bold');

    subplot(3,2,2)
    plot(position, cell2mat(output(2:length(output(:,1)), ...
                    find(strncmp(output(1,:),'R1',2)))),'b-o', 'linewidth', 2);
    xlabel('Position', 'fontweight', 'bold');
    ylabel('R_1 (A)', 'fontweight', 'bold');

    subplot(3,2,3)
    plot(position, 180/pi*cell2mat(output(2:length(output(:,1)), ...
                    find(strncmp(output(1,:),'w0',2)))),'b-o', 'linewidth', 2);
    xlabel('Position', 'fontweight', 'bold');
    ylabel('w_0 (degrees)', 'fontweight', 'bold');

    subplot(3,2,4)
    plot(position, 180/pi*cell2mat(output(2:length(output(:,1)), ...
                    find(strncmp(output(1,:),'w1',2)))),'b-o', 'linewidth', 2);
    xlabel('Position', 'fontweight', 'bold');
    ylabel('w_1 (degrees)', 'fontweight', 'bold');

    subplot(3,2,5)
    plot(position, 180/pi*cell2mat(output(2:length(output(:,1)), ...
                    find(strncmp(output(1,:),'alpha',5)))),'b-o', 'linewidth', 2);
    xlabel('Position', 'fontweight', 'bold');
    ylabel('Alpha (degrees)', 'fontweight', 'bold');

    subplot(3,2,6)
    plot(position, cell2mat(output(2:length(output(:,1)), ...
                    find(strncmp(output(1,:),'rise',4)))),'b-o', 'linewidth', 2);
    xlabel('Position', 'fontweight', 'bold');
    ylabel('Rise/Residue (A)', 'fontweight', 'bold');

    fig1 = gcf;
    set(gcf, 'paperunits', 'inches')
    set(gcf, 'paperposition', [0.5 0.5 7.5 10])
    saveas(fig1, [outDir, '/CrickParam_plots'], 'pdf')
    print(fig1, [outDir, '/CrickParam_plots.jpg'], '-djpeg', '-r300')
    print(fig1, [outDir, '/CrickParam_plots.small.jpg'], '-djpeg', '-r72')

    figure('visible', 'off');
    subplot(2,1,1)
    plot(position,AI,'b-o');
    for i = 1:size(AI, 2)
        plot(position, AI(:, i), 'ro', position, erfnorm(xfit(i,:), position), 'k-', 'linewidth', 2);
    end
    xlabel('Position', 'fontweight', 'bold');
    ylabel('AI', 'fontweight', 'bold');

    subplot(2,1,2)
    plot(position, delPhase,'b-o', 'linewidth', 2);
    xlabel('Position', 'fontweight', 'bold');
    ylabel('Phase Difference (degrees)', 'fontweight', 'bold');

    fig2 = gcf;
    set(gcf, 'paperunits', 'inches')
    set(gcf, 'paperposition', [0.5 0.5 7.5 10])
    saveas(fig2, [outDir, '/AI_profiles'], 'pdf')
    print(fig2, [outDir, '/AI_profiles.jpg'], '-djpeg', '-r300')
    print(fig2, [outDir, '/AI_profiles.small.jpg'], '-djpeg', '-r72')

    coilfits = {header; fitdata};
    parameters = output;
    AIphases = Phases;
    close all
    disp(sprintf('\nfinished:REPORT'));

end

% ERFNORM fits the AI profile using an erf(x) 'error' function corresponding
% to the integral of the normalized Gaussian function from -inf to the point x.
function [F, varargout] = erfnorm(x, xdata)
    % The insertion index; Iindex  = 0, & 0.5, 1, -0.5 for canonical coils, and
    % ones containing 4, 1, and 3 residue insertions, respectively.
    Iindex = x(1);
    % Simple constant background.
    Back = x(2);
    % Midpoint of accommodation region where AI profile has reached its halfway
    % point. Defines the structural position of the insertion.
    mu = x(3);
    % Gaussian width, which is proportional to the accommodation length La by
    % La = 2*sqrt(2)*sig
    sig = x(4);

    F = (Iindex/2)*(1 + erf((xdata - mu)/(sig*sqrt(2)))) + Back;

    % now the Jacobian
    if (nargout > 1)
      if (size(xdata, 2) ~= 1), xdata = xdata'; end
      J(:, 1) = (1/2)*(1 + erf((xdata - mu)/(sig*sqrt(2))));
      J(:, 2) = 1;
      J(:, 3) = -(Iindex/sqrt(2*pi()))*exp(-(xdata - mu).^2/(2*sig^2))*1/sig;
      J(:, 4) = J(:, 3).*(xdata - mu)/sig;
      varargout{1} = J;
     end
end

