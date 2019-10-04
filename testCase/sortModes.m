%% Sort fields
% Author: Marc Rovira

%% Preamble

clear; close all; clc;

%% Sort outputed modes

flag_U = 1; % 0: scalar, 1: velocity

%N = 5; % Number of modes
N = 200;

coord = 'xyz';

mkdir output/sortedModes
mkdir output/sortedModes/s
mkdir output/sortedModes/U

f = waitbar(0, 'Sorting fields...');

for i = 1:N

    orig = dlmread('pointCloud_reference.xy'); % x y z and scalar from boxUniform output

    if flag_U == 0% if Scalar
        mode = dlmread(strcat('output/mode/mode_s_', num2str(i - 1), '.dat'));

        orig(:, end) = mode; % substitute 4th column (scalar) for disordered scalar

        A = sortrows(orig, [3 1 2]); % Sort orig such that: 1st:y 2nd:x 3rd:z

        fileID = fopen(strcat('output/sortedModes/s/mode_s_', num2str(i - 1), '.dat'), 'w');
        fprintf(fileID, '%5.6e\n', A(1:end - 1, end)');
        fprintf(fileID, '%5.6e\n', A(end, end)');
        fclose(fileID);

    elseif flag_U == 1% if Velocity

        for j = 1:3

            mode = dlmread(strcat('output/mode/mode_U', coord(j), '_', num2str(i - 1), '.dat'));
            %mode = dlmread(strcat('output/data/U', coord(j), '_', num2str(i - 1), '.dat'));

            orig(:, end) = mode; % substitute 4th-6th columns for velocity

            A = sortrows(orig, [3 1 2]); % Sort orig such that: 1st:y 2nd:x 3rd:z

            fileID = fopen(strcat('output/sortedModes/U/mode_U', coord(j), '_', num2str(i - 1), '.dat'), 'w');
            fprintf(fileID, '%5.6e\n', A(1:end - 1, end)');
            fprintf(fileID, '%5.6e', A(end, end));
            fclose(fileID);

            waitbar(i / N, f, 'Sorting fields...');

        end

    end

end

close(f)