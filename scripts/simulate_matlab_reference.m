function simulate_matlab_reference(script_name, therapy, days, output_path)
%SIMULATE_MATLAB_REFERENCE Run a tutorial scenario and persist the trajectories.
%
% script_name - MATLAB script configuring the model (example1 / example2)
% therapy     - 'none' or 'anti_pd1'
% days        - simulation horizon in days
% output_path - destination CSV path

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'model')));
addpath(genpath(fullfile(repo_root, 'utils')));
addpath(fullfile(repo_root, 'scripts'));
addpath(repo_root);

prev_dir = pwd;
cleanup = onCleanup(@() cd(prev_dir));
cd(repo_root);

run(fullfile(repo_root, 'scripts', script_name + ".m"));
[model_ic, success] = initial_conditions(model);
if ~success
    error('Initial conditions failed for %s', script_name);
end

config = getconfigset(model_ic);
times = 0:1:days;
set(config, 'StopTime', days);
set(config.SolverOptions, 'OutputTimes', times);

switch therapy
    case "none"
        simData = sbiosimulate(model_ic);
    case "anti_pd1"
        if ~exist('dose_schedule', 'var')
            error('Dose schedule unavailable in script %s for therapy export.', script_name);
        end
        simData = sbiosimulate(model_ic, dose_schedule);
    otherwise
        error('Unsupported therapy option: %s', therapy);
end

frame = simdata_to_table(simData);
output_dir = fileparts(output_path);
if ~isempty(output_dir) && ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
writetable(frame, output_path);
end

function frame = simdata_to_table(simData)
times = simData.Time(:);

cancer_cells = get_series(simData, 'C_total');
dead_cells = get_series(simData, 'C_x', 'V_T');
t_cells = tumour_t_cells(simData);
volume_ul = get_series(simData, 'V_T');
volume_l = volume_ul * 1e-6;

diameter_cm = 2.0 * ((3.0 * (volume_l * 1e3)) ./ (4.0 * pi)).^(1.0 / 3.0);
pd1 = get_series(simData, 'H_PD1_C1');
density = t_cells ./ max(volume_l * 1e6, eps);

frame = table(times, cancer_cells, dead_cells, t_cells, volume_l, diameter_cm, pd1, density, ...
    'VariableNames', {'time_days', 'cancer_cells', 'dead_cells', 't_cells', ...
    'tumour_volume_l', 'tumour_diameter_cm', 'pd1_occupancy', 'tcell_density_per_ul'});
end

function series = get_series(simData, name, varargin)
indices = strcmp(simData.DataNames, name);
if ~any(indices)
    error('Quantity %s not found in SimData object.', name);
end

if ~isempty(varargin)
    compartment = varargin{1};
    keep = false(size(indices));
    for i = 1:numel(indices)
        info = simData.DataInfo{i};
        if isfield(info, 'Compartment') && strcmp(info.Compartment, compartment)
            keep(i) = true;
        end
    end
    indices = indices & keep;
end

if ~any(indices)
    error('Quantity %s not available in the requested compartment.', name);
end

series = sum(simData.Data(:, indices), 2);
end

function series = tumour_t_cells(simData)
series = zeros(size(simData.Time(:)));
for i = 1:numel(simData.DataNames)
    name = simData.DataNames{i};
    if startsWith(name, 'T') && numel(name) >= 2 && all(isstrprop(name(2:end), 'digit')) && name(2) ~= '0'
        info = simData.DataInfo{i};
        if isfield(info, 'Compartment') && strcmp(info.Compartment, 'V_T')
            series = series + simData.Data(:, i);
        end
    end
end
end
