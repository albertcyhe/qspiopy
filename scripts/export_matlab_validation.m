function export_matlab_validation(output_dir)
%EXPORT_MATLAB_VALIDATION Run SimBiology scenarios and export trajectories.
%
% This helper aligns the MATLAB tutorial models with the Python validation
% harness by emitting CSV artefacts for Example 1 (control and anti-PD1) and
% Example 2 (anti-PD1).  The output schema mirrors the Python reference
% results so downstream scripts can compute error metrics directly.

if nargin < 1 || isempty(output_dir)
    output_dir = 'artifacts/matlab_validation';
end

repo_root = fileparts(fileparts(mfilename('fullpath')));
scripts_dir = fullfile(repo_root, 'scripts');

addpath(genpath(fullfile(repo_root, 'model')));
addpath(genpath(fullfile(repo_root, 'utils')));
addpath(scripts_dir);

scenarios = { ...
    struct('name', 'example1_control', 'script', 'example1.m', 'therapy', 'none', 'horizon', 400), ...
    struct('name', 'example1_treated', 'script', 'example1.m', 'therapy', 'anti_pd1', 'horizon', 400), ...
    struct('name', 'example2_treated', 'script', 'example2.m', 'therapy', 'anti_pd1', 'horizon', 400) ...
    };
benchmark_replicates = 25;

output_path = fullfile(repo_root, output_dir);
if ~exist(output_path, 'dir')
    mkdir(output_path);
end

for k = 1:numel(scenarios)
    scenario = scenarios{k};
    sim_data = run_scenario(repo_root, scripts_dir, scenario);
    frame = simdata_to_table(sim_data);
    writetable(frame, fullfile(output_path, [scenario.name, '_matlab.csv']));
end

benchmark_seconds = benchmark_matlab_reference(repo_root, scripts_dir, benchmark_replicates);

performance_path = fullfile(output_path, 'performance.json');
payload = struct( ...
    'replicates', benchmark_replicates, ...
    'matlab_reference_seconds', benchmark_seconds ...
    );
json_text = jsonencode(payload, 'PrettyPrint', true);
fid = fopen(performance_path, 'w');
if fid == -1
    error('Unable to open %s for writing.', performance_path);
end
clean = onCleanup(@() fclose(fid));
fwrite(fid, json_text);

end

function sim_data = run_scenario(repo_root, scripts_dir, scenario)
prev_dir = pwd;
cleanup = onCleanup(@() cd(prev_dir));
cd(repo_root);

script_name = erase(scenario.script, '.m');
if isempty(which(script_name))
    addpath(scripts_dir);
end
eval(script_name);

[model_ic, success] = initial_conditions(model);
if ~success
    error('Initial conditions failed to converge for %s', scenario.name);
end

config = getconfigset(model_ic);
times = 0:1:scenario.horizon;
set(config, 'StopTime', scenario.horizon);
set(config.SolverOptions, 'OutputTimes', times);

switch scenario.therapy
    case 'none'
        sim_data = sbiosimulate(model_ic);
    case 'anti_pd1'
        sim_data = sbiosimulate(model_ic, dose_schedule);
    otherwise
        error('Unsupported therapy: %s', scenario.therapy);
end
end

function frame = simdata_to_table(sim_data)
times = sim_data.Time(:);

cancer_cells = get_series(sim_data, 'C_total');
dead_cells = get_series(sim_data, 'C_x', 'V_T');
t_cells = tumour_t_cells(sim_data);
volume_ul = get_series(sim_data, 'V_T');
volume_l = volume_ul * 1e-6;

diameter_cm = 2.0 * ((3.0 * (volume_l * 1e3)) ./ (4.0 * pi)).^(1.0 / 3.0);
pd1 = get_series(sim_data, 'H_PD1_C1');
density = t_cells ./ max(volume_l * 1e6, eps);

frame = table(times, cancer_cells, dead_cells, t_cells, volume_l, diameter_cm, pd1, density, ...
    'VariableNames', {'time_days', 'cancer_cells', 'dead_cells', 't_cells', ...
    'tumour_volume_l', 'tumour_diameter_cm', 'pd1_occupancy', 'tcell_density_per_ul'});
end

function series = get_series(sim_data, name, varargin)
indices = find(strcmp(sim_data.DataNames, name));
if isempty(indices)
    error('Quantity %s not found in SimData object.', name);
end

if ~isempty(varargin)
    compartment = varargin{1};
    keep = false(size(indices));
    for i = 1:numel(indices)
        info = sim_data.DataInfo{indices(i)};
        if isfield(info, 'Compartment') && strcmp(info.Compartment, compartment)
            keep(i) = true;
        end
    end
    indices = indices(keep);
end

if isempty(indices)
    error('Quantity %s not available in the requested compartment.', name);
end

series = sum(sim_data.Data(:, indices), 2);
end

function series = tumour_t_cells(sim_data)
series = zeros(size(sim_data.Time(:)));
for i = 1:numel(sim_data.DataNames)
    name = sim_data.DataNames{i};
    if startsWith(name, 'T') && numel(name) >= 2 && all(isstrprop(name(2:end), 'digit')) && name(2) ~= '0'
        info = sim_data.DataInfo{i};
        if isfield(info, 'Compartment') && strcmp(info.Compartment, 'V_T')
            series = series + sim_data.Data(:, i);
        end
    end
end
end

function seconds = benchmark_matlab_reference(repo_root, scripts_dir, replicates)
prev_dir = pwd;
cleanup = onCleanup(@() cd(prev_dir));
cd(repo_root);

script_name = 'example1';
if isempty(which(script_name))
    addpath(scripts_dir);
end
eval(script_name);

[model_ic, success] = initial_conditions(model);
if ~success
    error('Initial conditions failed during benchmark setup.');
end

config = getconfigset(model_ic);
times = 0:1:400;
set(config, 'StopTime', 400);
set(config.SolverOptions, 'OutputTimes', times);

% warm-up run
sbiosimulate(model_ic, dose_schedule);

tic;
for i = 1:replicates
    sbiosimulate(model_ic, dose_schedule);
end
seconds = toc;
end
