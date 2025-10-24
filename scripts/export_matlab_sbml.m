function export_matlab_sbml(script_name, output_file)
%EXPORT_MATLAB_SBML Export a SimBiology tutorial model to SBML.
%
% export_matlab_sbml('example1', 'artifacts/sbml/example1.xml')

if nargin < 2
    error('Usage: export_matlab_sbml(script_name, output_file)');
end

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

output_dir = fileparts(output_file);
if ~isempty(output_dir) && ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

sbmlexport(model_ic, output_file);
end
