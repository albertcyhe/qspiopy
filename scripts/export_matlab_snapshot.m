function export_matlab_snapshot(script_name, output_dir)
%EXPORT_MATLAB_SNAPSHOT Freeze a SimBiology tutorial model to disk.
%
% export_matlab_snapshot('example1', 'artifacts/matlab_frozen_model')
%
% Generates a directory named after the script (example1, example2, ...)
% containing:
%   - odes.txt              : getequations(model) output
%   - fluxes.txt            : reaction flux listing
%   - species.csv           : species metadata and initial amounts
%   - parameters.csv        : parameter values
%   - compartments.csv      : compartment volumes
%   - reactions.csv         : stoichiometry-aware reaction metadata
%   - stoichiometry.csv     : long-form stoichiometry matrix
%   - rules.csv             : assignment/rate rules
%   - events.csv            : model events

if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('artifacts', 'matlab_frozen_model');
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

target_root = fullfile(repo_root, output_dir, script_name);
if ~exist(target_root, 'dir')
    mkdir(target_root);
end

eq_text = string(getequations(model_ic));
write_text(fullfile(target_root, 'equations.txt'), eq_text);

config = getconfigset(model_ic);

export_configset(config, fullfile(target_root, 'configset.json'));
export_species(model_ic, fullfile(target_root, 'species.csv'), config);
export_parameters(model_ic, fullfile(target_root, 'parameters.csv'));
export_compartments(model_ic, fullfile(target_root, 'compartments.csv'));
export_reactions(model_ic, fullfile(target_root, 'reactions.csv'));
export_stoichiometry(model_ic, fullfile(target_root, 'stoichiometry.csv'));
export_reaction_parameters(model_ic, fullfile(target_root, 'reaction_parameters.csv'));
export_rules(model_ic, fullfile(target_root, 'rules.csv'));
export_events(model_ic, fullfile(target_root, 'events.csv'));
export_doses(model_ic, fullfile(target_root, 'doses.csv'));
export_variants(model_ic, fullfile(target_root, 'variants.csv'));

end

function export_configset(config, destination)
info = struct();
info.SolverType = char(config.SolverType);
info.StopTime = config.StopTime;
info.TimeUnits = char(config.TimeUnits);

solver_opts = config.SolverOptions;
info.SolverOptions = struct( ...
    'AbsoluteTolerance', solver_opts.AbsoluteTolerance, ...
    'RelativeTolerance', solver_opts.RelativeTolerance, ...
    'MaxStep', solver_opts.MaxStep, ...
    'InitialStep', solver_opts.InitialStep, ...
    'MaxOrder', solver_opts.MaxOrder ...
);

compile_opts = config.CompileOptions;
info.CompileOptions = struct( ...
    'DimensionalAnalysis', compile_opts.DimensionalAnalysis, ...
    'UnitConversion', compile_opts.UnitConversion, ...
    'DefaultSpeciesDimension', char(compile_opts.DefaultSpeciesDimension) ...
);

write_json(destination, info);
end


function export_species(model, destination, config)
species = model.Species;
n = numel(species);
data = cell(n, 11);
default_dim = char(config.CompileOptions.DefaultSpeciesDimension);
for i = 1:n
    s = species(i);
    compartment = char(s.Parent.Name);
    data{i, 1} = sprintf('%s.%s', compartment, s.Name);
    data{i, 2} = s.Name;
    data{i, 3} = compartment;
    data{i, 4} = s.InitialAmount;
    data{i, 5} = char(s.InitialAmountUnits);
    data{i, 6} = s.ConstantAmount;
    data{i, 7} = s.Value;
    data{i, 8} = s.BoundaryCondition;
    if isprop(s, 'NonNegative')
        data{i, 9} = s.NonNegative;
    else
        data{i, 9} = false;
    end
    units = char(s.Units);
    data{i, 10} = units;
    if isempty(units)
        data{i, 11} = default_dim;
    elseif contains(units, '/') || contains(units, 'per')
        data{i, 11} = 'concentration';
    else
        data{i, 11} = 'amount';
    end
end
headers = {'identifier', 'name', 'compartment', 'initial_amount', 'initial_units', 'constant_amount', 'current_value', 'boundary_condition', 'nonnegative', 'units', 'interpreted_dimension'};
write_cell_table(headers, data, destination);
end

function export_parameters(model, destination)
params = model.Parameters;
n = numel(params);
data = cell(n, 5);
for i = 1:n
    p = params(i);
    data{i, 1} = p.Name;
    data{i, 2} = p.Value;
    data{i, 3} = char(p.ValueUnits);
    data{i, 4} = p.ConstantValue;
    if isprop(p, 'InitialValue')
        data{i, 5} = p.InitialValue;
    else
        data{i, 5} = '';
    end
end
headers = {'name', 'value', 'units', 'constant_value', 'initial_value'};
write_cell_table(headers, data, destination);
end

function export_compartments(model, destination)
comps = model.Compartments;
n = numel(comps);
data = cell(n, 4);
for i = 1:n
    c = comps(i);
    data{i, 1} = c.Name;
    data{i, 2} = c.Capacity;
    data{i, 3} = char(c.CapacityUnits);
    data{i, 4} = c.ConstantCapacity;
end
headers = {'name', 'capacity', 'units', 'constant'};
write_cell_table(headers, data, destination);
end

function export_reactions(model, destination)
rxns = model.Reactions;
n = numel(rxns);
data = cell(n, 6);
for i = 1:n
    r = rxns(i);
    data{i, 1} = r.Name;
    data{i, 2} = strjoin(arrayfun(@(x) x.Name, r.Reactants, 'UniformOutput', false), ';');
    data{i, 3} = strjoin(arrayfun(@(x) x.Name, r.Products, 'UniformOutput', false), ';');
    data{i, 4} = char(r.ReactionRate);
    if isempty(r.KineticLaw)
        data{i, 5} = '';
        data{i, 6} = '';
    else
        data{i, 5} = r.KineticLaw.KineticLawName;
        data{i, 6} = char(r.KineticLaw.Expression);
    end
end
headers = {'name', 'reactants', 'products', 'reaction_rate', 'kinetic_law', 'kinetic_expression'};
write_cell_table(headers, data, destination);
end

function export_stoichiometry(model, destination)
[S, species, reactions] = getstoichmatrix(model);
if iscell(species)
    species_names = species;
else
    species_names = arrayfun(@(s) s.Name, species, 'UniformOutput', false);
end
if iscell(reactions)
    reaction_names = reactions;
else
    reaction_names = arrayfun(@(r) r.Name, reactions, 'UniformOutput', false);
end
rows = {};
for i = 1:size(S, 1)
    for j = 1:size(S, 2)
        value = S(i, j);
        if value ~= 0
            rows(end + 1, :) = {species_names{i}, reaction_names{j}, value}; %#ok<AGROW>
        end
    end
end
headers = {'species', 'reaction', 'stoichiometry'};
write_cell_table(headers, rows, destination);
end

function export_reaction_parameters(model, destination)
rxns = model.Reactions;
rows = {};
for i = 1:numel(rxns)
    r = rxns(i);
    if isempty(r.KineticLaw)
        continue;
    end
    params = r.KineticLaw.Parameters;
    for j = 1:numel(params)
        p = params(j);
        rows(end + 1, :) = {r.Name, p.Name, p.Value, char(p.ValueUnits)}; %#ok<AGROW>
    end
end
headers = {'reaction', 'name', 'value', 'units'};
write_cell_table(headers, rows, destination);
end

function export_rules(model, destination)
rules = model.Rules;
n = numel(rules);
data = cell(n, 5);
for i = 1:n
    rule = rules(i);
    rule_type = char(rule.RuleType);
    expr = char(rule.Rule);
    target = '';
    if contains(expr, '=')
        parts = strsplit(expr, '=');
        target = strtrim(parts{1});
        expr = strtrim(strjoin(parts(2:end), '='));
    elseif isprop(rule, 'Variable')
        target = char(rule.Variable);
    end
data{i, 1} = i;
data{i, 2} = rule.Name;
data{i, 3} = rule_type;
data{i, 4} = target;
data{i, 5} = expr;
end
headers = {'rule_index', 'name', 'type', 'target', 'expression'};
write_cell_table(headers, data, destination);
end

function export_events(model, destination)
events = model.Events;
n = numel(events);
data = cell(n, 6);
for i = 1:n
    event = events(i);
    data{i, 1} = i;
    data{i, 2} = event.Name;
    data{i, 3} = char(event.Trigger);
    if isprop(event, 'DelayValue')
        data{i, 4} = event.DelayValue;
    else
        data{i, 4} = 0;
    end
    if isprop(event, 'DelayType')
        data{i, 5} = char(event.DelayType);
    else
        data{i, 5} = '';
    end
    fcns = event.EventFcns;
    if iscell(fcns)
        data{i, 6} = strjoin(string(fcns), ';');
    else
        data{i, 6} = char(fcns);
    end
end
headers = {'event_index', 'name', 'trigger', 'delay', 'delay_type', 'assignments'};
write_cell_table(headers, data, destination);
end

function export_doses(model, destination)
doses = model.Doses;
n = numel(doses);
data = cell(n, 9);
for i = 1:n
    dose = doses(i);
    data{i, 1} = i;
    data{i, 2} = dose.Name;
    data{i, 3} = dose.Type;
    data{i, 4} = char(dose.TargetName);
    data{i, 5} = dose.Amount;
    data{i, 6} = char(dose.AmountUnits);
    data{i, 7} = dose.StartTime;
    data{i, 8} = dose.Interval;
    data{i, 9} = dose.RepeatCount;
end
headers = {'dose_index', 'name', 'type', 'target', 'amount', 'amount_units', 'start_time', 'interval', 'repeat_count'};
write_cell_table(headers, data, destination);
end

function export_variants(model, destination)
variants = model.Variants;
n = numel(variants);
rows = cell(n, 4);
for i = 1:n
    variant = variants(i);
    rows{i, 1} = i;
    rows{i, 2} = variant.Name;
    rows{i, 3} = variant.Active;
    rows{i, 4} = format_variant_content(variant.Content);
end
headers = {'variant_index', 'name', 'active', 'content'};
write_cell_table(headers, rows, destination);
end

function text = format_variant_content(content)
if isempty(content)
    text = '';
    return;
end
pairs = cellfun(@(entry) strjoin(string(entry), ':'), content, 'UniformOutput', false);
text = strjoin(pairs, ';');
end

function write_text(path, content)
fid = fopen(path, 'w');
if fid == -1
    error('Unable to write %s', path);
end
clean = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strtrim(content));
end

function write_json(path, data)
fid = fopen(path, 'w');
if fid == -1
    error('Unable to write %s', path);
end
clean = onCleanup(@() fclose(fid));
json_text = jsonencode(data, 'PrettyPrint', true);
fprintf(fid, '%s', json_text);
end

function write_cell_table(headers, data, destination)
fid = fopen(destination, 'w');
if fid == -1
    error('Unable to write %s', destination);
end
clean = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(headers, ','));
for i = 1:size(data, 1)
    row = cellfun(@sanitize_field, data(i, :), 'UniformOutput', false);
    fprintf(fid, '%s\n', strjoin(row, ','));
end
end

function out = sanitize_field(value)
if isnumeric(value) && isscalar(value) && ~isnan(value)
    out = num2str(value, '%.15g');
elseif islogical(value)
    out = char(string(value));
elseif isempty(value)
    out = '';
else
    out = char(strrep(string(value), ',', ';'));
end
out = strrep(out, sprintf('\n'), ' ');
out = strrep(out, sprintf('\r'), ' ');
end
