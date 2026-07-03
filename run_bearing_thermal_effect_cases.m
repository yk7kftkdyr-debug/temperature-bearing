clear; clc;

projectRoot = fileparts(mfilename('fullpath'));
desktopOut = fullfile('/Users/bailingguan/Desktop', '轴承热效应');
if exist(desktopOut, 'dir') ~= 7
    mkdir(desktopOut);
end

addpath(projectRoot);
addpath(fullfile(projectRoot, '球轴承程序'));
addpath(fullfile(projectRoot, '滚子轴承程序'));

% Table supplied by user:
% ambient, outer ring, inner ring, rolling element, shaft, housing.
thermalTable = [
    20,  40,  60,  55,  65,  35;
    50,  75, 105,  95, 110,  70;
    80, 115, 150, 138, 155, 105;
   100, 140, 180, 165, 185, 130];

rho_oil = 970; % kg/m^3, from uploaded thermal-effect document
reference_T = 20;
reference_mu = oil_dynamic_viscosity(reference_T, rho_oil);
viscosity_alpha = 0.02815;
records = repmat(empty_record(), 0, 1);

for row = 1:size(thermalTable, 1)
    ambient_T = thermalTable(row, 1);
    outer_T = thermalTable(row, 2);
    inner_T = thermalTable(row, 3);
    rolling_T = thermalTable(row, 4);
    shaft_T = thermalTable(row, 5);
    housing_T = thermalTable(row, 6);

    nu_T = oil_kinematic_viscosity(ambient_T);
    eta_T = oil_dynamic_viscosity(ambient_T, rho_oil);
    alpha_p_T = pressure_viscosity_coefficient(ambient_T);

    cfg = make_micro_interface_config();
    cfg.thermal.enabled = true;
    cfg.thermal.Ct1 = 1;
    cfg.thermal.Ct2 = 1;
    cfg.thermal.reference_temperature_C = reference_T;
    cfg.thermal.reference_viscosity_Pa_s = reference_mu;
    cfg.thermal.viscosity_temperature_coefficient = viscosity_alpha;

    tempDir = fullfile(desktopOut, sprintf('ambient_%03dC', ambient_T));
    ballOut = fullfile(tempDir, 'ball');
    rollerOut = fullfile(tempDir, 'roller');
    ensure_dir(ballOut);
    ensure_dir(rollerOut);

    ballData = default_ball_input();
    ballData(16) = 20000;
    ballData(17) = 0;
    ballData(18) = 20000;
    ballData(26) = rho_oil;
    ballData(27) = ambient_T;
    ballData(29) = eta_T;
    ballData(30) = alpha_p_T;
    % Ball program mapping: 46 rolling element, 47 outer, 48 inner,
    % 49 ambient, 59 shaft, 60 housing.
    ballData(46) = rolling_T;
    ballData(47) = outer_T;
    ballData(48) = inner_T;
    ballData(49) = ambient_T;
    ballData(59) = shaft_T;
    ballData(60) = housing_T;
    records(end + 1, 1) = run_ball_case(projectRoot, ballData, cfg, ballOut, ...
        ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, ...
        nu_T, eta_T, alpha_p_T); %#ok<SAGROW>

    rollerData = defaultBearingInput();
    rollerData(19) = 0;
    rollerData(20) = 80000;
    rollerData(26) = rho_oil;
    rollerData(27) = ambient_T;
    rollerData(29) = eta_T;
    rollerData(30) = alpha_p_T;
    % Roller program mapping: 44 roller, 45 outer, 46 inner,
    % 47 ambient, 57 shaft, 58 housing.
    rollerData(44) = rolling_T;
    rollerData(45) = outer_T;
    rollerData(46) = inner_T;
    rollerData(47) = ambient_T;
    rollerData(57) = shaft_T;
    rollerData(58) = housing_T;
    records(end + 1, 1) = run_roller_case(projectRoot, rollerData, cfg, rollerOut, ...
        ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, ...
        nu_T, eta_T, alpha_p_T); %#ok<SAGROW>

    writetable(struct2table(records), fullfile(desktopOut, 'bearing_thermal_effect_summary.csv'), 'Encoding', 'UTF-8');
end

writetable(struct2table(records), fullfile(desktopOut, 'bearing_thermal_effect_summary.csv'), 'Encoding', 'UTF-8');
disp('BEARING_THERMAL_EFFECT_CASES_DONE');
disp(desktopOut);

function nu = oil_kinematic_viscosity(T)
nu = 27.6 * exp(-0.02815 * (T - 40)); % cSt
end

function eta = oil_dynamic_viscosity(T, rho)
eta = rho * oil_kinematic_viscosity(T) * 1e-6; % Pa*s
end

function alpha_p = pressure_viscosity_coefficient(T)
alpha_p = 1.28e-8 * exp(-0.010 * (T - 20)); % Pa^-1
end

function ensure_dir(pathName)
if exist(pathName, 'dir') ~= 7
    mkdir(pathName);
end
end

function record = run_ball_case(projectRoot, datafromvb, cfg, outDir, ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T)
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(projectRoot, '球轴承程序'));
record = base_record('ball', ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T);
record.Fx_N = datafromvb(16);
record.Fy_N = datafromvb(17);
record.Fz_N = datafromvb(18);
clear_case_artifacts("ball");
lastwarn('');
try
    runLog = evalc('[loadj, returndata] = qiujieend(datafromvb, cfg);');
    [record.warning_flag, record.comment] = classify_warning(runLog);
    record = fill_metrics(record, "ball", loadj, returndata);
    record = validate_record(record);
    copyfile(fullfile(projectRoot, 'results_bearing_base', 'ball', '*'), outDir);
catch ME
    record.warning_flag = 'error';
    record.comment = sprintf('%s at %s:%d', ME.message, error_file(ME), error_line(ME));
end
end

function record = run_roller_case(projectRoot, datafromvb, cfg, outDir, ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T)
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir));
cd(fullfile(projectRoot, '滚子轴承程序'));
record = base_record('roller', ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T);
record.Fx_N = NaN;
record.Fy_N = datafromvb(19);
record.Fz_N = datafromvb(20);
clear_case_artifacts("roller");
lastwarn('');
try
    runLog = evalc('[loadj, returndata] = qiujieall(datafromvb, cfg);');
    [record.warning_flag, record.comment] = classify_warning(runLog);
    record = fill_metrics(record, "roller", loadj, returndata);
    record = validate_record(record);
    copyfile(fullfile(projectRoot, 'results_bearing_base', 'roller', '*'), outDir);
catch ME
    record.warning_flag = 'error';
    record.comment = sprintf('%s at %s:%d', ME.message, error_file(ME), error_line(ME));
end
end

function record = base_record(bearingType, ambient_T, outer_T, inner_T, rolling_T, shaft_T, housing_T, nu_T, eta_T, alpha_p_T)
record = empty_record();
record.bearing_type = bearingType;
record.ambient_temperature_C = ambient_T;
record.outer_ring_temperature_C = outer_T;
record.inner_ring_temperature_C = inner_T;
record.rolling_element_temperature_C = rolling_T;
record.shaft_temperature_C = shaft_T;
record.housing_temperature_C = housing_T;
record.kinematic_viscosity_cSt = nu_T;
record.dynamic_viscosity_Pa_s = eta_T;
record.pressure_viscosity_Pa_inv = alpha_p_T;
end

function [warningFlag, comment] = classify_warning(runLog)
[lastMsg, lastId] = lastwarn;
warningFlag = 'none';
comment = '运行完成';
hasWarningText = contains(runLog, 'Warning', 'IgnoreCase', true) || contains(runLog, '警告');
hasSingular = contains(runLog, 'singular', 'IgnoreCase', true) || contains(lastMsg, 'singular', 'IgnoreCase', true) || ...
    contains(runLog, '奇异') || contains(lastMsg, '奇异');
if hasSingular
    warningFlag = 'singular_matrix';
    comment = '运行完成，但出现奇异矩阵相关警告';
elseif hasWarningText || strlength(string(lastMsg)) > 0 || strlength(string(lastId)) > 0
    warningFlag = 'warning';
    if strlength(string(lastMsg)) > 0
        comment = char(lastMsg);
    else
        comment = '运行完成，但命令输出包含 warning';
    end
end
end

function record = fill_metrics(record, bearingType, loadj, returndata)
record.loaded_element_count = loadj;
if bearingType == "ball"
    contact = load('q1q2a1a2.mat', 'q1q2a1a2');
    contact = contact.q1q2a1a2;
    record.max_contact_load = max(contact(:, 1:2), [], 'all');
    record.final_residual = NaN;
    residualLimit = NaN;
    record.residual_status = 'not_available_after_stiffness_perturbation';
else
    record.max_contact_load = max_named_from_mats( ...
        {'Q1.mat', 'Q2.mat'}, {'Q1', 'Q2'});
    record.final_residual = scalar_from_mat('result111.mat', 'result111');
    residualLimit = 0.1;
    record.residual_status = 'available';
end
record.max_contact_stress = max_named_from_mats( ...
    {'Ph1.mat', 'Ph2.mat'}, {'Ph1', 'Ph2'});
record.min_oil_film_thickness = min_named_from_mats( ...
    {'oilh1.mat', 'oilh2.mat'}, {'oilh1', 'oilh2'});
if bearingType == "ball"
    record.PV_value = max_named_from_mats( ...
        {'pvzhi1.mat', 'pvzhi2.mat', 'pvzhi1nonload.mat'}, ...
        {'pvzhi1', 'pvzhi2', 'pvzhi1nonload'});
    record.cage_slip_ratio = value_or_nan(returndata, 6);
    record.equivalent_stiffness_x = value_or_nan(returndata, 1);
    record.equivalent_stiffness_y = value_or_nan(returndata, 2);
    record.life_or_damage_index = value_or_nan(returndata, 9);
else
    record.PV_value = max_named_from_mats( ...
        {'pvzhi1.mat', 'pvzhi2.mat', 'pvzhinonload.mat'}, ...
        {'pvzhi1', 'pvzhi2', 'pvzhinonload'});
    record.cage_slip_ratio = value_or_nan(returndata, 5);
    record.equivalent_stiffness_x = value_or_nan(returndata, 1);
    record.equivalent_stiffness_y = value_or_nan(returndata, 2);
    record.life_or_damage_index = value_or_nan(returndata, 6);
end
record.residual_limit = residualLimit;
record.stiffness_source = 'contact_finite_difference';
record.damping_model_active = false;
end

function value = max_named_from_mats(fileNames, variableNames)
values = read_named_values(fileNames, variableNames);
value = max(values, [], 'omitnan');
end

function value = min_named_from_mats(fileNames, variableNames)
values = read_named_values(fileNames, variableNames);
value = min(values, [], 'omitnan');
end

function values = read_named_values(fileNames, variableNames)
assert(numel(fileNames) == numel(variableNames), ...
    'ThermalSummary:VariableMapping', ...
    'MAT file and variable-name lists must have equal length.');
values = [];
for i = 1:numel(fileNames)
    if exist(fileNames{i}, 'file') ~= 2
        continue;
    end
    data = load(fileNames{i}, variableNames{i});
    if ~isfield(data, variableNames{i})
        continue;
    end
    candidate = data.(variableNames{i});
    if isnumeric(candidate) && ~isempty(candidate)
        values = [values; candidate(:)]; %#ok<AGROW>
    end
end
if isempty(values)
value = NaN;
values = value;
end
end

function record = validate_record(record)
criticalValues = [record.max_contact_load, record.max_contact_stress, ...
    record.min_oil_film_thickness, record.PV_value, ...
    record.cage_slip_ratio, record.equivalent_stiffness_x, ...
    record.equivalent_stiffness_y, record.life_or_damage_index];
residualValid = ~isfinite(record.residual_limit) || ...
    (isfinite(record.final_residual) && ...
    record.final_residual <= record.residual_limit);
record.valid_result = all(isfinite(criticalValues)) && ...
    record.min_oil_film_thickness > 0 && residualValid;
if record.valid_result && strcmp(record.warning_flag, 'singular_matrix')
    record.comment = ['completed_with_warning：最终关键输出有限且油膜为正；' ...
        '迭代过程中出现奇异矩阵相关警告'];
elseif ~record.valid_result && ~strcmp(record.warning_flag, 'error')
    record.warning_flag = 'invalid_result';
    record.comment = '结果包含NaN、Inf或非正油膜厚度';
end
end

function value = scalar_from_mat(fileName, variableName)
data = load(fileName, variableName);
value = data.(variableName);
value = value(1);
end

function clear_case_artifacts(bearingType)
common = {'Ph1.mat','Ph2.mat','oilh1.mat','oilh2.mat', ...
    'pvzhi1.mat','pvzhi2.mat'};
if bearingType == "ball"
    files = [common, {'q1q2a1a2.mat','pvzhi1nonload.mat','result333.mat'}];
else
    files = [common, {'Q1.mat','Q2.mat','pvzhinonload.mat','result111.mat'}];
end
for i = 1:numel(files)
    if exist(files{i}, 'file') == 2
        delete(files{i});
    end
end
end

function value = value_or_nan(values, index)
if numel(values) >= index
    value = values(index);
else
    value = NaN;
end
end

function record = empty_record()
record = struct( ...
    'bearing_type', '', ...
    'ambient_temperature_C', NaN, ...
    'outer_ring_temperature_C', NaN, ...
    'inner_ring_temperature_C', NaN, ...
    'rolling_element_temperature_C', NaN, ...
    'shaft_temperature_C', NaN, ...
    'housing_temperature_C', NaN, ...
    'kinematic_viscosity_cSt', NaN, ...
    'dynamic_viscosity_Pa_s', NaN, ...
    'pressure_viscosity_Pa_inv', NaN, ...
    'Fx_N', NaN, ...
    'Fy_N', NaN, ...
    'Fz_N', NaN, ...
    'loaded_element_count', NaN, ...
    'max_contact_load', NaN, ...
    'max_contact_stress', NaN, ...
    'min_oil_film_thickness', NaN, ...
    'PV_value', NaN, ...
    'cage_slip_ratio', NaN, ...
    'equivalent_stiffness_x', NaN, ...
    'equivalent_stiffness_y', NaN, ...
    'life_or_damage_index', NaN, ...
    'final_residual', NaN, ...
    'residual_limit', NaN, ...
    'residual_status', '', ...
    'stiffness_source', '', ...
    'damping_model_active', false, ...
    'valid_result', false, ...
    'warning_flag', '', ...
    'comment', '');
end

function fileName = error_file(ME)
if isempty(ME.stack)
    fileName = 'unknown';
else
    fileName = ME.stack(1).file;
end
end

function lineNo = error_line(ME)
if isempty(ME.stack)
    lineNo = 0;
else
    lineNo = ME.stack(1).line;
end
end

function datafromvb = default_ball_input()
datafromvb = zeros(75, 1);
datafromvb(1)=15; datafromvb(2)=0.02223; datafromvb(3)=0.1253; datafromvb(4)=0.5232; datafromvb(5)=0.5232; datafromvb(6)=40; datafromvb(7)=0.0115; datafromvb(35)=5;
datafromvb(61)=7800; datafromvb(68)=7800; datafromvb(51)=7800; datafromvb(50)=7800; datafromvb(9)=2.06e+11; datafromvb(10)=2.06e+11; datafromvb(11)=2.06e+11;
datafromvb(65)=2.06e+11; datafromvb(12)=0.3; datafromvb(13)=0.3; datafromvb(14)=0.3; datafromvb(66)=0.3; datafromvb(69)=0; datafromvb(15)=10000; datafromvb(16)=20000; datafromvb(17)=0;
datafromvb(18)=3000; datafromvb(19)=0; datafromvb(20)=0; datafromvb(21)=0.275; datafromvb(22)=0.225; datafromvb(23)=1e-07; datafromvb(24)=1e-07; datafromvb(25)=1e-07; datafromvb(26)=970;
datafromvb(27)=20; datafromvb(28)=0.0966; datafromvb(29)=0.0318; datafromvb(30)=1.28e-008; datafromvb(31)=3.2e-002; datafromvb(32)=1; datafromvb(33)=0.4e-03; datafromvb(34)=3.8;
datafromvb(74)=2; datafromvb(36)=370e-03; datafromvb(37)=220e-03; datafromvb(38)=120.0e-03; datafromvb(39)=225e-03; datafromvb(40)=1.96e11; datafromvb(41)=2.18e11; datafromvb(42)=0.3; datafromvb(43)=0.3;
datafromvb(44)=0.5e-03; datafromvb(45)=-0.5e-03; datafromvb(54)=11.6e-06; datafromvb(55)=11.8e-06; datafromvb(56)=11.8e-06; datafromvb(57)=11.8e-06; datafromvb(58)=11.8e-06; datafromvb(46)=180;
datafromvb(47)=170; datafromvb(48)=190; datafromvb(49)=27; datafromvb(59)=195; datafromvb(60)=160; datafromvb(50)=7870; datafromvb(51)=7870; datafromvb(52)=7860; datafromvb(53)=8360;
datafromvb(70)=0; datafromvb(71)=0.2; datafromvb(72)=0.2; datafromvb(73)=0.2;
end