function result = bearing_thermo(params)
%BEARING_THERMO 滚动轴承低风险粘温等效参数模型。
% =========================
% 模块名称：滚动轴承热-粘等效参数模块
% 功能说明：计算 mu(T)、h(T)、Kb(T)、Cb(T) 与可选等效轴承力
% 输入输出：输入 params 结构，输出 result 结构
% 物理意义：球/滚子轴承兼容的降阶热粘参数映射
% =========================
%
% 本模型属于：
% Rolling bearing thermo-viscous equivalent parameter model
% (ball & roller bearing compatible reduced-order model)
%
% Thermo-viscous reduced-order rolling bearing model for ball and roller
% bearing systems; no fluid dynamic instability modeling included.
%
% 不涉及流体动压滑动轴承失稳机理，不包含 oil whirl、oil whip、
% 特征值稳定性判据、0.45 倍频判据或三次非线性刚度。

validateattributes(params, {'struct'}, {'scalar'}, mfilename, 'params');
require_fields(params, {'temp', 'T0', 'mu0', 'alpha', 'h_HD'});
temp = finite_scalar(params.temp, 'temp');
T0 = finite_scalar(params.T0, 'T0');
mu0 = positive_scalar(params.mu0, 'mu0');
alpha = nonnegative_scalar(params.alpha, 'alpha');
validateattributes(params.h_HD, {'numeric'}, ...
    {'real', 'finite', 'positive'}, mfilename, 'h_HD');
h_HD = params.h_HD;

temperature_mode = char(get_option(params, 'temperature_mode', 'inlet'));
if ~strcmpi(temperature_mode, 'inlet')
    error('bearing_thermo:UnsupportedTemperatureMode', ...
        'Only temperature_mode=''inlet'' is implemented in this phase.');
end

beta = positive_scalar(get_option(params, 'stiffness_exponent', 1), ...
    'stiffness_exponent');
n = finite_scalar(get_option(params, 'damping_exponent', 0.75), ...
    'damping_exponent');
if n < 0.3 || n > 1.0
    error('bearing_thermo:InvalidDampingExponent', ...
        'damping_exponent must be within [0.3, 1.0].');
end

% =========================
% 1. 粘度计算 mu(T)
% =========================
% 来源：Reid / Hamrock viscosity-temperature relation。
% 温度变化时缓存键同步变化，因此每个新温度状态都会重新计算。
mu_T = cached_andrade(temp, T0, mu0, alpha);
viscosity_ratio = mu_T / mu0;

% =========================
% 2. Reynolds 一致性说明
% =========================
% 本模块不求解完整 Reynolds PDE，仅保留温度经粘度影响润滑性能的
% 降阶映射，不声称获得完整压力场或流体动压稳定性结论。

% =========================
% 3. 油膜厚度 h(T)
% =========================
% 来源：Hamrock-Dowson EHL empirical relation。
h_T = h_HD .* viscosity_ratio.^0.68;

% =========================
% 4. 刚度 Kb(T)（EHL 降阶）
% =========================
% Kb 的严格定义为 dF_b/dq，且 F_b 为压力场面积积分。
% 当前表达以 h(T) 近似压力场变化，是 EHL 降阶等效映射。
% 来源：Hamrock-Dowson + Childs 线性化轴承动力学。
% 不代表真实压力场求导结果，也不是普适的直接反比物理定律。
Kb_T = [];
if isfield(params, 'Kb0') && ~isempty(params.Kb0)
    stiffness_scale = (h_HD ./ h_T).^beta;
    Kb_T = scale_equivalent(params.Kb0, stiffness_scale, 'Kb0');
end

% =========================
% 5. 阻尼 Cb(T)（Reynolds 扰动工程近似）
% =========================
% 阻尼为 Reynolds 扰动等效形式的工程近似，不构成完整 Reynolds PDE
% 求解结果；C0 必须由用户提供，禁止隐式硬编码经验常数。
Cb_T = [];
if isfield(params, 'C0') && ~isempty(params.C0)
    damping_scale = viscosity_ratio.^n;
    Cb_T = scale_equivalent(params.C0, damping_scale, 'C0');
end

% =========================
% 6. 轴承力 F_b
% =========================
% 按用户接口返回等效内力；系统平衡方程中的恢复力负号由调用方处理。
F_b = [];
has_q = isfield(params, 'delta_q') && ~isempty(params.delta_q);
has_qdot = isfield(params, 'delta_qdot') && ~isempty(params.delta_qdot);
if has_q || has_qdot
    if isempty(Kb_T) || isempty(Cb_T) || ~has_q || ~has_qdot
        error('bearing_thermo:IncompleteForceInput', ...
            'Kb0, C0, delta_q and delta_qdot are all required for F_b.');
    end
    F_b = apply_operator(Kb_T, params.delta_q, 'Kb_T', 'delta_q') + ...
        apply_operator(Cb_T, params.delta_qdot, 'Cb_T', 'delta_qdot');
end

[temp_next, Q_friction] = optional_thermal_update(params, temp, mu_T);

result = struct( ...
    'mu_T', mu_T, ...
    'h_T', h_T, ...
    'Kb_T', Kb_T, ...
    'Cb_T', Cb_T, ...
    'F_b', F_b, ...
    'temp_next', temp_next, ...
    'Q_friction', Q_friction, ...
    'temperature_mode', 'inlet', ...
    'model_name', ['Thermo-viscous reduced-order rolling bearing model ' ...
        'for ball and roller bearing systems']);
end

function require_fields(params, names)
for index = 1:numel(names)
    if ~isfield(params, names{index})
        error('bearing_thermo:MissingField', ...
            'Required field "%s" is missing.', names{index});
    end
end
end

function value = get_option(params, name, default_value)
if isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
else
    value = default_value;
end
end

function value = finite_scalar(value, name)
validateattributes(value, {'numeric'}, ...
    {'real', 'finite', 'scalar'}, mfilename, name);
end

function value = positive_scalar(value, name)
validateattributes(value, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, name);
end

function value = nonnegative_scalar(value, name)
validateattributes(value, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, name);
end

function mu_T = cached_andrade(temp, T0, mu0, alpha)
persistent cached_key cached_mu
key = [temp, T0, mu0, alpha];
if isempty(cached_key) || ~isequaln(key, cached_key)
    cached_mu = mu0 * exp(-alpha * (temp - T0));
    if ~isfinite(cached_mu) || cached_mu <= 0
        error('bearing_thermo:InvalidViscosity', ...
            ['Andrade relation produced a nonpositive or nonfinite ' ...
            'viscosity.']);
    end
    cached_key = key;
end
mu_T = cached_mu;
end

function scaled = scale_equivalent(base, scale, name)
validateattributes(base, {'numeric'}, {'real', 'finite'}, mfilename, name);
if strcmp(name, 'C0') && any(base(:) < 0)
    error('bearing_thermo:NegativeDamping', 'C0 must be nonnegative.');
end
if isscalar(scale)
    scaled = base .* scale;
elseif isequal(size(base), size(scale))
    scaled = base .* scale;
else
    error('bearing_thermo:ScaleDimensionMismatch', ...
        ['%s and its thermal scale must be scalar-compatible or ' ...
        'equal-sized.'], name);
end
end

function force = apply_operator(operator, state, operator_name, state_name)
validateattributes(state, {'numeric'}, {'real', 'finite', 'vector'}, ...
    mfilename, state_name);
state = state(:);
if isscalar(operator)
    force = operator .* state;
elseif isvector(operator)
    if numel(operator) ~= numel(state)
        error('bearing_thermo:ForceDimensionMismatch', ...
            '%s and %s must have equal lengths.', operator_name, state_name);
    end
    force = operator(:) .* state;
elseif ismatrix(operator) && size(operator, 1) == size(operator, 2)
    if size(operator, 2) ~= numel(state)
        error('bearing_thermo:ForceDimensionMismatch', ...
            '%s column count must equal the length of %s.', ...
            operator_name, state_name);
    end
    force = operator * state;
else
    error('bearing_thermo:InvalidOperator', ...
        '%s must be a scalar, vector, or square matrix.', operator_name);
end
end

function [temp_next, Q_friction] = optional_thermal_update(params, temp, mu_T)
temp_next = [];
Q_friction = [];
if ~isfield(params, 'thermal_update') || isempty(params.thermal_update)
    return;
end
cfg = params.thermal_update;
validateattributes(cfg, {'struct'}, {'scalar'}, mfilename, 'thermal_update');
if ~isfield(cfg, 'enabled') || ~cfg.enabled
    return;
end
required = {'dt', 'm_eff', 'cp', 'UA', 'T_ambient', 'omega', ...
    'loss_factor'};
require_fields(cfg, required);
dt = positive_scalar(cfg.dt, 'thermal_update.dt');
m_eff = positive_scalar(cfg.m_eff, 'thermal_update.m_eff');
cp = positive_scalar(cfg.cp, 'thermal_update.cp');
UA = nonnegative_scalar(cfg.UA, 'thermal_update.UA');
T_ambient = finite_scalar(cfg.T_ambient, 'thermal_update.T_ambient');
omega = finite_scalar(cfg.omega, 'thermal_update.omega');
loss_factor = nonnegative_scalar(cfg.loss_factor, ...
    'thermal_update.loss_factor');
if dt * UA / (m_eff * cp) > 1
    error('bearing_thermo:UnstableThermalStep', ...
        'Require dt*UA/(m_eff*cp) <= 1 for the explicit thermal update.');
end
% loss_factor 为有量纲标定系数，使 Q_friction 的单位为 W。
Q_friction = mu_T * omega^2 * loss_factor;
temp_next = temp + dt * ...
    (Q_friction - UA * (temp - T_ambient)) / (m_eff * cp);
if ~isfinite(temp_next)
    error('bearing_thermo:InvalidTemperatureUpdate', ...
        'The optional thermal update produced a nonfinite temperature.');
end
end
