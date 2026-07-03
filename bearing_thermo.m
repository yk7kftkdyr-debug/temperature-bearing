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
%
% 当前工程接入状态（论文与结果解释必须保持一致）：
% ffLOAD/ff2 仅使用本模块的 mu_T 与 h_T，通过接触变形重新求解载荷，
% 等效刚度来自原接触模型的有限差分；显式 Kb_T、Cb_T、F_b 为预留
% 接口，当前主求解不激活阻尼修正，也不构成完整 THD 温度场模型。

validateattributes(params, {'struct'}, {'scalar'}, mfilename, 'params');
require_fields(params, {'T0', 'mu0', 'alpha', 'h_HD'});
T0 = finite_scalar(params.T0, 'T0');
mu0 = positive_scalar(params.mu0, 'mu0');
alpha = nonnegative_scalar(params.alpha, 'alpha');
validateattributes(params.h_HD, {'numeric'}, ...
    {'real', 'finite', 'positive'}, mfilename, 'h_HD');
h_HD = params.h_HD;

temperature_mode = lower(char(get_option(params, ...
    'temperature_mode', 'inlet')));
film_temperature_weight = bounded_unit_scalar(get_option(params, ...
    'film_temperature_weight', 0.5), 'film_temperature_weight');
switch temperature_mode
    case 'inlet'
        require_fields(params, {'temp'});
        temp_effective = finite_scalar(params.temp, 'temp');
    case 'effective_contact'
        require_fields(params, {'temp_outer', 'temp_inner'});
        temp_effective = [ ...
            finite_scalar(params.temp_outer, 'temp_outer'), ...
            finite_scalar(params.temp_inner, 'temp_inner')];
        if numel(h_HD) ~= 2
            error('bearing_thermo:InvalidFilmReferenceSize', ...
                ['effective_contact requires h_HD=[outer,inner] ' ...
                'with exactly two elements.']);
        end
        h_HD = reshape(h_HD, 1, 2);
    otherwise
        error('bearing_thermo:UnsupportedTemperatureMode', ...
            ['temperature_mode must be ''inlet'' or ' ...
            '''effective_contact''.']);
end

beta = positive_scalar(get_option(params, 'stiffness_exponent', 1), ...
    'stiffness_exponent');
n = finite_scalar(get_option(params, 'damping_exponent', 0.75), ...
    'damping_exponent');
if n < 0.3 || n > 1.0
    error('bearing_thermo:InvalidDampingExponent', ...
        'damping_exponent must be within [0.3, 1.0].');
end
film_viscosity_exponent = positive_scalar(get_option(params, ...
    'film_viscosity_exponent', 0.68), 'film_viscosity_exponent');
if film_viscosity_exponent > 1.0
    error('bearing_thermo:InvalidFilmViscosityExponent', ...
        'film_viscosity_exponent must be within (0, 1.0].');
end

% =========================
% 1. 粘度计算 mu(T)
% =========================
% 来源：Reid / Hamrock viscosity-temperature relation。
% 温度变化时缓存键同步变化，因此每个新温度状态都会重新计算。
mu_T = cached_andrade(temp_effective, T0, mu0, alpha);
viscosity_ratio = mu_T ./ mu0;

% =========================
% 2. Reynolds 一致性说明
% =========================
% 本模块不求解完整 Reynolds PDE，仅保留温度经粘度影响润滑性能的
% 降阶映射，不声称获得完整压力场或流体动压稳定性结论。

% =========================
% 3. 油膜厚度 h(T)
% =========================
% 来源：Hamrock-Dowson EHL empirical relation。
% 默认指数 0.68 用于点接触；滚子线接触可由调用方显式传入 0.70。
h_T = h_HD .* viscosity_ratio.^film_viscosity_exponent;

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

[temp_next, Q_shear, Q_friction, Q_gen, thermal_update_status] = ...
    optional_thermal_update(params, temp_effective, mu_T);
Cb_T_available = ~isempty(Cb_T);
Cb_T_used_in_equilibrium = false;
damping_model_active = false;
thermal_feedback_active = false;

result = struct( ...
    'mu_T', mu_T, ...
    'h_T', h_T, ...
    'Kb_T', Kb_T, ...
    'Cb_T', Cb_T, ...
    'F_b', F_b, ...
    'temp_next', temp_next, ...
    'Q_shear', Q_shear, ...
    'Q_friction', Q_friction, ...
    'Q_gen', Q_gen, ...
    'thermal_update_status', thermal_update_status, ...
    'film_viscosity_exponent', film_viscosity_exponent, ...
    'film_temperature_weight', film_temperature_weight, ...
    'temperature_mode', temperature_mode, ...
    'temp_effective', temp_effective, ...
    'temperature_fallback_used', false, ...
    'Cb_T_available', Cb_T_available, ...
    'Cb_T_used_in_equilibrium', Cb_T_used_in_equilibrium, ...
    'damping_model_active', damping_model_active, ...
    'thermal_feedback_active', thermal_feedback_active, ...
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

function value = bounded_unit_scalar(value, name)
value = finite_scalar(value, name);
if value < 0 || value > 1
    error('bearing_thermo:InvalidUnitIntervalValue', ...
        '%s must be within [0, 1].', name);
end
end

function mu_T = cached_andrade(temp, T0, mu0, alpha)
persistent cached_key cached_mu
key = [temp, T0, mu0, alpha];
if isempty(cached_key) || ~isequaln(key, cached_key)
    cached_mu = mu0 .* exp(-alpha .* (temp - T0));
    if any(~isfinite(cached_mu)) || any(cached_mu <= 0)
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
if isscalar(base) || isscalar(scale)
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

function [temp_next, Q_shear, Q_friction, Q_gen, status] = ...
        optional_thermal_update(params, temp_effective, mu_T)
temp_next = [];
Q_shear = [];
Q_friction = [];
Q_gen = [];
status = 'disabled';
if ~isfield(params, 'thermal_update') || isempty(params.thermal_update)
    return;
end
cfg = params.thermal_update;
validateattributes(cfg, {'struct'}, {'scalar'}, mfilename, 'thermal_update');
if ~isfield(cfg, 'enabled') || ~cfg.enabled
    return;
end
if ~isfield(cfg, 'Q_contact') || isempty(cfg.Q_contact) || ...
        ~isfield(cfg, 'v_slip') || isempty(cfg.v_slip)
    status = 'missing_dissipation_input';
    return;
end
required = {'dt', 'm_eff', 'cp', 'UA', 'T_ambient', 'omega', ...
    'xi_mu', 'eta_f'};
require_fields(cfg, required);
dt = positive_scalar(cfg.dt, 'thermal_update.dt');
m_eff = positive_scalar(cfg.m_eff, 'thermal_update.m_eff');
cp = positive_scalar(cfg.cp, 'thermal_update.cp');
UA = nonnegative_scalar(cfg.UA, 'thermal_update.UA');
T_ambient = finite_scalar(cfg.T_ambient, 'thermal_update.T_ambient');
omega = finite_scalar(cfg.omega, 'thermal_update.omega');
xi_mu = nonnegative_scalar(cfg.xi_mu, 'thermal_update.xi_mu');
eta_f = nonnegative_scalar(cfg.eta_f, 'thermal_update.eta_f');
validateattributes(cfg.Q_contact, {'numeric'}, ...
    {'real', 'finite', 'vector'}, mfilename, 'thermal_update.Q_contact');
validateattributes(cfg.v_slip, {'numeric'}, ...
    {'real', 'finite', 'vector'}, mfilename, 'thermal_update.v_slip');
if numel(cfg.Q_contact) ~= numel(cfg.v_slip)
    error('bearing_thermo:DissipationDimensionMismatch', ...
        'Q_contact and v_slip must have equal lengths.');
end
if dt * UA / (m_eff * cp) > 1
    error('bearing_thermo:UnstableThermalStep', ...
        'Require dt*UA/(m_eff*cp) <= 1 for the explicit thermal update.');
end
% Q_shear 是集总模型的标定型剪切损失接口，其中 xi_mu 负责量纲闭合；
% 它不等价于完整能量方程的局部剪切耗散项，也不包含速度梯度与空间积分。
Q_shear = sum(xi_mu .* mu_T .* omega.^2);
Q_friction = eta_f .* sum(abs(cfg.Q_contact(:) .* cfg.v_slip(:)));
Q_gen = Q_shear + Q_friction;
% effective_contact 使用内外圈温度均值作为单节点集总温度；这是工程聚合，
% 不表示已经求解内外圈空间温度场。
temp_reference = mean(temp_effective);
temp_next = temp_reference + dt * ...
    (Q_gen - UA * (temp_reference - T_ambient)) / (m_eff * cp);
if ~isfinite(temp_next)
    error('bearing_thermo:InvalidTemperatureUpdate', ...
        'The optional thermal update produced a nonfinite temperature.');
end
status = 'computed_open_loop';
end
