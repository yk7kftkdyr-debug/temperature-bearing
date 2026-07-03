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

contact_result = contact_diagnostics(params, mu_T);
if contact_result.present
    Cb_T = contact_result.Cb_T_total;
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

[temp_next, Q_shear, Q_friction, Q_gen, tau_th, ...
    Q_shear_contact, Q_friction_contact, thermal_update_status] = ...
    optional_thermal_update(params, contact_result);
Cb_T_available = ~isempty(Cb_T);
Cb_T_used_in_equilibrium = false;
damping_model_active = contact_result.damping_model_active;
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
    'tau_th', tau_th, ...
    'Q_shear_contact', Q_shear_contact, ...
    'Q_friction_contact', Q_friction_contact, ...
    'thermal_update_status', thermal_update_status, ...
    'film_viscosity_exponent', film_viscosity_exponent, ...
    'film_temperature_weight', film_temperature_weight, ...
    'temperature_mode', temperature_mode, ...
    'temp_effective', temp_effective, ...
    'temperature_fallback_used', false, ...
    'Cb_T_available', Cb_T_available, ...
    'Cb_T_used_in_equilibrium', Cb_T_used_in_equilibrium, ...
    'damping_model_active', damping_model_active, ...
    'pressure_valid', contact_result.pressure_valid, ...
    'damping_valid', contact_result.damping_valid, ...
    'heat_valid', contact_result.heat_valid, ...
    'p_eq_contact', contact_result.p_eq_contact, ...
    'Cb_T_contact', contact_result.Cb_T_contact, ...
    'Cb_T_total', contact_result.Cb_T_total, ...
    'contact_diagnostics_status', contact_result.status, ...
    'thermal_feedback_active', thermal_feedback_active, ...
    'model_name', ['Thermo-viscous reduced-order rolling bearing model ' ...
        'for ball and roller bearing systems']);
end

function output = contact_diagnostics(params, mu_T)
output = struct('present',false,'pressure_valid',[], ...
    'damping_valid',[],'heat_valid',[],'p_eq_contact',[], ...
    'Cb_T_contact',[],'Cb_T_total',[], ...
    'damping_model_active',false,'status','not_requested', ...
    'Q_contact',[],'A_eff',[],'h_T',[],'U',[],'v_slip',[], ...
    'mu_contact',[]);
if ~isfield(params,'contact') || isempty(params.contact)
    return;
end
output.present = true;
contact = params.contact;
required = {'Q_contact','A_eff','h_T','U','v_slip'};
for index = 1:numel(required)
    if ~isfield(contact,required{index}) || isempty(contact.(required{index}))
        output.status = 'missing_contact_input';
        return;
    end
end
Q_contact = real_matrix(contact.Q_contact,'contact.Q_contact');
A_eff = real_matrix(contact.A_eff,'contact.A_eff');
h_contact = real_matrix(contact.h_T,'contact.h_T');
U = real_matrix(contact.U,'contact.U');
v_slip = real_matrix(contact.v_slip,'contact.v_slip');
expected_size = size(Q_contact);
if size(Q_contact,2) ~= 2 || ~isequal(size(A_eff),expected_size) || ...
        ~isequal(size(h_contact),expected_size) || ...
        ~isequal(size(U),expected_size) || ...
        ~isequal(size(v_slip),expected_size)
    error('bearing_thermo:ContactDimensionMismatch', ...
        'All contact arrays must be equal-sized N-by-2 matrices.');
end
mu_nodes = two_node_value(mu_T,'mu_T');
mu_contact = repmat(mu_nodes,size(Q_contact,1),1);
pressure_valid = Q_contact>0 & A_eff>0 & ...
    isfinite(Q_contact) & isfinite(A_eff);
damping_valid = pressure_valid & h_contact>0 & ...
    isfinite(h_contact) & isfinite(mu_contact);
heat_valid = damping_valid & isfinite(U) & isfinite(v_slip);
p_eq_contact = NaN(expected_size);
p_eq_contact(pressure_valid) = ...
    Q_contact(pressure_valid)./A_eff(pressure_valid);

Cb_T_contact = [];
Cb_T_total = [];
damping_active = false;
if isfield(params,'damping_calibration') && ...
        ~isempty(params.damping_calibration)
    gamma_c = nonnegative_scalar(params.damping_calibration, ...
        'damping_calibration');
    Cb_T_contact = NaN(expected_size);
    Cb_T_contact(damping_valid) = gamma_c .* ...
        mu_contact(damping_valid) .* A_eff(damping_valid).^2 ./ ...
        h_contact(damping_valid).^3;
    if any(~isfinite(Cb_T_contact(damping_valid)))
        error('bearing_thermo:InvalidEquivalentDamping', ...
            'Equivalent damping is nonfinite.');
    end
    if any(damping_valid(:))
        Cb_T_total = sum(Cb_T_contact,1,'omitnan');
        damping_active = true;
    end
end
output = struct('present',true,'pressure_valid',pressure_valid, ...
    'damping_valid',damping_valid,'heat_valid',heat_valid, ...
    'p_eq_contact',p_eq_contact,'Cb_T_contact',Cb_T_contact, ...
    'Cb_T_total',Cb_T_total,'damping_model_active',damping_active, ...
    'status','computed','Q_contact',Q_contact,'A_eff',A_eff, ...
    'h_T',h_contact,'U',U,'v_slip',v_slip, ...
    'mu_contact',mu_contact);
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

function value = real_matrix(value, name)
validateattributes(value, {'numeric'}, {'real','2d'}, mfilename, name);
end

function value = two_node_value(value, name)
validateattributes(value, {'numeric'}, ...
    {'real','finite','vector'}, mfilename, name);
if isscalar(value)
    value = [value,value];
elseif numel(value) == 2
    value = reshape(value,1,2);
else
    error('bearing_thermo:InvalidTwoNodeValue', ...
        '%s must be a scalar or two-element vector.',name);
end
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

function [temp_next, Q_shear, Q_friction, Q_gen, tau_th, ...
        Q_shear_contact, Q_friction_contact, status] = ...
        optional_thermal_update(params, contact)
temp_next = [];
Q_shear = [];
Q_friction = [];
Q_gen = [];
tau_th = [];
Q_shear_contact = [];
Q_friction_contact = [];
status = 'disabled';
if ~isfield(params, 'thermal_update') || isempty(params.thermal_update)
    return;
end
cfg = params.thermal_update;
validateattributes(cfg, {'struct'}, {'scalar'}, mfilename, 'thermal_update');
if ~isfield(cfg, 'enabled') || ~cfg.enabled
    return;
end
if ~contact.present || ~strcmp(contact.status,'computed')
    status = 'missing_contact_input';
    return;
end
required = {'eta_s','eta_f','m_eff','cp','UA','dt', ...
    'T_ambient','temp_state'};
for index = 1:numel(required)
    if ~isfield(cfg,required{index}) || isempty(cfg.(required{index}))
        status = 'missing_calibration_input';
        return;
    end
end
dt = positive_scalar(cfg.dt, 'thermal_update.dt');
m_eff = positive_two_node(cfg.m_eff,'thermal_update.m_eff');
cp = positive_two_node(cfg.cp,'thermal_update.cp');
UA = nonnegative_two_node(cfg.UA,'thermal_update.UA');
T_ambient = two_node_value(cfg.T_ambient,'thermal_update.T_ambient');
temp_state = two_node_value(cfg.temp_state,'thermal_update.temp_state');
eta_s = nonnegative_scalar(cfg.eta_s, 'thermal_update.eta_s');
eta_f = nonnegative_scalar(cfg.eta_f, 'thermal_update.eta_f');
relax = bounded_positive_unit_scalar(get_option(cfg,'relax',0.3), ...
    'thermal_update.relax');
if any(~any(contact.heat_valid,1))
    status = 'missing_contact_heat_input';
    return;
end
Q_shear_contact = NaN(size(contact.Q_contact));
Q_friction_contact = NaN(size(contact.Q_contact));
valid = contact.heat_valid;
Q_shear_contact(valid) = eta_s .* contact.mu_contact(valid) .* ...
    contact.U(valid).^2 ./ contact.h_T(valid) .* contact.A_eff(valid);
Q_friction_contact(valid) = eta_f .* ...
    abs(contact.Q_contact(valid).*contact.v_slip(valid));
Q_shear = sum(Q_shear_contact,1,'omitnan');
Q_friction = sum(Q_friction_contact,1,'omitnan');
Q_gen = Q_shear + Q_friction;
tau_th = Inf(1,2);
positive_UA = UA>0;
tau_th(positive_UA) = ...
    m_eff(positive_UA).*cp(positive_UA)./UA(positive_UA);
if any(dt>0.2.*tau_th(isfinite(tau_th)))
    warning('bearing_thermo:LargeThermalTimeStep', ...
        'dt exceeds 0.2 times a finite thermal time constant.');
end
T_raw = temp_state + dt.* ...
    (Q_gen-UA.*(temp_state-T_ambient))./(m_eff.*cp);
temp_next = relax.*T_raw + (1-relax).*temp_state;
has_min = isfield(cfg,'T_min') && ~isempty(cfg.T_min);
has_max = isfield(cfg,'T_max') && ~isempty(cfg.T_max);
if xor(has_min,has_max)
    error('bearing_thermo:IncompleteTemperatureBounds', ...
        'T_min and T_max must be provided together.');
end
if has_min
    T_min = two_node_value(cfg.T_min,'thermal_update.T_min');
    T_max = two_node_value(cfg.T_max,'thermal_update.T_max');
    if any(T_min>=T_max)
        error('bearing_thermo:InvalidTemperatureBounds', ...
            'T_min must be less than T_max at both nodes.');
    end
    if any(temp_next<T_min | temp_next>T_max)
        temp_next = [];
        status = 'temperature_out_of_range';
        return;
    end
end
if any(~isfinite(temp_next))
    error('bearing_thermo:InvalidTemperatureUpdate', ...
        'The optional thermal update produced a nonfinite temperature.');
end
status = 'computed_open_loop';
end

function value = positive_two_node(value, name)
value = two_node_value(value,name);
if any(value<=0)
    error('bearing_thermo:InvalidPositiveTwoNodeValue', ...
        '%s must be positive.',name);
end
end

function value = nonnegative_two_node(value, name)
value = two_node_value(value,name);
if any(value<0)
    error('bearing_thermo:InvalidNonnegativeTwoNodeValue', ...
        '%s must be nonnegative.',name);
end
end

function value = bounded_positive_unit_scalar(value, name)
value = finite_scalar(value,name);
if value<=0 || value>1
    error('bearing_thermo:InvalidRelaxation', ...
        '%s must be within (0,1].',name);
end
end
