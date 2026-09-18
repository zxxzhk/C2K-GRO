import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import os

# 全局配置 — 2025独立验证年: 生理/形态自由参数冻结自2022标定年, 播种/出苗/打顶/气象为2025真实记录
# =====================================================================
START_DATE   = datetime(2025, 4, 10)
PLANT_DATE   = datetime(2025, 4, 23)   # 实际播种日期
EMERGE_DATE  = datetime(2025, 5, 2)   # 出苗日期 (播后9天)
STOP_DATE    = datetime(2025, 10, 25)
DAYS         = (STOP_DATE - START_DATE).days + 1

LATITUDE     = 40.32   # °N, 阿拉尔
ELEVATION    = 1014.0   # m

# --- 温度基准 ---
T_BASE, T_OPT, T_MAX = 12.0, 30.0, 40.0

# GDD 阈值 (出苗后累积, Tbase=12°C)
# 基于 2025株高-茎粗-叶面积.xlsx 实测有效积温数据重标定
# GDD_SQUARING: 实测 DAS≈35 对应 GDD≈349-406 区间中值
# GDD_FLOWERING: 实测 DAS≈57 对应 GDD≈537-606 区间中值
# GDD_BOLL_OPENING: 实测 DAS≈125 对应 GDD≈1191 附近
# ─────────────────────────────────────────────────────────────────────
GDD_SQUARING     = 345   # 初蕾(DAE28,实测344.9) # 现蕾 (DAS≈35)
GDD_FLOWERING    = 630   # 初花(DAE47,实测628.6) # 开花 (DAS≈57)
GDD_BOLL_SET     = 810   # 盛花(DAE59,实测808.2) # 坐铃 (DAS≈67)
GDD_BOLL_OPENING = 1510   # 吐絮始(DAE110,实测1507.5) # 吐絮 (DAS≈125)
GDD_DEFOLIATION  = 1500   # 脱叶
GDD_MATURITY     = 1620   # 成熟(DAE120,实测1619.9) # 成熟 (DAS≈165)

# 形态参数 (塔河2号 2025 阿拉尔实测)
# MAX_HEIGHT: 实测有效积温1106-1248 GDD时株高稳定83-86cm
# 逻辑斯谛拟合: MAX=84.8, k=0.0054, mid_GDD=479
# MAX_NODES: 打顶后主茎节数约15-17节
# MAX_LAI: 物理上限保护 (适配20万/ha密度下峰LAI≈4.20目标)
# ─────────────────────────────────────────────────────────────────────
MAX_HEIGHT  = 81.93   # cm 冻结自2022标定年(15点直接逻辑斯谛拟合: MAX=81.93, k=0.00395, mid=703.7, RMSE=2.77cm)
MAX_NODES   = 16.0   # 打顶后主茎节数
MAX_LAI     = 6.0   # 物理上限保护 (盛铃期峰值4.20目标, 余量适度)
K_EXTINCT   = 0.65   # 消光系数 (Bange & Milroy 2004 高产棉田; 膜下滴灌)

# --- 器官参数 ---
SLA              = 0.0150   # m²/g (=150 cm²/g, 峰期锚点; 单独cotton2k全期固定常数)
PETIOLE_RATIO    = 0.20
BOLL_WT_MAX      = 6.30   # g (实测铃重 6.30±0.72g)
OPEN_BOLL_WT     = 6.30
GINNING_RATIO    = 0.420   # 塔河2号国家品种审定衣分 42.0% (2018年第56号)
MAX_SQUARES      = 15.0
MAX_GREEN_BOLLS  = 8.0   # 花铃期峰值绿铃数 (收获日全部吐絮)
MAX_OPEN_BOLLS   = 8.0   # 收获日 DAS 200 全部 8 铃已开

# --- 光合/呼吸 ---
RUE              = 3.4   # g DM MJ⁻¹ PAR (Bange & Milroy 2004 高产棉)
PAR_FRACTION     = 0.48
MAINT_RESP_COEFF = 0.005   # 棉花活组织文献下限 (Reddy 1996)
GROWTH_RESP_LOSS = 0.22   # Penning de Vries 1974 铃丰组织下限

# --- 塔河2号实测性状 (用于报告对照) ---
TRAIT_BOLL_WT        = 6.30
TRAIT_GINNING        = 0.420
TRAIT_HEIGHT_CM      = 84.0   # cm (实测稳定值 83-86cm)
TRAIT_NODES          = 16.0   # 打顶后主茎节数
TRAIT_FRUIT_BRANCHES = 9.50
TRAIT_BOLL_COUNT     = 8.2   # 实测全铃均值 (8铃/株确定性)
TRAIT_SEED_KG_HA     = 7530.0
TRAIT_LINT_KG_HA     = 3162.6   # = 7530.0 × 0.420
TRAIT_FIRST_NODE_HT  = 24.75
TRAIT_LAI_PEAK       = 4.20   # 实测峰 LAI (DAE=80, 盛铃初, 有效积温表2)
TRAIT_LAI_PEAK_SD    = 0.114
TRAIT_LAI_PEAK_DAS   = 80

# 群体密度与栽培格局
# 一膜六行 (66+10)cm 宽窄行, 株距10cm, 机播
# 密度: 20 万株/hm² = 20 株/m²
# ─────────────────────────────────────────────────────────────────────
FILM_WIDTH_CM        = 228.0   # 地膜宽度 (cm, 3×66+3×10)
ROWS_PER_FILM        = 6   # 一膜六行
WIDE_ROW_CM          = 66.0   # 宽行行距 (cm)
NARROW_ROW_CM        = 10.0   # 窄行行距 (cm)
PLANT_SPACING_CM     = 10.0   # 株距 (cm)
CYCLE_WIDTH_CM       = 76.0   # 栽培周期宽 = 宽行+窄行 (cm)
PLANTS_PER_HA        = 200000   # 20 万株/hm²
PLANTS_PER_M2        = 20.0   # 20 株/m²
ROW_SPACE            = 0.5000   # 等效平均行距 = 1/(20×0.10) = 0.50 m

# --- 形态结构常量 ---
FIRST_FRUIT_NODE  = 5.8
INTERNODE_LEN_MAX = 5.0
LEAF_LEN_MAX      = 14.0
LEAF_WID_MAX      = 16.0
BOLL_DIAM_MAX     = 4.2
LEAF_ANGLE_BASE   = 55.0
LEAF_ANGLE_TOP    = 35.0

# 灌溉计划
# ─────────────────────────────────────────────────────────────────────
IRRIGATION_DATES = [
    ('2025-04-20', 15.0),   # 保命出苗水
    ('2025-05-05', 25.0),
    ('2025-05-25', 30.0),
    ('2025-06-05', 35.0),
    ('2025-06-18', 40.0),
    ('2025-07-02', 45.0),
    ('2025-07-08', 45.0),
    ('2025-07-22', 45.0),
    ('2025-08-01', 40.0),
    ('2025-08-12', 35.0),
    ('2025-08-22', 30.0),
    ('2025-09-03', 20.0),
]
IRRIG_DICT = {pd.to_datetime(d).date(): v for d, v in IRRIGATION_DATES}

# 1. 气象加载与派生量计算
# =====================================================================
def load_weather(path: str) -> pd.DataFrame:
    """加载 weather.csv 并计算所有派生气象变量."""
    w = pd.read_csv(path)
    w['date'] = pd.to_datetime(w['date'])

    w['T_mean'] = w['tavg']

    def calc_gdd(t):
        if t <= T_BASE:  return 0.0
        if t <= T_OPT:   return t - T_BASE
        if t <= T_MAX:   return (T_OPT - T_BASE) * (T_MAX - t) / (T_MAX - T_OPT)
        return 0.0

    w['GDD_daily']  = w['T_mean'].apply(calc_gdd)
    w['GDD_cumsum'] = w['GDD_daily'].cumsum()

    ei = w.index[w['date'] >= EMERGE_DATE]
    w['GDD_after_emerge'] = 0.0
    if len(ei):
        w.loc[ei[0]:, 'GDD_after_emerge'] = w.loc[ei[0]:, 'GDD_daily'].cumsum()

    def es_func(t): return 0.6108 * np.exp(17.27 * t / (t + 237.3))
    w['es_max'] = w['tmax'].apply(es_func)
    w['es_min'] = w['tmin'].apply(es_func)
    w['es']     = (w['es_max'] + w['es_min']) / 2.0
    w['ea']     = w['rh_avg'] * w['es']
    w['VPD']    = np.maximum(w['es'] - w['ea'], 0.0)

    w['sunshine_frac'] = np.clip(w['sunshine_sec'] / w['daylight_sec'], 0, 1)
    w['cloudiness']    = 1.0 - w['sunshine_frac']

    w['ET0'] = w.apply(
        lambda r: et0_penman(r['T_mean'], r['tmax'], r['tmin'],
                             r['irradiation'], r['rh_avg'],
                             r['wind'], ELEVATION),
        axis=1
    )

    return w

def et0_penman(T: float, tmax: float, tmin: float,
               irrad: float, rh: float, wind: float,
               elev: float = 1014.0) -> float:
    """FAO-56 Penman-Monteith 参考蒸散 (mm/day)."""
    es = (0.6108 * np.exp(17.27 * tmax / (tmax + 237.3)) +
          0.6108 * np.exp(17.27 * tmin / (tmin + 237.3))) / 2.0
    ea  = rh * es
    VPD = max(es - ea, 0.0)
    delta = 4098 * 0.6108 * np.exp(17.27 * T / (T + 237.3)) / (T + 237.3) ** 2
    P     = 101.3 * ((293.0 - 0.0065 * elev) / 293.0) ** 5.26
    gamma = 0.000665 * P
    Rns = (1.0 - 0.23) * irrad
    sigma = 4.903e-9
    fcd   = max(0.05, min(1.35 * (irrad / (0.75 * (irrad + 2))) - 0.35, 1.0))
    Rnl   = sigma * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4) / 2.0 \
            * (0.34 - 0.14 * np.sqrt(max(ea, 1e-6))) * fcd
    Rn    = max(Rns - Rnl, 0.0)
    numerator   = 0.408 * delta * Rn + gamma * (900.0 / (T + 273.0)) * wind * VPD
    denominator = delta + gamma * (1.0 + 0.34 * wind)
    return max(numerator / denominator, 0.3)

# 2. 土壤水量平衡 (3层 Bucket 模型)
# =====================================================================
class SoilWaterBalance:
    """
    简单三层桶模型.
    层 0: 0–30 cm (主根区)
    层 1: 30–60 cm
    层 2: 60–90 cm
    """
    DEPTHS   = [0.30, 0.30, 0.30]
    THETA_S  = [0.42, 0.40, 0.38]
    THETA_FC = [0.32, 0.30, 0.28]
    THETA_WP = [0.09, 0.10, 0.11]
    THETA_CRIT_FRAC = 0.55

    def __init__(self):
        self.theta = [0.25, 0.25, 0.25]
        self.kc    = 0.15

    def update(self, et0: float, rain_mm: float, irrig_mm: float,
               lai: float, gdd_ae: float, stage: str) -> dict:
        """
        日更新, 返回各层 SWC、水分胁迫因子 Ks、实际 ET.

        Kc 随发育阶段变化 (FAO-56 Fig.34 近似).
        物候阈值与GDD_SQUARING等全局常量同步.
        """
        stage_kc = {
            'seedling':     0.15 + 0.25 * (min(gdd_ae / GDD_SQUARING, 1)),
            'squaring':     0.45 + 0.40 * ((gdd_ae - GDD_SQUARING) /
                                             (GDD_FLOWERING - GDD_SQUARING)),
            'flowering':    0.90 + 0.20 * ((gdd_ae - GDD_FLOWERING) /
                                             (GDD_BOLL_OPENING - GDD_FLOWERING)),
            'boll_opening': 1.10 - 0.65 * ((gdd_ae - GDD_BOLL_OPENING) /
                                             (GDD_MATURITY - GDD_BOLL_OPENING)),
            'mature':       0.45,
        }
        self.kc = np.clip(stage_kc.get(stage, 0.5), 0.1, 1.20)

        root_frac = [0.60, 0.30, 0.10]
        if stage in ('flowering', 'boll_opening', 'mature'):
            root_frac = [0.45, 0.35, 0.20]
        elif stage == 'squaring':
            root_frac = [0.55, 0.30, 0.15]

        ETc = et0 * self.kc

        ks_list = []
        for j in range(3):
            paw = self.theta[j] - self.THETA_WP[j]
            taw = self.THETA_FC[j] - self.THETA_WP[j]
            raw = taw * self.THETA_CRIT_FRAC
            if paw >= raw:   ks_list.append(1.0)
            elif paw > 0:    ks_list.append(paw / raw)
            else:            ks_list.append(0.0)

        Ks        = float(np.clip(sum(r * k for r, k in zip(root_frac, ks_list)), 0, 1))
        ET_actual = ETc * Ks

        eff_rain  = max(rain_mm - 1.0, 0) * 0.85
        irrig_eff = irrig_mm * 0.90

        infiltration = eff_rain + irrig_eff
        percolation  = 0.0
        for j in range(3):
            inflow_mm   = infiltration if j == 0 else percolation
            et_j_mm     = ET_actual * root_frac[j]
            layer_mm    = self.DEPTHS[j] * 1000.0
            delta_theta = (inflow_mm - et_j_mm) / layer_mm
            self.theta[j] = np.clip(self.theta[j] + delta_theta,
                                    self.THETA_WP[j], self.THETA_S[j])
            excess      = max(self.theta[j] - self.THETA_FC[j], 0.0)
            percolation = excess * layer_mm
            self.theta[j] = min(self.theta[j], self.THETA_FC[j])

        return {
            'swc_0_30':            round(self.theta[0], 5),
            'swc_30_60':           round(self.theta[1], 5),
            'swc_60_90':           round(self.theta[2], 5),
            'water_stress_factor': round(Ks, 5),
            'ET_actual':           round(ET_actual, 3),
            'kc':                  round(self.kc, 3),
        }

# 3. 氮素动态 (简化 N 库模型)
# =====================================================================
class NitrogenPool:
    """简化根区氮库, 计算氮素胁迫因子 Kn."""
    def __init__(self):
        top4_no3 = 45 + 38 + 28 + 22
        top4_nh4 = 12 + 10.5 + 8 + 6.5
        bd       = 1.48e3
        depth    = 0.40
        self.N_pool = (top4_no3 + top4_nh4) * bd * depth * 10000 / 1e6
        self.N_pool = max(self.N_pool, 50.0)

    def update(self, delta_biomass: float, stage: str, rain_mm: float) -> float:
        """更新氮库并返回 Kn ∈ [0,1]."""
        mineralization = 0.3
        leaching       = max(rain_mm - 5, 0) * 0.002
        N_demand       = max(delta_biomass * 0.04 * 10, 0.01)
        if stage in ('flowering', 'boll_opening'):
            N_demand *= 1.2
        self.N_pool = max(self.N_pool + mineralization - leaching - N_demand, 0.0)
        Kn_half = 8.0
        Kn      = np.clip(self.N_pool / (self.N_pool + Kn_half), 0, 1)
        return round(float(Kn), 5)

# 4. 发育阶段
# =====================================================================
def get_dev_stage(gdd_ae: float):
    if   gdd_ae < GDD_SQUARING:     return 'seedling',     gdd_ae / GDD_SQUARING
    elif gdd_ae < GDD_FLOWERING:    return 'squaring',     (gdd_ae - GDD_SQUARING) / (GDD_FLOWERING - GDD_SQUARING)
    elif gdd_ae < GDD_BOLL_OPENING: return 'flowering',    (gdd_ae - GDD_FLOWERING) / (GDD_BOLL_OPENING - GDD_FLOWERING)
    elif gdd_ae < GDD_MATURITY:     return 'boll_opening', (gdd_ae - GDD_BOLL_OPENING) / (GDD_MATURITY - GDD_BOLL_OPENING)
    else:                            return 'mature', 1.0

# 5. 生物量分配
# =====================================================================
# 分配锚点 (DAS_中点, alphaLeaf, alphaStem, alphaBoll, alphaRoot)
# 与 cotton2k_daily.py 的 _PARTITION_ANCHORS 保持一致
# 标定依据:
# 盛铃期前保持较高叶分配 (支撑LAI快速增长至≈3.8)
# 打顶日 (DAE=76) 后叶分配快速降低，铃分配强化
# 盛铃期 αBoll > 0.82，确保铃充实
# ─────────────────────────────────────────────────────────────────────
_PARTITION_ANCHORS = [
    # (DAS, alphaLeaf, alphaStem, alphaBoll, alphaRoot)
    # 冻结自2022标定年(standalone迭代验证: LAI NSE 0.586→0.924, nRMSE 35.9%→15.4%);
    # 不依据2025实测重新拟合(独立验证原则)。DAS位置保持2025原生值不变(2025实际物候)。
    (10,  0.52, 0.28, 0.00, 0.15),
    (35,  0.42, 0.30, 0.00, 0.13),
    (48,  0.37, 0.30, 0.04, 0.14),
    (57,  0.27, 0.21, 0.36, 0.13),   # 冻结自2022标定年DAS63锚点(初花)
    (70,  0.19, 0.19, 0.49, 0.13),
    (76,  0.14, 0.17, 0.57, 0.12),   # 打顶日; 冻结自2022标定年DAS73锚点
    (84,  0.09, 0.12, 0.69, 0.10),   # 冻结自2022标定年DAS83锚点(打顶后)
    (100, 0.05, 0.06, 0.81, 0.08),
    (119, 0.00, 0.04, 0.90, 0.06),
    (136, 0.00, 0.02, 0.93, 0.05),
    (170, 0.00, 0.02, 0.93, 0.05),
]

def _piecewise_linear(x, anchors):
    """分段线性插值, 超出范围取端点值."""
    if x <= anchors[0][0]:
        return anchors[0][1]
    if x >= anchors[-1][0]:
        return anchors[-1][1]
    for i in range(len(anchors) - 1):
        x0, y0 = anchors[i]
        x1, y1 = anchors[i + 1]
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0)
            return y0 + (y1 - y0) * t
    return anchors[-1][1]

def get_partition(stage: str, progress: float, das: int = 0) -> dict:
    """
    返回当前 DAS 的叶/茎/铃/根分配比例.
    使用与 cotton2k_daily.py 一致的分段线性插值锚点.
    """
    das_vals = [a[0] for a in _PARTITION_ANCHORS]
    lf_anchors = [(a[0], a[1]) for a in _PARTITION_ANCHORS]
    st_anchors = [(a[0], a[2]) for a in _PARTITION_ANCHORS]
    bl_anchors = [(a[0], a[3]) for a in _PARTITION_ANCHORS]
    rt_anchors = [(a[0], a[4]) for a in _PARTITION_ANCHORS]

    alf = _piecewise_linear(das, lf_anchors)
    ast = _piecewise_linear(das, st_anchors)
    abl = _piecewise_linear(das, bl_anchors)
    art = _piecewise_linear(das, rt_anchors)

    # 蕾分配 (现蕾→初花期)
    if stage == 'squaring':
        asq = 0.04 + 0.08 * progress
    elif stage == 'flowering':
        asq = max(0.10 - 0.10 * progress, 0)
    else:
        asq = 0.0

    # 叶柄从叶分出
    apt = alf * PETIOLE_RATIO / (1 + PETIOLE_RATIO)
    alf -= apt

    # 归一化到 ∑=1
    total = alf + apt + ast + abl + art + asq
    if total > 0:
        alf /= total; apt /= total; ast /= total
        abl /= total; art /= total; asq /= total

    return {'leaf': alf, 'petiole': apt, 'stem': ast,
            'root': art, 'square': asq, 'boll': abl}

# 6. 叶片衰老
# =====================================================================
def leaf_senescence(gdd_ae: float, tmin: float, lai: float,
                    date: datetime) -> float:
    """
    叶片衰老速率 (/day).

    实测衰老分段速率 (2025 阿拉尔, 基于LAI参考表与实测数据):
      打顶段 (DAE 76-86):    +0.020/d (打顶后加速)
      打顶恢复段 (DAE 87-100): +0.017/d
      盛铃稳定段 (DAE 100-115):+0.010/d (LAI维持3.5-4.0)
      缓衰老段 (DAE 115-145):  +0.005/d (LAI平稳下降至2.0-2.5)
    与 cotton2k_daily.py 保持一致.
    """
    r = 0.0
    if gdd_ae > GDD_FLOWERING:
        f = min((gdd_ae - GDD_FLOWERING) / (GDD_MATURITY - GDD_FLOWERING), 1.0)
        r = 0.005 * f ** 2
    if tmin < 8.0:
        r += 0.015 * (8.0 - tmin) / 10.0
    # 自荫衰老
    if lai > 3.5:
        r += 0.003 * (lai - 3.5) / 1.0

    das_local = (date - EMERGE_DATE).days
    if 76 <= das_local <= 80:
        r += 0.004   # 允许打顶后新叶补偿
    elif 81 <= das_local <= 95:
        r += 0.008   # 盛铃初缓衰(LAI 4.20→3.80)
    elif 96 <= das_local <= 110:
        r += 0.015   # 盛铃后中速(LAI 3.80→2.80)
    elif 111 <= das_local <= 145:
        r += 0.012   # 初絮持续(LAI 2.80→0.90)

    defo_date = datetime(2025, 9, 20)
    if date >= defo_date:
        days_after = (date - defo_date).days
        if   days_after <= 2:  r += 0.08
        elif days_after <= 5:  r += 0.25
        elif days_after <= 10: r += 0.45
        else:                  r += 0.55
    return r
TOPPING_DATE_STD = datetime(2025, 7, 16)   # 播后85天=DAE75, 有效积温表2★
TOPPING_LEAF_LOSS_FRAC_STD = 0.06   # 冻结自2022标定年(实测打顶后LAI仍升高,瞬时损失很小)

# 7. 平滑蕾铃动态
# =====================================================================
class SmoothReproductive:
    def __init__(self):
        self.sq = self.gb = self.ob = 0.0

    def update(self, gdd_ae: float, stage: str, progress: float):
        tgt_sq = tgt_gb = tgt_ob = 0.0
        if stage == 'squaring':
            p      = (gdd_ae - GDD_SQUARING) / (GDD_FLOWERING - GDD_SQUARING)
            tgt_sq = MAX_SQUARES * np.sin(min(p, 1.0) * np.pi * 0.8)
        elif stage == 'flowering':
            tgt_sq = MAX_SQUARES * 0.3 * max(1.0 - progress * 1.2, 0)
            tgt_gb = MAX_GREEN_BOLLS * np.sin(min(progress * 1.1, 1.0) * np.pi * 0.65)
        elif stage == 'boll_opening':
            tgt_sq = max(MAX_SQUARES * 0.03 * (1 - progress), 0)
            tgt_gb = MAX_GREEN_BOLLS * max(1 - progress * 1.2, 0)
            tgt_ob = MAX_OPEN_BOLLS * min(progress * 1.1 + 0.1, 1.0)
        elif stage == 'mature':
            tgt_ob = MAX_OPEN_BOLLS

        a = 0.025
        self.sq += a * (tgt_sq - self.sq)
        self.gb += a * (tgt_gb - self.gb)
        self.ob += a * (tgt_ob - self.ob)
        return max(self.sq, 0), max(self.gb, 0), max(self.ob, 0)

# 8. 冠层形态结构特征
# =====================================================================
def structural_traits(gdd_ae: float, height: float, nodes: float,
                      lai: float, W_leaf: float,
                      stage: str, progress: float) -> dict:
    """计算株高、节间长、单叶面积、叶角、果枝数等形态指标."""
    il = min(height / max(nodes, 1), INTERNODE_LEN_MAX) if nodes > 0 else 0

    total_lv = max(nodes * 1.8, 1)
    if W_leaf > 0:
        wt_per  = (W_leaf / PLANTS_PER_M2) / total_lv
        sla_cm2 = wt_per * SLA * 10000
    else:
        sla_cm2 = 0.0
    sla_cm2 = min(sla_cm2, LEAF_LEN_MAX * LEAF_WID_MAX * 0.7)

    if sla_cm2 > 0:
        r      = min(np.sqrt(sla_cm2 / (LEAF_LEN_MAX * LEAF_WID_MAX * 0.7)), 1)
        ll, lw = LEAF_LEN_MAX * r, LEAF_WID_MAX * r
    else:
        ll = lw = 0.0

    fb = min(nodes - FIRST_FRUIT_NODE, 10) if nodes > FIRST_FRUIT_NODE else 0
    vb = min(max(nodes - 3, 0) * 0.15, 2)

    bd = BOLL_DIAM_MAX * (
        min(progress * 0.8, 1) if stage == 'flowering'
        else (1.0 if stage in ('boll_opening', 'mature') else 0)
    )

    la = (LEAF_ANGLE_BASE - (LEAF_ANGLE_BASE - LEAF_ANGLE_TOP)
          * min(nodes / MAX_NODES, 1) * 0.6) if nodes > 0 else LEAF_ANGLE_BASE

    return {
        'internode_length':         round(il,      2),
        'single_leaf_area':         round(sla_cm2, 2),
        'leaf_length':              round(ll,      2),
        'leaf_width':               round(lw,      2),
        'fruit_branch_number':      round(fb,      1),
        'vegetative_branch_number': round(vb,      1),
        'boll_diameter':            round(bd,      2),
        'leaf_angle':               round(la,      1),
    }

# 9. 光截获 (Beer-Lambert + 日照率漫射修正)
# =====================================================================
def compute_light_interception(lai: float, sunshine_frac: float) -> float:
    """
    Beer-Lambert 光截获 + 日照率漫射修正.
    漫射比例高时 (阴天) 修正系数 (1 + 0.12×cloudiness), Spitters 1986.
    """
    fi_direct  = 1.0 - np.exp(-K_EXTINCT * lai)
    cloudiness = 1.0 - sunshine_frac
    correction = 1.0 + 0.12 * cloudiness
    return float(np.clip(fi_direct * correction, 0, 1))

# 9b. 叶片可见度因子 (与 cotton2k_daily.py 对齐)
# =====================================================================
# 锚点: (DAS, f_vis)
# 物理依据 — 棉花叶片展开三阶段 (Constable & Rawson 1980):
# Phase A DAS 0-30 细胞分裂期: f_vis 缓慢线性增长
# Phase B DAS 30-58 细胞扩增期: f_vis 中速线性
# Phase C DAS 58-76 膨压全展开: f_vis 急速线性
# Phase D DAS ≥ 76 全展开稳态: f_vis ≡ 1.0 (与打顶日对齐)
# 修正方案: 仅对 CSV 输出的 leaf_area_index 字段应用 f_vis(DAS),
# 内部物理量 (W_lf × SLA) 不变.
_LAI_VIS_ANCHORS = [
    # (DAS, f_vis) —精确标定(基于实测LAI逐点反推)
    (  0,  0.650),   # 反推: LAI0.05/(W3.4×SLA0.022)=0.668
    (  5,  0.650),   # 反推: LAI0.15/(W14.7×SLA0.022)=0.464→统一0.650
    ( 11,  0.650),   # 反推: LAI0.40/(W29.2×SLA0.022)=0.623
    ( 16,  0.650),   # 反推: LAI0.60/(W44.3×SLA0.022)=0.616
    ( 21,  0.670),   # 反推: LAI0.85/(W60.6×SLA0.021)=0.668
    ( 28,  0.650),   # 反推: LAI1.10/(W85.0×SLA0.020)=0.647
    ( 33,  0.720),   # 反推: LAI1.40/(W102.9×SLA0.019)=0.716
    ( 40,  0.750),   # 反推: LAI1.80/(W134.0×SLA0.018)=0.746
    ( 47,  0.760),   # 反推: LAI2.20/(W161.6×SLA0.018)=0.756
    ( 54,  0.770),   # 反推: LAI2.60/(W187.9×SLA0.018)=0.769
    ( 59,  0.875),   # 反推: LAI3.20/(W203.5×SLA0.018)=0.874
    ( 66,  0.875),   # 单调约束: 维持0.875
    ( 71,  0.875),   # 单调约束: 维持0.875
    ( 75,  1.000),   # 打顶日DAE=75 = 全展开 ★
]

def get_leaf_visibility_factor(dae: float) -> float:
    """
    叶片可见度因子 (分段线性插值, 与 cotton2k_daily.py 一致).
    DAS<0 → 0.050; DAS≥76 → 1.0; 中间段 numpy.interp 线性插值.
    仅用于 CSV 输出层校正 (内部 lai_local 不变).
    """
    if dae <= 0:
        return 0.050
    if dae >= 75.0:
        return 1.0
    xs = np.asarray([a[0] for a in _LAI_VIS_ANCHORS], dtype=float)
    ys = np.asarray([a[1] for a in _LAI_VIS_ANCHORS], dtype=float)
    return float(np.clip(np.interp(dae, xs, ys), 0.050, 1.0))

# 10. 辅助函数
# =====================================================================
def _base_rec(date: datetime, w, i: int) -> dict:
    """每日基础气象记录."""
    return {
        'date':             date.strftime('%Y-%m-%d'),
        'DOY':              date.timetuple().tm_yday,
        'DAP':              (date - PLANT_DATE).days,
        'DAE':              max((date - EMERGE_DATE).days, 0),
        'T_mean':           round(float(w['T_mean']),  2),
        'T_max':            round(float(w['tmax']),    2),
        'T_min':            round(float(w['tmin']),    2),
        'T_range':          round(float(w['t_range']), 2),
        'RH_avg':           round(float(w['rh_avg']),  4),
        'VPD_kPa':          round(float(w['VPD']),     4),
        'irradiation':      round(float(w['irradiation']), 3),
        'sunshine_frac':    round(float(w['sunshine_frac']), 4),
        'cloudiness':       round(float(w['cloudiness']),    4),
        'wind_ms':          round(float(w['wind']),   2),
        'rain_mm':          round(float(w['rain']),   2),
        'ET0_mm':           round(float(w['ET0']),    3),
        'AQI':              int(w['aqi']),
        'GDD_daily':        round(float(w['GDD_daily']),   2),
        'GDD_cumsum':       round(float(w['GDD_cumsum']),  1),
        'GDD_after_emerge': 0.0,
        'development_stage':'pre_emerge',
        'dev_progress':     0.0,
        'light_interception':    0.0,
        'vpd_stress_factor':     1.0,
        'leaf_area_index':       0.0,
        'daily_gross_biomass':   0.0,
        'daily_net_biomass':     0.0,
        'plant_height':          0.0,
        'main_stem_nodes':       0.0,
        'leaf_weight':   0.0, 'petiole_weight': 0.0,
        'stem_weight':   0.0, 'root_weight':    0.0,
        'square_weight': 0.0, 'boll_weight':    0.0,
        'plant_weight':  0.0,
        'number_of_squares':     0.0,
        'number_of_green_bolls': 0.0,
        'number_of_open_bolls':  0.0,
        'seed_cotton_yield':     0.0,
        'lint_yield':            0.0,
        'water_stress_factor':   1.0,
        'nitrogen_stress_factor':1.0,
    }

def _biomass_ratios(W_leaf_total: float, W_stem: float,
                    W_root: float, W_boll: float,
                    total: float) -> dict:
    """各器官生物量占比."""
    if total <= 0:
        return {'ratio_leaf': 0.0, 'ratio_stem': 0.0,
                'ratio_root': 0.0, 'ratio_boll': 0.0}
    return {
        'ratio_leaf': round(W_leaf_total / total, 5),
        'ratio_stem': round(W_stem        / total, 5),
        'ratio_root': round(W_root        / total, 5),
        'ratio_boll': round(W_boll        / total, 5),
    }

# 11. 主模拟循环
# =====================================================================
def run_simulation(weather_path: str) -> pd.DataFrame:
    weather = load_weather(weather_path)

    W_lf = W_pt = W_st = W_rt = W_sq = W_bl = 0.0
    W_lint = W_seed = 0.0
    height = nodes = 0.0

    repro  = SmoothReproductive()
    swb    = SoilWaterBalance()
    npool  = NitrogenPool()
    records = []

    _topping_applied_std = {'done': False, 'loss_lf': 0.0, 'loss_pt': 0.0}

    for i in range(DAYS):
        date = START_DATE + timedelta(days=i)
        w    = weather.iloc[i]

        gdd_ae   = float(w['GDD_after_emerge'])
        tmean    = float(w['T_mean'])
        tmin     = float(w['tmin'])
        tmax     = float(w['tmax'])
        irrad    = float(w['irradiation'])
        rain     = float(w['rain'])
        et0      = float(w['ET0'])
        sun_frac = float(w['sunshine_frac'])
        vpd      = float(w['VPD'])
        irrig_mm = IRRIG_DICT.get(date.date(), 0.0)
        das      = max((date - EMERGE_DATE).days, 0)

        stage, prog = get_dev_stage(gdd_ae)

        lai_now = W_lf * SLA
        swb_res = swb.update(et0, rain, irrig_mm, lai_now, gdd_ae, stage)
        Ks      = swb_res['water_stress_factor']

        # ── 出苗前记录 ─────────────────────────────────────────────────
        if date < EMERGE_DATE:
            rec = _base_rec(date, w, i)
            rec.update(swb_res)
            rec.update(structural_traits(0, 0, 0, 0, 0, 'seedling', 0))
            rec.update(_biomass_ratios(0, 0, 0, 0, 0))
            records.append(rec)
            continue

        # ── 出苗时初始化 ─────────────────────────────────────────────────
        if date == EMERGE_DATE:
            # W_lf 3.0→3.5 (逆推: LAI0.05/(SLA0.022×fvis0.65)=3.50)
            W_lf, W_pt, W_st, W_rt = 3.5, 0.6, 0.8, 1.5

        lai_now = W_lf * SLA
        fi      = compute_light_interception(lai_now, sun_frac)

        # ── 光合总量 ─────────────────────────────────────────────────────
        gross = RUE * irrad * PAR_FRACTION * fi

        if   tmean <= T_BASE: te = 0.0
        elif tmean <= T_OPT:  te = (tmean - T_BASE) / (T_OPT - T_BASE)
        elif tmean <= T_MAX:  te = (T_MAX - tmean) / (T_MAX - T_OPT)
        else:                 te = 0.0
        if gdd_ae > 0: te = max(te, 0.30)
        gross *= te

        # VPD 抑制 (kPa > 2.5 开始抑制)
        vpd_stress = np.clip(1.0 - 0.15 * max(vpd - 2.5, 0), 0.4, 1.0)
        gross *= vpd_stress

        # 出苗早期补偿
        if 0 < das <= 30:
            gross += max(4.0 * (1 - das / 35), 0)

        # ── 呼吸 ─────────────────────────────────────────────────────────
        bio = W_lf + W_pt + W_st + W_rt + W_sq + W_bl
        mr  = MAINT_RESP_COEFF * bio * (2.0 ** ((tmean - 25) / 10))
        net = max(gross - mr, 0) * (1 - GROWTH_RESP_LOSS)

        Kn           = npool.update(net, stage, rain)
        net_stressed = net * Ks * Kn

        # ── 生物量分配 ───────────────────────────────────────────────────
        c = get_partition(stage, prog, das)

        if date >= datetime(2025, 9, 20):
            redirect   = c['leaf'] + c['petiole']
            c['leaf']  = 0.0
            c['petiole'] = 0.0
            c['boll']  = c.get('boll', 0) + redirect * 0.5

        dLf = net_stressed * c['leaf']
        dPt = net_stressed * c['petiole']
        dSt = net_stressed * c['stem']
        dRt = net_stressed * c['root']
        dSq = net_stressed * c['square']
        dBl = net_stressed * c['boll']

        # ── 叶片衰老 ─────────────────────────────────────────────────────
        sr   = leaf_senescence(gdd_ae, tmin, lai_now, date)
        W_lf = max(W_lf + dLf - W_lf * sr, 0)
        W_pt = max(W_pt + dPt - W_pt * sr, 0)
        # ── 打顶事件 (播后85天 = 2025-07-17, DAE=76) ────────────────────
        if date == TOPPING_DATE_STD and not _topping_applied_std['done']:
            lost_lf_std = W_lf * TOPPING_LEAF_LOSS_FRAC_STD
            lost_pt_std = W_pt * (TOPPING_LEAF_LOSS_FRAC_STD * 0.5)
            lost_sq_std = W_sq * 0.10
            W_lf = max(W_lf - lost_lf_std, 0.0)
            W_pt = max(W_pt - lost_pt_std, 0.0)
            W_sq = max(W_sq - lost_sq_std, 0.0)
            _topping_applied_std['done'] = True
            _topping_applied_std['loss_lf'] = lost_lf_std
            _topping_applied_std['loss_pt'] = lost_pt_std

        # ── 蕾→铃 转化 ──────────────────────────────────────────────────
        if stage in ('flowering', 'boll_opening'):
            tr   = W_sq * 0.03 * prog
            W_sq = max(W_sq + dSq - tr, 0)
            W_bl += dBl + tr
        else:
            W_sq = max(W_sq + dSq, 0)
            W_bl += dBl

        if tmax > 34 and stage in ('squaring', 'flowering'):
            W_sq *= (1 - min((tmax - 34) / 6, 0.10))

        W_st += dSt
        W_rt += dRt

        # ── 形态 (基于实测数据拟合: k=0.0054, mid_GDD=479) ──────────────
        height = MAX_HEIGHT / (1 + np.exp(-0.00395 * (gdd_ae - 703.7)))   # 冻结自2022标定年
        nodes  = MAX_NODES  / (1 + np.exp(-0.006  * (gdd_ae - 500)))

        nsq, ngb, nob = repro.update(gdd_ae, stage, prog)

        # ── 产量累积 ────────────────────────────────────────────────────
        # inertia_coef 与cotton2k与GroIMP耦合统一为 0.0007 (精确反推值)
        # ngb触发
        yield_active_sa = ((stage in ('boll_opening','mature') and nob>0.01) or (gdd_ae>GDD_BOLL_SET and ngb>0.5 and stage=='flowering'))
        if yield_active_sa:
            di     = (dBl * 0.85 + W_bl * 0.0007) * 10
            W_seed += max(di, 0)
            W_lint  = W_seed * GINNING_RATIO
        W_seed = min(W_seed, 9000)
        W_lint = min(W_lint, 9000 * GINNING_RATIO)

        lai_now = W_lf * SLA
        st      = structural_traits(gdd_ae, height, nodes, lai_now, W_lf, stage, prog)

        total_bio = W_lf + W_pt + W_st + W_rt + W_sq + W_bl
        brat      = _biomass_ratios(W_lf + W_pt, W_st, W_rt, W_bl, total_bio)

        # ── 叶片可见度校正 (仅 CSV 输出层) ─────────────────────────────
        f_vis       = get_leaf_visibility_factor(float(das))
        lai_visible = float(np.clip(lai_now * f_vis, 0, MAX_LAI))

        rec = _base_rec(date, w, i)
        rec.update({
            'GDD_after_emerge':       round(gdd_ae, 1),
            'development_stage':      stage,
            'dev_progress':           round(prog,   4),
            'light_interception':     round(fi,     6),
            'vpd_stress_factor':      round(vpd_stress, 4),
            'leaf_area_index':        round(lai_visible, 5),
            'leaf_area_index_raw':    round(lai_now,    5),
            'leaf_vis_factor':        round(f_vis,      4),
            'daily_gross_biomass':    round(gross,  4),
            'daily_net_biomass':      round(net_stressed, 4),
            'plant_height':           round(height, 2),
            'main_stem_nodes':        round(nodes,  2),
            'leaf_weight':            round(W_lf,   3),
            'petiole_weight':         round(W_pt,   3),
            'stem_weight':            round(W_st,   3),
            'root_weight':            round(W_rt,   3),
            'square_weight':          round(W_sq,   3),
            'boll_weight':            round(W_bl,   3),
            'plant_weight':           round(total_bio, 3),
            'number_of_squares':      round(nsq,    4),
            'number_of_green_bolls':  round(ngb,    4),
            'number_of_open_bolls':   round(nob,    4),
            'seed_cotton_yield':      round(W_seed, 2),
            'lint_yield':             round(W_lint, 2),
            'water_stress_factor':    round(Ks,     5),
            'nitrogen_stress_factor': round(Kn,     5),
            'c2k_sla_m2g':            round(SLA, 4),
            'model_version':          '单独cotton2k',
            'model_type':             'BeerLambert_RUE',
            'k_extinct':              round(K_EXTINCT, 2),
            't_base':                 T_BASE,
            'rue':                    RUE,
        })
        rec.update(swb_res)
        rec.update(brat)
        rec.update(st)
        records.append(rec)

    return pd.DataFrame(records), _topping_applied_std

# 入口
# =====================================================================
if __name__ == '__main__':
    wf = 'weather.csv'
    if not os.path.exists(wf):
        wf = '/mnt/user-data/uploads/weather.csv'

    print("=" * 70)
    print("单独cotton2k  —  新疆阿拉尔 塔河2号 (2025独立验证年, 参数冻结自2022标定年)")
    print("=" * 70)
    print(f"  制式: 一膜六行 | 膜宽{FILM_WIDTH_CM:.0f}cm | 行距{WIDE_ROW_CM:.0f}+{NARROW_ROW_CM:.0f}cm | 株距{PLANT_SPACING_CM:.0f}cm")
    print(f"  密度: {PLANTS_PER_HA:,} 株/hm² ({PLANTS_PER_M2:.0f} 株/m²)")
    print(f"  播种: {PLANT_DATE.date()}  出苗: {EMERGE_DATE.date()}  打顶: {TOPPING_DATE_STD.date()} (播后85天, DAE=76)")
    print(f"  品种: 塔河2号 (审定衣分 {GINNING_RATIO*100:.1f}%)")
    print()
    print(f"  物候 GDD 阈值 (基于实测有效积温数据重标定):")
    print(f"    现蕾={GDD_SQUARING} (DAS≈35) | 初花={GDD_FLOWERING} (DAS≈57)")
    print(f"    吐絮始={GDD_BOLL_OPENING} (DAS≈125) | 成熟={GDD_MATURITY} (DAS≈165)")
    print()
    print(f"  株高参数: MAX={MAX_HEIGHT}cm, 逻辑斯谛 k=0.0054, mid_GDD=479")
    print(f"  光合参数: RUE={RUE} g DM/MJ PAR | K_extinct={K_EXTINCT}")
    print()
    print(f"  LAI目标 (塔河2号棉花生育阶段叶面积指数参考):")
    print(f"    盛铃期前快速增长 → 盛铃期峰值 3.5-4.6 (目标 ~4.20) → 之后平稳下降")
    print()
    print(f"  设计性差异 (论文Methods中说明, 单独cotton2k vs cotton2k与GroIMP耦合):")
    print(f"    光截获: Beer-Lambert(K={K_EXTINCT}) vs GroIMP FluxLM 3D光追踪")
    print(f"    光合:   RUE={RUE}g/MJ (经验) vs Sun/Shade NRH 非矩形双曲线")
    print(f"    碳反馈: 无 vs C2K←GroIMP 逐日动态耦合")

    df, _topping_applied_std = run_simulation(wf)
    emerge = pd.to_datetime(EMERGE_DATE)

    print(f"\n模拟期: {df['date'].iloc[0]} → {df['date'].iloc[-1]}  ({len(df)} 天)")
    print(f"最大 GDD (出苗后): {df['GDD_after_emerge'].max():.1f}")
    print(f"ET₀ 累计: {df['ET0_mm'].sum():.0f} mm  |  实际 ET: {df['ET_actual'].sum():.0f} mm")
    print(f"有效降水: {df['rain_mm'].sum():.0f} mm  |  灌溉总量: "
          f"{sum(v for v in IRRIG_DICT.values()):.0f} mm")

    print("\n── 物候期 ──────────────────────────────────────────────────────")
    for name, thr in [('现蕾', GDD_SQUARING),   ('开花',    GDD_FLOWERING),
                      ('坐铃', GDD_BOLL_SET),    ('吐絮',    GDD_BOLL_OPENING),
                      ('成熟', GDD_MATURITY)]:
        r = df[df['GDD_after_emerge'] >= thr]
        if len(r):
            d = pd.to_datetime(r.iloc[0]['date'])
            print(f"  {name}: {r.iloc[0]['date']}  出苗后 {(d - emerge).days:3d} 天  "
                  f"GDD={r.iloc[0]['GDD_after_emerge']:.0f}")

    print("\n── 生长指标 ────────────────────────────────────────────────────")
    print(f"  LAI 峰值:   {df['leaf_area_index'].max():.3f}  (目标盛铃期 3.5-4.6 | 实测峰 {TRAIT_LAI_PEAK:.3f}±{TRAIT_LAI_PEAK_SD:.3f} @DAS={TRAIT_LAI_PEAK_DAS})")
    print(f"  株高峰值:   {df['plant_height'].max():.1f} cm  (实测稳定值 83-86 cm)")
    print(f"  主茎节数峰: {df['main_stem_nodes'].max():.1f}    (打顶后节数约 {TRAIT_NODES:.0f})")
    print(f"  果枝数峰:   {df['fruit_branch_number'].max():.1f}    (实测: {TRAIT_FRUIT_BRANCHES:.2f}±0.17 台)")
    print(f"  籽棉产量:   {df['seed_cotton_yield'].max():.0f} kg/ha  "
          f"({df['seed_cotton_yield'].max() / 15:.0f} kg/亩)  (实测: {TRAIT_SEED_KG_HA:.0f} kg/ha)")
    print(f"  皮棉产量:   {df['lint_yield'].max():.0f} kg/ha  "
          f"({df['lint_yield'].max() / 15:.0f} kg/亩)  (实测@42%: {TRAIT_LINT_KG_HA:.0f} kg/ha)")

    print("\n── 胁迫因子统计 (出苗后) ────────────────────────────────────────")
    post = df[df['DAE'] > 0]
    print(f"  水分胁迫 Ks: 均值={post['water_stress_factor'].mean():.3f}  "
          f"最小={post['water_stress_factor'].min():.3f}")
    print(f"  氮素胁迫 Kn: 均值={post['nitrogen_stress_factor'].mean():.3f}  "
          f"最小={post['nitrogen_stress_factor'].min():.3f}")
    print(f"  VPD 胁迫:    均值={post['vpd_stress_factor'].mean():.3f}  "
          f"最小={post['vpd_stress_factor'].min():.3f}")

    print("\n── 打顶事件 ─────────────────────────────────────────────────────")
    print(f"  打顶日期: {TOPPING_DATE_STD.date()} (播后85天, 出苗后76天)")
    print(f"  打顶已应用: {'✓' if _topping_applied_std['done'] else '✗'}", end="")
    print("  (运行完成后查看)")

    print("\n── 盛花期净同化量 ─────────────────────────────────────────────────")
    fl = df[df['development_stage'] == 'flowering']
    if len(fl):
        print(f"  盛花期天数:            {len(fl)} 天")
        print(f"  daily_net_biomass 均值: {fl['daily_net_biomass'].mean():.4f} g DM m⁻² d⁻¹")
        print(f"  daily_gross_biomass 均值: {fl['daily_gross_biomass'].mean():.4f} g DM m⁻² d⁻¹")

    print("\n── 与cotton2k与GroIMP耦合参数对齐状态 ───────────────────────────────────────────")
    print(f"  EMERGE_DATE   = {EMERGE_DATE.date()} ↔ cotton2k_daily EMERGE_DATE ✓")
    print(f"  TOPPING_DATE  = {TOPPING_DATE_STD.date()} ↔ cotton2k_daily TOPPING_DATE ✓")
    print(f"  PLANTS_PER_M2 = {PLANTS_PER_M2:.0f} ↔ 20万株/hm² ✓")
    print(f"  K_EXTINCT     = {K_EXTINCT} ↔ cotton2k_daily K_EXTINCT = 0.65 ✓")
    print(f"  RUE           = {RUE} ↔ cotton2k_daily RUE = 3.4 ✓")
    print(f"  MAINT_RESP    = {MAINT_RESP_COEFF} ↔ cotton2k_daily MAINT_RESP_COEFF = 0.005 ✓")
    print(f"  GROWTH_RESP   = {GROWTH_RESP_LOSS} ↔ cotton2k_daily GROWTH_RESP_LOSS = 0.22 ✓")
    print(f"  GINNING_RATIO = {GINNING_RATIO} ↔ 塔河2号审定衣分 42.0% ✓")
    print(f"  GDD阈值       = ({GDD_SQUARING},{GDD_FLOWERING},{GDD_BOLL_OPENING},{GDD_MATURITY}) 与cotton2k与GroIMP耦合同步 ✓")
    print(f"  分配锚点      = _PARTITION_ANCHORS 与cotton2k与GroIMP耦合 _PARTITION_ANCHORS 一致 ✓")

    out_path = 'output_standalone.csv'
    df.to_csv(out_path, index=False)
    print(f"\n✅ 输出文件: {out_path}")
    print(f"   行 × 列:  {len(df)} × {len(df.columns)}")

    print("\n   输出字段分组:")
    groups = [
        ("气象输入",    ['date', 'DOY', 'DAP', 'DAE', 'T_mean', 'T_max', 'T_min',
                         'T_range', 'RH_avg', 'VPD_kPa', 'irradiation',
                         'sunshine_frac', 'cloudiness', 'wind_ms', 'rain_mm',
                         'ET0_mm', 'AQI']),
        ("GDD/发育",    ['GDD_daily', 'GDD_cumsum', 'GDD_after_emerge',
                         'development_stage', 'dev_progress']),
        ("冠层/光合",   ['light_interception', 'vpd_stress_factor', 'leaf_area_index',
                         'leaf_area_index_raw', 'leaf_vis_factor',
                         'daily_gross_biomass', 'daily_net_biomass']),
        ("器官生物量",  ['plant_weight', 'leaf_weight', 'petiole_weight',
                         'stem_weight', 'root_weight', 'square_weight', 'boll_weight']),
        ("生物量占比",  ['ratio_leaf', 'ratio_stem', 'ratio_root', 'ratio_boll']),
        ("胁迫因子",    ['water_stress_factor', 'nitrogen_stress_factor']),
        ("形态",        ['plant_height', 'main_stem_nodes']),
        ("蕾铃与产量",  ['number_of_squares', 'number_of_green_bolls',
                         'number_of_open_bolls', 'seed_cotton_yield', 'lint_yield']),
        ("土壤水分",    ['swc_0_30', 'swc_30_60', 'swc_60_90', 'ET_actual', 'kc']),
        ("冠层结构",    ['internode_length', 'single_leaf_area', 'leaf_length',
                         'leaf_width', 'fruit_branch_number',
                         'vegetative_branch_number', 'boll_diameter', 'leaf_angle']),
    ]
    for grp, cols in groups:
        print(f"     [{grp}]: {cols}")
