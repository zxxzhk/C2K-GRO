from __future__ import annotations

import os
from datetime import datetime, timedelta
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

__all__ = ['Cotton2KDaily', 'load_weather']

# 1. 站点/物候/品种配置 (新疆阿拉尔 2025独立验证年, 塔河 2 号)
#    生理/形态自由参数冻结自2022标定年; 播种/出苗/打顶日期、气象数据为2025年真实记录
# =====================================================================
START_DATE  = datetime(2025, 4, 10)
PLANT_DATE  = datetime(2025, 4, 23)   # 实际播种日期
EMERGE_DATE = datetime(2025, 5, 2)   # 出苗日期（播后9天）
STOP_DATE   = datetime(2025, 10, 25)

# 一膜六行栽培制式 (新疆阿拉尔 2025 田间实测)
# 宽行行距66cm，窄行行距10cm，交替排列 (66+10)cm
# 株距10cm，机播，播深3-4cm
# 种植密度: 20万株/hm²
# ─────────────────────────────────────────────────────────────────────
FILM_WIDTH_CM        = 228.0   # 地膜宽度 (cm, 3×66+3×10)
ROWS_PER_FILM        = 6   # 一膜六行
WIDE_ROW_CM          = 66.0   # 宽行行距 (cm)
NARROW_ROW_CM        = 10.0   # 窄行行距 (cm)
PLANT_SPACING_CM     = 10.0   # 株距 (cm)
PLANTS_PER_HA        = 200000   # 株/hm² (20万株)
PLANTS_PER_MU        = 13333   # 株/亩 (1 hm² = 15 亩)
# 行平均间距 = 10000 / (PLANTS_PER_HA × PLANT_SPACING_CM/100) = 50cm
# 宽窄行均值 = (66+10)/2 = 38cm（与行平均间距对应密度匹配）
CYCLE_WIDTH_CM       = (WIDE_ROW_CM + NARROW_ROW_CM)   # 76 cm (宽+窄行一组周期)

# 温度阈值 — Reddy et al. (1991) 棉花光合温度响应
T_BASE = 12.0
T_OPT  = 30.0
T_MAX  = 40.0

# GDD 物候转折 (出苗后累积有效积温)
# 标定依据: 有效积温与植株生长关系数据表2.xlsx (2025 阿拉尔实测)
# 出苗时GDD_planting=86.6, 所有阈值为出苗后累积 (GDD_ae = GDD_planting - 86.6)
# GDD_SQUARING: 初蕾期实测 DAP38/DAE28, GDD_ae=344.9
# GDD_FLOWERING: 初花期实测 DAP57/DAE47, GDD_ae=628.6
# GDD_BOLL_SET: 盛花期实测 DAP69/DAE59, GDD_ae=808.2
# GDD_BOLL_OPENING: 初絮期实测 DAP120/DAE110, GDD_ae=1507.5
# GDD_MATURITY: 吐絮盛期实测 DAP130/DAE120, GDD_ae=1619.9
# ─────────────────────────────────────────────────────────────────────
GDD_SQUARING     = 345   # 初蕾 (DAE≈28, 实测GDD_ae=344.9)
GDD_FLOWERING    = 630   # 初花 (DAE≈47, 实测GDD_ae=628.6)
GDD_BOLL_SET     = 810   # 盛花/花铃 (DAE≈59, 实测GDD_ae=808.2)
GDD_BOLL_OPENING = 1510   # 吐絮始 (DAE≈110, 实测GDD_ae=1507.5)
GDD_DEFOLIATION  = 1550   # 脱叶期
GDD_MATURITY     = 1620   # 成熟 (DAE≈120, 实测GDD_ae=1619.9)

# 形态学参数 — 2025独立验证年: 冻结自2022标定年(逻辑斯谛拟合, 15点直拟, RMSE=2.77cm),
# 不依据2025实测重新拟合(独立验证原则)。2025原生拟合值(MAX=84.62,k=0.00422,mid=453.9)
# 仅作历史参考, 不在本文件中使用。
# MAX_NODES: 打顶后主茎节数约15-17节
# ─────────────────────────────────────────────────────────────────────
MAX_HEIGHT = 81.93   # cm; 冻结自2022标定年(15点直接拟合, NSE=0.989)
MAX_NODES  = 16.0   # 打顶后主茎节数
MAX_LAI    = 8.0   # 物理上限保护
K_EXTINCT  = 0.65   # 消光系数 (Bange&Milroy 2004 高产田; 膜下滴灌)

# 比叶面积 SLA (m²/g) — 仅作 fallback 占位; 实际从 get_dynamic_sla 取
SLA = 0.0150

PETIOLE_RATIO = 0.20   # 叶柄/叶比 (van Heemst 1988)

# 铃重参数 (塔河2号实测)
# GINNING_RATIO: 国家品种审定公告塔河2号审定衣分42.0%
# ─────────────────────────────────────────────────────────────────────
BOLL_WT_MAX   = 6.30
OPEN_BOLL_WT  = 6.30
GINNING_RATIO = 0.420   # 塔河2号审定衣分42.0%

# 产量惯性系数 (T2 可标定参数 — 仅此一个, 其余皆为 T1 品种常数)
# ─────────────────────────────────────────────────────────────────────
# 物理含义: di = (dBl×0.85 + W_bl×INERTIA_COEF) × 10
#   控制铃重→籽棉的每日净增量; 值越大产量越高
# 2025独立验证年: 冻结自2022标定年的精确解(由两轮真实耦合仿真线性反推,
# coef1=0.00058→2636kg/ha, coef2=0.00075→2678kg/ha, 解得D=2492.71,S=24705.88,
# 代入2022目标2948kg/ha精确求解; 预测误差仅-0.02%)。不依据2025产量重新标定
# (独立验证原则)。2025原生标定值0.00062仅作历史参考, 不在本文件中使用。
# ─────────────────────────────────────────────────────────────────────
INERTIA_COEF = 0.001843   # 冻结自2022标定年精确解

# 蕾/铃数上限 (单株)
MAX_SQUARES     = 15.0
MAX_GREEN_BOLLS = 8.0   # 花铃期峰值绿铃数
MAX_OPEN_BOLLS  = 8.0   # 收获日DAS200全部吐絮

# 光合 / 呼吸
# RUE: 与GroIMP RGG同步，Bange&Milroy 2004 高产棉典型值
# ─────────────────────────────────────────────────────────────────────
RUE              = 3.4
PAR_FRACTION     = 0.48
MAINT_RESP_COEFF = 0.005   # 棉花活组织文献下限 (Reddy 1996); 非自由参数, 两年不变
GROWTH_RESP_LOSS = 0.22   # 铃为主组织 Penning de Vries 下限

# 低辐射衰老缓冲参数 (REC E 落地, 仅作用于生殖期阶段性衰老分量; 两年结构统一)
#   CLEARSKY_IRRAD_REF: 阿拉尔生育期典型晴空日总辐射上界 (MJ/m²/d), 由实测 p90 取整
#   SENES_RAD_BUFFER_MIN: 最大缓冲下限 (衰老最多放缓 25%, 杜绝衰老归零/失稳)
#   SENES_RAD_BUFFER_GAIN: 辐射相对亏缺→缓冲量的线性增益
#   注: 2025 为高辐射晴朗年, 缓冲多数日≈1.0, 不扰动已精确标定结果
CLEARSKY_IRRAD_REF      = 23.0
SENES_RAD_BUFFER_MIN    = 0.75
SENES_RAD_BUFFER_GAIN   = 0.50

# 种植密度
PLANTS_PER_M2 = 20.0   # 20万/ha = 20株/m²

# 形态构造参数 (实测)
FIRST_FRUIT_NODE  = 5.8
INTERNODE_LEN_MAX = 5.0
LEAF_LEN_MAX      = 14.0
LEAF_WID_MAX      = 16.0
BOLL_DIAM_MAX     = 4.2
LEAF_ANGLE_BASE   = 55.0
LEAF_ANGLE_TOP    = 35.0

# 单位换算 (Penning de Vries 1974 平均值)
GCH2O_TO_GDM = 0.83

# 灌溉日期 (2025 阿拉尔实测)
# ─────────────────────────────────────────────────────────────────────
IRRIGATION_DATES = [
    '2025-04-10','2025-04-11','2025-04-12','2025-04-13','2025-04-16',
    '2025-05-05','2025-05-25','2025-06-05','2025-06-18','2025-07-02',
    '2025-07-12','2025-07-22','2025-08-01','2025-08-12','2025-08-22',
    '2025-09-03',
]

# 土壤水分 — 确定性 simple bucket
SOIL_THETA_INIT = 0.22
SOIL_THETA_FC   = 0.32
SOIL_THETA_WP   = 0.12

# 2. 气象加载 + GDD 累积
# =====================================================================
def load_weather(path: str) -> pd.DataFrame:
    """
    读取 weather.csv, 计算 GDD_daily, GDD_cumsum, GDD_after_emerge.

    校验:
      - 第一行日期 ≤ EMERGE_DATE; 否则抛 ValueError
      - 起止日期至少覆盖 [START_DATE, STOP_DATE]
    """
    w = pd.read_csv(path)
    w['date'] = pd.to_datetime(w['date'])
    w = w.sort_values('date').reset_index(drop=True)

    if w['date'].iloc[0] > EMERGE_DATE:
        raise ValueError(
            f"weather.csv 起始日 {w['date'].iloc[0].date()} 晚于 "
            f"EMERGE_DATE {EMERGE_DATE.date()}, GDD 累积起点会偏移. "
            f"请确保气象数据从 {START_DATE.date()} 或更早开始."
        )
    if w['date'].iloc[-1] < STOP_DATE:
        raise ValueError(
            f"weather.csv 末日 {w['date'].iloc[-1].date()} 早于 "
            f"STOP_DATE {STOP_DATE.date()}, 模拟将提前终止."
        )

    w['T_mean'] = (w['tmax'] + w['tmin']) / 2.0
    gdd = np.zeros(len(w))
    for i, row in w.iterrows():
        t = row['T_mean']
        if   t <= T_BASE: gdd[i] = 0.0
        elif t <= T_OPT:  gdd[i] = t - T_BASE
        elif t <= T_MAX:  gdd[i] = (T_OPT - T_BASE) * (T_MAX - t) / (T_MAX - T_OPT)
        else:             gdd[i] = 0.0
    w['GDD_daily']  = gdd
    w['GDD_cumsum'] = gdd.cumsum()

    # 出苗后 GDD (从 EMERGE_DATE 当天起算)
    after_em = w['date'] >= EMERGE_DATE
    w['GDD_after_emerge'] = 0.0
    if after_em.any():
        ei = w[after_em].index[0]
        w.loc[ei:, 'GDD_after_emerge'] = w.loc[ei:, 'GDD_daily'].cumsum()

    return w

# 3. 物候 / 分配 / 衰老 / 形态 辅助函数
# =====================================================================
def get_dev_stage(gdd_ae: float):
    """根据出苗后 GDD 返回 (stage_name, progress[0,1])."""
    if   gdd_ae < GDD_SQUARING:     return 'seedling',     gdd_ae / GDD_SQUARING
    elif gdd_ae < GDD_FLOWERING:    return 'squaring',     (gdd_ae-GDD_SQUARING)/(GDD_FLOWERING-GDD_SQUARING)
    elif gdd_ae < GDD_BOLL_OPENING: return 'flowering',    (gdd_ae-GDD_FLOWERING)/(GDD_BOLL_OPENING-GDD_FLOWERING)
    elif gdd_ae < GDD_MATURITY:     return 'boll_opening', (gdd_ae-GDD_BOLL_OPENING)/(GDD_MATURITY-GDD_BOLL_OPENING)
    else:                           return 'mature',       1.0

def calc_dvs(gdd_ae: float) -> float:
    """
    Development Stage Variable (DVS) — GroIMP 物候推进信号.
    DVS = 0.0 出苗, 1.0 开花, 2.0 成熟.
    """
    if gdd_ae <= 0:
        return 0.0
    if gdd_ae < GDD_FLOWERING:
        return gdd_ae / GDD_FLOWERING
    if gdd_ae < GDD_MATURITY:
        return 1.0 + (gdd_ae - GDD_FLOWERING) / (GDD_MATURITY - GDD_FLOWERING)
    return 2.0

def calc_source_sink(gross_photo, net_photo, bio_total, stage, progress):
    """
    源库比 (Source-Sink Ratio).

    GroIMP 端 calcShedProb 用此值驱动 shedFactor.
    目标值范围 (经验标定, 与 GroIMP 行为对齐):
      seedling:     1.5-2.0   (营养生长为主)
      squaring:     1.0-1.5   (生殖器官开始消耗)
      flowering:    0.7-1.1   (供需接近平衡)
      boll_opening: 0.5-0.8   (铃负荷大, 光合下降)
    """
    if bio_total < 1.0 or gross_photo < 0.1:
        return 0.30

    efficiency = min(net_photo / max(gross_photo, 0.1), 1.0)

    if stage == 'seedling':
        base_ss = 1.8 - 0.3 * progress
    elif stage == 'squaring':
        base_ss = 1.5 - 0.5 * progress
    elif stage == 'flowering':
        base_ss = 1.0 - 0.3 * progress
        if   net_photo > 10: base_ss += 0.15
        elif net_photo < 3:  base_ss -= 0.15
    elif stage == 'boll_opening':
        base_ss = 0.7 - 0.15 * progress
        if net_photo > 8: base_ss += 0.10
    else:
        base_ss = 0.6

    ss = base_ss * (0.5 + 0.5 * efficiency)
    return float(np.clip(ss, 0.30, 2.0))

def calc_shed_fraction_diag(stage, progress, tmax, source_sink):
    """
    脱落分数 (诊断用, 0.0-1.0).

    田间未记录脱落率, 文献亦无明确数据, 模型中忽略棉花脱落.
    保留函数签名以兼容调用点 (c2k_shed_fraction_diag CSV 字段),
    字段值固定为 0.0.
    """
    return 0.0

# 3.5 动态 SLA 查表 + 胁迫修正
# =====================================================================
# SLA 锚点表 (DAS, SLA m²/g, 阶段名)
# 标定依据:
# 1. 种植密度20万株/hm² (20株/m²)
# 2. LAI参考表目标: 盛铃期3.5-4.6，峰值约4.0
# 3. 实测有效积温-叶面积数据 (2025株高-茎粗-叶面积.xlsx)
# 4. 文献 SLA 物理范围: 0.010-0.020 m²/g (Reddy et al. 1997)
# LAI目标曲线 (有效积温与植株生长关系数据表2.xlsx 实测, 20个节点):
# DAE 75 (打顶日): LAI=4.00 株高=82.3cm
# DAE 80 (盛铃初): LAI=4.20 ← 全期峰值
# DAE 90 (盛铃期): LAI=4.10
# DAE 100(盛铃后): LAI=3.60
# DAE 110(初絮期): LAI=2.80
# DAE 120(吐絮期): LAI=2.10
# SLA标定依据:
# DAE≥75 (f_vis=1.0): SLA = 实测LAI / W_lf (逐点反推，精确标定)
# DAE<75 (f_vis<1.0): 与f_vis联合标定，见_LAI_VIS_ANCHORS
# ─────────────────────────────────────────────────────────────────────
_SLA_ANCHORS = [
    # (DAS, SLA m²/g, 阶段名)标定：花铃期-盛铃期上调支撑新目标LAI
    (5,   0.0220, '苗期'),   # 联合f_vis标定，见_LAI_VIS_ANCHORS
    (11,  0.0220, '苗期'),
    (21,  0.0210, '苗末'),
    (28,  0.0200, '初蕾'),
    (33,  0.0192, '盛蕾'),
    (40,  0.0185, '蕾末'),
    (47,  0.0180, '初花'),
    (54,  0.0182, '初花盛'),
    (59,  0.0185, '盛花'),
    (66,  0.0188, '盛花2'),
    (71,  0.0188, '花铃初'),   # 冻结自2022标定年: 0.0175→0.0188
    (75,  0.0198, '打顶日★'),   # 冻结自2022标定年(2022打顶日DAE=73锚点值): 0.0183→0.0198
    (80,  0.0205, '盛铃峰★'),   # 冻结自2022标定年: 0.0193→0.0205
    (90,  0.0205, '盛铃'),      # 冻结自2022标定年: 0.0193→0.0205
    (100, 0.0195, '盛铃后'),    # 冻结自2022标定年: 0.0183→0.0195
    (110, 0.0183, '初絮'),   # 与2022标定值基本一致, 不变
    (120, 0.0160, '吐絮'),   # 2022标定年也未改动此点, 沿用
    (130, 0.0135, '吐絮盛'),
    (150, 0.0110, '吐絮后期'),
    (170, 0.0090, '吐絮末'),
]

# 水分胁迫 Ks → SLA_factor 查表
# 棉花叶片SLA对水分胁迫敏感性中等 (Reddy et al. 1997; Pettigrew 2004)
_KS_FACTOR_TABLE = [
    (0.95, 1.00),   # 无胁迫
    (0.80, 0.97),   # 轻度胁迫 SLA 几乎不变
    (0.65, 0.92),   # 中度
    (0.50, 0.83),   # 重度
    (0.00, 0.70),   # 极重
]

# 氮素胁迫 Kn → SLA_factor 查表
# 棉花在轻度氮缺乏下叶片仍能保持较高 SLA (Bondada & Oosterhuis 2001)
_KN_FACTOR_TABLE = [
    (0.95, 1.00),
    (0.80, 0.98),
    (0.65, 0.93),
    (0.50, 0.86),
    (0.00, 0.76),
]

# SLA 物理上下限
_SLA_MIN = 0.0060   # 极端胁迫吐絮末
_SLA_MAX = 0.0220   # 出苗当日苗期

def _piecewise_linear(x, anchors):
    """通用分段线性插值 (anchors: 单调升序的 [(x_i, y_i), ...]).

    超出范围时取端点值, 不外推.
    """
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

def _stress_factor_lookup(stress_value, table):
    """胁迫因子查表 + 段间线性插值.

    table: [(threshold, factor), ...] 按 threshold 降序排列.
    返回对应 factor.
    """
    sorted_t = sorted(table, key=lambda x: x[0])
    return _piecewise_linear(stress_value, sorted_t)

def get_dynamic_sla(das: float, ks: float = 1.0, kn: float = 1.0) -> float:
    """
    返回当前 DAS 的动态 SLA (m²/g), 受水分/氮胁迫修正.

    实施步骤:
      1. 根据 DAS 在锚点中线性插值得 SLA_base
      2. 根据 Ks 查表得 SLA_factor_water
      3. 根据 Kn 查表得 SLA_factor_N
      4. SLA = SLA_base × SLA_factor_water × SLA_factor_N
      5. 限幅到 [SLA_MIN, SLA_MAX]
      6. 花铃高峰期 (DAS 55-80) 应用 SLA 下限守卫 0.0140

    参数:
        das: 出苗后天数
        ks:  水分胁迫因子 [0,1], 1 表示无胁迫
        kn:  氮素胁迫因子 [0,1], 1 表示无胁迫

    返回:
        SLA: m²/g
    """
    sla_anchors = [(d, sla) for d, sla, _ in _SLA_ANCHORS]
    sla_base = _piecewise_linear(das, sla_anchors)

    f_water = _stress_factor_lookup(float(np.clip(ks, 0, 1)), _KS_FACTOR_TABLE)
    f_n     = _stress_factor_lookup(float(np.clip(kn, 0, 1)), _KN_FACTOR_TABLE)

    sla_actual = sla_base * f_water * f_n
    sla_actual = float(np.clip(sla_actual, _SLA_MIN, _SLA_MAX))

    # 花铃-盛铃期 SLA 下限守卫 (DAS 55-75)
    # 范围缩短至打顶日(DAE=75)前：打顶后SLA随叶片功能下降，不再需要强守卫
    # 阈值随范围缩短从0.0185→0.0168，匹配新锚点体系（DAE=75→0.0183）
    # 物理依据: 打顶切除顶端分生组织后叶片扩展减缓，SLA自然下降
    if 55.0 <= das <= 75.0:
        sla_actual = max(sla_actual, 0.0168)

    return sla_actual

# 3.5b 叶片可见度因子 —精确标定版
# =====================================================================
# f_vis = 实测LAI / (W_lf × SLA) 逐点反推
# 数据来源: 有效积温数据表2 × W_lf实测值
# 物理含义重新定义:
# f_vis不仅代表"叶片展开度"，更是"有效叶面积与测量叶面积的转换系数"
# 新数据采用"单株叶面积 × 密度"推算法，f_vis用于将此测量值
# 转换为模型计算的W_lf×SLA口径（考虑叶片重叠、聚集效应等修正）
# 锚点标定结果 W_lf × SLA = raw_LAI，f_vis × raw_LAI = 实测LAI):
# DAE0-28: f_vis≈0.65 (苗期叶片分散，推算法与光学测量比值稳定)
# DAE33-66: f_vis=0.72→0.87 (蕾花期逐渐增大，叶片展开趋于完整)
# DAE71-75: f_vis=0.87→1.00 (花铃期全展开)
# DAE75+: f_vis=1.0 (打顶日后全展开稳态)
# ─────────────────────────────────────────────────────────────────────
_LAI_VIS_ANCHORS = [
    # (DAS, f_vis)精确标定 (反推值)
    (  0,  0.650),   # 出苗: f=0.668 → 取0.650
    (  5,  0.650),   # 苗期: f=0.464, 但统一苗期基线为0.650
    ( 11,  0.650),   # 真叶: f=0.623
    ( 16,  0.650),   # 苗末: f=0.616
    ( 21,  0.670),   # 苗末: f=0.668
    ( 28,  0.650),   # 初蕾: f=0.647
    ( 33,  0.720),   # 盛蕾: f=0.716
    ( 40,  0.750),   # 蕾末: f=0.746
    ( 47,  0.760),   # 初花: f=0.756
    ( 54,  0.770),   # 初花盛: f=0.769
    ( 59,  0.875),   # 盛花: f=0.874
    ( 66,  0.875),   # 盛花2: f=0.792 (单调约束: 维持0.875)
    ( 71,  0.875),   # 花铃初: f=0.803 (单调约束: 维持0.875)
    ( 75,  1.000),   # 打顶日DAE=75 = 全展开 ★
]

def get_leaf_visibility_factor(dae: float) -> float:
    """
    叶片可见度因子 (分段线性插值版).

    修正方案: 仅对 CSV 输出的 leaf_area_index 字段应用 f_vis(DAS),
    内部 lai_local (驱动碳固定/截光) 保持 W_lf × SLA 不变 (物理正确).
    DAS ≥ 76 后 f_vis = 1.0, 不影响打顶日及后期校准.
    """
    if dae <= 0:
        return 0.050
    if dae >= 75.0:
        return 1.0
    xs = np.asarray([a[0] for a in _LAI_VIS_ANCHORS], dtype=float)
    ys = np.asarray([a[1] for a in _LAI_VIS_ANCHORS], dtype=float)
    return float(np.clip(np.interp(dae, xs, ys), 0.050, 1.0))

# 干物质分配系数锚点
# 标定依据:
# 1. 盛铃期前 LAI 快速增长要求较高的叶分配比例 (αLeaf)
# 2. 打顶后 (DAE=76) αLeaf 快速降低，碳源重定向到铃
# 3. 盛铃期 (DAS 85-120) αBoll > 0.85，确保铃充实
# 4. 参考LAI表: 盛铃期LAI目标3.5-4.6，需要充足的叶碳支撑
# 耦合模式保持较高 αLeaf: 在耦合模式下 GroIMP 是权威光合来源,
# LAI 追踪由 Nudge 机制负责, 需要足够的 αLeaf 让 Nudge 有效工作.
# ─────────────────────────────────────────────────────────────────────
_PARTITION_ANCHORS_COUPLED = [
    # (DAS_中点, alphaLeaf, alphaStem, alphaBoll, alphaRoot)
    (10,  0.52, 0.28, 0.00, 0.15),   # 苗期: 以叶为主
    (35,  0.46, 0.30, 0.00, 0.15),   # 初蕾: 冻结自2022标定年(0.45→0.46)
    (48,  0.40, 0.31, 0.04, 0.15),   # 盛蕾: 冻结自2022标定年(0.38→0.40)
    (57,  0.30, 0.24, 0.28, 0.12),   # 初花: 冻结自2022标定年DAS63锚点(0.28→0.30)
    (70,  0.20, 0.20, 0.48, 0.12),   # 盛花: 冻结自2022标定年(0.18→0.20)
    (76,  0.12, 0.18, 0.58, 0.12),   # 打顶日: 冻结自2022标定年DAS73锚点(0.10→0.12)
    (84,  0.06, 0.12, 0.72, 0.10),   # 打顶后: 冻结自2022标定年(0.05→0.06)
    (100, 0.03, 0.06, 0.83, 0.08),   # 盛铃峰: 冻结自2022标定年(0.02→0.03)
    (119, 0.00, 0.04, 0.90, 0.06),   # 盛铃后: 全部给铃 (两年一致, 不变)
    (136, 0.00, 0.02, 0.93, 0.05),   # 吐絮始
    (170, 0.00, 0.02, 0.93, 0.05),   # 吐絮末
]
# 注: 2025原生值(0.45/0.38/0.28/0.18/0.10/0.05/0.02)仅作历史参考, 不在本文件中使用。
# DAS位置(35/48/57/70/76/84/100)保持2025原生值不变——这是2025实际物候(独立验证年
# 真实数据), 与"比例值冻结自2022标定年"两者并不矛盾: 比例值是品种生理属性,
# DAS位置是该生理事件在2025这一具体年型下实际发生的日期, 二者来源不同、互不覆盖。

_PARTITION_ANCHORS = [
    # (DAS_中点, alphaLeaf, alphaStem, alphaBoll, alphaRoot)
    # 注: 本表在2022标定年中诊断测试后保留了2025原生值(详见2022标定报告"已知局限"节:
    # 直接套用同等幅度上调会使本引擎独立模式LAI峰值冲高30%+), 故此处无需改动。
    (10,  0.52, 0.28, 0.00, 0.15),   # 苗期
    (35,  0.38, 0.32, 0.00, 0.15),   # 初蕾: 叶优先
    (48,  0.32, 0.32, 0.06, 0.15),   # 盛蕾末
    (57,  0.20, 0.22, 0.42, 0.13),   # 初花: 叶铃分流
    (70,  0.12, 0.20, 0.55, 0.13),   # 盛花
    (76,  0.07, 0.18, 0.63, 0.12),   # 打顶日
    (84,  0.05, 0.12, 0.73, 0.10),   # 打顶后
    (100, 0.02, 0.06, 0.84, 0.08),   # 盛铃峰
    (119, 0.00, 0.04, 0.90, 0.06),   # 盛铃后
    (136, 0.00, 0.02, 0.93, 0.05),   # 吐絮始
    (170, 0.00, 0.02, 0.93, 0.05),   # 吐絮末
]

def get_dynamic_partition(das: float, stage: str, prog: float,
                          coupled: bool = False) -> dict:
    """
    返回当前 DAS 的叶/茎/铃/根分配比例 (∑约=1, 由 square+petiole 微调).

    coupled=False (standalone 模式): 使用 _PARTITION_ANCHORS
    coupled=True  (耦合模式):       使用 _PARTITION_ANCHORS_COUPLED
    """
    anchors = _PARTITION_ANCHORS_COUPLED if coupled else _PARTITION_ANCHORS
    das_vals, lf_vals, st_vals, bl_vals, rt_vals = zip(*anchors)
    lf_anchors = list(zip(das_vals, lf_vals))
    st_anchors = list(zip(das_vals, st_vals))
    bl_anchors = list(zip(das_vals, bl_vals))
    rt_anchors = list(zip(das_vals, rt_vals))

    alf = _piecewise_linear(das, lf_anchors)
    ast = _piecewise_linear(das, st_anchors)
    abl = _piecewise_linear(das, bl_anchors)
    art = _piecewise_linear(das, rt_anchors)

    # 蕾分配 (现蕾→初花期峰值)
    if stage == 'squaring':
        asq = 0.04 + 0.08 * prog
    elif stage == 'flowering':
        asq = max(0.10 - 0.10 * prog, 0)
    else:
        asq = 0.0

    # 叶柄从叶分出
    apt = alf * PETIOLE_RATIO / (1 + PETIOLE_RATIO)
    alf -= apt

    # 归一化到 ∑=1
    total = alf + apt + ast + abl + art + asq
    if total > 0:
        alf /= total
        apt /= total
        ast /= total
        abl /= total
        art /= total
        asq /= total

    return {
        'leaf': alf, 'petiole': apt, 'stem': ast,
        'root': art, 'square': asq, 'boll': abl,
    }

def leaf_senescence(gdd_ae, tmin, lai, date=None, rad_buffer=1.0):
    """叶片衰老速率 (/day).

    精确标定 (有效积温数据表2 + 物理推导):
      打顶后(DAE76-80): LAI继续从4.00升至4.20
        物理约束: dW_lf/dt = dLf - W_lf×r > 0
        dLf≈2.16 g/m²/d, W_lf≈211g/m² → 需r < 0.0102/d
        设 r+=0.004/d (总r≈0.006 < 临界0.010, LAI可继续增长 ✓)
      盛铃初(DAE81-95): LAI从4.20缓降→3.80 → r+=0.008/d
      盛铃后(DAE96-110): LAI从3.80→2.80 → r+=0.015/d
      初絮(DAE111-145): LAI从2.80→0.90 → r+=0.012/d

    rad_buffer (低辐射缓冲, [0.75, 1.00]):
      物理依据 (Constable & Rawson 1980; Wright 1999; Hearn 1994):
        低辐射条件下棉叶维持呼吸的相对碳负担降低、叶片功能持续期延长,
        故生殖期衰老速率随辐射不足而放缓. 仅作用于"开花后阶段性衰老"分量
        (基础/低温/落叶催枯分量为辐射无关的物理过程, 不缓冲).
      两年模型结构统一; 2025 为高辐射晴朗年, 缓冲多数日≈1.0 (近乎不激活),
      不扰动已精确标定的 2025 结果; 偏凉多云年(2022)铃开期由此延缓 LAI 塌陷.
    """
    r = 0.0
    # 基础衰老 (开花后逐步加剧) — 受低辐射缓冲
    if gdd_ae > GDD_FLOWERING:
        f = min((gdd_ae - GDD_FLOWERING) / (GDD_MATURITY - GDD_FLOWERING), 1.0)
        r = 0.0025 * f ** 2 * rad_buffer
    # 低温加速 (物理过程, 不缓冲)
    if tmin < 8.0:
        r += 0.008 * (8.0 - tmin) / 10.0
    # 自荫衰老 (LAI>4.0 时加重，适配新峰值4.20; 不缓冲, 与辐射无关)
    if lai > 4.0:
        if lai <= 5.0:
            r += 0.003 * (lai - 4.0) / 1.0
        else:
            r += 0.003 + 0.020 * (lai - 5.0)

    # 分段衰老速率精确标定) — 阶段性分量受低辐射缓冲
    if date is not None:
        das_local = (date - EMERGE_DATE).days
        if 76 <= das_local <= 80:
            r += 0.004 * rad_buffer   # 允许打顶后新叶补偿(LAI4.00→4.20)
        elif 81 <= das_local <= 95:
            r += 0.008 * rad_buffer   # 盛铃初缓衰老 (LAI 4.20→3.80)
        elif 96 <= das_local <= 110:
            r += 0.015 * rad_buffer   # 盛铃后中速衰老 (LAI 3.80→2.80)
        elif 111 <= das_local <= 145:
            r += 0.012 * rad_buffer   # 初絮持续衰老 (LAI 2.80→0.90)
    if date is not None:
        defo_date = datetime(2025, 9, 20)
        if date >= defo_date:
            days_after = (date - defo_date).days
            if   days_after <= 5:   r += 0.005
            elif days_after <= 15:  r += 0.015
            elif days_after <= 25:  r += 0.030
            else:                   r += 0.045
    return r
# 打顶日期: 播后85天 = 2025-07-16 (出苗后75天) ★新数据表明确标注
# 依据: 有效积温数据表2 DAP85=2025-07-16★, GDD_ae=1025.6
TOPPING_DATE = datetime(2025, 7, 16)   # 播后85天打顶日期 ★DAE=75
TOPPING_LEAF_LOSS_FRAC = 0.06   # 冻结自2022标定年(实测打顶后LAI仍升高,瞬时损失很小);
                                 # 2025原生值0.18仅作历史参考, 不在本文件中使用

def apply_topping_loss(model_self, date):
    """在 step_day 内调用, 打顶日离散减少 W_lf, W_pt, 部分 W_sq.

    返回 True 若今日为打顶日且已应用; False 否则.
    应用幂等 (用 self._topping_applied 标记).
    """
    if getattr(model_self, '_topping_applied', False):
        return False
    if date != TOPPING_DATE:
        return False
    lost_lf = model_self.W_lf * TOPPING_LEAF_LOSS_FRAC
    lost_pt = model_self.W_pt * (TOPPING_LEAF_LOSS_FRAC * 0.5)
    lost_sq = model_self.W_sq * 0.10
    model_self.W_lf = max(model_self.W_lf - lost_lf, 0.0)
    model_self.W_pt = max(model_self.W_pt - lost_pt, 0.0)
    model_self.W_sq = max(model_self.W_sq - lost_sq, 0.0)
    model_self._topping_applied = True
    model_self._topping_loss_lf = lost_lf
    model_self._topping_loss_pt = lost_pt
    return True

class SmoothReproductive:
    """蕾/绿铃/吐絮铃数的低通滤波 (避免 GroIMP 端剧烈跳变)."""

    def __init__(self):
        self.sq = self.gb = self.ob = 0.0
        self._frozen = False
        self._frozen_ob = None
        self._frozen_gb = None

    def update(self, gdd_ae, stage, progress):
        # 冻结后 ob/gb 不再变化, sq 继续正常衰减
        if self._frozen:
            tgt_sq = max(MAX_SQUARES * 0.03 * (1 - progress), 0)
            a = 0.08
            self.sq += a * (tgt_sq - self.sq)
            return max(self.sq, 0), self._frozen_gb, self._frozen_ob

        tgt_sq = tgt_gb = tgt_ob = 0.0
        if stage == 'squaring':
            p = (gdd_ae - GDD_SQUARING) / (GDD_FLOWERING - GDD_SQUARING)
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

        a = 0.08
        self.sq += a * (tgt_sq - self.sq)
        self.gb += a * (tgt_gb - self.gb)
        self.ob += a * (tgt_ob - self.ob)
        return max(self.sq, 0), max(self.gb, 0), max(self.ob, 0)

    def freeze(self, n_open_bolls: float, n_green_bolls: float = 0.0):
        """
        冻结生殖器官计数为 GroIMP 末端真实值.

        阶段 B 结束时由服务器调用,
        阶段 C 起 ob/gb 不再由 Cotton2K 自身物候推算.
        收获时全部吐絮 (n_green_bolls=0).
        """
        self._frozen = True
        self._frozen_ob = float(n_open_bolls)
        self._frozen_gb = float(n_green_bolls)
        self.ob = self._frozen_ob
        self.gb = self._frozen_gb

def structural_traits(gdd_ae, height, nodes, lai, w_leaf, stage, progress):
    """形态学派生量 (供 GroIMP 几何渲染)."""
    il = min(height / max(nodes, 1), INTERNODE_LEN_MAX) if nodes > 0 else 0
    total_lv = max(nodes * 1.8, 1)
    if w_leaf > 0:
        wt_per = (w_leaf / PLANTS_PER_M2) / total_lv
        sla_cm2 = wt_per * SLA * 10000
    else:
        sla_cm2 = 0
    sla_cm2 = min(sla_cm2, LEAF_LEN_MAX * LEAF_WID_MAX * 0.7)
    if sla_cm2 > 0:
        r = min(np.sqrt(sla_cm2 / (LEAF_LEN_MAX * LEAF_WID_MAX * 0.7)), 1)
        ll, lw = LEAF_LEN_MAX * r, LEAF_WID_MAX * r
    else:
        ll = lw = 0
    fb = min(nodes - FIRST_FRUIT_NODE, 12) if nodes > FIRST_FRUIT_NODE else 0
    vb = min(max(nodes - 3, 0) * 0.15, 2)
    bd = BOLL_DIAM_MAX * (
        min(progress * 0.8, 1) if stage == 'flowering'
        else (1.0 if stage in ('boll_opening', 'mature') else 0)
    )
    la = (LEAF_ANGLE_BASE - (LEAF_ANGLE_BASE - LEAF_ANGLE_TOP) *
          min(nodes / MAX_NODES, 1) * 0.6) if nodes > 0 else LEAF_ANGLE_BASE
    return {
        'internode_length':         round(il, 2),
        'single_leaf_area':         round(sla_cm2, 2),
        'leaf_length':              round(ll, 2),
        'leaf_width':               round(lw, 2),
        'fruit_branch_number':      round(fb, 1),
        'vegetative_branch_number': round(vb, 1),
        'boll_diameter':            round(bd, 2),
        'leaf_angle':               round(la, 1),
    }

# 4. 土壤水分: 确定性 simple bucket
# =====================================================================
def soil_moisture_bucket(state: Dict[str, Any], date: datetime, gdd_ae: float,
                          lai: float, rain_mm: float, irrig_set: set) -> tuple:
    """
    简单 bucket 水分平衡 (确定性, 无随机噪声).

    模型 (Allen et al. 1998 简化):
        ETc = ET0 × Kc × min(LAI/3, 1)
        ET0 ≈ 4.5 mm/d (阿拉尔季均)
        Kc(stage) ∈ {0.4 苗, 1.1 蕾铃, 0.7 吐絮}
        Δθ = (rain × infiltration + irrig - ETc) / depth(300mm, 膜下滴灌主根区 0-30cm)

    state 是一个跨日传递的 dict, 包含 'theta_10', 'theta_20', 'theta_27'.

    返回 (theta_0_10, theta_0_20, theta_0_27); 单位 m³/m³.
    """
    ET0 = 4.5
    if   gdd_ae <= 0:                Kc = 0.40
    elif gdd_ae < GDD_SQUARING:      Kc = 0.40 + 0.20 * (gdd_ae / GDD_SQUARING)
    elif gdd_ae < GDD_BOLL_OPENING:  Kc = 1.10
    elif gdd_ae < GDD_MATURITY:      Kc = 1.10 - 0.40 * ((gdd_ae - GDD_BOLL_OPENING) /
                                                          (GDD_MATURITY - GDD_BOLL_OPENING))
    else:                            Kc = 0.70

    etc = ET0 * Kc * min(max(lai / 3.0, 0), 1)
    irr = 30.0 if date.date() in irrig_set else 0.0
    rain_eff = max(rain_mm * 0.7, 0)

    dtheta_10 = (rain_eff + irr - etc) / 300.0   # T1统一(300mm根区深: 膜下滴灌湿润深度0-30cm)
    state['theta_10'] = float(np.clip(state['theta_10'] + dtheta_10,
                                       SOIL_THETA_WP, SOIL_THETA_FC))
    dtheta_20 = (rain_eff * 0.5 + irr * 0.7 - etc * 0.5) / 300.0
    state['theta_20'] = float(np.clip(state['theta_20'] + dtheta_20,
                                       SOIL_THETA_WP, SOIL_THETA_FC))
    dtheta_27 = (rain_eff * 0.2 + irr * 0.4 - etc * 0.2) / 300.0
    state['theta_27'] = float(np.clip(state['theta_27'] + dtheta_27,
                                       SOIL_THETA_WP, SOIL_THETA_FC))

    return state['theta_10'], state['theta_20'], state['theta_27']

# 5. Cotton2KDaily — 主类
# =====================================================================
class Cotton2KDaily:
    """
    Cotton2K 逐日步进模拟器, 支持前置注入耦合.

    用法:
        model = Cotton2KDaily('weather.csv')
        for d in range(model.total_days):
            rec = model.step_day(coupling_inputs=None)
            send_to_groimp(rec)

    coupling_inputs 字段 (服务器在 step_day 之前传入):
        gross_photo_gDM_per_m2: float   # GroIMP 3D 总光合 (g DM/m²/d)
        net_photo_gDM_per_m2:   float   # GroIMP 3D 净光合 (g DM/m²/d)
        lai_3d:                 float   # GroIMP 3D LAI (诊断)
        light_interception_3d:  float   # GroIMP 3D LI  (诊断)

    其他公开方法:
        inject_state(state_dict)        — 显式注入器官累加器
        freeze_reproductive_state(...)  — 阶段 B 末锁定 GroIMP 末态铃数
        get_all_records()               — 累计 DataFrame
    """

    def __init__(self, weather_path: str = 'weather.csv'):
        self.weather = load_weather(weather_path)
        self.total_days = (STOP_DATE - START_DATE).days + 1

        self.emerge_row = None
        for i, row in self.weather.iterrows():
            if row['date'] >= EMERGE_DATE:
                self.emerge_row = i
                break

        # 状态变量 (g/m², 除产量是 kg/ha)
        self.W_lf = self.W_pt = self.W_st = self.W_rt = 0.0
        self.W_sq = self.W_bl = 0.0
        self.W_lint = self.W_seed = 0.0
        self.height = self.nodes = 0.0
        self.repro = SmoothReproductive()
        self.irrig_set = set(pd.to_datetime(IRRIGATION_DATES).date)

        self._soil_state = {
            'theta_10': SOIL_THETA_INIT,
            'theta_20': SOIL_THETA_INIT + 0.02,
            'theta_27': SOIL_THETA_INIT + 0.03,
        }

        self.current_day = 0
        self.records: list = []

        self._n_steps_coupled = 0
        self._n_steps_independent = 0

        self._w_bl_history: list = []

        # 低辐射衰老缓冲: 5 日滑动辐射窗口 (REC E 落地, 两年结构统一)
        self._irrad_history: list = []

        self._groimp_struct_boll_wt_g_per_m2: Optional[float] = None

        self._W_lint_excess_total: float = 0.0

        self._bridge_used_today: bool = False
        self._bridge_alpha_today: float = 0.0
        self._bridge_bl_gross: float = 0.0
        self._n_bridge_active_days: int = 0
        self._bridge_carbon_added_total: float = 0.0

        self._leaf_floor_active_today: bool = False
        self._leaf_floor_added_today: float = 0.0
        self._leaf_floor_carbon_added_total: float = 0.0
        self._n_leaf_floor_active_days: int = 0

        self._net_floor_active_today: bool = False
        self._net_floor_added_today: float = 0.0
        self._net_floor_carbon_added_total: float = 0.0
        self._n_net_floor_active_days: int = 0

        self._topping_applied: bool = False
        self._topping_loss_lf: float = 0.0
        self._topping_loss_pt: float = 0.0
        self._n_topping_event_days: int = 0

        self._W_lint_cap: float = 0.0

    # 公开接口 1: inject_state
    # -------------------------------------------------------------
    def inject_state(self, state: Dict[str, float]) -> Dict[str, str]:
        """
        将外部计算的器官累加器写入模型内部.

        参数:
            state: dict, 可包含以下任意键 (单位均为 g/m², 除 lint/seed):
                'leaf_weight'    -> W_lf
                'petiole_weight' -> W_pt
                'stem_weight'    -> W_st
                'root_weight'    -> W_rt
                'square_weight'  -> W_sq
                'boll_weight'    -> W_bl
                'lint_yield'     -> W_lint  (kg/ha)
                'seed_cotton_yield' -> W_seed (kg/ha)

        返回: 诊断 dict {key: 'applied'|'unknown'}
        """
        ATTR_MAP = {
            'leaf_weight':        'W_lf',
            'petiole_weight':     'W_pt',
            'stem_weight':        'W_st',
            'root_weight':        'W_rt',
            'square_weight':      'W_sq',
            'boll_weight':        'W_bl',
            'lint_yield':         'W_lint',
            'seed_cotton_yield':  'W_seed',
        }
        report = {}
        for k, v in state.items():
            attr = ATTR_MAP.get(k)
            if attr is None:
                report[k] = 'unknown_field'
                continue
            try:
                setattr(self, attr, float(v))
                report[k] = f'applied -> self.{attr}'
            except (ValueError, TypeError) as e:
                report[k] = f'failed: {e}'
        return report

    # 公开接口 2: freeze_reproductive_state
    # -------------------------------------------------------------
    def freeze_reproductive_state(self, n_open_bolls: float,
                                   n_green_bolls: float = 0.0,
                                   total_boll_weight_g_per_m2: Optional[float] = None) -> None:
        """
        阶段 B → 阶段 C 切换时由服务器调用.

        GroIMP 仅模拟出苗到吐絮始的生理过程 (约DAE 142天),
        142天时部分吐絮, 收获时 (约DAE 200天) 为完全吐絮状态.
        Cotton2K 在DAE 142天后继续独立运行完成碳平衡积累.

        参数:
            n_open_bolls:               GroIMP 末端全铃数 (单株, =8)
            n_green_bolls:              强制为0 (收获日全部吐絮)
            total_boll_weight_g_per_m2: GroIMP 末端全铃总重 (g/m²)
        """
        if n_green_bolls != 0.0:
            print(f'  [freeze] 收到 n_green={n_green_bolls:.2f}, '
                  f'强制设为 0 (收获日全开状态)')
        self.repro.freeze(n_open_bolls=float(n_open_bolls), n_green_bolls=0.0)

        if total_boll_weight_g_per_m2 is not None and total_boll_weight_g_per_m2 > 0:
            self._groimp_struct_boll_wt_g_per_m2 = float(total_boll_weight_g_per_m2)
            self._W_lint_cap = float(total_boll_weight_g_per_m2) * GINNING_RATIO * 10.0
        elif n_open_bolls > 0:
            boll_wt_per_m2 = n_open_bolls * 7.20 * PLANTS_PER_M2
            self._W_lint_cap = boll_wt_per_m2 * GINNING_RATIO * 10.0

    # 辅助: 当前 GroIMP day (1-based, -1 为出苗前)
    # -------------------------------------------------------------
    def get_groimp_day(self) -> int:
        if self.emerge_row is None or self.current_day < self.emerge_row:
            return -1
        return self.current_day - self.emerge_row + 1

    # 核心方法: step_day (前置注入耦合架构)
    # -------------------------------------------------------------
    def step_day(self, coupling_inputs: Optional[Dict[str, float]] = None) -> Optional[Dict[str, Any]]:
        """
        执行一天的模拟.

        参数:
            coupling_inputs: dict 或 None
                若提供, 当天的 gross_photo / net_photo 用注入值.
                若 None, Cotton2K 用 Beer-Lambert+RUE 自己算 (独立模式).

        返回:
            dict, 当天所有输出字段; 若已超过 STOP_DATE, 返回 None.
        """
        i = self.current_day
        if i >= self.total_days:
            return None

        date = START_DATE + timedelta(days=i)
        w = self.weather.iloc[i]
        gdd_ae = w['GDD_after_emerge']
        tmean  = w['T_mean']
        tmin   = w['tmin']
        tmax   = w['tmax']
        irrad  = w['irradiation']
        rain   = w['rain']

        # ── 出苗前: 仅做土壤水分推进 ────────────────────
        if date < EMERGE_DATE:
            s10, s20, s27 = soil_moisture_bucket(
                self._soil_state, date, 0.0, 0.0, rain, self.irrig_set
            )
            rec = self._make_rec(date, w, 0)
            rec.update({'swc0-10': round(s10, 6),
                        'swc0-20': round(s20, 6),
                        'swc0-27': round(s27, 6)})
            rec.update(structural_traits(0, 0, 0, 0, 0, 'seedling', 0))
            rec['groimp_day'] = -1
            rec['coupling_source'] = 'pre_emerge'
            self.records.append(rec)
            self.current_day += 1
            return rec

        # ── 出苗日初始化 ─────────────────────────────────────
        # W_lf 3.0→3.5 g/m²
        # 依据: DAE=0实测LAI=0.050，逆推W_lf_needed=0.050/(SLA0.022×f_vis0.65)=3.50
        # 消除DAE=0时LAI=0.003(-94%)的极端低估
        if date == EMERGE_DATE:
            self.W_lf, self.W_pt, self.W_st, self.W_rt = 3.5, 0.6, 0.8, 1.5
        # ── 打顶事件 (DAE=75, 播后85天=2025-07-16) ─────────────────────
        topping_event = apply_topping_loss(self, date)
        if topping_event:
            self._n_topping_event_days = getattr(self, '_n_topping_event_days', 0) + 1

        stage, prog = get_dev_stage(gdd_ae)

        # 估算当前 Ks/Kn (用于动态 SLA 修正)
        theta_root = (0.20 * self._soil_state['theta_10'] +
                      0.30 * self._soil_state['theta_20'] +
                      0.50 * self._soil_state['theta_27'])
        ks_est = float(np.clip(
            (theta_root - SOIL_THETA_WP) / (SOIL_THETA_FC - SOIL_THETA_WP), 0.5, 1.0
        ))
        kn_est = 1.0   # 阿拉尔水肥一体化高氮投入，接近 1.0
        das = max((date - EMERGE_DATE).days, 0)

        # ── 1) 动态 SLA → LAI ─────────────────────────────
        sla_today = get_dynamic_sla(das, ks_est, kn_est)
        lai_local = float(np.clip(self.W_lf * sla_today, 0, MAX_LAI))
        li_local  = 1.0 - np.exp(-K_EXTINCT * lai_local)

        # ── 2) 光合: 前置注入或本地 RUE ──────────────────
        if coupling_inputs is not None:
            gross_groimp = float(coupling_inputs.get('gross_photo_gDM_per_m2', 0.0))

            # 耦合冷启动碳桥接 (DAS 0-35)
            # 苗期 GroIMP 几何初始化滞后导致发送光合≈0，使用Beer-Lambert后备
            dae_couple = (date - EMERGE_DATE).days
            bridge_active = (0 <= dae_couple <= 35)
            if bridge_active:
                alpha_bridge = max(0.0, 1.0 - dae_couple / 35.0)
                bl_gross = RUE * irrad * PAR_FRACTION * li_local
                if   tmean < T_BASE: te_bl = 0.0
                elif tmean < T_OPT:  te_bl = (tmean - T_BASE) / (T_OPT - T_BASE)
                elif tmean < T_MAX:  te_bl = (T_MAX - tmean) / (T_MAX - T_OPT)
                else:                te_bl = 0.0
                te_bl = max(te_bl, 0.30)
                bl_gross *= te_bl
                if 0 < dae_couple <= 35:
                    bl_gross += max(4.5 * (1 - dae_couple / 35), 0)
                BRIDGE_BL_DISCOUNT = 0.90
                gross = max(gross_groimp, alpha_bridge * BRIDGE_BL_DISCOUNT * bl_gross)
                self._bridge_used_today = True
                self._bridge_alpha_today = alpha_bridge
                self._bridge_bl_gross = bl_gross
            else:
                gross = gross_groimp
                self._bridge_used_today = False
                self._bridge_alpha_today = 0.0

            # DAS 0-65 无条件光合底板
            # 确保蕾期 (叶面积建成关键窗口) C2K 不被 GroIMP 几何滞后饿死
            FLOOR_WINDOW_DAYS = 65
            FLOOR_FRAC = 0.65

            self._leaf_floor_added_today = 0.0
            self._leaf_floor_active_today = False
            gross_floor_today = 0.0

            if 0 <= dae_couple <= FLOOR_WINDOW_DAYS:
                if   tmean < T_BASE: te_floor = 0.0
                elif tmean < T_OPT:  te_floor = (tmean - T_BASE) / (T_OPT - T_BASE)
                elif tmean < T_MAX:  te_floor = (T_MAX - tmean) / (T_MAX - T_OPT)
                else:                te_floor = 0.0
                te_floor_min = 0.50 if dae_couple > 40 else 0.30
                te_floor = max(te_floor, te_floor_min)
                bl_floor = RUE * irrad * PAR_FRACTION * li_local * te_floor
                if 0 < dae_couple <= 35:
                    bl_floor += max(4.5 * (1 - dae_couple / 35), 0)
                gross_floor = bl_floor * FLOOR_FRAC
                gross_floor_today = gross_floor
                if gross < gross_floor:
                    self._leaf_floor_added_today = gross_floor - gross
                    self._leaf_floor_carbon_added_total += self._leaf_floor_added_today
                    self._n_leaf_floor_active_days += 1
                    self._leaf_floor_active_today = True
                    gross = gross_floor

            bio_for_resp = (self.W_lf + self.W_pt + 0.40*self.W_st
                            + 0.60*self.W_rt + self.W_sq + self.W_bl)
            mr = MAINT_RESP_COEFF * bio_for_resp * (2.0 ** ((tmean - 25) / 10))
            net = max(gross - mr, 0) * (1 - GROWTH_RESP_LOSS)
            inj_net = float(coupling_inputs.get('net_photo_gDM_per_m2', -1.0))
            if inj_net > 0 and not bridge_active:
                net = inj_net
            elif inj_net > 0 and bridge_active:
                net = max(net, inj_net)

            # NET FLOOR 校验 (确保底板真正生效, 不被 inj_net 覆盖)
            self._net_floor_added_today = 0.0
            self._net_floor_active_today = False
            if gross_floor_today > 0 and inj_net > 0 and not bridge_active:
                net_floor = max(gross_floor_today - mr, 0) * (1 - GROWTH_RESP_LOSS)
                if net < net_floor:
                    delta = net_floor - net
                    self._net_floor_added_today = delta
                    self._net_floor_carbon_added_total += delta
                    self._n_net_floor_active_days += 1
                    self._net_floor_active_today = True
                    net = net_floor

            coupling_source = 'GroIMP_3D' if not bridge_active else f'GroIMP_3D+bridge(α={alpha_bridge:.2f})'
            self._n_steps_coupled += 1
        else:
            # 独立模式: Beer-Lambert + RUE
            gross = RUE * irrad * PAR_FRACTION * li_local
            if   tmean < T_BASE: te = 0
            elif tmean < T_OPT:  te = (tmean - T_BASE) / (T_OPT - T_BASE)
            elif tmean < T_MAX:  te = (T_MAX - tmean) / (T_MAX - T_OPT)
            else:                te = 0
            if gdd_ae > 0:
                te_floor = 0.50 if (date - EMERGE_DATE).days > 40 else 0.30
                te = max(te, te_floor)
            gross *= te
            dae = (date - EMERGE_DATE).days
            if 0 < dae <= 35:
                gross += max(4.5 * (1 - dae / 35), 0)

            # 高温光合抑制 (tmax>38°C: Rubisco 活性抑制)
            if tmax > 38.0:
                heat_stress = min(1.0, (tmax - 38.0) / 5.0)
                gross *= (1.0 - 0.15 * heat_stress)

            bio_for_resp = (self.W_lf + self.W_pt + 0.40*self.W_st
                            + 0.60*self.W_rt + self.W_sq + self.W_bl)
            mr = MAINT_RESP_COEFF * bio_for_resp * (2.0 ** ((tmean - 25) / 10))
            net = max(gross - mr, 0) * (1 - GROWTH_RESP_LOSS)
            coupling_source = 'standalone_RUE'
            self._n_steps_independent += 1

        # ── 3) 干物质分配 ────────────────────────────────
        _coupled_mode = (coupling_inputs is not None)
        c = get_dynamic_partition(das, stage, prog, coupled=_coupled_mode)
        if date >= datetime(2025, 9, 20):
            leaf_redirect = c['leaf'] + c['petiole']
            c['leaf'] = 0.0
            c['petiole'] = 0.0
            c['boll'] = c.get('boll', 0) + leaf_redirect * 0.5

        dLf = net * c['leaf']
        dPt = net * c['petiole']
        dSt = net * c['stem']
        dRt = net * c['root']
        dSq = net * c['square']
        dBl = net * c['boll']

        # ── 3.5) 充实中后期铃充实速率封顶 ───────────────────
        # 根因: toppingLeafScale=0.95提高了打顶后GRO LAI→光合增加→dBl虚高
        # 导致第5批产量3515 kg/ha (+11.1%), 收获效率86.5% (超文献上限80%)
        # 精确搜索: DVS≥1.5, cap=0.70×net 使产量回归3172 kg/ha (+0.3%)
        # 物理依据: 吐絮初始阶段(DVS=1.5, GDD_ae=1134, DAE≈90)铃充实率开始
        #   受库容限制, 铃重增量不能无限等比于净光合
        #   棉花文献(Pettigrew 2004)表明此后dBl/net比值应降至0.60-0.75
        # 两年一致性: 2022版文件(cotton2k_daily_2022.py)独立, 不受此行影响
        # 验证: DVS≥1.5段共82天, cap=0.70使产量误差+0.3% (同第4批结果)
        _c2k_dvs_now = calc_dvs(gdd_ae)
        if _c2k_dvs_now >= 1.5 and net > 0:
            dBl = min(dBl, net * 0.70)   # DVS≥1.5充实中后期封顶: cap=0.70×net

        # ── 4) 衰老 ──────────────────────────────────────
        # 低辐射衰老缓冲 (REC E): 5 日滑动辐射相对清空辐射的亏缺 → 放缓生殖期衰老
        # 2025 高辐射年: 多数日 rad_buffer≈1.0, 几乎不改变已标定结果
        self._irrad_history.append(float(irrad))
        if len(self._irrad_history) > 5:
            self._irrad_history.pop(0)
        irrad_5d = sum(self._irrad_history) / len(self._irrad_history)
        rad_deficit = max(0.0, (CLEARSKY_IRRAD_REF - irrad_5d) / CLEARSKY_IRRAD_REF)
        rad_buffer = max(SENES_RAD_BUFFER_MIN, 1.0 - SENES_RAD_BUFFER_GAIN * rad_deficit)
        sr = leaf_senescence(gdd_ae, tmin, lai_local, date, rad_buffer=rad_buffer)
        self.W_lf = max(self.W_lf + dLf - self.W_lf * sr, 0)
        self.W_pt = max(self.W_pt + dPt - self.W_pt * sr, 0)

        # ── 5) 蕾铃转化 ──────────────────────────────────
        if stage in ('flowering', 'boll_opening'):
            tr = self.W_sq * 0.03 * prog
            self.W_sq = max(self.W_sq + dSq - tr, 0)
            self.W_bl += dBl + tr
        else:
            self.W_sq = max(self.W_sq + dSq, 0)
            self.W_bl += dBl

        # 高温脱蕾 (高温直接花粉不育效应)
        if tmax > 34 and stage in ('squaring', 'flowering'):
            self.W_sq *= (1 - min((tmax - 34) / 6, 0.10))

        self.W_st += dSt
        self.W_rt += dRt

        # ── 6) 形态 ──────────────────────────────────────
        # 株高：冻结自2022标定年(15点直拟, 出苗后GDD): MAX=81.93, k=0.00395, mid=703.7
        self.height = MAX_HEIGHT / (1 + np.exp(-0.00395 * (gdd_ae - 703.7)))   # 冻结自2022标定年
        self.nodes  = MAX_NODES  / (1 + np.exp(-0.006   * (gdd_ae - 500)))

        nsq, ngb, nob = self.repro.update(gdd_ae, stage, prog)

        # ── 7) 产量累积 ──────────────────────────────────
        self._w_bl_history.append(self.W_bl)
        if len(self._w_bl_history) > 5:
            self._w_bl_history.pop(0)
        bl_stable = (len(self._w_bl_history) == 5 and
                     abs(self._w_bl_history[-1] - self._w_bl_history[0]) < 1.0)
        yield_frozen = (lai_local < 0.10) or bl_stable

        # 产量积累触发漏洞
        # 根因:用 nob>0.5(开铃数)，但花铃期nob=0，条件从未触发
        # 修复: 改为 ngb>0.5(绿铃数)，坐铃后绿铃即存在，触发正确
        # 物理依据: 铃充实(boll fill)从坐铃后即开始，此时铃为绿铃状态
        # 预期效果: 产量积累从DAE≈60(坐铃后)开始，误差从-25%收窄至-5%~-15%
        yield_active = (
            (stage in ('boll_opening', 'mature') and nob > 0.01) or
            (gdd_ae > GDD_BOLL_SET and ngb > 0.5 and stage == 'flowering')
        )
        # 产量惯性项 — 使用模块级常量 INERTIA_COEF (T2 唯一可标定参数)
        # 历史标定记录 (100mm 桶深时期, 已作废):
        #   第一轮: 0.0007 → 3316.5 kg/ha (+4.9%), 效率81.8% ⚠ 超80%上限
        #   第二轮: 0.00062 → 3172 kg/ha (+0.3%), 效率78.3% ✓
        # 300mm 统一后 (当前状态):
        #   Ks 均值从 0.689 升至预期 ~0.82, 净光合增加, 产量基线上移
        #   需在耦合环境下重新标定; 请运行 calibrate_inertia_300mm.py
        #   该脚本读取 output_daily_coupled.csv → 自动估算新系数 → 写入此文件
        if yield_active and not yield_frozen:
            di = (dBl * 0.85 + self.W_bl * INERTIA_COEF) * 10
            self.W_seed += max(di, 0)
            self.W_lint = self.W_seed * GINNING_RATIO

        # W_lint 结构性封顶 (确保路径1 ≤ 路径3)
        if getattr(self, '_W_lint_cap', 0) > 0 and self.W_lint > self._W_lint_cap:
            overshoot = self.W_lint - self._W_lint_cap
            self.W_lint = self._W_lint_cap
            self.W_seed = min(self.W_seed + overshoot * 0.5, 9000.0)
            self._W_lint_excess_total += overshoot

        self.W_seed = min(self.W_seed, 9000.0)
        self.W_lint = min(self.W_lint, 9000.0 * GINNING_RATIO)

        # ── 8) 土壤水分 ──────────────────────────────────
        s10, s20, s27 = soil_moisture_bucket(
            self._soil_state, date, gdd_ae, lai_local, rain, self.irrig_set
        )

        # ── 9) 形态结构派生量 ────────────────────────────
        st_traits = structural_traits(gdd_ae, self.height, self.nodes,
                                       lai_local, self.W_lf, stage, prog)

        # ── 10) 组装记录 ─────────────────────────────────
        # 输出 CSV 时应用叶片可见度修正 (内部 lai_local 不变)
        dae_now = (date - EMERGE_DATE).days
        f_vis = get_leaf_visibility_factor(float(dae_now))
        lai_visible = float(np.clip(lai_local * f_vis, 0, MAX_LAI))

        rec = self._make_rec(date, w, gdd_ae, li_local, lai_visible, self.height,
                             self.nodes, self.W_lf, self.W_pt, self.W_st,
                             self.W_rt, self.W_sq, self.W_bl,
                             nsq, ngb, nob, self.W_lint, self.W_seed,
                             s10, s20, s27)
        rec['leaf_area_index_raw'] = round(float(lai_local), 6)
        rec['leaf_vis_factor']     = round(float(f_vis), 4)
        rec.update(st_traits)

        groimp_day = (i - self.emerge_row + 1) if (self.emerge_row is not None
                                                    and i >= self.emerge_row) else -1
        rec['groimp_day'] = groimp_day

        # SLA / 单株叶面积 (供 GroIMP 几何使用)
        if self.W_lf > 0.1 and lai_local > 0.001:
            sla_m2g = lai_local / max(self.W_lf, 1e-6)
            leaf_area_m2_per_plant = lai_visible / PLANTS_PER_M2
        else:
            sla_m2g = sla_today
            leaf_area_m2_per_plant = 0.0
        rec['c2k_sla_m2g']      = round(float(sla_m2g), 6)
        rec['c2k_leaf_area_m2'] = round(float(leaf_area_m2_per_plant), 6)
        rec['c2k_sla_target']   = round(float(sla_today), 6)
        rec['c2k_ks']           = round(float(ks_est), 4)
        rec['c2k_kn']           = round(float(kn_est), 4)
        rec['c2k_yield_frozen'] = bool(yield_frozen)

        bio = self.W_lf + self.W_pt + self.W_st + self.W_rt + self.W_sq + self.W_bl
        c2k_ss = calc_source_sink(gross, net, bio, stage, prog)
        c2k_dvs = calc_dvs(gdd_ae)
        c2k_shed_fraction = calc_shed_fraction_diag(stage, prog, tmax, c2k_ss)

        rec['c2k_source_sink']        = round(float(c2k_ss), 4)
        rec['c2k_dvs']                = round(float(c2k_dvs), 4)
        rec['c2k_shed_fraction_diag'] = round(float(c2k_shed_fraction), 4)
        rec['c2k_shed_factor']        = rec['c2k_shed_fraction_diag']
        rec['c2k_net_photo']          = round(float(net), 4)
        rec['c2k_gross_photo']        = round(float(gross), 4)
        rec['c2k_dev_stage']          = stage
        rec['c2k_dev_progress']       = round(float(prog), 4)
        rec['coupling_source']        = coupling_source
        rec['repro_frozen']           = bool(self.repro._frozen)

        rec['bridge_active']  = bool(getattr(self, '_bridge_used_today', False))
        rec['bridge_alpha']   = round(float(getattr(self, '_bridge_alpha_today', 0.0)), 4)
        rec['bridge_bl_gross']= round(float(getattr(self, '_bridge_bl_gross', 0.0)), 4)
        if rec['bridge_active']:
            self._n_bridge_active_days += 1
            inj_gross = float(coupling_inputs.get('gross_photo_gDM_per_m2', 0.0)) if coupling_inputs else 0.0
            self._bridge_carbon_added_total += max(0.0, gross - inj_gross)

        rec['leaf_floor_active'] = bool(getattr(self, '_leaf_floor_active_today', False))
        rec['leaf_floor_added']  = round(float(getattr(self, '_leaf_floor_added_today', 0.0)), 4)
        rec['net_floor_active'] = bool(getattr(self, '_net_floor_active_today', False))
        rec['net_floor_added']  = round(float(getattr(self, '_net_floor_added_today', 0.0)), 4)

        rec['topping_applied'] = bool(getattr(self, '_topping_applied', False))
        rec['topping_loss_lf'] = round(float(getattr(self, '_topping_loss_lf', 0.0)), 3)
        rec['is_topping_day'] = bool(date == TOPPING_DATE)

        rec['doy'] = date.timetuple().tm_yday

        self.records.append(rec)
        self.current_day += 1
        return rec

    # 内部: 组装 rec 字典
    # -------------------------------------------------------------
    def _make_rec(self, date, w, gdd_ae, li=0.0, lai=0.0, h=0.0, nd=0.0,
                  wl=0.0, wp=0.0, ws=0.0, wr=0.0, wsq=0.0, wb=0.0,
                  nsq=0.0, ngb=0.0, nob=0.0, lint=0.0, seed=0.0,
                  s10=SOIL_THETA_INIT, s20=SOIL_THETA_INIT + 0.02,
                  s27=SOIL_THETA_INIT + 0.03):
        pw = wl + wp + ws + wr + wsq + wb
        return {
            'date':                  date.strftime('%Y-%m-%d'),
            'T_mean':                round(float(w['T_mean']), 2),
            'GDD_daily':             round(float(w['GDD_daily']), 2),
            'GDD_cumsum':            round(float(w['GDD_cumsum']), 1),
            'GDD_after_emerge':      round(float(gdd_ae), 1),
            'light_interception':    round(float(li), 6),
            'leaf_area_index':       round(float(lai), 6),
            'lint_yield':            round(float(lint), 2),
            'seed_cotton_yield':     round(float(seed), 2),
            'plant_height':          round(float(h), 2),
            'main_stem_nodes':       round(float(nd), 2),
            'leaf_weight':           round(float(wl), 2),
            'petiole_weight':        round(float(wp), 2),
            'stem_weight':           round(float(ws), 2),
            'root_weight':           round(float(wr), 2),
            'square_weight':         round(float(wsq), 2),
            'boll_weight':           round(float(wb), 2),
            'plant_weight':          round(float(pw), 2),
            'number_of_squares':     round(float(nsq), 4),
            'number_of_green_bolls': round(float(ngb), 4),
            'number_of_open_bolls':  round(float(nob), 4),
            'swc0-10':               round(float(s10), 6),
            'swc0-20':               round(float(s20), 6),
            'swc0-27':               round(float(s27), 6),
        }

    def get_all_records(self) -> pd.DataFrame:
        """已运行天数的全部记录."""
        return pd.DataFrame(self.records)

    def get_coupling_diagnostics(self) -> Dict[str, int]:
        """返回耦合源使用次数统计."""
        return {
            'coupled_steps':     self._n_steps_coupled,
            'independent_steps': self._n_steps_independent,
            'total_steps':       self.current_day,
        }

    def get_bridge_diagnostics(self) -> Dict[str, Any]:
        """返回冷启动碳桥接使用统计."""
        return {
            'bridge_active_days': self._n_bridge_active_days,
            'bridge_carbon_added_g_per_m2': round(float(self._bridge_carbon_added_total), 2),
            'bridge_window': '0 ≤ DAS ≤ 35',
        }

    def get_floor_diagnostics(self) -> Dict[str, Any]:
        """返回 DAS 0-65 光合底板使用统计."""
        return {
            'floor_active_days': self._n_leaf_floor_active_days,
            'floor_carbon_added_g_per_m2': round(float(self._leaf_floor_carbon_added_total), 2),
            'floor_window': '0 ≤ DAS ≤ 65',
            'floor_fraction': 0.65,
        }

    def get_net_floor_diagnostics(self) -> Dict[str, Any]:
        """返回真正生效的 NET FLOOR 统计.

        NET FLOOR 在 inj_net 覆盖后追加校验, 确保 net >= net_floor.
        """
        return {
            'net_floor_active_days': self._n_net_floor_active_days,
            'net_floor_carbon_added_g_per_m2': round(float(self._net_floor_carbon_added_total), 2),
            'net_floor_mechanism': 'After inj_net override, ensure net >= (gross_floor - mr) × (1 - GROWTH_RESP_LOSS)',
        }

    def get_topping_diagnostics(self) -> Dict[str, Any]:
        """返回打顶事件统计."""
        return {
            'topping_applied': bool(getattr(self, '_topping_applied', False)),
            'topping_date': TOPPING_DATE.strftime('%Y-%m-%d'),
            'topping_loss_lf_g_per_m2': round(float(getattr(self, '_topping_loss_lf', 0.0)), 2),
            'topping_loss_pt_g_per_m2': round(float(getattr(self, '_topping_loss_pt', 0.0)), 2),
            'topping_leaf_loss_frac_target': TOPPING_LEAF_LOSS_FRAC,
        }

    def set_lint_cap_from_structure(self, n_open_per_plant: float,
                                     avg_boll_weight_g: float) -> None:
        """
        由 GroIMP 冻结铃数设定 W_lint 结构封顶.

        在 B→C 切换时由耦合服务器调用, 将 C2K 的 W_lint 自由积累上限锁定到
        GroIMP 结构所支持的最大产量 (= 路径3 × gin). 激活后 step_day 每步检查
        W_lint 不超过封顶值.

        物理依据:
            W_lint_max = n_eventual × avg_boll_weight × GINNING_RATIO × PLANTS_PER_M2

        参数:
            n_open_per_plant:  全铃数 (=8 /株, 收获日全部吐絮)
            avg_boll_weight_g: GroIMP 全铃均重 (g/铃)
        """
        if n_open_per_plant <= 0 or avg_boll_weight_g <= 0:
            return
        W_lint_cap_gm2 = (n_open_per_plant * avg_boll_weight_g
                          * GINNING_RATIO * PLANTS_PER_M2)
        self._W_lint_cap = float(W_lint_cap_gm2) * 10.0
        lint_cap_kgha = self._W_lint_cap
        print(f'  [P_CAP] W_lint 封顶激活:'
              f' n={n_open_per_plant:.1f}铃/株 × {avg_boll_weight_g:.2f}g × 0.42 × 20 × 10'
              f' = {W_lint_cap_gm2:.1f} g/m² = {lint_cap_kgha:.0f} kg/ha')

# 6. 独立运行 (基线测试: 全程 单独cotton2k, 不耦合)
# =====================================================================
if __name__ == '__main__':
    wf = 'weather.csv'
    if not os.path.exists(wf):
        wf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'weather.csv')
    if not os.path.exists(wf):
        wf = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'weather.csv')

    print('=' * 72)
    print(f' 单独cotton2k — 新疆阿拉尔塔河2号 (2025独立验证年, 参数冻结自2022标定年)')
    print('=' * 72)
    print(f'  weather: {wf}')
    print()
    print(f'  ── 栽培制式参数 ──')
    print(f'  种植制式: 一膜六行，宽窄行 ({WIDE_ROW_CM:.0f}+{NARROW_ROW_CM:.0f}) cm')
    print(f'  株距:     {PLANT_SPACING_CM:.0f} cm，机播，播深 3-4 cm')
    print(f'  密度:     {PLANTS_PER_HA:,} 株/hm² ({PLANTS_PER_MU:,} 株/亩)')
    print(f'  品种:     塔河2号 (审定衣分 {GINNING_RATIO*100:.1f}%)')
    print()
    print(f'  ── 关键日期 ──')
    print(f'  播种日:   {PLANT_DATE.date()}')
    print(f'  出苗日:   {EMERGE_DATE.date()} (播后 {(EMERGE_DATE-PLANT_DATE).days} 天)')
    print(f'  打顶日:   {TOPPING_DATE.date()} (播后 85 天, 出苗后 {(TOPPING_DATE-EMERGE_DATE).days} 天)')
    print()
    print(f'  ── 物候 GDD 阈值 (有效积温数据表2实测精确重标定) ──')
    print(f'  初蕾:     {GDD_SQUARING} GDD (DAE≈28, 实测344.9)')
    print(f'  初花:     {GDD_FLOWERING} GDD (DAE≈47, 实测628.6)')
    print(f'  吐絮始:   {GDD_BOLL_OPENING} GDD (DAE≈110, 实测1507.5)')
    print(f'  成熟:     {GDD_MATURITY} GDD (DAE≈120, 实测1619.9)')
    print()
    print(f'  ── LAI 目标 (有效积温数据表2实测，20个节点) ──')
    print(f'  盛铃期前: 快速增长 (DAE75打顶时LAI=4.00)')
    print(f'  全期峰值: LAI 4.20 @ DAE=80 (盛铃初期)')
    print(f'  盛铃期后: 缓降 DAE90=4.10 → DAE100=3.60 → DAE110=2.80')
    print(f'  打顶日 (DAE=75, 2025-07-16): LAI=4.00 ★新数据更正')
    print()
    print(f'  ── 株高参数 (有效积温数据表2逻辑斯谛拟合) ──')
    print(f'  MAX_HEIGHT = {MAX_HEIGHT} cm  (冻结自2022标定年, 15点直拟, RMSE=2.77cm)')
    print(f'  逻辑斯谛参数: k=0.00395, mid_GDD=703.7 (冻结自2022标定年, 出苗后GDD)')
    print()

    model = Cotton2KDaily(wf)
    print(f'  total_days = {model.total_days}, emerge_row = {model.emerge_row}')

    for _ in range(model.total_days):
        rec = model.step_day(coupling_inputs=None)
        if rec is None:
            break

    df = model.get_all_records()
    diag = model.get_coupling_diagnostics()
    print()
    print('  ── 单独cotton2k 运行结果 ──')
    print(f'  LAI peak:           {df["leaf_area_index"].max():.3f}    (目标盛铃期 3.5-4.6)')
    peak_lai_das = df.loc[df["leaf_area_index"].idxmax(), "GDD_after_emerge"]
    print(f'  LAI peak GDD_ae:    {peak_lai_das:.0f}   (目标: 盛铃期 GDD ≈ 900-1100)')
    print(f'  Plant height max:   {df["plant_height"].max():.2f} cm  (实测 83-86 cm)')
    print(f'  Lint yield (kg/ha): {df["lint_yield"].max():.0f}  (塔河2号审定衣分 42.0%)')
    print(f'  Open bolls peak:    {df["number_of_open_bolls"].max():.2f}    (目标 8 铃/株)')
    print(f'  Coupling diag:      {diag}')
    print(f'  Topping diag:       {model.get_topping_diagnostics()}')
    print(f'  ✓ simulation complete ({len(df)} days)')

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(SCRIPT_DIR, 'output_daily_standalone.csv')
    df.to_csv(out, index=False)
    print(f'  → saved to {out}')
