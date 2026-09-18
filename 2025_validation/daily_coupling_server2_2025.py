# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
daily_coupling_server.py — cotton2k与GroIMP耦合 逐日耦合服务器
==================================================================
站点:    新疆阿拉尔 (40.55°N, 海拔 ~1012 m)
品种:    塔河 2 号 (国家品种审定公告 2018 年第 56 号; 审定衣分 42.0%)
栽培制式: 一膜六行 (膜宽 228cm, 宽行 66cm + 窄行 10cm, 株距 10cm)
密度:    20 万株/hm² = 20 株/m²
播种:    2025-04-23  出苗: 2025-05-02  打顶: 2025-07-16 (播后85天, DAE=75) [2025独立验证年]
身份:    生理/形态自由参数冻结自2022标定年(不依据2025实测重新拟合); 播种/出苗/打顶/气象/
        田间实测对照值(MEASURED_*)均为2025年真实记录, 用于盲验证评估
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import socket
import sys
import tempfile
import time
from datetime import datetime
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from cotton2k_daily import Cotton2KDaily, GCH2O_TO_GDM, GINNING_RATIO, INERTIA_COEF, \
    get_dynamic_sla, get_leaf_visibility_factor, SOIL_THETA_FC, SOIL_THETA_WP


# 1. 常量与单位换算
# =======================================================================
DEFAULT_HOST           = '127.0.0.1'
DEFAULT_PORT           = 9527
GROWTH_PERIOD          = 142   # 与 RGG 的 GROWTH_PERIOD 一致
SOCKET_RECV_CHUNK      = 65536
SOCKET_LINE_TERMINATOR = b'\n'

# 地块面积换算 (中国法定 GB/T 28415-2012)
MU_TO_M2 = 666.667
HA_TO_M2 = 10000.0
HA_TO_MU = HA_TO_M2 / MU_TO_M2   # ≈ 15.0

# 种植格局参数 (一膜六行栽培制式, 新疆阿拉尔 2025)
# 宽行 66cm + 窄行 10cm 交替排列, 株距 10cm
# 膜宽 = 3×66 + 3×10 = 228cm
# 密度 = 20 万株/hm² = 20 株/m²
# ─────────────────────────────────────────────────────────────────────
FILM_WIDTH_CM        = 228.0   # 地膜宽度 (cm, 3×66+3×10)
ROWS_PER_FILM        = 6   # 一膜六行
WIDE_ROW_CM          = 66.0   # 宽行行距 (cm)
NARROW_ROW_CM        = 10.0   # 窄行行距 (cm)
PLANT_SPACING_CM     = 10.0   # 株距 (cm)
CYCLE_WIDTH_CM       = 76.0   # 栽培周期宽 = 宽行+窄行 (cm)
GROUND_AREA_PER_PLANT_CM2 = 500.0   # 单株地面占地 = 1/20 m² = 500 cm²

# 塔河 2 号品种常数
PLANTS_PER_MU      = 13333
LINT_RATIO_DEFAULT = 0.420   # 国家审定衣分 42.0% (塔河2号 2018第56号)
PLANTS_PER_HA      = 200000   # 20 万株/hm²
BOLL_WT_TARGET     = (5.58, 7.02)   # 实测铃重 6.30±0.72 g

# 铃数目标 (8 铃/株确定性, 田间均值 8.2 铃/株)
BOLL_COUNT_TARGET  = 8.2

# LAI 目标范围 (塔河2号棉花生育阶段叶面积指数参考)
# 盛铃期前快速增长, 盛铃期达峰 (3.5-4.6), 之后平稳下降
# ─────────────────────────────────────────────────────────────────────
LAI_PEAK_TARGET    = (3.5, 4.6)   # 盛铃期 LAI 目标区间

# 净同化量单位换算 (Penning de Vries 1974)
NET_PHOTO_UNIT  = 'gCH2O'
NET_PHOTO_TO_DM = GCH2O_TO_GDM   # 0.83

# FAO-56 Tetens 常数
TETENS_A = 17.27
TETENS_B = 237.3

# 物候阶段划分 (出苗起天数, 基于实测有效积温数据重标定)
# seedling: DAS 1-35 (GDD 0-380, 苗期)
# squaring: DAS 36-57 (GDD 380-570, 现蕾-盛蕾)
# flowering: DAS 58-100 (GDD 570-1050, 初花-花铃)
# boll_opening: DAS 101-142(GDD 1050-1200+, 盛铃-吐絮)
# ─────────────────────────────────────────────────────────────────────
PHASES: Dict[str, Tuple[int, int]] = {
    'seedling':     (1,   35),
    'squaring':     (36,  57),
    'flowering':    (58,  100),
    'boll_opening': (101, 142),
}

# Nudge 消融实验开关 (由 --nudge-off 参数控制)
# True  → NudgeON:  Python 端在 step_day() 后对次日 W_lf 施加软约束修正 (方案A v2)
# False → NudgeOFF: 跳过修正; 同时向 RGG 发 nudge_enabled=false (RGG 端硬开关)
# 实现说明 (方案A v2):
#   step_day() 后执行, 用当天 visible LAI (W_lf×SLA×f_vis) 对比当天 GroIMP LAI
#   门槛 DAS>=45 规避 GroIMP 冷启动期 (ratio<0.90), W_lf 下限 0.5 g/m²
# 参数与 RGG 端 c2kNudgeLAI() 保持一致:
LAI_NUDGE_THRESHOLD = 0.15   # 相对偏差阈值 (15%), 低于此值不修正
LAI_NUDGE_TAU       = 6.0    # 指数收敛时间常数 (天)
LAI_NUDGE_MAX_RATE  = 0.20   # 单日修正上限 (gap 的 20%)
NUDGE_ENABLED = True   # 默认 NudgeON; --nudge-off 时改为 False

MEASURED_LINT_KGHA  = 3162.6   # = 7530.0 × 0.420
MEASURED_SEED_KGHA  = 7530.0
MEASURED_PLANT_HT   = 84.0   # cm (实测稳定值 83-86 cm)
MEASURED_NODES      = 16.0   # 打顶后主茎节数

# 实测 LAI 关键节点
# 更新依据: 有效积温与植株生长关系数据表2.xlsx (2025 阿拉尔, 20个实测节点)
# 峰值LAI=4.20@DAE=80(盛铃初期), 打顶日(DAE=75)LAI=4.00
MEASURED_LAI_PEAK     = 4.20   # 全期峰值LAI (DAE=80, 盛铃初期)
MEASURED_LAI_PEAK_SD  = 0.10   # 估算标准差
MEASURED_LAI_PEAK_DAS = 80   # 出苗后峰值天(盛铃初期DAP90=DAE80)

# 完整 20 个实测 LAI 锚点 (有效积温数据表2.xlsx, DAE 0-120)
MEASURED_LAI_SERIES = [
    (0,   0.05),  (5,   0.15),  (8,   0.22),  (11,  0.40),
    (16,  0.60),  (21,  0.85),  (28,  1.10),  (33,  1.40),
    (40,  1.80),  (47,  2.20),  (54,  2.60),  (59,  3.20),
    (66,  3.50),  (71,  3.80),  (75,  4.00),  (80,  4.20),
    (90,  4.10),  (100, 3.60),  (110, 2.80),  (120, 2.10),
]

# 实测产量 + 单株
MEASURED_OPEN_BOLLS_PER_PLANT    = 4.8   # 区一+区二均值: (5.4+4.2)/2
MEASURED_BOLL_WEIGHT_G           = 6.76   # 区一+区二均值: (6.84+6.68)/2
MEASURED_SEED_COTTON_G_PER_PLANT = 32.4   # 区一+区二均值: (36.82+28.05)/2

# 2. 气象预处理 (FAO-56 标准化)
# =======================================================================
def _safe_float(s: Any, default: float = 0.0) -> float:
    try:
        if s is None:
            return default
        s = str(s).strip()
        if s == '':
            return default
        return float(s)
    except (ValueError, TypeError):
        return default

def normalize_date(date_str: str) -> str:
    date_str = date_str.strip()
    if date_str.startswith('\ufeff'):
        date_str = date_str[1:]
    for fmt in ('%Y-%m-%d', '%Y/%m/%d', '%Y%m%d', '%Y.%m.%d'):
        try:
            return datetime.strptime(date_str, fmt).strftime('%Y-%m-%d')
        except ValueError:
            continue
    try:
        parts = date_str.replace('-', '/').split('/')
        if len(parts) != 3:
            raise ValueError
        return f"{int(parts[0]):04d}-{int(parts[1]):02d}-{int(parts[2]):02d}"
    except Exception:
        raise ValueError(f"无法解析日期格式: {date_str!r}")

def tetens_es(t_c: float) -> float:
    """Tetens-Murray 饱和水汽压 (kPa). FAO-56 Eq.11."""
    return 0.6108 * math.exp(TETENS_A * t_c / (t_c + TETENS_B))

def fao56_vpd(tmax: float, tmin: float, rh_avg: float) -> float:
    """FAO-56 Eq.12 日均 VPD (kPa)."""
    es_mean = 0.5 * (tetens_es(tmax) + tetens_es(tmin))
    ea = rh_avg * es_mean
    return max(es_mean - ea, 0.0)

def detect_rh_scale(rh_values: List[float]) -> str:
    vals = [v for v in rh_values if v > 0]
    if not vals:
        return 'invalid'
    vmax = max(vals)
    if vmax <= 1.05:
        return 'fraction'
    if vmax <= 101.0:
        return 'percent'
    return 'invalid'

def preprocess_weather(src_path: str) -> Tuple[str, Dict[str, Dict[str, Any]]]:
    with open(src_path, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames
    if not fieldnames:
        raise ValueError(f"气象文件无表头: {src_path}")

    rh_col = 'rh_avg' if 'rh_avg' in fieldnames else None
    rh_scale = 'fraction'
    if rh_col:
        sample_rh = [_safe_float(r.get(rh_col, 0)) for r in rows[:60]]
        rh_scale = detect_rh_scale(sample_rh)
        if rh_scale == 'invalid':
            print(f"  [气象] ⚠ RH 列含异常值, 按 0-1 处理")
            rh_scale = 'fraction'

    tmp = tempfile.NamedTemporaryFile(
        mode='w', suffix='.csv', delete=False, newline='', encoding='utf-8'
    )
    writer = csv.DictWriter(tmp, fieldnames=fieldnames)
    writer.writeheader()

    weather_index: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        try:
            norm_date = normalize_date(row['date'])
        except Exception as e:
            print(f"  [气象] 跳过异常行 {row.get('date','?')}: {e}")
            continue
        row['date'] = norm_date
        writer.writerow(row)

        try:
            tmax = _safe_float(row.get('tmax'))
            tmin = _safe_float(row.get('tmin'))
            tavg = _safe_float(row.get('tavg'), default=(tmax + tmin) / 2.0)
            rh_raw = _safe_float(row.get('rh_avg'))
            rh_avg = rh_raw / 100.0 if rh_scale == 'percent' else rh_raw
            rh_avg = max(0.0, min(1.0, rh_avg))
            tdew = _safe_float(row.get('tdew'), default=tmin - 4.0)
            daylight_sec = _safe_float(row.get('daylight_sec'))
            sunshine_sec = _safe_float(row.get('sunshine_sec'))
            irradiation = _safe_float(row.get('irradiation'))
            aqi = int(_safe_float(row.get('aqi')))

            vpd = round(fao56_vpd(tmax, tmin, rh_avg), 4)
            sunshine_frac = (round(sunshine_sec / daylight_sec, 4)
                             if daylight_sec > 0 else 0.0)

            weather_index[norm_date] = {
                'tmax': tmax, 'tmin': tmin, 'tavg': tavg,
                'rh_avg': rh_avg, 'tdew': tdew, 'vpd': vpd,
                'sunshine_frac': sunshine_frac,
                'daylight_sec': daylight_sec, 'sunshine_sec': sunshine_sec,
                'irradiation': irradiation, 'aqi': aqi,
            }
        except Exception as e:
            print(f"  [气象] 衍生字段失败 ({norm_date}): {e}")
            weather_index[norm_date] = {}

    tmp.close()
    print(f"  [气象] RH 单位识别: {rh_scale}")
    return tmp.name, weather_index

# 3. Socket I/O
# =======================================================================
def make_line_reader(conn: socket.socket) -> Callable[[], Dict[str, Any]]:
    buf = bytearray()

    def read_one() -> Dict[str, Any]:
        nonlocal buf
        while SOCKET_LINE_TERMINATOR not in buf:
            chunk = conn.recv(SOCKET_RECV_CHUNK)
            if not chunk:
                raise ConnectionError("GroIMP 断开连接 (对端关闭)")
            buf.extend(chunk)
        idx = buf.index(SOCKET_LINE_TERMINATOR)
        line = bytes(buf[:idx])
        del buf[:idx + 1]
        if not line:
            return {}
        try:
            return json.loads(line.decode('utf-8'))
        except json.JSONDecodeError as e:
            raise ConnectionError(f"JSON 解析失败: {e} ({line[:120]!r}...)")

    return read_one

def send_json(conn: socket.socket, data: Dict[str, Any]) -> None:
    msg = json.dumps(data, ensure_ascii=False, separators=(',', ':')) + '\n'
    conn.sendall(msg.encode('utf-8'))

# 4. GroIMP 反馈 → coupling_inputs 转换
# =======================================================================
_coupling_stats = {
    'attempts': 0,
    'used_RGG_3D_direct': 0,
    'used_python_rebuild': 0,
    'fallback_invalid':   0,
}

def build_coupling_inputs(
    groimp_feedback: Optional[Dict[str, Any]],
    last_rec: Optional[Dict[str, Any]],
    wx: Dict[str, Any],
) -> Tuple[Optional[Dict[str, float]], Dict[str, Any]]:
    """
    将 GroIMP 上一日反馈转换为 cotton2k_daily.step_day 的 coupling_inputs.

    优先级:
      1. 若 groimp_feedback 中含 RGG 的 3D 权威字段
         (groimp_net_photo_gCH2O / groimp_gross_photo_gCH2O), 直接换算
         g DM/m²/d 后填入 coupling_inputs.
      2. 若仅含 LI/LAI 没有 photo 字段, 用 RUE 在 Python 端基于 GroIMP-LI
         重算 gross/net (fallback, 兼容老版 RGG).
      3. 若 GroIMP 反馈无效, 返回 None (Cotton2K 走 standalone 路径).

    返回 (coupling_inputs_or_None, report_dict).
    """
    _coupling_stats['attempts'] += 1

    if groimp_feedback is None:
        return None, {'source': 'no_feedback', 'reason': 'first_day'}

    gli  = _safe_float(groimp_feedback.get('groimp_light_interception'), -1.0)
    glai = _safe_float(groimp_feedback.get('groimp_lai_3d'), -1.0)
    gpar = _safe_float(groimp_feedback.get('groimp_intercepted_par'), -1.0)

    if gli < 0.0 or gli > 1.0 or glai < 0.01:
        _coupling_stats['fallback_invalid'] += 1
        return None, {'source': 'invalid', 'reason': f'gli={gli:.3f},glai={glai:.3f}'}

    # 路径 1: RGG+ 直接发来的 3D 权威光合 (g CH2O/m²/d)
    gross_gch2o = _safe_float(groimp_feedback.get('groimp_gross_photo_gCH2O'), -1.0)
    net_gch2o   = _safe_float(groimp_feedback.get('groimp_net_photo_gCH2O'),   -1.0)

    if net_gch2o > 0 and gross_gch2o > 0:
        gross_dm = gross_gch2o * NET_PHOTO_TO_DM
        net_dm   = net_gch2o   * NET_PHOTO_TO_DM
        _coupling_stats['used_RGG_3D_direct'] += 1
        return {
            'gross_photo_gDM_per_m2': gross_dm,
            'net_photo_gDM_per_m2':   net_dm,
            'lai_3d':                 glai,
            'light_interception_3d':  gli,
        }, {
            'source':    'RGG_3D_direct',
            'gross_dm':  round(gross_dm, 3),
            'net_dm':    round(net_dm, 3),
            'gli':       round(gli, 4),
        }

    # 路径 2: Python 端 RUE 回退 (老版 RGG 不发 3D photo)
    RUE_DM = 3.4
    PAR_FRAC = 0.48
    MAINT_COEF = 0.005
    GROWTH_FRAC = 0.22

    if gpar > 0:
        intercept_par = gpar
    else:
        solar = _safe_float(wx.get('irradiation'), 0.0)
        intercept_par = solar * PAR_FRAC * gli

    tavg = _safe_float(wx.get('tavg'), 25.0)
    t_opt = 28.0
    temp_resp = max(0.0, 1.0 - ((tavg - t_opt) / 15.0) ** 2)

    gross_dm = intercept_par * RUE_DM * temp_resp
    pw_g_per_m2 = _safe_float((last_rec or {}).get('plant_weight'), 0.0)
    q10f = 2.0 ** ((tavg - 25.0) / 10.0)
    maint = pw_g_per_m2 * MAINT_COEF * q10f
    net_dm = max(0.0, gross_dm - maint)
    net_dm *= (1.0 - GROWTH_FRAC)

    _coupling_stats['used_python_rebuild'] += 1
    return {
        'gross_photo_gDM_per_m2': gross_dm,
        'net_photo_gDM_per_m2':   net_dm,
        'lai_3d':                 glai,
        'light_interception_3d':  gli,
    }, {
        'source':    'python_RUE_rebuild',
        'gross_dm':  round(gross_dm, 3),
        'net_dm':    round(net_dm, 3),
        'gli':       round(gli, 4),
    }

# 5. LAI 一致性软约束 (方案A v2)
# =======================================================================
# GroIMP 有效性门槛 (数据驱动标定):
# DAS<45 GroIMP 几何冷启动期 ratio<0.90, 不可信; DAS>=45 起连续稳定可信.
# 2025年门槛待首次运行后用数据验证, 暂与2022年保持一致 (DAS=45).
NUDGE_GROIMP_START_DAY = 46    # groimp_day >= 此值才允许触发 (对应 DAS=45)
NUDGE_GROIMP_MIN_LAI   = 0.10  # GroIMP LAI 低于此值时几何结构尚未建成


def apply_wlf_nudge(
    model: 'Cotton2KDaily',
    groimp_lai_today: float,
    c2k_lai_today: float,
    groimp_day: int,
) -> Tuple[float, float, bool]:
    """
    LAI 一致性软约束 v2: 在 step_day() 之后执行, 修正量次日生效.

    修正口径:
      C2K LAI  = rec['leaf_area_index'] = W_lf x SLA x f_vis  (visible, 与实测同口径)
      GRO LAI  = groimp_lai_3d                                  (GroIMP 3D, 与实测同口径)
      两者均为当天同步值.

    有效性门槛:
      groimp_day >= 46 (DAS=45): GroIMP 几何可信起点 (GRO/C2K ratio>=0.90 连续稳定).

    参数:
        model            : Cotton2KDaily 实例 (直接修改 model.W_lf)
        groimp_lai_today : 当日 GroIMP 3D LAI (groimp_feedback['groimp_lai_3d'])
        c2k_lai_today    : 当日 C2K visible LAI (rec['leaf_area_index'])
        groimp_day       : 当前生育天 (1-based)

    返回:
        (nudge_wlf_delta, c2k_lai_used, fired)
    """
    # 门槛1: GroIMP 几何建成期 (DAS<45) 不可信
    if groimp_day < NUDGE_GROIMP_START_DAY:
        return 0.0, c2k_lai_today, False

    # 门槛2: GroIMP LAI 过低表示几何尚未建成
    if groimp_lai_today < NUDGE_GROIMP_MIN_LAI:
        return 0.0, c2k_lai_today, False

    # 门槛3: C2K LAI 过低无意义
    if c2k_lai_today < 0.05:
        return 0.0, c2k_lai_today, False

    # 相对偏差 (均为 visible 口径)
    max_lai = max(groimp_lai_today, c2k_lai_today)
    rel_diff = abs(groimp_lai_today - c2k_lai_today) / max_lai
    if rel_diff <= LAI_NUDGE_THRESHOLD:
        return 0.0, c2k_lai_today, False

    # DVS 衰减因子 (与 RGG 端 c2kNudgeLAI() 完全一致)
    if groimp_day <= 80:
        dvs_factor = 1.0
    elif groimp_day <= 110:
        dvs_factor = 1.0 - (groimp_day - 80) / 30.0 * 0.90
    else:
        dvs_factor = 0.10

    # 计算当日 SLA x f_vis (visible SLA)
    # das: step_day 执行后 current_day 已推进一天, 需回退一天
    das = max((model.current_day - model.emerge_row) - 1, 0) \
        if model.emerge_row is not None else 0
    theta_root = (
        0.20 * model._soil_state.get('theta_10', SOIL_THETA_FC) +
        0.30 * model._soil_state.get('theta_20', SOIL_THETA_FC) +
        0.50 * model._soil_state.get('theta_27', SOIL_THETA_FC)
    )
    ks_est = float(np.clip(
        (theta_root - SOIL_THETA_WP) / (SOIL_THETA_FC - SOIL_THETA_WP),
        0.5, 1.0
    ))
    sla    = get_dynamic_sla(das, ks_est, 1.0)
    f_vis  = get_leaf_visibility_factor(float(das))
    sla_vis = sla * f_vis

    if sla_vis < 1e-6:
        return 0.0, c2k_lai_today, False

    # 目标 W_lf: 使得次日 visible LAI ≈ GroIMP LAI
    target_wlf = groimp_lai_today / sla_vis
    gap = target_wlf - model.W_lf

    daily_rate = 1.0 - math.exp(-1.0 / LAI_NUDGE_TAU)
    nudge = gap * daily_rate * dvs_factor

    # 限幅: 单日不超过 gap 的 MAX_RATE
    max_nudge = abs(gap) * LAI_NUDGE_MAX_RATE * dvs_factor
    nudge = max(-max_nudge, min(max_nudge, nudge))

    # 应用修正; W_lf 下限 0.5 g/m²
    model.W_lf = max(0.5, model.W_lf + nudge)

    return float(nudge), c2k_lai_today, True


# 6. Cotton2K → GroIMP 数据准备
# =======================================================================
def prepare_c2k_data_for_groimp(
    rec: Dict[str, Any],
    groimp_day: int,
    wx: Dict[str, Any],
) -> Dict[str, Any]:
    """
    构造发往 GroIMP 的当日 JSON 字段.

    设计原则:
      1. 字段名向后兼容 RGG (避免 RGG 端 parseJSONFloat 读到 -1).
      2. 单位与 RGG 一致: 干重 g/m², LI 0-1, 产量 kg/ha.
      3. 慢变字段每 5 天发送 1 次 (减带宽).
      4. lint_yield 单位严格为 kg/ha; 回填用 boll_weight × 衣分 × 10.
      5. c2k_shed_factor 字段保留以兼容 RGG, 语义为诊断值 ∈ [0,1].
    """
    def rr(key: str, nd: int = 2, default: float = 0.0) -> float:
        return round(_safe_float(rec.get(key, default), default), nd)

    c2k_net_raw = rr('c2k_net_photo', 4)
    c2k_net_dm  = round(c2k_net_raw * NET_PHOTO_TO_DM, 4)

    data: Dict[str, Any] = {
        'groimp_day':            rec.get('groimp_day', groimp_day),
        'date':                  rec.get('date', ''),
        'doy':                   rec.get('doy', 0),

        'T_mean':                rr('T_mean', 1),
        'GDD_daily':             rr('GDD_daily', 2),
        'GDD_cumsum':            rr('GDD_cumsum', 1),
        'GDD_after_emerge':      rr('GDD_after_emerge', 1),

        'plant_height':          rr('plant_height', 2),
        'leaf_area_index':       rr('leaf_area_index', 4),
        'light_interception':    rr('light_interception', 4),

        'leaf_weight':           rr('leaf_weight', 2),
        'stem_weight':           rr('stem_weight', 2),
        'boll_weight':           rr('boll_weight', 2),
        'plant_weight':          rr('plant_weight', 2),
        'root_weight':           rr('root_weight', 2),

        'number_of_squares':     rr('number_of_squares', 1),
        'number_of_green_bolls': rr('number_of_green_bolls', 1),
        'number_of_open_bolls':  rr('number_of_open_bolls', 1),

        'c2k_sourceSink':        round(_safe_float(rec.get('c2k_source_sink'), -1.0), 3),
        'c2k_dvs':               round(_safe_float(rec.get('c2k_dvs'), 0.0), 3),
        'c2k_shedFactor':        round(_safe_float(rec.get('c2k_shed_fraction_diag'), 0.0), 4),
        'c2k_shed_fraction_diag':round(_safe_float(rec.get('c2k_shed_fraction_diag'), 0.0), 4),
        'c2k_devStage':          rec.get('c2k_dev_stage', 'seedling'),
        'c2k_devProgress':       round(_safe_float(rec.get('c2k_dev_progress'), 0.0), 3),

        'c2k_netPhoto':          c2k_net_raw,
        'c2k_netPhoto_gDM':      c2k_net_dm,
        'c2k_netPhoto_unit':     NET_PHOTO_UNIT,
        'c2k_grossPhoto':        rr('c2k_gross_photo', 3),

        'c2k_sla_m2g':           round(_safe_float(rec.get('c2k_sla_m2g'), 0.013), 4),
        'c2k_leaf_area_m2':      rr('c2k_leaf_area_m2', 4),

        'main_stem_nodes':       rr('main_stem_nodes', 1),
        'fruit_branch_number':   rr('fruit_branch_number', 1),
        'single_leaf_area':      rr('single_leaf_area', 2),
        'boll_diameter':         rr('boll_diameter', 2),
        'leaf_angle':            rr('leaf_angle', 1),

        'tavg':          round(wx.get('tavg', rec.get('T_mean', 0)), 1),
        'rh_avg':        round(wx.get('rh_avg', 0.0), 4),
        'tdew':          round(wx.get('tdew', 0.0), 1),
        'vpd':           round(wx.get('vpd', 0.0), 4),
        'sunshine_frac': round(wx.get('sunshine_frac', 0.0), 4),
        'daylight_sec':  wx.get('daylight_sec', 0.0),

        'coupling_source': rec.get('coupling_source', 'unknown'),
        # nudge_enabled: 方案A v2 由 Python 端执行, 此字段仍发送以保持 RGG 端诊断一致性
        'nudge_enabled':   NUDGE_ENABLED,
    }

    # 慢变字段 (每 5 天 + 前 5 天发送)
    if groimp_day % 5 == 0 or groimp_day <= 5:
        data.update({
            'leaf_length':       rr('leaf_length', 2),
            'leaf_width':        rr('leaf_width', 2),
            'internode_length':  rr('internode_length', 2),
            'petiole_weight':    rr('petiole_weight', 2),
            'seed_cotton_yield': rr('seed_cotton_yield', 2),
            'square_weight':     rr('square_weight', 2),
        })

        # lint_yield 单位严格为 kg/ha
        # boll_weight 单位是 g/m², 换算: 1 g/m² × 10 = 10 kg/ha
        lint_direct_kgha = rr('lint_yield', 2)
        if lint_direct_kgha <= 0:
            boll_w_g_per_m2 = rr('boll_weight', 3)
            if boll_w_g_per_m2 > 0:
                data['lint_yield'] = round(
                    boll_w_g_per_m2 * LINT_RATIO_DEFAULT * 10.0, 2
                )
                data['lint_yield_estimated'] = True
                data['lint_yield_unit_note'] = 'estimated_from_boll_weight_x_lint_x_10'
            else:
                data['lint_yield'] = 0.0
                data['lint_yield_estimated'] = False
        else:
            data['lint_yield'] = lint_direct_kgha
            data['lint_yield_estimated'] = False

    return data

# 6. 评估指标 (NSE / r / r² / RMSE / MBE / Willmott d)
# =======================================================================
def model_metrics(obs: List[float], sim: List[float]) -> Dict[str, float]:
    """
    SCI 评估指标: NSE / r / r² / RMSE / MBE / Willmott d

    三个指标的严格定义 (Moriasi et al. 2007; Willmott & Matsuura 2005):
      r   (Pearson相关系数):  衡量模拟值与观测值的线性趋势同步程度 [−1, 1]
      r²  (Pearson², 趋势一致性): 仅表征趋势协变, 非预测精度 [0, 1]
                                   注: r² ≠ CoD; 论文中禁止混用
      NSE (Nash-Sutcliffe效率, 即模型评估CoD):
          NSE = 1 - Σ(obs-sim)² / Σ(obs-obs_mean)²  [可为负]
          当存在系统偏差时 NSE << r², 两者差距正是偏差量的定量体现.
    """
    if len(obs) < 3 or len(obs) != len(sim):
        return {}
    o = np.asarray(obs, dtype=float)
    s = np.asarray(sim, dtype=float)
    n = len(o)
    o_mean = float(np.mean(o))
    ss_res = float(np.sum((o - s) ** 2))
    ss_tot = float(np.sum((o - o_mean) ** 2))
    rmse = math.sqrt(ss_res / n)
    mbe = float(np.mean(s - o))
    nrmse = rmse / o_mean * 100.0 if o_mean != 0 else float('nan')
    nse = (1.0 - ss_res / ss_tot) if ss_tot > 0 else float('nan')
    if np.std(o) > 0 and np.std(s) > 0:
        r_pearson = float(np.corrcoef(o, s)[0, 1])
        r_sq = r_pearson ** 2
    else:
        r_pearson = float('nan')
        r_sq = float('nan')
    denom_d = float(np.sum((np.abs(s - o_mean) + np.abs(o - o_mean)) ** 2))
    d_idx = (1.0 - ss_res / denom_d) if denom_d > 0 else float('nan')
    return {
        'n': n, 'obs_mean': o_mean,
        'RMSE': rmse, 'nRMSE_%': nrmse, 'MBE': mbe,
        'r_sq': r_sq,
        'NSE': nse,
        'd_index': d_idx, 'r_Pearson': r_pearson,
    }

def rate_nrmse(nrmse: float) -> str:
    """Jamieson et al. (1991) 分级."""
    if math.isnan(nrmse):
        return 'N/A'
    if nrmse < 10: return 'Excellent'
    if nrmse < 20: return 'Good'
    if nrmse < 30: return 'Fair'
    return 'Poor'

# 7. 主耦合循环
# =======================================================================
def run_daily_coupling(
    weather_path: str,
    host: str,
    port: int,
    output_dir: str,
    handshake_timeout_s: float = 60.0,
    nudge_enabled: bool = True,
) -> None:
    _coupling_stats['attempts'] = 0
    _coupling_stats['used_RGG_3D_direct'] = 0
    _coupling_stats['used_python_rebuild'] = 0
    _coupling_stats['fallback_invalid'] = 0

    global NUDGE_ENABLED
    NUDGE_ENABLED = nudge_enabled

    print('=' * 70)
    print(f'  cotton2k与GroIMP耦合 — 逐日耦合服务器 [2025独立验证年, 参数冻结自2022标定年]')
    nudge_tag = ('NudgeON  (方案A v2: visible口径+step后同步比对+DAS>=45门槛, THRESHOLD=15% tau=6d MAX=20%)'
                 if NUDGE_ENABLED else
                 'NudgeOFF ★ 消融实验 — LAI一致性软约束关闭')
    print(f'  Nudge 模式: {nudge_tag}')
    print(f'  新疆阿拉尔 | 塔河2号 (审定衣分 42.0%)')
    print(f'  制式: 一膜六行 膜宽{FILM_WIDTH_CM:.0f}cm | 行距 {WIDE_ROW_CM:.0f}+{NARROW_ROW_CM:.0f}cm | 株距 {PLANT_SPACING_CM:.0f}cm')
    print(f'  密度: {PLANTS_PER_HA:,} 株/hm² ({PLANTS_PER_MU:,} 株/亩, {PLANTS_PER_HA/10000:.0f} 株/m²)')
    print(f'  播种: 2025-04-23  出苗: 2025-05-02  打顶: 2025-07-16 (播后85天, DAE=75)')
    print('=' * 70)
    print(f'  LAI 目标: 盛铃期 {LAI_PEAK_TARGET[0]:.1f}-{LAI_PEAK_TARGET[1]:.1f} (峰值约 4.0)')
    print(f'  物候 GDD 阈值 (出苗后): 初蕾=345, 初花=630, 吐絮始=1510, 成熟=1620')
    print(f'  实测峰LAI: {MEASURED_LAI_PEAK:.3f}±{MEASURED_LAI_PEAK_SD:.3f} (DAS={MEASURED_LAI_PEAK_DAS})')
    print(f'  实测籽棉: {MEASURED_SEED_KGHA:.0f} kg/ha → 皮棉: {MEASURED_LINT_KGHA:.0f} kg/ha (衣分42%)')
    print('=' * 70)
    print(f'  [土壤] 根区桶深 = 300 mm (膜下滴灌 0-30cm 主根区, 两年统一)')
    print(f'  [T2]  INERTIA_COEF = {INERTIA_COEF:.5f}  ← 产量唯一可标定参数')
    print(f'  [消融] Nudge OFF 模式 — LAI 一致性软约束已禁用 (THRESHOLD=999)')
    if abs(INERTIA_COEF - 0.00062) < 1e-7:
        print(f'  [!]  INERTIA_COEF 仍为 100mm 旧值; 请先重标定再对比消融结果')
    print('=' * 70)

    if not os.path.exists(weather_path):
        raise FileNotFoundError(f"气象文件不存在: {weather_path}")
    print(f'\n[气象] 加载: {weather_path}')
    norm_path, weather_index = preprocess_weather(weather_path)
    print(f'[气象] 归一化完成, {len(weather_index)} 天')
    if not weather_index:
        raise RuntimeError('气象文件无有效数据')

    vpd_series = [v.get('vpd', 0) for v in weather_index.values() if v]
    sf_series  = [v.get('sunshine_frac', 0) for v in weather_index.values() if v]
    if vpd_series:
        print(f'[气象] 季均 VPD = {sum(vpd_series)/len(vpd_series):.3f} kPa, '
              f'季均日照率 = {sum(sf_series)/len(sf_series):.1%}')

    print(f'\n[C2K] 加载模型: {norm_path}')
    model = Cotton2KDaily(norm_path)
    print(f'[C2K] total_days = {model.total_days}, emerge_row = {model.emerge_row}')

    # ── Socket 监听 ───────────────────────────────────────────────
    server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    try:
        server.bind((host, port))
    except OSError as e:
        raise RuntimeError(f'端口 {host}:{port} 占用: {e}')
    server.listen(1)
    server.settimeout(handshake_timeout_s)
    print(f'\n[Socket] 监听 {host}:{port} (超时 {handshake_timeout_s}s)')
    print('[Socket] 等待 GroIMP 连接 (在 GroIMP 中执行 init())...')

    conn: Optional[socket.socket] = None
    try:
        conn, addr = server.accept()
        conn.settimeout(None)
        print(f'[Socket] ✓ 已连接: {addr}')

        read_line = make_line_reader(conn)
        try:
            handshake = read_line()
            print(f'[Socket] 握手: {handshake.get("msg", "OK")}')
        except Exception as e:
            print(f'[Socket] 握手失败: {e}')
            return

        # === 阶段 A: 出苗前 (Cotton2K 独立) ===
        print(f'\n{"="*60}\n  阶段 A: 播种 → 出苗 (C2K 独立 {model.emerge_row} 天)\n{"="*60}')
        for _ in range(model.emerge_row):
            model.step_day(coupling_inputs=None)
        print(f'  完成 {model.emerge_row} 天')

        # === 阶段 B: 双向耦合 142 天 ===
        print(f'\n{"="*60}\n  阶段 B: 出苗 → 吐絮始 ({GROWTH_PERIOD} 天前置注入耦合)\n{"="*60}')
        print(f"  {'Day':>4}|{'C2K LAI':>8}|{'GRO LAI':>8}|{'Δ%':>6}|"
              f"{'VPD':>5}|{'Stage':>12}|{'Source':>14}|{'BollW':>6}|{'ms':>5}")
        print(f"  {'-'*4}+{'-'*8}+{'-'*8}+{'-'*6}+{'-'*5}+{'-'*12}+{'-'*14}+{'-'*6}+{'-'*5}")

        coupling_log: List[Dict[str, Any]] = []
        groimp_feedback: Optional[Dict[str, Any]] = None
        last_rec: Dict[str, Any] = {}
        # Nudge 追踪 (方案A v2: Python 端 W_lf 修正)
        nudge_log: List[Dict[str, Any]] = []
        nudge_total_delta = 0.0

        for groimp_day in range(1, GROWTH_PERIOD + 1):
            t0 = time.time()

            cur_date_idx = model.current_day
            if cur_date_idx >= len(model.weather):
                print(f'  [警告] Cotton2K 已耗尽气象数据, day={groimp_day}')
                break
            cur_date = model.weather.iloc[cur_date_idx]['date'].strftime('%Y-%m-%d')
            wx = weather_index.get(cur_date, {})

            # 真·前置注入: 上一日 GroIMP 反馈 → coupling_inputs → step_day
            coupling, cinp_report = build_coupling_inputs(groimp_feedback, last_rec, wx)
            rec = model.step_day(coupling_inputs=coupling)
            if rec is None:
                print(f'  [警告] step_day 返回 None, day={groimp_day}')
                break

            c2k_data = prepare_c2k_data_for_groimp(rec, groimp_day, wx)
            try:
                send_json(conn, {'day': groimp_day, 'data': c2k_data})
            except (OSError, BrokenPipeError) as e:
                print(f'\n  [错误] 发送失败 day={groimp_day}: {e}')
                break

            try:
                response = read_line()
            except ConnectionError as e:
                print(f'\n  [错误] GroIMP 断开 day={groimp_day}: {e}')
                break
            groimp_feedback = response.get('feedback', {}) if response else {}

            elapsed_ms = (time.time() - t0) * 1000.0

            groimp_lai = _safe_float(groimp_feedback.get('groimp_lai_3d'))
            groimp_li  = _safe_float(groimp_feedback.get('groimp_light_interception'))
            c2k_lai = _safe_float(rec.get('leaf_area_index'))
            ss      = _safe_float(rec.get('c2k_source_sink'))

            lai_dev_pct = 0.0
            if c2k_lai > 0.01 and groimp_lai > 0.01:
                lai_dev_pct = (groimp_lai - c2k_lai) / max(c2k_lai, groimp_lai) * 100.0

            # ── 方案A v2 Nudge: 在 step_day() 后执行, 修正次日生效 ──────
            nudge_delta = 0.0
            nudge_fired = False
            if NUDGE_ENABLED and groimp_lai >= 0.0 and c2k_lai >= 0.0:
                nudge_delta, _c2k_used, nudge_fired = apply_wlf_nudge(
                    model, groimp_lai, c2k_lai, groimp_day
                )
                if nudge_fired:
                    nudge_total_delta += nudge_delta
                    nudge_log.append({
                        'groimp_day':       groimp_day,
                        'date':             rec.get('date', ''),
                        'groimp_lai_today': round(groimp_lai, 4),
                        'c2k_lai_today':    round(c2k_lai, 4),
                        'rel_diff_pct':     round(abs(groimp_lai - c2k_lai) /
                                                  max(groimp_lai, c2k_lai) * 100, 2),
                        'nudge_wlf_delta':  round(nudge_delta, 4),
                        'wlf_post':         round(model.W_lf, 4),
                        'c2k_dev_stage':    rec.get('c2k_dev_stage', ''),
                    })
            # ─────────────────────────────────────────────────────────────

            log_entry = {
                'date': rec.get('date', ''),
                'groimp_day': groimp_day,
                'c2k_lai': round(c2k_lai, 4),
                'groimp_lai': round(groimp_lai, 4),
                'lai_dev_pct': round(lai_dev_pct, 2),
                'c2k_li': round(_safe_float(rec.get('light_interception')), 4),
                'groimp_li': round(groimp_li, 4),
                'groimp_n_leaves': int(_safe_float(groimp_feedback.get('groimp_n_leaves'))),
                'c2k_source_sink': round(ss, 4),
                'c2k_dvs': round(_safe_float(rec.get('c2k_dvs')), 4),
                'c2k_dev_stage': rec.get('c2k_dev_stage', ''),
                'groimp_n_bolls':          int(_safe_float(groimp_feedback.get('groimp_n_bolls'))),
                'groimp_n_open_bolls':     int(_safe_float(groimp_feedback.get('groimp_n_open_bolls'))),
                'groimp_n_eventual_bolls': int(_safe_float(groimp_feedback.get('groimp_n_eventual_bolls'))),
                'groimp_n_squares':        int(_safe_float(groimp_feedback.get('groimp_n_squares'))),
                'groimp_n_fruit_branches': int(_safe_float(groimp_feedback.get('groimp_n_fruit_branches'))),
                'groimp_avg_boll_wt':      round(_safe_float(groimp_feedback.get('groimp_avg_boll_weight')), 2),
                'groimp_total_boll_wt':    round(_safe_float(groimp_feedback.get('groimp_total_boll_weight')), 1),
                'groimp_total_boll_wt_all':  round(_safe_float(
                    groimp_feedback.get('groimp_total_boll_weight_all',
                                         groimp_feedback.get('groimp_total_boll_weight', 0))), 2),
                'groimp_total_boll_wt_open': round(_safe_float(
                    groimp_feedback.get('groimp_total_boll_weight_open', 0)), 2),
                'groimp_avg_boll_wt_all':    round(_safe_float(
                    groimp_feedback.get('groimp_avg_boll_weight_all',
                                         groimp_feedback.get('groimp_avg_boll_weight', 0))), 2),
                'groimp_avg_boll_wt_open':   round(_safe_float(
                    groimp_feedback.get('groimp_avg_boll_weight_open', 0)), 2),
                'groimp_current_green_brown_wt': round(_safe_float(
                    groimp_feedback.get('groimp_current_green_brown_weight', 0)), 2),
                'groimp_shed_today':       int(_safe_float(groimp_feedback.get('groimp_n_shed_today'))),
                'groimp_shed_total':       int(_safe_float(groimp_feedback.get('groimp_n_shed_total'))),
                'groimp_sq_shed':          int(_safe_float(groimp_feedback.get('groimp_n_square_shed'))),
                'groimp_bl_shed':          int(_safe_float(groimp_feedback.get('groimp_n_boll_shed'))),
                'groimp_pos0_shed':        int(_safe_float(groimp_feedback.get('groimp_n_pos0_shed'))),
                'groimp_pos1_shed':        int(_safe_float(groimp_feedback.get('groimp_n_pos1_shed'))),
                'groimp_max_fruit_branches': int(_safe_float(groimp_feedback.get('groimp_max_fruit_branches', 0))),
                'groimp_max_fruit_pos':      int(_safe_float(groimp_feedback.get('groimp_max_fruit_pos', 1))),
                'tavg': round(wx.get('tavg', 0), 1),
                'rh_avg': round(wx.get('rh_avg', 0), 4),
                'vpd': round(wx.get('vpd', 0), 4),
                'sunshine_frac': round(wx.get('sunshine_frac', 0), 4),
                'c2k_netPhoto_gCH2O': c2k_data.get('c2k_netPhoto', 0.0),
                'c2k_netPhoto_gDM':   c2k_data.get('c2k_netPhoto_gDM', 0.0),
                'coupling_source':    rec.get('coupling_source', 'unknown'),
                'coupling_input_used':bool(coupling is not None),
                'coupling_path':      cinp_report.get('source', 'none'),
                'coupling_gross_dm':  cinp_report.get('gross_dm', 0.0),
                'coupling_net_dm':    cinp_report.get('net_dm', 0.0),
                'nudge_wlf_delta':    round(nudge_delta, 4),
                'nudge_fired':        nudge_fired,
                'elapsed_ms':         round(elapsed_ms, 1),
            }
            coupling_log.append(log_entry)
            last_rec = {k: v for k, v in rec.items() if not k.startswith('_')}

            if groimp_day % 10 == 0 or groimp_day <= 5 or groimp_day >= GROWTH_PERIOD - 1:
                stage = rec.get('c2k_dev_stage', '')
                src   = rec.get('coupling_source', '?')[:13]
                bollw = _safe_float(groimp_feedback.get('groimp_avg_boll_weight'))
                vpd   = wx.get('vpd', 0.0)
                print(f'  {groimp_day:>4}|{c2k_lai:>8.3f}|{groimp_lai:>8.3f}|'
                      f'{lai_dev_pct:>+5.1f}%|{vpd:>5.3f}|'
                      f'{stage:>12}|{src:>14}|{bollw:>5.1f}g|{elapsed_ms:>5.0f}')

        # === 阶段 B 结束: 冻结 GroIMP 铃数 ===
        # 设计原则: freeze 只锁定形态学量 (n_open_bolls, n_green_bolls),
        # 不修改碳平衡量 (W_bl/W_lint/W_seed).
        # 阶段 C 41 天 Cotton2K 继续基于真实碳平衡积累,
        # 路径1 反映真实可达产量, 路径3 反映理论最大潜力.
        # 收获时 (DAS 200) 全部 8 铃吐絮, n_green_bolls 强制为 0.
        # 铃数取法: 取 142 天中 n_eventual_bolls 的峰值,
        # 避免 L-system 末端物候推进导致部分铃被标记为已收获.
        if coupling_log:
            last_log = coupling_log[-1]

            peak_open = 0.0
            peak_eventual = 0.0
            peak_idx  = len(coupling_log) - 1
            for idx, entry in enumerate(coupling_log):
                v  = _safe_float(entry.get('groimp_n_open_bolls'))
                ev = _safe_float(entry.get('groimp_n_eventual_bolls', 0))
                if ev > peak_eventual:
                    peak_eventual = ev
                    peak_idx = idx
                if v > peak_open:
                    peak_open = v
            peak_log = coupling_log[peak_idx]

            n_eventual_groimp = peak_eventual if peak_eventual > 0 else peak_open
            n_open_at_d142    = _safe_float(last_log.get('groimp_n_open_bolls'))
            n_green_at_d142   = _safe_float(last_log.get('groimp_n_bolls'))

            peak_boll_wt_all   = _safe_float(peak_log.get('groimp_total_boll_wt_all'))
            last_boll_wt_all   = _safe_float(last_log.get('groimp_total_boll_wt_all'))
            if peak_boll_wt_all <= 0 and last_boll_wt_all <= 0:
                peak_boll_wt_all = _safe_float(peak_log.get('groimp_total_boll_wt'))
                last_boll_wt_all = _safe_float(last_log.get('groimp_total_boll_wt'))
            total_boll_wt_per_plant = max(peak_boll_wt_all, last_boll_wt_all)
            total_boll_wt_per_m2    = total_boll_wt_per_plant * (PLANTS_PER_HA / HA_TO_M2)

            print(f'\n{"="*70}')
            print(f'  阶段 B → C 切换: 冻结到收获日 (DAS 200) 全开状态')
            print(f'{"="*70}')
            print(f'  142 天末态 (DAS 142):')
            print(f'    已吐絮铃:  {n_open_at_d142:.2f}/株')
            print(f'    绿/褐铃:   {n_green_at_d142:.2f}/株')
            print(f'    全铃合计:  {n_open_at_d142+n_green_at_d142:.2f}/株')
            print(f'    全铃重量:  {last_boll_wt_all:.2f} g/株 (含绿+褐+开)')
            print(f'  Eventual 峰值 (day {peak_log.get("groimp_day","?")}):'
                  f' {peak_eventual:.2f}/株')
            print(f'')
            print(f'  冻结状态 (反映 DAS 200 收获日全开):')
            print(f'    n_open_bolls (传入 freeze) = {n_eventual_groimp:.2f} 个/株 (全铃数)')
            print(f'    n_green_bolls (强制)       = 0.00 个/株 (收获日全部已开)')
            print(f'    全铃总重 (诊断)            = {total_boll_wt_per_m2:.1f} g/m²')
            print(f'  ✓ freeze 不覆盖 W_bl/W_lint/W_seed (碳平衡自治)')

            model.freeze_reproductive_state(
                n_open_bolls=n_eventual_groimp,
                n_green_bolls=0.0,
                total_boll_weight_g_per_m2=total_boll_wt_per_m2,
            )

            # 设置 W_lint 结构封顶 (确保路径1 ≤ 路径3)
            # 使用全铃均重 (groimp_avg_boll_wt_all)
            if hasattr(model, 'set_lint_cap_from_structure'):
                boll_wts_valid = [e['groimp_avg_boll_wt_all'] for e in coupling_log
                                  if e.get('groimp_avg_boll_wt_all', 0) > 0]
                if not boll_wts_valid:
                    boll_wts_valid = [e['groimp_avg_boll_wt'] for e in coupling_log
                                      if e.get('groimp_avg_boll_wt', 0) > 0]
                avg_boll_for_cap = (np.mean(boll_wts_valid[-5:])
                                    if len(boll_wts_valid) >= 5
                                    else (boll_wts_valid[-1] if boll_wts_valid else 6.3))
                model.set_lint_cap_from_structure(
                    n_open_per_plant=n_eventual_groimp,
                    avg_boll_weight_g=float(avg_boll_for_cap),
                )

            print(f'  ✓ Cotton2K 阶段 C 基于真实碳平衡继续积累干物质')
            print(f'    路径1 = Cotton2K W_lint (真实可达产量)')
            print(f'    路径3 = 全铃 × 全铃均重 × 衣分 × 密度 (理论上限)')
            print(f'    收获效率 = 路径1/路径3, 文献 60-80% (Bange & Milroy 2004)')

        # === 阶段 C: 收获期 (Cotton2K 独立 + repro 冻结) ===
        print(f'\n{"="*60}\n  阶段 C: 收获期 (C2K 独立 + repro 冻结)\n{"="*60}')
        remaining = max(0, model.total_days - model.current_day)
        for _ in range(remaining):
            r = model.step_day(coupling_inputs=None)
            if r is None:
                break
        print(f'  完成 {remaining} 天')

        try:
            send_json(conn, {'day': -1, 'cmd': 'FINISH'})
        except Exception:
            pass

        save_results(model, coupling_log, output_dir,
                     nudge_log=nudge_log, nudge_total_delta=nudge_total_delta)

    finally:
        if conn is not None:
            try: conn.close()
            except Exception: pass
        try: server.close()
        except Exception: pass
        try: os.unlink(norm_path)
        except Exception: pass

# 8. 结果保存 + SCI 风格统计报告
# =======================================================================
def save_results(model, coupling_log: List[Dict[str, Any]], output_dir: str,
                 nudge_log: Optional[List[Dict[str, Any]]] = None,
                 nudge_total_delta: float = 0.0) -> None:
    os.makedirs(output_dir, exist_ok=True)
    suffix = 'ON' if NUDGE_ENABLED else 'OFF'

    df = model.get_all_records()
    out_csv = os.path.join(output_dir, f'output_nudge_{suffix}.csv')
    df.to_csv(out_csv, index=False)
    print(f'\n  [输出] {out_csv}  ({len(df)} × {len(df.columns)})')

    log_csv = os.path.join(output_dir, f'coupling_log_{suffix}.csv')
    if coupling_log:
        with open(log_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=coupling_log[0].keys())
            writer.writeheader()
            writer.writerows(coupling_log)
        print(f'  [输出] {log_csv}  ({len(coupling_log)} 天)')

    # Nudge 明细日志
    if nudge_log:
        nudge_csv = os.path.join(output_dir, f'nudge_detail_log_{suffix}.csv')
        with open(nudge_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=nudge_log[0].keys())
            writer.writeheader()
            writer.writerows(nudge_log)
        print(f'  [输出] {nudge_csv}  ({len(nudge_log)} 次触发)')

    if not coupling_log:
        return

    print(f'\n{"="*70}')
    print(f'  cotton2k与GroIMP耦合 统计报告')
    print(f'{"="*70}')

    valid = [e for e in coupling_log
             if e['groimp_lai'] > 0.01 and e['c2k_lai'] > 0.01]

    # ── A. LAI 全期与分阶段评估 ──────────────────────────────
    if valid:
        obs = [e['c2k_lai']    for e in valid]
        sim = [e['groimp_lai'] for e in valid]
        m = model_metrics(obs, sim)

        c2k_peak = max(obs) if obs else 0.0
        gro_peak = max(sim) if sim else 0.0
        peak_in_target_c2k = LAI_PEAK_TARGET[0] <= c2k_peak <= LAI_PEAK_TARGET[1]
        peak_in_target_gro = LAI_PEAK_TARGET[0] <= gro_peak <= LAI_PEAK_TARGET[1]
        print(f'\n  === LAI 峰值 vs 田间实测 ===')
        print(f'  实测峰 LAI:        {MEASURED_LAI_PEAK:.3f} ± {MEASURED_LAI_PEAK_SD:.3f}'
              f'  (DAS={MEASURED_LAI_PEAK_DAS}, 4区均值)')
        print(f'  Cotton2K 模拟峰:    {c2k_peak:.3f}    误差 {(c2k_peak-MEASURED_LAI_PEAK)/MEASURED_LAI_PEAK*100:+.1f}%'
              f'  {"✓" if peak_in_target_c2k else "⚠"} (目标 {LAI_PEAK_TARGET[0]:.1f}-{LAI_PEAK_TARGET[1]:.1f})')
        print(f'  GroIMP 模拟峰:     {gro_peak:.3f}    误差 {(gro_peak-MEASURED_LAI_PEAK)/MEASURED_LAI_PEAK*100:+.1f}%'
              f'  {"✓" if peak_in_target_gro else "⚠"} (目标 {LAI_PEAK_TARGET[0]:.1f}-{LAI_PEAK_TARGET[1]:.1f})')

        print(f'\n  === 全期 LAI: GroIMP vs Cotton2K (n={m.get("n", 0)}) ===')
        if m:
            print(f'  ── 趋势相关 (Pearson) ──────────────────────────────────')
            print(f'  r  (Pearson相关系数):  {m["r_Pearson"]:>8.4f}   ← 论文主引, 衡量趋势同步')
            print(f'  r² (Pearson², 趋势):   {m["r_sq"]:>8.4f}   ← 仅表征趋势; ≠ 决定系数(CoD)')
            print(f'  ── 预测精度 (模型评估CoD) ──────────────────────────────')
            print(f'  NSE (=CoD, Nash-Sutcliffe): {m["NSE"]:>7.4f}   ← 论文主引, 衡量预测偏差')
            print(f'  Willmott d:           {m["d_index"]:>8.4f}')
            print(f'  ── 误差量 ──────────────────────────────────────────────')
            print(f'  RMSE:                 {m["RMSE"]:>8.4f}')
            print(f'  nRMSE:                {m["nRMSE_%"]:>7.2f}%   [{rate_nrmse(m["nRMSE_%"])}]')
            print(f'  MBE:                  {m["MBE"]:>+8.4f}   '
                  f'{"(GroIMP > C2K)" if m["MBE"]>0 else "(GroIMP < C2K)"}')
            print(f'  [!] r²≠CoD; CoD=NSE={m["NSE"]:.4f}; 论文中禁止将 r² 标注为 R²(决定系数)')

        print(f'\n  === 分阶段 LAI 评估 (Jamieson et al. 1991 分级) ===')
        print(f"  {'阶段':<14}|{'n':>4}|{'r²(P)':>7}|{'RMSE':>7}|{'MBE':>8}|"
              f"{'nRMSE':>7}|{'d':>6}|{'等级':>10}")
        print(f"  {'-'*14}+{'-'*4}+{'-'*7}+{'-'*7}+{'-'*8}+{'-'*7}+{'-'*6}+{'-'*10}")
        for pname, (d0, d1) in PHASES.items():
            pv = [e for e in valid if d0 <= e['groimp_day'] <= d1]
            if len(pv) < 3:
                continue
            po = [e['c2k_lai']    for e in pv]
            ps = [e['groimp_lai'] for e in pv]
            pm = model_metrics(po, ps)
            if not pm:
                continue
            print(f'  {pname:<14}|{pm["n"]:>4}|{pm["r_sq"]:>+7.3f}|'
                  f'{pm["RMSE"]:>7.3f}|{pm["MBE"]:>+8.3f}|'
                  f'{pm["nRMSE_%"]:>6.1f}%|{pm["d_index"]:>6.3f}|'
                  f'{rate_nrmse(pm["nRMSE_%"]):>10}')

    # ── B. 气象统计 ──────────────────────────────────────
    vpds = [e['vpd'] for e in coupling_log if e.get('vpd', 0) > 0]
    sunfs = [e['sunshine_frac'] for e in coupling_log if e.get('sunshine_frac', -1) >= 0]
    if vpds:
        print(f'\n  === 气象驱动 ===')
        print(f'  VPD (FAO-56) 均值/最大: {np.mean(vpds):.3f} / {max(vpds):.3f} kPa')
        print(f'  日照率 均值/最小:        {np.mean(sunfs):.1%} / {min(sunfs):.1%}')

    # ── C. 铃数统计 ──────────────────────────────────────
    last = coupling_log[-1]
    n_bl   = int(last.get('groimp_n_bolls', 0))
    n_open = int(last.get('groimp_n_open_bolls', 0))
    n_eventual = n_bl + n_open
    max_branches = int(last.get('groimp_max_fruit_branches', 0))
    max_pos = int(last.get('groimp_max_fruit_pos', 1))

    print(f'\n  === 铃数统计 (脱落模块已删除, 确定性 8 铃/株) ===')
    print(f'  RGG 参数: MAX_FRUIT_BRANCHES={max_branches}, MAX_FRUIT_POS={max_pos}')
    print(f'  DAS 142 末态: 吐絮 {n_open} 铃 + 绿/褐 {n_bl} 铃 = {n_eventual} 铃/株')
    print(f'  目标: {BOLL_COUNT_TARGET:.1f} 铃/株 (田间均值)')
    print(f'  收获日 (DAS 200): 全部 {n_eventual} 铃吐絮')

    # ── D. 铃重统计 (全铃口径) ───────────────────────────
    boll_wts_all_log = [e['groimp_avg_boll_wt_all'] for e in coupling_log
                        if e.get('groimp_avg_boll_wt_all', 0) > 0]
    boll_wts_open_log = [e['groimp_avg_boll_wt'] for e in coupling_log
                         if e['groimp_avg_boll_wt'] > 0]
    boll_wts = boll_wts_all_log if boll_wts_all_log else boll_wts_open_log

    total_boll_wt_all   = float(last.get('groimp_total_boll_wt_all',
                                          last.get('groimp_total_boll_wt', 0)))
    total_boll_wt_open  = float(last.get('groimp_total_boll_wt_open', 0))
    current_gb_brown_wt = float(last.get('groimp_current_green_brown_wt', 0))

    print(f'\n  === 铃重统计 (GroIMP 3D 反馈, 全铃口径) ===')
    if boll_wts:
        print(f'  全铃均重: {boll_wts[-1]:.2f} g/铃'
              f'  (目标 {BOLL_WT_TARGET[0]:.2f}–{BOLL_WT_TARGET[1]:.2f} g)')
    print(f'  全铃总重 (含绿+褐+开): {total_boll_wt_all:.1f} g/株')
    if total_boll_wt_open > 0:
        print(f'    其中: 已开铃: {total_boll_wt_open:.1f} g/株'
              f'   绿/褐铃: {current_gb_brown_wt:.1f} g/株')

    # ── E_pre. Nudge 消融诊断 (方案A v2) ────────────────────────────────
    print(f'\n  === Nudge 消融实验状态 (方案A v2: visible LAI 同步比对) ===')
    _nudge_mode_str = ('NudgeON (v2: visible口径, step后, DAS>=45门槛)'
                       if NUDGE_ENABLED else
                       'NudgeOFF ★ 消融 — LAI一致性软约束已关闭')
    print(f'  运行模式:  {_nudge_mode_str}')
    _n_nudge_days   = sum(1 for e in coupling_log if e.get('nudge_fired', False))
    _total_nudge    = sum(_safe_float(e.get('nudge_wlf_delta', 0)) for e in coupling_log)
    if not NUDGE_ENABLED:
        print(f'  Nudge 触发次数: 0 次 (已消融, 碳平衡与3D结构自然分歧)')
        print(f'  论文表述: "LAI一致性软约束(Nudge)已关闭; Cotton2K碳平衡')
        print(f'            与GroIMP三维结构独立演化, LAI偏差直接报告。"')
    else:
        _pos = sum(1 for e in coupling_log
                   if e.get('nudge_fired') and e.get('nudge_wlf_delta', 0) > 0)
        _neg = sum(1 for e in coupling_log
                   if e.get('nudge_fired') and e.get('nudge_wlf_delta', 0) < 0)
        print(f'  Nudge 有效窗口:     DAS>=45 (groimp_day>={NUDGE_GROIMP_START_DAY}) '
              f'且 GroIMP LAI>={NUDGE_GROIMP_MIN_LAI}')
        print(f'  Nudge 触发次数:     {_n_nudge_days} 天'
              f'  (上推 {_pos} 天 / 下压 {_neg} 天)')
        print(f'  累积 W_lf 修正量:   {_total_nudge:+.2f} g/m²')
        print(f'  参数: THRESHOLD={LAI_NUDGE_THRESHOLD:.0%}'
              f'  tau={LAI_NUDGE_TAU:.1f}天'
              f'  MAX_RATE={LAI_NUDGE_MAX_RATE:.0%}')
    print(f'  配套RGG: CottonModel_2025_{"Coupled" if NUDGE_ENABLED else "NudgeOFF"}.rgg')

    # ── E. 通信延迟 ──────────────────────────────────────
    lat = np.array([e['elapsed_ms'] for e in coupling_log])
    print(f'\n  === 通信效率 (n={len(lat)} 天) ===')
    print(f'  均值/中位/标准差: {lat.mean():.0f} / {np.median(lat):.0f} / {lat.std():.0f} ms')
    print(f'  P95/P99/最大值:   {np.percentile(lat,95):.0f} / '
          f'{np.percentile(lat,99):.0f} / {lat.max():.0f} ms')
    print(f'  总耗时 Σ:         {lat.sum()/1000.0:.1f} s')

    # ── F. 光截获物理性 ─────────────────────────────────
    if valid:
        lis = [e['groimp_li'] for e in valid if e['groimp_li'] > 0]
        if lis:
            flowering_li = [e['groimp_li'] for e in valid
                            if PHASES['flowering'][0] <= e['groimp_day'] <= PHASES['flowering'][1]
                            and e['groimp_li'] > 0]
            li_mean = float(np.mean(lis))
            fli_mean = float(np.mean(flowering_li)) if flowering_li else 0.0
            li_max = float(np.max(lis))
            print(f'\n  === 光截获率 ===')
            print(f'  全期均值/花铃期均值/最大: {li_mean*100:.1f}% / {fli_mean*100:.1f}% / {li_max*100:.1f}%')
            if li_max > 0.97:
                print(f'  ⚠ 最大 LI > 97%, 超物理合理上限')
            elif fli_mean > 0.95:
                print(f'  ⚠ 花铃期 LI 偏高, 建议检查 FluxLM 上限')
            else:
                print(f'  ✓ 物理合理区间')

    # ── G. 产量三路核算 ────────────────────────────────────────
    print(f'\n  === 产量三路交叉核算 (全铃收获日 DAS 200 口径) ===')
    print(f'  制式: 一膜六行 | 膜宽{FILM_WIDTH_CM:.0f}cm | 行距{WIDE_ROW_CM:.0f}+{NARROW_ROW_CM:.0f}cm | 株距{PLANT_SPACING_CM:.0f}cm')
    print(f'  密度: {PLANTS_PER_HA:,} 株/ha | 衣分: {LINT_RATIO_DEFAULT*100:.1f}% (塔河2号审定)')

    # 路径 1: Cotton2K 内部 W_lint (真实碳平衡输出)
    c2k_lint_kg_ha = float(getattr(model, 'W_lint', 0.0))
    c2k_seed_kg_ha = float(getattr(model, 'W_seed', 0.0))
    if c2k_lint_kg_ha <= 0 and 'lint_yield' in df.columns:
        c2k_lint_kg_ha = float(df['lint_yield'].max())

    # 路径 2: GroIMP 全铃重量 × 衣分 × 密度 (全铃口径)
    gro_lint_kg_ha = (total_boll_wt_all * LINT_RATIO_DEFAULT * PLANTS_PER_HA) / 1000.0
    gro_lint_kg_ha_open = (total_boll_wt_open * LINT_RATIO_DEFAULT * PLANTS_PER_HA) / 1000.0

    # 路径 3: 全铃数 × 全铃均重 × 衣分 × 密度 (收获日 DAS 200 理论上限)
    if boll_wts and n_eventual > 0:
        recon_lint_kg_ha = (n_eventual * boll_wts[-1] * LINT_RATIO_DEFAULT *
                            PLANTS_PER_HA) / 1000.0
    else:
        recon_lint_kg_ha = 0.0
    if boll_wts_open_log and n_open > 0:
        recon_lint_kg_ha_open = (n_open * boll_wts_open_log[-1] * LINT_RATIO_DEFAULT *
                                  PLANTS_PER_HA) / 1000.0
    else:
        recon_lint_kg_ha_open = 0.0

    print(f'  ──────────────────────────────────────────────────────────────────')
    print(f'  路径 1 (Cotton2K W_lint, 真实碳平衡):          '
          f'{c2k_lint_kg_ha:>8,.0f} kg/ha')
    print(f'  路径 2 (GroIMP 全铃重量 × 衣分 × 密度):       '
          f'{gro_lint_kg_ha:>8,.0f} kg/ha')
    print(f'  路径 3 (全铃数 × 全铃均重 × 衣分 × 密度):      '
          f'{recon_lint_kg_ha:>8,.0f} kg/ha  ← 收获日 DAS 200 理论上限')
    print(f'  ──────────────────────────────────────────────────────────────────')
    print(f'  路径3 = {n_eventual} 铃 × {boll_wts[-1] if boll_wts else 0:.2f}g/铃'
          f' × {LINT_RATIO_DEFAULT:.3f} × {PLANTS_PER_HA:,} / 1000 = {recon_lint_kg_ha:.0f} kg/ha')

    if c2k_lint_kg_ha > 0 and recon_lint_kg_ha > 0:
        diff_pct = (recon_lint_kg_ha - c2k_lint_kg_ha) / c2k_lint_kg_ha * 100.0
        harvest_eff = c2k_lint_kg_ha / max(recon_lint_kg_ha, 1) * 100
        print(f'  路径1 vs 路径3 偏差: {diff_pct:+.1f}%  '
              f'(收获效率 = 路径1/路径3 = {harvest_eff:.1f}%)')
        if 60.0 <= harvest_eff <= 80.0:
            print(f'  ✓ 收获效率落入文献区间 [60%, 80%] (Bange & Milroy 2004)')
        elif 50.0 <= harvest_eff < 60.0:
            print(f'  ✓ 收获效率 ~50-60% (取消脱落后偏低端常态, 路径3 为理论上限)')
        elif harvest_eff < 50.0:
            print(f'  ⚠ 收获效率 < 50%: 检查路径1 碳平衡或路径3 铃数/均重是否过估')
        else:
            print(f'  ⚠ 收获效率 > 80%: 检查 W_seed 是否达 cap 或路径3 偏低')
    print(f'\n  === 与田间实测对照 (实测皮棉 = {MEASURED_LINT_KGHA:.0f} kg/ha) ===')
    if c2k_lint_kg_ha > 0:
        err = (c2k_lint_kg_ha - MEASURED_LINT_KGHA) / MEASURED_LINT_KGHA * 100
        rating = ('Excellent' if abs(err) < 10 else
                  'Good'      if abs(err) < 20 else
                  'Fair'      if abs(err) < 30 else 'Poor')
        print(f'  路径 1 误差:  {err:+.1f}%  [{rating}, Jamieson 1991]')
    if recon_lint_kg_ha > 0:
        err = (recon_lint_kg_ha - MEASURED_LINT_KGHA) / MEASURED_LINT_KGHA * 100
        print(f'  路径 3 误差:  {err:+.1f}%  (理论最大值, 通常高估)')

    # ── I. 耦合源使用统计 ────────────────────────────────
    print(f'\n  === 耦合源诊断 ===')
    n_attempts = _coupling_stats['attempts']
    n_rgg3d    = _coupling_stats['used_RGG_3D_direct']
    n_pyfb     = _coupling_stats['used_python_rebuild']
    n_invalid  = _coupling_stats['fallback_invalid']
    if n_attempts > 0:
        print(f'  尝试次数:                {n_attempts}')
        print(f'  RGG 3D 直传:        {n_rgg3d}  ({n_rgg3d/n_attempts*100:.1f}%)')
        print(f'  Python RUE 回退:         {n_pyfb}  ({n_pyfb/n_attempts*100:.1f}%)')
        print(f'  GroIMP 反馈无效 → 独立:   {n_invalid}  ({n_invalid/n_attempts*100:.1f}%)')
        if n_rgg3d == 0:
            print(f'  ⚠ 全部走 Python 回退. 检查 RGG 端 socketSendFeedback 是否含 '
                  f'groimp_net_photo_gCH2O 字段')
        elif n_pyfb > n_attempts * 0.2:
            print(f'  ⚠ 超过 20% 走 Python 回退')

    diag = model.get_coupling_diagnostics()
    print(f'  Cotton2KDaily 内部步进诊断: {diag}')

    # ── J. 冷启动碳桥接诊断 ──────────────────────────────
    if hasattr(model, 'get_bridge_diagnostics'):
        bdiag = model.get_bridge_diagnostics()
        print(f'\n  === 冷启动碳桥接诊断 ===')
        print(f'  桥接窗口:        {bdiag.get("bridge_window", "DAS 0-35")}')
        print(f'  桥接激活天数:     {bdiag.get("bridge_active_days", 0)}')
        print(f'  累计桥接补碳:     {bdiag.get("bridge_carbon_added_g_per_m2", 0):.1f} g DM/m²')
        if bdiag.get("bridge_active_days", 0) >= 30:
            print(f'  ✓ 桥接机制正常激活, 苗期碳缺口已补齐')
        else:
            print(f'  ⚠ 桥接激活天数 < 30, 可能 GroIMP 苗期已提供充足光合')

    # ── J2. 光合 GROSS FLOOR 诊断 ────────────────────────
    if hasattr(model, 'get_floor_diagnostics'):
        fdiag = model.get_floor_diagnostics()
        print(f'\n  === GROSS FLOOR 诊断 ===')
        print(f'  底板窗口:        {fdiag.get("floor_window", "DAS 0-65")}')
        print(f'  底板比例:        {fdiag.get("floor_fraction", 0.65):.2f} × Beer-Lambert')
        print(f'  GROSS 校正天数:   {fdiag.get("floor_active_days", 0)}')
        print(f'  累计 GROSS 校正:  {fdiag.get("floor_carbon_added_g_per_m2", 0):.1f} g DM/m²')

    # ── J3. NET FLOOR 诊断 ────────────────────────────────
    if hasattr(model, 'get_net_floor_diagnostics'):
        ndiag = model.get_net_floor_diagnostics()
        n_active = ndiag.get("net_floor_active_days", 0)
        c_added = ndiag.get("net_floor_carbon_added_g_per_m2", 0)
        print(f'\n  === NET FLOOR 诊断 ===')
        print(f'  NET FLOOR 激活天数: {n_active}')
        print(f'  NET FLOOR 累计补碳: {c_added:.1f} g DM/m²')
        if n_active >= 10:
            print(f'  ✓ NET FLOOR 已生效')
        elif n_active >= 1:
            print(f'  ✓ NET FLOOR 偶发激活')
        else:
            print(f'  ⚠ NET FLOOR 未激活: 检查 inj_net 是否正常')

    # ── J4. 打顶事件诊断 ──────────────────────────────────
    if hasattr(model, 'get_topping_diagnostics'):
        tdiag = model.get_topping_diagnostics()
        print(f'\n  === 打顶事件诊断 ===')
        print(f'  打顶日期:        {tdiag.get("topping_date", "2025-07-16")}')
        print(f'  (播后85天, 出苗后75天)')
        print(f'  打顶已应用:      {"✓" if tdiag.get("topping_applied") else "✗ 未应用 (检查日期是否在模拟期内)"}')
        print(f'  W_lf 离散损失:    {tdiag.get("topping_loss_lf_g_per_m2", 0):.1f} g/m² '
              f'(目标 {tdiag.get("topping_leaf_loss_frac_target", 0.18)*100:.0f}% 顶端切除)')

    # ── K. 全期 LAI vs 20 实测锚点 (论文级 KPI) ─────────
    if coupling_log:
        c2k_peak_lai = max([e.get('c2k_lai', 0) for e in coupling_log], default=0)
        gro_peak_lai = max([e.get('groimp_lai', 0) for e in coupling_log], default=0)
        obs_peak = MEASURED_LAI_PEAK
        print(f'\n  === 全期 LAI vs 20 实测锚点 (论文级核心 KPI) ===')
        print(f'  田间实测峰 LAI:    {obs_peak:.3f} ± {MEASURED_LAI_PEAK_SD:.3f}'
              f'  (DAS={MEASURED_LAI_PEAK_DAS}, 4区均值)')
        c2k_err = (c2k_peak_lai - obs_peak) / obs_peak * 100 if obs_peak > 0 else 0
        gro_err = (gro_peak_lai - obs_peak) / obs_peak * 100 if obs_peak > 0 else 0
        print(f'  Cotton2K 模拟峰:    {c2k_peak_lai:.3f}  误差 {c2k_err:+.1f}%  '
              f'{"✓" if abs(c2k_err)<15 else "⚠" if abs(c2k_err)<30 else "✗"}')
        print(f'  GroIMP 模拟峰:     {gro_peak_lai:.3f}  误差 {gro_err:+.1f}%  '
              f'{"✓" if abs(gro_err)<15 else "⚠" if abs(gro_err)<30 else "✗"}')

        print(f'\n  --- 全期 20 实测锚点对照 ---')
        print(f'  {"DAS":>4} {"Obs":>7} {"C2K":>7} {"GRO":>7} {"C2K err":>9} {"GRO err":>9}')

        if hasattr(model, 'get_all_records'):
            df_rec = model.get_all_records()
            df_rec['date'] = pd.to_datetime(df_rec['date'])
            try:
                from cotton2k_daily import EMERGE_DATE as _ED
                emerge_dt = pd.Timestamp(_ED)
            except Exception:
                emerge_dt = df_rec['date'].iloc[0]
            df_rec['DAS'] = (df_rec['date'] - emerge_dt).dt.days

            obs_x = np.array([d for d,_ in MEASURED_LAI_SERIES])
            obs_y = np.array([v for _,v in MEASURED_LAI_SERIES])
            sim_c2k_y = []
            sim_gro_y = []
            for das, obs in MEASURED_LAI_SERIES:
                r = df_rec[df_rec['DAS']==das]
                cv = float(r['leaf_area_index'].iloc[0]) if not r.empty else np.nan
                gv = np.nan
                for e in coupling_log:
                    if e.get('groimp_day') == das + 1:
                        gv = e.get('groimp_lai', np.nan)
                        break
                sim_c2k_y.append(cv)
                sim_gro_y.append(gv)
                c2k_e = (cv-obs)/obs*100 if not (np.isnan(cv) or obs==0) else np.nan
                gro_e = (gv-obs)/obs*100 if not (np.isnan(gv) or obs==0) else np.nan
                print(f'  {das:>4} {obs:>7.3f} {cv:>7.3f} {gv:>7.3f} '
                      f'{c2k_e:>+8.1f}% {gro_e:>+8.1f}%')

            sim_c2k_y = np.array(sim_c2k_y, dtype=float)
            sim_gro_y = np.array(sim_gro_y, dtype=float)
            obs_mean = obs_y.mean()
            ss_res_c2k = np.nansum((sim_c2k_y - obs_y) ** 2)
            ss_res_gro = np.nansum((sim_gro_y - obs_y) ** 2)
            ss_tot = np.sum((obs_y - obs_mean) ** 2)
            nse_c2k = 1 - ss_res_c2k / ss_tot if ss_tot > 0 else np.nan
            nse_gro = 1 - ss_res_gro / ss_tot if ss_tot > 0 else np.nan
            rmse_c2k = np.sqrt(np.nanmean((sim_c2k_y - obs_y) ** 2))
            rmse_gro = np.sqrt(np.nanmean((sim_gro_y - obs_y) ** 2))
            mbe_c2k = np.nanmean(sim_c2k_y - obs_y)
            mbe_gro = np.nanmean(sim_gro_y - obs_y)
            nrmse_c2k = rmse_c2k / obs_mean * 100
            nrmse_gro = rmse_gro / obs_mean * 100
            print(f'\n  --- 模型 vs 实测 NSE/RMSE/MBE (n=20) ---')
            print(f'  {"":>8}{"NSE":>10}{"RMSE":>10}{"MBE":>10}{"nRMSE%":>10}')
            grade_c = "Good" if nse_c2k>=0.65 else "Fair" if nse_c2k>=0.50 else "Poor"
            grade_g = "Good" if nse_gro>=0.65 else "Fair" if nse_gro>=0.50 else "Poor"
            print(f'  Cotton2K:{nse_c2k:>10.4f}{rmse_c2k:>10.3f}{mbe_c2k:>+10.3f}{nrmse_c2k:>10.1f}  [{grade_c}]')
            print(f'  GroIMP:  {nse_gro:>10.4f}{rmse_gro:>10.3f}{mbe_gro:>+10.3f}{nrmse_gro:>10.1f}  [{grade_g}]')

    print(f'\n{"="*70}\n  cotton2k与GroIMP耦合 仿真完成\n{"="*70}')

# 9. 入口
# =======================================================================
def _make_default_output_dir() -> str:
    """在脚本所在目录下建 runs/时间戳子目录, 防止覆盖."""
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(SCRIPT_DIR, 'runs', datetime.now().strftime('%Y%m%d_%H%M%S'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'cotton2k与GroIMP耦合 — 逐日耦合服务器 [2025独立验证年, 参数冻结自2022标定年] '
                    f'(配合 CottonModel_2025_Coupled.rgg / CottonModel_2025_NudgeOFF.rgg)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--weather', default='weather2025.csv', help='气象 CSV')
    parser.add_argument('--host',    default=DEFAULT_HOST)
    parser.add_argument('--port',    type=int, default=DEFAULT_PORT)
    parser.add_argument('--output',  default=None,
                        help='输出目录 (默认 runs/YYYYMMDD_HHMMSS/)')
    parser.add_argument('--timeout', type=float, default=60.0,
                        help='GroIMP 握手超时秒数 (默认 60s)')
    parser.add_argument('--nudge-off', action='store_true', default=False,
                        help='消融实验: 关闭LAI一致性软约束(Nudge); '
                             '配合 CottonModel_2025_NudgeOFF.rgg 使用')
    args = parser.parse_args()

    output_dir = args.output or _make_default_output_dir()

    try:
        nudge_on = not args.nudge_off
        if not nudge_on:
            print('[消融实验] --nudge-off 激活: 配合 CottonModel_2025_NudgeOFF.rgg')
        run_daily_coupling(args.weather, args.host, args.port,
                           output_dir, args.timeout, nudge_enabled=nudge_on)
    except KeyboardInterrupt:
        print('\n[中断] 用户 Ctrl+C 退出')
        sys.exit(130)
    except Exception as e:
        print(f'\n[致命错误] {type(e).__name__}: {e}')
        raise
