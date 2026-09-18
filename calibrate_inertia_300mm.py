#!/usr/bin/env python3
"""
calibrate_inertia_300mm.py
──────────────────────────
Cotton2K 300mm 根区 INERTIA_COEF 自动重标定脚本

用途
────
土壤桶深从 100mm 统一为 300mm 后，Ks 均值升高（虚假胁迫减少），
净光合偏高，产量基线上移，需重新找到使耦合产量 ≈ 3162.6 kg/ha 的
INERTIA_COEF 值并写入 cotton2k_daily.py。

使用方法
────────
# 第一步：运行完整 2025 耦合实验
#   python daily_coupling_server.py --weather weather.csv

# 第二步：运行本脚本（指向上一步输出的 CSV）
#   python calibrate_inertia_300mm.py --csv runs/YYYYMMDD_HHMMSS/output_daily_coupled.csv

# 若 --csv 省略，脚本自动搜索最新 runs/ 子目录
#   python calibrate_inertia_300mm.py

# --dry-run：只显示建议，不写文件
#   python calibrate_inertia_300mm.py --csv XXX.csv --dry-run

输出
────
• 控制台：当前产量 / 目标 / 误差 / 建议新系数 / 预计轮次
• 自动 patch cotton2k_daily.py（覆盖 INERTIA_COEF = ... 行）
• 输出标定历史文件 calibration_history.csv（追加写）

标定算法
────────
线性近似（产量与 INERTIA_COEF 在小范围内近似线性）：
    new_coef = old_coef × (target_lint / current_lint)

收敛判据：|误差| < 2%（约等于 Excellent 级精度）
通常 1–3 轮耦合实验即可收敛。
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from datetime import datetime
from pathlib import Path

# ─────────────────────── 可配置常量 ────────────────────────────────────
TARGET_LINT_KGHA    = 3162.6    # 2025 实测皮棉产量 kg/ha
TARGET_ERR_PCT      = 2.0       # 收敛判据：误差绝对值 < 2% 视为达标
COEF_MIN            = 0.00030   # 物理下限（防止跑飞）
COEF_MAX            = 0.00120   # 物理上限
COTTON2K_DAILY      = 'cotton2k_daily.py'  # 相对或绝对路径
HISTORY_CSV         = 'calibration_history.csv'


# ─────────────────────── 辅助函数 ──────────────────────────────────────

def find_latest_csv() -> Path | None:
    """在 runs/ 子目录中搜索最新的 output_daily_coupled.csv。"""
    runs = Path('runs')
    if not runs.exists():
        return None
    candidates = sorted(runs.glob('*/output_daily_coupled.csv'), reverse=True)
    return candidates[0] if candidates else None


def read_lint_yield(csv_path: Path) -> float:
    """从 output_daily_coupled.csv 读取最终皮棉产量（lint_yield 列最大值）。"""
    best = 0.0
    with open(csv_path, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            v = float(row.get('lint_yield', 0) or 0)
            if v > best:
                best = v
    return best


def read_current_coef(py_path: Path) -> float:
    """从 cotton2k_daily.py 读取当前 INERTIA_COEF 值。"""
    text = py_path.read_text(encoding='utf-8')
    m = re.search(r'^INERTIA_COEF\s*=\s*([0-9.eE+-]+)', text, re.MULTILINE)
    if not m:
        raise ValueError(f'在 {py_path} 中未找到 INERTIA_COEF 定义行')
    return float(m.group(1))


def patch_coef(py_path: Path, new_coef: float, dry_run: bool = False) -> None:
    """将 INERTIA_COEF 更新写入 cotton2k_daily.py。"""
    text = py_path.read_text(encoding='utf-8')
    # 替换模式：INERTIA_COEF = <数字>  # 可选注释
    pattern = r'^(INERTIA_COEF\s*=\s*)[0-9.eE+-]+(.*)'
    new_line = rf'\g<1>{new_coef:.6f}   # ← 300mm 标定 ({datetime.now().strftime("%Y-%m-%d")})\2'
    new_text, n = re.subn(pattern, new_line, text, flags=re.MULTILINE)
    if n == 0:
        raise RuntimeError('patch 失败：未找到 INERTIA_COEF 行，请检查文件格式')
    if n > 1:
        raise RuntimeError(f'patch 失败：找到 {n} 处 INERTIA_COEF 定义，应恰好为 1 处')
    if dry_run:
        print(f'  [dry-run] 将写入: INERTIA_COEF = {new_coef:.6f}')
        return
    py_path.write_text(new_text, encoding='utf-8')
    print(f'  ✓ 已写入 {py_path}: INERTIA_COEF = {new_coef:.6f}')


def append_history(record: dict) -> None:
    """追加写入标定历史 CSV。"""
    fieldnames = ['timestamp', 'old_coef', 'current_lint', 'target_lint',
                  'error_pct', 'new_coef', 'converged', 'csv_source']
    write_header = not Path(HISTORY_CSV).exists()
    with open(HISTORY_CSV, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            w.writeheader()
        w.writerow(record)


def estimate_rounds(err_pct: float, damping: float = 0.8) -> int:
    """粗估收敛轮次（线性近似 + 保守阻尼）。"""
    import math
    if abs(err_pct) <= TARGET_ERR_PCT:
        return 0
    n = math.ceil(math.log(TARGET_ERR_PCT / abs(err_pct)) / math.log(damping))
    return max(1, n)


# ─────────────────────── 主程序 ────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Cotton2K 300mm INERTIA_COEF 自动重标定',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--csv', default=None,
        help='耦合实验输出 CSV 路径 (默认: 自动搜索 runs/ 最新目录)',
    )
    parser.add_argument(
        '--target', type=float, default=TARGET_LINT_KGHA,
        help=f'目标皮棉产量 kg/ha (默认: {TARGET_LINT_KGHA})',
    )
    parser.add_argument(
        '--c2k', default=COTTON2K_DAILY,
        help=f'cotton2k_daily.py 路径 (默认: {COTTON2K_DAILY})',
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='只显示建议，不修改文件',
    )
    parser.add_argument(
        '--force-coef', type=float, default=None,
        help='直接指定新系数（跳过自动估算，强制写入）',
    )
    args = parser.parse_args()

    # ── 1. 定位 CSV ──────────────────────────────────────────────────
    print('=' * 70)
    print('  Cotton2K 300mm INERTIA_COEF 重标定')
    print(f'  目标皮棉产量: {args.target:.1f} kg/ha  |  收敛判据: |误差| < {TARGET_ERR_PCT}%')
    print('=' * 70)

    csv_path = Path(args.csv) if args.csv else find_latest_csv()
    if csv_path is None or not csv_path.exists():
        print('\n[错误] 未找到耦合输出 CSV。请先运行耦合实验，或用 --csv 指定路径。')
        print('  示例:')
        print('    python daily_coupling_server.py --weather weather.csv')
        print('    python calibrate_inertia_300mm.py --csv runs/xxx/output_daily_coupled.csv')
        sys.exit(1)
    print(f'\n  [读取] {csv_path}')

    # ── 2. 读取当前产量 ──────────────────────────────────────────────
    current_lint = read_lint_yield(csv_path)
    if current_lint <= 0:
        print('[错误] 未能从 CSV 读到有效 lint_yield，请检查文件列名（应含 lint_yield 列）。')
        sys.exit(1)

    # ── 3. 读取当前系数 ──────────────────────────────────────────────
    py_path = Path(args.c2k)
    if not py_path.exists():
        print(f'[错误] 未找到 {py_path}，请用 --c2k 指定正确路径。')
        sys.exit(1)

    old_coef = read_current_coef(py_path)
    err_pct  = (current_lint - args.target) / args.target * 100.0

    print(f'\n  当前 INERTIA_COEF : {old_coef:.6f}')
    print(f'  当前模拟产量      : {current_lint:,.1f} kg/ha')
    print(f'  目标产量          : {args.target:,.1f} kg/ha')
    print(f'  误差              : {err_pct:+.2f}%')

    # ── 4. 判断是否收敛 ──────────────────────────────────────────────
    if abs(err_pct) <= TARGET_ERR_PCT and args.force_coef is None:
        print(f'\n  ✓ 误差 {err_pct:+.2f}% 已在 ±{TARGET_ERR_PCT}% 以内，无需重标定！')
        print(f'  INERTIA_COEF = {old_coef:.6f} 即为 300mm 最终标定值。')
        append_history(dict(
            timestamp=datetime.now().isoformat(timespec='seconds'),
            old_coef=old_coef, current_lint=round(current_lint, 1),
            target_lint=args.target, error_pct=round(err_pct, 3),
            new_coef=old_coef, converged=True,
            csv_source=str(csv_path),
        ))
        print('=' * 70)
        return

    # ── 5. 计算建议新系数 ────────────────────────────────────────────
    if args.force_coef is not None:
        new_coef = args.force_coef
        print(f'\n  [强制写入] INERTIA_COEF = {new_coef:.6f}')
    else:
        # 线性近似：产量 ∝ INERTIA_COEF（在小范围内成立）
        raw = old_coef * (args.target / current_lint)
        # 保守阻尼（防止过冲）
        damping = 0.80
        new_coef = old_coef + (raw - old_coef) * damping
        new_coef = float(max(COEF_MIN, min(COEF_MAX, new_coef)))

        est_new_lint = current_lint * (new_coef / old_coef)
        est_new_err  = (est_new_lint - args.target) / args.target * 100.0
        rounds = estimate_rounds(err_pct)

        print(f'\n  [线性估算]')
        print(f'   原始估算        : {raw:.6f}')
        print(f'   保守阻尼 (×{damping:.2f}) : {new_coef:.6f}')
        print(f'   预计下轮产量    : {est_new_lint:,.0f} kg/ha  (误差 {est_new_err:+.1f}%)')
        print(f'   预计收敛轮次    : {rounds} 轮')
        print(f'\n  建议写入: INERTIA_COEF = {new_coef:.6f}')

    # ── 6. 安全检查 ──────────────────────────────────────────────────
    delta_pct = abs(new_coef - old_coef) / old_coef * 100
    if delta_pct > 30:
        print(f'\n  ⚠ 警告: 系数变化 {delta_pct:.1f}% > 30%，可能不稳定。')
        print(f'    建议检查耦合结果是否正常（Ks 均值 > 0.75，LAI 峰值 3.5–4.8）')
        if not args.dry_run:
            ans = input('    是否仍要写入? [y/N] ').strip().lower()
            if ans != 'y':
                print('  已取消。')
                sys.exit(0)

    # ── 7. 写入文件 ──────────────────────────────────────────────────
    patch_coef(py_path, new_coef, dry_run=args.dry_run)

    # ── 8. 记录历史 ──────────────────────────────────────────────────
    append_history(dict(
        timestamp=datetime.now().isoformat(timespec='seconds'),
        old_coef=old_coef, current_lint=round(current_lint, 1),
        target_lint=args.target, error_pct=round(err_pct, 3),
        new_coef=round(new_coef, 7), converged=False,
        csv_source=str(csv_path),
    ))
    print(f'  📋 标定历史已追加: {HISTORY_CSV}')

    # ── 9. 下一步引导 ────────────────────────────────────────────────
    print('\n' + '=' * 70)
    if args.dry_run:
        print('  [dry-run 完成] 未写入任何文件。去掉 --dry-run 后重新运行以应用更改。')
    else:
        print('  下一步:')
        print('  1. 重新启动 GroIMP (加载 CottonModel_2025_Coupled.rgg)')
        print('  2. python daily_coupling_server.py --weather weather.csv')
        print('  3. python calibrate_inertia_300mm.py --csv runs/新目录/output_daily_coupled.csv')
        print(f'  目标: |误差| < {TARGET_ERR_PCT}%  (当前 {err_pct:+.2f}%)')
        print()
        print('  收敛后同步更新 NudgeOFF 耦合服务器并重跑 --nudge-off 实验。')
    print('=' * 70)


if __name__ == '__main__':
    main()
