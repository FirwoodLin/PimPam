import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from datetime import datetime
from collections import defaultdict

def parse_tasklet_cycles(file_path):
    pattern = re.compile(r'DPU:\s*(\d+),\s*tasklet:\s*(\d+),\s*cycle:\s*(\d+),\s*root_num:\s*(\d+)')
    dpu_tasklet_data = defaultdict(lambda: defaultdict(int))

    with open(file_path, 'r') as f:
        for line in f:
            match = pattern.search(line)
            if match:
                dpu, tasklet, cycle, _ = map(int, match.groups())
                dpu_tasklet_data[dpu][tasklet] = cycle

    return dpu_tasklet_data

def get_unique_filename(base_path):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_name, ext = os.path.splitext(base_path)
    counter = 1
    while True:
        new_path = f"{base_name}_{timestamp}" + (f"_{counter}" if counter > 1 else "") + ext
        if not os.path.exists(new_path):
            return new_path
        counter += 1

def plot_min_range_stacked(dpu_tasklet_data, output_path):
    dpu_ids = sorted(dpu_tasklet_data.keys())
    x = np.arange(len(dpu_ids))

    min_cycles = []
    ranges = []
    max_cycles = []

    for dpu in dpu_ids:
        tasklet_cycles = [dpu_tasklet_data[dpu].get(tid, 0) for tid in range(16)]
        min_val = min(tasklet_cycles)
        max_val = max(tasklet_cycles)
        min_cycles.append(min_val)
        ranges.append(max_val - min_val)
        max_cycles.append(max_val)

    max_runtime_cycle = max(max_cycles)
    avg_min = np.mean(min_cycles)
    avg_range = np.mean(ranges)
    runtime_ms = max_runtime_cycle * 2.85e-6

    # 绘图
    plt.figure(figsize=(14, 6))
    ax = plt.gca()

    # 堆叠柱状图：下层为 min，顶层为 range
    plt.bar(x, min_cycles, color='skyblue', label='Min Cycle')
    plt.bar(x, ranges, bottom=min_cycles, color='salmon', label='Range (Max - Min)')

    # 样式设置
    plt.title("DPU Thread Cycle: Min and Range", fontsize=14)
    plt.xlabel("DPU ID", fontsize=12)
    plt.ylabel("Cycle Count", fontsize=12)

    # 控制 x 轴标签显示数量
    n_labels = 7
    step = max(1, len(dpu_ids) // n_labels)
    xticks_labels = [str(dpu_ids[i]) if i % step == 0 else "" for i in range(len(dpu_ids))]
    plt.xticks(x, xticks_labels, rotation=45, fontsize=8)

    # 去掉底部边框线
    ax.spines['bottom'].set_visible(False)

    # 图例和注释
    plt.legend(fontsize=10)
    text = f"Max Cycle ≈ {max_runtime_cycle} → Runtime ≈ {runtime_ms:.2f} ms\n" \
           f"Avg Min = {avg_min:.0f}, Avg Range = {avg_range:.0f}"
    plt.text(0.02, 0.95, text, transform=ax.transAxes, fontsize=10,
             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8))

    plt.tight_layout()
    out_path = get_unique_filename(output_path)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"图像已保存至: {out_path}")

# === 主程序入口 ===
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python plot_min_range.py input.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    if not os.path.isfile(input_file):
        print(f"Error: file '{input_file}' not found.")
        sys.exit(1)

    dpu_tasklet_data = parse_tasklet_cycles(input_file)
    output_base = os.path.splitext(input_file)[0] + "_min_range.png"
    plot_min_range_stacked(dpu_tasklet_data, output_base)
