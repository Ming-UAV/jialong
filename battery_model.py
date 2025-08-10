import numpy as np
import matplotlib.pyplot as plt

# 参数（Type A）
beta0 = 0.0107
beta1 = 0.4163
b = 0.3959

# 定义公式
def NDC(C_rate, N_c):
    return 1 - beta0 * np.exp(beta1 * C_rate) * (N_c ** b)

# 数据
N_vals = np.linspace(0, 1000, 200)  # 循环次数
C_rates = [0.2, 0.5, 1.0, 2.0]     # 不同C_rate
colors = ['b', 'k', 'r', 'g']
linestyles = ['dotted', 'solid', 'dashed', 'dashdot']

# 绘图
plt.figure(figsize=(8,6))
for C, color, ls in zip(C_rates, colors, linestyles):
    ndc_vals = NDC(C, N_vals)
    plt.plot(N_vals, ndc_vals, color=color, linestyle=ls, label=f'{C}C Data')

    # 找到第一个小于等于 0.8 的循环次数
    below_80 = np.where(ndc_vals <= 0.8)[0]
    if below_80.size > 0:
        N_at_80 = N_vals[below_80[0]]
        print(f"{C}C: 容量降到 80% 时的循环次数 ≈ {N_at_80:.1f}")
    else:
        print(f"{C}C: 在 250 循环内未降到 80%")

plt.xlabel('Cycles')
plt.ylabel('Normalized Discharge Capacity')
plt.title('Capacity Fade at Different C-rates')
plt.legend()
plt.grid(True)
import os

out_dir = r"D:\OneDrive\Desktop\dissertation-GA-method\results"
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "battery_fade_model.png")
plt.savefig(out_path, dpi=300, bbox_inches='tight')
plt.show()


print(f"图像已保存到: {out_path}")