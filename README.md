# 项目说明

这是一个基于遗传算法和贪心覆盖策略的电动货车路径规划系统，支持客户服务时间窗、充电站充电时间窗、车辆容量约束等多种约束条件，并提供数据预处理、聚类选站、路径生成及评估功能。

---

## 目录结构

```
├── algorithm/                # 核心算法包
│   ├── core.py               # 核心遗传算法实现（run_gavrptw, eval_vrptw, 带站点版本等）
│   ├── greedy_cover.py       # 聚类与贪心覆盖算法（cluster_customers, assign_customers_to_centers 等）
│   └── utils.py              # 实例加载、文本转 JSON 等工具函数
│
├── data/
│   └── json_customize/       # 定制化测试与生产数据 JSON
│
├── tests/                    # 各种单元测试与集成测试脚本
│   ├── test_kmeans.py        # 训练簇中心并验证簇分配完整性
│   ├── test_compute_distance.py  # 验证站点到兴趣点距离计算正确性
│   ├── test_customize.py     # 基本 GA 路径生成不抛错测试
│   ├── test_timewindow.py    # 带时间窗与带站点版本的运行测试
│   ├── test_logisic.py       # 极简示例测试: 检查是否经过充电站
│   └── test_ind_station_setect.py  # 最近充电站选择逻辑测试
│
├── text2json_customize.py    # 文本格式数据转 JSON 脚本入口
├── assign_customer.py        # 综合示例：聚类 + 覆盖 + 站点选择演示脚本
├── main.py                   # PyCharm 示例入口脚本
├── requirements.txt          # Python 依赖列表
└── README.md                 # 项目说明文档（本文件）
```

---

## 安装与依赖

1. 克隆仓库：

   ```bash
   git clone <仓库地址>.git
   cd <项目根目录>
   ```
2. 创建并激活虚拟环境（可选）：

   ```bash
   python3 -m venv venv
   source venv/bin/activate   # Linux/macOS
   venv\\Scripts\\activate  # Windows
   ```
3. 安装依赖：

   ```bash
   pip install -r requirements.txt
   ```

---

## 快速开始

1. **数据预处理**：将原始文本数据转换为 JSON 格式

   ```bash
   python text2json_customize.py
   ```

   生成的 JSON 文件会保存到 `data/json_customize/` 目录。

2. **聚类与站点覆盖**：示例脚本

   ```bash
   python assign_customer.py
   ```

   输出聚类中心、客户覆盖映射以及贪心选出的充电站列表。

3. **运行路径规划（带充电站）**：

   ```bash
   python -c "from algorithm.core import run_gavrptw_with_stations; \
   run_gavrptw_with_stations(instance_name='c105', unit_cost=8.0, init_cost=100.0, \
   wait_cost=50, delay_cost=50, ind_size=18, pop_size=400, cx_pb=0.85, mut_pb=0.02, \
   n_gen=100, energy_per_km=1, export_csv=False, customize_data=True)"
   ```

   结果会打印最佳路径及相关成本。

---

## 测试

项目使用 `pytest` 进行单元与集成测试，覆盖聚类、距离计算、时间窗、充电站逻辑等。运行所有测试：

```bash
pytest -q
```

---

## 开发与贡献

* **代码规范**：PEP8，函数注释齐全。
* **新增功能**：请先在 `tests/` 下添加对应的测试用例，确保逻辑无误后再提交 PR。
* **文档更新**：如功能或接口有改动，请同步更新 `README.md` 与代码内文档字符串。

欢迎 issue 和 PR，共同完善此电动货车路径优化项目！
