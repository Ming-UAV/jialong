import random
from algorithm.utils import load_instance
from algorithm.core import run_gavrptw_with_stations
#测试是否会选择最近的充电站
def test_run_gavrptw_executes():
    random.seed(64)
    instance_name = 'extrme_simple_test'

    # 先加载实例，自动计算客户数
    from algorithm.utils import load_instance
    inst = load_instance(f"data/json_customize/{instance_name}.json")
    cust_keys = [k for k in inst if k.startswith('customer_')]
    ind_size = len(cust_keys)

    run_gavrptw_with_stations(
        instance_name=instance_name,
        unit_cost=8.0,
        init_cost=100.0,
        wait_cost=0,
        delay_cost=0,
        ind_size=ind_size,   # <-- 动态取
        pop_size=400,
        cx_pb=0.85,
        mut_pb=0.02,
        n_gen=50,
        energy_per_km=1,
        export_csv=False,
        customize_data=True
    )
