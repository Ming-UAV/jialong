import os
import json
import numpy as np
import pytest
from algorithm.greedy_cover import cluster_customers

@pytest.mark.parametrize("K", [1, 3, 5])
def test_cluster_customers_smoke(K):
    # 1. 项目根目录就是这个测试文件的目录
    project_root = os.path.dirname(__file__)
    json_file = os.path.join(
        project_root, "data", "json_customize", "Customized_Data.json"
    )
    assert os.path.exists(json_file), f"文件不存在：{json_file}"

    # 2. 调用聚类
    clusters, centers = cluster_customers(
        instance_name="Customized_Data",
        n_clusters=K,
        customize_data=True,
        random_state=42
    )

    # 3. 断言中心格式
    assert isinstance(centers, np.ndarray)
    assert centers.ndim == 2 and centers.shape[1] == 2

    # 4. 断言所有 customer 都被赋到某个簇
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    all_ids = {k for k in data.keys() if k.startswith("customer_")}
    assigned = {cid for members in clusters.values() for cid in members}
    assert all_ids == assigned

    # 5. 簇编号连续，每个簇至少一个成员
    assert set(clusters.keys()) == set(range(len(clusters)))
    for memb in clusters.values():
        assert len(memb) >= 1


# 只弹图,需要修改调用的json文件名称
clusters, centers = cluster_customers("RC204", 5, customize_data=True, visualize=True)


