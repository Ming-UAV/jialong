
import numpy as np
import json

from sklearn.cluster import KMeans
#随机化仓库位置，
def load_pois(n=10, seed=42):
    np.random.seed(seed)
    return np.random.rand(n, 2) * 100  # 生成 [0, 100] 范围内的二维兴趣点

#计算充电站和兴趣点的距离，一个充电站和每个兴趣点的距离
def compute_station_poi_distances(station: np.ndarray,
                                  pois: np.ndarray
                                 ) -> np.ndarray:
    """
    :param station: 单个充电站坐标，shape=(2,) 或 (1,2)
    :param pois: 兴趣点坐标数组，shape=(n,2)
    :return: numpy 数组，shape=(n,)，第 i 个元素是 station 到 POI[i] 的欧氏距离
    """
    station = np.asarray(station).reshape(1, 2)
    pois    = np.asarray(pois).reshape(-1, 2)
    # pois - station 会广播成 (n,2)，然后按行计算 L2 范数
    distances = np.linalg.norm(pois - station, axis=1)
    return distances

#k-means初步聚类

import os
import json
import numpy as np
from sklearn.cluster import KMeans
from .utils import load_instance
from . import BASE_DIR

import os
import numpy as np
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from .utils import load_instance
from . import BASE_DIR

def cluster_customers(instance_name: str,
                      n_clusters: int = None,
                      customize_data: bool = False,
                      random_state: int = 0,
                      visualize: bool = False,
                      save_plot: bool = False,
                      plot_path: str = None,
                      auto_k: bool = False,
                      k_range: range = range(2, 11)
                      ):
    """
    对 JSON 文件中的客户坐标进行 K-Means 聚类。

    :param instance_name: JSON 文件名（不含 .json）
    :param n_clusters: 指定簇数，如果 auto_k=True 则此参数会被忽略
    :param auto_k: 是否基于 Silhouette Score 自动选簇数
    :param k_range: 候选簇数区间，仅在 auto_k=True 时使用
    :param customize_data: 是否使用 data/json_customize 目录
    :param random_state: 随机种子，保证结果可复现
    :param visualize: 是否弹出绘图窗口
    :param save_plot: 是否保存绘图结果
    :param plot_path: 保存图像路径（含文件名），仅在 save_plot=True 时有效
    :return:
      - clusters: dict, key 为簇编号（0..K-1），value 为该簇包含的 customer_id 列表
      - centers: numpy.ndarray, shape=(K,2)，每行是聚类中心坐标
    """
    # 1. 定位并读取 JSON
    json_data_dir = (os.path.join(BASE_DIR, 'data', 'json_customize')
                     if customize_data
                     else os.path.join(BASE_DIR, 'data', 'json'))
    json_file = os.path.join(json_data_dir, f'{instance_name}.json')
    instance = load_instance(json_file=json_file)
    if instance is None:
        return {}, np.empty((0, 2))

    # 2. 提取 customer 坐标
    coords = []
    customer_ids = []
    for key, val in instance.items():
        if key.startswith('customer_'):
            coords.append([val['coordinates']['x'], val['coordinates']['y']])
            customer_ids.append(key)
    coords = np.array(coords)
    if coords.size == 0:
        return {}, np.empty((0, 2))

    # 3. 自动选最优 K
    if auto_k:
        best_k = None
        best_score = -1.0
        for k in k_range:
            km = KMeans(n_clusters=k, random_state=random_state).fit(coords)
            labels = km.labels_
            # 只有 k>=2 时才能计算 silhouette
            if k > 1:
                score = silhouette_score(coords, labels)
                # print(f"[auto_k] K={k}, silhouette={score:.4f}")
                if score > best_score:
                    best_score, best_k = score, k
        n_clusters = best_k
        # print(f"[auto_k] 选择最佳簇数：K={n_clusters}, 得分={best_score:.4f}")

    # 4. 最终 KMeans
    kmeans = KMeans(n_clusters=n_clusters,
                    random_state=random_state).fit(coords)
    labels = kmeans.labels_
    centers = kmeans.cluster_centers_

    # 5. 根据 labels 构建簇
    clusters = {i: [] for i in range(n_clusters)}
    for cid, label in zip(customer_ids, labels):
        clusters[label].append(cid)

    # 6. 可视化
    if visualize or save_plot:
        plt.figure(figsize=(8, 6))
        # 客户：不同簇用不同颜色
        for i in range(n_clusters):
            mask = labels == i
            plt.scatter(coords[mask, 0],
                        coords[mask, 1],
                        label=f'Cluster {i}',
                        alpha=0.6)
        # 聚类中心
        plt.scatter(centers[:, 0],
                    centers[:, 1],
                    c='red',
                    marker='x',
                    s=100,
                    label='Centers')
        plt.title(f'K-Means Clustering (K={n_clusters})')
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.legend()
        plt.grid(True)

        if save_plot and plot_path:
            os.makedirs(os.path.dirname(plot_path), exist_ok=True)
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        if visualize:
            plt.show()
        plt.close()

    return clusters, centers

from .utils import load_instance


def assign_customers_to_centers(instance_name: str,
                                centers: np.ndarray,
                                radius: float,
                                customize_data: bool = False
                               ) -> dict[int, list[str]]:
    """
    给定聚类中心和覆盖半径，判断哪些客户点被哪个中心覆盖。

    :param instance_name: JSON 文件名（不含 .json）
    :param centers:       shape=(K,2) 的聚类中心坐标数组
    :param radius:        覆盖半径
    :param customize_data: 是否从 data/json_customize 读取
    :return:
      cover_map: dict
        key   = center_index (0 <= idx < K)
        value = List of customer keys (e.g. 'customer_1', 'customer_23') within radius
    """
    # 1. 定位并加载 JSON 实例
    json_dir = (os.path.join(BASE_DIR, 'data', 'json_customize')
                if customize_data
                else os.path.join(BASE_DIR, 'data', 'json'))
    json_file = os.path.join(json_dir, f'{instance_name}.json')
    instance = load_instance(json_file=json_file)
    if instance is None:
        # 文件不存在，所有站都空列表
        return {i: [] for i in range(centers.shape[0])}

    # 2. 提取所有 customer 坐标和对应 key
    pois = []
    keys = []
    for key, val in instance.items():
        if key.startswith('customer_'):
            keys.append(key)
            pois.append([val['coordinates']['x'], val['coordinates']['y']])
    pois = np.array(pois)  # shape=(n,2)
    if pois.size == 0:
        return {i: [] for i in range(centers.shape[0])}

    # 3. 为每个中心计算覆盖客户
    cover_map: dict[int, list[str]] = {}
    for idx, center in enumerate(centers):
        # 广播 pois - center => (n,2)，然后按行 L2 范数
        dists = np.linalg.norm(pois - center.reshape(1,2), axis=1)
        # 找出所有 dist <= radius 的索引
        covered_idxs = np.nonzero(dists <= radius)[0]
        # 用这些索引映射回 keys
        cover_map[idx] = [keys[i] for i in covered_idxs.tolist()]

    return cover_map

def greedy_select_stations(universe: set[str],
                           cover_map: dict[int, list[str]]
                          ) -> tuple[list[int], set[str]]:
    """
    用贪心算法选出最少的站点，使其覆盖 universe 中的所有元素（客户 key）。

    :param universe:   需要覆盖的客户 key 的集合，例如 {'customer_1', 'customer_2', ...}
    :param cover_map:  站点覆盖映射，key=站点索引, value=list of customer key
    :return:
      - selected: List[int] 被选中的站点索引（按选取顺序）
      - uncovered: set[str] 最终仍未被覆盖到的客户 key 集合（可能为空）
    """
    # 尚未被覆盖的客户集合
    uncovered = set(universe)
    # 最终选中的站点索引列表
    selected = []

    # 为了计算方便，将 cover_map 的值都转换为 set
    station_sets = {idx: set(custs) for idx, custs in cover_map.items()}

    # 只要还有客户未被覆盖，就继续选
    while uncovered:
        # 从所有站点中选一个，能覆盖最多“新”客户的那一个
        best_station, best_cov = max(
            station_sets.items(),
            key=lambda item: len(item[1] & uncovered),
            default=(None, set())
        )
        # 如果没有任何站可以再覆盖新客户，就跳出
        if best_station is None or len(best_cov & uncovered) == 0:
            break

        # 记录该站并更新 uncovered
        selected.append(best_station)
        uncovered -= station_sets[best_station]

    return selected

def find_uncovered_customers(covered_customers: set[str],
                             max_customers: int = 100
                            ) -> set[str]:
    """
    比较 covered_customers 和假定的全集（customer_1~customer_{max_customers}），
    返回那些未被覆盖的客户编号。

    :param covered_customers: 实际被覆盖到的客户 key 集合
    :param max_customers:     客户编号最大值（默认 100），全集为 customer_1...customer_max_customers
    :return: 未被覆盖的客户 key 集合
    """
    # 构造全集
    full_set = {f"customer_{i}" for i in range(1, max_customers + 1)}
    # 差集运算，得到未被覆盖的
    return full_set - covered_customers



