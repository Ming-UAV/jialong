import numpy as np

from algorithm.greedy_cover import compute_station_poi_distances


def test_compute_station_poi_distances_basic():
    # 单个充电站
    station = [1, 1]
    # 三个兴趣点
    pois = [
        [1, 3],
        [2, 1],
        [3, 1],
    ]
    # 预期距离
    expected = np.array([1.0, 1.0, 2.0])
    result = compute_station_poi_distances(station, pois)
    assert isinstance(result, np.ndarray)
    assert result.shape == expected.shape
    assert np.allclose(result, expected)