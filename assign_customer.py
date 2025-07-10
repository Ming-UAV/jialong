from algorithm.greedy_cover import cluster_customers, assign_customers_to_centers, greedy_select_stations,find_uncovered_customers
# 先聚类
clusters, centers = cluster_customers(
    instance_name="c105",
    n_clusters=range(2,15),
    auto_k= True,
    customize_data=True,
    random_state=42,
    visualize = True,
)

# 再根据半径 radius生成覆盖映射
cover_map = assign_customers_to_centers(
    instance_name="c105",
    centers=centers,
    radius=30,#station的覆盖半径，根据电池电量改
    customize_data=True
)

for c_idx, customers in cover_map.items():
    print(f"中心 {c_idx} 覆盖客户：", customers)
# 准备所有客户 key,遍历每个集合点的值，作并集，从而把能去的兴趣点存入allcustomer数组,但是这样会导致，没有空集
all_customers = set()
for custs in cover_map.values():
    all_customers.update(custs)

# 贪心选站
best_stations = greedy_select_stations(all_customers, cover_map)
print("Selected station indices:", best_stations)
uncovered = find_uncovered_customers(all_customers, max_customers=100)

print("以下客户未被覆盖到：", sorted(uncovered))
#打印集群坐标
for idx, (cx, cy) in enumerate(centers, start=1):
    print(f"station {idx}  {cx:.4f}   {cy:.4f}")

print("\n>>> 下面是被选中的站点：")
filtered = [oid for oid in best_stations if oid != 0]
for new_id, orig_id in enumerate(filtered, start=1):
    # orig_id-1 因为 centers 是 0-based 的
    cx, cy = centers[orig_id-1]
    print(f"  station {new_id}  {cx:.4f}   {cy:.4f}")

