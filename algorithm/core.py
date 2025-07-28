import os
import io
import random
import matplotlib.pyplot as plt
import json
from datetime import datetime
from csv import DictWriter
from deap import base, creator, tools
from . import BASE_DIR
from .utils import make_dirs_for_file, exist, load_instance, merge_rules

#编码流程，通过判断客户能否加入路线，从而将客户号码编入路线数组中
def ind2route(individual, instance):
    '''gavrptw.core.ind2route(individual, instance)'''
    #定义参数变量
    route = []
    vehicle_capacity = instance['vehicle_capacity']
    depart_due_time = instance['depart']['due_time']
    # Initialize a sub-route
    sub_route = []
    vehicle_load = 0
    elapsed_time = 0
    last_customer_id = 0
    #循环遍历
    for customer_id in individual:
        # Update vehicle load将id文本化，提取出需求，更新车辆载荷
        demand = instance[f'customer_{customer_id}']['demand']
        updated_vehicle_load = vehicle_load + demand
        # Update elapsed time提取车辆过去时间和返回仓库时间，更新总消耗时间
        service_time = instance[f'customer_{customer_id}']['service_time']
        return_time = instance['distance_matrix'][customer_id][0]
        #更新总消耗时间为之前消耗时间加上，从上个客户到目前客户时间和服务时间和返回时间
        updated_elapsed_time = elapsed_time + \
            instance['distance_matrix'][last_customer_id][customer_id] + service_time + return_time
        # Validate vehicle load and elapsed time载重不超以及总时间小于客户的要求时间情况下，将客户加入路线
        if (updated_vehicle_load <= vehicle_capacity) and (updated_elapsed_time <= depart_due_time):
            # Add to current sub-route
            sub_route.append(customer_id)
            vehicle_load = updated_vehicle_load
            elapsed_time = updated_elapsed_time - return_time
        else:
            # Save current sub-route
            route.append(sub_route)
            # Initialize a new sub-route and add to it
            sub_route = [customer_id]
            vehicle_load = demand
            elapsed_time = instance['distance_matrix'][0][customer_id] + service_time
        # Update last customer ID
        last_customer_id = customer_id
    if sub_route != []:
        # Save current sub-route before return if not empty
        route.append(sub_route)
    return route
def print_route(route, merge=False):
    """
    打印每条子航线，遇到 station_x 则用 Sx 的形式标注。
    """
    for vid, sub_route in enumerate(route, start=1):
        # 从仓库（0）出发
        labels = ['0']
        for node in sub_route:
            if isinstance(node, str) and node.startswith('station_'):
                # 站点：station_3 -> S3
                sid = node.split('_', 1)[1]
                labels.append(f'S{sid}')
            else:
                # 客户节点，直接是数字
                labels.append(str(node))
        labels.append('0')  # 回到仓库
        line = ' - '.join(labels)
        print(f"  Vehicle {vid}'s route: {line}")

#计算成本函数，先遍历子路径，再遍历子路径中每个客户的成本
def eval_vrptw(individual, instance, unit_cost=1.0, init_cost=0, wait_cost=0, delay_cost=0):
    '''gavrptw.core.eval_vrptw(individual, instance, unit_cost=1.0, init_cost=0, wait_cost=0,
        delay_cost=0)'''
    total_cost = 0
    route = ind2route(individual, instance)
    total_cost = 0
    for sub_route in route:
        sub_route_time_cost = 0
        sub_route_distance = 0
        elapsed_time = 0
        last_customer_id = 0
        for customer_id in sub_route:
            # Calculate section distance
            distance = instance['distance_matrix'][last_customer_id][customer_id]
            # Update sub-route distance
            sub_route_distance = sub_route_distance + distance
            # Calculate time cost
            arrival_time = elapsed_time + distance
            time_cost = wait_cost * max(instance[f'customer_{customer_id}']['ready_time'] - \
                arrival_time, 0) + delay_cost * max(arrival_time - \
                instance[f'customer_{customer_id}']['due_time'], 0)
            # Update sub-route time cost
            sub_route_time_cost = sub_route_time_cost + time_cost
            # Update elapsed time
            elapsed_time = arrival_time + instance[f'customer_{customer_id}']['service_time']
            # Update last customer ID
            last_customer_id = customer_id
        # Calculate transport cost
        sub_route_distance = sub_route_distance + instance['distance_matrix'][last_customer_id][0]
        sub_route_transport_cost = init_cost + unit_cost * sub_route_distance
        # Obtain sub-route cost
        sub_route_cost = sub_route_time_cost + sub_route_transport_cost
        # Update total cost
        total_cost = total_cost + sub_route_cost
    fitness = 1.0 / total_cost
    return (fitness, )
#交换部分基因
def cx_partially_matched(ind1, ind2):
    '''gavrptw.core.cx_partially_matched(ind1, ind2)'''
    cxpoint1, cxpoint2 = sorted(random.sample(range(min(len(ind1), len(ind2))), 2))
    part1 = ind2[cxpoint1:cxpoint2+1]
    part2 = ind1[cxpoint1:cxpoint2+1]
    rule1to2 = list(zip(part1, part2))
    is_fully_merged = False
    while not is_fully_merged:
        rule1to2, is_fully_merged = merge_rules(rules=rule1to2)
    rule2to1 = {rule[1]: rule[0] for rule in rule1to2}
    rule1to2 = dict(rule1to2)
    ind1 = [gene if gene not in part2 else rule2to1[gene] for gene in ind1[:cxpoint1]] + part2 + \
        [gene if gene not in part2 else rule2to1[gene] for gene in ind1[cxpoint2+1:]]
    ind2 = [gene if gene not in part1 else rule1to2[gene] for gene in ind2[:cxpoint1]] + part1 + \
        [gene if gene not in part1 else rule1to2[gene] for gene in ind2[cxpoint2+1:]]
    return ind1, ind2

#变异：逆序选中的部分基因片段
def mut_inverse_indexes(individual):
    '''gavrptw.core.mut_inverse_indexes(individual)'''
    start, stop = sorted(random.sample(range(len(individual)), 2))
    temp = individual[start:stop+1]
    temp.reverse()
    individual[start:stop+1] = temp
    return (individual, )


def run_gavrptw(instance_name, unit_cost, init_cost, wait_cost, delay_cost, ind_size, pop_size, \
    cx_pb, mut_pb, n_gen, export_csv=False, customize_data=False):
    '''gavrptw.core.run_gavrptw(instance_name, unit_cost, init_cost, wait_cost, delay_cost,
        ind_size, pop_size, cx_pb, mut_pb, n_gen, export_csv=False, customize_data=False)'''
    if customize_data:
        json_data_dir = os.path.join(BASE_DIR, 'data', 'json_customize')
    else:
        json_data_dir = os.path.join(BASE_DIR, 'data', 'json')
    json_file = os.path.join(json_data_dir, f'{instance_name}.json')
    instance = load_instance(json_file=json_file)
    if instance is None:
        return
    creator.create('FitnessMax', base.Fitness, weights=(1.0, ))
    creator.create('Individual', list, fitness=creator.FitnessMax)
    toolbox = base.Toolbox()
    # Attribute generator
    toolbox.register('indexes', random.sample, range(1, ind_size + 1), ind_size)
    # Structure initializers
    toolbox.register('individual', tools.initIterate, creator.Individual, toolbox.indexes)
    toolbox.register('population', tools.initRepeat, list, toolbox.individual)
    # Operator registering
    toolbox.register('evaluate', eval_vrptw, instance=instance, unit_cost=unit_cost, \
        init_cost=init_cost, wait_cost=wait_cost, delay_cost=delay_cost)
    toolbox.register('select', tools.selRoulette)
    toolbox.register('mate', cx_partially_matched)
    toolbox.register('mutate', mut_inverse_indexes)
    pop = toolbox.population(n=pop_size)
    # Results holders for exporting results to CSV file
    csv_data = []
    print('Start of evolution')
    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    print(f'  Evaluated {len(pop)} individuals')
    # Begin the evolution
    for gen in range(n_gen):
        print(f'-- Generation {gen} --')
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < cx_pb:
                toolbox.mate(child1, child2)
                del child1.fitness.values
                del child2.fitness.values
        for mutant in offspring:
            if random.random() < mut_pb:
                toolbox.mutate(mutant)
                del mutant.fitness.values
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        print(f'  Evaluated {len(invalid_ind)} individuals')
        # The population is entirely replaced by the offspring
        pop[:] = offspring
        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum([x**2 for x in fits])
        std = abs(sum2 / length - mean**2)**0.5
        print(f'  Min {min(fits)}')
        print(f'  Max {max(fits)}')
        print(f'  Avg {mean}')
        print(f'  Std {std}')
        # Write data to holders for exporting results to CSV file
        if export_csv:
            csv_row = {
                'generation': gen,
                'evaluated_individuals': len(invalid_ind),
                'min_fitness': min(fits),
                'max_fitness': max(fits),
                'avg_fitness': mean,
                'std_fitness': std,
            }
            csv_data.append(csv_row)
    print('-- End of (successful) evolution --')
    best_ind = tools.selBest(pop, 1)[0]
    print(f'Best individual: {best_ind}')
    print(f'Fitness: {best_ind.fitness.values[0]}')
    print_route(ind2route(best_ind, instance))
    print(f'Total cost: {1 / best_ind.fitness.values[0]}')
    if export_csv:
        csv_file_name = f'{instance_name}_uC{unit_cost}_iC{init_cost}_wC{wait_cost}' \
            f'_dC{delay_cost}_iS{ind_size}_pS{pop_size}_cP{cx_pb}_mP{mut_pb}_nG{n_gen}.csv'
        csv_file = os.path.join(BASE_DIR, 'results', csv_file_name)
        print(f'Write to file: {csv_file}')
        make_dirs_for_file(path=csv_file)
        if not exist(path=csv_file, overwrite=True):
            with io.open(csv_file, 'wt', encoding='utf-8', newline='') as file_object:
                fieldnames = [
                    'generation',
                    'evaluated_individuals',
                    'min_fitness',
                    'max_fitness',
                    'avg_fitness',
                    'std_fitness',
                ]
                writer = DictWriter(file_object, fieldnames=fieldnames, dialect='excel')
                writer.writeheader()
                for csv_row in csv_data:
                    writer.writerow(csv_row)
    est_ind = tools.selBest(pop, 1)[0]
    routes = ind2route(best_ind, instance)

    # 调用可视化路径函数
    plot_vrptw_routes(instance_name, routes, customize=customize_data)


import os
import json
import matplotlib.pyplot as plt
from datetime import datetime

def plot_vrptw_instance(instance_name, customize=True, save_image=True, timestamp=True):
    """
    可视化 VRPTW 问题的客户、仓库和充电站坐标分布图，并可选择保存图像。

    :param instance_name: JSON 文件的实例名（不包含 .json 后缀）
    :param customize:     是否从 json_customize 文件夹读取，默认 True
    :param save_image:    是否保存图片到 results 文件夹，默认 True
    :param timestamp:     是否在文件名中添加时间戳，避免覆盖历史记录，默认 True
    """
    BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    json_folder = 'json_customize' if customize else 'json'
    json_data_dir = os.path.join(BASE_DIR, 'data', json_folder)
    json_file = os.path.join(json_data_dir, f'{instance_name}.json')

    if not os.path.exists(json_file):
        print(f"❌ 文件不存在：{json_file}")
        return

    # 加载 JSON 数据
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)

    # 提取客户坐标
    customers_x, customers_y = [], []
    for key, val in data.items():
        if key.startswith('customer_'):
            customers_x.append(val['coordinates']['x'])
            customers_y.append(val['coordinates']['y'])

    # 提取充电站坐标
    stations_x, stations_y = [], []
    for key, val in data.items():
        if key.startswith('station_'):
            stations_x.append(val['coordinates']['x'])
            stations_y.append(val['coordinates']['y'])

    # 仓库坐标
    depot_x = data['depart']['coordinates']['x']
    depot_y = data['depart']['coordinates']['y']

    # 绘图
    plt.figure(figsize=(8, 6))
    # 客户
    plt.scatter(customers_x, customers_y,
                color='blue', marker='o', label='Customers')
    # 充电站
    if stations_x:
        plt.scatter(stations_x, stations_y,
                    color='green', marker='^', s=80, label='Stations')
    # 仓库
    plt.scatter(depot_x, depot_y,
                color='red', marker='s', s=100, label='Depot')

    # 添加客户编号标注
    for key, val in data.items():
        if key.startswith('customer_'):
            x, y = val['coordinates']['x'], val['coordinates']['y']
            plt.text(x + 0.3, y + 0.3, key.split('_')[1],
                     fontsize=9, color='blue')

    # 添加充电站编号标注
    for key, val in data.items():
        if key.startswith('station_'):
            x, y = val['coordinates']['x'], val['coordinates']['y']
            plt.text(x + 0.3, y + 0.3, key.split('_')[1],
                     fontsize=9, color='green')

    # 仓库标注
    plt.text(depot_x + 0.5, depot_y + 0.5,
             'Depot', fontsize=10, fontweight='bold', color='red')

    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title(f'VRPTW Instance: {instance_name}')
    plt.legend()
    plt.grid(True)

    # 保存图片
    if save_image:
        result_dir = os.path.join(BASE_DIR, 'results')
        os.makedirs(result_dir, exist_ok=True)

        if timestamp:
            time_str = datetime.now().strftime('%Y%m%d_%H%M%S')
            image_file = os.path.join(
                result_dir, f'{instance_name}_{time_str}.png')
        else:
            image_file = os.path.join(result_dir, f'{instance_name}.png')

        plt.savefig(image_file, dpi=300, bbox_inches='tight')
        print(f"✅ 图片已保存到：{image_file}")

    plt.show()


import os
import json
import matplotlib.pyplot as plt
from datetime import datetime


def plot_vrptw_routes(instance_name, routes, customize=True, save_image=True, timestamp=False):
    """
    可视化 VRPTW 问题的客户、仓库和车辆路径。

    :param instance_name: JSON 文件的实例名（不带 .json 后缀）
    :param routes: 路径数据，例如 [[3, 7, 1], [5, 9, 6], ...]
    :param customize: 是否使用 json_customize 目录，默认 True
    :param save_image: 是否保存图像，默认 True
    :param timestamp: 保存图像时是否添加时间戳，默认 False
    """
    BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    json_folder = 'json_customize' if customize else 'json'
    json_data_dir = os.path.join(BASE_DIR, 'data', json_folder)
    json_file = os.path.join(json_data_dir, f'{instance_name}.json')

    if not os.path.exists(json_file):
        print(f"❌ 文件不存在：{json_file}")
        return

    with open(json_file, 'r') as f:
        data = json.load(f)

    # 仓库坐标
    depot_x = data['depart']['coordinates']['x']
    depot_y = data['depart']['coordinates']['y']

    # 提取客户坐标
    customers_x, customers_y = [], []
    for key in data.keys():
        if key.startswith('customer_'):
            customers_x.append(data[key]['coordinates']['x'])
            customers_y.append(data[key]['coordinates']['y'])

    # 开始绘图
    plt.figure(figsize=(8, 6))
    plt.scatter(customers_x, customers_y, color='blue', marker='o', label='Customers')
    plt.scatter(depot_x, depot_y, color='red', marker='s', s=100, label='Depot')

    # 添加客户编号标注
    for key in data.keys():
        if key.startswith('customer_'):
            x = data[key]['coordinates']['x']
            y = data[key]['coordinates']['y']
            plt.text(x + 0.3, y + 0.3, key.split('_')[1], fontsize=9)

    plt.text(depot_x + 0.5, depot_y + 0.5, 'Depot', fontsize=10, fontweight='bold')

    # 路径颜色轮换
    colors = ['orange', 'green', 'purple', 'cyan', 'magenta', 'brown', 'pink', 'gray']

    # 绘制路径
    for idx, route in enumerate(routes):
        color = colors[idx % len(colors)]
        path_x = [depot_x]  # 从仓库出发
        path_y = [depot_y]
        for customer_id in route:
            key = f'customer_{customer_id}'
            path_x.append(data[key]['coordinates']['x'])
            path_y.append(data[key]['coordinates']['y'])
        path_x.append(depot_x)  # 回到仓库
        path_y.append(depot_y)

        plt.plot(path_x, path_y, color=color, linestyle='-', marker='>', label=f'Vehicle {idx + 1}')

    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title(f'VRPTW Routes: {instance_name}')
    plt.legend()
    plt.grid(True)

    # 保存图片
    if save_image:
        result_dir = os.path.join(BASE_DIR, 'results')
        os.makedirs(result_dir, exist_ok=True)
        if timestamp:
            time_str = datetime.now().strftime('%Y%m%d_%H%M%S')
            image_file = os.path.join(result_dir, f'{instance_name}_routes_{time_str}.png')
        else:
            image_file = os.path.join(result_dir, f'{instance_name}_routes.png')
        plt.savefig(image_file, dpi=300, bbox_inches='tight')
        print(f"✅ 路径图已保存到：{image_file}")

    plt.show()


import os
import io
import random
from csv import DictWriter
from deap import base, creator, tools
import matplotlib.pyplot as plt
import numpy as np

from .utils import load_instance, make_dirs_for_file, calculate_distance
from .core import print_route, plot_vrptw_routes

# 1) 带充电站的解码逻辑
def ind2route_with_stations(individual, instance, energy_per_km: float = 1.0,):
    """
    把一个访问顺序解码成多条子航线，允许途中去充电站充电。
    子航线里如果访问了 station，会在列表中插入 'station_<idx>' 字符串。
    """
    # 1. 恢复 key -> 索引映射
    cust_keys = sorted([k for k in instance if k.startswith('customer_')],
                       key=lambda x: int(x.split('_')[1]))
    stat_keys = sorted([k for k in instance if k.startswith('station_')],
                       key=lambda x: int(x.split('_')[1]))
    nodes = ['depart'] + cust_keys + stat_keys
    key2idx = {k:i for i,k in enumerate(nodes)}
    D = instance['distance_matrix']
    speed = instance['vehicle_speed']  # 直接读 JSON 里定义的 speed
    max_e = instance['vehicle_capacity']
    charge_rt = instance['vehicle_charge_rate']
    degr_rt = instance['vehicle_degradation_rate']

    e_left = max_e
    t_elapsed = 0.0
    last_idx = 0   # depot == index 0

    route = []
    sub = []
    # —— 0) 先算出每个 customer 返回“最近充电站”的最小距离 —— #
    # stat_idxs 是所有 station 在矩阵里的索引
    stat_idxs = [key2idx[s] for s in stat_keys]
    # nearest_station_dist[cidx] = min distance from cidx 到任意 station
    nearest_station_dist = {}
    for ckey in cust_keys:
        cidx = key2idx[ckey]
        nearest_station_dist[cidx] = min(D[cidx][sidx] for sidx in stat_idxs)
    for cid in individual:
        ckey = f'customer_{cid}'
        cidx = key2idx[ckey]

        # —— 1) 尝试直飞 last→cust→depot —— #
        d1 = D[last_idx][cidx]
        d2 = D[cidx][0]
        d2_station = nearest_station_dist[cidx]
        need_dir = (d1 + d2_station) * energy_per_km
        arr_t = t_elapsed + d1/speed
        ret_t = arr_t + instance[ckey]['service_time'] + d2/speed

        if e_left >= need_dir :
            sub.append(cid)
            e_left -= need_dir
            t_elapsed = arr_t + instance[ckey]['service_time']
            last_idx = cidx
            continue
        stations_sorted = sorted(
            stat_keys,
            key=lambda skey: (
                    D[last_idx][key2idx[skey]] +
                    D[key2idx[skey]][cidx]
            )
        )
        # —— 2) 直飞不行，尝试各充电站 —— #
        charged = False
        for skey in stations_sorted:
            sidx = key2idx[skey]
            d_ls = D[last_idx][sidx]
            d_sc = D[sidx][cidx]
            need_via = (d_ls ) * energy_per_km
            arr_v = t_elapsed + d_ls/speed + d_sc/speed
            ret_v = arr_v + instance[ckey]['service_time'] + d2/speed

            if e_left >= need_via:
                # 1) 飞到充电站
                e_left -= need_via
                t_elapsed += d_ls / speed
                last_idx = sidx

                # 2) 动态充电：把电充满
                need_charge = max_e - e_left
                t_charge = need_charge / charge_rt
                t_elapsed += t_charge
                e_left = max_e

                # 3) 插入站点标记
                sub.append(skey)

                # 4) 飞到客户并服务
                e_left -= d_sc * energy_per_km
                t_elapsed += d_sc / speed
                t_elapsed += instance[ckey]['service_time']
                sub.append(cid)
                last_idx = cidx

                charged = True
                break

        if charged:
            continue

        # —— 3) 都不行，则收尾当前子航线，换电池新航线 —— #
        if sub:
            route.append(sub)
        sub = [cid]
        d0c = D[0][cidx]
        e_left = max_e - d0c*energy_per_km
        t_elapsed = d0c/speed + instance[ckey]['service_time']
        last_idx = cidx

    if sub:
        route.append(sub)
    return route

# 2) 对应的评价函数
def eval_vrptw_with_stations(individual, instance,
                             unit_cost=1.0, init_cost=0.0,
                             wait_cost=0.0, delay_cost=0.0,
                             energy_per_km=1.0):
    speed = instance.get('vehicle_speed', 1.0)
    # 1) 先重建 key->idx 的映射
    cust_keys = sorted([k for k in instance if k.startswith('customer_')],
                       key=lambda x: int(x.split('_')[1]))
    stat_keys = sorted([k for k in instance if k.startswith('station_')],
                       key=lambda x: int(x.split('_')[1]))
    nodes = ['depart'] + cust_keys + stat_keys
    key2idx = {k: i for i, k in enumerate(nodes)}

    routes = ind2route_with_stations(individual, instance, energy_per_km)
    total = 0.0

    for sub in routes:
        dist = 0.0
        time_pen = 0.0
        elapsed = 0.0
        last_idx = key2idx['depart']

        for node in sub:
            # 2) 根据 node 的类型，查映射表拿 idx
            if isinstance(node, str) and node.startswith('station_'):
                idx = key2idx[node]
            else:
                # node 是客户的数字，比如 3，就对应 'customer_3'
                idx = key2idx[f'customer_{node}']

            # 3) 用 idx 去 distance_matrix 里取距离
            d = instance['distance_matrix'][last_idx][idx]
            dist += d
            t_fly = d / speed
            arr = elapsed + t_fly


            # 4) 只有 customer 节点才有时间窗惩罚/服务时间
            if not isinstance(node, str):
                key = f'customer_{node}'
                time_pen += (
                    wait_cost * max(instance[key]['ready_time'] - arr, 0)
                    + delay_cost * max(arr - instance[key]['due_time'], 0)
                )
                elapsed = arr + instance[key]['service_time']
            else:
                # 充电站只算飞行时间，不算服务
                elapsed = arr

            last_idx = idx

        # 回到 depot
        dist += instance['distance_matrix'][last_idx][key2idx['depart']]
        total += init_cost + unit_cost * dist + time_pen

    return (1.0 / total, )




def eval_vrptw_with_stations_chargingtime (individual, instance,
                             unit_cost=1.0, init_cost=0.0,
                             wait_cost=0.0, delay_cost=0.0,
                             energy_per_km=1.0):
    """
    对带充电站路径进行成本评估。
    与 eval_vrptw 不同之处在于：
      - 路径中已包含充电站节点（由 ind2route_with_stations 插入）
      - 飞行段时间 = distance / vehicle_speed
      - 充电站节点需计算充电时间并累加到 elapsed
      - 仅对客户节点计算软时间窗惩罚
      - 最终成本 = init_cost * (#子航线) + unit_cost * 总距离 + time_penalty
    返回： (fitness,) ，fitness = 1/total_cost
    """

    # 1. 全局参数读取
    speed     = instance['vehicle_speed']               # 车辆飞行速度
    bat_cap   = instance['vehicle_capacity']         # 电池满电量
    charge_rt = instance['vehicle_charge_rate']      # 充电速率（电量单位/时间单位）

    # 2. 构建节点到距离矩阵索引的映射
    cust_keys = sorted([k for k in instance if k.startswith('customer_')],
                       key=lambda x: int(x.split('_')[1]))
    stat_keys = sorted([k for k in instance if k.startswith('station_')],
                       key=lambda x: int(x.split('_')[1]))
    nodes     = ['depart'] + cust_keys + stat_keys
    key2idx   = {k: i for i, k in enumerate(nodes)}
    D         = instance['distance_matrix']

    # 3. 解码：获得含充电站的子航线列表
    routes = ind2route_with_stations(individual, instance, energy_per_km)

    total_cost = 0.0

    # 4. 遍历每条子航线进行成本计算
    for sub in routes:
        dist      = 0.0     # 累计飞行距离成本部分
        time_pen  = 0.0     # 累计客户时间窗惩罚
        elapsed   = 0.0     # 子航线内部的时序进度
        remaining = bat_cap # 当前剩余电量
        last_idx  = key2idx['depart']


        for node in sub:
            # —— 飞行段 —— #
            # 找到下一个节点在矩阵中的索引
            if isinstance(node, str):
                # node 本身就是 'station_X'
                idx = key2idx[node]
            else:
                # node 是客户编号的整数，要拼出 'customer_<node>'
                idx = key2idx[f'customer_{node}']
            # 读取两点距离
            d   = D[last_idx][idx]
            dist += d
            # 耗电 = 距离 * 单位耗电率
            remaining -= d * energy_per_km
            # 飞行时间 = 距离 / 速度
            t_fly   = d / speed
            arrival = elapsed + t_fly

            if isinstance(node, str):
                # —— 充电站节点 —— #
                # 1) （可选）时间窗惩罚
                ready = instance[node].get('ready_time', 0.0)
                due   = instance[node].get('due_time', float('inf'))
                time_pen += wait_cost * max(ready - arrival, 0) \
                          + delay_cost * max(arrival - due, 0)

                # 2) 动态充电：补满电所需时间
                need_charge = bat_cap - remaining
                t_charge    = need_charge / charge_rt if charge_rt > 0 else 0.0

                # 累加充电时间
                elapsed   = arrival + t_charge
                # 电量回满
                remaining = bat_cap

            else:
                # —— 客户节点 —— #
                # 1) 时间窗软惩罚
                key = f'customer_{node}'
                ready = instance[key]['ready_time']
                due   = instance[key]['due_time']
                time_pen += wait_cost * max(ready - arrival, 0) \
                          + delay_cost * max(arrival - due, 0)

                # 2) 服务时间累加
                service = instance[key]['service_time']
                elapsed = arrival + service

            # 更新当前位置
            last_idx = idx

        # —— 回仓段 —— #
        # 回到仓库的距离成本
        dist += D[last_idx][ key2idx['depart'] ]
        # 子航线总成本 = init_cost + unit_cost*dist + time_pen
        total_cost += init_cost + unit_cost * dist + time_pen

    # 防除零
    if total_cost <= 0:
        return (float('inf'),)
    # 适应度 = 1 / 总成本
    return (1.0 / total_cost,)


# 3) GA 主流程改造
def run_gavrptw_with_stations(instance_name,
                              unit_cost, init_cost, wait_cost, delay_cost,
                              ind_size, pop_size, cx_pb, mut_pb, n_gen,
                              energy_per_km=1.0,
                              export_csv=False, customize_data=False):
    # load JSON
    json_dir = os.path.join(BASE_DIR, 'data',
                            'json_customize' if customize_data else 'json')
    inst_file = os.path.join(json_dir, f'{instance_name}.json')
    instance = load_instance(inst_file)
    if instance is None:
        print(f"❌ 找不到 {inst_file}")
        return

    # DEAP setup
    creator.create('FitnessMax', base.Fitness, weights=(1.0,))
    creator.create('Individual', list, fitness=creator.FitnessMax)
    tb = base.Toolbox()
    tb.register('idx', random.sample, range(1, ind_size+1), ind_size)
    tb.register('ind', tools.initIterate, creator.Individual, tb.idx)
    tb.register('pop', tools.initRepeat, list, tb.ind)

    tb.register('evaluate', eval_vrptw_with_stations_chargingtime,
                instance=instance,
                unit_cost=unit_cost,
                init_cost=init_cost,
                wait_cost=wait_cost,
                delay_cost=delay_cost,
                energy_per_km=energy_per_km)
    tb.register('select', tools.selRoulette)
    tb.register('mate', cx_partially_matched)
    tb.register('mutate', mut_inverse_indexes)

    pop = tb.pop(n=pop_size)
    # 初评
    for ind in pop:
        ind.fitness.values = tb.evaluate(ind)

    stats = []
    for gen in range(n_gen):
        off = tb.select(pop, len(pop))
        off = list(map(tb.clone, off))
        # crossover + mutate
        for c1, c2 in zip(off[::2], off[1::2]):
            if random.random() < cx_pb:
                tb.mate(c1, c2)
                del c1.fitness.values; del c2.fitness.values
        for m in off:
            if random.random() < mut_pb:
                tb.mutate(m)
                del m.fitness.values
        # re-evaluate
        invalid = [i for i in off if not i.fitness.valid]
        for i in invalid:
            i.fitness.values = tb.evaluate(i)
        pop[:] = off

        if export_csv:
            fits = [i.fitness.values[0] for i in pop]
            stats.append({
                'gen':gen,'min':min(fits),'max':max(fits),
                'avg':sum(fits)/len(fits),
                'std':(sum(x*x for x in fits)/len(fits)- (sum(fits)/len(fits))**2)**0.5
            })

    best = tools.selBest(pop, 1)[0]
    print(f"最佳染色体: {best}")
    print(f"适应度: {best.fitness.values[0]}")
    best_routes = ind2route_with_stations(best, instance, energy_per_km)
    print_route(best_routes)
    print_route_with_times_chargingtime(best_routes, instance,energy_per_km)
    plot_vrptw_routes_with_stations(instance_name, best_routes, customize=customize_data,save_image=True)


    # 导出 CSV
    if export_csv and stats:
        out = os.path.join(BASE_DIR,'results',
                           f"{instance_name}_stations_ec{energy_per_km}.csv")
        make_dirs_for_file(out)
        with io.open(out,'wt',encoding='utf-8',newline='') as f:
            w = DictWriter(f, fieldnames=list(stats[0].keys()))
            w.writeheader()
            for r in stats:
                w.writerow(r)
import os, json
import matplotlib.pyplot as plt
from datetime import datetime

def plot_vrptw_routes_with_stations(instance_name,
                                    routes,
                                    customize=True,
                                    save_image=True,
                                    timestamp=False):
    BASE_DIR = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    folder = 'json_customize' if customize else 'json'
    json_file = os.path.join(BASE_DIR, 'data', folder, f"{instance_name}.json")
    with open(json_file, 'r') as f:
        data = json.load(f)

    # 提取坐标
    depot = data['depart']['coordinates']
    # 客户坐标，key 用 int(customer_id)
    customers = {
        int(k.split('_')[1]): v['coordinates']
        for k, v in data.items() if k.startswith('customer_')
    }
    # 充电站坐标，同理
    stations = {
        int(k.split('_')[1]): v['coordinates']
        for k, v in data.items() if k.startswith('station_')
    }

    plt.figure(figsize=(8,6))

    # 仓库
    plt.scatter(depot['x'], depot['y'], c='red', marker='s', s=100, label='Depot')
    plt.text(depot['x']+0.2, depot['y']+0.2, 'Depot')

    # 客户
    for cid, coord in customers.items():
        plt.scatter(coord['x'], coord['y'], c='blue', marker='o')
        plt.text(coord['x']+0.2, coord['y']+0.2, str(cid), fontsize=8)
    plt.scatter([], [], c='blue', marker='o', label='Customer')

    # 充电站
    for sid, coord in stations.items():
        plt.scatter(coord['x'], coord['y'], c='green', marker='^', s=80)
        plt.text(coord['x']+0.2, coord['y']+0.2, f"S{sid}", fontsize=8)
    plt.scatter([], [], c='green', marker='^', label='Station')

    # 绘制每条子航线
    colors = ['orange','purple','cyan','magenta','brown','pink','gray']
    for i, sub in enumerate(routes):
        path_x = []
        path_y = []
        # 从 depot 出发
        path_x.append(depot['x']); path_y.append(depot['y'])
        for node in sub:
            if isinstance(node, str) and node.startswith('station_'):
                sid = int(node.split('_')[1])
                coord = stations[sid]
            else:
                cid = int(node)
                coord = customers[cid]
            path_x.append(coord['x']); path_y.append(coord['y'])
        # 回 depot
        path_x.append(depot['x']); path_y.append(depot['y'])
        plt.plot(path_x, path_y,
                 color=colors[i % len(colors)],
                 linestyle='-',
                 marker='>',
                 label=f'Vehicle {i+1}')

    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.title(f'Routes w/ Stations (K={len(routes)})')
    plt.legend(loc='best')
    plt.grid(True)

    if save_image:
        outdir = os.path.join(BASE_DIR, 'results')
        os.makedirs(outdir, exist_ok=True)
        # 始终生成时间戳，格式 YYYYmmdd_HHMMSS
        time_str = datetime.now().strftime('%Y%m%d_%H%M%S')
        # 在文件名里加个下划线分隔
        filename = f"{instance_name}_routes_{time_str}.png"
        filepath = os.path.join(outdir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"✅ 路径图已保存到：{filepath}")
    plt.close()



def print_route_with_times(route, instance):
    """
    打印每条子航线，并显示到达每个客户/站点时刻，
    并且正确区分客户编号和站点编号在 distance_matrix 中的行列位置。
    同时将距离除以 instance['vehicle_speed'] 得到的时间用来计算 arrival。
    """
    D = instance['distance_matrix']
    speed = instance.get('vehicle_speed', 1.0)  # 默认 1

    # 1) 构建节点顺序：0（仓库） + 所有客户升序 + 所有站点升序
    cust_ids    = sorted(int(k.split('_')[1]) for k in instance if k.startswith('customer_'))
    station_ids = sorted(int(k.split('_')[1]) for k in instance if k.startswith('station_'))
    node2idx = {0: 0}
    # 客户从 1 开始
    for i, cid in enumerate(cust_ids, start=1):
        node2idx[cid] = i
    # 站点接着客户编号后面
    for j, sid in enumerate(station_ids, start=1 + len(cust_ids)):
        node2idx[f'station_{sid}'] = j

    for vid, sub in enumerate(route, start=1):
        print(f"Vehicle {vid}:")
        elapsed   = 0.0
        last_idx  = node2idx[0]  # depot 的矩阵索引
        print(f"  Depart depot at t=0.00")
        for node in sub:
            idx      = node2idx[node]
            distance = D[last_idx][idx]
            t_travel = distance / speed      # ← 用速度换算时间
            arrival  = elapsed + t_travel

            if isinstance(node, int):
                # 客户节点
                print(f"    -> Customer {node} at t={arrival:.2f}")
                service = instance[f'customer_{node}']['service_time']
                elapsed = arrival + service
            else:
                # 站点节点 'station_X'
                sid = node.split('_', 1)[1]
                print(f"    -> Station {sid} at t={arrival:.2f}")
                # 如果站点存在固定充电时间，则加上：
                # charge = instance[f'station_{sid}'].get('charge_time', 0.0)
                # elapsed = arrival + charge
                elapsed = arrival

            last_idx = idx

        # 回到 depot
        back_dist = D[last_idx][ node2idx[0] ]
        t_back    = back_dist / speed        # ← 同样要除以速度
        print(f"    -> Return depot at t={elapsed + t_back:.2f}\n")


def print_route_with_times_chargingtime (route, instance,energy_per_km):
    """
    打印每条子航线，并显示到达每个客户/站点时刻，
    正确区分客户编号和站点编号，并动态加上充电时间。
    """
    D          = instance['distance_matrix']
    speed      = instance['vehicle_speed']
    bat_cap    = instance['vehicle_capacity']
    charge_rt  = instance['vehicle_charge_rate']
    energy_km = energy_per_km

    # 构建矩阵索引映射：0=depot, 1..=customers, then stations
    cust_ids    = sorted(int(k.split('_')[1]) for k in instance if k.startswith('customer_'))
    station_ids = sorted(int(k.split('_')[1]) for k in instance if k.startswith('station_'))
    node2idx    = {0: 0}
    for i, cid in enumerate(cust_ids, start=1):
        node2idx[cid] = i
    for j, sid in enumerate(station_ids, start=1 + len(cust_ids)):
        node2idx[f'station_{sid}'] = j

    for vid, sub in enumerate(route, start=1):
        print(f"Vehicle {vid}:")
        elapsed   = 0.0
        remaining = bat_cap
        last_idx  = node2idx[0]  # depot index
        print(f"  Depart depot at t=0.00")

        for node in sub:
            idx      = node2idx[node]
            # 1) 飞行到下一个节点
            dist     = D[last_idx][idx]
            t_fly    = dist / speed
            arrival  = elapsed + t_fly
            remaining -= dist * energy_km

            if isinstance(node, int):
                # 客户节点：打印到达时间 + 服务
                print(f"    -> Customer {node} at t={arrival:.2f}")
                service = instance[f'customer_{node}']['service_time']
                elapsed = arrival + service
            else:
                # 充电站节点：打印到达时间
                sid = node.split('_', 1)[1]
                print(f"    -> Station {sid} at t={arrival:.2f}", end='')

                # 2) 动态计算充电时间
                need_charge = bat_cap - remaining
                t_charge    = need_charge / charge_rt if charge_rt > 0 else 0.0


                # 打印充电时长
                print(f", charging for {t_charge:.2f}")
                # 累加充电时间并恢复电量
                elapsed   = arrival + t_charge
                remaining = bat_cap

            last_idx = idx

        # 返回 depot
        back_dist = D[last_idx][ node2idx[0] ]
        t_back    = back_dist / speed
        print(f"    -> Return depot at t={elapsed + t_back:.2f}\n")
