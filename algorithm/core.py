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
def ind2route_with_stations(individual, instance, energy_per_km: float = 1.0):
    """
    把一个访问顺序解码成多条子航线，允许途中去充电站充电。
    子航线里如果访问了 station，会在列表中插入 'station_<idx>' 字符串。
    支持：每连续飞行里程达到 max_cont_km(默认 3.6km) 时，自动“落地+起飞”，扣固定能量。
    """
    # —— 安全转换为 float 的小工具 —— #
    def f(x, default=0.0):
        try:
            return float(x)
        except (TypeError, ValueError):
            return default

    # 1) key 映射
    cust_keys = sorted([k for k in instance if k.startswith('customer_')],
                       key=lambda x: int(x.split('_')[1]))
    stat_keys = sorted([k for k in instance if k.startswith('station_')],
                       key=lambda x: int(x.split('_')[1]))
    nodes   = ['depart'] + cust_keys + stat_keys
    key2idx = {k: i for i, k in enumerate(nodes)}
    D       = instance['distance_matrix']
    speed        = f(instance.get('vehicle_speed', 0.0))
    max_e        = f(instance.get('vehicle_capacity', 0.0))
    charge_rt    = f(instance.get('vehicle_charge_rate', 0.0))
    landing_energy = f(instance.get('vehicle_vertical_energy_consume', 0.0))
    pay_cap      = f(instance.get('payload_capacity', instance.get('vehicle_initial_payload', 0.0)))
    max_cont_km  = f(instance.get('max_continuous_km', 3.6))  # 可从 JSON 覆盖

    k_payload = 0.3

    # 2) 单位里程能耗（随载荷）
    def rate(payload_now: float) -> float:
        frac = 0.0 if pay_cap <= 0 else max(0.0, payload_now) / pay_cap
        return energy_per_km * (1.0 + k_payload * frac)

    # 3) 每个 customer 到最近充电站的最短距离（若无站，则用回仓距离代替）
    stat_idxs = [key2idx[s] for s in stat_keys]
    nearest_station_dist = {}
    for ckey in cust_keys:
        cidx = key2idx[ckey]
        if stat_idxs:
            nearest_station_dist[cidx] = min(D[cidx][sidx] for sidx in stat_idxs)
        else:
            nearest_station_dist[cidx] = D[cidx][0]  # 没有充电站时，至少能回仓

    # 4) 计算一段距离内需要几次强制落地
    def forced_landings_needed(segment_len: float, since_touch: float) -> int:
        if max_cont_km <= 0:
            return 0
        return int((since_touch + max(0.0, segment_len)) // max_cont_km)

    # 5) 主循环状态
    e_left         = max_e
    t_elapsed      = 0.0
    last_idx       = 0      # depot
    payload_left   = pay_cap
    route, sub     = [], []
    since_touchdown = 0.0   # 自上次“落地”起累计的连续飞行里程

    for cid in individual:
        ckey = f'customer_{cid}'
        cidx = key2idx[ckey]

        # 尝试直飞 last -> customer
        d1         = D[last_idx][cidx]
        d2_depot   = D[cidx][0]                  # 时间参考
        d2_station = nearest_station_dist[cidx]  # 能量安全参考
        demand     = f(instance[ckey].get('demand', 0.0))
        svc        = f(instance[ckey].get('service_time', 0.0))

        # 直飞段需要的强制落地次数
        fl_d1 = forced_landings_needed(d1, since_touchdown)
        # 客户点服务后（已落地），再去最近站的强制落地次数
        fl_to_station_after = forced_landings_needed(d2_station, 0.0)

        # 直飞方案最低能量需求：保证“能到客户且还能到最近站”
        need_dir = (
            d1 * rate(payload_left) + fl_d1 * landing_energy +
            d2_station * rate(payload_left - demand) + fl_to_station_after * landing_energy
        )

        # 到达/返回时间（用于时间窗惩罚时打印用，这里只是保留，不参与能量判定）
        arr_t = t_elapsed + (d1 / speed if speed > 0 else 0.0)
        _ret_t = arr_t + svc + (d2_depot / speed if speed > 0 else 0.0)

        if (e_left >= need_dir) and (payload_left - demand >= 0):
            # 直飞可行：扣除当前段能量（含途中强制落地），落地服务
            sub.append(cid)
            e_left -= d1 * rate(payload_left) + fl_d1 * landing_energy

            if max_cont_km > 0:
                since_touchdown = (since_touchdown + d1) % max_cont_km
            # 到客户视为落地
            since_touchdown = 0.0

            t_elapsed   = arr_t + svc
            last_idx    = cidx
            payload_left -= demand
            continue

        # 直飞不行：尝试 last -> station -> customer 最短绕行
        stations_sorted = sorted(
            stat_keys,
            key=lambda skey: (D[last_idx][key2idx[skey]] + D[key2idx[skey]][cidx])
        )

        charged = False
        for skey in stations_sorted:
            sidx = key2idx[skey]
            d_ls = D[last_idx][sidx]
            d_sc = D[sidx][cidx]

            fl_ls = forced_landings_needed(d_ls, since_touchdown)
            fl_sc = forced_landings_needed(d_sc, 0.0)

            # 先保证能到站（到站就能充满）
            need_via_first_leg = d_ls * rate(payload_left) + fl_ls * landing_energy
            if (e_left >= need_via_first_leg) and (payload_left - demand >= 0):
                # last -> station
                e_left    -= d_ls * rate(payload_left) + fl_ls * landing_energy
                t_elapsed += (d_ls / speed if speed > 0 else 0.0)
                if max_cont_km > 0:
                    since_touchdown = (since_touchdown + d_ls) % max_cont_km
                since_touchdown = 0.0  # 到站视为落地
                last_idx = sidx
                sub.append(skey)

                # 充电到满
                need_charge = max(0.0, max_e - e_left)
                t_charge    = (need_charge / charge_rt) if (charge_rt > 0 and need_charge > 0) else 0.0
                t_elapsed  += t_charge
                e_left      = max_e

                # station -> customer
                e_left    -= d_sc * rate(payload_left) + fl_sc * landing_energy
                t_elapsed += (d_sc / speed if speed > 0 else 0.0)
                since_touchdown = 0.0  # 落地

                # 服务
                t_elapsed   += svc
                payload_left -= demand
                sub.append(cid)
                last_idx = cidx

                charged = True
                break

        if charged:
            continue

        # 仍然不行：收尾当前子航线，从仓库重新起飞一条
        if sub:
            route.append(sub)
        sub = [cid]

        d0c   = D[0][cidx]
        fl_d0c = forced_landings_needed(d0c, 0.0)
        e_left       = max_e - (d0c * rate(pay_cap) + fl_d0c * landing_energy)
        payload_left = pay_cap - demand
        t_elapsed    = (d0c / speed if speed > 0 else 0.0) + svc
        last_idx     = cidx
        since_touchdown = 0.0

    if sub:
        route.append(sub)
    return route

# 2) 对应的评价函数



def eval_vrptw_with_stations_chargingtime(individual, instance,
                             unit_cost=1.0, init_cost=0.0,
                             wait_cost=0.0, delay_cost=0.0,
                             energy_per_km=1.0,
                             decoder_fn=None):
    """
    评估带充电站的成本（能耗随载荷变化 + 强制降落/起飞规则）：
      - rate(payload) = energy_per_km * (1 + k_payload * (payload / payload_cap))
      - 每连续飞行 max_cont_km 需要执行 1 次“降落+起飞”：
          * 额外能量: landing_energy
          * 额外时间: landing_time
      - 到达 station/客户 都视为一次“落地”（连续里程清零）
    返回 (fitness,) ，fitness = 1 / total_cost
    """

    # —— 安全取 float —— #
    def f(x, default=0.0):
        try:
            return float(x)
        except (TypeError, ValueError):
            return default

    # 1) 全局参数
    speed        = f(instance.get('vehicle_speed', 0.0))
    bat_cap      = f(instance.get('battery_capacity', instance.get('vehicle_capacity', 0.0)))
    charge_rt    = f(instance.get('vehicle_charge_rate', 0.0))
    payload_cap  = f(instance.get('vehicle_initial_payload', 0.0))
    k_payload    = 0.3                                 # 满载比空载多的能耗比例
    energy_km    = f(energy_per_km, 1.0)

    # —— 强制降落参数 —— #
    max_cont_km  = f(instance.get('max_continuous_km', 3.6))          # 连续飞行阈值 (km)
    landing_energy = f(instance.get('vehicle_vertical_energy_consume', 0.0))  # 每次降落+起飞能量
    landing_time   = f(instance.get('vehicle_vertical_time_consume', 0.0))    # 每次降落+起飞时间(与speed同时间单位)

    # 2) 节点映射
    cust_keys = sorted([k for k in instance if k.startswith('customer_')],
                       key=lambda x: int(x.split('_')[1]))
    stat_keys = sorted([k for k in instance if k.startswith('station_')],
                       key=lambda x: int(x.split('_')[1]))
    nodes   = ['depart'] + cust_keys + stat_keys
    key2idx = {k: i for i, k in enumerate(nodes)}
    D       = instance['distance_matrix']

    # 3) 解码（支持自定义解码器）
    if decoder_fn is None:
        routes = ind2route_with_stations(individual, instance, energy_km)
    else:
        routes = decoder_fn(individual, instance, energy_km)

    # 4) 载荷相关单位里程能耗
    def rate(payload_now: float) -> float:
        frac = 0.0 if payload_cap <= 0 else max(0.0, payload_now) / payload_cap
        return energy_km * (1.0 + k_payload * frac)

    # —— 计算一段内需要触发几次“强制落地” —— #
    def forced_landings_needed(segment_len: float, since_touch: float) -> int:
        if max_cont_km <= 0:
            return 0
        return int((since_touch + max(0.0, segment_len)) // max_cont_km)

    total_cost = 0.0

    # 5) 逐条子航线评估
    for sub in routes:
        dist           = 0.0        # 距离成本
        time_pen       = 0.0        # 时间窗惩罚
        elapsed        = 0.0        # 子航线内累计时间
        remaining      = bat_cap    # 剩余电量
        payload_left   = payload_cap# 剩余载荷（子航线起始满载）
        last_idx       = key2idx['depart']
        since_touchdown = 0.0       # 自上次“落地”以来的连续飞行里程

        for node in sub:
            # 下一节点索引
            if isinstance(node, str):      # 'station_X'
                idx = key2idx[node]
                is_station = True
            else:                          # 客户 int
                idx = key2idx[f'customer_{node}']
                is_station = False

            # —— 飞行段 —— #
            d      = D[last_idx][idx]
            dist  += d

            # 计算本段需要的强制降落次数
            fl = forced_landings_needed(d, since_touchdown)

            # 按本段起飞时的载荷扣能量 + 强制降落能量
            remaining -= d * rate(payload_left) + fl * landing_energy

            # 飞行时间 + 强制降落时间
            t_fly   = (d / speed) if speed > 0 else 0.0
            arrival = elapsed + t_fly + fl * landing_time

            # 到达节点即视为“落地”，连续里程清零
            since_touchdown = 0.0

            if is_station:
                # —— 充电站（可选时间窗） —— #
                ready = f(instance[node].get('ready_time', 0.0))
                due   = f(instance[node].get('due_time', float('inf')))
                time_pen += wait_cost * max(ready - arrival, 0.0) \
                          + delay_cost * max(arrival - due, 0.0)

                # 充电至满（补满所需电量）
                need_charge = max(0.0, bat_cap - remaining)
                t_charge    = (need_charge / charge_rt) if (charge_rt > 0 and need_charge > 0) else 0.0
                elapsed     = arrival + t_charge
                remaining   = bat_cap      # 充满

                # 载荷不变
            else:
                # —— 客户 —— #
                key     = f'customer_{node}'
                ready   = f(instance[key].get('ready_time', 0.0))
                due     = f(instance[key].get('due_time', float('inf')))
                service = f(instance[key].get('service_time', 0.0))
                demand  = f(instance[key].get('demand', 0.0))

                time_pen += wait_cost * max(ready - arrival, 0.0) \
                          + delay_cost * max(arrival - due, 0.0)

                elapsed      = arrival + service
                payload_left = max(0.0, payload_left - demand)  # 服务后投放，载荷减少

            last_idx = idx

        # —— 回仓：按你的原逻辑，仅计距离成本（不计时间/电量/强制降落） —— #
        dist += D[last_idx][ key2idx['depart'] ]
        total_cost += init_cost + unit_cost * dist + time_pen

    if total_cost <= 0:
        return (float('inf'),)
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
                energy_per_km=energy_per_km
                )
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




def print_route_with_times_chargingtime(route, instance, energy_per_km):
    """
    打印每条子航线：到达时间 / 充电时间 / 剩余载荷 / SOC，
    并考虑单位里程能耗随载荷变化：
      rate = energy_per_km * (1 + k_payload * (payload_left / payload_cap))

    同时加入“最大连续飞行里程”强制落地逻辑：
      - 阈值：instance['max_continuous_km']（缺省 3.6）
      - 每触发一次落地，扣除 instance['vehicle_vertical_energy_consume'] 能量
      - 不额外计时（与解码器一致）

    额外汇总：
      - Total energy consumed
      - Mission completion time / makespan
    """
    # —— 安全取 float —— #
    def f(x, default=0.0):
        try:
            return float(x)
        except (TypeError, ValueError):
            return default

    D            = instance['distance_matrix']
    speed        = f(instance.get('vehicle_speed', 0.0))
    bat_cap      = f(instance.get('vehicle_capacity', 0.0))           # 电池容量
    charge_rt    = f(instance.get('vehicle_charge_rate', 0.0))
    payload_cap  = f(instance.get('vehicle_initial_payload', 0.0))    # 载荷上限
    landing_energy = f(instance.get('vehicle_vertical_energy_consume', 0.0))
    max_cont_km  = f(instance.get('max_continuous_km', 3.6))          # 最大连续飞行里程
    k_payload    = 0.3                                                 # 载荷能耗系数
    energy_km    = float(energy_per_km)                                # 空载/基准每公里能耗

    # 单位里程能耗率：随当前载荷线性变化
    def rate(payload_now: float) -> float:
        frac = 0.0 if payload_cap <= 0 else max(0.0, payload_now) / payload_cap
        return energy_km * (1.0 + k_payload * frac)

    # 计算一段距离内需要几次“强制落地”
    def forced_landings_needed(segment_len: float, since_touch: float) -> int:
        if max_cont_km <= 0:
            return 0
        return int((since_touch + max(0.0, segment_len)) // max_cont_km)

    # 0=depot, 1..=customers, then stations
    cust_ids    = sorted(int(k.split('_')[1]) for k in instance if k.startswith('customer_'))
    station_ids = sorted(int(k.split('_')[1]) for k in instance if k.startswith('station_'))
    node2idx    = {0: 0}
    for i, cid in enumerate(cust_ids, start=1):
        node2idx[cid] = i
    for j, sid in enumerate(station_ids, start=1 + len(cust_ids)):
        node2idx[f'station_{sid}'] = j

    total_charge_count = 0
    per_vehicle_counts = []
    total_energy       = 0.0
    mission_makespan   = 0.0

    for vid, sub in enumerate(route, start=1):
        print(f"Vehicle {vid}:")
        elapsed        = 0.0
        remaining      = bat_cap
        last_idx       = node2idx[0]
        payload_left   = payload_cap
        veh_energy     = 0.0
        charge_count   = 0
        since_touch    = 0.0  # 自上次“落地”以来已连续飞行里程

        print(f"  Depart depot at t=0.00, payload={payload_left}, SOC={remaining:.2f}/{bat_cap:.2f}")

        for node in sub:
            idx   = node2idx[node]
            dist  = D[last_idx][idx]
            t_fly = dist / speed if speed > 0 else 0.0

            # —— 本段强制落地次数 —— #
            fl = forced_landings_needed(dist, since_touch)
            extra_energy = fl * landing_energy

            # —— 本段水平飞行能量 —— #
            hop_energy  = dist * rate(payload_left)

            # —— 能量结算 —— #
            remaining   -= (hop_energy + extra_energy)
            total_energy += (hop_energy + extra_energy)
            veh_energy  += (hop_energy + extra_energy)

            # —— 连续里程更新（到达后先更新，再因“落地”清零）—— #
            if max_cont_km > 0:
                since_touch = (since_touch + dist) % max_cont_km

            arrival = elapsed + t_fly
            if remaining < -1e-9:
                print(f"    [WARN] energy below 0 after hop: {remaining:.2f}")

            # —— 打印强制落地信息 —— #
            if fl > 0:
                print(f"    [Forced landings: {fl}, extra_energy={extra_energy:.2f}]")

            if isinstance(node, int):
                # 客户节点：服务后投放，载荷减少；到达即“落地”
                demand  = f(instance[f'customer_{node}'].get('demand', 0.0))
                service = f(instance[f'customer_{node}'].get('service_time', 0.0))
                new_payload = payload_left - demand

                print(f"    -> Customer {node} at t={arrival:.2f}, "
                      f"drop={demand}, payload_left={max(new_payload,0):.2f}, "
                      f"SOC={remaining:.2f}/{bat_cap:.2f}")

                payload_left = new_payload
                elapsed      = arrival + service

                # 到达客户视为“落地”，清零连续里程
                since_touch = 0.0

            else:
                # 充电站节点：只充电，不改载荷；到达即“落地”
                sid = node.split('_', 1)[1]
                print(f"    -> Station {sid} at t={arrival:.2f}, "
                      f"payload_left={payload_left:.2f}, SOC={remaining:.2f}/{bat_cap:.2f}", end='')

                need_charge = max(0.0, bat_cap - remaining)
                if need_charge > 0 and charge_rt > 0:
                    t_charge = need_charge / charge_rt
                    charge_count += 1
                else:
                    t_charge = 0.0

                print(f", charging for {t_charge:.2f}")
                elapsed   = arrival + t_charge
                remaining = bat_cap

                # 到站视为“落地”，清零连续里程
                since_touch = 0.0

            last_idx = idx

        # —— 回仓 —— #
        back_dist = D[last_idx][ node2idx[0] ]
        t_back    = back_dist / speed if speed > 0 else 0.0

        # 回程强制落地
        fl_back        = forced_landings_needed(back_dist, since_touch)
        extra_back     = fl_back * landing_energy
        back_energy    = back_dist * rate(payload_left)

        remaining  -= (back_energy + extra_back)
        total_energy += (back_energy + extra_back)
        veh_energy   += (back_energy + extra_back)

        finish_t  = elapsed + t_back
        mission_makespan = max(mission_makespan, finish_t)

        if fl_back > 0:
            print(f"    [Forced landings on return: {fl_back}, extra_energy={extra_back:.2f}]")

        print(f"    -> Return depot at t={finish_t:.2f}, payload_left={payload_left:.2f}, "
              f"SOC={remaining:.2f}/{bat_cap:.2f}  |  energy_used={veh_energy:.2f}\n")

        per_vehicle_counts.append(charge_count)
        total_charge_count += charge_count

    vehicle_num = len(route)
    avg_charge_per_vehicle = total_charge_count / vehicle_num if vehicle_num > 0 else 0.0

    print(f"Total charging events: {total_charge_count}")
    print(f"Average charging events per vehicle: {avg_charge_per_vehicle:.2f}")
    print(f"Total energy consumed: {total_energy:.2f}")
    print(f"Mission completion time (makespan): {mission_makespan:.2f}")

    return {
        "total_energy": total_energy,
        "mission_makespan": mission_makespan,
        "total_charging_events": total_charge_count,
        "avg_charging_per_vehicle": avg_charge_per_vehicle,
    }


def ind2route_always_charge(individual, instance, energy_per_km: float = 1.0):
    """
    解码：每服务完一个客户，必定飞到“最近充电站”并充满，再继续下一个客户。
    若从当前位置无法满足“到客户 + 客户到最近站”的能量需求，或载荷不足，
    则结束当前子航线，开启新子航线（从仓库满载、满电出发）。

    路由表示：客户用 int，充电站用 'station_<id>' 字符串。
    能耗：rate(payload) = energy_per_km * (1 + k_payload * (payload/pay_cap)).
    """
    # -------- 安全取值并转成 float，避免 JSON 里是字符串导致类型错误 --------
    def f(x, default=0.0):
        try:
            return float(x)
        except (TypeError, ValueError):
            return default

    # 1) 节点映射
    cust_keys = sorted([k for k in instance if k.startswith('customer_')],
                       key=lambda x: int(x.split('_')[1]))
    stat_keys = sorted([k for k in instance if k.startswith('station_')],
                       key=lambda x: int(x.split('_')[1]))
    if not stat_keys:
        raise ValueError("ind2route_always_charge 需要至少 1 个 station_* 节点。")

    nodes   = ['depart'] + cust_keys + stat_keys
    key2idx = {k: i for i, k in enumerate(nodes)}

    D         = instance['distance_matrix']
    speed     = f(instance.get('vehicle_speed', 0.0))
    max_e     = f(instance.get('vehicle_capacity', instance.get('battery_capacity', 0.0)))
    charge_rt = f(instance.get('vehicle_charge_rate', 0.0))
    pay_cap   = f(instance.get('payload_capacity',
                     instance.get('vehicle_initial_payload', instance.get('payload', 0.0))))
    k_payload = 0.3
    base_rate = f(energy_per_km, 1.0)

    # 2) 载荷相关单位里程能耗
    def rate(payload_now: float) -> float:
        frac = 0.0 if pay_cap <= 0 else max(0.0, payload_now) / pay_cap
        return base_rate * (1.0 + k_payload * frac)

    # 3) 预计算：每个客户最近的充电站（索引与键名）
    stat_idxs = [key2idx[s] for s in stat_keys]
    nearest_station_idx = {}
    nearest_station_key = {}
    for ckey in cust_keys:
        cidx = key2idx[ckey]
        # 选距离最小的站
        sidx = min(stat_idxs, key=lambda si: D[cidx][si])
        nearest_station_idx[cidx] = sidx
        nearest_station_key[cidx] = nodes[sidx]  # 形如 'station_3'

    # 4) 主循环
    route, sub = [], []
    last_idx = 0        # 从仓库出发
    e_left   = max_e
    payload_left = pay_cap
    t_elapsed = 0.0     # 目前只累积，不做时间窗约束

    for cid in individual:
        ckey = f'customer_{cid}'
        cidx = key2idx[ckey]
        demand = f(instance[ckey]['demand'])

        if demand > pay_cap:
            raise ValueError(f"客户 {cid} 需求({demand}) 超过载荷上限({pay_cap})。")

        sidx = nearest_station_idx[cidx]        # 该客户最近站
        skey = nearest_station_key[cidx]

        # 定义一个尝试函数：从当前 last_idx 出发，去客户并再去最近站
        def can_do_from_here(e_now: float, payload_now: float) -> bool:
            d_last_c = D[last_idx][cidx]     # last -> customer
            d_c_s    = D[cidx][sidx]         # customer -> nearest station
            need     = d_last_c * rate(payload_now) + d_c_s * rate(payload_now - demand)
            return (e_now >= need) and (payload_now - demand >= 0)

        # 如果从当前基点/剩余能量无法完成“到客户+到站”，则开新子航线从仓库出发
        if not can_do_from_here(e_left, payload_left):
            if sub:
                route.append(sub)
            # 重置为新子航线（新车满载、满电、从仓库起飞）
            sub = []
            last_idx = 0
            e_left   = max_e
            payload_left = pay_cap
            t_elapsed = 0.0
            # 再次检查从仓库是否可行，不可行则认为参数设置不合理
            if not can_do_from_here(e_left, payload_left):
                d0c = D[0][cidx]; dcs = D[cidx][sidx]
                need_dbg = d0c * rate(pay_cap) + dcs * rate(pay_cap - demand)
                raise ValueError(
                    f"从仓库出发也无法完成 客户{cid} -> 最近站 的能量需求，"
                    f"需能量≈{need_dbg:.3f} > 电池容量 {max_e:.3f}。"
                )

        # —— 执行：last -> customer
        d_last_c = D[last_idx][cidx]
        e_left  -= d_last_c * rate(payload_left)
        t_elapsed += d_last_c / speed
        # 记录客户
        sub.append(cid)
        # 服务并投放
        service_time = f(instance[ckey]['service_time'])
        t_elapsed += service_time
        payload_left -= demand

        # —— 必定去最近充电站
        d_c_s = D[cidx][sidx]
        e_left  -= d_c_s * rate(payload_left)     # 注意：已投放后的载荷
        t_elapsed += d_c_s / speed
        # 插入站点标记
        sub.append(skey)
        last_idx = sidx

        # 充电至满
        need_charge = max_e - e_left
        t_charge    = need_charge / charge_rt if charge_rt > 0 else 0.0
        t_elapsed  += t_charge
        e_left      = max_e

        # 注意：不在站点补载荷；载荷只会递减到 0，如不足将触发新子航线

    # 收尾
    if sub:
        route.append(sub)
    return route

def run_gavrptw_always_charge(instance_name,
                              unit_cost, init_cost, wait_cost, delay_cost,
                              ind_size, pop_size, cx_pb, mut_pb, n_gen,
                              energy_per_km=1.0,
                              export_csv=False, customize_data=False):
    """
    GA 主流程（每服务完一个客户必去最近充电站充满的解码逻辑）
    - 评估函数：eval_vrptw_with_stations_chargingtime（能耗随载荷变化）
    - 解码器：ind2route_always_charge（客户 -> 最近站 -> 客户 -> 最近站 ...）
    """

    # --- 加载实例 ---
    json_dir  = os.path.join(BASE_DIR, 'data', 'json_customize' if customize_data else 'json')
    inst_file = os.path.join(json_dir, f'{instance_name}.json')
    instance  = load_instance(inst_file)
    if instance is None:
        print(f"❌ 找不到 {inst_file}")
        return

    # --- DEAP 注册（避免重复注册报错） ---
    try:
        creator.create('FitnessMax', base.Fitness, weights=(1.0,))
    except Exception:
        pass
    try:
        creator.create('Individual', list, fitness=creator.FitnessMax)
    except Exception:
        pass

    tb = base.Toolbox()
    tb.register('idx', random.sample, range(1, ind_size + 1), ind_size)
    tb.register('ind', tools.initIterate, creator.Individual, tb.idx)
    tb.register('pop', tools.initRepeat, list, tb.ind)

    # 评估：使用“必到站”的解码器，并考虑载荷影响能耗
    tb.register('evaluate', eval_vrptw_with_stations_chargingtime,
                instance=instance,
                unit_cost=unit_cost,
                init_cost=init_cost,
                wait_cost=wait_cost,
                delay_cost=delay_cost,
                energy_per_km=float(energy_per_km),
                decoder_fn=ind2route_always_charge)

    tb.register('select', tools.selRoulette)
    tb.register('mate', cx_partially_matched)
    tb.register('mutate', mut_inverse_indexes)

    # --- 初始种群与打分 ---
    pop = tb.pop(n=pop_size)
    for ind in pop:
        ind.fitness.values = tb.evaluate(ind)

    stats = []
    # --- 演化 ---
    for gen in range(n_gen):
        off = tb.select(pop, len(pop))
        off = list(map(tb.clone, off))

        # 交叉
        for c1, c2 in zip(off[::2], off[1::2]):
            if random.random() < cx_pb:
                tb.mate(c1, c2)
                del c1.fitness.values; del c2.fitness.values

        # 变异
        for m in off:
            if random.random() < mut_pb:
                tb.mutate(m)
                del m.fitness.values

        # 重新评估
        invalid = [i for i in off if not i.fitness.valid]
        for i in invalid:
            i.fitness.values = tb.evaluate(i)
        pop[:] = off

        if export_csv:
            fits = [i.fitness.values[0] for i in pop]
            avg  = sum(fits)/len(fits)
            std  = (sum(x*x for x in fits)/len(fits) - avg*avg)**0.5
            stats.append({'gen': gen, 'min': min(fits), 'max': max(fits), 'avg': avg, 'std': std})

    # --- 最优个体与展示（同一个解码器，保证一致） ---
    best = tools.selBest(pop, 1)[0]
    print(f"最佳染色体: {best}")
    print(f"适应度: {best.fitness.values[0]}")

    best_routes = ind2route_always_charge(best, instance, float(energy_per_km))
    print_route(best_routes)
    print_route_with_times_chargingtime(best_routes, instance, float(energy_per_km))
    plot_vrptw_routes_with_stations(instance_name, best_routes,
                                    customize=customize_data, save_image=True)

    # --- 导出 CSV ---
    if export_csv and stats:
        out = os.path.join(BASE_DIR, 'results',
                           f"{instance_name}_always_charge_ec{energy_per_km}.csv")
        make_dirs_for_file(out)
        with io.open(out, 'wt', encoding='utf-8', newline='') as f:
            w = DictWriter(f, fieldnames=list(stats[0].keys()))
            w.writeheader()
            for r in stats:
                w.writerow(r)

    return best, best_routes
