
import os
import io
import fnmatch
from json import load, dump
from . import BASE_DIR
def make_dirs_for_file(path):
    '''gavrptw.uitls.make_dirs_for_file(path)'''
    try:
        os.makedirs(os.path.dirname(path))
    except OSError:
        pass


def guess_path_type(path):
    '''gavrptw.uitls.guess_path_type(path)'''
    if os.path.isfile(path):
        return 'File'
    if os.path.isdir(path):
        return 'Directory'
    if os.path.islink(path):
        return 'Symbolic Link'
    if os.path.ismount(path):
        return 'Mount Point'
    return 'Path'


def exist(path, overwrite=False, display_info=True):
    '''gavrptw.uitls.exist(path, overwrite=False, display_info=True)'''
    if os.path.exists(path):
        if overwrite:
            if display_info:
                print(f'{guess_path_type(path)}: {path} exists. Overwrite.')
            os.remove(path)
            return False
        if display_info:
            print(f'{guess_path_type(path)}: {path} exists.')
        return True
    if display_info:
        print(f'{guess_path_type(path)}: {path} does not exist.')
    return False


def load_instance(json_file):
    '''gavrptw.uitls.load_instance(json_file)'''

    if exist(path=json_file, overwrite=False, display_info=True):
        with io.open(json_file, 'rt', encoding='utf-8', newline='') as file_object:
            return load(file_object)
    return None


def merge_rules(rules):
    '''gavrptw.uitls.merge_rules(rules)'''
    is_fully_merged = True
    for round1 in rules:
        if round1[0] == round1[1]:
            rules.remove(round1)
            is_fully_merged = False
        else:
            for round2 in rules:
                if round2[0] == round1[1]:
                    rules.append((round1[0], round2[1]))
                    rules.remove(round1)
                    rules.remove(round2)
                    is_fully_merged = False
    return rules, is_fully_merged


def calculate_distance(customer1, customer2):
    '''gavrptw.uitls.calculate_distance(customer1, customer2)'''
    return ((customer1['coordinates']['x'] - customer2['coordinates']['x'])**2 + \
        (customer1['coordinates']['y'] - customer2['coordinates']['y'])**2)**0.5
def text2json(customize=False):
    '''gavrptw.uitls.text2json(customize=False)'''
    text_data_dir = os.path.join(BASE_DIR, 'data', 'text_customize' if customize else 'text')
    json_data_dir = os.path.join(BASE_DIR, 'data', 'json_customize' if customize else 'json')
    for text_file in map(lambda text_filename: os.path.join(text_data_dir, text_filename), \
        fnmatch.filter(os.listdir(text_data_dir), '*.txt')):
        json_data = {}
        with io.open(text_file, 'rt', encoding='utf-8', newline='') as file_object:
            for line_count, line in enumerate(file_object, start=1):
                if line_count in [2, 3, 4, 6, 7, 8, 9]:
                    pass
                elif line_count == 1:
                    # <Instance name>
                    json_data['instance_name'] = line.strip()
                elif line_count == 5:
                    # <Maximum vehicle number>, <Vehicle capacity>
                    values = line.strip().split()
                    json_data['max_vehicle_number'] = int(values[0])
                    json_data['vehicle_capacity'] = float(values[1])

                elif line_count == 10:
                    # Custom number = 0, depart
                    # <Custom number>, <X coordinate>, <Y coordinate>,
                    # ... <Demand>, <Ready time>, <Due date>, <Service time>
                    values = line.strip().split()
                    json_data['depart'] = {
                        'coordinates': {
                            'x': float(values[1]),
                            'y': float(values[2]),
                        },
                        'demand': float(values[3]),
                        'ready_time': float(values[4]),
                        'due_time': float(values[5]),
                        'service_time': float(values[6]),
                    }
                else:
                    # <Custom number>, <X coordinate>, <Y coordinate>,
                    # ... <Demand>, <Ready time>, <Due date>, <Service time>
                    values = line.strip().split()
                    json_data[f'customer_{values[0]}'] = {
                        'coordinates': {
                            'x': float(values[1]),
                            'y': float(values[2]),
                        },
                        'demand': float(values[3]),
                        'ready_time': float(values[4]),
                        'due_time': float(values[5]),
                        'service_time': float(values[6]),
                    }


        # 目前只包含 depart + customer
        customer_keys = [k for k in json_data if k.startswith('customer_')]
        # 需要把所有 station_* 也加进来
        station_keys = [k for k in json_data if k.startswith('station_')]

        # 然后 node_list = depart + customers + stations
        nodes = ['depart'] + sorted(customer_keys, key=lambda x: int(x.split('_')[1])) \
                + sorted(station_keys, key=lambda x: int(x.split('_')[1]))
        json_data['distance_matrix'] = [
            [calculate_distance(json_data[i], json_data[j]) for i in nodes]
            for j in nodes
        ]

        json_file_name = f"{json_data['instance_name']}.json"
        json_file = os.path.join(json_data_dir, json_file_name)
        print(f'Write to file: {json_file}')
        make_dirs_for_file(path=json_file)
        with io.open(json_file, 'wt', encoding='utf-8', newline='') as file_object:
            dump(json_data, file_object, sort_keys=True, indent=4, separators=(',', ': '))




def text3json(customize=False):
    '''gavrptw.utils.text2json(customize=False)'''
    text_data_dir = os.path.join(BASE_DIR, 'data', 'text_customize' if customize else 'text')
    json_data_dir = os.path.join(BASE_DIR, 'data', 'json_customize' if customize else 'json')

    for txt in fnmatch.filter(os.listdir(text_data_dir), '*.txt'):
        fp_txt = os.path.join(text_data_dir, txt)
        json_data = {}
        with io.open(fp_txt, 'rt', encoding='utf-8') as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                parts = line.split()
                head = parts[0].lower()

                # 1) Instance name
                if 'instance_name' not in json_data:
                    json_data['instance_name'] = line
                    continue

                # 2) Vehicle + battery info (两列数字)
                if head.isdigit() and len(parts) == 6 and 'max_vehicle_number' not in json_data:
                    json_data['max_vehicle_number'] = int(parts[0])
                    json_data['vehicle_capacity']   = float(parts[1])
                    json_data['vehicle_speed'] = float(parts[2])
                    json_data['vehicle_charge_rate'] = float(parts[3])
                    json_data['vehicle_vertical_energy_consume'] = float(parts[4])
                    json_data['vehicle_initial_payload'] =float(parts[5])


                    continue

                # 3) 充电站行： station ID X Y
                if head == 'station' and len(parts) >= 4:
                    sid = parts[1]
                    json_data[f'station_{sid}'] = {
                        'coordinates': {'x': float(parts[2]), 'y': float(parts[3])}
                    }
                    continue

                # 4) 仓库（depart）：customer编号0
                if head == '0' and len(parts) >= 7:
                    _, x, y, demand, ready, due, service = parts[:7]
                    json_data['depart'] = {
                        'coordinates': {'x': float(x), 'y': float(y)},
                        'demand':       float(demand),
                        'ready_time':   float(ready),
                        'due_time':     float(due),
                        'service_time': float(service),
                    }
                    continue

                # 5) 普通客户行：编号 1,2,3…
                if head.isdigit() and head != '0' and len(parts) >= 7:
                    cid = head
                    _, x, y, demand, ready, due, service = parts[:7]
                    json_data[f'customer_{cid}'] = {
                        'coordinates': {'x': float(x), 'y': float(y)},
                        'demand':       float(demand),
                        'ready_time':   float(ready),
                        'due_time':     float(due),
                        'service_time': float(service),
                    }
                    continue

                # 其余行一律跳过

        # 构造节点顺序：depart → customers → stations
        cust_keys    = sorted([k for k in json_data if k.startswith('customer_')],
                              key=lambda x: int(x.split('_')[1]))
        station_keys = sorted([k for k in json_data if k.startswith('station_')],
                              key=lambda x: int(x.split('_')[1]))
        nodes = ['depart'] + cust_keys + station_keys

        # 重新计算 distance_matrix（现在包含 station）
        json_data['distance_matrix'] = [
            [calculate_distance(json_data[i], json_data[j]) for i in nodes]
            for j in nodes
        ]

        # 写 JSON
        out_name = json_data['instance_name'] + '.json'
        fp_json  = os.path.join(json_data_dir, out_name)
        print(f'Write to file: {fp_json}')
        make_dirs_for_file(fp_json)
        with io.open(fp_json, 'wt', encoding='utf-8', newline='') as fo:
            dump(json_data, fo, sort_keys=True, indent=4, separators=(',',': '))
