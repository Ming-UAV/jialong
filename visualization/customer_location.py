import os
import json
import matplotlib.pyplot as plt

from algorithm.core import plot_vrptw_instance

# 自动获取项目根目录
def main():
    plot_vrptw_instance('c105', customize=True,save_image=True)

if __name__ == '__main__':
    main()