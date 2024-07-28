import csv
import random
from datetime import datetime

import numpy as np

import dem_data
from main import analysis_by_dpderl_simplified, analysis_by_r3, judge_is_out_bound

data = []
terrain_type_dict = {
    "data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif": "山区",
    "data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif": "平原",
    "data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif": "丘陵"
}


def analysis_accuracy():
    print(f"精度测试运行起始时间:{datetime.now()}")
    for i in range(0, 5000 + 1, 10):
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif", 30)
        analysis_accuracy_repeat(dem, i)
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif", 30)
        analysis_accuracy_repeat(dem, i)
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif", 30)
        analysis_accuracy_repeat(dem, i)


"""
dem: 选择的地形文件
r: 实际半径,因此需转化为格网半径
h: 观察点距离海平面高度
limit: 同一高度测多少次
"""


def analysis_accuracy_repeat(dem, i):
    global data

    data_line = []
    # 生成一个从10到20的随机整数
    lon = 0
    lat = 0
    is_out_bound = True
    r_distance = 2000
    while is_out_bound:
        random_integer_x = random.randint(1, 3600)
        random_integer_y = random.randint(1, 3600)
        lon = dem.start_x + random_integer_x * dem.dx
        lat = dem.start_y + random_integer_y * dem.dy
        is_out_bound = judge_is_out_bound(lon, lat, r_distance, dem)
    start_angle = random.uniform(0, 360)  # 生成起始角度，范围在 0 到 360 之间
    end_angle = random.uniform(start_angle, 360)  # 生成终止角度，确保大于起始角度且在 0 到 360 之间
    h_center = dem.height[int((lon - dem.start_x) / dem.dx)][int((lat - dem.start_y) / dem.dy)]
    # SPDERL
    # if dem.file_path == "data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif":
    #     print("----------0.0")
    #     lon = 114.2806944
    #     lat = 34.29569444
    #     h_center = dem.height[int((lon - dem.start_x) / dem.dx)][int((lat - dem.start_y) / dem.dy)]
    #     start_angle = 72.02715637
    #     end_angle = 338.661188
    #     i = 10
    print(
        f"第{i + 1}次测试:\ndem.file_path:{dem.file_path}\nlon:{lon}\nlat:{lat}\nh_center:{h_center}\nh:{i}\nstart_angle:{start_angle}\nend_angle:{end_angle}\nr:{r_distance}")
    result_spderl, actually_count, result_opt_rl_spderl = analysis_by_dpderl_simplified(lon, lat, h_center,
                                                                                        r_distance, i,
                                                                                        start_angle,
                                                                                        end_angle, dem)

    # R3
    result_r3 = analysis_by_r3(lon, lat, h_center, r_distance, i,
                               start_angle,
                               end_angle, dem)
    r3_count = np.count_nonzero(result_r3 == 1)
    visibility_rate = r3_count / actually_count
    # 计算不同值的数量
    diff_count_pderl_optrl = np.count_nonzero(result_spderl != result_opt_rl_spderl)
    diff_count_pderl_r3 = np.count_nonzero(result_spderl != result_r3)
    diff_count_optrl_r3 = np.count_nonzero(result_opt_rl_spderl != result_r3)
    # 错误率
    actually_error_rate_pderl = diff_count_pderl_r3 / actually_count
    actually_error_rate_optrl = diff_count_optrl_r3 / actually_count
    diff_accuacy = actually_error_rate_pderl - actually_error_rate_optrl
    print(
        f"actually_count:{actually_count}\ndiff_count_pderl_r3:{diff_count_pderl_r3}\ndiff_count_optrl_r3:{diff_count_optrl_r3}\ndiff_count_pderl_optrl :{diff_count_pderl_optrl}\nactually_error_rate_pderl:{actually_error_rate_pderl}\nactually_error_rate_optrl:{actually_error_rate_optrl}\nvisibility_rate:{visibility_rate}\n")
    data_line.extend(
        [dem.file_path, lon, lat, start_angle, end_angle, i, actually_error_rate_pderl, actually_error_rate_optrl, actually_count, diff_count_pderl_r3, diff_count_optrl_r3, visibility_rate, diff_accuacy])
    data.append(data_line)


analysis_accuracy()
print(f"精度测试运行结束时间:{datetime.now()}")
# 写入CSV文件
with open('output/参考线优化0127_9.csv', 'w', newline='') as csvfile:
    # 列名
    headers = ['DEM数据文件', '经度', '纬度', "起始角度", "终止角度", "视点高度", "SPDERL算法整体错误率",
               "参考线改进后算法整体错误率", "总个数", "pderl误差数量", "optrl误差数量", "可视率", "精度增量"]
    writer = csv.writer(csvfile)
    writer.writerow(headers)  # 写入列名
    writer.writerows(data)  # 写入数据行
print("CSV文件生成完毕！")
