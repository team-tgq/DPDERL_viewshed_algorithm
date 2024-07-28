import csv
import math
import random
from datetime import datetime

import numpy as np

import dem_data
from main import analysis_by_dpderl_simplified, analysis_by_r3
from output_veiwshed_grid import generate_viewable_raster
from output_viewshed_grid_nochange import generate_viewable_raster_nochange

data = []
terrain_type_dict = {
    "data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif": "山区",
    "data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif": "平原",
    "data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif": "丘陵"
}


def analysis_accuracy():
    # # 起始高度为1 accuracy
    # dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif", 30)
    # # 生成一个0-1之间的浮点数
    # step = random.random() * 0.6
    # analysis_accuracy_repeat(dem, 097.1, 28.3 + step, 097.9, 28.1 + step, 1, 5)
    # dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif", 30)
    # step = random.random() * 0.6
    # analysis_accuracy_repeat(dem, 114.1, 34.3 + step, 114.9, 34.1 + step, 1, 5)
    #
    # dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif", 30)
    # step = random.random() * 0.6
    # analysis_accuracy_repeat(dem, 119.1, 41.3 + step, 119.9, 41.1 + step, 1, 5)
    print(f"精度测试运行起始时间:{datetime.now()}")
    start_height = 1
    end_height = 5000

    for height in range(start_height, end_height + 1, 50):
        # 起始高度为1 accuracy
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif", 30)
        # 生成一个0-1之间的浮点数
        step = random.random() * 0.6
        analysis_accuracy_repeat(dem, 097.1, 28.3 + step, 097.9, 28.1 + step, height, 50)
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif", 30)
        step = random.random() * 0.6
        analysis_accuracy_repeat(dem, 114.1, 34.3 + step, 114.9, 34.1 + step, height, 50)
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif", 30)
        step = random.random() * 0.6
        analysis_accuracy_repeat(dem, 119.1, 41.3 + step, 119.9, 41.1 + step, height, 50)


"""
dem: 选择的地形文件
r: 实际半径,因此需转化为格网半径
h: 观察点距离海平面高度
limit: 同一高度测多少次
"""


def analysis_accuracy_repeat(dem, min_lon, max_lat, max_lon, min_lat, h, limit):
    global data
    count = 0
    while count < limit:
        count += 1
        data_line = []
        lon_random_1 = random.random() * (max_lon - min_lon) + min_lon
        lon_random_2 = random.random() * (max_lon - min_lon) + min_lon
        lon_r = abs(lon_random_1 - lon_random_2) / 2
        if lon_r < 0.0027:
            count -= 1
            continue
        lon = (min(lon_random_1, lon_random_2) + max(lon_random_1, lon_random_2)) / 2

        lat_random_1 = random.random() * (max_lat - min_lat) + min_lat
        lat_random_2 = random.random() * (max_lat - min_lat) + min_lat
        lat_r = abs(lat_random_1 - lat_random_2) / 2
        if lat_r < 0.0027:
            count -= 1
            continue

        lat = (min(lat_random_1, lat_random_2) + max(lat_random_1, lat_random_2)) / 2
        to_lon = 0
        to_lat = 0
        if lat_r > lon_r:
            to_lon = lon + lon_r
            to_lat = lat
        else:
            to_lon = lon
            to_lat = lat + lat_r

        # 计算半径
        r = math.sqrt((lon - to_lon) ** 2 + (lat - to_lat) ** 2)
        r_distance = r * dem.rdx

        start_angle = random.uniform(0, 360)  # 生成起始角度，范围在 0 到 360 之间
        end_angle = random.uniform(start_angle, 360)  # 生成终止角度，确保大于起始角度且在 0 到 360 之间
        h_center = dem.height[int((lon - dem.start_x) / dem.dx)][int((lat - dem.start_y) / dem.dy)]
        print(
            f"第{count}次测试:\ndem.file_path:{dem.file_path}\nlon:{lon}\nlat:{lat}\nh_center:{h_center}\nh:{h}\nstart_angle:{start_angle}\nend_angle:{end_angle}\nr:{r}")
        # SPDERL
        result_spderl, actually_count,_ = analysis_by_dpderl_simplified(lon, lat, h_center, r_distance, h,
                                                                        start_angle,
                                                                        end_angle, dem)
        file_name = f"Spderl_{terrain_type_dict[dem.file_path]}_{lon}_{lat}.tif"
        x_grid_observe, y_grid_observe, x_grid_center, y_grid_center, a, b, _ = dem.get_point_location(lon, lat)

        # generate_viewable_raster(result_spderl, lon, lat, r_distance, h,
        #                          start_angle,
        #                          end_angle, dem.file_path, file_name, x_grid_center, y_grid_center, x_grid_observe, y_grid_observe)

        # R3
        file_name = f"R3_{terrain_type_dict[dem.file_path]}_{lon}_{lat}.tif"
        result_r3 = analysis_by_r3(lon, lat, h_center, r_distance, h,
                                   start_angle,
                                   end_angle, dem)

        # generate_viewable_raster(result_r3, lon, lat,  r_distance, h,
        #                          start_angle,
        #                          end_angle, dem.file_path, file_name, x_grid_center, y_grid_center, x_grid_observe, y_grid_observe)
        # 结果数组总数量 并不等于实际分析的数量
        sum_point = result_spderl.size
        # 计算不同值的数量
        diff_count = np.count_nonzero(result_spderl != result_r3)

        # 错误率
        actually_error_rate = diff_count / actually_count

        sum_error_rate = diff_count / sum_point

        print(
            f"diff_count:{diff_count}\nsum_point:{sum_point}\nactually_count:{actually_count}\nactually_error_rate:{actually_error_rate}\n")
        data_line.extend(
            [dem.file_path, lon, lat, start_angle, end_angle, h, actually_error_rate])
        data.append(data_line)


analysis_accuracy()
print(f"精度测试运行结束时间:{datetime.now()}")
# 写入CSV文件
with open('output/精度测试_0109.csv', 'w', newline='') as csvfile:
    # 列名
    headers = ['DEM数据文件', '经度', '纬度', "起始角度", "终止角度", "视点高度", "SPDERL算法整体错误率"]
    writer = csv.writer(csvfile)
    writer.writerow(headers)  # 写入列名
    writer.writerows(data)  # 写入数据行
print("CSV文件生成完毕！")
