import csv
import math
from datetime import datetime
import random

import numpy as np

import dem_data
from main import analysis_by_dpderl_simplified, analysis_by_pderl_with_angle, analysis_by_xdraw, \
    analysis_by_xpderl_with_angle

data = []


def analysis_aggregation_of_error_points():
    print(f"聚合点误差测试运行起始时间:{datetime.now()}")
    for h in range(3, 76, 3):
        # dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif", 30)
        # analysis_get_neighbor_error_with_pderl(0.1, h, 100, dem)
        # dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif", 30)
        # analysis_get_neighbor_error_with_pderl(0.1, h, 100, dem)
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif", 30)
        analysis_get_neighbor_error_with_pderl(0.1, h, 100, dem)


def analysis_get_neighbor_error_with_pderl(r, h, max_count, dem):
    global data

    data_line = []
    if max_count < 1:
        return

    count = 0
    while count <= max_count:
        data_line = []
        radom = random.random()
        dlon = dem.dx
        dlat = dem.dy
        start_lon = dem.start_x + r + dlon + 0.00181
        start_lat = dem.start_y + r + dlat + 0.00181
        max_lon = dem.max_lon - r - dlon - 0.00181
        max_lat = dem.max_lat - r - dlat - 0.00181
        start_angle = random.uniform(0, 360)  # 生成起始角度，范围在 0 到 360 之间
        end_angle = random.uniform(start_angle, 360)  # 生成终止角度，确保大于起始角度且在 0 到 360 之间

        lon = start_lon + radom * (max_lon - start_lon)
        lat = start_lat + radom * (max_lat - start_lat)
        h_center = dem.height[int((lon - dem.start_x) / dem.dx)][int((lat - dem.start_y) / dem.dy)]
        print(f"result_our_{h}_{count}:\n{lon} {lat} {h_center} {r * dem.rdx} {h} {start_angle} {end_angle}")
        _, _, result_our = analysis_by_dpderl_simplified(lon, lat, h_center, r * dem.rdx, h, start_angle, end_angle,
                                                         dem)
        result_pderl = analysis_by_pderl_with_angle(lon, lat, h_center, r * dem.rdx, h, start_angle, end_angle, dem)
        result_xdraw = analysis_by_xdraw(lon, lat, h_center, r * dem.rdx, h, start_angle, end_angle, dem)
        result_xpderl = analysis_by_xpderl_with_angle(lon, lat, h_center, r * dem.rdx, h, start_angle, end_angle, dem)

        neighbor_our = get_neighbor_err(result_our, result_pderl, 10, 0.9)
        neighbor_xdraw = get_neighbor_err(result_xdraw, result_pderl, 10, 0.9)
        neighbor_xpderl = get_neighbor_err(result_xpderl, result_pderl, 10, 0.9)
        sum_point = result_pderl.size
        count += 1
        data_line.extend(
            [dem.file_path, lon, lat, start_angle, end_angle, h, *(neighbor_our[1:11] / sum_point),
             *(neighbor_xdraw[1:11] / sum_point),
             *(neighbor_xpderl[1:11] / sum_point)])
        print(data_line)

        data.append(data_line)


def get_neighbor_err(result_test, result_pderl, max_neighbor, rate):
    neighbor_err_count = np.zeros(15)
    lon_count, lat_count = result_pderl.shape
    i = 0
    while i < lon_count:
        j = 0
        while j < lat_count:
            if result_pderl[i, j] != result_test[i, j]:
                count = 1
                tmp_neighbor = 1
                while tmp_neighbor <= max_neighbor:
                    lons = i - tmp_neighbor
                    lone = i + tmp_neighbor
                    lats = j - tmp_neighbor
                    late = j + tmp_neighbor
                    if lons < 0 or lats < 0 or lone >= lon_count or late >= lat_count:
                        tmp_neighbor += 1
                        continue
                    tmp_lat = lats
                    while tmp_lat <= late:
                        if result_pderl[lons, tmp_lat] != result_test[lons, tmp_lat]:
                            count += 1
                        if result_pderl[lone, tmp_lat] != result_test[lone, tmp_lat]:
                            count += 1
                        tmp_lat += 1
                    tmp_lon = lons + 1
                    while tmp_lon < lone:
                        if result_pderl[tmp_lon, lats] != result_test[tmp_lon, lats]:
                            count += 1
                        if result_pderl[tmp_lon, late] != result_test[tmp_lon, late]:
                            count += 1
                        tmp_lon += 1
                    if count >= math.ceil(rate * (tmp_neighbor * 2 + 1) * (tmp_neighbor * 2 + 1)):
                        neighbor_err_count[tmp_neighbor] += 1
                    tmp_neighbor += 1
            j += 1
        i += 1
    return neighbor_err_count


analysis_aggregation_of_error_points()
print(f"误差点聚合运行结束时间:{datetime.now()}")
# 写入CSV文件
with open('output/误差点聚合0315_丘陵.csv', 'w', newline='') as csvfile:
    # 列名
    headers = ['DEM数据文件', '经度', '纬度', "起始角度", "终止角度", "视点高度",
               "1邻域（XPDERL）", "2邻域（XPDERL）", "3邻域（XPDERL）", "4邻域（XPDERL）", "5邻域（XPDERL）", "6邻域（XPDERL）",
               "7邻域（XPDERL）", "X8邻域（XPDERL）", "9邻域（XPDERL）", "10邻域（XPDERL）", "1邻域（Our）", "2邻域（Our）",
               "3邻域（Our）", "4邻域（Our）", "5邻域（Our）", "6邻域（Our）",
               "7邻域（Our）", "X8邻域（Our）", "9邻域（Our）", "10邻域（Our）", "1邻域（XDraw）", "2邻域（XDraw）", "3邻域（XDraw）",
               "4邻域（XDraw）", "5邻域（XDraw）", "6邻域（XDraw）",
               "7邻域（XDraw）", "8邻域（XDraw）", "9邻域（XDraw）", "10邻域（XDraw）"]
    writer = csv.writer(csvfile)
    writer.writerow(headers)  # 写入列名
    writer.writerows(data)  # 写入数据行
print("CSV文件生成完毕！")
