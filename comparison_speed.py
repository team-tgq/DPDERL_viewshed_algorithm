import csv
import random
from datetime import datetime
import dem_data
from main import analysis_by_dpderl_simplified, analysis_by_r3, analysis_by_xdraw

data = []


def analysis_speed():
    print(f"速度测试运行起始时间:{datetime.now()}")
    start_height = 1
    end_height = 5000
    for height in range(start_height, end_height + 1, 10):
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif", 30)
        analysis_speed_repeat(dem, 10000, height, 50)
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif", 30)
        analysis_speed_repeat(dem, 10000, height, 50)
        dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif", 30)
        analysis_speed_repeat(dem, 10000, height, 50)


"""
dem: 选择的地形文件
r: 实际半径,因此需转化为格网半径
h: 观察点距离海平面高度
limit: 同一高度测多少次
"""


def analysis_speed_repeat(dem, r, h, limit):
    global data
    r_grid = r / dem.rdx
    start_lon = dem.start_x + r_grid + dem.dx
    start_lat = dem.start_y + r_grid + dem.dy
    max_lon = dem.max_lon - r_grid - dem.dx
    max_lat = dem.max_lat - r_grid - dem.dy
    count = 0
    while count < limit:
        data_line = []
        lon = random.uniform(start_lon, max_lon)
        lat = random.uniform(start_lat, max_lat)
        start_angle = random.uniform(0, 360)  # 生成起始角度，范围在 0 到 360 之间
        end_angle = random.uniform(start_angle, 360)  # 生成终止角度，确保大于起始角度且在 0 到 360 之间
        date = datetime.now()
        h_center = dem.height[int((lon - dem.start_x) / dem.dx)][int((lat - dem.start_y) / dem.dy)]
        print(
            f"第{count + 1}次测试:\ndem.file_path:{dem.file_path}\nlon:{lon}\nlat:{lat}\nh_center:{h_center}\nh:{h}\nstart_angle:{start_angle}\nend_angle:{end_angle}\n")
        # 获取观察点底部高程值
        # SPDERL
        analysis_by_dpderl_simplified(lon, lat, h_center, r, h,
                                      start_angle,
                                      end_angle, dem)
        spderl_time = datetime.now() - date
        date = datetime.now()

        # xdraw
        analysis_by_xdraw(lon, lat, h_center, r, h,
                          start_angle,
                          end_angle, dem)
        xdraw_time = datetime.now() - date
        #  spderl_time.total_seconds() 时间间隔总秒数
        print(f"spderl_time:{spderl_time.total_seconds()}\nxdraw_time:{xdraw_time.total_seconds()}\n")
        data_line.extend(
            [dem.file_path, lon, lat, start_angle, end_angle, h, spderl_time.total_seconds(),
             xdraw_time.total_seconds(),
             spderl_time.total_seconds() / xdraw_time.total_seconds()])
        data.append(data_line)
        count += 1


analysis_speed()
print(f"速度测试运行结束时间:{datetime.now()}")
# 写入CSV文件
with open('output/速度测试_0109.csv', 'w', newline='') as csvfile:
    # 列名
    headers = ['DEM数据文件', '经度', '纬度', "起始角度", "终止角度", "视点高度", "SPDERL算法耗时(秒)",
               "xdraw算法耗时(秒)", "时间消耗比"]
    writer = csv.writer(csvfile)
    writer.writerow(headers)  # 写入列名
    writer.writerows(data)  # 写入数据行

print("CSV文件生成完毕！")
