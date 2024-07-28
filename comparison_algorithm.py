from datetime import datetime

import numpy as np

import dem_data
from main import analysis_by_r3, analysis_by_dpderl_simplified, get_initial_param_spderl, \
    get_start_loc, judge_is_out_bound
import random

from output_viewshed_grid_nochange import generate_viewable_raster_nochange

high_error_array = []
dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif", 30)

terrain_type_dict = {
    "data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif": "山区",
    "data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif": "平原",
    "data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif": "丘陵"
}

def comparison_with_r3(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle, horizontal_end_angle,
                       cnt, dem=dem):
    global high_error_array
    # 获取当前时间
    time = datetime.now()

    print(
        f"第{cnt}次测试:\nhorizontal_start_angle:{horizontal_start_angle}\nhorizontal_end_angle:{horizontal_end_angle}\n")

    result_spderl, count, result_opt_rl_spderl = analysis_by_dpderl_simplified(x_center, y_center, h_center, radius,
                                                                               h_stand,
                                                                               horizontal_start_angle,
                                                                               horizontal_end_angle, dem)
    spderl_time = datetime.now() - time

    time = datetime.now()

    result_r3 = analysis_by_r3(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle,
                               horizontal_end_angle, dem)
    r3_time = datetime.now() - time
    # 结果数组总数量 并不等于实际分析的数量
    sum_point = result_spderl.size
    # 计算不同值的数量
    diff_count = np.count_nonzero(result_spderl != result_r3)

    diff_array = np.where(result_spderl != result_r3, 1, 0)
    r3_count = np.count_nonzero(result_r3 == 1)
    visibility_rate = r3_count / count
    # 错误率
    diff_count_pderl_optrl = np.count_nonzero(result_spderl != result_opt_rl_spderl)
    diff_count_pderl_r3 = np.count_nonzero(result_spderl != result_r3)
    diff_count_optrl_r3 = np.count_nonzero(result_opt_rl_spderl != result_r3)
    # 错误率
    actually_error_rate_pderl = diff_count_pderl_r3 / count
    actually_error_rate_optrl = diff_count_optrl_r3 / count
    print(
        f"actually_count:{count}\ndiff_count_pderl_optrl :{diff_count_pderl_optrl}\ndiff_count_pderl_r3:{diff_count_pderl_r3}\ndiff_count_optrl_r3:{diff_count_optrl_r3}\nactually_error_rate_pderl:{actually_error_rate_pderl}\nactually_error_rate_optrl:{actually_error_rate_optrl}\nvisibility_rate:{visibility_rate}\n")
    return diff_count, count


# 山区
# dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N28_00_E097_00_DEM.tif", 30)

# 平原
dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N34_00_E114_00_DEM.tif", 30)

# 丘陵
# dem = dem_data.Dem("data/Copernicus_DSM_COG_10_N41_00_E119_00_DEM.tif", 30)

# 夹角越小 错误越多

x_center = 114.64375000000001
y_center = 34.52513888888889
h_center = dem.height[int((x_center - dem.start_x) / dem.dx)][int((y_center - dem.start_y) / dem.dy)]
# to_x = 97.421364
# to_y = 28.706327
h_stand = 110
radius = 2000

horizontal_start_angle = 128.91810081731498
horizontal_end_angle = 326.4852929095307


def visualization():
    result_spderl, count, result_opt_rl_spderl = analysis_by_dpderl_simplified(x_center, y_center, h_center, radius,
                                                                               h_stand,
                                                                               horizontal_start_angle,
                                                                               horizontal_end_angle, dem)

    result_spderl = np.rot90(result_spderl)
    start_lon, start_lat = get_start_loc(x_center, y_center, radius)
    file_name = f"spderl_{terrain_type_dict[dem.file_path]}_{x_center}_{y_center}.tif"
    generate_viewable_raster_nochange(result_spderl, x_center, y_center, radius, h_stand, horizontal_start_angle,
                                      horizontal_end_angle, dem.file_path, file_name, start_lon, start_lat)

    result_opt_rl_spderl = np.rot90(result_opt_rl_spderl)
    file_name = f"opt_rl_spderl_{terrain_type_dict[dem.file_path]}_{x_center}_{y_center}.tif"
    generate_viewable_raster_nochange(result_opt_rl_spderl, x_center, y_center, radius, h_stand, horizontal_start_angle,
                                      horizontal_end_angle, dem.file_path, file_name, start_lon, start_lat)

    result_r3 = analysis_by_r3(x_center, y_center, h_center, radius, h_stand, horizontal_start_angle,
                               horizontal_end_angle, dem)
    result_r3 = np.rot90(result_r3)
    file_name = f"r3_{terrain_type_dict[dem.file_path]}_{x_center}_{y_center}.tif"
    generate_viewable_raster_nochange(result_r3, x_center, y_center, radius, h_stand, horizontal_start_angle,
                                      horizontal_end_angle, dem.file_path, file_name, start_lon, start_lat)


def repeat_test():
    all_count = 0
    all_diff_count = 0
    for i in range(10, 3000 + 1, 10):
        start_angle = random.uniform(0, 360)  # 生成起始角度，范围在 0 到 360 之间
        end_angle = random.uniform(start_angle, 360)  # 生成终止角度，确保大于起始角度且在 0 到 360 之间
        is_out_bound = True
        r_distance = 2000
        lon = 0
        lat = 0
        while is_out_bound:
            random_integer_x = random.randint(1, 3600)
            random_integer_y = random.randint(1, 3600)
            lon = dem.start_x + random_integer_x * dem.dx
            lat = dem.start_y + random_integer_y * dem.dy
            is_out_bound = judge_is_out_bound(lon, lat, r_distance, dem)
        x_center = dem.start_x + 561 * dem.dx
        y_center = dem.start_y + 550 * dem.dy
        h_center = dem.height[int((x_center - dem.start_x) / dem.dx)][int((y_center - dem.start_y) / dem.dy)]
        diff_count, count = comparison_with_r3(x_center, y_center, h_center, radius, h_stand, start_angle, end_angle,
                                               i + 1, dem)
        all_count += count
        all_diff_count += diff_count
    all_error_rate = all_diff_count / all_count
    print(f"all_error_rate:{all_error_rate}")
    print(high_error_array)


comparison_with_r3(x_center, y_center, h_center, 2000, h_stand, horizontal_start_angle,
                   horizontal_end_angle, 1,
                   dem)

# repeat_test()
# visualization()
