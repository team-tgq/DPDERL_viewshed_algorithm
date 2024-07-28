import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams['font.sans-serif'] = ['SimHei']  # 替换为您安装的中文字体名称
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示异常的问题
# 读取CSV文件
df = pd.read_csv('output/返稿_速度.csv', encoding='gbk')

# 按照地形和高度分组，并计算SPDERL和xdraw算法耗时的平均值
grouped_df = df.groupby(['DEM数据文件', '视点高度']).agg({
    'SPDERL算法耗时(秒)': 'mean',
    'xdraw算法耗时(秒)': 'mean'
}).reset_index()

# 打印结果
print(grouped_df)

# 如果需要将结果保存为新的CSV文件
grouped_df.to_csv('output/refurbishments/effective_execution_times.csv', index=False)
