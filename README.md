# 项目名称
本研究提出一种分区策略实现任意视场角所形成可视区域均可采用DPDERL算法进行可视性分析且解决分区策略所引起边缘点、搜索冗余等问题，并通过优化构建初始参考线避免观察点位于网格点时出现大量误差，最终将视频监控物理特性视场角与PDERL算法通过分区策略结合实现视频监控的快速近似视域分析。本研究算法为视频监控等需分析范围根据物理特性动态变化的设备，在视域分析领域中实现了速度与精度之间的有效平衡，同时对于需要快速响应且对精度有较高要求的视域分析场景有重要的参考价值。
###### **安装步骤**
```sh
git clone https://github.com/shaojintian/Best_README_template.git
```
## 文件目录说明

```
filetree
├── /data/
│  ├── /output/
│  │  ├── grid
│  │  └── image
│  │  └── shape
├── aggregation_of_error_points.py
├── comparison_accuracy.py
├── comparison_algorithm.py
├── comparison_optrl_with_peder.py
├── comparison_speed.py
├── dem_data.py
├── linked_line.py
├── main.py
├── partition_algorithm.py
├── partition_optimal_reference_lineal_gorithm.py
```
## 代码简介

### 主程序
在main.py下包含了本研究算法DPDERL、XPDERL、PDERL、XDraw、R3算法主体框架及对输入的参数数据预处理方法。
### DPDERL算法具体实现 
在partition_optimal_reference_lineal_gorithm.py中包含一个DPDERL算法的实现类，其中定义了算法所需要的参数(视场角、视距、视域半径等)及针对不同视场角所采取不同分区策略对应算法的具体实现。
### 速度测试（论文原标题：speed）
运行comparison_speed.py即可测试随机参数DPDERL算法与XDraw算法在山区、丘陵、平原地形下的耗时，并生成相应的csv文件。
### 精度对比测试（论文原标题：Accuracy Comparison Experiment）
运行comparison_optrl_with_peder.py即可测试当观察点位于网格点时DPDERL算法与PDERL算法在山区、丘陵、平原地形下的误差值，并生成相应的csv文件。
### 整体精度测试（论文原标题：Overall Accuracy）
运行comparison_accuracy.py即可测试以R3算法为基准DPDERL算法的整体误差率。
### 误差点聚合（论文原标题：Aggregation of error points）
运行aggregation_of_error_points.py即可测试以PDERL算法为基准DPDERL算法、XPDERL算法、XDraw算法的聚合点误差程度。
### 数据
在/data/文件目录下包含了本研究所用的所有DSM数据
### 输出文件
/output/grid目录下存储了缺陷分析所用的示例生成的栅格数据
/output/shape目录下存储了缺陷分析所用的示例生成的shapefile数据
/output/下包含了各实验对应测试csv文件

## 版本控制

该项目使用Git进行版本管理。您可以在repository参看当前可用版本。
