# Language

- [English](#english)
- [中文](#中文)

---

## English
## Project Introduction
This study proposes a partitioning strategy that enables visibility analysis using the DPDERL algorithm for any field of view, addressing issues such as edge points and search redundancies caused by the partitioning strategy. By optimizing the construction of initial reference lines, the method avoids significant errors when observation points are located at grid points. Ultimately, the study integrates the physical characteristics of video surveillance field of view with the PDERL algorithm through a partitioning strategy to achieve rapid approximate visibility analysis for video surveillance. The algorithm developed in this research provides a dynamic balance between speed and accuracy in visibility analysis for devices whose operational scope changes based on physical characteristics, such as video surveillance systems. It holds significant referential value for visibility analysis scenarios that require rapid response and high precision.

### main program
Under main.py contains the main framework of the algorithms DPDERL, XPDERL, PDERL, XDraw, R3 algorithms of this study and the method of preprocessing the input parameter data.
### Implementation of the DPDERL Algorithm
In the file partition_optimal_reference_lineal_gorithm.py, an implementation class for the DPDERL (Dynamic Range proximity-direction-elevation Reference Line) algorithm is included. This class is designed to encapsulate the essential parameters required by the algorithm such as field of view angle, visibility radius, and observer height. It also details the specific implementations of the algorithm that vary according to different field of view angles, employing distinct partitioning strategies accordingly.
### Speed Test
To test the execution time of the DPDERL algorithm compared to the XDraw algorithm under varying terrain conditions such as mountainous, hilly, and plain areas, run comparison_speed.py. This script will evaluate the performance of the algorithms with random parameters and generate a corresponding CSV file detailing the time consumption for each scenario.
### Accuracy Comparison Experiment
To conduct an error evaluation when observation points are situated on grid points, run comparison_optrl_with_peder.py. This script facilitates the comparison of the DPDERL algorithm and the PDERL algorithm in terms of their error values across various terrains such as mountainous, hilly, and plain areas. It will also generate a corresponding CSV file that documents the discrepancies observed for each scenario.
### Overall Accuracy
To evaluate the overall error rate of the DPDERL algorithm with the R3 algorithm as the benchmark, execute comparison_accuracy.py. This script is designed to systematically measure and compare the accuracy discrepancies between the DPDERL algorithm and the established R3 algorithm.
### Aggregation of error points
To assess the aggregation of error points for the DPDERL, XPDERL, and XDraw algorithms using the PDERL algorithm as a benchmark, execute aggregation_of_error_points.py. This script is intended to measure and compare the extent of error aggregation where the calculations of these algorithms differ from those of the PDERL algorithm.
### data
The /data/ directory contains all the Digital Surface Model (DSM) data utilized in this study.
### output
The /output/grid directory contains the raster data generated from examples used for defect analysis. The /output/shape directory stores the shapefile data generated from examples used for defect analysis. Additionally, the /output/ directory includes the CSV files corresponding to tests conducted in various experiments.

---

## 中文

## 项目介绍
本研究提出一种分区策略实现任意视场角所形成可视区域均可采用DPDERL算法进行可视性分析且解决分区策略所引起边缘点、搜索冗余等问题，并通过优化构建初始参考线避免观察点位于网格点时出现大量误差，最终将视频监控物理特性视场角与PDERL算法通过分区策略结合实现视频监控的快速近似视域分析。本研究算法为视频监控等需分析范围根据物理特性动态变化的设备，在视域分析领域中实现了速度与精度之间的有效平衡，同时对于需要快速响应且对精度有较高要求的视域分析场景有重要的参考价值。

## 文件目录说明

```
filetree
├── /data/
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
在partition_optimal_reference_lineal_gorithm.py中包含一个DPDERL算法的实现类，其中定义了算法所需要的参数(视场角、视域半径、视点高度等)及针对不同视场角所采取不同分区策略对应算法的具体实现。
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

## 版本控制

该项目使用Git进行版本管理。您可以在repository参看当前可用版本。
