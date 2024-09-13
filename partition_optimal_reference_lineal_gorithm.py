import math

from linked_line import LinkedLinePDE


class PartitionOptimalRLAlgorithm:
    def __init__(self, horizontal_start_angle, horizontal_end_angle,
                 x_grid_observe, y_grid_observe, x_grid_count, y_grid_count, see_height, x_grid_center,
                 y_grid_center, min_x, max_x, min_y, max_y, x_distance, y_distance, dem,
                 is_observe_not_in_left_top_diagonal, is_observe_not_in_right_bottom_diagonal, result, quadrant,
                 is_in_grid_point):
        self.horizontal_start_angle = horizontal_start_angle
        self.horizontal_end_angle = horizontal_end_angle
        self.x_grid_observe = x_grid_observe
        self.y_grid_observe = y_grid_observe
        self.x_grid_count = x_grid_count
        self.y_grid_count = y_grid_count
        self.see_height = see_height
        self.x_grid_center = x_grid_center
        self.y_grid_center = y_grid_center
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y
        self.x_distance = x_distance
        self.y_distance = y_distance
        self.dem = dem
        self.is_observe_not_in_left_top_diagonal = is_observe_not_in_left_top_diagonal
        self.is_observe_not_in_right_bottom_diagonal = is_observe_not_in_right_bottom_diagonal
        self.result = result
        self.quadrant = quadrant
        self.is_in_grid_point = is_in_grid_point

    # Analysis area on the right --> Vertical grid lines (south to north and west to east)
    def analysis_right(self):
        count = 0
        # e->height*p for the first point
        start_base_e = 0.0
        # Start segment of the reference line
        start_base_line = None
        # The current point against the old baseline
        base_line = None
        # Part of the line that will be connected to the new reference line is temporarily stored in this variable when the scope of the line is uncertain
        current_new_line = None

        u = 1 / (self.dem.dy * self.dem.rdy)
        # Absolute equation direction based on the direction of analysis: here is an example where north is above south
        # Inverse of vertical grid spacing field distance use this to calculate segment coefficients a-->u = 1 / (dem.dy * dem.rdy)
        north_k, north_b = self.calculate_line_equation(self.horizontal_start_angle)
        south_k, south_b = self.calculate_line_equation(self.horizontal_end_angle)

        # No need to expand on a grid
        if south_k == 0:
            self.is_observe_not_in_right_bottom_diagonal = 0
        # Assign the regular latitude to find the other latitude from the equation
        # x_gird_center is a grid point to the left of the observation point if it is not on a grid point, and x_gird_center+1 is the first column to analyze (not considered on a grid point)
        x_grid_index = self.x_grid_center + 1

        # calculate start boundary search start and end index start round down and end round up
        start_y_bound = int(x_grid_index * south_k + south_b)
        start_y_index = start_y_bound - self.is_observe_not_in_right_bottom_diagonal
        end_y_bound = math.ceil(x_grid_index * north_k + north_b)
        end_y_index = end_y_bound + self.is_observe_not_in_left_top_diagonal

        # We need to think about the Y-axis
        if north_k == math.inf:
            end_y_index = self.max_y
        else:
            end_y_index = self.max_y if end_y_index > self.max_y else end_y_index

        start_y_index = self.min_y if start_y_index < self.min_y else start_y_index
        # The starting grid point is shifted down a few
        diff_y = self.y_grid_center - start_y_index

        start_distance_y_index = self.y_grid_count - diff_y - 1

        y_grid_index = start_y_index
        x_distance_index = self.x_grid_count
        y_distance_index = start_distance_y_index
        last_height = self.dem.height[x_grid_index][y_grid_index - 1] - self.see_height
        p = 1 / self.x_distance[x_distance_index]

        # The observation point is not a grid point and the algorithm just treats it as a grid point above or below the central grid point and does not evaluate it
        # Construct the right half initial reference line
        if self.is_in_grid_point is True:
            self.point_in_line_constx(end_y_index, y_grid_index, x_grid_index, x_distance_index)
            # while y_grid_index <= end_y_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     y_grid_index += 1
            #     y_distance_index += 1
            #     count += 1
            #  The target point continues to move to the right
            x_grid_index += 1
            x_distance_index += 1

            const_increase_y = 0 if north_k == math.inf else math.ceil(north_k)
            const_decrease_y = 0 if south_k == math.inf else math.floor(south_k)

            # The y-coordinate of the first observation point in the region
            start_distance_y_index += const_decrease_y

            # Expand the boundary up and down in the y direction
            start_y_index += const_decrease_y

            end_y_index += const_increase_y
            if start_distance_y_index < 0:
                adjust_y = abs(start_distance_y_index)
                start_distance_y_index = 0
                start_y_index += adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index
            p = 1 / self.x_distance[x_distance_index]
        while y_grid_index <= end_y_index and y_distance_index < 2 * self.y_grid_count - 1:
            # When the height of the target point is less than the height of the observation point, there will be a negative number
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # By recording e for the first segment, a for each segment, and d for the end of each segment, the entire terrain line e=ad+b can be expressed and recorded in the PDE equation space
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next
            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            y_grid_index += 1
            y_distance_index += 1
            count += 1

        # The target point keeps moving to the right
        x_grid_index += 1
        x_distance_index += 1
        # Since dy varies linearly and x is incremented uniformly each time, the increase in y is fixed but not by a multiple of the grid points so there may be redundant areas
        const_increase_y = 0 if north_k == math.inf else math.ceil(north_k)
        const_decrease_y = 0 if south_k == math.inf else math.floor(south_k)

        # The y-coordinate of the first observation point in the region
        start_distance_y_index += const_decrease_y

        # Expand the boundary up and down in the y direction
        start_y_index += const_decrease_y

        end_y_index += const_increase_y

        # On the right, the explicitness and implicitness of reference lines and points are calculated layer by layer
        while x_grid_index <= self.max_x and start_distance_y_index < 2 * self.y_grid_count - 1:
            # To prevent exceeding the bounds of the result array, the index is callback
            if start_distance_y_index < 0:
                adjust_y = abs(start_distance_y_index)
                start_distance_y_index = 0
                start_y_index += adjust_y

            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index

            # It's a little bit outside the analysis area
            last_height = self.dem.height[x_grid_index, y_grid_index - 1] - self.see_height

            p = 1 / self.x_distance[x_distance_index]
            base_line = start_base_line

            # Determine whether the starting point is explicit or implicit
            # Initialize d, a of the starting point
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            """
                Why do we need next: because the first point is aided by a point outside the scope of the analysis, if you go out and the first segment of the aided segment itself contains the next column d you take its own reference segment as the first segment, go out to the aided reference segment and copy e0
                Why do we want to compare K: because the reference line for the survey line is e0 for the reference line is not its first corresponding point so we need to find the reference line e for the corresponding point
                If the value of d at the end of a line segment below the reference line is less than the current value, find the corresponding value of △Ei
                k: direction d beseLint.EndK: direction d of the reference line The latter must cover the former
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di − di − 1) Here, only the E of the reference line is computed through the linked list iteration, and the survey line is not involved in the operation
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # The previous step amounts to sieving out the reference line segments in which the extended point participates
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
                This is the difference between the slope of the current starting point axis and the reference line
                a*(base_line.end_d-d) means that k is the k of the current target point, so the obtained value is the e of the corresponding reference line in the same direction D. 
                It is not necessarily equal to ei, it is the value in the range [ei,ei+1]
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # Visibility check
            if e_diff >= 0:
                # Update the guideline if it's not blocked because e is bigger
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 1 if y_grid_index >= self.y_grid_center else 4
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and x_distance_index >= self.x_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # Process from south to north
            y_grid_index += 1
            y_distance_index += 1
            count += 1
            # Determine whether the next point is explicit or implicit
            while y_grid_index <= self.max_y and y_grid_index <= end_y_index:
                # Form a new survey line segment
                count += 1
                d = self.y_distance[y_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # When the d and EndK of the reference line are less than the current k, the reference line iteration operation is needed to update the reference line and judge the visibility of the corresponding point
                while base_line and base_line.end_d < d:
                    """
                        Find the increment of the difference in the axis slope
                        The reason why a is fixed here is that it is accumulated for only one segment of the survey line and its slope is fixed, so a is unchanged. However, a segment of the survey line may contain multiple segments of the reference line, so the reference line will be updated constantly
                        △Ei=△Ei-1+(di-di-1)(ALi-ARi) Where E is the E between the reference line and the survey line
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # If the front and back visibility are different, then there is an intersection point
                    # it's visual
                    if e_diff >= 0:
                        # Become invisible
                        if e_diff_last < 0:
                            # d Not accurate at all.
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # The d value of the intersection
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # Insert and move the segment to the intersection point and connect the following reference lines
                            # If previously visible but now not visible, record the crossing point and insert the old reference line to the back
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # This also updates the guide line
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # The d-value of the intersection
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # Become invisible
                    if e_diff_last < 0:
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # Insert and move the segment to the intersection point and connect the following reference lines
                        # If previously visible but now not visible, record the crossing point and insert the old reference line to the back
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # Visibility check
                if e_diff >= 0:
                    # Update the guideline if it's not blocked because e is bigger
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 1 if y_grid_index >= self.y_grid_center else 4
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and x_distance_index >= self.x_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # Process from south to north
                y_grid_index += 1
                y_distance_index += 1

            # Outermost loop corresponds to grid point abscissa plus 1: west to east process
            x_grid_index += 1
            x_distance_index += 1
            start_y_index += const_decrease_y
            start_distance_y_index += const_decrease_y
            end_y_index += const_increase_y
        return self.result, count

    # on the underside --> Horizontal grid lines (north to south, west to east)
    def analysis_bottom(self):
        count = 0
        # 第一个点对应的e->height*p
        start_base_e = 0.0
        # 参考线的起始线段
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None
        u = 1 / (self.dem.dx * self.dem.rdx)
        east_k, east_b = self.calculate_line_equation(self.horizontal_start_angle)
        west_k, west_b = self.calculate_line_equation(self.horizontal_end_angle)

        # 若观察点不在网格点上则x_gird_center是其左边的一个网格点,x_gird_center+1为要分析的第一列（没考虑到在格网点上的情况）
        y_grid_index = self.y_grid_center

        # 计算起始边界搜索起点与终点索引 需考虑斜率为零的情况
        if west_k == 0:
            start_x_index = self.min_x
        else:
            start_x_index = int((y_grid_index - west_b) / west_k) - self.is_observe_not_in_left_top_diagonal
        if east_k == 0:
            end_x_index = self.max_x
        else:
            end_x_index = math.ceil(
                (y_grid_index - east_b) / east_k) + self.is_observe_not_in_right_bottom_diagonal

        # 起始终止角度位于y轴
        if west_k == math.inf:
            start_x_index = self.x_grid_center
        else:
            start_x_index = self.min_x if start_x_index < self.min_x else start_x_index

        if east_k == math.inf:
            end_x_index = self.x_grid_center
        else:
            end_x_index = self.max_x if end_x_index > self.max_x else end_x_index

        diff_x = self.x_grid_center - start_x_index

        start_distance_x_index = self.x_grid_count - diff_x - 1

        x_grid_index = start_x_index
        x_distance_index = start_distance_x_index
        y_distance_index = self.y_grid_count - 1
        last_height = self.dem.height[x_grid_index - 1, y_grid_index] - self.see_height
        p = -1 / self.y_distance[y_distance_index]

        if self.is_in_grid_point is True:
            self.point_in_line_consty(end_x_index, x_grid_index, y_grid_index, y_distance_index)
            # while x_grid_index <= end_x_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     x_grid_index += 1
            #     x_distance_index += 1
            #     count += 1
            y_grid_index -= 1
            y_distance_index -= 1

            const_increase_x = 0 if east_k == 0 or east_k == math.inf else math.ceil(-1 / east_k)
            const_decrease_x = 0 if west_k == 0 or west_k == math.inf else math.floor(-1 / west_k)

            start_x_index += const_decrease_x
            start_distance_x_index += const_decrease_x
            end_x_index += const_increase_x
            if start_distance_x_index < 0:
                adjust_x = abs(start_distance_x_index)
                start_distance_x_index = 0
                start_x_index += adjust_x
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index
            p = -1 / self.y_distance[y_distance_index]
        while x_grid_index <= end_x_index and x_distance_index < 2 * self.x_grid_count - 1:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            x_grid_index += 1
            x_distance_index += 1
            count += 1

        y_grid_index -= 1
        y_distance_index -= 1

        const_increase_x = 0 if east_k == 0 or east_k == math.inf else math.ceil(-1 / east_k)
        const_decrease_x = 0 if west_k == 0 or west_k == math.inf else math.floor(-1 / west_k)

        start_x_index += const_decrease_x
        start_distance_x_index += const_decrease_x
        end_x_index += const_increase_x

        # 下边逐层计算参考线与点的显隐性 从左往右分析如果起始x已经大于最大x则停止
        while y_grid_index >= self.min_y and start_distance_x_index < 2 * self.x_grid_count - 1:
            if start_distance_x_index < 0:
                adjust_x = abs(start_distance_x_index)
                start_distance_x_index = 0
                start_x_index += adjust_x
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index - 1, y_grid_index] - self.see_height

            p = -1 / self.y_distance[y_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e
            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 4 if x_grid_index >= self.x_grid_center else 3
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and y_distance_index <= self.y_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从西向东过程
            x_grid_index += 1
            x_distance_index += 1

            count += 1

            # 判断后续点的显隐性别
            while x_grid_index <= self.max_x and x_grid_index <= end_x_index:
                # 形成新的调查线段
                count += 1
                d = self.x_distance[x_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line and base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 4 if x_grid_index >= self.x_grid_center else 3
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and y_distance_index <= self.y_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从西向东过程
                x_grid_index += 1
                x_distance_index += 1

            # 最外层循环 对应网格点横坐标加1：从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

            start_x_index += const_decrease_x
            start_distance_x_index += const_decrease_x
            end_x_index += const_increase_x
        return self.result, count

    # On the left --> vertical grid lines (north to south, east to west)
    def analysis_left(self):
        # 第一个点对应的e->height*p
        start_base_e = 0.0
        # 参考线的起始线段
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        u = 1 / (self.dem.dy * self.dem.rdy)
        count = 0
        north_k, north_b = self.calculate_line_equation(self.horizontal_end_angle)
        south_k, south_b = self.calculate_line_equation(self.horizontal_start_angle)

        x_grid_index = self.x_grid_center

        # 计算起始边界搜索起点与终点索引
        end_y_bound = int(x_grid_index * south_k + south_b)
        end_y_index = end_y_bound - self.is_observe_not_in_left_top_diagonal
        start_y_bound = int(x_grid_index * north_k + north_b)
        start_y_index = start_y_bound + self.is_observe_not_in_right_bottom_diagonal + 1

        # 需要考虑在y轴上的情况 此时不能直接用方程去表示
        if north_k == math.inf:
            start_y_index = self.max_y
        else:
            start_y_index = self.max_y if start_y_index > self.max_y else start_y_index
        if south_k == math.inf:
            end_y_index = self.min_y
        else:
            end_y_index = self.min_y if end_y_index < self.min_y else end_y_index

        # 计算起始格网点与中心格网点的偏移量
        diff_y = self.y_grid_center - start_y_index

        start_distance_y_index = self.y_grid_count - diff_y - 1

        y_grid_index = start_y_index
        x_distance_index = self.x_grid_count - 1
        y_distance_index = start_distance_y_index
        last_height = self.dem.height[x_grid_index, y_grid_index + 1] - self.see_height

        p = -1 / self.x_distance[x_distance_index]

        const_increase_y = 0 if north_k == math.inf else math.ceil(-north_k)
        const_decrease_y = 0 if south_k == math.inf else math.floor(-south_k)
        if const_increase_y + start_y_index >= self.max_y:
            y_grid_index = self.max_y
            start_y_index = self.max_y
            start_distance_y_index = 2 * self.y_grid_count - 1
            y_distance_index = 2 * self.y_grid_count - 1
            const_increase_y = 0

        if self.is_in_grid_point is True:
            self.point_in_line_constx(end_y_index, y_grid_index, x_grid_index, x_distance_index)
            # while y_grid_index >= end_y_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     y_grid_index -= 1
            #     y_distance_index -= 1
            #     count += 1
            x_grid_index -= 1
            x_distance_index -= 1
            max_y_index = 2 * self.y_grid_count - 1
            start_distance_y_index += const_increase_y
            start_y_index += const_increase_y
            end_y_index += const_decrease_y
            if start_distance_y_index > max_y_index:
                adjust_y = abs(start_distance_y_index - max_y_index)
                start_distance_y_index = max_y_index
                start_y_index -= adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index
            p = -1 / self.x_distance[x_distance_index]
        while y_grid_index >= end_y_index and y_distance_index > 0:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)
            #  print(f"p:{p},y_distance_index:{y_distance_index},d:{d}")
            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            self.init_reference_line_by_r3(x_distance_index,y_distance_index)
            last_height = current_height

            y_grid_index -= 1
            y_distance_index -= 1
            count += 1

        x_grid_index -= 1
        x_distance_index -= 1

        start_distance_y_index += const_increase_y
        start_y_index += const_increase_y
        end_y_index += const_decrease_y
        # 左边逐层计算参考线与点的显隐性
        max_y_index = 2 * self.y_grid_count - 1
        count += 1
        while x_grid_index >= self.min_x and start_distance_y_index > 0:
            if start_distance_y_index > max_y_index:
                adjust_y = abs(start_distance_y_index - max_y_index)
                start_distance_y_index = max_y_index
                start_y_index -= adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index, y_grid_index + 1] - self.see_height

            p = -1 / self.x_distance[x_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 2 if y_grid_index >= self.y_grid_center else 3
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and x_distance_index <= self.x_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

            # 判断后续点的显隐性别
            while y_grid_index >= self.min_y and y_grid_index >= end_y_index:
                # 形成新的调查线段
                # print(f"p:{p},y_distance_index:{y_distance_index},d:{d}")
                count += 1
                d = -self.y_distance[y_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 2 if y_grid_index >= self.y_grid_center else 3
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and x_distance_index <= self.x_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从北向南过程
                y_grid_index -= 1
                y_distance_index -= 1

            # 最外层循环 对应网格点横坐标加1：从东向西过程
            x_grid_index -= 1
            x_distance_index -= 1

            end_y_index += const_decrease_y
            start_distance_y_index += const_increase_y
            start_y_index += const_increase_y
        return self.result, count

    # on top --> Horizontal grid lines (south to north, east to west)
    def analysis_top(self):
        # 第一个点对应的e->height*p
        start_base_e = 0.0
        # 参考线的起始线段
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        u = 1 / (self.dem.dx * self.dem.rdx)
        count = 0
        west_k, west_b = self.calculate_line_equation(self.horizontal_start_angle)
        east_k, east_b = self.calculate_line_equation(self.horizontal_end_angle)

        # 若观察点不在网格点上则x_gird_center是其左边的一个网格点,x_gird_center+1为要分析的第一列（没考虑到在格网点上的情况）
        y_grid_index = self.y_grid_center + 1

        # 计算起始边界搜索起点与终点索引 需考虑斜率为零的情况
        if east_k == 0:
            start_x_index = self.max_x
        else:
            start_x_index = math.ceil(
                (y_grid_index - east_b) / east_k) + self.is_observe_not_in_left_top_diagonal

        if west_k == 0:
            end_x_index = self.min_x
        else:
            end_x_index = int((y_grid_index - west_b) / west_k) - self.is_observe_not_in_right_bottom_diagonal

        # 起始终止角度位于y轴
        if west_k == math.inf:
            end_x_index = self.x_grid_center
        else:
            end_x_index = self.min_x if end_x_index < self.min_x else end_x_index

        if east_k == math.inf:
            start_x_index = self.x_grid_center
        else:
            start_x_index = self.max_x if start_x_index > self.max_x else start_x_index


        diff_x = self.x_grid_center - start_x_index

        start_distance_x_index = self.x_grid_count - diff_x - 1

        x_grid_index = start_x_index
        x_distance_index = start_distance_x_index
        y_distance_index = self.y_grid_count
        last_height = self.dem.height[x_grid_index + 1, y_grid_index] - self.see_height
        p = 1 / self.y_distance[y_distance_index]

        if self.is_in_grid_point is True > 0:
            self.point_in_line_consty(end_x_index, x_grid_index, y_grid_index, y_distance_index)
            # while x_grid_index >= end_x_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     x_grid_index -= 1
            #     x_distance_index -= 1
            #     count += 1
            y_grid_index += 1
            y_distance_index += 1

            # 直线线性变换,y均匀增加,x增加也为固定值.(目前未考虑斜率为0的情况)
            const_increase_x = 0 if east_k == 0 or east_k == math.inf else math.ceil(1 / east_k)
            const_decrease_x = 0 if west_k == 0 or west_k == math.inf else math.floor(1 / west_k)
            max_right_index = 2 * self.x_grid_count - 1
            # 拓展 x 方向
            start_x_index += const_increase_x
            start_distance_x_index += const_increase_x
            end_x_index += const_decrease_x
            if start_distance_x_index > max_right_index:
                adjust_x = abs(start_distance_x_index - max_right_index)
                start_distance_x_index = max_right_index
                start_x_index -= adjust_x
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index
            p = 1 / self.y_distance[y_distance_index]
        # 构造上半边初始参考线
        while x_grid_index >= end_x_index and x_distance_index > 0:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next
            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            x_grid_index -= 1
            x_distance_index -= 1
            count += 1

        y_grid_index += 1
        y_distance_index += 1

        # 直线线性变换,y均匀增加,x增加也为固定值.(目前未考虑斜率为0的情况)
        const_increase_x = 0 if east_k == 0 or east_k == math.inf else math.ceil(1 / east_k)
        const_decrease_x = 0 if west_k == 0 or west_k == math.inf else math.floor(1 / west_k)

        # 拓展 x 方向
        start_x_index += const_increase_x
        start_distance_x_index += const_increase_x
        end_x_index += const_decrease_x

        # 上边逐层计算参考线与点的显隐性
        # 防止超出x边界
        max_right_index = 2 * self.x_grid_count - 1
        while y_grid_index <= self.max_y and start_distance_x_index > 0:
            if start_distance_x_index > max_right_index:
                adjust_x = abs(start_distance_x_index - max_right_index)
                start_distance_x_index = max_right_index
                start_x_index -= adjust_x

            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index + 1, y_grid_index] - self.see_height

            p = 1 / self.y_distance[y_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 1 if x_grid_index >= self.x_grid_center else 2
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and y_distance_index >= self.y_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从东向西过程
            x_grid_index -= 1
            x_distance_index -= 1
            count += 1
            # 判断后续点的显隐性别
            while x_grid_index >= self.min_x and x_grid_index >= end_x_index:
                # 形成新的调查线段
                count += 1
                d = -self.x_distance[x_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 1 if x_grid_index >= self.x_grid_center else 2
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and y_distance_index >= self.y_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从东向西过程
                x_grid_index -= 1
                x_distance_index -= 1

            # 最外层循环 对应网格点横坐标加1：从南向北过程
            # 拓展 x 边界
            start_distance_x_index += const_increase_x
            start_x_index += const_increase_x
            end_x_index += const_decrease_x

            y_grid_index += 1
            y_distance_index += 1
        return self.result, count

    # is divided into two categories:
    # ① Parts larger than 180 are located to the left of the observation point --> vertical grid line
    # ② The part less than 180 is located to the right of the observation point --> the vertical grid line
    # Three line equations are used to divide the region into left region and right region: l1-> left region, l2-> middle region, l3-> right region
    def analysis_left_right(self):
        # 第一个点对应的e->height*p
        start_base_e = 0.0
        # 参考线的起始线段
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        # 左侧
        # 位于左侧 --> 纵向网格线(从北向南,从东向西)
        u = 1 / (self.dem.dy * self.dem.rdy)
        count = 0
        north_k, north_b = self.calculate_line_equation(self.horizontal_end_angle)

        # 首列，观察点不在格网点上则在左侧
        x_grid_index = self.x_grid_center

        # 计算起始边界搜索起点与终点索引
        end_y_index = self.min_y
        start_y_bound = int(x_grid_index * north_k + north_b)
        start_y_index = start_y_bound + self.is_observe_not_in_right_bottom_diagonal

        start_y_index = self.max_y if start_y_index > self.max_y else start_y_index
        end_y_index = self.min_y if end_y_index < self.min_y else end_y_index
        # 计算起始格网点向上移动了几个
        diff_y = self.y_grid_center - start_y_index

        start_distance_y_index = self.y_grid_count - diff_y - 1

        y_grid_index = start_y_index
        x_distance_index = self.x_grid_count - 1
        y_distance_index = start_distance_y_index
        last_height = self.dem.height[x_grid_index, y_grid_index + 1] - self.see_height

        p = -1 / self.x_distance[x_distance_index]

        # 若观察点位于格网点上
        if self.is_in_grid_point is True:
            self.point_in_line_constx(end_y_index, y_grid_index, x_grid_index, x_distance_index)
            # while y_grid_index >= end_y_index:
            #     # 第一列可见性通过R3算法计算
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     y_grid_index -= 1
            #     y_distance_index -= 1
            #     count += 1
            max_y_index = 2 * self.y_grid_count - 1
            # 向分析方向移动一格,起始y保持不变
            x_grid_index -= 1
            x_distance_index -= 1
            const_increase_y = 0 if north_k == math.inf else math.ceil(-north_k)
            start_distance_y_index += const_increase_y
            start_y_index += const_increase_y
            if start_distance_y_index > max_y_index:
                adjust_y = abs(start_distance_y_index - max_y_index)
                start_distance_y_index = max_y_index
                start_y_index -= adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index
            p = -1 / self.x_distance[x_distance_index]

        # 观察点不在格网点上时，通过R3算法利用首列构造左半边初始参考线
        while y_grid_index >= end_y_index and y_distance_index > 0:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            # 第一列可见性通过R3算法计算
            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            y_grid_index -= 1
            y_distance_index -= 1
            count += 1

        # 因为是dy线性变化每次x都是均匀递增,这样y增加的值也是固定的
        x_grid_index -= 1
        x_distance_index -= 1

        const_increase_y = 0 if north_k == math.inf else math.ceil(-north_k)

        start_distance_y_index += const_increase_y
        start_y_index += const_increase_y

        # 左边逐层计算参考线与点的显隐性
        max_y_index = 2 * self.y_grid_count - 1
        while x_grid_index >= self.min_x and start_distance_y_index > 0:
            if start_distance_y_index > max_y_index:
                adjust_y = abs(start_distance_y_index - max_y_index)
                start_distance_y_index = max_y_index
                start_y_index -= adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index, y_grid_index + 1] - self.see_height

            p = -1 / self.x_distance[x_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 2 if y_grid_index >= self.y_grid_center else 3
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and x_distance_index <= self.x_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

            # 判断后续点的显隐性别
            while y_grid_index >= self.min_y and y_grid_index >= end_y_index:
                # 形成新的调查线段
                d = -self.y_distance[y_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)
                count += 1
                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 2 if y_grid_index >= self.y_grid_center else 3
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and x_distance_index <= self.x_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从北向南过程
                y_grid_index -= 1
                y_distance_index -= 1

            # 最外层循环 对应网格点横坐标加1：从东向西过程
            start_distance_y_index += const_increase_y
            start_y_index += const_increase_y
            x_grid_index -= 1
            x_distance_index -= 1

        # 右侧
        # 分析区域位于右侧 --> 纵向网格线(从南到北 从西到东)
        # 根据分析的方向来绝对方程所在方向：以下为例，其中north 代表的方程在 south的上方
        # 参考线的起始线段
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        u = 1 / (self.dem.dx * self.dem.rdx)
        north_k, north_b = self.calculate_line_equation(self.horizontal_start_angle)

        # 先将规律变化的纬度进行赋值 从而根据方程求出另一纬
        # 若观察点不在网格点上则x_gird_center是其左边的一个网格点,x_gird_center+1为要分析的第一列（没考虑到在格网点上的情况）
        x_grid_index = self.x_grid_center + 1

        # 计算起始边界搜索起点与终点索引
        start_y_index = self.min_y
        end_y_bound = math.ceil(x_grid_index * north_k + north_b)
        end_y_index = end_y_bound + self.is_observe_not_in_left_top_diagonal
        if north_k == math.inf:
            end_y_index = self.max_y
        else:
            end_y_index = self.max_y if end_y_index > self.max_y else end_y_index

        start_y_index = self.min_y if start_y_index < self.min_y else start_y_index
        # 计算起始格网点向下移动了几个
        diff_y = self.y_grid_center - start_y_index

        start_distance_y_index = self.y_grid_count - diff_y - 1

        y_grid_index = start_y_index
        x_distance_index = self.x_grid_count
        y_distance_index = start_distance_y_index
        last_height = self.dem.height[x_grid_index, y_grid_index - 1] - self.see_height
        p = 1 / self.x_distance[x_distance_index]

        if self.is_in_grid_point is True:
            self.point_in_line_constx(end_y_index, y_grid_index, x_grid_index, x_distance_index)
            # 首列可见性由R3决定,从下一列开始构造参考线
            # while y_grid_index <= end_y_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     y_grid_index += 1
            #     y_distance_index += 1
            #     count += 1

            # 还原y值并扩展

            const_increase_y = 0 if north_k == math.inf else math.ceil(north_k)

            # 目标点继续向右移动
            x_grid_index += 1
            x_distance_index += 1
            end_y_index += const_increase_y

            if start_distance_y_index < 0:
                adjust_y = abs(start_distance_y_index)
                start_distance_y_index = 0
                start_y_index += adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index
            p = 1 / self.x_distance[x_distance_index]
        # 观察点不是网格点该算法也只是将其视为中心网格点上方或者下方的网格点 并没有对其求值
        # 构造右半边初始参考线
        while y_grid_index <= end_y_index and y_distance_index < 2 * self.y_grid_count - 1:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            # 首列可见性由R3决定
            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            y_grid_index += 1
            y_distance_index += 1
            count += 1

        # 因为是dy线性变化每次x都是均匀递增,这样y增加的值也是固定的
        const_increase_y = 0 if north_k == math.inf else math.ceil(north_k)

        # 目标点继续向右移动
        x_grid_index += 1
        x_distance_index += 1

        end_y_index += const_increase_y

        # 右边逐层计算参考线与点的显隐性
        while x_grid_index <= self.max_x and start_distance_y_index < 2 * self.y_grid_count - 1:
            # 防止超出结果数组边界,回调索引
            if start_distance_y_index < 0:
                adjust_y = abs(start_distance_y_index)
                start_distance_y_index = 0
                start_y_index += adjust_y

            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index, y_grid_index - 1] - self.see_height

            p = 1 / self.x_distance[x_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 1 if y_grid_index >= self.y_grid_center else 4
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                # 判断当前象限的点不在范围内及不处于当前象限点
                if is_in_range and x_distance_index >= self.x_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从南向北过程
            y_grid_index += 1
            y_distance_index += 1

            # 判断后续点的显隐性别
            while y_grid_index <= self.max_y and y_grid_index <= end_y_index:
                # 形成新的调查线段
                count += 1
                d = self.y_distance[y_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                # print(f"{base_line.end_d},{d},{y_grid_index},{x_grid_index}")
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 1 if y_grid_index >= self.y_grid_center else 4
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    # 判断当前象限的点不在范围内及不处于当前象限点
                    if is_in_range and x_distance_index >= self.x_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从南向北过程
                y_grid_index += 1
                y_distance_index += 1

            # 最外层循环 对应网格点横坐标加1：从西向东过程
            x_grid_index += 1
            x_distance_index += 1
            end_y_index += const_increase_y
        return self.result, count

    # Divided into two categories:
    # ① The part greater than 180 is located below the observation point --> lateral grid line
    # ② The part less than 180 is located to the right of the observation point --> the vertical grid line
    # tune the intermediate index
    # Divide the region into lower and right regions by three line equations: l1-> lower left region, l2-> middle region, l3-> right region
    def analysis_right_bottom(self):
        # 第一个点对应的e->height*p
        start_base_e = 0.0
        # 参考线的起始线段
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        # if horizontal_end_angle != 180:
        # 下侧
        u = 1 / (self.dem.dx * self.dem.rdx)
        count = 0
        west_k, west_b = self.calculate_line_equation(self.horizontal_end_angle)
        mid_x = self.x_grid_center + 1

        y_grid_index = self.y_grid_center

        # 计算起始边界搜索起点与终点索引 需考虑斜率为零的情况
        if west_k == 0:
            start_x_index = self.min_x
        else:
            start_x_index = int((y_grid_index - west_b) / west_k) - self.is_observe_not_in_left_top_diagonal
        end_x_index = mid_x

        # 起始终止角度位于y轴
        if west_k == math.inf:
            start_x_index = self.x_grid_center
        else:
            start_x_index = self.min_x if start_x_index < self.min_x else start_x_index

        diff_x = self.x_grid_center - start_x_index

        start_distance_x_index = self.x_grid_count - diff_x - 1
        x_grid_index = start_x_index
        x_distance_index = start_distance_x_index
        y_distance_index = self.y_grid_count - 1
        last_height = self.dem.height[x_grid_index - 1, y_grid_index] - self.see_height
        p = -1 / self.y_distance[y_distance_index]

        if self.is_in_grid_point is True:
            self.point_in_line_consty(end_x_index, x_grid_index, y_grid_index, y_distance_index)
            # while x_grid_index <= end_x_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     x_grid_index += 1
            #     x_distance_index += 1
            #     count += 1
            const_decrease_x = 0 if west_k == 0 or west_k == math.inf else math.floor(-1 / west_k)

            y_grid_index -= 1
            y_distance_index -= 1

            start_x_index += const_decrease_x
            start_distance_x_index += const_decrease_x
            if start_distance_x_index < 0:
                adjust_x = abs(start_distance_x_index)
                start_distance_x_index = 0
                start_x_index += adjust_x
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index
            p = -1 / self.y_distance[y_distance_index]
        while x_grid_index <= end_x_index and x_distance_index < 2 * self.x_grid_count - 1:
            count += 1
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            x_grid_index += 1
            x_distance_index += 1

            # 直线线性变换,y均匀增加,x增加也为固定值.
        const_decrease_x = 0 if west_k == 0 or west_k == math.inf else math.floor(-1 / west_k)

        y_grid_index -= 1
        y_distance_index -= 1

        start_x_index += const_decrease_x
        start_distance_x_index += const_decrease_x

        # 下边逐层计算参考线与点的显隐性
        while y_grid_index >= self.min_y and start_distance_x_index < 2 * self.x_grid_count - 1:
            if start_distance_x_index < 0:
                adjust_x = abs(start_distance_x_index)
                start_distance_x_index = 0
                start_x_index += adjust_x
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index - 1, y_grid_index] - self.see_height

            p = -1 / self.y_distance[y_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 4 if x_grid_index >= self.x_grid_center else 3
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and y_distance_index <= self.y_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从西向东过程
            x_grid_index += 1
            x_distance_index += 1
            count += 1
            # 判断后续点的显隐性
            while x_grid_index <= self.max_x and x_grid_index <= end_x_index:
                count += 1
                # 形成新的调查线段
                d = self.x_distance[x_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 4 if x_grid_index >= self.x_grid_center else 3
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and y_distance_index <= self.y_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从西向东过程
                x_grid_index += 1
                x_distance_index += 1

            # 最外层循环 对应网格点横坐标加1：从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

            start_x_index += const_decrease_x
            start_distance_x_index += const_decrease_x

        # 右侧
        # 参考线的起始线段
        u = 1 / (self.dem.dy * self.dem.rdy)
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        north_k, north_b = self.calculate_line_equation(self.horizontal_start_angle)

        # 先将规律变化的纬度进行赋值 从而根据方程求出另一纬
        # 若观察点不在网格点上则x_gird_center是其左边的一个网格点,x_gird_center+1为要分析的第一列（没考虑到在格网点上的情况）
        x_grid_index = self.x_grid_center + 1

        # 计算起始边界搜索起点与终点索引
        start_y_index = self.min_y
        end_y_bound = math.ceil(x_grid_index * north_k + north_b)
        end_y_index = end_y_bound + self.is_observe_not_in_left_top_diagonal
        # 需要考虑在y轴上的情况 此时不能直接用方程去表示
        if north_k == math.inf:
            end_y_index = self.max_y
        else:
            end_y_index = self.max_y if end_y_index > self.max_y else end_y_index

        start_y_index = self.min_y if start_y_index < self.min_y else start_y_index
        start_distance_y_index = 0
        y_grid_index = start_y_index
        x_distance_index = self.x_grid_count
        y_distance_index = start_distance_y_index
        last_height = self.dem.height[x_grid_index, y_grid_index - 1] - self.see_height
        p = 1 / self.x_distance[x_distance_index]

        # 观察点不是网格点该算法也只是将其视为中心网格点上方或者下方的网格点 并没有对其求值

        if self.is_in_grid_point is True:
            self.point_in_line_constx(end_y_index, y_grid_index, x_grid_index, x_distance_index)
            # while y_grid_index <= end_y_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     y_grid_index += 1
            #     y_distance_index += 1
            #     count += 1
            # 因为是dy线性变化每次x都是均匀递增,这样y增加的值也是固定的
            const_increase_y = 0 if north_k == math.inf else math.ceil(north_k)

            # 目标点继续向右移动
            x_grid_index += 1
            x_distance_index += 1
            end_y_index += const_increase_y
            if start_distance_y_index < 0:
                adjust_y = abs(start_distance_y_index)
                start_distance_y_index = 0
                start_y_index += adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index
            p = 1 / self.x_distance[x_distance_index]
        # 构造右半边初始参考线
        while y_grid_index <= end_y_index and y_distance_index < 2 * self.y_grid_count - 1:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            y_grid_index += 1
            y_distance_index += 1
            count += 1

        # 因为是dy线性变化每次x都是均匀递增,这样y增加的值也是固定的
        const_increase_y = 0 if north_k == math.inf else math.ceil(north_k)

        # 目标点继续向右移动
        x_grid_index += 1
        x_distance_index += 1

        end_y_index += const_increase_y

        # 右边逐层计算参考线与点的显隐性
        while x_grid_index <= self.max_x and start_distance_y_index < 2 * self.y_grid_count - 1:
            count += 1
            # 防止超出结果数组边界,回调索引
            if start_distance_y_index < 0:
                adjust_y = abs(start_distance_y_index)
                start_distance_y_index = 0
                start_y_index += adjust_y

            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index, y_grid_index - 1] - self.see_height

            p = 1 / self.x_distance[x_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 1 if y_grid_index >= self.y_grid_center else 4
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and x_distance_index >= self.x_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从南向北过程
            y_grid_index += 1
            y_distance_index += 1
            count += 1
            # 判断后续点的显隐性别
            while y_grid_index <= self.max_y and y_grid_index <= end_y_index:
                count += 1
                # 形成新的调查线段
                d = self.y_distance[y_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 1 if y_grid_index >= self.y_grid_center else 4
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and x_distance_index >= self.x_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从南向北过程
                y_grid_index += 1
                y_distance_index += 1

            # 最外层循环 对应网格点横坐标加1：从西向东过程
            x_grid_index += 1
            x_distance_index += 1
            end_y_index += const_increase_y
        return self.result, count

    # Divided into two categories:
    # ① The part less than 180 is located below the observation point --> lateral grid line
    # ② The part greater than 180 is to the left of the observation point --> the vertical grid line
    # Divide regions into lower and left regions by three line equations: l1-> upper left region, l2-> middle region, l3-> lower right region
    def analysis_left_bottom(self):
        # 第一个点对应的e->height*p
        start_base_e = 0.0
        # 参考线的起始线段
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        # 位于左侧 --> 纵向网格线(从北向南,从东向西)
        u = 1 / (self.dem.dy * self.dem.rdy)
        count = 0
        north_k, north_b = self.calculate_line_equation(self.horizontal_end_angle)

        x_grid_index = self.x_grid_center

        # 需要考虑在y轴上的情况 此时不能直接用方程去表示
        if north_k == math.inf:
            start_y_index = self.max_y
        else:
            # 计算起始边界搜索起点与终点索引
            start_y_bound = int(x_grid_index * north_k + north_b)
            start_y_index = start_y_bound + self.is_observe_not_in_right_bottom_diagonal + 1

        end_y_index = self.min_y
        start_y_index = self.max_y if start_y_index > self.max_y else start_y_index

        # 计算起始格网点向上移动了几个
        diff_y = self.y_grid_center - start_y_index

        start_distance_y_index = self.y_grid_count - diff_y - 1

        y_grid_index = start_y_index
        x_distance_index = self.x_grid_count - 1
        y_distance_index = start_distance_y_index
        last_height = self.dem.height[x_grid_index, y_grid_index + 1] - self.see_height

        p = -1 / self.x_distance[x_distance_index]

        if self.is_in_grid_point is True:
            self.point_in_line_constx(end_y_index, y_grid_index, x_grid_index, x_distance_index)
            # while y_grid_index >= end_y_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     y_grid_index -= 1
            #     y_distance_index -= 1
            #     count += 1
            x_grid_index -= 1
            x_distance_index -= 1
            const_increase_y = 0 if north_k == math.inf else math.ceil(-north_k)

            start_distance_y_index += const_increase_y
            start_y_index += const_increase_y

            max_y_index = 2 * self.y_grid_count - 1
            if start_distance_y_index > max_y_index:
                adjust_y = abs(start_distance_y_index - max_y_index)
                start_distance_y_index = max_y_index
                start_y_index -= adjust_y

            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index
            p = -1 / self.x_distance[x_distance_index]
        while y_grid_index >= end_y_index and y_distance_index > 0:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)
            #  print(f"p:{p},y_distance_index:{y_distance_index},d:{d}")
            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next
            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            y_grid_index -= 1
            y_distance_index -= 1
            count += 1

        # 因为是dy线性变化每次x都是均匀递增,这样y增加的值也是固定的
        # 后续工作-->根据起始与终止角度 判断边界是增加还是缩小
        # const_increase_y = 0 if north_k == -1 else self.adjust_bound(north_k * self.dem.dx / self.dem.dy)

        x_grid_index -= 1
        x_distance_index -= 1
        const_increase_y = 0 if north_k == math.inf else math.ceil(-north_k)

        start_distance_y_index += const_increase_y
        start_y_index += const_increase_y
        # 左边逐层计算参考线与点的显隐性
        max_y_index = 2 * self.y_grid_count - 1
        count += 1
        while x_grid_index >= self.min_x and start_distance_y_index > 0:
            if start_distance_y_index > max_y_index:
                adjust_y = abs(start_distance_y_index - max_y_index)
                start_distance_y_index = max_y_index
                start_y_index -= adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index, y_grid_index + 1] - self.see_height

            p = -1 / self.x_distance[x_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去 负责角度最小范围
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 2 if y_grid_index >= self.y_grid_center else 3
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and x_distance_index <= self.x_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

            # 判断后续点的显隐性别
            while y_grid_index >= self.min_y and y_grid_index >= end_y_index:
                # 形成新的调查线段
                count += 1
                d = -self.y_distance[y_distance_index] * p
                # print(f"dx,dy:{x_distance_index},{y_distance_index}\nx,y{x_grid_index, y_grid_index}\nd:{d}")
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                # 负责角度最大范围 若当前d大于所有d则参考线构建错误
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 2 if y_grid_index >= self.y_grid_center else 3
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and x_distance_index <= self.x_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从北向南过程
                y_grid_index -= 1
                y_distance_index -= 1

            # 最外层循环 对应网格点横坐标加1：从东向西过程
            x_grid_index -= 1
            x_distance_index -= 1

            start_distance_y_index += const_increase_y
            start_y_index += const_increase_y
        # 下侧
        # 位于下侧 --> 横向网格线(从北向南，从西到东)
        # 参考线的起始线段
        u = 1 / (self.dem.dx * self.dem.rdx)

        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        east_k, east_b = self.calculate_line_equation(self.horizontal_start_angle)

        # 若观察点不在网格点上则x_gird_center是其左边的一个网格点,x_gird_center+1为要分析的第一列（没考虑到在格网点上的情况）
        y_grid_index = self.y_grid_center

        # 计算起始边界搜索起点与终点索引 需考虑斜率为零的情况
        start_x_index = self.x_grid_center - 1

        if east_k == 0:
            end_x_index = self.max_x
        else:
            end_x_index = math.ceil(
                (y_grid_index - east_b) / east_k) + self.is_observe_not_in_right_bottom_diagonal

        if east_k == math.inf:
            end_x_index = self.x_grid_center
        else:
            end_x_index = self.max_x if end_x_index > self.max_x else end_x_index

        start_x_index = self.min_x if start_x_index < self.min_x else start_x_index

        diff_x = self.x_grid_center - start_x_index

        start_distance_x_index = self.x_grid_count - diff_x - 1

        x_grid_index = start_x_index
        x_distance_index = start_distance_x_index
        y_distance_index = self.y_grid_count - 1
        last_height = self.dem.height[x_grid_index - 1, y_grid_index] - self.see_height
        p = -1 / self.y_distance[y_distance_index]

        if self.is_in_grid_point is True:
            self.point_in_line_consty(end_x_index, x_grid_index, y_grid_index, y_distance_index)
            # while x_grid_index <= end_x_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     x_grid_index += 1
            #     x_distance_index += 1
            #     count += 1
            y_grid_index -= 1
            y_distance_index -= 1
            # 直线线性变换,y均匀增加,x增加也为固定值.(需要根据角度来对扩展边界值进行优化)
            const_increase_x = 0 if east_k == 0 or east_k == math.inf else math.ceil(-1 / east_k)

            end_x_index += const_increase_x

            if start_distance_x_index < 0:
                adjust_x = abs(start_distance_x_index)
                start_distance_x_index = 0
                start_x_index += adjust_x

            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index
            p = -1 / self.y_distance[y_distance_index]
        while x_grid_index <= end_x_index and x_distance_index < 2 * self.x_grid_count - 1:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            last_height = current_height

            x_grid_index += 1
            x_distance_index += 1
            count += 1

        y_grid_index -= 1
        y_distance_index -= 1
        # 直线线性变换,y均匀增加,x增加也为固定值.(需要根据角度来对扩展边界值进行优化)
        const_increase_x = 0 if east_k == 0 or east_k == math.inf else math.ceil(-1 / east_k)

        end_x_index += const_increase_x

        # 下边逐层计算参考线与点的显隐性
        while y_grid_index >= self.min_y and start_distance_x_index < 2 * self.x_grid_count - 1:
            if start_distance_x_index < 0:
                adjust_x = abs(start_distance_x_index)
                start_distance_x_index = 0
                start_x_index += adjust_x
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index - 1, y_grid_index] - self.see_height

            p = -1 / self.y_distance[y_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                if base_line.next is None:
                    break
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                current_quadrant = 4 if x_grid_index >= self.x_grid_center else 3
                is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                if is_in_range and y_distance_index <= self.y_grid_count:
                    self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从西向东过程
            x_grid_index += 1
            x_distance_index += 1

            count += 1

            # 判断后续点的显隐性别
            while x_grid_index <= self.max_x and x_grid_index <= end_x_index:
                # 形成新的调查线段
                count += 1
                d = self.x_distance[x_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    current_quadrant = 4 if x_grid_index >= self.x_grid_center else 3
                    is_in_range = self.point_is_in_range(x_grid_index, y_grid_index, current_quadrant)
                    if is_in_range and y_distance_index <= self.y_grid_count:
                        self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从西向东过程
                x_grid_index += 1
                x_distance_index += 1

            # 最外层循环 对应网格点横坐标加1：从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1
            end_x_index += const_increase_x
        return self.result, count

    # Execute Xpderl in 360 degrees
    def analysis_by_xpderl(self):

        # 纵向格网间距实地距离的倒数 用它来计算线段系数a
        u = 1 / (self.dem.dy * self.dem.rdy)

        # p,d根据不同区域会有所变化，但均符合以下两点：①p始终保持正值 ②d值随着调查线的构建应该逐渐变大
        # 邻近值p->1/x
        p = 0.0

        # 当前经纬度对应的高程值
        current_height = 0.0
        last_height = 0.0

        # 方向值d->y/x
        d = 0.0
        # 用于对比旧基准线上一次的当前点
        last_d = 0.0
        # 交点方向值
        cross_d = 0.0
        # 记录一次求交运算中交点区间的最小值
        min_d = 0.0

        # 线段对应系数a->u*(current_height-last_height)
        a = 0.0

        # 第一个点对应的e->height*p
        start_base_e = 0.0

        current_base_e = 0.0

        # 参考线的起始线段
        start_base_line = None
        # 用于对比旧基准线的当前点
        base_line = None
        # 即将接入新参考线的一部分线段,在不确定该线段定义域范围时暂存于该变量
        current_new_line = None

        # 边界索引值
        x_right_index = 2 * self.x_grid_count - 1
        y_top_index = 2 * self.y_grid_count - 1

        # 当前点对应参考线和调查线的e差值 dif_e=current_e(调查线)-current_base_e(参考线)
        e_diff = 0.0
        e_diff_last = 0.0

        # 右半面
        # 以下算法认为横向间隔和纵向间隔相同
        # 右半圆: 从南往北,从西往东算
        # 通过约束纵向y索引来限制判断区域
        # 拓展边界
        start_index_adjust = -1 if self.is_observe_not_in_right_bottom_diagonal else 0
        end_index_adjust = 1 if self.is_observe_not_in_left_top_diagonal else 0

        # 约束循环条件
        start_y_index = self.y_grid_center + start_index_adjust
        end_y_index = self.y_grid_center + end_index_adjust + 1

        # 调整对应的d边界
        start_distance_y_index = self.y_grid_count + start_index_adjust - 1

        x_grid_index = self.x_grid_center + 1
        y_grid_index = start_y_index
        x_distance_index = self.x_grid_count
        y_distance_index = start_distance_y_index
        last_height = self.dem.height[x_grid_index, y_grid_index - 1] - self.see_height

        if self.is_in_grid_point is True:
            self.point_in_line_constx(end_y_index, y_grid_index, x_grid_index, x_distance_index)
            # while y_grid_index <= end_y_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     y_grid_index += 1
            #     y_distance_index += 1
            # 目标点继续向右移动
            x_grid_index += 1
            x_distance_index += 1

            # 区域内第一个观察点的y坐标
            start_distance_y_index -= 1

            # 向y方向上下拓展边界
            start_y_index -= 1

            end_y_index += 1
            if start_distance_y_index < 0:
                adjust_y = abs(start_distance_y_index)
                start_distance_y_index = 0
                start_y_index += adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index
            p = 1 / self.x_distance[x_distance_index]

        # 观察点不是网格点该算法也只是将其视为中心网格点上方或者下方的网格点 并没有对其求值
        # 构造右半边初始参考线
        while y_grid_index <= end_y_index:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            #  第一列是可见的
            self.result[x_distance_index][y_distance_index] = 1
            last_height = current_height

            y_grid_index += 1
            y_distance_index += 1

        # 区域内第一个观察点的y坐标
        start_distance_y_index -= 1

        # 向y方向上下拓展边界 x型分区
        start_y_index -= 1
        end_y_index += 1

        # 目标点继续向右移动
        x_grid_index += 1
        x_distance_index += 1

        # 右边逐层计算参考线与点的显隐性
        while x_grid_index <= self.max_x and start_distance_y_index < 2 * self.y_grid_count - 1:
            # 防止超出结果数组边界,回调索引
            if start_distance_y_index < 0:
                start_y_index += 1
                start_distance_y_index = 0

            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index, y_grid_index - 1] - self.see_height

            p = 1 / self.x_distance[x_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从南向北过程
            y_grid_index += 1
            y_distance_index += 1

            # 判断后续点的显隐性别
            while y_grid_index <= self.max_y and y_grid_index <= end_y_index:
                # 形成新的调查线段
                d = self.y_distance[y_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从南向北过程
                y_grid_index += 1
                y_distance_index += 1

            # 最外层循环 对应网格点横坐标加1：从西向东过程
            x_grid_index += 1
            x_distance_index += 1
            start_y_index -= 1
            start_distance_y_index -= 1
            end_y_index += 1
        # print("右面：\n", result)

        # print("右面是否相同：\n", judge_with_true.are_arrays_equal(result, true_result))
        # 上半面及下半面->横向格网间距实地距离的倒数 用它来计算线段系数a
        u = 1 / (self.dem.dx * self.dem.rdx)

        # 上半面
        # 上半圆:从南往北，从东往西算

        start_index_adjust = 1 if self.is_observe_not_in_right_bottom_diagonal else 0
        end_index_adjust = -1 if self.is_observe_not_in_left_top_diagonal else 0

        start_x_index = self.x_grid_center + start_index_adjust + 1
        start_distance_x_index = self.x_grid_count + start_index_adjust
        end_x_index = self.x_grid_center + end_index_adjust

        x_grid_index = start_x_index
        y_grid_index = self.y_grid_center + 1
        x_distance_index = start_distance_x_index
        y_distance_index = self.y_grid_count
        last_height = self.dem.height[x_grid_index + 1, y_grid_index] - self.see_height
        start_base_line = None

        if self.is_in_grid_point is True:
            self.point_in_line_consty(end_x_index, x_grid_index, y_grid_index, y_distance_index)
            # while x_grid_index >= end_x_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     x_grid_index -= 1
            #     x_distance_index -= 1
            y_grid_index += 1
            y_distance_index += 1

            max_right_index = 2 * self.x_grid_count - 1
            # 拓展 x 方向
            start_x_index += 1
            start_distance_x_index += 1
            end_x_index -= 1
            if start_distance_x_index > max_right_index:
                adjust_x = abs(start_distance_x_index - max_right_index)
                start_distance_x_index = max_right_index
                start_x_index -= adjust_x
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index
            p = 1 / self.y_distance[y_distance_index]
        # 构造上半边初始参考线
        while x_grid_index >= end_x_index:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            #  第一列是可见的
            # result[x_distance_index][y_distance_index] = 1
            last_height = current_height

            x_grid_index -= 1
            x_distance_index -= 1

        # 拓展 x 方向
        start_x_index += 1
        start_distance_x_index += 1
        end_x_index -= 1

        y_grid_index += 1
        y_distance_index += 1

        # 上边逐层计算参考线与点的显隐性
        # 防止超出x边界
        max_right_index = 2 * self.x_grid_count - 1
        while y_grid_index <= self.max_y:
            if start_distance_x_index > max_right_index:
                start_x_index -= 1
                start_distance_x_index = max_right_index
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index + 1, y_grid_index] - self.see_height

            p = 1 / self.y_distance[y_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                # result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从东向西过程
            x_grid_index -= 1
            x_distance_index -= 1

            # 判断后续点的显隐性别
            while x_grid_index >= self.min_x and x_grid_index >= end_x_index:
                # 形成新的调查线段
                d = -self.x_distance[x_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next

                    # 上下半面 并没有赋值为1
                    self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从东向西过程
                x_grid_index -= 1
                x_distance_index -= 1

            # 最外层循环 对应网格点横坐标加1：从南向北过程
            # 拓展 x 边界
            start_distance_x_index += 1
            start_x_index += 1
            end_x_index -= 1

            y_grid_index += 1
            y_distance_index += 1
        # print("上面：\n", result)

        u = 1 / (self.dem.dy * self.dem.rdy)

        # 左半面
        # 左半圆:从北往南，从东往西算
        start_index_adjust = 1 if self.is_observe_not_in_left_top_diagonal else 0
        end_index_adjust = -1 if self.is_observe_not_in_right_bottom_diagonal else 0

        start_y_index = self.y_grid_center + start_index_adjust + 1
        end_y_index = self.y_grid_center + end_index_adjust
        start_distance_y_index = self.y_grid_count + start_index_adjust

        x_grid_index = self.x_grid_center
        y_grid_index = start_y_index
        x_distance_index = self.x_grid_count - 1
        y_distance_index = start_distance_y_index
        last_height = self.dem.height[x_grid_index, y_grid_index + 1] - self.see_height
        start_base_line = None

        if self.is_in_grid_point is True:
            self.point_in_line_constx(end_y_index, y_grid_index, x_grid_index, x_distance_index)
            # while y_grid_index >= end_y_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     y_grid_index -= 1
            #     y_distance_index -= 1
            x_grid_index -= 1
            x_distance_index -= 1
            max_y_index = 2 * self.y_grid_count - 1
            start_distance_y_index += 1
            start_y_index += 1
            end_y_index -= 1
            if start_distance_y_index > max_y_index:
                adjust_y = abs(start_distance_y_index - max_y_index)
                start_distance_y_index = max_y_index
                start_y_index -= adjust_y
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index
            p = -1 / self.x_distance[x_distance_index]

        # 构造左半边初始参考线
        while y_grid_index >= end_y_index:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            #  第一列是可见的
            self.result[x_distance_index][y_distance_index] = 1
            last_height = current_height

            y_grid_index -= 1
            y_distance_index -= 1

        start_distance_y_index += 1
        start_y_index += 1
        end_y_index -= 1
        x_grid_index -= 1
        x_distance_index -= 1

        # 左边逐层计算参考线与点的显隐性
        max_y_index = 2 * self.y_grid_count - 1
        while x_grid_index >= self.min_x:
            if start_distance_y_index > max_y_index:
                start_distance_y_index = max_y_index
                start_y_index -= 1
            y_grid_index = start_y_index
            y_distance_index = start_distance_y_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index, y_grid_index + 1] - self.see_height

            p = -1 / self.x_distance[x_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = -self.y_distance[y_distance_index] * p
            a = u * (current_height - last_height)

            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                self.result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0
            last_d = d
            last_height = current_height

            # 从北向南过程
            y_grid_index -= 1
            y_distance_index -= 1

            # 判断后续点的显隐性别
            while y_grid_index >= self.min_y and y_grid_index >= end_y_index:
                # 形成新的调查线段
                d = -self.y_distance[y_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next
                    self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从北向南过程
                y_grid_index -= 1
                y_distance_index -= 1

            # 最外层循环 对应网格点横坐标加1：从东向西过程
            end_y_index -= 1
            start_distance_y_index += 1
            start_y_index += 1
            x_grid_index -= 1
            x_distance_index -= 1
        # print("左面：\n", result)

        # 上半面及下半面->横向格网间距实地距离的倒数 用它来计算线段系数a
        u = 1 / (self.dem.dx * self.dem.rdx)

        # 下半面
        # 下半圆：从北往南，从西往东算
        start_index_adjust = -1 if self.is_observe_not_in_left_top_diagonal else 0
        end_index_adjust = 1 if self.is_observe_not_in_right_bottom_diagonal else 0

        start_x_index = self.x_grid_center + start_index_adjust
        end_x_index = self.x_grid_center + end_index_adjust + 1
        start_distance_x_index = self.x_grid_count + start_index_adjust - 1

        x_grid_index = start_x_index
        y_grid_index = self.y_grid_center
        x_distance_index = start_distance_x_index
        y_distance_index = self.y_grid_count - 1
        last_height = self.dem.height[x_grid_index - 1, y_grid_index] - self.see_height
        start_base_line = None

        if self.is_in_grid_point is True:
            self.point_in_line_consty(end_x_index, x_grid_index, y_grid_index, y_distance_index)
            # while x_grid_index <= end_x_index:
            #     self.init_reference_line_by_r3(x_distance_index, y_distance_index)
            #     x_grid_index += 1
            #     x_distance_index += 1
            y_grid_index -= 1
            y_distance_index -= 1

            start_x_index -= 1
            start_distance_x_index -= 1
            end_x_index += 1
            if start_distance_x_index < 0:
                adjust_x = abs(start_distance_x_index)
                start_distance_x_index = 0
                start_x_index += adjust_x
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index
            p = -1 / self.y_distance[y_distance_index]
        # 构造下半边初始参考线
        while x_grid_index <= end_x_index:
            # 当目标点的高度小于观察点的高度会存在负数的情况
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            if start_base_line is None:
                # 通过记录第一段的e、每段的a和每段的结尾d,可以在PDE方程空间中表达和记录整个地形线 e=ad+b
                start_base_line = LinkedLinePDE(d, a)
                base_line = start_base_line
                start_base_e = current_height * p
            else:
                base_line.link_forward(LinkedLinePDE(d, a))
                base_line = base_line.next

            #  第一列是可见的
            # result[x_distance_index][y_distance_index] = 1
            last_height = current_height

            x_grid_index += 1
            x_distance_index += 1

        start_x_index -= 1
        start_distance_x_index -= 1
        end_x_index += 1
        y_grid_index -= 1
        y_distance_index -= 1

        # 下边逐层计算参考线与点的显隐性
        while y_grid_index >= self.min_y:
            if start_distance_x_index < 0:
                start_distance_x_index = 0
                start_x_index += 1
            x_grid_index = start_x_index
            x_distance_index = start_distance_x_index

            # 分析区域外拓展了一点
            last_height = self.dem.height[x_grid_index - 1, y_grid_index] - self.see_height

            p = -1 / self.y_distance[y_distance_index]
            base_line = start_base_line

            # 判断起点的显隐性
            # 初始化起点的d、a
            current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
            d = self.x_distance[x_distance_index] * p
            a = u * (current_height - last_height)

            """
            为什么需要next呢:因为第一个点借助了分析范围外的一点,如果出去借助线段本身的第一个线段就包含了下一列的d就把自身的参考线段作为第一个线段,出去借助的参考线段并复制e0
            为什么要比较K呢：因为对于调查线为e0对于的参考线并不是其对应的第一个点 因此需要找到对应点的参考线e
            如果参考线下一个线段末端的d值小于当前值,找到对应的△Ei的值
            k：后续点的方向d  beseLint.EndK:参考线的方向d 后者必须覆盖前者
            """
            while base_line and base_line.next and base_line.next.end_d < d:
                base_line = base_line.next
                # △Ei=△Ei-1+Ai(di-di-1) 这里仅通过链表迭代计算参考线的E,调查线未参与运算
                start_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            start_base_line = base_line
            current_base_e = start_base_e

            # 上一步相当于把延伸的点所参与的参考线段筛出去
            while base_line and base_line.end_d < d:
                if base_line.next is None:
                    break
                base_line = base_line.next
                current_base_e += (base_line.end_d - base_line.pre.end_d) * base_line.a

            """
            这是当前起点轴斜率与参考线的差值 e 此处相当于△E0
            current_base_e值等于ei+1、 base_line.a*(base_line.end_d-d)就表示其中k用的是当前目标点的k,因此获取值即为同一方向d上对应参考线的e它与ei不一定相等,它是处于[ei,ei+1]区间内的值
            """
            e_diff = current_height * p - current_base_e + base_line.a * (base_line.end_d - d)

            # 可视性判断
            if e_diff >= 0:
                # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                start_base_line = current_new_line = LinkedLinePDE(d, a)
                start_base_e = current_height * p
                # result[x_distance_index][y_distance_index] = 1
            else:
                self.result[x_distance_index][y_distance_index] = 0

            last_d = d
            last_height = current_height

            # 从西向东过程
            x_grid_index += 1
            x_distance_index += 1

            # 判断后续点的显隐性别
            while x_grid_index <= self.max_x and x_grid_index <= end_x_index:
                # 形成新的调查线段
                d = self.x_distance[x_distance_index] * p
                current_height = self.dem.height[x_grid_index, y_grid_index] - self.see_height
                a = u * (current_height - last_height)

                min_d = last_d
                # 当参考线的d及EndK小于当前k时需要进行参考线迭代运算 来实现参考线更新及对应点的可视性判断
                while base_line.end_d < d:
                    """
                    求出轴斜率差值的增量
                    这里a固定的原因是仅针对调查线中的一条线段来进行累计其斜率固定因此a不变而调查线的一段线段可能会包含多条参考线的线段，因此参考线会不停的更新
                    △Ei=△Ei-1+(di-di-1)(ALi-ARi) 这里的E就是参考线与调查线之间的E
                    """
                    e_diff_last = e_diff + (base_line.end_d - min_d) * (a - base_line.a)

                    # 判断是否为交点 如果前后可视性不同则比存在交点
                    # 当前是可视的
                    if e_diff >= 0:
                        # 变为不可视
                        if e_diff_last < 0:
                            # d很小时算不准
                            if e_diff < 5e-15:
                                cross_d = base_line.end_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                            # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                            current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                            current_new_line = current_new_line.next

                            # 这条语句同时还更新了参考线
                            current_new_line.link_forward(base_line)
                    else:
                        if e_diff_last >= 0:
                            if e_diff > -5e-15:
                                cross_d = min_d
                            else:
                                # 交点的d值
                                cross_d = min_d + e_diff / (base_line.a - a)

                            current_new_line = LinkedLinePDE(cross_d, base_line.a)
                            if base_line.pre is not None:
                                base_line.pre.link_forward(current_new_line)
                            else:
                                start_base_line = current_new_line
                                start_base_e = current_height * p - a * (d - cross_d)
                    e_diff = e_diff_last
                    min_d = base_line.end_d
                    if base_line.next is None:
                        break
                    base_line = base_line.next

                e_diff_last = e_diff + (d - min_d) * (a - base_line.a)
                if e_diff >= 0:
                    # 变为不可视
                    if e_diff_last < 0:
                        # d很小时算不准
                        if e_diff < 5e-15:
                            cross_d = d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        # 将交点所含线段插入，并移至该线段，并连接后续的参考线
                        # 之前可见现在不可见，则记录交叉点并将旧参考线插入到后面
                        current_new_line.link_forward(LinkedLinePDE(cross_d, a))
                        current_new_line = current_new_line.next

                        # 这条语句同时还更新了参考线
                        current_new_line.link_forward(base_line)
                else:
                    if e_diff_last >= 0:
                        if e_diff > -5e-15:
                            cross_d = min_d
                        else:
                            # 交点的d值
                            cross_d = min_d + e_diff / (base_line.a - a)

                        current_new_line = LinkedLinePDE(cross_d, base_line.a)
                        if base_line.pre is not None:
                            base_line.pre.link_forward(current_new_line)
                        else:
                            start_base_line = current_new_line
                            start_base_e = current_height * p - a * (d - cross_d)
                e_diff = e_diff_last

                # 可视性判断
                if e_diff >= 0:
                    # 如果未被挡 则需要更新参考线 因为它e更大所以才没被挡
                    current_new_line.link_forward(LinkedLinePDE(d, a))
                    current_new_line = current_new_line.next

                    self.result[x_distance_index][y_distance_index] = 1
                else:
                    self.result[x_distance_index][y_distance_index] = 0

                last_d = d
                last_height = current_height

                # 从西向东过程
                x_grid_index += 1
                x_distance_index += 1

            # 最外层循环 对应网格点横坐标加1：从北向南过程
            start_x_index -= 1
            start_distance_x_index -= 1
            end_x_index += 1
            y_grid_index -= 1
            y_distance_index -= 1

        return self.result, self.result.size

    # Generate the boundary line equation according to the Angle value and the position in the grid corresponding to the observation point
    def calculate_line_equation(self, angle):
        quadrant = self.judge_angle_at_quadrant(angle)
        # 根据输入角度转为对应区间与x夹角值 ×会出现错误
        # 输入角度转为对应与 x 正方向夹角值
        if quadrant == 1:
            angle = 90 - angle
        elif quadrant == 2:
            angle = 450 - angle
        elif quadrant == 3:
            angle = 270 - angle
        else:
            angle = 270 - angle
        # 将角度转换为弧度
        if angle == 90 or angle == 270:
            return math.inf, math.inf
        angle_rad = math.radians(angle)

        # 计算斜率
        slope = math.tan(angle_rad)
        # 设置一个阈值，如果计算结果接近于0，则视为0
        threshold = 1e-10  # 例如, 1e-10
        if abs(slope) < threshold:
            slope = 0

        # 解方程 y = mx + b，计算截距 b
        b = self.y_grid_observe - slope * self.x_grid_observe

        # 返回斜率和截距
        return slope, b

    # According to the current grid index (not the result array) and the Angle range, determine whether it is within the visible distance, need to consider in combination with the quadrant
    def point_is_in_range(self, x, y, current_quadrant):
        angle = 0
        # 分象限考虑 先计算当前坐标与靠近北方向的夹角

        # 位于1、3象限 需计算与y轴夹角
        if current_quadrant == 1 or current_quadrant == 3:
            # 返回弧度制
            # 考虑同y轴的情况
            if y != self.y_grid_observe:
                atan_value = math.atan(abs((x - self.x_grid_observe) / (y - self.y_grid_observe))) + self.quadrant[
                    current_quadrant - 1]
                angle = math.degrees(atan_value) % 360
            else:
                angle = 90 if x > self.x_grid_observe else 270
        # 转为角度值 可能存在负数的情况

        # 位于2、4象限 需计算与x夹角
        elif current_quadrant == 2 or current_quadrant == 4:
            if x != self.x_grid_observe:
                atan_value = math.atan(abs((y - self.y_grid_observe) / (x - self.x_grid_observe))) + self.quadrant[
                    current_quadrant - 1]
                angle = math.degrees(atan_value) % 360
            else:
                angle = 0 if y > self.x_grid_observe else 180
        if self.horizontal_start_angle <= angle <= self.horizontal_end_angle:
            return True
        else:
            return False

    # Return Cartesian quadrant based on Angle without calling properties of class itself
    @staticmethod
    def judge_angle_at_quadrant(angle):
        if 0 <= angle < 90:
            quadrant = 1
        elif 90 <= angle < 180:
            quadrant = 4
        elif 180 <= angle < 270:
            quadrant = 3
        else:
            quadrant = 2
        return quadrant

    # Construct initial reference lines using R3
    # Each time the grid coordinates are passed in, the r3 algorithm is used to determine whether the grid point is visible from the observation point
    def init_reference_line_by_r3(self, x_distance_index, y_distance_index):
        real_lon = self.dem.real_distance
        real_lat = self.dem.real_distance
        current_quadrant = 1
        # 判断格网点象限位置-第一象限
        if x_distance_index >= self.x_grid_count and y_distance_index >= self.y_grid_count:
            current_quadrant = 1
            i = x_distance_index
            j = y_distance_index
            # 观察点到目标点的距离
            r = math.sqrt(self.x_distance[i] ** 2 + self.y_distance[j] ** 2)
            # 高度比距离，观察点至目标点仰角
            s = (self.dem.height[self.min_x + i][self.min_y + j] - self.see_height) / r
            max_s = s
            # 当前观察点至目标点与横向格网线的交点数
            x_intersection = math.floor(self.x_distance[i] / real_lon)

            # 获取当前目标点单个网格 LOS对应的dx[i] 的dy[j]
            dp_lat = self.y_distance[j] / self.x_distance[i] * real_lon

            # 格网点对角线长度
            dp_l = math.sqrt(dp_lat ** 2 + real_lon ** 2)

            # 循环遍历x交点，看纵向网格线，从目标点到观察点进行计算的

            k = 1
            while k <= x_intersection:
                # 第k个交点时对应的y长度
                dy = dp_lat * k

                # 交点至观察点线的长度
                length = r - dp_l * k

                lon_index = i - k + self.min_x

                # 表示目标点在当前纬度格点之间的相对位置
                inner_y = dy % real_lat
                # 向下取整 负数的时候会与int强转有差异 -2.1 --> -3 而不是-2
                lat_index = j - math.floor(dy / real_lat) + self.min_y

                # 线性插值
                h = inner_y * (self.dem.height[lon_index][lat_index - 1] - self.dem.height[lon_index][
                    lat_index]) / real_lat + \
                    self.dem.height[lon_index][lat_index] - self.see_height

                current_s = h / length

                if current_s > max_s:
                    max_s = current_s
                k += 1

            y_intersection = math.floor(self.y_distance[j] / real_lat)
            dp_lon = self.x_distance[i] / self.y_distance[j] * real_lat
            dp_l = math.sqrt(dp_lon ** 2 + real_lat ** 2)
            k = 1
            while k <= y_intersection:
                dx = dp_lon * k
                length = r - dp_l * k

                lat_index = j - k + self.min_y
                inner_x = dx % real_lon
                lon_index = i - math.floor(dx / real_lon) + self.min_x

                h = inner_x * (self.dem.height[lon_index - 1][lat_index] - self.dem.height[lon_index][
                    lat_index]) / real_lon + \
                    self.dem.height[lon_index][lat_index] - self.see_height

                current_s = h / length
                if current_s > max_s:
                    max_s = current_s
                k += 1

            is_in_range = self.point_is_in_range(self.min_x + i, self.min_y + j, current_quadrant)
            if s >= max_s and is_in_range:
                self.result[i][j] = 1
            else:
                self.result[i][j] = 0

        # 第二象限
        if x_distance_index < self.x_grid_count and y_distance_index >= self.y_grid_count:
            current_quadrant = 2
            i = x_distance_index
            j = y_distance_index
            # 观察点到目标点的距离
            r = math.sqrt(self.x_distance[i] ** 2 + self.y_distance[j] ** 2)
            # 高度比距离，观察点至目标点仰角
            s = (self.dem.height[self.min_x + i][self.min_y + j] - self.see_height) / r
            max_s = s

            # 当前观察点至目标点与横向格网线的交点数
            x_intersection = math.floor(-self.x_distance[i] / real_lon)

            # 获取当前目标点单个网格 LOS对应的dx[i] 的dy[j]
            dp_lat = self.y_distance[j] / -self.x_distance[i] * real_lon

            # 格网点对角线长度
            dp_l = math.sqrt(dp_lat ** 2 + real_lon ** 2)

            # 循环遍历x交点，看纵向网格线，从目标点到观察点进行计算的

            k = 1
            while k <= x_intersection:
                # 第k个交点时对应的y长度
                dy = dp_lat * k

                # 交点至观察点线的长度
                length = r - dp_l * k

                lon_index = i + k + self.min_x

                # 表示目标点在当前纬度格点之间的相对位置
                inner_y = dy % real_lat
                lat_index = j - math.floor(dy / real_lat) + self.min_y

                # 线性插值
                h = inner_y * (self.dem.height[lon_index][lat_index - 1] - self.dem.height[lon_index][
                    lat_index]) / real_lat + \
                    self.dem.height[lon_index][lat_index] - self.see_height

                current_s = h / length

                if current_s > max_s:
                    max_s = current_s
                k += 1

            y_intersection = math.floor(self.y_distance[j] / real_lat)
            dp_lon = -self.x_distance[i] / self.y_distance[j] * real_lat
            dp_l = math.sqrt(dp_lon ** 2 + real_lat ** 2)
            k = 1
            while k <= y_intersection:
                dx = dp_lon * k
                length = r - dp_l * k

                lat_index = j - k + self.min_y
                inner_x = dx % real_lon
                lon_index = i + math.floor(dx / real_lon) + self.min_x

                h = inner_x * (self.dem.height[lon_index + 1][lat_index] - self.dem.height[lon_index][
                    lat_index]) / real_lon + \
                    self.dem.height[lon_index][lat_index] - self.see_height

                current_s = h / length
                if current_s > max_s:
                    max_s = current_s
                k += 1

            is_in_range = self.point_is_in_range(self.min_x + i, self.min_y + j, current_quadrant)
            if s >= max_s and is_in_range:
                self.result[i][j] = 1
            else:
                self.result[i][j] = 0

        # 第三象限
        if x_distance_index < self.x_grid_count and y_distance_index < self.y_grid_count:
            # 观察点到目标点的距离
            current_quadrant = 3
            i = x_distance_index
            j = y_distance_index
            r = math.sqrt(self.x_distance[i] ** 2 + self.y_distance[j] ** 2)
            # 高度比距离，观察点至目标点仰角
            s = (self.dem.height[self.min_x + i][self.min_y + j] - self.see_height) / r
            max_s = s

            # 当前观察点至目标点与横向格网线的交点数
            x_intersection = math.floor(-self.x_distance[i] / real_lon)

            # 获取当前目标点单个网格 LOS对应的dx[i] 的dy[j]
            dp_lat = -self.y_distance[j] / -self.x_distance[i] * real_lon

            # 格网点对角线长度
            dp_l = math.sqrt(dp_lat ** 2 + real_lon ** 2)

            # 循环遍历x交点，看纵向网格线，从目标点到观察点进行计算的

            k = 1
            while k <= x_intersection:
                # 第k个交点时对应的y长度
                dy = dp_lat * k

                # 交点至观察点线的长度
                length = r - dp_l * k

                lon_index = i + k + self.min_x

                # 表示目标点在当前纬度格点之间的相对位置
                inner_y = dy % real_lat
                lat_index = j + math.floor(dy / real_lat) + self.min_y

                # 线性插值
                h = inner_y * (self.dem.height[lon_index][lat_index + 1] - self.dem.height[lon_index][
                    lat_index]) / real_lat + \
                    self.dem.height[lon_index][lat_index] - self.see_height

                current_s = h / length

                if current_s > max_s:
                    max_s = current_s
                k += 1

            y_intersection = math.floor(-self.y_distance[j] / real_lat)
            dp_lon = -self.x_distance[i] / -self.y_distance[j] * real_lat
            dp_l = math.sqrt(dp_lon ** 2 + real_lat ** 2)
            k = 1
            while k <= y_intersection:
                dx = dp_lon * k
                length = r - dp_l * k

                lat_index = j + k + self.min_y
                inner_x = dx % real_lon
                lon_index = i + int(dx / real_lon) + self.min_x

                h = inner_x * (self.dem.height[lon_index + 1][lat_index] - self.dem.height[lon_index][
                    lat_index]) / real_lon + \
                    self.dem.height[lon_index][lat_index] - self.see_height

                current_s = h / length
                if current_s > max_s:
                    max_s = current_s
                k += 1

            is_in_range = self.point_is_in_range(self.min_x + i, self.min_y + j, current_quadrant)
            if s >= max_s and is_in_range:
                self.result[i][j] = 1
            else:
                self.result[i][j] = 0

        # 第四象限
        if x_distance_index >= self.x_grid_count and y_distance_index < self.y_grid_count:
            current_quadrant = 4
            i = x_distance_index
            j = y_distance_index
            # 观察点到目标点的距离
            r = math.sqrt(self.x_distance[i] ** 2 + self.y_distance[j] ** 2)
            # 高度比距离，观察点至目标点仰角
            s = (self.dem.height[self.min_x + i][self.min_y + j] - self.see_height) / r
            max_s = s
            # 当前观察点至目标点与横向格网线的交点数
            x_intersection = math.floor(self.x_distance[i] / real_lon)

            # 获取当前目标点单个网格 LOS对应的dx[i] 的dy[j]
            dp_lat = -self.y_distance[j] / self.x_distance[i] * real_lon

            # 格网点对角线长度
            dp_l = math.sqrt(dp_lat ** 2 + real_lon ** 2)

            # 循环遍历x交点，看纵向网格线，从目标点到观察点进行计算的
            k = 1
            while k <= x_intersection:
                # 第k个交点时对应的y长度
                dy = dp_lat * k

                # 交点至观察点线的长度
                length = r - dp_l * k

                lon_index = i - k + self.min_x

                # 表示目标点在当前纬度格点之间的相对位置
                inner_y = dy % real_lat
                lat_index = j + math.floor(dy / real_lat) + self.min_y

                # 线性插值
                h = inner_y * (self.dem.height[lon_index][lat_index + 1] - self.dem.height[lon_index][
                    lat_index]) / real_lat + \
                    self.dem.height[lon_index][lat_index] - self.see_height

                current_s = h / length

                if current_s > max_s:
                    max_s = current_s
                k += 1

            y_intersection = math.floor(-self.y_distance[j] / real_lat)
            dp_lon = self.x_distance[i] / -self.y_distance[j] * real_lat
            dp_l = math.sqrt(dp_lon ** 2 + real_lat ** 2)
            k = 1

            while k <= y_intersection:
                dx = dp_lon * k
                length = r - dp_l * k

                lat_index = j + k + self.min_y
                inner_x = dx % real_lon
                lon_index = i - math.floor(dx / real_lon) + self.min_x

                h = inner_x * (self.dem.height[lon_index - 1][lat_index] - self.dem.height[lon_index][
                    lat_index]) / real_lon + \
                    self.dem.height[lon_index][lat_index] - self.see_height

                current_s = h / length
                if current_s > max_s:
                    max_s = current_s
                k += 1

            is_in_range = self.point_is_in_range(self.min_x + i, self.min_y + j, current_quadrant)
            if s >= max_s and is_in_range:
                self.result[i][j] = 1
            else:
                self.result[i][j] = 0

    # For x fixed, determine the visibility of the same y grid point
    def point_in_line_constx(self, end_index, start_index, const_x_map_index, const_x_vis_index):
        max_h_up = 0
        max_h_down = 0
        index_up = self.y_grid_count + 1
        index_down = self.y_grid_count - 1
        map_index_up = self.y_grid_center + 1
        map_index_down = self.y_grid_center - 1
        max_k = -math.inf
        while map_index_up <= end_index and index_up < self.y_grid_count:
            cur_h = self.dem.height[const_x_map_index][map_index_up]
            distance = (map_index_up - self.y_grid_center) * self.dem.dy
            k = (cur_h - self.see_height) / distance
            # if cur_h > max_h_up:
            #     max_h_up = cur_h
            #     self.result[const_x_vis_index][index_up] = 1
            if k > max_k:
                max_k = k
                self.result[const_x_vis_index][index_up] = 1
            map_index_up += 1
            index_up += 1

        max_k = -math.inf
        while map_index_down >= start_index and index_down >= 0:
            cur_h = self.dem.height[const_x_map_index][map_index_down]
            distance = (map_index_down - self.y_grid_center) * self.dem.dy
            k = (cur_h - self.see_height) / distance
            if k > max_k:
                max_k = k
                self.result[const_x_vis_index][index_down] = 1
            # if cur_h > max_h_down:
            #     max_h_down = cur_h
            #     self.result[const_x_vis_index][index_down] = 1
            map_index_down -= 1
            index_down -= 1

    def point_in_line_consty(self, end_index, start_index, const_y_map_index, const_y_vis_index):
        max_h_up = 0
        max_h_down = 0
        index_up = self.x_grid_count + 1
        index_down = self.x_grid_count - 1
        map_index_up = self.x_grid_center + 1
        map_index_down = self.x_grid_center - 1
        max_k = -math.inf
        while map_index_up <= end_index and index_up < self.x_grid_count:
            cur_h = self.dem.height[map_index_up][const_y_map_index]
            distance = (map_index_up - self.x_grid_center) * self.dem.dx
            k = (cur_h - self.see_height) / distance
            if k > max_k:
                max_k = k
                self.result[index_up][const_y_vis_index] = 1
            # if cur_h > max_h_up:
            #     max_h_up = cur_h
            #     self.result[index_up][const_y_vis_index] = 1
            map_index_up += 1
            index_up += 1

        max_k = -math.inf
        while map_index_down >= start_index and index_down >= 0:
            cur_h = self.dem.height[map_index_down][const_y_map_index]
            distance = (map_index_down - self.x_grid_center) * self.dem.dx
            k = (cur_h - self.see_height) / distance
            if k > max_k:
                max_k = k
                self.result[index_down][const_y_vis_index] = 1
            # if cur_h > max_h_down:
            #     max_h_down = cur_h
            #     self.result[index_down][const_y_vis_index] = 1
            map_index_down -= 1
            index_down -= 1