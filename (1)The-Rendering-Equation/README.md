![](media/pbr-equation.jpg)

# 从渲染方程开始...
渲染方程描述了着色点沿一特定观察方向出射的光的总量与入射光和BRDF的关系。它来源于物理中的能量守恒定律，形式简洁，含义直观。1986年James Kajiya将之引入计算机图形学，概括了当时存在的一些渲染算法，例如分布式光线追踪和辐射度方法，并为后来的路径追踪等算法提供了直接的理论基础。在今天，渲染方程已经成为基于物理的渲染中各种复杂光照算法的出发点，让我们从学习渲染方程开始，一窥光照传输的门径。
## 辐射度量学
在基于物理的渲染中，我们需要一种方法对光进行量化。辐射度量学研究的是对电磁辐射的测量，与光线的物理传输过程紧密相关，渲染方程中出现的各个物理量就属于辐射度量学。这些物理量还需要光度学和色度学的知识，才会转化成我们熟悉的RGB颜色值。<br>渲染中常用的辐射度量有energy、flux/power、intensity、irradiance和radiance（这几个物理量基本都保留英文，不做翻译）。有关它们的详细定义请见其它材料，下面我用一个例子展示它们的作用和区别。
### 例1
![](media/the-life-of-light.png)
<p align="center">图1 一条由点光源、着色点、眼睛构成的简单光路</p>

考虑图1中的情况：对于点光源，我们无法定义它的出射radiance（因为没有面积），但是它的flux（记为 $\Phi$）均匀分布在半径为 $d$的球面上，球面面积为 $4\pi d^2$。因此，球面上单位面积的光通量（理想情况的irradiance）为 $\frac{\Phi}{4\pi d^2}$。再引入入射角的影响，着色点 $p$接收到的irradiance为 $\frac{\Phi}{4\pi d^2} \cdot (n \cdot \omega_i)$。表面的反射作用由BRDF（记为 $f_r(p,\omega_i,\omega_o)$）建模。最终出射radiance为 $\frac{\Phi}{4\pi d^2} \cdot (n \cdot \omega_i) \cdot f_r(p,\omega_i,\omega_o)$，这也是眼睛接收到的radiance。<br>几个需要注意的问题：<br>1.点光源无法定义出射radiance，对于某一方向，可以用intensity，值为 $I=\frac{\Phi}{4\pi}$。<br>2.BRDF的定义是出射的radiance除以入射的irradiance，为了得到后者我们需要将入射光乘以 $cos\theta$。<br>3.radiance是一根光线的特性，着色点出射的radiance与到达眼睛的radiance是相等的，没有平方衰减，距离增加的效果是在屏幕中对应的像素减少。
## 渲染方程