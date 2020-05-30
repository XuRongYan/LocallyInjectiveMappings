# Locally Injective Mappings

### 1 Motivation

在该论文之前的许多变形方法都没有考虑出现模型面片翻转的问题，即没有保证映射$f$的单射性质，在变形过程中出现翻转就说明存在点$p$与点$q$，使得$f(p) = f(q), p \neq q$。该论文通过在原本的能量函数中添加一个障碍函数来惩罚体积或者面积为0的元素来避免翻转，并提出了一个优化方法来优化该能量，使得能量能够快速收敛。

### 2 Problem Statement

定义映射$ f: \mathbb{R}^{3} \rightarrow \mathbb{R}^{3} $，该映射必须保证为单射，即没有发生翻转的四面体

### 3 Method

对于一般变形问题的通用模型可以用以下优化问题来表示：
$$
\underset{\mathbf{v}}{\operatorname{argmin}} E(\mathbf{v})+\alpha\|\mathbf{C} \mathbf{v}-\mathbf{d}\|^{2} \tag{1}
$$
其中$E(\mathbf{v})$为任意的能量，第二项为位置约束，$ \alpha \geqslant 0 $ 为位置约束的权重,，$\mathbf{C}$取单位矩阵，而该方法在该模型的基础上加上一个障碍函数，障碍函数在一个元素的面积或者体积逼近0时变得无穷大来惩罚翻转：
$$
\underset{\mathbf{v}}{\operatorname{argmin}} \quad E(\mathbf{v})+\alpha\|\mathbf{C} \mathbf{v}-\mathbf{d}\|^{2}+\beta \sum_{j \in \mathcal{E}} \phi_{j}\left(c_{j}(\mathbf{v})\right) \tag{2}
$$
其中$\beta > 0$是于控制障碍函数的惩罚力度的，$c_{j}(\mathbf{v})$为对第j个元素体积/面积的约束：
$$
c_{j}(\mathbf{v})=\lambda_{j}(\mathbf{v})-\varepsilon \tag{3}
$$
$ \lambda_{j}(\mathbf{v}) $为元素的体积/面积，而$ \varepsilon = 10^{-5} \cdot \min_{j}{\lambda_{j}(\mathbf{v}_{0})} $，$ \lambda_{j}(\mathbf{v}_{0}) $为初始状态下模型第j个元素的体积/面积。障碍函数$\phi_j$是由三次多项式$g_{j}$，该式子满足$g_{j}(0)=0, \quad g_{j}\left(s_{j}\right)=1, \quad g_{j}^{\prime}\left(s_{j}\right)=0, \quad g_{j}^{\prime \prime}\left(s_{j}\right)=0$；$g_{j}$可表示为：
$$
g_{j}(x)=\frac{1}{s_{j}^{3}} x^{3}-\frac{3}{s_{j}^{2}} x^{2}+\frac{3}{s_{j}} x \tag{4}
$$


参数$s_{j}$表示面积在触发障碍之前可以缩小到多小，本方法中取$s_{j} \in (0.1, 1)$。



### 4 Optimization

该方法通过LM方法对公式(2)进行优化，优化过程可表示为以下两步：

1. $\mathbf{p}_{i}=\left(\nabla^{2} E_{b a r}\left(\mathbf{v}_{i}\right)+\mu_{i} \mathbf{I}\right)^{-1} \cdot \nabla E_{b a r}\left(\mathbf{v}_{i}\right)$
2. $\mathbf{v}_{i+1}=\mathbf{v}_{i}-\sigma_{i} \mathbf{p}_{i}$

其中$\sigma_{i}$与$\mu_{i}$运用线性搜索的方法进行确定（将会在**附录**中进行说明）。$\sigma_{i}$代表迭代的步长，而$\mu_{i}$在尽可能小的同时，保证Hessian矩阵可逆。

### 5 Substepping

由于公式(2)对于一些比较极端的情况不能给出理想的解，所以该方法通过松弛位置约束的方法来使的求解更加稳定：
$$
E_{b a r}(\mathbf{v})=E(\mathbf{v})+\alpha\left\|\mathbf{C} \mathbf{v}-\mathbf{t}_{i}\right\|^{2}+\beta \sum_{j \in \mathcal{E}} \phi_{j}\left(c_{j}(\mathbf{v})\right) \tag{5}
$$
$\mathbf{t_{i}}$将在每一迭代中进行更新，如果Hessian矩阵为可逆的，那么$\mathbf{t}_{i} = \mathbf{d}$，$\mathbf{t}_{i}$可表示为：
$$
\mathbf{t}_{i}=\mathbf{C} \mathbf{v}_{i-1}+\frac{1}{1+\mu_{i}^{2}}\left(\mathbf{d}-\mathbf{C} \mathbf{v}_{i-1}\right) \tag{6}
$$

### 6 Computation of $\alpha$ In Equation 1

$$
\gamma = \frac{r}{\| \mathbf{Cv - d} \|^{2}} (E(\mathbf{v}) + \beta \sum_{j \in \epsilon}{\phi_{j}(c_{j}(\mathbf{v}))}) \tag{7}
$$

$$
\alpha_{i + 1} = \min{(\max(\alpha_{i}, \gamma), t)} \tag{8}
$$

其中$r$ 为控制软约束影响力的参数，而取$r = 1000, t = 10^{16}$能得到比较好的效果。



### 附录

### 1 Line Search of $\mu_{i}$

![image-20200530144256474](/Users/xry/Library/Application Support/typora-user-images/image-20200530144256474.png)

### 2 Line Search of $\sigma_{i}$

![image-20200530141717549](/Users/xry/Library/Application Support/typora-user-images/image-20200530141717549.png)

### 3 Hessian and Gradient of Barries function

$$
\nabla\phi_{j}(\mathbf{v}) = -\frac{g_{j}^{\prime}}{g_{j}^{2}} \cdot \nabla \lambda_{j}(\mathbf{v}), \tag{1}
$$

$$
\nabla^{2} \phi_{j}(\mathbf{v}) = \frac{2(g_{j}^{\prime})^{2} - g_{j}^{\prime \prime}g_{j}^{\prime}}{g_{j}^{3}} \cdot \nabla \lambda_{j}(\mathbf{v}) \cdot \nabla \lambda_{j}(\mathbf{v})^{T} - \frac{g_{j}^{\prime}}{g_{j}^{2}} \cdot \nabla^{2} \lambda_{j}(\mathbf{v}). \tag{2}
$$

$$
\nabla \lambda_{j}(\mathbf{v}) = \frac{1}{2}(\mathbf{M} + \mathbf{M}^{T}) \mathbf{v} \tag{3}
$$

$$
\nabla^{2}\lambda_{j}(\mathbf{v}) = \frac{1}{2}(\mathbf{M} + \mathbf{M}^{T}) \tag{4}
$$

$\mathbf{v}$为网格的顶点向量，可以表示为：
$$
\mathbf{v} = (v_{1x}, v_{1y}, v_{1z}, v_{2x}, v_{2y}, v_{2z}, \cdots, v_{nx}, v_{ny}, v_{nz})^{T} \tag{5}
$$
若已知三角形面片三点$v_{1}, v_{2}, v_{3}$则三角形面积可以表示为：
$$
S_{\triangle} = \frac{1}{2}（\mathbf{v_{2} - v_{1}}） \cdot (\mathbf{v_{3} - v_{1}}) \tag{6}
$$
则可将其整理为二次型形式：
$$
S_{\triangle} = \frac{1}{2}\mathbf{v_{j}}^{T} \mathbf{M}_{j} \mathbf{v_{j}} \tag{7}
$$
则将每个面片的$\mathbf{M}_{j}$拼接可得到$\mathbf{M}$

### 4 Hessian and Gradient of Soft Constrain

$$
f_{1} = \|\mathbf{Cv - t_{i}}\|^{2} \tag{8}
$$

$$
\nabla f_{1} = 2\mathbf{C^{T}Cv - 2C^{T}t_{i}} \tag{9}
$$

$$
\nabla^{2}f_{1} = 2\mathbf{C^{T}C} \tag{10}
$$



### 5 Hessian and Gradient of $E(\mathbf{v})$

拟采用ARAP能量作为$E(\mathbf{v})$：
$$
E(\mathbf{v}) = \sum_{j \in \epsilon} \| \mathbf{G_{j}v - R_{j}} \|^{2} \tag{11}
$$

$$
\nabla E(\mathbf{v}) = 2\mathbf{G^{T}Gv} - 2\mathbf{G^{T}R} \tag{12}
$$

$$
\nabla^{2}E(\mathbf{v}) = 2\mathbf{G^{T}G} \tag{13}
$$


