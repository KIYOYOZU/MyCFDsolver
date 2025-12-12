# myCFD

<div align="center">

![MATLAB](https://img.shields.io/badge/MATLAB-R2018b+-orange.svg)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey.svg)
![Status](https://img.shields.io/badge/status-Active-success.svg)

**面向 CFD 教学和研究的 MATLAB 数值模拟框架**

*A MATLAB-based CFD Framework for Teaching and Research*

[特性](#主要特性) • [快速开始](#快速开始) • [算例](#算例说明) • [文档](#文档)

</div>

---

## 项目简介 | Project Overview

**myCFD** 是一个专为计算流体力学（CFD）教学和研究设计的 MATLAB 代码库。项目采用**分层架构设计**，将通用求解器框架（`src/`）与具体算例（`2dCases/`）分离，具有良好的可扩展性和可维护性。

当前实现了 **Couette 流启动过程**的数值模拟，包含完整的解析解验证、误差分析和可视化功能，可作为 CFD 入门的优秀教学案例。

**myCFD** is a MATLAB-based computational fluid dynamics (CFD) framework designed for teaching and research. It features a **layered architecture** that separates the generic solver framework (`src/`) from specific test cases (`2dCases/`), ensuring excellent extensibility and maintainability.

Currently implements numerical simulation of **Couette flow startup process**, with complete analytical solution verification, error analysis, and visualization capabilities, serving as an excellent teaching case for CFD beginners.

---

## 主要特性 | Key Features

- **模块化架构**: 通用求解器框架与算例分离，易于扩展新的物理问题
- **面向对象设计**: 采用抽象基类和继承机制，支持多种求解器和边界条件
- **高性能计算**: 充分利用 MATLAB 向量化运算，计算速度达 ~157,000 步/秒
- **完整验证**: 与解析解对比，L2 相对误差 < 1e-5，精度达 99.9993%
- **丰富可视化**: 自动生成演化图、对比图、云图和动画，支持高质量输出（300 DPI）
- **教学友好**: 提供单文件简化版本和完整 LaTeX 教程文档

---

## 目录结构 | Directory Structure

```
myCFD/
├── src/                    # 通用求解器框架（核心库）
│   └── 2d/                 # 2D 求解器模块
│       ├── solvers/        # 求解器抽象基类和具体实现
│       │   ├── BaseSolver.m              # 抽象基类
│       │   └── ExplicitEulerSolver.m     # 显式欧拉求解器
│       ├── boundary/       # 边界条件模块
│       │   ├── BoundaryCondition.m       # 抽象基类
│       │   ├── DirichletBC.m             # Dirichlet 边界条件
│       │   ├── PeriodicBC.m              # 周期边界条件
│       │   └── CompositeBoundaryCondition.m  # 组合边界条件
│       ├── mesh/           # 网格生成模块
│       │   └── Mesh2D.m                  # 2D 结构化网格
│       └── utils/          # 工具函数
│           └── compute_error.m           # 误差分析
│
└── 2dCases/                # 算例目录
    └── coutte_flow/        # Couette 流算例
        ├── run_couette.m           # 主计算脚本（入口）
        ├── postprocess_couette.m   # 后处理脚本（可视化）
        ├── CouetteAnalytical.m     # 解析解类
        ├── simplified_couette.m    # 单文件教学版
        ├── couette_tutorial.tex    # LaTeX 教程文档
        └── results/                # 计算结果（自动生成）
```

---

## 快速开始 | Quick Start

### 环境要求

- **MATLAB**: R2018b 或更高版本
- **操作系统**: Windows / Linux / macOS
- **内存**: 建议 4GB 以上

### 运行算例

#### 方法 1: 模块化版本（推荐）

```bash
# 计算阶段
matlab -batch "cd('2dCases/coutte_flow'); run_couette"

# 后处理阶段（生成图形和动画）
matlab -batch "cd('2dCases/coutte_flow'); postprocess_couette"
```

#### 方法 2: 单文件版本（教学演示）

```bash
matlab -batch "cd('2dCases/coutte_flow'); simplified_couette"
```

#### 方法 3: MATLAB GUI

```matlab
% 在 MATLAB 命令窗口中执行
cd 2dCases/coutte_flow
run_couette          % 运行计算
postprocess_couette  % 生成可视化
```

### 预期输出

运行成功后，将在 `2dCases/coutte_flow/results/` 目录下生成：

- **演化图**: `velocity_evolution.png` - 速度场随时间演化
- **对比图**: `comparison_with_analytical.png` - 数值解与解析解对比
- **云图**: `velocity_contour.png` - 速度场云图
- **动画**: `couette_animation.mp4` - 流场演化动画
- **数据文件**: `simulation_results.mat` - 完整计算结果

---

## 算例说明 | Test Cases

### Couette 流

**物理问题**: 两平行平板间的剪切流动，上板以恒定速度 U 运动，下板静止。

**控制方程**:
```
∂u/∂t = ν · ∂²u/∂y²
```

**边界条件**:
- 下壁面 (y=0): u = 0 (无滑移)
- 上壁面 (y=h): u = U (运动壁面)
- X 方向: 周期性边界条件

**数值方法**:
- **时间离散**: 显式欧拉格式 (O(Δt))
- **空间离散**: 中心差分法 (O(Δy²))
- **稳定性条件**: α_y = νΔt/Δy² ≤ 0.5

**解析解**: Fourier 级数展开（100 项求和）
```
u(y,t) = U·(y/h) + (2U/π)·Σ[(-1)^n/n · sin(nπy/h) · exp(-n²π²νt/h²)]
```

**性能指标**:

| 指标 | 数值 |
|------|------|
| 总时间步数 | 200,001 步 |
| 计算耗时 | ~1.3 秒 |
| 计算速度 | ~157,000 步/秒 |
| 最终 L2 误差 | 6.85×10⁻⁶ |
| 相对精度 | 99.9993% |
| 内存占用 | ~50 MB |

---

## 技术栈 | Tech Stack

- **编程语言**: MATLAB (R2018b+)
- **架构模式**: 面向对象 + 模块化分层
- **数值方法**: 有限差分法（显式欧拉 + 中心差分）
- **可视化**: MATLAB Graphics + VideoWriter
- **文档**: LaTeX + Markdown

---

## 代码示例 | Code Examples

### 创建求解器

```matlab
% 创建网格
mesh = Mesh2D(Nx, Ny, Lx, Ly);

% 配置物理参数
physical_params = struct('nu', 0.1);

% 配置时间参数
time_params = struct('dt', 0.0001, 't_end', 20.0);

% 创建求解器
solver = ExplicitEulerSolver(mesh, physical_params, time_params);

% 检查稳定性
solver.check_stability();
```

### 设置边界条件

```matlab
% 创建组合边界条件
bc = CompositeBoundaryCondition();

% 添加 Dirichlet 边界条件（上下壁面）
bc.add_condition(DirichletBC('bottom', 0.0));  % 下壁面 u=0
bc.add_condition(DirichletBC('top', U));       % 上壁面 u=U

% 添加周期边界条件（左右边界）
bc.add_condition(PeriodicBC('x'));
```

### 运行模拟

```matlab
% 初始化速度场
u = zeros(mesh.Nx, mesh.Ny);

% 运行模拟
[u_final, u_history, t_snapshots] = solver.run_simulation(u, bc);

% 计算误差
analytical = CouetteAnalytical(U, nu, h);
u_exact = analytical.compute(mesh.y, t_end);
error_L2 = compute_error(u_final(:, :), u_exact, 'L2');
```

---

## 扩展开发 | Extension Guide

### 添加新求解器

1. 继承 `BaseSolver` 抽象基类
2. 实现 `time_step()` 和 `check_stability()` 方法
3. 放置在 `src/2d/solvers/` 目录

```matlab
classdef NewSolver < BaseSolver
    methods
        function u_new = time_step(obj, u, dt)
            % 实现时间推进逻辑
        end

        function check_stability(obj)
            % 实现稳定性检查
        end
    end
end
```

### 添加新算例

1. 在 `2dCases/` 下创建新目录
2. 复制 `coutte_flow/` 作为模板
3. 修改物理参数和边界条件
4. 编写算例说明文档 `README.md`

---

## 文档 | Documentation

### 现有文档资源

| 文档 | 位置 | 说明 |
|------|------|------|
| **项目指南** | `CLAUDE.md` | 完整的开发指南和代码规范 |
| **算例说明** | `2dCases/coutte_flow/README.md` | Couette 流算例详细文档 |
| **理论推导** | `2dCases/coutte_flow/couette_tutorial.tex` | LaTeX 教学文档（理论+数值方法） |

### 编译 LaTeX 文档

```bash
cd 2dCases/coutte_flow
xelatex couette_tutorial.tex
xelatex couette_tutorial.tex  # 第二次生成完整目录和交叉引用
```

---

## 性能优化 | Performance

### 向量化运算

```matlab
% ✓ 推荐：向量化运算
u_new(2:Nx-1, 2:Ny-1) = u(2:Nx-1, 2:Ny-1) + alpha_y * (...
    u(2:Nx-1, 3:Ny) - 2*u(2:Nx-1, 2:Ny-1) + u(2:Nx-1, 1:Ny-2));

% ✗ 避免：显式循环
for i = 2:Nx-1
    for j = 2:Ny-1
        u_new(i,j) = ...;
    end
end
```

### 句柄类使用

```matlab
% 大型数据结构使用 handle 基类，避免深拷贝
classdef Mesh2D < handle
    properties
        x, y, Nx, Ny
    end
end
```

---

## 未来规划 | Roadmap

### 高优先级

- [ ] 1D 优化版本（性能提升 20 倍）
- [ ] 单元测试框架（MATLAB Unit Testing）
- [ ] 版本控制配置（`.gitignore`）

### 中优先级

- [ ] Poiseuille-Couette 流（添加压力梯度）
- [ ] 高阶时间格式（RK2 / Crank-Nicolson）
- [ ] 自适应时间步长（CFL 条件）

### 长期规划

- [ ] 3D 求解器模块
- [ ] 非结构化网格支持
- [ ] 并行计算（Parallel Computing Toolbox）

---

## 贡献指南 | Contributing

欢迎贡献代码、报告问题或提出建议！

1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启 Pull Request

### 提交前检查清单

- [ ] 代码通过稳定性检查
- [ ] 精度验证通过（L2 误差 < 1e-5）
- [ ] 无硬编码路径
- [ ] 注释覆盖率 > 60%
- [ ] 生成的图形质量良好
- [ ] 未提交 `results/` 目录

---

## 许可证 | License

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 致谢 | Acknowledgments

- 感谢所有为 CFD 教学和开源社区做出贡献的研究者
- 本项目受经典 CFD 教材和开源项目启发

---

## 联系方式 | Contact

如有问题或建议，欢迎通过以下方式联系：

- **Issues**: [GitHub Issues](https://github.com/yourusername/myCFD/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/myCFD/discussions)

---

<div align="center">

**⭐ 如果这个项目对您有帮助，请给个 Star！**

**⭐ If this project helps you, please give it a Star!**

Made with ❤️ for CFD Education

</div>
