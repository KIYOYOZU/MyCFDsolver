# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述 | Project Overview

**myCFD** 是一个面向 CFD（计算流体力学）教学和研究的 MATLAB 代码库。当前实现了 Couette 流启动过程的数值模拟，采用**分层架构设计**：通用求解器框架（src/）与具体算例（2dCases/）分离，具有良好的可扩展性。

- **主要语言**: MATLAB (R2018b+)
- **架构模式**: 面向对象 + 模块化分层
- **开源许可**: MIT License

---

## 目录结构与职责 | Directory Structure

```
myCFD/
├── src/                    # 通用求解器框架（核心库）
│   └── 2d/                 # 2D 求解器模块
│       ├── solvers/        # 求解器抽象基类和具体实现
│       ├── boundary/       # 边界条件模块（抽象基类+具体实现）
│       ├── mesh/           # 网格生成模块
│       └── utils/          # 工具函数（误差分析等）
│
└── 2dCases/                # 算例目录（独立的物理问题）
    └── coutte_flow/        # Couette 流算例
        ├── run_couette.m           # 主计算脚本（入口）
        ├── postprocess_couette.m   # 后处理脚本（可视化）
        ├── CouetteAnalytical.m     # 解析解类（算例专用）
        ├── simplified_couette.m    # 单文件教学版（备份）
        └── results/                # 计算结果（自动生成，不纳入版本控制）
```

### 职责划分

| 目录 | 职责 | 修改原则 |
|------|------|----------|
| **src/** | 通用求解器框架，可复用于多个算例 | 不可破坏向后兼容性和抽象接口 |
| **2dCases/** | 具体物理问题实现 | 可自由添加新算例子目录 |
| **results/** | 自动生成的计算结果 | 不纳入版本控制（建议添加到 .gitignore） |

---

## 代码规范 | Coding Standards

### 性能优化原则

1. **向量化优先**：避免显式 for 循环，充分利用 MATLAB 矩阵运算
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

2. **句柄类使用**：大型数据结构使用 `handle` 基类，避免深拷贝
   ```matlab
   classdef Mesh2D < handle
       % 避免每次传递参数时拷贝整个网格数据
   end
   ```

### 面向对象设计规范

1. **抽象基类与继承**
   - 抽象基类使用 `Abstract` 关键字定义接口
   - 子类必须实现所有抽象方法
   ```matlab
   % 抽象基类示例
   classdef (Abstract) BaseSolver < handle
       methods (Abstract)
           u_new = time_step(obj, u, dt)
           check_stability(obj)
       end
   end
   ```

2. **类命名与文件组织**
   - 类文件名必须与类名一致（如 `ExplicitEulerSolver.m`）
   - 使用驼峰命名法（PascalCase）

### 注释与文档规范

- **单行注释**: 使用 `%`
- **章节注释**: 使用 `%%` 分隔代码块
- **函数文档**: 在函数声明后添加 Help Text
  ```matlab
  function result = my_function(param)
      % MY_FUNCTION - 函数简要描述
      %
      % 输入:
      %   param - 参数描述
      %
      % 输出:
      %   result - 返回值描述
  ```

---

## 开发工作流 | Development Workflow

### 运行与测试流程

```bash
# 方法1：模块化版本（推荐用于开发）
matlab -batch "cd('2dCases/coutte_flow'); run_couette"         # 计算阶段
matlab -batch "cd('2dCases/coutte_flow'); postprocess_couette"  # 后处理阶段

# 方法2：单文件版本（教学演示）
matlab -batch "cd('2dCases/coutte_flow'); simplified_couette"
```

### 验证要求

任何修改后必须通过以下验证：

1. **稳定性检查**（自动）
   ```matlab
   alpha_y = nu * dt / dy^2;
   if alpha_y > 0.5
       error('稳定性条件违反：α_y > 0.5');
   end
   ```

2. **精度验证**（与解析解对比）
   - L2 相对误差应 < 1e-5
   - 稳态速度剖面应为线性分布

3. **可视化输出检查**
   - 确保生成 4 类图形：演化图、对比图、云图、动画
   - 检查图形质量（300 DPI，无警告）

### 路径管理规范

算例脚本必须自动添加 src/ 到搜索路径：
```matlab
% 在 run_couette.m 开头
case_dir = fileparts(mfilename('fullpath'));
src_path = fullfile(case_dir, '..', '..', 'src');
addpath(genpath(src_path));
```

---

## 扩展指南 | Extension Guide

### 添加新求解器

1. **创建求解器类**（继承 `BaseSolver`）
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

2. **文件位置**: `src/2d/solvers/NewSolver.m`

3. **测试要求**: 使用 Couette 流算例验证精度和稳定性

### 添加新边界条件

1. **创建边界条件类**（继承 `BoundaryCondition`）
   ```matlab
   classdef NeumannBC < BoundaryCondition
       methods
           function u = apply(obj, u)
               % 实现边界条件施加逻辑
           end
       end
   end
   ```

2. **文件位置**: `src/2d/boundary/NeumannBC.m`

3. **集成方式**: 通过 `CompositeBoundaryCondition` 组合多个边界条件

### 创建新算例

1. **目录结构**（复制 `2dCases/coutte_flow` 作为模板）
   ```
   2dCases/new_case/
   ├── run_new_case.m          # 主计算脚本
   ├── postprocess_new_case.m  # 后处理脚本
   ├── NewCaseAnalytical.m     # 解析解或基准解（如适用）
   └── README.md               # 算例说明文档
   ```

2. **必须包含内容**
   - 物理问题描述
   - 控制方程和边界条件
   - 数值方法和参数配置
   - 验证结果（与解析解或文献对比）

3. **文档要求**
   - 每个算例必须有 README.md
   - 关键理论推导建议写入 .tex 文件

---

## 常见开发任务 | Common Tasks

### 修改物理参数

在 `run_couette.m` 中修改：
```matlab
% 物理参数
U  = 1.0;   % 上壁速度 (m/s)
nu = 0.1;   % 运动粘度 (m²/s)

% 时间参数
dt = 0.0001;  % 时间步长 (s)
t_end = 20.0; % 终止时间 (s)
```

**注意**: 修改 `dt` 或 `nu` 后，必须验证稳定性条件 `alpha_y ≤ 0.5`

### 调整网格分辨率

```matlab
% 网格参数
Nx = 32;   % X方向网格点数（Couette流中保持32即可）
Ny = 51;   % Y方向网格点数（影响精度，建议 51-101）
```

**网格加密原则**:
- Y方向加密可显著提高精度
- 同时需减小 `dt` 以保持 `alpha_y ≤ 0.5`

### 添加新的可视化

在 `postprocess_couette.m` 中添加绘图代码：
```matlab
% 绘制新图形
figure;
plot(y, u_final, 'b-', 'LineWidth', 2);
xlabel('Y Position (m)');
ylabel('Velocity (m/s)');
title('Final Velocity Profile');
saveas(gcf, fullfile(results_dir, 'new_figure.png'));
```

---

## 核心架构说明 | Core Architecture

### 求解器设计模式

```
BaseSolver (抽象基类)
├── 通用属性: grid_info, physical_params, time_params
├── 抽象方法: time_step(), check_stability()
└── 通用方法: run_simulation()

ExplicitEulerSolver (具体实现)
├── 时间推进: 显式欧拉格式 u^{n+1} = u^n + Δt·L(u^n)
├── 空间离散: 中心差分法（二阶精度）
└── 稳定性: α_y = νΔt/Δy² ≤ 0.5
```

### 边界条件管理

```
CompositeBoundaryCondition (组合管理器)
├── 按顺序应用多个边界条件
└── 支持的边界条件类型:
    ├── DirichletBC (固定值): u = g(x,y,t)
    ├── PeriodicBC (周期性): u(0) = u(L)
    └── [可扩展]: NeumannBC, RobinBC 等
```

### 数值方法细节

| 项 | 方法 | 精度 |
|---|------|------|
| **时间离散** | 显式欧拉 | O(Δt) |
| **空间离散** | 中心差分 | O(Δy²) |
| **稳定性条件** | α_y = νΔt/Δy² ≤ 0.5 | CFL 条件 |

---

## 故障排查 | Troubleshooting

### 1. 稳定性问题

**症状**: 程序崩溃或数值爆炸，提示"稳定性条件违反"

**解决方案**:
```matlab
% 减小时间步长
dt_max = 0.5 * dy^2 / nu;
dt = 0.8 * dt_max;  % 使用 80% 的最大允许值
```

### 2. 精度问题

**症状**: L2 误差 > 1e-5，数值解与解析解偏差大

**原因**:
- 网格太粗（Ny < 51）
- 时间步长太大
- 计算时间不足（未达到稳态）

**解决方案**:
```matlab
Ny = 101;            % 加密网格
dt = 0.00005;        % 减小时间步长
t_end = 40.0;        % 延长计算时间（2倍特征扩散时间）
```

### 3. 性能问题

**症状**: 计算速度 < 50,000 步/秒

**检查项**:
- 是否使用了显式循环（应使用向量化）
- 是否频繁分配内存（预分配数组）
- 是否在循环中进行了不必要的计算

**优化建议**:
```matlab
% ✓ 预分配数组
u_history = cell(n_snapshots, 1);

% ✓ 向量化运算
diffusion = u(:, 3:Ny) - 2*u(:, 2:Ny-1) + u(:, 1:Ny-2);

% ✗ 避免在循环中执行耗时操作
% 不要在时间步循环中调用 disp() 或 fprintf()
```

---

## 禁止事项 | Restrictions

### 严格禁止

1. **破坏抽象接口**: 不得修改 `BaseSolver` 或 `BoundaryCondition` 的抽象方法签名
2. **硬编码路径**: 算例代码中不得包含绝对路径（如 `C:\Users\...`）
3. **提交结果文件**: 不得将 `results/` 目录纳入版本控制
4. **违反稳定性条件**: 不得使用 `alpha_y > 0.5` 的参数组合

### 强烈不推荐

1. **在 src/ 中引入算例特定代码**: 保持框架的通用性
2. **使用全局变量**: 所有参数应通过类属性或函数参数传递
3. **忽略误差分析**: 任何新算例必须提供验证结果
4. **过度优化**: 保持代码可读性，避免过早优化

---

## 数值方法说明 | Numerical Methods

### 控制方程

```matlab
% 一维非定常扩散方程（含源项）
∂u/∂t = ν · ∂²u/∂y² + S
```

### 离散格式

```matlab
% 显式欧拉时间推进
u^{n+1} = u^n + Δt · [ν · (u^n_{j+1} - 2u^n_j + u^n_{j-1})/Δy² + S]

% 稳定性条件（Von Neumann 分析）
α_y = ν·Δt/Δy² ≤ 0.5
```

### 解析解（Couette 流）

```matlab
% Fourier 级数展开（100 项求和）
u(y,t) = U·(y/h) + (2U/π)·Σ[(-1)^n/n · sin(nπy/h) · exp(-n²π²νt/h²)]
         ￣￣￣￣￣   ￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣￣
         稳态分量              瞬态分量（指数衰减）
```

---

## 性能基准 | Performance Benchmarks

当前算例（Couette 流）性能指标：

| 指标 | 数值 |
|------|------|
| 总时间步数 | 200,001 步 |
| 计算耗时 | ~1.3 秒 |
| 计算速度 | ~157,000 步/秒 |
| 最终 L2 误差 | 6.85×10⁻⁶ |
| 相对精度 | 99.9993% |
| 内存占用 | ~50 MB |

---

## 文档与教程 | Documentation

### 现有文档资源

| 文档 | 位置 | 说明 |
|------|------|------|
| **算例说明** | `2dCases/coutte_flow/README.md` | 完整的 Couette 流算例文档 |
| **理论推导** | `2dCases/coutte_flow/couette_tutorial.tex` | LaTeX 教学文档（理论+数值方法） |
| **LaTeX 指南** | `2dCases/coutte_flow/document/LaTeX算法块编写指南.md` | 算法环境使用指南 |

### 编译 LaTeX 文档

```bash
cd 2dCases/coutte_flow
xelatex couette_tutorial.tex
xelatex couette_tutorial.tex  # 第二次生成完整目录和交叉引用
```

---

## 扩展方向建议 | Future Extensions

### 高优先级

1. **1D 优化版本**: 移除冗余的 X 方向，性能提升 20 倍
2. **单元测试框架**: 使用 MATLAB Unit Testing Framework
3. **版本控制配置**: 添加 `.gitignore` 文件

### 中优先级

4. **Poiseuille-Couette 流**: 添加压力梯度源项
5. **高阶时间格式**: 实现 RK2 或 Crank-Nicolson 格式
6. **自适应时间步长**: 根据 CFL 条件动态调整 dt

### 长期规划

7. **3D 求解器模块**: 实现 `src/3d/` 目录
8. **非结构化网格**: 支持三角形/四面体网格
9. **并行计算**: 使用 MATLAB Parallel Computing Toolbox

---

## Git 工作流建议 | Git Workflow

### .gitignore 配置（建议添加）

```gitignore
# MATLAB 自动生成文件
*.asv
*.mex*

# 计算结果（不纳入版本控制）
results/
*.mat
*.png
*.mp4

# LaTeX 临时文件
*.aux
*.log
*.out
*.toc
*.synctex.gz

# 系统文件
.DS_Store
Thumbs.db
```

### 提交前检查清单

- [ ] 代码通过稳定性检查
- [ ] 精度验证通过（L2 误差 < 1e-5）
- [ ] 无硬编码路径
- [ ] 注释覆盖率 > 60%
- [ ] 生成的图形质量良好
- [ ] 未提交 `results/` 目录

---

**最后更新**: 2025-12-12
**适用版本**: myCFD v1.0
**MATLAB 版本**: R2018b+
