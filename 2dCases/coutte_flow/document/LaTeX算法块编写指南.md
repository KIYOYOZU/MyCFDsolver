# LaTeX 算法块编写指南与常见错误

## 1. 宏包选择

### ✅ 推荐配置
```latex
\usepackage{algorithm}
\usepackage{algpseudocode}  % algorithmicx 家族，命令首字母大写
```

### ❌ 避免使用
```latex
\usepackage{algorithmic}  % 旧包，命令全大写如 \STATE, \IF
```

**原因：** `algpseudocode` 提供更现代的语法，命令更直观（`\State`, `\If`, `\While`）

---

## 2. 浮动位置设置

### ✅ 推荐
```latex
\begin{algorithm}[!htbp]  % 允许LaTeX自动调整位置
\caption{算法名称}
\label{alg:mylabel}
\begin{algorithmic}[1]
...
\end{algorithmic}
\end{algorithm}
```

### ❌ 易出错
```latex
\begin{algorithm}[H]  % 强制当前位置
```

**问题：** `[H]` 强制算法在当前位置显示，若内容过长会超出页面导致显示不全

**浮动参数说明：**
- `h`: 尽量在当前位置（here）
- `t`: 页面顶部（top）
- `b`: 页面底部（bottom）
- `p`: 独立页面（page）
- `!`: 忽略LaTeX的严格限制

---

## 3. 内容过长的处理方法

### 方法一：缩小字体
```latex
\begin{algorithm}[!htbp]
\caption{算法名称}
\footnotesize  % 或 \small
\begin{algorithmic}[1]
...
\end{algorithmic}
\end{algorithm}
```

### 方法二：拆分算法
```latex
\begin{algorithm}[!htbp]
\caption{算法1：初始化与计算}
...
\end{algorithm}

\begin{algorithm}[!htbp]
\caption{算法2：后处理与输出}
...
\end{algorithm}
```

---

## 4. 编号与标题设置

### ✅ 推荐写法
```latex
\section{算法设计与实现}

\begin{algorithm}[!htbp]
\caption{算法1：主求解器}  % 编号直接在caption中
\label{alg:solver}
\begin{algorithmic}[1]
...
\end{algorithmic}
\end{algorithm}
```

### ❌ 易出错写法
```latex
\section{算法设计与实现}

\subsection{算法1：主求解器}  % 独立的subsection标题

\begin{algorithm}[!htbp]
\caption{主求解器}
...
\end{algorithm}
```

**问题：** 使用浮动参数时，算法块可能浮动到subsection标题之前，导致标题出现在算法块之后

---

## 5. 算法风格：注重思想而非细节

### ✅ 好的写法（高层次抽象）
```latex
\begin{algorithmic}[1]
\State \textbf{初始化阶段：}
\State 计算网格间距 $\Delta x, \Delta y$ 和扩散数 $\alpha_y = \nu \Delta t / (\Delta y)^2$；
\State 检查稳定性条件 $\alpha_y \leq 0.5$，若违反则终止；
\State
\State \textbf{时间推进循环：}
\While{$t < t_{\text{end}}$}
    \State 对所有内部点应用显式格式：$u_j^{n+1} = u_j^n + \alpha_y(u_{j+1}^n - 2u_j^n + u_{j-1}^n)$；
    \State 应用边界条件；
    \State 更新时间 $t \gets t + \Delta t$；
\EndWhile
\end{algorithmic}
```

### ❌ 不推荐的写法（过度细节）
```latex
\begin{algorithmic}[1]
\State $\Delta x \gets L_x/(N_x-1)$
\State $\Delta y \gets L_y/(N_y-1)$
\State $\alpha_y \gets \nu \cdot \Delta t / (\Delta y)^2$
\For{$i=1$ \textbf{to} $N_x$}
    \For{$j=2$ \textbf{to} $N_y-1$}
        \State $u\_new(i,j) \gets u(i,j) + \alpha_y \cdot [u(i,j+1) - 2u(i,j) + u(i,j-1)]$
    \EndFor
\EndFor
\end{algorithmic}
```

**原因：** 学术论文中的算法应展示核心思想和数学公式，而非编程实现细节

---

## 6. 中文化设置

```latex
% 在导言区添加
\algrenewcommand\algorithmicrequire{\textbf{输入:}}
\algrenewcommand\algorithmicensure{\textbf{输出:}}
\floatname{algorithm}{算法}
```

**使用效果：**
```latex
\begin{algorithmic}[1]
\Require 参数 $x, y, z$
\Ensure 结果 $result$
...
\end{algorithmic}
```
显示为：
- **输入:** 参数 $x, y, z$
- **输出:** 结果 $result$

---

## 7. 常用命令对照表

| algpseudocode（推荐） | algorithmic（旧） | 说明 |
|---------------------|------------------|------|
| `\State` | `\STATE` | 单行语句 |
| `\If{条件}...\EndIf` | `\IF{条件}...\ENDIF` | 条件语句 |
| `\While{条件}...\EndWhile` | `\WHILE{条件}...\ENDWHILE` | 循环 |
| `\For{条件}...\EndFor` | `\FOR{条件}...\ENDFOR` | For循环 |
| `\Return` | `\RETURN` | 返回 |
| `\Require` | `\REQUIRE` | 输入 |
| `\Ensure` | `\ENSURE` | 输出 |

---

## 8. 完整示例模板

```latex
\documentclass[12pt,a4paper]{ctexart}

\usepackage{amsmath, amssymb}
\usepackage{algorithm}
\usepackage{algpseudocode}

% 中文化
\algrenewcommand\algorithmicrequire{\textbf{输入:}}
\algrenewcommand\algorithmicensure{\textbf{输出:}}
\floatname{algorithm}{算法}

\begin{document}

\section{算法设计}

\begin{algorithm}[!htbp]
\caption{算法1：显式求解器}
\label{alg:explicit_solver}
\begin{algorithmic}[1]
\Require 网格参数 $N_x, N_y$；物性参数 $\nu$；时间步长 $\Delta t$
\Ensure 速度场 $u(x,y,t)$
\State
\State \textbf{初始化：}
\State 计算扩散数 $\alpha = \nu \Delta t / (\Delta x)^2$；
\State 检查稳定性条件 $\alpha \leq 0.5$；
\State
\State \textbf{时间推进：}
\While{$t < t_{\text{end}}$}
    \State 对内部点应用离散格式；
    \State 应用边界条件；
    \State 更新时间 $t \gets t + \Delta t$；
\EndWhile
\State
\Return 速度场历史数据
\end{algorithmic}
\end{algorithm}

\end{document}
```

---

## 9. 常见编译错误及解决

| 错误信息 | 原因 | 解决方法 |
|---------|------|---------|
| `Undefined control sequence \State` | 使用了 `algorithmic` 包但写了 `\State` | 改用 `algpseudocode` 包 |
| `Float too large for page` | 算法内容过长超出页面 | 使用 `\footnotesize` 或拆分算法 |
| 算法标题出现在算法块之后 | 使用了独立的 `\subsection` 与浮动参数 | 将编号整合到 `\caption` 中 |
| `Missing $ inserted` | 数学符号未用 `$...$` 包裹 | 检查所有变量是否在数学模式中 |

---

## 10. 最佳实践总结

1. ✅ 使用 `algpseudocode` 包（命令首字母大写）
2. ✅ 使用 `[!htbp]` 浮动参数而非 `[H]`
3. ✅ 内容过长时用 `\footnotesize` 或拆分算法
4. ✅ 编号直接写在 `\caption` 中，避免独立 `\subsection`
5. ✅ 注重算法思想，避免过度实现细节
6. ✅ 合理使用分组（`\textbf{阶段名：}`）提升可读性
7. ✅ 数学公式用 `$...$` 包裹，复杂公式用 `\qquad` 缩进
8. ✅ 中文文档记得设置中文化命令
9. ✅ 多次编译以更新交叉引用
10. ✅ 参考优秀论文的算法块格式

---

**创建日期：** 2025-12-04
**适用场景：** 学术论文、技术报告、毕业论文
