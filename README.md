# Cubic SplineInterpolation

## clamped boundary condition

The $m_0$ and $m_n​$ are given and meet:
$$
S'(x_0)=m_0, S'(x_n)=m_n
$$

## natural boundary conditions

The $M_0$ and $M_n​$ are given and meet:
$$
S''(x_0)=M_0,S''(x_n)=M_n \\
2m_0+m_1=3f[x_0,x_1]-\frac{h_1}{2}M_0 \\
m_{n-1}+2m_n=3f[x_{n-1},x_n]+\frac{h_n}{2}M_n
$$

## periodic boundary conditions(unfinished)

$$
S''(x_0)=S''(x_n),S'(x_0)=S'(x_n)
$$

### Test LaTeX
<img src="https://latex.codecogs.com/png.latex?\dpi{150}&space;S'(x_0)=m_0,&space;S'(x_n)=m_n"/>