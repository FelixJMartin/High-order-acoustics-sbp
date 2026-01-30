# High-order-acoustics-sbp

High-order **SBP–Projection finite difference methods** for the acoustic wave equation in **1D and 2D**, with emphasis on stability, well-posed boundary conditions, interfaces, convergence studies, and RK4 time integration.

---

## Acoustic Wave Equation – SBP-Projection Methods

- **Course:** 1TD354 – Scientific Computing for Partial Differential Equations  
- **Institution:** Uppsala University  
- **Language:** MATLAB  

---

## 📌 Project Overview

This project investigates the numerical solution of the **acoustic wave equation** in both **one and two spatial dimensions** using **provably stable high-order SBP–Projection finite difference methods**.

The main focus areas are:
- Well-posed boundary and interface conditions  
- Stability via energy estimates and eigenvalue analysis  
- High-order convergence studies  
- Wave propagation in heterogeneous media  

Time integration is performed using the **classical 4th-order Runge–Kutta method (RK4)**.

<p align="center">
  <img src="Figures/Assignement_2/Ass_03.png" width="70%">
</p>

---

## 📁 Repository Structure

```text
theory/
├── Assignement 1
├── Assignement 3

Code/
├── Assignement 1
├── Assignement 3

Figures/
├── Assignement 1
├── Assignement 3

```


---

## 🧪 Assignment 1 — 1D Acoustic Wave Equation

### Tasks
1. **Well-posedness analysis**
   - Derive admissible boundary conditions
   - Compare Dirichlet and characteristic BCs
   - Present continuous energy estimates

2. **SBP-Projection discretization**
   - Central SBP operators (D₁)
   - Upwind SBP operators (D⁺ / D⁻)
   - Discrete energy stability proof

3. **Numerical stability**
   - Eigenvalue analysis of discretization matrix `M`
   - Plot eigenvalues of `hM`
   - Compare Dirichlet vs characteristic BCs
   - Determine CFL numbers for RK4

4. **Convergence study**
   - Validate against analytic solution
   - Compute discrete L₂-errors
   - Measure convergence rates (6th & 7th order)

### Numerical stability study of eigenvalues:

![Eigenvalues](Figures/Assignement_1/task3_eigvals_Mh.jpg)


## 🧊 Assignment 3 — Two-Dimensional Acoustic Wave Equation

This assignment extends the SBP-Projection framework to **two spatial dimensions**. The focus is on **well-posedness**, **efficient implementation using Kronecker products**, and **wave propagation in both homogeneous and heterogeneous media**.

---

### 🔹 Task 1: 2D Acoustic Wave Equation on a Square Domain

We solve the two-dimensional acoustic wave equation on the square domain:
\[
\Omega = [-1,\,1] \times [-1,\,1]
\]

The governing equations are:
\[
C u_t + A u_x + B u_y + D u = 0
\]

with the solution vector:
\[
u = \begin{bmatrix} p & v & w \end{bmatrix}^T
\]

#### Boundary Conditions
Normal velocity is set to zero on all boundaries:
- \( v = 0 \) on \( \partial\Omega_E \cup \partial\Omega_W \)
- \( w = 0 \) on \( \partial\Omega_S \cup \partial\Omega_N \)

These boundary conditions are shown to yield a **well-posed initial boundary-value problem**.

#### Initial Conditions
\[
p(x,y,0) = e^{-100(x^2 + y^2)}, \quad
v(x,y,0) = 0, \quad
w(x,y,0) = 0
\]

Material parameters are set to:
- \( \rho = 1 \)
- \( c = 1 \)
- \( \beta = 0 \)

---

### 🔹 Discretization and Numerical Method

The domain is discretized using an **equidistant \( m \times m \) grid** with spacing:
\[
h = \frac{2}{m-1}
\]

To efficiently construct the 2D operators, **Kronecker products** are used:
\[
D_x = D \otimes I, \qquad D_y = I \otimes D
\]

The semi-discrete SBP-Projection formulation is:
\[
u_t = - P C^{-1} (D_x + D_y + D) P u
\]

All matrices are stored in **sparse format** to reduce memory usage and computational cost.

---

### 🔹 Numerical Setup
- SBP operators: **7th order accurate upwind**
- Grid size: \( m = 200 \)
- Time integration: **4th order Runge–Kutta (RK4)**
- CFL number: **0.05**

---

### 🔹 Plots to Present

The pressure field \( p(x,y,t) \) is visualized at selected time instances:

```text
figures/assignment3/p_t0.png
figures/assignment3/p_t1.png
figures/assignment3/p_t2.png



