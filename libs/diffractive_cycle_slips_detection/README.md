# Cycle-Slip Detection Using Total Variation Denoising (TVD)

This document explains the mathematical development behind our TVD-based cycle-slip detection approach and shows how to determine the number of cycle slips from the optimization output.

---

## 1. Problem Overview

Cycle slips appear as abrupt jumps in a phase signal, for instance, in GNSS or phase tracking applications. Our goal is to estimate a denoised phase signal \(x\) from a noisy measurement \(y\) such that \(x\) is close to \(y\) while being piecewise constant. The discontinuities in \(x\) correspond to cycle slips.

---

## 2. Total Variation Denoising (TVD)

We formulate the problem as:

\[
\min_{x} \; \frac{1}{2}\|x - y\|_2^2 + \lambda \sum_{i=1}^{n-1} |x_{i+1} - x_i|
\]

- **Data Fidelity:** \(\frac{1}{2}\|x - y\|_2^2\) ensures \(x\) remains close to \(y\).
- **Regularization (TV):** \(\lambda \sum_{i=1}^{n-1} |x_{i+1} - x_i|\) encourages a piecewise constant solution, highlighting abrupt jumps (cycle slips).

---

## 3. Variable Splitting

Because the absolute value is nonsmooth, we introduce auxiliary variables \(z\) such that:

\[
z_i \geq |x_{i+1} - x_i|, \quad i=1,\dots,n-1.
\]

The reformulated problem becomes:

\[
\begin{array}{ll}
\min\limits_{x,z} & \frac{1}{2}\|x - y\|_2^2 + \lambda \sum_{i=1}^{n-1} z_i \\
\text{s.t.} & -z_i \leq x_{i+1} - x_i \leq z_i, \quad i=1,\dots,n-1, \\
            & z_i \geq 0.
\end{array}
\]

At the optimum, we expect \(z_i = |x_{i+1} - x_i|\).

---

## 4. Building the Constraints Using a Difference Operator

### 4.1 Difference Operator

We define a matrix \(E\) that computes differences:
\[
E \, x = \begin{bmatrix} x_2 - x_1 \\ x_3 - x_2 \\ \vdots \\ x_n - x_{n-1} \end{bmatrix}
\]
For a signal \(x \in \mathbb{R}^n\), \(E\) is an \((n-1) \times n\) matrix with:
- \(-1\) on the main diagonal,
- \(+1\) on the first superdiagonal.

In MATLAB, we create \(E\) as:
```matlab
m = n - 1;
E = spdiags([-ones(m,1), ones(m,1)], [0,1], m, n);
```

### 4.2 Expressing the Constraints

We combine \(x\) and \(z\) into one variable:
\[
u = \begin{bmatrix} x \\ z \end{bmatrix}.
\]
The constraints are:
- \(x(i+1) - x(i) - z_i \le 0\) → in matrix form: \([E \quad -I]\,u \le 0\),
- \(-\bigl(x(i+1)-x(i)\bigr) - z_i \le 0\) → in matrix form: \([-E \quad -I]\,u \le 0\).

Stacking these gives:
\[
A = \begin{bmatrix} E & -I \\ -E & -I \end{bmatrix}, \quad b = 0.
\]

---

## 5. Downsampling

When working with high-frequency signals (e.g., 30,000 samples for a 5-minute 100 Hz signal), downsampling reduces computational load. If a downsampling factor \(dsFactor > 1\) is provided, the signal is downsampled using MATLAB’s `downsample` function.

---

## 6. The Complete Optimization Problem

The final optimization problem solved by the function is:

\[
\begin{array}{ll}
\min\limits_{x,z} & \frac{1}{2}\|x - y_{ds}\|_2^2 + \lambda \sum_{i=1}^{n-1} z_i \\
\text{s.t.} & [E \quad -I]\,u \le 0, \\
            & [-E \quad -I]\,u \le 0, \\
            & z \ge 0.
\end{array}
\]

Here, \(y_{ds}\) is the (optionally downsampled) observed signal.

---

## 7. Cycle-Slip Detection

After solving for \(u = [x; z]\):
- The vector \(z\) gives the absolute differences \(|x(i+1)-x(i)|\).
- To detect cycle slips, set a threshold. If \(z_i\) exceeds the threshold, a cycle slip is deemed to have occurred at that index.

For example, if cycle slips correspond to significant jumps (e.g., near multiples of \(2\pi\)), you might use:
```matlab
threshold = 1.5; % adjust as appropriate
cycle_slip_count = sum(z_est > threshold);
```

---

## 8. Summary

- **TVD Objective:** Balances fidelity and smoothness, enforcing piecewise constant behavior.
- **Variable Splitting:** Introduces \(z\) to handle the nonsmooth \(|x(i+1)-x(i)|\) term.
- **Constraints:** Built using a sparse difference operator \(E\) and identity matrices.
- **Downsampling:** Reduces problem size for high-frequency signals.
- **Cycle Slip Detection:** Significant \(z\) values (above a threshold) indicate cycle slips.

This theoretical framework underpins the MATLAB implementation for detecting cycle slips using TVD. Adjust the threshold based on your specific signal characteristics to accurately count the cycle slips.

