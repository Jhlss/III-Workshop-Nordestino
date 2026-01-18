# III-Workshop-Nordestino
Notebook related to the poster presented at the III Nordestino Workshop.

## We are developing a package to calculate the discretized Koopman operator to obtain a computational approach of error for 
diffusion coefficients of systems with noise from Adjoint Methods.

```julia
using RigorousAdjointKoopman#This is a development version. To install the latest stable version, use:
# using Pkg
# Pkg.add("RigorousAdjointKoopman")
``` 
## Example:

Let consider $D:[0,1]\to [0,1]$, where $T(x) = 2x + 0.01\sin(2\pi x) \mod(1)$.
```julia
function D(x) 
    if 0<= x < 1/2
        return 2*x + 0.01*sin(2*pi*x)
    elseif 1/2 <= x <= 1
       return (2*x - 1) + 0.01*sin(2*pi*x)
    end
end
```
```julia
using Plots
plot(D, 0:0.01:1, label="D(x)", title="Perturbed Doubling Map", xlabel="x", ylabel="D(x)")
```
## Let’s declare the basis of Hat functions.

$$
h_i(x) =
\begin{cases}
N({x - x_{i-1}}), & x \in [x_{i-1}, x_i], \\
N({x_{i+1} - x}), & x \in [x_i, x_{i+1}], \\
0, & \text{otherwise}.
\end{cases}
$$

```julia
N = 256
B = HatNP(N)
```
## Using the discretization via the projection operator $\Pi_{\delta} \sum_{i=1}^N \phi(x_i) h_i(x)$, we can approximate $K$


```julia
K_delta = assemble_Koopman_boundary_condition(B, D, 1/4)#Here when we get a big noise the matrix is not markov. Arrange the assemble.
```

#### Assume there exists a stationary probability measure $\mu$ such that $\mu K=\mu$.

For an observable $\phi:[-1,1]\to\R$, where $\phi(x) = x^2$ we have:


```julia
ϕ(x) = x^2
```
 Using the discretization via the projection operator again.

```julia
psi_delta = [ϕ(i/N) for i in 0:N]
```

The invariance $\mu K=\mu$ implies, for any test function $g$,
$$
\mu\big((Id-K)g\big)=0.
$$

### Poisson equation for the Koopman operator

Let $K$ be the Koopman operator and  $\phi$ an observable.
We want to find a function $g$ and a constant $c$ such that

$$
(Id - K) g = \phi - c,
$$

together with the normalization condition

$$
\int g \, d\mu = 0.
$$

Integrating both sides of the Poisson equation and using the invariance of
the measure $\mu$, we obtain

$$
c = \int \phi \, d\mu.
$$

#### Solving the Poisson equation allows you to decompose:


$$
\phi - \int \phi \, d\mu  = g -g\circ T,
$$

which is fundamental for:

- Ergodic sums

- Central limit theorems

- Variance formulas

- Rigorous error bounds


#### Approximate solution and residual

Suppose
- $\tilde g$ is an approximate solution.

- $\tilde c$ is a approximate mean.


Define the residual

$$
r := (Id - K)\tilde g - (\phi - \tilde c).
$$


This measures how badly the Poisson equation is violated.


Then, integrating with respect to $\mu$,

$$
\int r \, d\mu
=
\int (Id - K)\tilde g \, d\mu
-
\int (\phi - \tilde c) \, d\mu.
$$

Since $\mu$ is invariant under $K$,


- $\int (Id - K)\tilde g \, d\mu = 0$,


and therefore

$$
c - \tilde c = \int r \, d\mu.
$$
The error in the approximation $\tilde c$ of the true average $c$
is exactly the integral of the residual:

$$
|c - \tilde c| \le \int |r_{\delta}| \, d\mu + \int|r-r_{\delta}|.
$$
Thus, controlling the residual gives a rigorous error bound on the computed  by $\tilde r_\delta$.

### We estimate $\tilde{c}$ by using Birkhoff averages:

```julia
sum_c = 0
M = 100000
X = rand()
for i in 1:M
    X = π_B_C(D(X) + 0.1*(2*rand()-1))
    sum_c += ϕ(D(X))
end
tilde_c  = sum_c / M
```

### Now let’s calculate a candidate for $g$. A good candidate $g$ satisfies:




If the system is mixing and  $\phi \in L^2(\mu)$, a formal solution of the
Poisson equation
$$
(Id - K) g = \phi - \int \phi \, d\mu
$$
is given by the Neumann series:

$$
g
=
\sum_{n=0}^{\infty} K^n \bigl(\phi - \int \phi \, d\mu).
$$

#### To use the computed we are estimating the values of were finite that we can work with the computer.


$$
\tilde g = \sum_{n=0}^{N} K^n \bigl( \phi - \tilde c \bigr)
$$


```julia
using IntervalArithmetic
tilde_g = zeros(N+1)
Kc = mid.(K_delta)
w = psi_delta .- tilde_c 
for i in 1:1000
    tilde_g += w
    w = Kc*w
end
return tilde_g
```

### Then we have that $\tilde r_{\delta}$:

```julia
tilde_r = maximum([sup(x) for x in abs.((psi_delta.-tilde_c)+(K_delta*tilde_g-tilde_g))])#this is the rigorous bound on the error infinity norm.
```
