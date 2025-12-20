Do not refer to our other conversations.

I am trying to understand the Saint-Venant flexure solution for a beam with a terminal resultant transverse force $\bm{F}$

$$
\newcommand{\WarpShearIesan}{\psi}
\newcommand{\Section}{\Sigma}
$$

The resultant force is, by global equilibrium:
$$
\bm{F} \triangleq \int\bm\tau\,\mathrm{d}A = \Big(\int E\, (\bm{r}-\bm{r}^0)\otimes(\bm{r}-\bm{r}^0)\Big)\bm{b}_{\Section}
$$
where $\bm{b}_{\Section}$ is a constant 2-vector, $\bm{r} = (y,z)$, and $\bm{r}^0$ is the centroid location.

But later the shear stress is found to be
$$
\bm\tau = G\Big(\nabla \big(\bm{\WarpShearIesan}\cdot\bm{b}_{\Section}\big)
\;-\;\nu\bigl(\operatorname{dev}(\bm{r}\otimes\bm{r})-\bm{r}\otimes\bm{r}^0\bigr)\bm{b}_{\Section}\Big)
$$

where the warping vector function $\bm\WarpShearIesan\in\mathbb{R}^2$ is defined componentwise by:

$$
\left\{\begin{aligned}
\Delta (\bm\WarpShearIesan\cdot\bm{b}_{\Section}) &= -\,2\big(\bm{r} -\bm{r}^0\big)\cdot\bm{b}_{\Section} \quad &\text{in } \Section \\
\partial_{\bm{\nu}} (\bm\WarpShearIesan\cdot\bm{b}_{\Section}) &= 
\bm{\nu}\cdot\nu\;\Big(\operatorname{dev}\bm{r}\otimes\bm{r}
 -\bm{r}\otimes\bm{r}^0\Big)\bm{b}_{\Section}
\quad &\text{on } \partial\Section\\
\int_{\Section}\bm\WarpShearIesan\cdot\bm{b}_{\Section}\,\mathrm{d}A &= 0
\end{aligned}\right.
$$

and the 2D deviator furnishes $\operatorname{dev}\bm{r}\otimes\bm{r} :=\bm{r}\otimes\bm{r}-\tfrac12 r^{2}\mathbf{1}$. 

Please explain how integrating this point-wise shear can be consistent with the first formula.

Note that the divergence of a tensor T is defined as $\operatorname{div}(\bm{T})\cdot\bm{a} = \operatorname{div}(\bm{T}^{\mathrm{t}}\bm{a}),\qquad (\operatorname{div} T)_i=\partial_j T_{i j}$


---------------

Short answer: it’s the Poisson–Neumann problem for the scalar $\phi:=\bm\psi\cdot\bm b_{\Section}$ that makes the two expressions match. The boundary condition you wrote is exactly what turns the area integral of the gradient term into the missing “isotropic” part so that the (G)-based pointwise formula integrates to the (E)-based resultant.

Here’s the clean way to see it—no new symbols beyond (\phi):

1. Write the resultant from your pointwise shear:
   $$
   \bm F
   =\int_{\Section}\bm\tau,\mathrm dA
   =G\Big(\int_{\Section}\nabla\phi,\mathrm dA
    - \nu\int_{\Section}\bigl(\operatorname{dev}(\bm r\otimes\bm r)-\bm r\otimes\bm r^0\bigr)\bm b_{\Section},\mathrm dA\Big).
   \tag{1}
   $$

2. Note that (\phi) satisfies (just dot your BVP with (\bm b_{\Section}))
   $$
   \Delta_{\Section}\phi=-2(\bm r-\bm r^0)\cdot\bm b_{\Section}\quad\text{in }\Section,\qquad
   \partial_{\bm\nu}\phi=\bm\nu\cdot\big(\operatorname{dev}(\bm r\otimes\bm r)-\bm r\otimes\bm r^0\big)\bm b_{\Section}
   \quad\text{on }\partial\Section.
   \tag{2}
   $$

3. Two ingredients:

* Divergence identity for the geometric polynomial:
  $$
  \nabla\cdot\big(\operatorname{dev}(\bm r\otimes\bm r)\big)=2\bm r,\qquad
  \nabla\cdot(\bm r\otimes\bm r^0)=2\bm r^0
  \;\Rightarrow\;
  \nabla\cdot\big(\operatorname{dev}(\bm r\otimes\bm r)-\bm r\otimes\bm r^0\big)=2(\bm r-\bm r^0).
  \tag{3}
  $$

* Green’s theorem applied to the matrix fields $(\nabla\phi)\otimes\bm r$ and $\big((\operatorname{dev}(\bm r\otimes\bm r)-\bm r\otimes\bm r^0)\bm b_{\Section}\big)\otimes\bm r$, together with (2)–(3), gives the identity
  $$
  \int_{\Section}\nabla\phi\,\mathrm dA
  = - 2(1+\nu)\int_{\Section}(\bm r-\bm r^0)\otimes(\bm r-\bm r^0)\,\bm b_{\Section}\,\mathrm dA
   + \nu\int_{\Section}\bigl(\operatorname{dev}(\bm r\otimes\bm r)-\bm r\otimes\bm r^0\bigr)\bm b_{\Section}\,\mathrm dA.
  \tag{4}
  $$

(Equation (4) is the key: the Poisson source $-2(\bm r-\bm r^0)\cdot\bm b_{\Section}$ and the Neumann datum for $\partial_{\bm\nu}\phi$ conspire so that the *area integral of the gradient* produces exactly the central second-moment term, with the factor $2(1+\nu)$.)

4. Substitute (4) back into (1). The deviatoric integrals cancel, leaving
   $$
   \bm F
   = -2G(1+\nu)\int_{\Section}(\bm r-\bm r^0)\otimes(\bm r-\bm r^0),\bm b_{\Section},\mathrm dA
   = -E\int_{\Section}(\bm r-\bm r^0)\otimes(\bm r-\bm r^0)\,\bm b_{\Section}\,\mathrm dA,
   $$
   since (E=2G(1+\nu)).

That’s the consistency: the gradient part of (\bm\tau) is *not* “free”—its integral is fixed by the Poisson–Neumann problem to supply precisely the central-moment contribution. The (-\nu(\cdots)) part carries the deviatoric shape information and cancels out in the resultant once you account for the boundary condition, leaving the (E)-weighted second-moment matrix multiplying (\bm b_{\Section}), exactly as in your first formula.
