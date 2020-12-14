
## Stability analysis



### Exercises

1. Let $X^\ast$ be the set of distinguishable points, i.e. the set of equivalence classes of distance-zero points, in a pseudometric space $X$. Show that $X^\ast$ is the terminal object in a (relatively small) category that also includes $X$ as an object.

2. Define the radius $r(X) = \min\{ D_X(x, X\wo\{x\}) \mid x\in X \}$.
    a) Prove that $r : \mathsf{FinPMet} \to \mathbb{R}_{\geq 0}$ is continuous with respect to a suitable pseudometric on $\mathsf{FinPMet}$.
    b) Rigorously restate the assertion that $\max\{d_X(\minmax(X),x) \mid x\in X\wo \{\minmax(X)\}\}$ varies continuously with perturbations in $X$. Prove that this rigorous statement is correct.

3. Let $\mathcal{D}$ be the set of finite decompositions of $[0,1]$ into intervals of non-negative length, i.e. $0 = a_0 \leq a_1 \leq a_2 \leq \cdots \leq a_n = 1$ for some $n \in \N$.
    a) Prove by construction that every $D\in\mathcal{D}$ is obtained as the (nonincreasing) sequence $m_i=\max\{\min\{d(x,\ell) \mid \ell \in  L_i\} \mid x \in X^\ast \wo L_i\}$ of a maxmin procedure on a finite pseudometric space $X$ seeded with $\ell_0$ so that $\max\{d(x,\ell_0) \mid x \in X^\ast \wo \{\ell_0\}\}=1$. Is it enough to consider only finite metric spaces?
    b) Define a reasonable metric or pseudometric on $\mathcal{D}$. Which is more appropriate?
    c) Prove that the map $(X,\ell_0) \mapsto D$ from part (a) is stable when $\abs{X} \leq 3$.

4. Let $\operatorname{Par}(n)$ denote the set of _partitions_ of $n$, $\operatorname{Par}=\bigsqcup_{n}{\operatorname{Par}(n)}$, and $\operatorname{Par}_{m \times n}$ denote the set of partitions whose Young diagrams are contained in the $m \times n$ rectangle; and let $\operatorname{Comp}(n)$ denote the set of _compositions_ of $n$.
    a) Prove that, for any $x \in X \in \mathsf{FinPMet}$ with $N = \abs{X}$, $Q^\pm(x, X) \in \operatorname{Par}_{N \times N}$.
    b) Furthermore prove that $Q^+(x, X) \in \operatorname{Comp}(N)$.

5. Define a _move_ in $\operatorname{Comp}(n)$ to be of one of the forms
    $$(\ldots, m, 1, n-1, \ldots) \leftrightarrow (\ldots, m, n, \ldots) \leftrightarrow (\ldots, m-1, 1, n, \ldots)$$
Define the "rank-radius" $s(X) = \min\{ Q^+(x, X) \mid x\in X\}$.
    a) Prove that a sufficiently small perturbation of any $x \in X\wo\fl(X)$ produces at most one move difference in $Q^+(x, X)$, and that there exists such a perturbation for each such $x$.
    b) Suppose that $X^\ast$ is in general position. Prove that then a sufficiently small perturbation of any $x \in X$ produces at most one move difference in $s(X)$, and that there exists such a perturbation for each such $x$.

\pagebreak
