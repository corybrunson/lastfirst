---
title: "Fixed and adaptive landmark sets for finite metric spaces"
author:
  - Jason Cory Brunson
  - Yara Skaf
output:
  pdf_document:
    fig_caption: yes
    keep_tex: yes
header-includes:
  \usepackage{format}
---

<!--Submission format
header-includes:
  \usepackage{lineno}
  \linenumbers
  \doublespacing
-->
<!--Print units
header-includes:
  \usepackage{lengthconvert}
  \Convertsetup{unit=cm}
-->

\pagebreak

# Introduction

Topological data analysis (TDA) is a maturing field in data science, at the interface of statistics, computer science, and mathematics.
Topology is the discipline at the interface of geometry (the study of shape) and analysis (the study of continuity) that focuses on geometric properties that are preserved under continuous transformations.
TDA consists in the use of computational theories of continuity to investigate the shape or structure of data.
While TDA is most commonly associated with persistent homology (PH) and mapper, it can be understood to encompass and generalize many conventional and even classical techniques, including cluster analysis, network analysis, and nearest neighbors prediction.

Two basic maneuvers in TDA are locality preservation and cardinality reduction. \emph{Locality preservation} is the property of some functions, such as projections or hashes, that nearby cases or points in the domain have nearby images in the codomain. (Continuity, as defined in analysis and topology, is a type of locality preservation.) This is the defining property of many non-linear dimensionality reduction techniques, including t-SNE and UMAP, but also an asymptotic property of nearest neighbors prediction and the locally-sensitive hashing used in its implementations. We use the term _cardinality reduction_[^cardinality-reduction], in contrast to dimension reduction, to describe techniques that produce more tractable or comprehensible representations of complex data by reducing the number of cases or points rather than the number of variables or dimensions used to represent them. Cardinality reduction describes cluster analysis and association (co-occurrence or correlation) network analysis, for example.

[^cardinality-reduction]: The term "cardinality reduction" takes different meanings in the data analysis literature, including the combining of different values of a categorical variable [@MicciBarreca2001; @Refaat2010] or of time series [@Hu2011] (also "numerosity reduction" [@Lin2007]). Our meaning is that of @ByczkowskaLipinska2017: an $n\times p$ data table of $n$ cases (rows) and $p$ variables (columns) can be dimension-reduced to an $n\times q$ table, where $q<p$, or cardinality-reduced to an $m\times p$ table, where $m<n$. The most common cardinality reduction method is data reduction.

More elaborate TDA techniques combine both maneuvers, often through the use of covering methods, as with the approximation of PH through witness complexes [@deSilva2004] or in the mapper construction [@Singh2007]. Covering methods, in turn, can be enhanced by strategic sampling from large data sets [@], as can other proximity-based techniques like nearest neighbors. The maxmin procedure is often used for this purpose, as it is deterministic, is computationally efficient, and generates more evenly distributed samples than random selection. However, maxmin comes with its own potential limitations in the analysis of data that vary greatly in density or have multiplicities. This is a frequent concern when sparse, heterogeneous, and incomplete data are modeled as finite pseudometric spaces.

In this paper, we develop a landmark sampling procedure complementary to maxmin, based on the rankings of points by a distance or similarity measure rather than on the raw values of such a measure. In the remainder of this section, we motivate the procedure, which we call lastfirst, as a counterpart to maxmin obtained by adapting this procedure from the use of fixed-radius balls to the use of fixed-cardinality neighborhoods. We then formally describe the procedure and prove some of its basic properties in Section\nbs\ref{sec:procedures}. In Section\nbs\ref{sec:implementation} we report the results of benchmark tests and robustness checks of lastfirst on simulated and empirical data sets. We describe some basic and novel applications to real-world data in Section\nbs\ref{sec:experiments}.

## Conventions

$(X, d_X)$ will refer to a finite pseudometric space with point set $X$ and pseudometric $d_X:X\times X\to\R_{\geq 0}$, which by definition satisfies all of the properties of a metric except that $d_X(x,y)=0$ implies $x=y$.
$(X,d_X)$ may be shortened to $X$, and $d_X$ to $d$, when clear from context.
If $x\neq y$ but $d(x,y)=0$ then $x$ and $y$ are said to be indistinguishable or co-located.
The cardinality of $Y\subseteq X$ (counting multiplicities) is denoted $\abs{Y}$, and the set of distinguishable points of $Y$---or of equivalence classes under co-location---is denoted $\supp{Y}$.
When $Y,Z\subseteq X$, let $Y\wo Z$ denote the set of points in $Y$ (with multiplicities) that are distinguishable from all points in $Z$. This means that, when defined, $\min_{y\in Y\wo Z,z\in Z}{d(y, z)}>0$.

We denote the diameter of $Y$ as $D(Y)=\max_{x,y\in Y}{d(x,y)}$, and we write:
\begin{align*}
d(Y,Z) &= \min_{y\in Y,z\in Z}{d(y,z)} & D(Y,Z) &= \max_{y\in Y,z\in Z}{(y,z)} \\
d(x,Y) &= d(\{x\},Y)                   & D(x,Y) &= D(\{x\},Y) \\
\end{align*}
If $d(x,y)=d(x,z)\implies y=z$, then $X$ is considered to be in _locally general position_, even if $d(x,y)=d(z,w)\nimplies \{x,y\}=\{z,w\}$. This is condition implies that $X$ is Hausdorff: $d(x,y)=0\implies x=y$.
$f:X \to Y$ will denote a morphism of pseudometric spaces, meaning a function (or morphism of sets) that preserves locality in the sense that $d_X(x,y)\geq d_Y(f(x),f(y))$.

We use the ball notation $B_{\eps}(x)$ for the set of points less than distance $\eps$ from a point $x$; that is, $B_{\eps}(x) = \{ y \mid d(x,y) < \eps \}$.
We use an overline to also include points exactly distance $\eps$ from $x$: $\cl{B_{\eps}}(x) = \{ y \mid d(x,y) \leq \eps \}$. (While these connote openness and closedness in the discrete topology, every such ball is both open and closed.)
If $\abs{\cl{B_\eps}(x)}\geq k$ and $\eps'<\eps\implies\abs{\cl{B_\eps'}(x)}<k$, then we say that $N^+_k(x)=\cl{B_\eps}(x)$ is the $k$-nearest neighborhood of $x$. When $X$ is in general position, therefore, $\abs{N_k(x)}=k$.

Throughout, let $N=\abs{X}$.
For convenience, we assume $0\in\N$.

## Background

We designed the lastfirst procedure to addresses an issue with the maxmin procedure that arises when, due to the use of binary or categorical variables or to limits on measurement resolution, a data set includes many duplicate or otherwise indistinguishable cases. These render the finite metric space representation of the data non-Hausdorff. While these issues may be negligible when such points are rare, they raise computational and interpretative concerns when they are common. Because our procedure is motivated by the same practical needs as the maxmin procedure, we begin with a discussion of those needs.

The earliest appearance of the maxmin[^maxmin] procedure of which we are aware is by @deSilva2004. The authors propose witness complexes, later generalized to alpha complexes [@], for the rapid approximation of persistent homology: Given a point cloud, a set of landmark points and their overlapping neighborhoods define a nerve, which stands in for the Vietoris--Rips complex at each scale. They use the maxmin procedure, which we define in Section\nbs\ref{sec:maxmin}, as an alternative to selecting landmark points uniformly at random. The procedure ensures that the landmarks are locally separated and roughly evenly distributed.
While the procedure improved little upon uniform random selection in most use cases, on some tasks it far outperformed.

[^maxmin]: This procedure is distinct from many other uses of "maxmin" and related terms.

Subsequent uses of maxmin include the selection of a sample of points from a computationally intractable point cloud for the purpose of downstream topological analysis, as when performing the mapper construction [@Singh2007]; and the optimization of a fixed-radius ball cover of a point cloud, in the sense of minimizing both the number of balls and their shared radius [@Dlotko2019]. In addition to approximating persistent homology [@deSilva2004; @Dlotko2019], maxmin has been used to reduce the sizes of simplicial complex models of point cloud data for the sake of visualization and exploration [@Singh2007; @Dlotko2019].

## Motivation

The ball covers mentioned above have been proposed as an alternative to mapper [@Dlotko2019], where they exchange complexity for computational cost, and the mapper construction itself relies on a crucial covering step that has received limited theoretical attention.
Conventionally, mappers rely on covers consisting either of overlapping intervals (when the lens is one-dimensional) or of their cartesian products (higher-dimensional).
For this purpose, we propose that ball covers, heuristically optimized using the maxmin procedure, have a potential advantage over conventional covers, alongside a potential disadvantage.

Conventionally, mappers use low-dimensional lens spaces $\mathbb{R}^p$ and one of two types of cover, based on overlapping intervals either of fixed length or at fixed quantiles [@Piekenbrock2020].
We think of this length or quantile as the _resolution_ of the cover.
When $m>1$, covers for $Y\subset\mathbb{R}^m$ can be obtained as the cartesian products of those of the coordinate projections of $Y$---so, if $\pi_1,\pi_2:\mathbb{R}^2\to\mathbb{R}$ are the coordinate projections, then interval covers $\{I_\alpha\}$ and $\{J_\beta\}$ of $\pi_1(Y)$ and $\pi_2(Y)$ give rise to the rectangle cover $\{I_\alpha \times J_\beta\}$ for $Y$.
While these cover types are manageable in very few dimensions, the number of sets scales geometrically with $p$, holding the resolution fixed.
Moreover, eventually most of the resulting cover sets will contain no points of $X$, and additional calculations will be needed to restrict to the non-empty sets.

In contrast, a cover obtained by centering balls at a subset of landmark points in $Y$ will have greater up-front computational cost but will be guaranteed to contain no empty sets, and the number of sets required to capture the topology of $Y$ will increase only with the geometric and topological complexity of $Y$, not with $p$. (We test this hypothesis in Section\nbs\ref{sec:experiments}.)

Nevertheless, the maxmin cover requires a meaningful distance metric: the dissimilarity of cases $x$ and $y$ is captured by their distance $d(x,y)$, regardless of where $x$ and $y$ are located in $X$, and the neighborhoods $B_r(x)$ and $B_r(y)$ about landmarks $x$ and $y$ play an equal role in the cover.^[Is there a term for this property, e.g. something being a "universal constant"?]
This means that cover sets centered at landmarks in sparse regions of $X$ will be more numerous and of lower cardinality than those centered in dense regions.
The assumption is violated in many real-world settings, including much of biomedicine. For example, in psychometric terms, this would mean that inter-case distance is an \emph{interval}, not only an \emph{ordinal}, variable, so that the distances between cases in a point cloud representation has a definite meaning independent of which cases are considered. This assumption is often made for convenience, but it generally does not follow from theory.

This motivates us to produce a counterpart to the ball cover that we call the _neighborhood cover_, each set of which may have a different radius but (nearly) the same cardinality. Especially in analyses of medical and healthcare data, underlying variables can often only be understood as ordinal. Other representations of high-dimensional data sets are commonly defined by similarity (or dissimilarity) measures such as cosine similarity rather than by coordinates and associated metrics. Furthermore, because measurements are coarse and often missing, such data often contain indistinguishable entries---cases all of whose measurements are equal and that are therefore represented as multiple instances of the same point. All of these attributes violate the assumptions of the ball cover approach and suggest the need for an ordinal counterpart.

# Procedures
\label{sec:procedures}

In this section we review the maxmin procedure and introduce the lastfirst procedure as a complement to it.

## Maxmin procedure
\label{sec:maxmin}

As defined for the purpose of landmark selection, maxmin takes as input a proper subset $L\subset X$ and returns as output a point $x\in X\wo L$. A related procedure, which we will call minmax, takes as input only $X$ and returns its \emph{Chebyshev center}, the point $x\in X$ for which $D(x,X\wo\{x\})$ is minimized. Whereas maxmin naturally increments a landmark set, minmax naturally initializes one: The Chebyshev center is the landmark that produces the single-set cover using the smallest radius. To help build intuition around these concepts, consider them together with two other procedures: Each procedure calculates either minmax or maxmin, and each does so either with reference to a proper subset of $X$ or with respect to $X$ itself.

Given $L\subset X$ and $x\in X\wo L$, define the \emph{minmax sets}
\begin{align*}
    \minmax(L) &= \{x\in X\wo L\mid D(x,L) = \min_{y\in X\wo L}{D(y,L)}\} \\
    \minmax(X) &= \{x\in X\mid D(x,X\wo\{x\}) = \min_{y\in X}{D(y,X\wo\{y\})}\}
\end{align*}
consisting of \emph{minmax points} and the \emph{maxmin sets}
\begin{align*}
    \maxmin(L) &= \{x\in X\wo L\mid d(x,L) = \max_{y\in X\wo L}{d(y,L)}\} \\
    \maxmin(X) &= \{x\in X\mid d(x,X\wo\{x\}) = \max_{y\in X}{d(y,X\wo\{y\})}\}
\end{align*}
consisting of \emph{maxmin points}.
Note that both $\minmax(\,\cdot\,)$ and $\maxmin(\,\cdot\,)$ are nonempty and that, when $X$ is in locally general position, they both have cardinality $1$.

\begin{algorithm}
\caption{Select a maxmin landmark set.}
\label{alg:maxmin}
\begin{algorithmic}
\REQUIRE finite metric space $(X,d)$
\REQUIRE at least one parameter $\eps\geq 0$ or $n\in\N$
\REQUIRE seed point $\ell_0 \in X$
\IF{only $\eps$ is given}
    \STATE $n \leftarrow 1$
\ENDIF
\IF{only $n$ is given}
    \STATE $\eps \leftarrow \infty$
\ENDIF
\STATE $L \leftarrow \varnothing$
\STATE $i \leftarrow 0$
\REPEAT
    \STATE $L \leftarrow L\cup\{\ell_i\}$
    \STATE $i \leftarrow i+1$
    \STATE $\ell_i \in \maxmin(L)$
    \STATE $d_{\operatorname{max}} \leftarrow d(\ell_i,L)$
\UNTIL $d_{\operatorname{max}} < \eps$ and $\abs{L} \geq n$
\RETURN maxmin landmark set $L$
\end{algorithmic}
\end{algorithm}

The \emph{maxmin procedure} for generating a landmark set $L\subseteq X$ proceeds as follows (see Algorithm\nbs\ref{alg:maxmin}):
First, choose a number $k\leq\uniq{X}$ of landmark points to generate or a radius $\eps\geq 0$ for which to require that the balls $\cl{B_{\eps}}(\ell)$ minimally cover $X$.
Choose a first landmark point $\ell_0\in X$. This choice may be arbitrary; we specifically consider three selection rules: the first point index in the object representing $X$, selection at random, and (random selection from) $\minmax(X)$.
Inductively over $i\in\N$, if ever $i\geq k$ or $d(L,X\wo L)\leq\eps$, then stop.
Otherwise, when $L=\{\ell_0,\ldots,\ell_{i-1}\}$, choose $\ell_i\in\maxmin(L)$, again according to a preferred rule.
If a fixed number $k$ of landmarks was prescribed, then set $\eps=d(L,X\wo L)$; if $\eps$ was prescribed, then set $k=\abs{L}$.

We will write the elements of landmark sets $L=\{\ell_0,\ldots,\ell_{k-1}\}$ in the order in which they were generated.
Note that, if $k=\uniq{X}$ or $\eps=0$, then $L=\supp{X}$.
When the procedure stops, $X=\bigcup_{i=0}^{k-1}{\cl{B_{\eps}}(\ell_i)}$, and this cover is minimal in two senses: The removal of any $\cl{B_\eps}(\ell_i)$ or any decrease in $\eps$ will obtain a collection of sets that fail to cover $X$. More generally, a non-minimal cover can be obtained by specifying both $k$ and $\eps$ in a compatible way. In Section\nbs\ref{sec:implementation}, we describe two adaptive parameters implemented in our software package that make these choices easier.

## Lastfirst procedure

The lastfirst procedure is defined analogously to the maxmin procedure, substituting "rank-distance" for the pseudometric $d_X$.

### Rank-distances

Rank-distance is an adaptive notion of distance with respect to nearest neighborhoods.[^rank-distance] It relies on the underlying pseudometric $d_X$ but takes values in $\N$ given by the ordinal of one point's distance from another.

[^rank-distance]: As used here, rank-distance is distinct from the _rank-distance_ between permutations, which is used to define rank correlation coefficients, and from the _ordinal distance_ proposed by Pattanaik and Xu (2008), a loosening of the concept of pseudometric that dispenses with the triangle inequality.

\begin{definition} (Rank-distance)
    For $x,y\in X$, define the \textit{rank-distance} $q_X : X \times X \longrightarrow \N$ as follows:
    \begin{equation*}
        q_X(x,y)=\abs{\{z\in X\mid d_X(x,z)<d_X(x,y)\}}+1%>
    \end{equation*}
\end{definition}

As with $d$, we allow ourselves to write $q=q_X$ when clear from context. Note that $q(x,x)=1$ and $q(x,y)\leq N$, and that $q$ is not, in general, symmetric. (It is therefore not a pseudometric.)^[These images might be made larger, with less or no transparency.]

\begin{example}\label{ex:rank-distance}
    Consider the simple case $X = \{a, b, c, d\}$, visualized below, equipped with the standard Euclidean metric:
    \begin{centeredTikz}
        [every label/.append style={text=black!60!blue, font=\scriptsize}]
        \draw[gray] (0,0) -- (5,0);

        \foreach \i in {0,...,3}
            \draw[gray] (\i,0.1) -- + (0,-0.25) node[font=\scriptsize, text=gray, below] {$\i$};
        \draw[gray] (3,0.1) -- + (0,-0.25) node[font=\scriptsize, text=gray, below] {};
        \draw[gray] (5,0.1) -- + (0,-0.25) node[font=\scriptsize, text=gray, below] {$5$};

        \node[circle, draw=blue!60, fill=blue!5, inner sep=0.5mm, label=above:{$a = 1$}] at (1,0) {};
        \node[circle, draw=blue!60, fill=blue!5, inner sep=0.5mm, label=above:{$b = 2$}] at (2,0) {};
        \node[circle, draw=blue!60, fill=blue!5, inner sep=0.5mm, label=above:{$c = 4$}] at (4,0.1) {};
        \node[circle, draw=blue!60, fill=blue!5, inner sep=0.5mm, label=below:{$d = 4$}] at (4,-0.1) {};
    \end{centeredTikz}

    Lack of symmetry of $q$ is shown by points $b$ and $c$ :
    \begin{align*}
        q(a,c) &= \abs{\{x \in X \mid \abs{x-a} < \abs{c-a} = 2\}} + 1 &
        q(c,a) &= \abs{\{x \in X \mid \abs{x-c} < \abs{a-c} = 2\}} + 1 \\
               &= \abs{\{a, b\}} + 1 &&= \abs{\{b, c, d\}} + 1 \\
               &= 3 &&= 4
    \end{align*}

    Observe in particular that $\max_{x\in X}{q(a,x)}=3<\abs{X}$; $q(x,\,\cdot\,)$ will not max out at $N$ when the most distant points from $x$ have multiplicity.

    However, $q(x,x) = 1$ whether $x$ has multiplicity or not:
    \begin{align*}
        q(a,a) &= \abs{\{x \in X \mid \abs{x-a} < \abs{a-a} = 0\}} + 1 &
        q(c,c) &= \abs{\{x \in X \mid \abs{x-c} < \abs{c-c} = 0\}} + 1 \\
               &= \abs{\varnothing} + 1 &&= \abs{\varnothing} + 1 \\
               &= 1 &&= 1
    \end{align*}
\end{example}

We term the unary rankings $q(x,\,\cdot\,)$ and $q(\,\cdot\,,x)$ the \emph{out- (from $x$)} and \emph{in- (to $x$) rankings} of $X$, respectively. These can then be used to define \textit{out-} and \textit{in-neighborhoods} of $x$.[^out-in]

[^out-in]: The terminology and notation are adapted from the theory of directed graphs. These definitions are the same as those for a complete directed graph on $X$ with directed arcs $x\to y$ weighted by rank-distance $q(x,y)$.

\begin{definition} ($k$-neighborhoods)
    For $x \in X$, define the \emph{$k$-out-neighborhoods} $N^+_k$ and \emph{$k$-in-neighborhoods} $N^-_k$ of $x$ as the sets
    \begin{align*}
        & N^+_k(x)=\{y\in X\mid q(x,y)\leq k\} \\
        & N^-_k(x)=\{y\in X\mid q(y,x)\leq k\}
    \end{align*}
\end{definition}

Note that $\varnothing \subseteq N^\pm_1(x) \subseteq \cdots \subseteq  N^\pm_N(x) = X$.
The $k$-out-neighborhoods of $x$ are the sets of points in $X$ that have rank-distance at most $k$ from $x$. This is equivalent to the $k$-nearest neighbors of $x$. The $k$-in-neighborhoods of $x$ are the sets of points in $X$ from which $x$ has rank-distance at most $k$.
These definitions can be adapted as follows to be relative to a subset $Y \sub X$, using the notation $q_X$ to emphasize that the points in $X\wo (Y\cup\{x\})$ are still involved in the calculation:

\begin{align*}
    & N^+_k(x,Y)=\{y\in Y\mid q_X(x,y)\leq k\} \\
    & N^-_k(x,Y)=\{y\in Y\mid q_X(y,x)\leq k\}
\end{align*}

\begin{example}\label{ex:rank-neighborhoods}
Consider the same $X$ as in Example\nbs\ref{ex:rank-distance}. Compute $N_k^+$ and $N_k^-$ for $a$ and $c$, using $k = 3$:
\begin{align*}
    N_3^+ (a) &= \{x \in X \mid q(a,x) \leq 3\} &
    N_3^+ (c) &= \{x \in X \mid q(c,x) \leq 3\} \\
    &= \{a, b, c, d\} &
    &= \{b, c, d\} \\
    \\
    N_3^- (a) &= \{x \in X \mid q(x,a) \leq 3\} &
    N_3^- (c) &= \{x \in X \mid q(x,c) \leq 3\} \\
    &= \{a, b\} &
    &= \{a, b, c, d\}
\end{align*}
\end{example}

Rank-distances are not as straightforward to compare among subsets of points. For example, for $y\neq x$, $q(x,y)$ takes integer values between $\cl{B_0}(x)+1$ and $N$. To explain and motivate their use, we now intuitively adapt the minmax and maxmin procedures to this setting, then state and prove a formal definition for their analogs.

### Rank sequences and landmark selection

Consider the case of choosing the next landmark point given a subset $L=\{\ell_0,\ldots,\ell_{i-1}\}\subset X$ of collected landmark points. For a ball cover, we choose $\ell_i\in X\wo L$ so that the \emph{minimum radius $\eps$} required for some $B_\eps(\ell_j)$ to contain $\ell_i$ is \emph{maximized}.
If $X$ is in general position, the choice of $\ell_i$ will be unique; otherwise, $\ell_i$ can be chosen at random from the subset of $X$ that satisfies this criterion.
Analogously, for a neighborhood cover, we want that the \emph{minimum cardinality $k$} required for some $N^+_k(\ell_j)$ to contain $\ell_i$ is \emph{maximized}. Switching perspective from out- to in-, reversing the roles of $L$ and $\ell_i$, we want $N^-_k(\ell_i,L)=1$ for the latest (largest) $k$ possible, say $k^-_1$. When indistinguishable points abound, this may still not uniquely determine $\ell_i$, so we may extend the principle: Among those $\ell\in X\wo L$ for which $N^-_{k^-_1}(\ell_i,L)=1$, choose $\ell_i$ for which $N^-_{k}(\ell_i,L)=2$ for the latest $k$ possible, say $k^-_2\geq k^-_1$. Continue this process until only one candidate $\ell$ remains (up to multiplicity), or until $N^-_{k}(\ell,L)=\abs{L}$, in which case we consider all remaining candidates equivalent.

Now consider the case of choosing an initial landmark point with respect to a subset $L\subset X$. For a ball cover, we above considered the Chebyshev center $\ell_0\in X$ so that the \emph{maximum radius $\eps$} required for $\cl{B_\eps}(\ell_0)$ to contain any $\ell\in L$ is \emph{minimized}.
This means that a one-set ball cover of $L$ centered at $\ell_0$ would have minimum radius.
Analogously, for a neighborhood cover, we would want the \emph{maximum cardinality $k$} required for $N^+_k(\ell_0)$ to contain any $\ell\in L$ to be \emph{minimized}. That is, we would want $N^+_k(\ell_0,L)=\abs{L}$ for the earliest (smallest) $k$ possible, say $k^+_{\abs{L}}$. To decide among several points $\ell$ with this property, we would then choose one for which $N^+_k(\ell,X)=\abs{L}-1$ for the earliest $k$ possible, say $k^+_{\abs{L}-1}\leq k^+_{\abs{L}}$. This process would likewise continue until only one candidate $\ell$ remains, or until $N^+_{k}(\ell,L)=0$ and all remaining candidates are considered equivalent.

We formalize these two procedures, and companion procedures without reference to $L$, by defining two sequences of neighborhood sizes that encode the optimized cardinalities $k_i$ and imposing on them two total orders that encode the preferences between them.

\begin{definition}\label{def:rank-sequence} (Rank sequences)
    For $x \in X$ and $L\subset X$, define the \emph{out-rank} ($Q^+$) and \emph{in-rank} ($Q^-$) \emph{sequences} of $k$-neighborhood cardinalities:
    \begin{align*}
        Q^\pm(x) &= (\abs{N^\pm_k(x)})_{k=1}^N &
        Q^\pm(x,L) &= (\abs{N^\pm_k(x,L)})_{k=1}^{\abs{L}}
    \end{align*}
\end{definition}

\begin{proposition}
Taking the $k^\pm_i$ as defined above,
\begin{align*}
    Q^+(x,L) &= (k^+_1-1,k^+_2-1,\ldots,k^+_{\abs{L}}-1,\abs{L}) \\
    Q^-(x,L) &= (1^{k^-_1},2^{k^-_2-k^-_1},\ldots,{(\abs{L}-1)}^{k^-_{\abs{L}-1}-k^-_{\abs{L}-2}},{\abs{L}}^{\abs{L}-k^-_{\abs{L}-1}})
\end{align*}
\end{proposition}

\begin{lemma}\label{lem:out-rank}
If $Q = Q^+(x,L)$, then $Q_i \geq i$.
\end{lemma}

\begin{proof}
For convenience, assume that $d(x,\ell_0) \leq d(x,\ell_1) \leq \cdots \leq d(x,\ell_{k-1})$, and consider $i>1$.
If $d(x,\ell) \leq d(x,\ell_{i-1})$, then $q(x,\ell) = \abs{\{ \ell' \in L \mid d(x,\ell') < d(x,\ell) \}} + 1 \leq \abs{\{ \ell' \in L \mid d(x,\ell') < d(x,\ell_{i-1}) \}} + 1 \leq \abs{\{ \ell_0, \ldots, \ell_{i-2} \}} + 1 = i$;
whereas, if $d(x,\ell) > d(x,\ell_{i-1})$, then $q(x,\ell) = \abs{\{ \ell' \in L \mid d(x,\ell') < d(x,\ell) \}} + 1 > \abs{\{ \ell' \in L \mid d(x,\ell') < d(x,\ell) \} \wo \{ \ell_{i-1} \}} + 1 \geq \abs{\{ \ell' \in L \mid d(x,\ell') < d(x,\ell_{i-1}) \}} + 1 = \abs{\{ \ell_0, \ldots, \ell_{i-2} \}} + 1 = i$.
(The reader may check that $q(x,\ell) \leq 1$ if and only if $d(x,\ell) \leq d(x,\ell_0)$.)
In the first instance, then, $Q_i = \abs{N^+_i(x,L)} = \abs{\{\ell \in L \mid q(x,\ell) \leq i\}} = \abs{\{\ell \in L \mid d(x,\ell) \leq d(x,\ell_{i-1})\}} \geq \abs{\{ \ell_0, \ldots, \ell_{i-1} \}} = i$.
\end{proof}

\begin{corollary}\label{cor:out-rank}
Given $L=\{ \ell_0, \ldots, \ell_{k-1} \}$, the maximum possible out-rank sequence $Q^+(x,L)$ is $Q=(1,2,\ldots,k)$.
\end{corollary}

An example of this computation is shown in Example\nbs\ref{ex:rank-sequence}.

\begin{example}\label{ex:rank-sequence}
Continuing Example\nbs\ref{ex:rank-neighborhoods}, we can compute the other $N_k^+$ and $N_k^-$ to obtain $Q^+$ and $Q^-$ for $a$ and $c$:
\begin{align*}
    Q^+ (a) &= (\abs{N^+_1(a)}, \abs{N^+_2(a)}, \abs{N^+_3(a)}, \abs{N^+_4(a)}) &
    Q^+ (c) &= (\abs{N^+_1(c)}, \abs{N^+_2(c)}, \abs{N^+_3(c)}, \abs{N^+_4(c)}) \\
    &= (\abs{\{a\}}, \abs{\{a, b\}}, \abs{\{a, b, c, d\}}, \abs{\{a, b, c, d\}}) &
    &= (\abs{\{c, d\}}, \abs{\{c, d\}}, \abs{\{b, c, d\}}, \abs{\{a, b, c, d\}}) \\
    &= (1, 2, 4, 4) &
    &= (2, 2, 3, 4) \\
    \\
    Q^- (a) &= (\abs{N^-_1(a)}, \abs{N^-_2(a)}, \abs{N^-_3(a)}, \abs{N^-_4(a)}) &
    Q^- (c) &= (\abs{N^-_1(c)}, \abs{N^-_2(c)}, \abs{N^-_3(c)}, \abs{N^-_4(c)}) \\
    &= (\abs{\{a\}}, \abs{\{a, b\}}, \abs{\{a, b\}}, \abs{\{a, b, c, d\}}) &
    &= (\abs{\{c, d\}}, \abs{\{c, d\}}, \abs{\{a, b, c, d\}}, \abs{\{a, b, c, d\}}) \\
    &= (1, 2, 2, 4) &
    &= (2, 2, 4, 4)
\end{align*}
\end{example}

\begin{definition} (Total orders on rank sequences)
    Let $a_n = (a_1,\ldots,a_M),b_n = (b_1,\ldots,b_M)\in\N^M$. Then
    \begin{itemize}
        \item $a_n < b_n$ in the reverse colexicographic (revcolex) order if $\exists i, a_i > b_i \wedge (\forall j>i, a_j = b_j)$.
        \item $a_n < b_n$ in the reverse lexicographic (revlex) order if $\exists i, a_i > b_i \wedge (\forall j<i, a_j = b_j)$.
    \end{itemize}
\end{definition}

Impose the revcolex order on the $Q^+$ to emphasize later out-distances. Impose the revlex order on the $Q^-$ to emphasize earlier in-distances. In both cases, sequences with more large values indicate points with lower rank-distances to or from more other points, which is why both orders are reversed.

We now define counterparts to the minmax and maxmin procedures using these totally ordered sequences.

\begin{definition}(Firstlast and lastfirst)
    Given $L\subset X$, define the \emph{firstlast sets}
    \begin{align*}
        \fl(L) &= \{x\in X\wo L\mid Q^+(x,L) = \min_{y\in X\wo L}{Q^+(y,L)}\} \\
        \fl(X) &= \{x\in X\mid Q^+(x,X\wo\{x\}) = \min_{y\in X}{Q^+(y,X\wo\{y\})}
    \end{align*}
    consisting of \emph{firstlast points} and the \emph{lastfirst sets}
    \begin{align*}
        \lf(L) &= \{x\in X\wo L\mid Q^-(x,L) = \max_{y\in X\wo L}{Q^-(y,L)}\} \\
        \lf(X) &= \{x\in X\mid Q^-(x,X\wo\{x\}) = \max_{y\in X}{Q^-(y,X\wo\{y\}}\}
    \end{align*}
    consisting of \emph{lastfirst points}.
\end{definition}

\begin{remark}
Replacing each $Q^\pm(x,X\wo\{x\})$ with $Q^\pm(x)$ would simply increase the integer in each position by $1$, so the latter could be used equivalently to the former.
\end{remark}

\begin{example}
    Return again to $X=\{a,b,c,d\}$ from Example\nbs\ref{ex:rank-distance}. We calculate an exhaustive lastfirst landmark set, seeded with a firstlast point:
    \begin{enumerate}
        \item We have
        \begin{align*}
            Q^+(a,\{b,c,d\}) &= (0,1,3,3) &
            Q^+(b,\{a,c,d\}) &= (0,1,3,3) \\
            Q^+(c,\{a,b,d\}) &= (1,1,2,3) &
            Q^+(d,\{a,b,c\}) &= (1,1,2,3)
        \end{align*}
        Under the revcolex order, $\min_{x\in X}{Q^+(x,X\wo\{x\})}=(0,1,3,3)$, and we randomly select $\ell_0=a$.
        \item We now have
        \begin{align*}
            Q^-(b,\{a\}) &= (0,1,1,1) &
            Q^-(c,\{a\}) &= (0,0,1,1) &
            Q^-(d,\{a\}) &= (0,0,1,1)
        \end{align*}
        Under the revlex order, $\max_{x\in X\wo\{a\}}{Q^-(x,\{a\})}=(0,0,1,1)$, and without loss of generality, we select $\ell_1=c$, so that now $L=\{a,c\}$.
        \item Only one point remains in $X\wo L=\{b\}$, so the exhaustive landmark set is $\{a,c,b\}$.
    \end{enumerate}
\end{example}

### Algorithms

For illustration, we provide a complete algorithm to identify the firstlast set of $X$ and a simplified algorithm to obtain a lastfirst landmark set using specified parameters. The former requires first calculating $Q^+(x,L)$ from $x$ and $L$, which is done in Algorithm\nbs\ref{alg:outnbhdseq}.
This algorithm uses the notation $[a,b]$ as shorthand for the arithmetic sequence $(a,a+1,\ldots,b)$ of indices.
Making use of these out-rank sequences, the recovery of a firstlast set proceeds as in Algorithm\nbs\ref{alg:firstlast}.

\begin{algorithm}
\caption{Compute the out-rank sequence of a point with respect to a proper subset.}
\label{alg:outnbhdseq}
\begin{algorithmic}
\REQUIRE finite pseudometric space $(X,d)$
\REQUIRE landmark set $L=\{\ell_0,\ldots,\ell_{k-1}\}\subset X$
\REQUIRE point $x\in X$
\STATE $D \leftarrow \verb|sort|(d(x,\ell_0),\ldots,d(x,\ell_{k-1}))=(d_1\leq\cdots\leq d_k)$
\STATE $\tilde{D} \leftarrow (\tilde{d}_0=-\infty, \tilde{d}_1=d_1, \ldots, \tilde{d}_k=d_k, \tilde{d}_{k+1}=+\infty)$
\STATE $\tilde{Q} \leftarrow (0,\ldots,0)\in\N^{[0,k+1]}$
\STATE $j \leftarrow 0$
\FOR{$i=1$ to $k+1$}
    \IF{$\tilde{d}_i > \tilde{d}_{i-1}$}
        \WHILE{$j<i$}
            \STATE $j \leftarrow j+1$
            \STATE $\tilde{Q}_j \leftarrow \tilde{Q}_{j-1}$
        \ENDWHILE
    \ENDIF
    \STATE $\tilde{Q}_j \leftarrow \tilde{Q}_j+1$
\ENDFOR
\STATE $Q \leftarrow \tilde{Q}_{[1,k]}$
\RETURN out-rank sequence $Q$
\end{algorithmic}
\end{algorithm}

\begin{lemma}
Algorithm\nbs\ref{alg:outnbhdseq} returns the out-rank sequence.
\end{lemma}

\begin{proof}
The algorithm takes $k+1$ steps, which we consider in three groups.
At any step $i$, $j$ never decreases and increases at most to $j=i$, so the value of $\tilde{Q}_{j'}$ is fixed once $j>j'$.

\begin{itemize}
\item[$i = 1$.]
In this step, $\tilde{d}_i = \tilde{d}_1 \geq 0 > -\infty = \tilde{d}_0 = \tilde{d}_{i-1}$, so the \verb|while| loop is run.
The loop begins with $j=0$ and ends with $j=1$, so it has only the vacuous effect $\tilde{Q}_1 \leftarrow \tilde{Q}_0$ on $\tilde{Q}$.
The only effect in this case is the final increment $\tilde{Q}_1 \leftarrow \tilde{Q}_1 + 1$, so that the step ends with $\tilde{Q} = (0,1,0,\ldots,0)$.
\item[$1 < i \leq k$.]
The intermediate steps proceed in two ways, depending on whether $d_i = \tilde{d}_i > \tilde{d}_{i-1} = d_{i-1}$---that is, whether the next-closest point in $L$ is farther from $x$ than the previous.
For convenience, let $i' \leq i-1$ denote the previous value of $i$ at which the distance from $x$ increased.
\begin{itemize}
\item[$d_i = d_{i-1}$.]
In this case, the only change to $\tilde{Q}$ is to increase $\tilde{Q}_{i'}$ by $1$.
Thus, at the end of the step, $\tilde{Q}$ equals the number of members of $L$ up to index $i$ at distance at most $d_{i'} = d_i$ from $x$.
\item[$d_i > d_{i-1}$.]
In this case, $j$ starts out at $i'$.
If $i' < i-1$, then the \verb|while| loop copies the value $\tilde{Q}_{i'}$ to the values of $\tilde{Q}_{[i'+1,i-1]}$.
By the same reasoning as in the previous case plus the fact that $d_i > d_{i-1}$, $\tilde{Q}_{i'}$ at this step equals $\abs{N^+_{i'}(x,L)} = \abs{N^+_{i-1}(x,L)}$.
Therefore, at the end of this step, $\tilde{Q}_{[i',i-1]}=(Q^+(x,L))_{[i',i-1]}$ and $j>i-1$, so that the returned values $Q_{[i',i-1]}$ will be correct.
Regardless, this value is also copied to $\tilde{Q}_i$, and $\tilde{Q}_i$ is incremented by $1$, so that $\tilde{Q}_i = \tilde{Q}_{i-1} + 1$.
\end{itemize}
Regardless of the case, at each step $\max\tilde{Q}$ increases by $1$, so that after step $i = k$ we have $\max\tilde{Q} = k$.
The last step wraps up the process in the case that $d_k = d_{k-1}$.
\item[$i = k + 1$.]
In this step, $\tilde{d}_i = \tilde{d}_{k+1} = +\infty > d_k = \tilde{d}_{i-1}$, so the \verb|while| loop is run.
If $d_k > d_{k-1}$, then this has no effect on $Q=\tilde{Q}_{[1,k]}$.
If, instead, $d_k = d_{k-1}$, then the effect is to copy $\tilde{Q}_{i'}$ to $\tilde{Q}_{[i'+1,k]}$.
As in the previous step, this makes $\tilde{Q}_{[i'+1,k]}=(Q^+(x,L))_{[i',k]}$, which ensures that the returned $Q_{[i',k]}$ are correct.
\end{itemize}
\end{proof}

\begin{algorithm}
\caption{Identify a firstlast set with respect to a proper subset.}
\label{alg:firstlast}
\begin{algorithmic}
\REQUIRE finite metric space $(X,d)$
\REQUIRE subset $L=\{\ell_0,\ldots,\ell_{k-1}\}\subset X$
\STATE $F \leftarrow \varnothing$
\STATE $Q \leftarrow (k,\ldots,k)\in\N^{k}$
\FOR{$x\in X$}
    \STATE $Q' \leftarrow Q^+(x,L)\in\N^{k}$ (Algorithm\nbs\ref{alg:outnbhdseq})
    \STATE $j \leftarrow k$
    \WHILE{$Q'_j = Q_j$}
        \STATE $j \leftarrow j-1$
    \ENDWHILE
    \IF{$j = 0$}
        \STATE $F \leftarrow F\cup\{x\}$
    \ELSE
        \IF{$Q'_j < Q_j$}
            \STATE $F \leftarrow \{x\}$
            \STATE $Q \leftarrow Q'$
        \ENDIF
    \ENDIF
\ENDFOR
\RETURN firstlast set $F$
\end{algorithmic}
\end{algorithm}

\begin{proposition}
Algorithm\nbs\ref{alg:firstlast} returns the firstlast set.
\end{proposition}

\begin{proof}
In the algorithm, $Q$ begins at what we know from Corollary\nbs\ref{cor:out-rank} to be the maximum possible out-rank sequence.
The algorithm then simply computes the out-rank sequence of every point $x \in X$ and either stores it in a list, if it is equal to the current minimum, or resets and replaces the stored list, if it is smaller.
Once all points have been considered, those in the stored list will constitute the firstlast set.
\end{proof}

The implementation of $\lf$ is more involved: Because an in-rank sequence $Q^-(x)$ cannot be calculated from the distances $\{d(x,y)\mid y\in X\}$ alone, much more computation is required to generate $Q^-(x,L)$ at each step.
A partially vectorized R implementation uses the combinatorial identity in Lemma\nbs\ref{lemma-revlex} to expedite this calculation.

\begin{lemma}\label{lemma-revlex}
Let $S$ be a set of nondecreasing sequences of length $N$ taking integer values in $[n]$.
Given $q = (q_1 \leq \cdots \leq q_N) \in S$, define $f(q) = ( \abs{\{ i \mid q_i = k \}} )_{k=1}^{n}$, and let $T = \{ f(q) \mid q \in S \}$.
Then the function $f$ is injective.
Impose the lexicographic (lex) order on $S$ and the revlex order on $T$.
Then, furthermore, for $p,q \in S$, $p < q$ if and only if $f(p) < f(q)$.
\end{lemma}

\begin{proof}
That $f$ is injective is shown via an inverse construction.
\textbf{This is basic and should be a citable fact.}
Given $t = (t_1, \ldots, t_n) \in T$, define $g(t) = (\, 1^{t_1}\, \cdots \, n^{t_n})$, where $a^b$ represents a constant sequence of value $a$ and length $b$.
By definition of $f$, $t_k$ is the number of $k$s in $q$; and, by the definition of $g$, $t_k$ is the number of $k$s in $g(t)$.
Since both $q$ and $g(t)$ are in ascending order and take values in $[n]$, it must be that $g(t) = q$.

Now suppose $p < q$ in $S$, which means that there exists $j \in [N]$ for which $p_j < q_j$ and if $i < j$ then $p_i = q_i$.
Take $k = p_j$, $s = f(p)$, and $t = f(q)$.
Then $s_i = t_i$ for all $i < k$, since these counts are determined by indices of $p$ and $q$ before $j$. Furthermore, $s_k > t_k$, since $p$ takes the value $k$ at least once among the remaining indices (namely $j$) while $q$ does not take this value among them.
Thus $s < t$ in the revlex order.
\end{proof}

The construction of a lastfirst landmark set in R proceeds as follows, as outlined in Algorithm\nbs\ref{alg:lastfirst-landmarks}.
Let $L$ denote the set of landmark points. At the start of the algorithm, a seed point $\ell_0$ is selected so that $L = \{\ell_0\}$.
At each iteration, $\lf(L)$ is computed and the next landmark point $\ell_i$ is selected from this set. The procedure terminates when $X=\bigcup_{j\leq i}N^+_1(\ell_j)$, hence $\supp{L}=\supp{X}$.

\begin{algorithm}
\caption{Calculate the lastfirst landmark sequence from a seed point.}
\label{alg:lastfirst-landmarks}
\begin{algorithmic}
\REQUIRE finite pseudometric space $(X,d)$
\REQUIRE at least one parameter $k>0$ (target cover set cardinality) or $n\in\N$ (number of landmarks)
\REQUIRE seed point $\ell_0 \in X$
\REQUIRE selection procedure \verb|pick|
\STATE $L \leftarrow \{ \ell_0 \}$ landmark set
\STATE $F \leftarrow L$ initial lastfirst set
\STATE $R \in \N^{N \times 0}$, a $0$-dimensional $\N$-valued matrix
\FOR{$i$ from $1$ to $\uniq{X}$}
    \STATE $\ell_i \leftarrow \verb|pick|(F)$
    \STATE $L \leftarrow L \cup \{\ell_i\}$
    \STATE $D_i \leftarrow (d_{i1},\ldots,d_{iN}) \in {\R_{\geq 0}}^N$, where $d_{ij} = d(\ell_i, x_j)$
    \STATE $Q_i \leftarrow \verb|rank|(D_i) \in {\N_{\geq 0}}^N$, so that $Q = (q_{i1},\ldots,q_{iN})$, where $q_{ij} = q(\ell_i, x_j)$
    \STATE $R \leftarrow [R, Q_i] \in \N^{N \times i}$
    \STATE $k_{\min} \leftarrow \max\{ \min_{j=1,i}{ R_{ij} } \}_{1 \leq i \leq N}$ (the minimum cardinality required for sets centered at $L$ to cover $X$)
    \IF{$D(L, X \wo L) = 0$}
        \STATE \textbf{break}
    \ENDIF
    \IF{$i \geq n$ and $k_{\min} \leq k$}
        \STATE \textbf{break}
    \ENDIF
    \STATE $R \leftarrow [ \verb|sort|({R_{1,\bullet}}^\top) \cdots \verb|sort|({R_{N,\bullet}}^\top) ]^\top \in \N^{N \times i}$
    \STATE $F \leftarrow X \wo \uniq{L}$
    \FOR{$j$ from $1$ to $i$}
        \STATE $F \leftarrow \{x_i \in F \mid R_{ij} = \max_{i'}{R_{i'j}}\}$
    \ENDFOR
\ENDFOR
\RETURN lastfirst landmark set $L$ with at least $n$ cover sets of cardinality at most $k$
\end{algorithmic}
\end{algorithm}

\begin{proposition}
Algorithm\nbs\ref{alg:lastfirst-landmarks} returns a lastfirst landmark set.
If $n$ is given as input and $k$ is not, $\abs{L} = n$. If $n$ and $k$ are both given, $\abs{L} \geq n$. Otherwise $L$ is minimal, meaning that no proper prefix (subset?) of $L$ gives a cover of $X$ by $k$-nearest neighborhoods.
\end{proposition}

\begin{proof}

\end{proof}

\textbf{Old versions below.}

\begin{proposition}
Let $L$ denote the landmark set returned by Algorithm\nbs\ref{alg:lastfirst-landmarks}. If $n$ is given as input and $k$ is not, $|L| = n$. If $n$ and $k$ are both given, $|L| \geq n$. Otherwise $L$ is minimal, meaning that no proper subset of $L$ gives a cover of $X$ by $k$-nearest neighborhoods.
\end{proposition}

\begin{proof}
Let $(X,d)$ be a finite metric space and $\ell_0 \in X$ be a seed point as required by Algorithm\nbs\ref{alg:lastfirst-landmarks}.

Recall that for the algorithm to terminate its loop and subsequently return the resulting landmark set $L$, both of the following exit conditions must hold: (1) $q_{\max} \leq k$ and (2) $|L| \geq n$.

Suppose first that $n$ is given.
Regardless of whether $k$ is provided, (2) must always hold for the algorithm to terminate, so $|L| \geq n$ for any $k$.
If $k$ is not given, it is set to $|X|$ before the loop, which means (1) holds from the first iteration onward since $|X|$ is by definition the maximum value $q$ can attain.
Then the algorithm terminates as soon as condition (2) is first met, which is when $|L| = n$.
\footnote{Note that $|L| > n$ if and only if $q_{\max} \leq k$ is not satisfied when $|L| = n$, meaning $X$ would not be covered by $k$-neighborhoods around $n$ landmark points, so more landmarks must be chosen to guarantee the algorithm produces a valid cover by neighborhoods of size $k$.}

Now suppose $n$ is not given. Then $k$ must be given since the algorithm requires at least one of the two parameters $n,k$ to be provided.
Therefore, $n$ is set to $0$ before the loop, meaning (2) $|L| \geq n = 0$ always holds, so the algorithm returns $L$ as soon as (1) is satisfied, i.e. when % $q_{\max} := \min_{\ell \in L} q(\ell, \ell_i) \leq k$ for any $\ell_i \in \mathrm{lf}(L)$.
\[
    q_{\max} := \min_{\ell \in L} q(\ell, \ell_i) \leq k \quad\A \ell_i \in \mathrm{lf}(L)
\]
This means the point $\ell_i \in \mathrm{lf}(L)$ that is farthest from any $\ell_j \in L$ by $q$ is still within a $k$-neighborhood of some landmark $\ell_j$.
Since this occurs and causes the loop to exit as soon as the space $X$ can be covered by $k$-neighborhoods of points in $L$, $|L|$ is as small as possible.
Further, $q_{\max} > k$ at any previous iteration before (1) is met, meaning there is at least one point in $X \wo L$ that would \textit{not} be covered by a $k$-neighborhood of any landmark in $L$. This smaller set $L$ is therefore insufficient cover the space, so no proper subset of $L$ can yield a cover of $X$ by $k$-nearest neighborhoods.

\end{proof}





### Tie handling

We might have defined two \textit{rank-distances} $\check{q}, \hat{q} : X \times X \longrightarrow \N$ ("$q$-check" and "$q$-hat") as follows:
\begin{align*}
& \check{q}(x,y)=\abs{\{z\in X\mid d(x,z)<d(x,y)\}}+1 \\%>
& \hat{q}(x,y)=\abs{\{z\in X\mid d(x,z)\leq d(x,y)\}}
\end{align*}
In this notation, $\check{q}=q$, while $\hat{q}(x,y)$ is the cardinality of the smallest ball centered at $x$ that contains $y$.
Letting $\dot{q}$ stand in for either $\check{q}$ or $\hat{q}$, we can define $\dot{N}^\pm_k(x)$ and $\dot{Q}^\pm(x)$ as before and arrive at corresponding notions of firstlast and lastfirst sets.

Note that then $\check{N}^\pm_1(x) \subseteq  \{x\} \subseteq \hat{N}^\pm_1(x)$, and that $\hat{q}(x,x)>1$ when $x$ has multiplicity.
These two rank-distances derive from two tie-handling schemes for calculating rankings of lists with duplicates: For example, if $a<b=c<d$ are the distances from $x$ to $y_1,y_2,y_3,y_4$, respectively, then $(\check{q}(x,y_1),\check{q}(x,y_2),\check{q}(x,y_3),\check{q}(x,y_4))=(1,2,2,4)$ and $(\hat{q}(x,y_1),\hat{q}(x,y_2),\hat{q}(x,y_3),\hat{q}(x,y_4))=(1,3,3,4)$.

Conceptually, these tools would produce landmark sets that yield neighborhood covers with smaller, rather than larger, neighborhoods in regions of high multiplicity.
While we do not use these ideas in this study, they may be suitable in some settings or for some purposes, for example when high multiplicity indicates a failure to distinguish similar but distinct cases that the analyst wishes to maintain greater separation between.
They may also play a role in stability analyses, e.g. by producing an interleaving sequence of nerves of covers.

\begin{example}\label{ex:rank-distance-max}

Recall $X=\{a,b,c,d\}$ from Example\nbs\ref{ex:rank-distance}. The rank-distance $\hat{q}$ is also asymmetric:
\begin{align*}
    \hat{q}(b,c) &= \abs{\{x \in X \mid \abs{x-b} \leq \abs{c-b} = 2\}} &
    \hat{q}(c,b) &= \abs{\{x \in X \mid \abs{x-c} \leq \abs{b-c} = 2\}} \\
        &= \abs{\{a, b, c, d\}} &&= \abs{\{b, c, d\}} \\
        &= 4 &&= 3
\end{align*}

Observe that $\hat{q}(x,x) = 1$ only for distinguishable points $x \in X$:
\begin{align*}
    \hat{q}(a,a) &= \abs{\{x \in X \mid \abs{x-a} \leq \abs{a-a} = 0\}} &
    \hat{q}(c,c) &= \abs{\{x \in X \mid \abs{x-c} \leq \abs{c-c} = 0\}} \\
       &= \abs{\{a\}} &&= \abs{\{c, d\}} \\
       &= 1 &&= 2
\end{align*}

Finally, observe that $\hat{q}(x,\,\cdot\,)$ always maxes out at $\abs{X}$: $\hat{q}(a,c) = \hat{q}(b,c) = \hat{q}(c,a) = \hat{q}(d,a) = \abs{X}$.

Continuing on as in Example\nbs\ref{ex:rank-neighborhoods}, we can compute $\hat{N}_2^+$ and $\hat{N}_2^-$ for $b$ and $c$:
\begin{align*}
    \hat{N}_2^+ (b) &= \{x \in X \mid \hat{q}(b,x) \leq 2\} &
    \hat{N}_2^+ (c) &= \{x \in X \mid \hat{q}(c,x) \leq 2\} \\
        &= \{a, b\} &&= \{c, d\} \\
        \\
    \hat{N}_2^- (b) &= \{x \in X \mid \hat{q}(x,b) \leq 2\} &
    \hat{N}_2^- (c) &= \{x \in X \mid \hat{q}(x,c) \leq 2\} \\
        &= \{a, b\} &&= \{c, d\}
\end{align*}

Similarly, we can compute the other $\hat{N}_k^+$ and $\hat{N}_k^-$ to obtain $\hat{Q}^+$ and $\hat{Q}^-$ for $b$ and $c$:
\begin{align*}
    \hat{Q}^+ (b) &= (\abs{\hat{N}^+_1(b)}, \abs{\hat{N}^+_2(b)}, \abs{\hat{N}^+_3(b)}, \abs{\hat{N}^+_4(b)}) &
    \hat{Q}^+ (c) &= (\abs{\hat{N}^+_1(c)}, \abs{\hat{N}^+_2(c)}, \abs{\hat{N}^+_3(c)}, \abs{\hat{N}^+_4(c)}) \\
        &= (\abs{\{b\}}, \abs{\{a, b\}}, \abs{\{a, b\}}, \abs{\{a, b, c, d\}}) &
        &= (\abs{\varnothing}, \abs{\{c, d\}}, \abs{\{b, c, d\}}, \abs{\{a, b, c, d\}}) \\
        &= (1, 2, 2, 4) &
        &= (0, 2, 3, 4) \\
        \\
    \hat{Q}^- (b) &= (\abs{\hat{N}^-_1(b)}, \abs{\hat{N}^-_2(b)}, \abs{\hat{N}^-_3(b)}, \abs{\hat{N}^-_4(b)}) &
    \hat{Q}^- (c) &= (\abs{\hat{N}^-_1(c)}, \abs{\hat{N}^-_2(c)}, \abs{\hat{N}^-_3(c)}, \abs{\hat{N}^-_4(c)}) \\
        &= (\abs{\{b\}}, \abs{\{a, b\}}, \abs{\{a, b, c, d\}}, \abs{\{a, b, c, d\}}) &
        &= (\abs{\varnothing}, \abs{\{c, d\}}, \abs{\{c, d\}}, \abs{\{a, b, c, d\}}) \\
        &= (1, 2, 4, 4) &&= (0, 2, 2, 4)
\end{align*}
\end{example}

# Implementation
\label{sec:implementation}

We have implemented the firstlast and lastfirst procedures, together with minmax and maxmin, in the R package landmark [@Brunson2021a]. Each procedure is implemented in C++ using Rcpp [@Eddelbuettel2011] and separately in R. For rank-distance--based procedures, the user can choose between $\check{q}$ and $\hat{q}$.
The landmark-generating procedures return the indices of the selected landmarks, optionally together with the sets of indices of the points in the cover set (ball or neighborhood) centered at each landmark.
In addition to the number of landmarks $m$ and either the radius $\eps$ of the balls or the cardinality $k$ of the neighborhoods, the user may also specify additive and multiplicative extension factors for $m$ and for $\eps$ or $k$. These will produce additional landmarks ($m$) and larger cover sets ($\eps$ or $k$) with increased overlaps, in order (for example) to construct more robust nerve complexes.

## Validation

We validated the firstlast and lastfirst procedures against several small example data sets, including that of Example\nbs\ref{ex:rank-distance}.
We also implemented each of the maxmin and lastfirst procedures, using C++ (for Euclidean distances only) and R (which uses the proxy package [@Meyer2021] to calculate distances other than Euclidean), and validated these against each other on several larger data sets, including as part of the benchmark tests reported in the next section.
We invite readers to install the package and experiment with new cases, as well as to request or write any desired additional features.

## Benchmark tests

We benchmarked the C++ and R implementations on three data sets: uniform samples from the unit circle $\Sph^1\subset\R^2$ convoluted with Gaussian noise, samples with duplication from the integer lattice $[0,23]\times[0,11]$ using the probability mass function $p(a,b) \propto 2^{-ab}$, and patients recorded at each critical care unit in MIMIC-III using the RT-similarity measure (Section\nbs\ref{sec:data}).
We conducted benchmarks using the bench package [@Hester2020] on the University of Florida high-performance cluster HiPerGator.

\begin{figure}
\includegraphics[width=.6666667\textwidth]{../figures/benchmark-circle-lattice}
\includegraphics[width=.3333333\textwidth]{../figures/benchmark-mimic}
\caption{
Benchmark results for computing landmarks on two families of artificial data (circle and lattice) and one collection of empirical data (RT-similarity space of critical care units in MIMIC-III). Some points are missing because benchmark tests did not complete within 1 hour.
\label{fig:benchmark}
}
\end{figure}

Benchmark results are reported in Figure\nbs\ref{fig:benchmark}.
R implementations used orders of magnitude more memory and took slightly longer. They appeared to scale slightly better in terms of time and slightly worse in terms of memory.
The additional calculations required for the lastfirst procedure increase runtimes by a median factor of 2.5 in our R implementations. The C++ implementation of lastfirst is based on combinatorial definitions and not optimized for speed, and as a result takes much longer---a median factor of almost 2000 relative to maxmin in C++---and failed to complete in many of our tests.

# Experiments
\label{sec:experiments}

## Empirical data
\label{sec:data}

### MIMIC-III

The open-access critical care database MIMIC-III ("Medical Information Mart for Intensive Care"), derived from the administrative and clinical records for 58,976 admissions of 46,520 patients over 12 years and maintained by the MIT Laboratory for Computational Physiology and collaborating groups, has been widely used for education and research [@Goldberger2000; @Johnson2016].
For our analyses we included data for patients admitted to five care units: coronary care (CCU), cardiac surgery recovery (CSRU), medical intensive care (MICU), surgical intensive care (SICU), and trauma/surgical intensive care (TSICU).[^mimic-units]
For each patient admission, we extracted the set of ICD-9/10 codes from the patient's record and several categorical demographic variables: age group (18–29, decades 30–39 through 70–79, and 80+), recorded gender (M or F), stated ethnicity (41 values),[^ethnicity] stated religion (Catholic, unspecified/unobtainable/missing, Protestant Quaker, Jewish, other, Episcopalian, Greek Orthodox, Christian Scientist, Buddhist, Muslim, Jehovah's Witness, Unitarian-Universalist, 7th Day Adventist, Hindu, Romanian Eastern Orthodox, Baptist, Hebrew, Methodist, or Lutheran), marital status (married, single, widowed, divorced, unknown/missing, separated, or life partner), and type of medical insurance (Medicare, private, Medicaid, povernment, or self pay).
Following @Zhong2020, we transformed these _relational-transaction (RT)_ data into a binary case-by-variable matrix $X \in \B^{n \times p}$ suitable for the cosine similarity measure, which was converted to a distance measure by subtraction from 1.
Because cosine similarity is monotonically related to the angle metric, our topological results will be the same up to this rescaling, so for simplicity we use cosine similarity in our experiments.

[^mimic-units]: <https://mimic.physionet.org/mimictables/transfers/>
[^ethnicity]: White, Black/African American, Unknown/Not Specified, Hispanic or Latino, Other, Unable to Obtain, Asian, Patient Declined to Answer, Asian – Chinese, Hispanic Latino – Puerto Rican, Black/Cape Verdean, White – Russian, Multi Race Ethnicity, Black/Haitian, Hispanic/Latino – Dominican, White – Other European, Asian – Asian Indian, Portuguese, White – Brazilian, Asian – Vietnamese, Black/African, Middle Eastern, Hispanic/Latino – Guatemalan, Hispanic/Latino – Cuban, Asian – Filipino, White – Eastern European, American Indian/Alaska Native, Hispanic/Latino – Salvadoran, Asian – Cambodian, Native Hawaiian or Other Pacific Islander, Asian – Korean, Asian – Other, Hispanic/Latino – Mexican, Hispanic/Latino – Central American (Other), Hispanic/Latino – Colombian, Caribbean Island, South American, Asian – Japanese, Hispanic/Latino – Honduran, Asian – Thai, American Indian/Alaska Native Federally Recognized Tribe

### Mexican Department of Health

The Mexican Department of Health (MXDH) has released official open-access data containing an assortment of patient-level clinical variables related to COVID-19 infection and outcomes. These data have been compiled into a database and made freely available on Kaggle[^kaggle], a collaborative data science platform, where they are maintained and updated regularly. The database includes information regarding over 724,000 patients confirmed to be COVID-positive via diagnostic laboratory testing. Two main types of information are present for each patient: (1) temporal data, and (2) categorical or binary variables. The temporal data consist of key dates associated with the clinical course of infection such as date of symptom onset, date of admission to a healthcare institution, and date of death (if applicable). The categorical or binary fields encode clinical factors likely to be associated with COVID-19 infection, severity, or outcome. These variables include information such as sex, state of patient residence, and intubation status, as well as binary fields encoding the presence or absence of a wide variety of comorbidities such as asthma, hypertension, cardiovascular disease. (For a full description of each field included in the data set, see Kaggle.*) Though these variables are categorical rather than continuous/numeric, there are sufficiently many of them (~50) to distinguish between many phenotypic subtypes of COVID-19 patients. Further, this data set is very complete in that every patient is required to contain a valid value for every field, which minimizes concerns around missing data.

[^kaggle]: <https://www.kaggle.com/lalish99/covid19-mx>

## Homology recovery

<!--
### Trefoil sample

...
We hypothesized that, as measured by the fraction of their parameter ranges, the maxmin and lastfirst covers would more consistently recover the homology of the circle than the regular interval and quantile covers.

### Spherical sample with variable density
-->

We compared the suitability of three landmarking procedures (uniformly random, maxmin, lastfirst) on datasets with varying density and duplication patterns by extending an example of @deSilva2004. Each expriment proceeded as follows: We sampled $n=540$ points from the sphere $\Sph^2\subset\R^3$ and selected $k=12$ landmark points. We then used the landmarks to compute PH and recorded the statistics $R_0,R_1,K_0,K_1$ as defined by @deSilva2004. From these statistics we compute the _relative dominance_ $(R_1 - R_0) / K_0$ and _absolute dominance_ $(R_1 - R_0) / K_1$ of the last interval over which all Betti numbers are correctly calculated.

The points $x=(r,\theta,\phi)$ were sampled using four procedures: uniform sampling, skewed sampling, uniform sampling with skewed boosting, and skewed sampling with skewed boosting. The first procedure was used by @deSilva2004 and here serves as a baseline case. For a sample $S$ (with multiplicities) generated from each of the other three procedures, the expected density $\lim_{\eps\to 0}\lim_{n\to\infty}\frac{1}{n}\abs{\{x=(r,\theta,\phi)\in S\mid \alpha-\eps<\phi<\alpha+\eps\}}$ of points near a given latitude $\alpha\in[0,\pi]$ is proportional to the quartic function $p:[0,1]\to[0,1]$ defined by $p(x)=(\frac{\phi}{\pi})^4$.
Skewed sampling is performed via rejection sampling: Points $x_i=(r_i,\theta_i,\phi_i)$ are sampled uniformly and rejected at random if a uniform random variable $t_i\in[0,1]$ satisfies $(\frac{\phi_i}{\pi})^\alpha<t_i$ until $n$ points have been kept [@Diaconis2013].
Skewed boosting is performed by first obtaining a (uniform or skewed) sample $T$ of size a fraction $\frac{n}{6}$ of the total, then sampling $n$ points (with replacement) from $T$ using the probability mass function satisfying $P(x_i)\propto(\frac{\phi_i}{\pi})^\beta$.
When performed separately, skewed sampling and skewed boosting use $\alpha=\beta=4$; when performed in sequence, they use $\alpha=\beta=2$.

The landmark points were selected in three ways: uniform random selection (without replacement), the maxmin procedure, and the lastfirst procedure.
We computed PH in Python GUDHI, using three implementations: Vietoris–Rips (VR) filtrations on the landmarks, alpha complexes on the landmarks, and witness complexes on the landmarks with the full sample as witnesses [@Maria2021; @Rouvreau2021; @Kachanovich2021].

The skewed data sets are dense (or multiple) at the south pole and sparse (or singular) at the north pole. We expect lastfirst to be more sensitive to this variation and place more landmarks toward the south. As measured by dominance, therefore, we hypothesized that lastfirst would be competitive with maxmin when samples are uniform and inferior to maxmin when samples are skewed.
Put differently, we expected lastfirst to better reject the homology of $\Sph^2$, i.e. to detect the statistical void at the north pole.

Figure\nbs\ref{fig:sphere} compares the relative dominance of the shperical homology groups in each case.
When PH is computed using VR or alpha complexes, maxmin better recovers the homology of the sphere except on uniform samples, while lastfirst and random selection better detect the void.
Random selection is usually better than lastfirst selection at detecting this void when samples are non-uniform, which indicates that lastfirst selection still oversamples from less dense regions.
Lastfirst and maxmin perform similarly when PH is computed using witness complexes.

\begin{figure}
\includegraphics[width=\textwidth]{../figures/homology-sphere-relative}
\caption{
Relative dominance of the spherical homology groups in the persistent homology of four samples from the sphere, using each of three landmark procedures and three persistence computations. Similar plots of absolute dominance (not shown) tell a consistent story, but the distributions are more skewed so the comparisons are less clear.
\label{fig:sphere}
}
\end{figure}

## Covers and nerves

Cardinality reduction techniques can be used to model a large number of cases represented by a large number of variables as a smaller number of clusters with similarity or overlap relations among them.
The deterministic sampling procedures maxmin and lastfirst provide clusters (cover sets) defined by proximity to the landmark cases and relations defined by their overlap.
The clusters obtained by these procedures occupy a middle ground between the regular intervals or quantiles commonly used to cover samples from Euclidean space and the emergent clusters obtained heuristically by penalizing between-cluster similarity and rewarding within-cluster similarity.
The maxmin procedure produces cover sets of (roughly) fixed radius, analogous to overlapping intervals of fixed length, while the lastfirst procedure produces cover sets of fixed size, analogous to the quantiles of an adaptive cover.
This makes them natural solutions to the task of covering an arbitrary finite metric space that may or may not contain important geometric or topological structure [@Singh2007].

As a practical test of this potential, we loosely followed the approach of @Dlotko2019 to construct covers and their nerves for each care unit of MIMIC-III, using maxmin and lastfirst.
We varied the number of landmarks (6, 12, 24, 36, 48, 60, 120) and the multiplicative extension of the cover sets' sizes (0, .1, .2).
We evaluated the procedures in three ways:

- **Clustering quality:** Both procedures yield _fuzzy_ clusters---that is to say, clusters that allow for some overlap. While clustering quality measures might be useful, most, including almost all that have been proposed for fuzzy clusterings, rely on coordinate-wise calculations, specifically data and cluster centroids [@Bouguessa2006; @Wang2007; @Falasconi2010]. To our knowledge, the sole exception to have appeared in a comprehensive comparison of such measures is the _modified partition coefficient_ [@Dave1996], defined as $$\operatorname{MPC}=1-\frac{k}{k-1}(1-\frac{1}{n}\sum_{i=1}^{n}{\sum_{j=1}^{k}{{u_{ij}}^2}})$$ where $U=(u_{ij})$ is the $n\times k$ fuzzy partition matrix: $u_{ij}$ encodes the extent of membership of point $x_i$ in cluster $c_j$, and $\sum_{j=1}^{k}{u_{ij}}=1$ for all $i$. When a point $x_i$ is contained in $m$ cover sets $c_j$, we equally distribute its membership so that $u_{ij}=\frac{1}{m}$ when $x_i\in c_j$ and $u_{ij}=0$ otherwise. Thus, the MPC quantifies the extent of overlap between all pairs of clusters. Like the partition coefficient from which it is adapted, the MPC takes the value $1$ on crisp partitions and is penalized by membership sharing, but it is standardized so that its range does not depend on $k$.
- **Discrimination of risk:** For purposes of clinical phenotyping, patient clusters are more useful that better discriminate between low- and high-risk subgroups. We calculate a cover-based risk estimate from individual outcomes $y_i$ as follows: For each cover set $c_j\subset X$, let $p_j=\frac{1}{\abs{c_j}}\sum_{x_i\in c_j}{y_i}$ be the incidence of the outcome in that set. Then compute the weighted sum $q_i=\sum_{x_i\in c_j}{u_{ij}p_j}$ of these incidence rates for each case. We measure how well the cover discriminates risk as the area under the receiver operating characteristic curve (AUROC).
<!--
- **Homological richness:** Independent of the "true" topological structure of a data set, a simplicial complex model is more useful when it is not too sparse to be connected (or nearly so) and not too dense to visualize or contain low-degree cycles. We compare the $0$- and total $>0$-degree Betti numbers of the simplicial complex models from the perspective of optimizing $-\beta_0+\sum_{i>0}{\beta_i}$.
-->

We hypothesized that lastfirst covers would exhibit less overlap than maxmin covers by virtue of their greater sensitivity to local density, and that they would outperform maxmin covers at risk prediction by reducing the sizes of cover sets in denser regions of the data (taking advantage of more homogeneous patient cohorts).

Figure\nbs\ref{fig:cover-mimic} presents the sizes of the nerves of the covers and the two evaluation statistics as functions of the number of landmarks.
The numbers of 1- and of 2-simplices grew at most roughly quadratically and roughly cubically, respectively. This suggests that the densities of the simplicial complex models were at most roughly constant, regardless of the number of landmarks.
Landmark covers grew fuzzier and generated more accurate predictions until the number of landmarks reached around 60, beyond which point most covers grew crisper while performance increased more slowly (and in one case decreased). This pattern held for covers with any fixed multiplicative extension.
Naturally, these extensions produced fuzzier clusters, but they also reduced the overall accuracy of the predictive model.
In addition to (i.e. independently of) these patterns, models fitted to smaller care units tended to outperform those fitted to larger care units.
Contrary to expectations, unextended maxmin covers were usually crisper than their lastfirst counterparts and yielded more accurate predictions, though extensions reduced the crispness of maxmin covers more dramatically than of lastfirst covers.
The same patterns were observed in the risk discrimination of maxmin versus lastfirst covers, with maxmin covers yielding the most accurate predictions when unextended but lastfirst covers retaining more accuracy after extension.

\begin{figure}
\includegraphics[width=.5\textwidth]{../figures/cover-simplices}
\includegraphics[width=.5\textwidth]{../figures/cover-evaluate}
\caption{
Summary and evaluation statistics versus number of 0-simplices (landmarks) for the covers generated using the maxmin and lastfirst procedures, with three multiplicative extensions in their size.
Left: the sizes of their nerves, as numbers of 1- and 2-simplices, using a square root--transformed vertical scale.
Right: the modified partition coefficient (MPC) and the c-statistic of the risk prediction model based on the cover sets (AUROC).
\label{fig:cover-mimic}
}
\end{figure}

Figure\nbs\ref{fig:cover-mx} presents the same evaluations for covers of the MXDH data.
In contrast to the MIMIC experiments, lastfirst-based nerves of the MXDH data grow sub-polynomially and are significantly sparser than maxmin-based nerves.
Lastfirst covers tend to be crisper, especially as the number of landmarks and the extension factors increase.
This indicates that the nearest neighborhoods form a more parsimonious cover of the data than the centered balls.
The predictive accuracies of the maxmin- and lastfirst-selected cover-based models converge with increasing landmarks^[We should increase the maximum number of landmarks to be sure.], though for smaller numbers different selection procedures perform best for different outcomes.

\begin{figure}
\includegraphics[width=.5\textwidth]{../figures/cover-simplices-mx}
\includegraphics[width=.5\textwidth]{../figures/cover-evaluate-mx}
\caption{
Summary and evaluation statistics versus number of 0-simplices (landmarks) for the covers generated using the maxmin and lastfirst procedures, with three multiplicative extensions in their size.
Left: the sizes of their nerves, as numbers of 1- and 2-simplices, using a square root--transformed vertical scale.
Right: the modified partition coefficient (MPC) and the c-statistics of the risk prediction models based on the cover sets (AUROC).
\label{fig:cover-mx}
}
\end{figure}

## Interpolative nearest neighbors prediction

Landmark points may also be used to trade accuracy for memory in neighborhood-based prediction modeling. Consider the following approach: A modeling process involves predictor data $X \in \R^{n \times p}$ and response data $y \in \R^{n \times 1}$, partitioned into training and testing sets $X_0,X_1$ and $y_0,y_1$ according to a partition $I_0 \sqcup I_1 = \{1,\ldots,n\}$ of the index set. Given $x \in X_1$, a nearest neighbors model computes the prediction $p(x) = \frac{1}{k}\sum_{q(x,x_i) \leq k}{y_i}$ by averaging the responses for the $k^\text{th}$ nearest neighbors of $x$ in $X_0$. By selecting a landmark set $L \subset X_0$, a researcher can reduce the computational cost of the model as follows: For each $\ell \in L$, calculate $p(\ell)$ as above. Then, for each $x \in X_1$, calculate $p_L(x) = \sum_{\ell \in L}{w(d(x,\ell)) p(\ell)} / \sum_{\ell \in L}{w(d(x,\ell))}$, where $w : \R_{\geq 0} \to \R_{\geq 0}$ is a weighting function (for example, $w(d)=d^{-1}$). The nearest neighbor predictions for $L$ thus serve as proxies for the responses associated with $X_0$.

We took this approach to the prediction of in-hospital mortality for patients with records in each critical care unit of MIMIC-III.
We then implemented the following procedure:

<!--Use time series nested cross-validation for the CTD data?-->
1. Determine a nested $6 \times 6$–fold split for train–fit–test cross-validation. That is, partition $[n] = \bigsqcup_{i=1}^{6}{I_i}$ into roughly equal parts, and partition each $[n] \wo I_i = \bigsqcup_{j=1}^{6}{J_{ij}}$ into roughly equal parts.
2. Iterate the following steps over each $i,j$:
    a) Generate a sequence $L$ of landmarks from the points $X_{([n] \wo I_i) \wo J_{ij}}$.
    b) Identify the $180$ nearest neighbors $N^+_{180}(\ell)$ of each landmark $\ell$. This was a fixed parameter, chosen for being slightly larger than the optimal neighborhood size in a previous study of individualized models [@Lee2015].
    c) Find the value of $k \in [180]$ and the weighting function $w$ (among those available) for which the predictions $p_L : X_{J_{ij}} \to [0,1]$ maximize the AUROC.
    d) Use the AUROC to evaluate the performance of the predictions $p_L : X_{I_i} \to [0,1]$ using these $k$ and $w$.

We replicated the experiment for each combination of procedure (random, maxmin, lastfirst) and number of landmarks ($\abs{L}=36,60,180,360$).
We hypothesized that, as measured by overall accuracy of the resulting predictive model, the maxmin and lastfirst procedures would outperform random selection, and that lastfirst would outperform maxmin, for similar reasons to those in the previous section.

Boxplots of the AUROCs for each cross-validation step are presented in Figure\nbs\ref{fig:knn-mimic}.
While both landmark procedures yielded stronger results than random selection, lastfirst performed on average slightly worse than maxmin on each data set.
Importantly, both landmark procedures also yielded more accurate predictions than a basic unweighted nearest-neighbors model, lending support to the modeling approach itself.
Interestingly, only on the largest data set (the MICU) did increasing the number of landmarks from 36 to 360 appreciably improve predictive accuracy (using all three selection procedures).

\begin{figure}
\includegraphics[width=\textwidth]{../figures/knn-auc}
\caption{
AUROCs of the interpolative predictive models of mortality in five MIMIC-III care units based on covers constructed using random, maxmin, and lastfirst procedures to generate landmarks.
Each boxplot summarizes AUROCs from $6 \times 6 = 36$ models, one for each combination of outer and inner fold.
AUROCs of simple nearest-neighbor predictive models are included for comparison.
\label{fig:knn-mimic}
}
\end{figure}

Over the course of the COVID-19 pandemic, hospitals and other facilities experienced periods of overburden and resource depletion, and best practices were continually learned and disseminated.
As a result, outcomes in the MXDH data reflect institutional- as well as population-level factors.
We took advantage of the rapid learning process in particular by adapting the nested CV approach above to a nested-longitudinal CV approach:
We partitioned the data by week, beginning with Week\nbs11 (March\nbs11--17) and ending with Week\nbs19 (May\nbs6--9, the last dates for which data were available).
For each week $i$, $11 < i \leq 19$, we trained prediction models on the data from Week $i-1$.
We then randomly partitioned Week\nbs$i$ into six roughly equal parts and optimized and evaluated the models as above.
(For this analysis, we only considered Gaussian weighting.)

Line plots of model performance are presented in Figure\nbs\ref{fig:knn-mx}, with one curve (across numbers of landmarks) per selection procedure, outcome, and week.
Again both landmark selection procedures yielded stronger results than random selection. Interestingly, this was more pronounced in later weeks, as the pandemic progressed, even as overall predictive accuracy declined.^[Do we know what major events in Brazil might have influenced these outcomes? We should at least also plot the number of cases, stratified by outcomes, per day over this period.]
Overall, performance improved slightly as the number of landmarks increased from 50 to 150 but either plateaued or declined from 150 to 250.

\begin{figure}
\includegraphics[width=\textwidth]{../figures/knn-auc-mx-gaussian}
\caption{
AUROCs of the sliding-window interpolative predictive models of intubation and mortality in the MXDH data based on covers constructed using random, maxmin, and lastfirst procedures to generate landmarks.
\label{fig:knn-mx}
}
\end{figure}

# Discussion

The definitions of our lastfirst and firstlast procedures are analogous to those of maxmin and minmax, substituting ranks in the role of distances.
In this way, lastfirst is an alternative to maxmin that is adaptive to the local density of the data, similar to the use of fixed quantiles in place of fixed-length intervals.
The maxmin and lastfirst procedures implicitly construct a minimal cover whose sets are centered at the selected landmarks, and the fixed-radius balls of maxmin correspond to the fixed-cardinality neighborhoods of lastfirst.
The rank-based procedures are more combinatorially complex and computationally expensive, primarily because rank-distances are asymmetric, which doubles (in the best case) or squares (in the worst case) the number of distances that must be calculated.
Nevertheless, the procedure can be performed in a reasonable time for many real-world uses.

We ran an experiment to compare maxmin and lastfirst landmarks at expediting the computation of persistent homology, extending an experiment of @deSilva2004. In addition to a uniform sample from a sphere, we drew skewed and bootstrapped samples in order to simulate data sets with variable density and multiplicity, in this case exhibiting a statistical void at one pole of the sphere, opposite a concentration at the other pole.
Whereas the maxmin procedure would sample more uniformly across the sphere despite this skew, the lastfirst procedure would concentrate landmarks toward the south.
Classical persistent homology is notoriously sensitive to outliers, and maxmin better recovered spherical homology than lastfirst, as expected. This is likely due in part to the filtration itself being based on distances rather than rank-distances between landmarks (and other points as witnesses).
Yet, compared to random sampling, lastfirst still oversampled from less dense regions. This is desirable for settings, such as healthcare, in which regions of missingness often indicate limitations of the data collection rather than rarity of case types in the population.

We also ran several experiments that used landmarks to obtain well-separated clusters of patients with common risk profiles and to more efficiently generate nearest neighbor predictions.
Because we designed lastfirst to produce cover sets of equal size despite variation in the density or multiplicity of the data, we expected it to outperform maxmin with respect to the crispness of clusters and to the accuracy of predictions.
In particular, we expected that the optimal neighborhood size for outcome prediction would be roughly equal across our data; as a result, by assigning each landmark case an equally-sized cohort of similar cases, we expected predictions based on these cohorts to outperform those based on cohorts using a fixed similarity threshold.

Contrary to expectations, maxmin produced crisper clusterings on average, and in the case of MIMIC-III more accurate predictions.
However, when the sets of these covers had their radii or cardinalities extended by a fixed proportion, those of lastfirst better preserved these qualities.
Additionally, in the case of MXDH, neither landmark selection procedure produced consistently more accurate predictions.

A possible explanation for the stronger performance of maxmin on MIMIC is that the data did not exhibit very strongly the patterns for which the lastfirst procedure is designed to account, namely variation in density and multiplicity.
As a result, the RT-similarity measure is in fact meaningful across the population: Whatever the baseline presentation of a patient, rather than a cohort of similar patients of some fixed size, their prognosis would be better guided by a cohort cut off at a fixed minimum similarity (or maximum distance).
This suggests that the use of personalized cohorts to improve predictive modeling, as employed by @Lee2015, may be strengthened by optimizing a fixed similarity threshold rather than a fixed cohort size.
It is worth noting that @Park2006 reached a similar conclusion.
In contrast, this stronger performance was not evident on MXDH, which contained fewer variables and as a result exhibited many more occurrences of duplication.

Another way to think about these results is in terms of a balance between relevance and pwoer, with fixed-radius balls (respectively, fixed-cardinality neighborhoods) providing training cohorts of roughly equal relevance (statistical power) to all test cases.
With sufficiently rich data, relevance can be more precisely measured and becomes more important to cohort definition, as with MIMIC.
When variables are fewer, as with MXDH, relevance is more difficult to measure, so that larger samples can improve performance even at the expense of such a measure.


# References

<!--
# To generate the LaTeX file, execute the following:
pandoc math-draft.md \
  -s \
  --number-sections \
  --bibliography=../lastfirst.bib \
  -o math-draft.tex
# To generate the PDF directly, execute the following:
pandoc math-draft.md \
  -t latex \
  --number-sections \
  --bibliography=../lastfirst.bib \
  --citeproc \
  -o math-draft.pdf
-->
