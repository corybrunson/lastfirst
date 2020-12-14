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

Topological data analysis (TDA) is a maturing field in data science, at the interface of statistics, computer science, and mathematics. Topology is often described as the mathematical study of shape but can also be defined as the study of continuity, and TDA consists in the use of computational theories of continuity to model the shape of data. While TDA is most commonly associated with persistent homology (PH) and mapper, it encompasses and generalizes many conventional and even classical techniques, including cluster analysis, association network analysis, and $k$-nearest neighbors (kNN) prediction.

Two basic maneuvers in TDA are locality preservation and cardinality reduction. \emph{Locality preservation} is the property of functions (such as annotations or hashes) that nearby cases or points in the domain have nearby images in the codomain. This is the defining property of many non-linear dimensionality reduction techniques, including t-SNE and UMAP, but also an asymptotic property of kNN prediction and the locally-sensitive hashing used in its implementations. We use the term \emph{cardinality reduction}\footnote{The term "cardinality reduction" takes different meanings in the data analysis literature, including the combining of different values of a categorical variable [@MicciBarreca2001; @Refaat2010] or of time series [@Hu2011] (also "numerosity reduction" [@Lin2007]). Our meaning is that of [@Byczkowska-Lipinska2017]: an $n\times p$ data table of $n$ cases (rows) and $p$ variables (columns) can be dimension-reduced to an $n\times q$ table, where $q<p$, or cardinality-reduced to an $m\times p$ table, where $m<n$. The most common cardinality reduction method is data reduction.}, in contrast to dimension reduction, to describe techniques that produce more tractable or comprehensible representations of complex data by reducing the number of cases or points rather than the number of variables or spatial dimensions used to represent them. Cardinality reduction describes the classical fields of cluster analysis and association (co-occurrence or correlation) network analysis, for example.

More elaborate TDA techniques combine locality preservation and cardinality reduction, often through the use of covering methods, as with the approximation of PH through witness complexes or in the mapper construction. Covering methods, in turn, can be enhanced by strategic sampling from large data sets, as can other spatial techniques like kNN. The maxmin procedure is often used for this purpose, as it is computationally efficient and generates more evenly distributed samples than random selection. However, maxmin comes with its own limitations, which may be more acute in the analysis of data that vary greatly in density or have multiplicities. This is a frequent concern when sparse, heterogeneous, and incomplete data are modeled as finite pseudometric spaces using similarity measures rather than vector space embeddings.

In this paper, we develop a complementary landmark sampling procedure to maxmin based on the rankings of points by a distance or similarity measure rather than on the raw values of such a measure. In the remainder of this section, we motivate the procedure, which we call lastfirst, as a counterpart to maxmin obtained by adapting this procedure from the use of fixed-radius balls to the use of fixed-cardinality neighborhoods. We then formally describe the procedure and prove some of its basic properties in Section\nbs\ref{sec:procedures}. In Section\nbs\ref{sec:evaluations} we report the results of benchmark tests and robustness analyses of lastfirst on simulated and empirical data sets. We describe some basic applications to real-world data in Section\nbs\ref{sec:applications}.

## Conventions

$(X, d_X)$ will refer to a finite pseudometric space with point set $X$ and pseudometric $d_X:X\times X\to\R_{\geq 0}$; $(X,d_X)$ may be shortened to $X$, and $d_X$ to $d$, when clear from context. The cardinality of $Y\subseteq X$ (including multiplicities) is denoted $\abs{Y}$, and the set of distinguishable points of $Y$ is denoted $\supp{Y}$. Throughout, we take $n=\abs{X}$. When $Y,Z\subseteq X$, let $Y\wo Z$ denote the set of points in $Y$ (with multiplicities) that are distinguishable from all points in $Z$. This means that, when nonempty, $\min_{y\in Y\wo Z,z\in Z}{d(y, z)}>0$. We denote the diameter of $Y$ as $D(Y)=\max_{x,y\in Y}{d(x,y)}$, and we write:
\begin{align*}
d(x,Y) &= \min_{y\in Y}{d(x,y)} & D(x,Y) &= \max_{y\in Y}{d(x,y)} \\
d(Y,Z) &= \min_{y\in Y,z\in Z}{d(y,z)} & D(Y,Z) &= \max_{y\in Y,z\in Z}{(y,z)} \\
\end{align*}

$f:X \to Y$ will denote a morphism of pseudometric spaces, meaning a function, or morphism of sets, that preserves locality in the sense that $d_X(x,y)\geq d_Y(f(x),f(y))$.
If $x\neq y$ but $d(x,y)=0$ then $x$ and $y$ are said to be co-located.
If $d(x,y)=d(x,z)\implies y=z$, then $X$ is considered to be in general position, even if $d(x,y)=d(z,w)\nimplies \{x,y\}=\{z,w\}$. This is condition implies that $X$ is Hausdorff: $d(x,y)=0\implies x=y$.

We use the ball notation $B_{\eps}(x)$ for the set of points less than distance $\eps$ from a central point $x$; that is, $B_{\eps}(x) = \{ y \mid d(x,y) < \eps \}$.
We use an overline to include points exactly distance $\eps$ from $x$: $\cl{B_{\eps}}(x) = \{ y \mid d(x,y) \leq \eps \}$. (While these connote openness and closedness, every such ball is both open and closed.)
When $\abs{\cl{B_\eps}(x)}\geq k$ and $\eps'<\eps\implies\abs{\cl{B_\eps'}(x)}<k$, then we say that $N^+_k(x)=\cl{B_\eps}(x)$ is the $k$-nearest neighborhood of $x$---though we will complicate this definition in the next section. When $X$ is in general position, $\abs{N_k(x)}=k$.

For convenience, we assume $0\in\N$.

## Background

The lastfirst procedure addresses a peculiar disadvantage of the maxmin procedure that arises when, due to the use of binary or categorical variables or to limits on measurement resolution, a data set includes many duplicate or otherwise indistinguishable cases. These render the finite metric space representation of the data non-Hausdorff: certain subsets of points, called \emph{points with multiplicity}, represent distinct cases but have zero distance. While these issues may be negligible when such points are rare, they raise computational and interpretative concerns when common. Because our procedure is motivated by the same practical needs that motivate the use of the maxmin procedure, we begin with a discussion of those needs.

The earliest appearance of the maxmin\footnote{Note that this procedure is distinct from many other uses of "maxmin" and related terms.} procedure of which we are aware is [@deSilva2004]. The authors propose witness complexes, later generalized to alpha complexes [@], for the rapid approximation of persistent homology: Given a point cloud, a set of landmark points and their overlapping neighborhoods define a nerve, which stands in for the Vietoris--Rips complex at each scale. They use the maxmin procedure, which we define in Section\nbs\ref{sec:maxmin}, as an alternative to selecting the landmark points uniformly at random. The procedure ensures that the landmarks are locally separated and roughly evenly dispersed throughout the point cloud.
While the procedure improved little upon uniform random selection in most use cases, on some tasks it far outperformed.

Subsequent uses of maxmin include the selection of a sample of points from a computationally intractable point cloud for the purpose of downstream topological analysis, as when performing the mapper construction [@Singh2007]; and the optimization of a fixed-radius ball cover of a point cloud, in the sense of minimizing both the number of balls and their shared radius [@Dlotko2019]. In addition to approximating persistent homology [@deSilva2004; @Dlotko2019], maxmin has been used to reduce the sizes of simplicial complex models of point cloud data for the sake of visualization and exploration [@Singh2007; @Dlotko2019].

## Motivation

The ball covers mentioned above have been proposed as an alternative to mapper [@Dlotko2019], where they exchange complexity for computational cost, and the mapper construction itself relies on a crucial covering step that has received limited theoretical attention.
Conventionally, mappers rely on covers consisting of overlapping intervals (when the lens is one-dimensional) or of their cartesian products (higher-dimensional).
For this purpose, we propose that ball covers, heuristically optimized or near-optimized using the maxmin procedure, have a potential advantage over conventional covers, alongside a potential disadvantage.

Conventionally, mappers use low-dimensional lens spaces $\mathbb{R}^m$ and one of two types of cover, based either on overlapping intervals of fixed length or on overlapping quantiles of fixed cardinality---or roughly fixed, in case of multiplicity.
We think of this length or cardinality as the resolution of the cover.
When $m>1$, covers for $Y\subset\mathbb{R}^m$ can be obtained as the cartesian products of those of the coordinate projections of $Y$---so, if $\pi_1,\pi_2:\mathbb{R}^2\to\mathbb{R}$ are the coordinate projections, then obtain interval covers, consisting of $I_\alpha$ and $J_\beta$, for $\pi_1(Y)$ and $\pi_2(Y)$, so that the rectangle $I_\alpha \times J_\beta$ contains $y\in Y$ when $\pi_1(y)\in I_\alpha$ and $\pi_2(y)\in J_\beta$.
While these cover types are manageable in very low dimensions, the number of sets scales exponentially with $m$, holding the resolution fixed.
Moreover, eventually most of the resulting cover sets will contain no points of $X$, and additional calculations will be needed to restrict to the non-empty sets.

In contrast, a cover obtained by centering balls at a subset of landmark points in $Y$ will have greater up-front computational cost but will be guaranteed to leave no empty sets, and the number of sets required to capture the topology of $Y$ will increase only with the topological complexity of $Y$, not with the dimension $m$. (We test this hypothesis using a point cloud sampled from a trefoil knot embedded in $\mathbb{R}^3$.)

Nevertheless, the maxmin cover relies on a meaningful distance metric: the dissimilarity of cases $x$ and $y$ is captured by their distance $d(x,y)$, regardless of where $x$ and $y$ are located in $X$, and the neighborhoods $B_r(x)$ and $B_r(y)$ about landmarks $x$ and $y$ play an equal role in the cover.
This means that cover sets centered at landmarks in sparse regions of $X$ will be more numerous and of lower cardinality than those centered in dense regions.
The assumption is violated in many real-world settings, including much of biomedicine. For example, in psychometric terms, this would mean that inter-case distance is an \emph{interval}, not only an \emph{ordinal}, variable, so that the distances between cases in a point cloud representation has a definite meaning independent of which cases are considered. This assumption is often made for convenience, but it generally does not follow from theory.

One of our motivations is to produce a counterpart to the ball cover, the (nearest) neighborhood cover, each set of which may have a different radius but (nearly) the same cardinality. Especially in analyses of medical and healthcare data, both underlying variables and multivariate similarity measures can often only be understood as ordinal. Other representations of high-dimensional data sets are commonly defined by similarity (or dissimilarity) measures such as cosine similarity rather than by coordinates and associated metrics. Furthermore, because measurements are coarse and often missing, such data often contain indistinguishable entries: cases all of whose measurements are equal and that are therefore represented as multiple instances of the same point. All of these attributes violate the assumptions of the ball cover approach and suggest the need for an ordinal counterpart.

\pagebreak

# Procedures

In this section we review the maxmin procedure as a prelude and introduce the lastfirst procedure as a complement to it.

## Maxmin procedure

As defined for the purpose of landmark selection, maxmin takes as input a proper subset $L\subset X$ and returns as output a point $x\in X\wo L$. A related procedure, which we will call minmax, takes as input only $X$ and returns its \emph{Chebyshev center}, the point $x\in X$ for which $D(x,X\wo\{x\})$ is minimized. Whereas maxmin naturally increments a landmark set, minmax naturally initializes one: The Chebyshev center is the landmark that produces the single-set cover using the smallest radius. To help build intuition around these concepts, we take a moment to see how they fit into a "two-by-two table" with two other procedures: Each procedure calculates either minmax or maxmin, and each does so either with reference to a proper subset of $X$ or with respect to $X$ itself. It will be this table that we adapt from the setting of fixed-radius balls to the setting of (roughly) fixed-cardinality neighborhoods.

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
Note that both $\minmax(\,\cdot\,)$ and $\maxmin(\,\cdot\,)$ are nonempty and that, when $X$ is in general position, they both have cardinality $1$.

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
First, choose a number $n\leq\uniq{X}$ of landmark points to generate or a radius $\eps\geq 0$ for which to require that the balls $\cl{B_{\eps}}(\ell)$ minimally cover $X$.
Choose a first landmark point $\ell_0\in X$. This choice may be arbitrary; we specifically consider three selection rules: the first point index in the object representing $X$, selection at random, and (random selection from) $\minmax(X)$.
Inductively over $i\in\N$, if ever $i\geq n$ or $d(L,X\wo L)\leq\eps$, then stop.
Otherwise, when $L=\{\ell_0,\ldots,\ell_{i-1}\}$, choose $\ell_i\in\maxmin(L)$, again according to a preferred rule.
If a fixed number $n$ of landmarks was prescribed, then set $\eps=d(L,X\wo L)$; if $\eps$ was prescribed, then set $n=\abs{L}$.

We will write the elements of landmark sets using set notation $L=\{\ell_0,\ldots,\ell_{n-1}\}$ but always in the order in which they were generated. **Justify using this notation rather than sequence notation.**
Note that, if $n=\uniq{X}$ or $\eps=0$, then $L=\supp{X}$.
When the procedure stops, $X=\bigcup_{i=0}^{n-1}{\cl{B_{\eps}}(\ell_i)}$, and this cover is minimal in two senses: The removal of any $\cl{B_\eps}(\ell_i)$ or any decrease in $\eps$ will obtain a collection of sets that fail to cover $X$. More generally, a non-minimal cover can be obtained by specifying both $n$ and $\eps$ in a compatible way. In Section\nbs\ref{sec:implementation}, we describe two adaptive parameters implemented in our software package that make these choices easier.

## Lastfirst procedure

The lastfirst procedure is defined analogously to the maxmin procedure, substituting "rank-distance" for the pseudometric $d_X$.

### Rank-distances

\textit{Rank-distance} is an adaptive notion of distance with respect to nearest neighborhoods.\footnote{Note that, as used here, rank-distance is distinct from the \emph{rank-distance} between permutations, which is used to define rank correlation coefficients, and from the \emph{ordinal distance} proposed by Pattanaik and Xu (2008), a loosening of the concept of pseudometric that dispenses with the triangle inequality.} It relies on the underlying pseudometric $d_X$ but takes values in $\N$ given by the number of nearest neighbors one point is away from another.

\begin{definition} (Rank-distance)
    For $x,y\in X$, define the \textit{rank-distance} $q_{X,d} : X \times X \longrightarrow \N$ as follows:
    \begin{equation*}
        q_{X,d}(x,y)=\abs{\{z\in X\mid d_X(x,z)<d_X(x,y)\}}+1%>
    \end{equation*}
\end{definition}

As with $d$, we allow ourselves to write $q=q_{X,d}$ when clear from context. Note that $q(x,x)=1$ and $q(x,y)\leq N$, and that $q$ is not, in general, symmetric. (It is therefore not a pseudometric.)

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

We term the unary rankings $q(x,\,\cdot\,)$ and $q(\,\cdot\,,x)$ the \emph{out- (from $x$)} and \emph{in- (to $x$) rankings} of $X$, respectively. These can then be used to define \textit{out-} and \textit{in-neighborhoods} of $x$.\footnote{The terminology and notation are adapted from the theory of directed graphs. These definitions are the same as those for a complete directed graph on $X$ with directed arcs $x\to y$ weighted by rank-distance $q(x,y)$.}

\begin{definition} ($k$-neighborhoods)
    For $x \in X$, define the \emph{$k$-out-neighborhoods} $N^+_k$ and \emph{$k$-in-neighborhoods} $N^-_k$ of $x$ as the sets
    \begin{align*}
        & N^+_k(x)=\{y\in X\mid q(x,y)\leq k\} \\
        & N^-_k(x)=\{y\in X\mid q(y,x)\leq k\}
    \end{align*}
\end{definition}

Note that $\varnothing \subseteq N^\pm_1(x) \subseteq \cdots \subseteq  N^\pm_N(x) = X$.
The $k$-out-neighborhoods of $x$ are the sets of points in $X$ that have rank-distance at most $k$ from $x$. This is equivalent to the $k$-nearest neighbors of $x$. The $k$-in-neighborhoods of $x$ are the sets of points in $X$ from which $x$ has rank-distance at most $k$.
These definitions can be adapted as follows to be relative to a subset $Y \sub X$, using the notation $q_{X,d}$ to emphasize that the points in $X\wo (Y\cup\{x\})$ are still involved in the calculation:

\begin{align*}
    & N^+_k(x,Y)=\{y\in Y\mid q_{X,d}(x,y)\leq k\} \\
    & N^-_k(x,Y)=\{y\in Y\mid q_{X,d}(y,x)\leq k\}
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

Consider the case of choosing the next landmark point given a subset $L=\{\ell_0,\ldots,\ell_{i-1}\}\subset X$ of collected landmark points. For a ball cover, we choose $\ell_i\in X\wo L$ so that the \emph{minimum radius $\eps$} required for some $B_\eps(\ell_j)$ to contain $\ell_i$ is \emph{maximized}. Analogously, for a neighborhood cover, we want that the \emph{minimum cardinality $k$} required for some $N^+_k(\ell_j)$ to contain $\ell_i$ is \emph{maximized}. Switching perspective from out- to in- and reversing the roles of $L$ and $\ell_i$, we want $N^-_k(\ell_i,L)=1$ for the latest (largest) $k$ possible, say $k^-_1$. When indistinguishable points abound, this may still not uniquely determine $\ell_i$, so we may extend the principle: Among those $\ell\in X\wo L$ for which $N^-_{k^-_1}(\ell_i,L)=1$, choose $\ell_i$ for which $N^-_{k}(\ell_i,L)=2$ for the latest $k$ possible, say $k^-_2\geq k^-_1$. Continue this process until only one candidate $\ell$ remains, or until $N^-_{k}(\ell,L)=\abs{L}$, in which case we consider all remaining candidates equivalent.

Now consider the case of choosing an initial landmark point with respect to a subset $L\subset X$. For a ball cover, we above considered the Chebyshev center $\ell_0\in X$ so that the \emph{maximum radius $\eps$} required for $\cl{B_\eps}(\ell_0)$ to contain any $\ell\in L$ is \emph{minimized}. Analogously, for a neighborhood cover, we would want the \emph{maximum cardinality $k$} required for $N^+_k(\ell_0)$ to contain any $\ell\in L$ to be \emph{minimized}. That is, we would want $N^+_k(\ell_0,L)=\abs{L}$ for the earliest (smallest) $k$ possible, say $k^+_{\abs{L}}$. To decide among several points $\ell$ with this property, we would then choose one for which $N^+_k(\ell,X)=\abs{L}-1$ for the earliest $k$ possible, say $k^+_{\abs{L}-1}\leq k^+_{\abs{L}}$. This process would likewise continue until only one candidate $\ell$ remains, or until $N^+_{k}(\ell,L)=0$ and all remaining candidates are considered equivalent.

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
        Under the revlex order, $\max_{x\in X\wo\{a\}}{Q^-(x,\{a\})}=(0,0,1,1)$, and without loss of generality (because $c$ and $d$ are indistinguishable) we select $\ell_1=c$, so that now $L=\{a,c\}$.
        \item Only one point remains in $X\wo L=\{b\}$, so the exhaustive landmark set is $\{a,c,b\}$.
    \end{enumerate}
\end{example}

\subsubsection{Algorithms}

For illustration, we provide complete algorithms to identify the firstlast set of $X$ and a simplified algorithm to obtain a lastfirst landmark set using specified parameters. The former requires first calculating $Q^+(x,L)$ from $x$ and $L$, which is done in Algorithm\nbs\ref{alg:outnbhdseq}.

\begin{algorithm}
\caption{Compute the out-rank sequence of a point with respect to a proper subset.}
\label{alg:outnbhdseq}
\begin{algorithmic}
\REQUIRE finite metric space $(X,d)$
\REQUIRE $L=\{\ell_0,\ldots,\ell_{k-1}\}\subset X$
\REQUIRE $x\in X$
\STATE $D \leftarrow (d(x,\ell_0),\ldots,d(x,\ell_{k-1}))$
\STATE $D \leftarrow \verb|sort|(D)=(d_1\leq\cdots\leq d_k)$
\STATE $d_0 \leftarrow -\infty$; $Q \leftarrow (0,\ldots,0)\in\N^k$; $j \leftarrow 0$
\FOR{$i=1$ to $k$}
    \IF{$d_i > d_{i-1}$}
        \WHILE{$j<i$}
            \STATE $j \leftarrow j+1$
            \STATE $Q_j \leftarrow Q_{j-1}$
        \ENDWHILE
    \ENDIF
    \STATE $Q_j \leftarrow Q_j+1$
\ENDFOR
\RETURN out-rank sequence $Q$
\end{algorithmic}
\end{algorithm}

Making use of these out-rank sequences, the recovery of a firstlast set proceeds as in Algorithm\nbs\ref{alg:firstlast}.

\begin{algorithm}
\caption{Identify a firstlast set with respect to a proper subset.}
\label{alg:firstlast}
\begin{algorithmic}
\REQUIRE finite metric space $(X,d)$
\REQUIRE subset $L=\{\ell_0,\ldots,\ell_{k-1}\}\subset X$
\STATE $F \leftarrow \varnothing$
\STATE $Q \leftarrow (k,\ldots,k)\in\N^{k}$
\FOR{$x\in X$}
    \STATE $Q' \leftarrow Q^+(x,L)\in\N^{k}$
    \STATE $j=k$
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
Algorithm\nbs\ref{alg:firstlast} returns the set of firstlast points.
\end{proposition}

The construction of a lastfirst landmark set proceeds as follows, and as outlined in Algorithm\nbs\ref{alg:lastfirst-landmarks}.
Let $L$ denote the set of landmark points. At the start of the algorithm, $L = \varnothing$.
An initial point $\ell_0$ is selected and added to the landmark set so that $L = \{\ell_0\}$. This point seeds the procedure.
At each iteration, $\lf(L)$ is computed and the next landmark point $\ell_i$ is selected from this set. The procedure terminates when $X=\bigcup_{j\leq i}N^+_1(\ell_j)$, meaning every $x\in X$ is co-located with a point in $L$.

\begin{algorithm}
\caption{Select a lastfirst landmark set.}
\label{alg:lastfirst-landmarks}
\begin{algorithmic}
\REQUIRE finite metric space $(X,d)$
\REQUIRE at least one parameter $k>0$ or $n\in\N$
\REQUIRE seed point $\ell_0 \in X$
\REQUIRE selection procedure \verb|pick|
\IF{only $k$ is given}
    \STATE $n \leftarrow 0$
\ENDIF
\IF{only $n$ is given}
    \STATE $k \leftarrow \abs{X}$
\ENDIF
\STATE $L \leftarrow \varnothing$
\STATE $i \leftarrow 0$
\REPEAT
    \STATE $L \leftarrow L\cup\{\ell_i\}$
    \STATE $i \leftarrow i+1$
    \STATE $F \leftarrow \lf(L)$
    \STATE $\ell_i \leftarrow \verb|pick|(F)$
    \STATE $q_{\operatorname{max}} \leftarrow q(\ell_i,L)$
\UNTIL $q_{\operatorname{max}} < k$ and $\abs{L} \geq n$
\RETURN lastfirst landmark set $L$
\end{algorithmic}
\end{algorithm}

The implementation of this lastfirst procedure is more involved: Because an in-rank sequence $Q^-(x)$ cannot be calculated from the distances $\{d(x,y)\mid y\in X\}$ alone, much more computation is required to generate $Q^-(x,L)$ at each step.

\begin{proposition}
Algorithm\nbs\ref{alg:lastfirst-landmarks} returns a sequence of lastfirst landmark points.
\end{proposition}

### Tie handling

We might have defined two \textit{rank-distances} $\check{q}, \hat{q} : X \times X \longrightarrow \N$ ("$q$-check" and "$q$-hat") as follows:
\begin{align*}
& \check{q}(x,y)=\abs{\{z\in X\mid d(x,z)<d(x,y)\}}+1 \\%>
& \hat{q}(x,y)=\abs{\{z\in X\mid d(x,z)\leq d(x,y)\}}
\end{align*}
In this notation, $\check{q}=q$, while $\hat{q}(x,y)$ is the cardinality of the smallest ball centered at $x$ that contains $y$. Note that then $\check{N}^\pm_1(x) \subseteq  \{x\} \subseteq \hat{N}^\pm_1(x)$, and that $\hat{q}(x,x)>1$ when $x$ has multiplicity.
These two rank-distances correspond to two tie-handling schemes for calculating rankings of lists with duplicates: For example, if $a<b=c<d$ are the distances from $x$ to $y_1,y_2,y_3,y_4$, respectively, then $(\check{q}(x,y_1),\check{q}(x,y_2),\check{q}(x,y_3),\check{q}(x,y_4))=(1,2,2,4)$ and $(\hat{q}(x,y_1),\hat{q}(x,y_2),\hat{q}(x,y_3),\hat{q}(x,y_4))=(1,3,3,4)$.

Letting $\dot{q}$ denote either $\check{q}$ or $\hat{q}$, we can define $\dot{N}^\pm_k(x)$ and $\dot{Q}^\pm(x)$ as before and arrive at corresponding notions of firstlast and lastfirst sets.
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

We have implemented the firstlast and lastfirst procedures, together with minmax and maxmin, in the R package landmark [@]. For rank-distance--based procedures, the user can choose between $\check{q}$ and $\hat{q}$.
The landmark-generating procedures return the indices of the selected landmarks, optionally together with the sets of indices of the points in the cover set (ball or neighborhood) centered at each landmark.
In addition to the number of landmarks $n$ and either the radius $\eps$ of the balls or the cardinality $k$ of the neighborhoods on which they are based, the user may also specify additive and multiplicative extension factors for $n$ and for $\eps$ or $k$. These will produce additional landmarks ($n$) and larger cover sets ($\eps$ or $k$) with increased overlaps, in order (for example) to construct more robust nerve complexes.
We report validation and benchmark tests in Section\nbs\ref{sec:evaluations}.

## Validation

We validated the firstlast and lastfirst procedures against several small example data sets, including that of Example\nbs\ref{ex:rank-distance}.
We also implemented each of the maxmin and lastfirst procedures, using C++ (for Euclidean distances only) and R (which calls the proxy function to calculate distances other than Euclidean), and validated these against each other on several larger data sets, including as part of the benchmark tests reported in the next section.

## Benchmark tests

We benchmarked the C++ and R implementations on three data sets: uniform samples from the unit circle $\Sph^1\subset\R^2$ convoluted with Gaussian noise, samples with duplication from the integer lattice $[0,11]\times[0,5]$ using the probability mass function $p(a,b) \propto 2^{-ab}$, and patients recorded at each critical care unit in MIMIC-III using the RT-similarity measure.
We conducted benchmarks using the bench package [@] on HiPerGator.

\begin{figure}
\includegraphics[width=.6666667\textwidth]{../figures/benchmark-circle-lattice}
\includegraphics[width=.3333333\textwidth]{../figures/benchmark-mimic}
\caption{
Benchmark results for computing landmarks on two families of artificial data (circle and lattice) and one collection of empirical data (RT-similarity space of critical care units in MIMIC-III). Some points are missing because benchmark tests did not complete within 1 hour.
\label{fig:benchmark}
}
\end{figure}

Benchmark results are reported in Figure\nbs\ref{fig:benchmark}.
As expected (**explain why this is expected above**), R implementations used orders of magnitude more memory and took slightly longer. They appeared to scale slightly better in terms of time and slightly worse in terms of memory.
The additional calculations required for the lastfirst procedure increase runtimes by a median factor of 2.5 in our R implementations. The C++ implementation of lastfirst is based on combinatorial definitions and not optimized for speed, and as a result takes much longer---a median factor of almost 2000 relative to maxmin in C++---and failed to complete in many of our tests.

# Experiments

## Data sets

### Simulated data

_Choose one (at least moderately) large data set for each combination of the following properties: high- versus low-dimensional, with and without high multiplicities. Have the simulated data sets exceed/bookend the empirical data sets on both sides of each property range._

### Empirical data
\label{sec:empirical-data}

* [MIMIC-III Cardiac Care Unit](https://mimic.physionet.org/mimictables/transfers/) under one or both similarity measures
* [México COVID-19 Comunicado Técnico Diario](https://www.gob.mx/salud/documentos/coronavirus-covid-19-comunicado-tecnico-diario-238449) under one similarity measure

## Robustness

We compared the suitability of three landmarking procedures (random, maxmin, lastfirst) on datasets with varying density and duplication patterns by extending an example of @deSilva2004. Each expriment proceeded as follows: We sampled $n=540$ points from the sphere $\Sph^2\subset\R^3$ and selected $k=12$ landmark points. We then used the landmarks to compute PH and recorded the statistics $R_0,R_1,K_0,K_1$ as defined by @deSilva2004. From these statistics we compute the _relative dominance_ $(R_1 - R_0) / K_0$ and _absolute dominance_ $(R_1 - R_0) / K_1$ of the last interval over which all Betti numbers are correctly calculated.

The points $x=(r,\theta,\phi)$ were sampled using four procedures: uniform sampling, skewed sampling, uniform sampling with skewed boosting, and skewed sampling with skewed boosting. The first procedure was used by @deSilva2004 and here serves as a baseline case. For a sample $S$ (with multiplicities) generated from each of the other three procedures, the expected density $\lim_{\eps\to 0}\lim_{n\to\infty}\frac{1}{n}\abs{\{x=(r,\theta,\phi)\in S\mid \alpha-\eps<\phi<\alpha+\eps\}}$ of points near a given latitude $\alpha\in[0,\pi]$ is proportional to the quartic function $p:[0,1]\to[0,1]$ defined by $p(x)=(\frac{\phi}{\pi})^4$.
Skewed sampling is performed via rejection sampling: Points $x_i=(r_i,\theta_i,\phi_i)$ are sampled uniformly and rejected at random if a uniform random variable $t_i\in[0,1]$ satisfies $(\frac{\phi_i}{\pi})^\alpha<t_i$ until $n$ points have been kept [@Diaconis2013].
Skewed boosting is performed by first obtaining a (uniform or skewed) sample $T$ of size a fraction $\frac{n}{6}$ of the total, then sampling $n$ points (with replacement) from $T$ using the probability mass function satisfying $P(x_i)\propto(\frac{\phi_i}{\pi})^\beta$.
When performed separately, skewed sampling and skewed boosting use $\alpha=\beta=4$; when performed in sequence, they use $\alpha=\beta=2$.

The landmark points were selected in three ways: uniform random selection (without replacement), the maxmin procedure, and the lastfirst procedure.
We computed PH in Python GUDHI [@], using three implementations: Vietoris–Rips filtrations on the landmarks, alpha complexes on the landmarks, and (Euclidean) witness complexes on the landmarks with the full sample as witnesses.

## Covers and nerves

Cardinality reduction techniques can be used to model a large number of cases represented by a large number of variables as a smaller number of clusters with similarity or overlap relations among them.
The deterministic sampling procedures maxmin and lastfirst provide clusters (cover sets) defined by proximity to the landmark cases and relations defined by the sharing of cases.
The clusters obtained by these procedures occupy a middle ground between the regular tiles or quantiles commonly used to cover samples from Euclidean space [@] and the emergent clusters obtained heuristically by penalizing between-cluster similarity and rewarding within-cluster similarity [@].
The maxmin procedure produces cover sets of fixed radius, analogous to the congruent tiles of overlapping tessellations, while the lastfirst procedure produces cover sets of fixed size, analogous to the quantiles of an adaptive cover [@].
This makes them natural solutions to the task of covering an arbitrary finite metric space that may or may not contain important geometric or topological structure [@Singh2007].

As a practical test of this potential, we loosely followed the approach of @Dlotko2019 to construct covers and their nerves for two real-world clinical data sets (see Section\nbs\ref{sec:empirical-data}), using maxmin and lastfirst.
We varied the number of landmarks (12, 24, 60) and the multiplicative extension of the cover sets (0, .1, .25).
We evaluated the procedures in three ways:

- **Clustering quality:** While not designed for clustering, both procedures yield _fuzzy_ clusters---that is to say, clusters that allow for some overlap. We do not expect these coverings to be competitive with clustering methods, but we think it reasonable to consider the two basic measures of clustering quality---compactness and separation---in our comparisons. A significant hindrance is that most clustering validation measures, including almost all that have been proposed for fuzzy clusterings, rely not only on inter-point distances but on coordinate-wise calculations (specifically, data and cluster centroids) [@Bouguessa2006; @Wang2007; @Falasconi2010]. To our knowledge, the sole exception to have appeared in a comprehensive comparison of such measures is the _modified partition coefficient_ [@Dave1996], defined as $$\operatorname{MPC}=1-\frac{k}{k-1}(1-\frac{1}{n}\sum_{i=1}^{n}{\sum_{j=1}^{k}{{u_{ij}}^2}})$$ where $U=(u_{ij})$ is the $n\times k$ fuzzy partition matrix: $u_{ij}$ encodes the extent of membership of point $x_i$ in cluster $c_j$, and $\sum_{j=1}^{k}{u_{ij}}=1$ for all $i$. When a point $x_i$ is contained in $m$ cover sets $c_j$, we equally distribute its membership so that $u_{ij}=\frac{1}{m}$ when $x_i\in c_j$ and $u_{ij}=0$ otherwise. Like the partition coefficient from which it is adapted, the MPC takes the value $1$ on crisp partitions and is penalized by membership sharing. The MPC is standardized so that its range does not depend on $k$.
- **Discrimination of risk:** For purposes of clinical phenotyping, patient clusters are more useful that discriminate between low- and high-risk subgroups. We calculate a cover-based risk estimate from individual outcomes $y_i$ as follows: For each cover set $c_j\subset X$, let $p_j=\frac{1}{\abs{c_j}}\sum_{x_i\in c_j}{y_i}$ be the incidence of the outcome in that set. Then compute the weighted sum $q_i=\sum_{x_i\in c_j}{u_{ij}p_j}$ of these incidence rates for each case. We measure how well the cover discriminates risk as the _area under the receiver operating curve_ (AUC).
<!--
- **Homological richness:** Independent of the "true" topological structure of a data set, a simplicial complex model is more useful when it is not too sparse to be connected (or nearly so) and not too dense to visualize or contain low-degree cycles. We compare the $0$- and total $>0$-degree Betti numbers of the simplicial complex models from the perspective of optimizing $-\beta_0+\sum_{i>0}{\beta_i}$.
-->

\begin{figure}
\includegraphics[width=.5\textwidth]{../figures/cover-simplices}
\includegraphics[width=.5\textwidth]{../figures/cover-evaluate}
\caption{
Summary and evaluation statistics versus number of 0-simplices (landmarks) for the covers generated using the maxmin and lastfirst procedures, with three multiplicative extensions in their size.
Left: the sizes of their nerves, as numbers of 1- and 2-simplices.
Right: the modified partition coefficient (MPC) and the c-statistic of the risk prediction model based on the cover sets (AUC).
\label{fig:cover}
}
\end{figure}

## Interpolative prediction

Landmark points may also be used to trade accuracy for memory in neighborhood-based prediction modeling. Consider the following approach: Suppose that a modeling process involves predictor data $X \in \R^{n \times p}$ and response data $y \in \R^{n \times 1}$, partitioned into training and testing sets $X_0,X_1$ and $y_0,y_1$ according to a partition $I_0 \sqcup I_1 = \{1,\ldots,n\}$ of the index set. Given $x \in X_1$, a traditional nearest neighbors model computes a prediction $p(x) = \frac{1}{k}\sum_{q(x,x_i) \leq k}{y_i}$ by averaging the responses for the $k^\text{th}$ nearest neighbors of $x$ in $X_0$. By selecting a landmark set $L \subset X_0$, a researcher can reduce the computational cost of the model as follows: For each $\ell \in L$, calculate $p(\ell)$ as above. Then, for $x \in X_1$, calculate $p'(x) = \sum_{\ell \in L}{w(d(x,\ell)) p(\ell)} / \sum_{\ell \in L}{w(d(x,\ell))}$, where $w : \R_{\geq 0} \to \R_{\geq 0}$ is a weighting function (for example, $w(d)=d^{-1}$). The nearest neighbor predictions for $L$ thus serve as a proxy for the responses associated with $X_0$.

We took this approach to the prediction of in-hospital mortality for cardiac care patients recorded in MIMIC-III and of Covid-19 positivity for patients tested by the Mexican Department of Health.
In both cases, we represented the point cloud as a binary case-by-variable matrix $X \in \B^{n \times p}$, suitable for using the cosine distance metric, using a procedure adapted from that of @Zhong2020, designed for demographic and diagnostic EHR data (there termed relational–transaction data).
We then implemented the following procedure:

1. Determine a nested $6 \times 6$–fold split for train–fit–test cross-validation [@???]. That is, partition $[n] = \bigsqcup_{i=1}^{6}{I_i}$ into roughly equal parts, and partition each $[n] \wo I_i = \bigsqcup_{j=1}^{6}{J_{ij}}$ into roughly equal parts. **Use time series nested cross-validation for the CTD data?**
2. Iterate the following steps over each $i,j$:
    a) Generate a sequence $L$ of landmarks from the points $X_{([n] \wo I_i) \wo J_{ij}}$.
    b) Identify the $180$ nearest neighbors $N^+_{180}(\ell)$ of each landmark $\ell$.
    c) Find the value of $k \in [180]$ and the weighting function $w$ (among those available) for which the predictions $p' : X_{J_{ij}} \to [0,1]$ maximize the AUROC.
    d) Use the AUROC to evaluate the performance of the predictions $p' : X_{I_i} \to [0,1]$ using these $k$ and $w$.

We replicated the procedure for each combination of procedure (random, maxmin, lastfirst) and number of landmarks ($\abs{L}=30,60,120$). To account for randomness, we performed $12$ replications using random landmarking.

\begin{figure}
\includegraphics[width=\textwidth]{../figures/knn-auc}
\caption{
C-statistics of the interpolative predictive models based on covers constructed using random, maxmin, and lastfirst procedures to generate landmarks.
For comparison, **describe the `NA` c-statistics**
\label{fig:knn}
}
\end{figure}

\pagebreak

# References

<!--
# To generate the LaTeX file, execute the following:
pandoc fixed-adaptive.md \
  -s \
  --number-sections \
  --bibliography=../lastfirst.bib \
  -o fixed-adaptive.tex
# To generate the PDF directly, execute the following:
pandoc fixed-adaptive.md \
  -t latex \
  --number-sections \
  --bibliography=../lastfirst.bib \
  --citeproc \
  -o fixed-adaptive.pdf
-->
