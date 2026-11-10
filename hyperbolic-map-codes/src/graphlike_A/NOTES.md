# graphlike_A: exact distance of a graph-like side

Implementation A of track T-A (prompt Section 8, C1'). Program: `graphlike_A.c`. Everything in this file is **proved** unless labelled otherwise.

## Setting

Let $H_X\in\mathbb F_2^{m_X\times n}$, $H_Z\in\mathbb F_2^{m_Z\times n}$ with $H_XH_Z^{\mathsf T}=0$, $k=n-\operatorname{rk}H_X-\operatorname{rk}H_Z$, and
$d_Z=\min\{|c|: c\in\ker H_X\setminus\operatorname{rowsp}H_Z\}$, $d_X=\min\{|c|: c\in\ker H_Z\setminus\operatorname{rowsp}H_X\}$.
We treat $d_Z$ under the hypothesis that every column of $H_X$ has weight $\le2$; $d_X$ is the same statement with $X$ and $Z$ exchanged.

$\Gamma=\Gamma^+_X$ is the multigraph with vertex set $V=\{1,\dots,m_X\}\sqcup\{\infty\}$ and edge set $E=\{1,\dots,n\}$: edge $j$ joins the two checks containing qubit $j$, or that check and $\infty$ if column $j$ has weight 1, and is a loop at $\infty$ if column $j$ is zero. Let $\partial:\mathbb F_2^E\to\mathbb F_2^V$, $\partial(e)=$ sum of the two endpoints of $e$ (so $\partial(e)=0$ for a loop), and $Z(\Gamma)=\ker\partial$, the cycle space. A *circuit* is the edge set of a closed trail without repeated vertices: a loop, two parallel edges, or a cycle of length $\ge3$.

## Lemma 1 (kernel and cycle space)

$\ker H_X=Z(\Gamma)$.

*Proof.* $H_X$ is $\partial$ followed by the projection forgetting the coordinate $\infty$. Every column of $\partial$ has even weight (2 or 0), so the $\infty$-coordinate of $\partial c$ is the sum of the other coordinates. Hence $\partial c=0$ iff $H_Xc=0$. $\square$

## Lemma 2 (logical test)

Let $L\in\mathbb F_2^{k\times n}$ with $H_ZL^{\mathsf T}=0$ and $\operatorname{rk}\binom{H_X}{L}=\operatorname{rk}H_X+k$. Then $\operatorname{rowsp}H_X+\operatorname{rowsp}L=\ker H_Z$, and for $c\in\ker H_X$: $c\in\operatorname{rowsp}H_Z$ iff $Lc=0$.

*Proof.* $\operatorname{rowsp}H_X\subseteq\ker H_Z$ by $H_XH_Z^{\mathsf T}=0$, and $\operatorname{rowsp}L\subseteq\ker H_Z$ by hypothesis; the sum has dimension $\operatorname{rk}H_X+k=n-\operatorname{rk}H_Z=\dim\ker H_Z$. Over any field $(\ker M)^\perp=\operatorname{rowsp}M$, so $c\in\operatorname{rowsp}H_Z$ iff $c\perp\ker H_Z$ iff $c\perp\operatorname{rowsp}H_X$ and $Lc=0$; the first condition holds since $c\in\ker H_X$. $\square$

The program computes $L$ as follows: a basis of $\ker H_Z$ from the reduced echelon form of $H_Z$, each vector reduced against an echelon basis of $\operatorname{rowsp}H_X$ and the previously accepted vectors, the nonzero remainders kept. It then *certifies* the hypotheses of Lemma 2 at run time: commutation $H_XH_Z^{\mathsf T}=0$, exactly $k$ remainders, $H_ZL^{\mathsf T}=0$ (sparse parity checks), and $\operatorname{rk}\binom{H_X}{L}=\operatorname{rk}H_X+k$ by a fresh elimination (`lcheck` in the output).

## Lemma 3 (minimum-weight logicals are circuits)

(a) Every $c\in Z(\Gamma)$ is a disjoint union of circuits.
(b) If $k\ge1$, every $c\in\ker H_X\setminus\operatorname{rowsp}H_Z$ with $|c|=d_Z$ is a circuit.

*Proof.* (a) Induction on $|c|$. If $c$ contains a loop $f$, then $f$ is a circuit and $c+f\in Z(\Gamma)$. Otherwise every vertex of the subgraph $\Gamma[c]$ has even degree $\ge2$. Take a longest path $u_0,\dots,u_t$ in $\Gamma[c]$ ($t\ge1$); $u_t$ has an edge $f$ of $c$ other than the last path edge, and by maximality its other end is some $u_s$, $s<t$ ($f$ is not a loop). Then $C=\{u_su_{s+1},\dots,u_{t-1}u_t,f\}$ is a circuit and $c+C\in Z(\Gamma)$ is shorter.
(b) As $k\ge1$, $\dim\ker H_X=n-\operatorname{rk}H_X>\operatorname{rk}H_Z$, so $\ker H_X\supsetneq\operatorname{rowsp}H_Z$ and $d_Z$ is defined. Write $c=c_1+\dots+c_r$ as in (a). By Lemma 2, $0\ne Lc=\sum_iLc_i$, so $Lc_i\ne0$ for some $i$, and $|c_i|\le|c|$ with equality iff $r=1$. Minimality forces $r=1$. $\square$

## Candidates

Fix a total order on $V$ with $\infty$ first. For $v\in V$ let $\Gamma_v$ be the subgraph induced on $\{u\ge v\}$ without loops, $T_v$ a breadth-first-search tree of the component of $v$ in $\Gamma_v$, $\delta_v(x)$ the depth of $x$ and $P_v(x)\subseteq E$ the tree path from $v$ to $x$. For an edge $e=xy$ of $\Gamma_v$ with $x,y\in T_v$ put
$$C(v,e)=P_v(x)+e+P_v(y)\in\mathbb F_2^E,\qquad \ell(v,e)=\delta_v(x)+1+\delta_v(y).$$
Annotations: $a(e)=Le\in\mathbb F_2^k$ (column $e$ of $L$), $A_v(x)=LP_v(x)$, so $LC(v,e)=A_v(x)+a(e)+A_v(y)$ by linearity, and $A_v(y)=A_v(x)+a(e)$ along a tree edge $e=xy$.

## Lemma 4

$C(v,e)\in Z(\Gamma)$ and $|C(v,e)|\le\ell(v,e)$.

*Proof.* $\partial P_v(x)=v+x$, $\partial e=x+y$, $\partial P_v(y)=v+y$ sum to $0$; $|a+b+c|\le|a|+|b|+|c|$. $\square$

## Proposition 5 (exactness)

Assume every column of $H_X$ has weight $\le2$, $k\ge1$, and $L$ as in Lemma 2. If some loop $e$ of $\Gamma$ has $a(e)\ne0$, then $d_Z=1$. Otherwise
$$d_Z=\min\{\ell(v,e):\ v\in V,\ e\in E(\Gamma_v),\ \text{both ends of }e\text{ in }T_v,\ A_v(x)+a(e)+A_v(y)\ne0\}.$$

*Proof.* A loop $e$ lies in $\ker H_X$ (Lemma 1), so a loop with $a(e)=Le\ne0$ is a $Z$-logical of weight 1 (Lemma 2); conversely a weight-one element of $Z(\Gamma)$ is a loop.
($\ge$) By Lemma 4 and Lemma 2 a candidate with nonzero annotation is in $\ker H_X\setminus\operatorname{rowsp}H_Z$, so $\ell(v,e)\ge|C(v,e)|\ge d_Z$.
($\le$) By Lemma 3(b) take a circuit $c$ with $|c|=d=d_Z$ and $Lc\ne0$; if no loop is nontrivial, $c$ is not a loop. Let $v$ be the least vertex of $c$ and write $c$ as the closed walk $v=u_0,e_1,u_1,\dots,e_d,u_d=v$. All $u_i\ge v$, so $c\subseteq\Gamma_v$ and all $u_i$ lie in $T_v$. Telescoping, with $P_v(u_0)=P_v(u_d)=\emptyset$ and every other $P_v(u_i)$ occurring twice,
$$\sum_{i=1}^dC(v,e_i)=\sum_{i=1}^d\big(P_v(u_{i-1})+e_i+P_v(u_i)\big)=c .$$
Applying $L$, some $i$ has $LC(v,e_i)\ne0$. Walking along $c$ in both directions gives $\delta_v(u_{i-1})\le i-1$ and $\delta_v(u_i)\le d-i$, hence $\ell(v,e_i)\le d$. $\square$

*Remarks.* (1) No uniqueness of shortest paths and no minimum cycle basis is needed: the candidate found need not be $c$, only some nontrivial candidate of nominal length $\le d$. This replaces the proof sketch in C1' (Horton's theorem with perturbed lengths plus the exchange property; there the phrase "neither is such a $b$" should read "not every such $b$ lies in $\operatorname{rowsp}H_Z$"). The argument is the usual one for shortest odd cycles with $\mathbb F_2$ replaced by $\mathbb F_2^k$: $\Gamma$ with the labels $a(e)$ is an $\mathbb F_2^k$-labelled graph and $d_Z$ is its shortest cycle of nonzero label (literature on shortest non-zero cycles in group-labelled graphs not checked).
(2) Restricting to $\Gamma_v$ (least vertex of the circuit) is a speed-up beyond C1'; with $\Gamma_v=\Gamma$ the statement holds for every vertex $v$ of $c$.
(3) Witness: for the minimizing $(v,e)$ the program stores $C(v,e)$ when it is found (parent pointers of a BFS never change after discovery) and rechecks it: weight $=d_Z$ (a smaller weight would contradict ($\ge$)), $H_Xc=0$, $Lc\ne0$, and $c\notin\operatorname{rowsp}H_Z$ by reduction against the echelon form of $H_Z$, independently of $L$ (`ver` in the output).

## Corollary 6 (B4, item 1)

If $\Gamma^+_X$ is connected and $k\ge1$, then $d_Z\le2\operatorname{rad}(\Gamma^+_X)+1$.

*Proof.* Take $\Gamma_v=\Gamma$ and $v$ a centre. The telescoping identity holds for every closed walk (its base point need not be $v$), so the $C(v,e)$ span $Z(\Gamma)$ by Lemma 3(a). As $L\ne0$ on $Z(\Gamma)$, some $LC(v,e)\ne0$, and $|C(v,e)|\le\ell(v,e)\le2\operatorname{ecc}(v)+1$. $\square$

## Lemma 7 (bookkeeping of the search)

Let $D$ be the current best value and $R=\lfloor(D-1)/2\rfloor$.
(a) If $\ell(v,e)\le D-1$ then $\delta_v(x),\delta_v(y)\le R$, since $|\delta_v(x)-\delta_v(y)|\le1$ for an edge of a BFS and $\delta_v(x)+\delta_v(y)\le D-2$. So the BFS discovers vertices only up to depth $R$, and stops at the first popped vertex of depth $>R$.
(b) For a tree edge $C(v,e)=0$; it is skipped.
(c) A non-tree edge with both ends discovered is examined once, when its later-popped end is scanned. If $D$ (hence $R$) drops during the BFS and the BFS stops before that end is popped, then that end has depth $\ge R_{\rm new}+1$ and the edge has $\ell\ge2(R_{\rm new}+1)\ge D_{\rm new}$, so nothing is lost.
(d) If no loop is nontrivial, every candidate has $\ell\ge2$, so the search stops as soon as $D=2$.
(e) The pruned search returns $d_Z$. Every value it records is some $\ell(v,e)\ge d_Z$ (Proposition 5, $\ge$). Let $c$, $v$, $i$ be as in the proof of Proposition 5 ($\le$), $d=d_Z$. If $D>d$ when the BFS from $v$ starts, then $R\ge\lfloor d/2\rfloor$ until $D$ drops to $d$; the vertices of $c$ have depth $\le\lfloor d/2\rfloor$, so they are discovered (R never increases, so a vertex of depth $\le R$ is discovered by the first scanned neighbour one level up, as in an unbounded BFS), and the edge $e_i$ is examined by (c) unless $D=d$ has been reached before. $\square$

## Proposition 8 (number of minimum-weight logicals, option `--count`)

Under the hypotheses of Proposition 5, with $d=d_Z$: the number of $c\in\ker H_X\setminus\operatorname{rowsp}H_Z$ with $|c|=d$ equals the number of circuits of $\Gamma$ of length $d$ with $Lc\ne0$ (Lemma 3(b), and circuits are determined by their edge sets). For $d=1$ these are the nontrivial loops. For $d\ge2$ each such circuit with least vertex $v$ corresponds to exactly two closed walks $v=u_0,\dots,u_d=v$ in $\Gamma_v$ with $u_0,\dots,u_{d-1}$ distinct and pairwise distinct edges (the two orientations), and its vertices satisfy $\delta_v(u_i)\le\min(i,d-i)$. The program enumerates these walks by depth-first search with the pruning $t+\delta_v(u_t)\le d$ and halves the total; it aborts (count `null`) after a given number of search nodes. $\square$

## Complexity (proved bounds, measured times in RESULTS.md)

Linear algebra: $O(n^3/64)$ word operations (eliminations of $H_X$, $H_Z$, of $\binom{H_X}{L}$, and the reduction of $n-\operatorname{rk}H_Z$ kernel vectors). Search: at most $|V|$ breadth-first searches, each $O((|V|+|E|)\lceil k/64\rceil)$; with the pruning of Lemma 7 each search only explores the ball of radius $\lfloor(d_Z-1)/2\rfloor$ in $\Gamma_v$ once a candidate of length $d_Z$ is known. Counting is exponential in $d$ in general (surface codes with $d\ge20$ exceed $10^8$ nodes).

## References

- J. D. Horton, A polynomial-time algorithm to find the shortest cycle basis of a graph, SIAM J. Comput. 16 (1987), 358-366 (the candidate set; only the idea is used, the proof above is self-contained).
- E. Dennis, A. Kitaev, A. Landahl, J. Preskill, Topological quantum memory, arXiv:quant-ph/0110143 (surface-code distance as a shortest homologically nontrivial cycle; standard).
- C. Gidney, Stim: a fast stabilizer circuit simulator, arXiv:2103.02202: `shortest_graphlike_error` addresses the analogous problem for detector error models; its algorithm was not checked here.
