# graphlike_B: exact distance on a graph-like side

Track T-A, implementation B (2026-11-10). Program: `graphlike_B.c` (wrapper `graphlike_B`).
Principle: shortest paths in double covers of the check graph, one cover per logical operator of the
opposite type, all covers searched at once (one bit per cover), with a vertex-elimination order.
No shortest-path trees and no Horton candidate set are used.

Files and usage:

* `./graphlike_B [--ref] NAME.hx NAME.hz`: one JSON line `{"n","k","dX","dZ","witX","witZ","secs"}`; a side that
  is not graph-like (or $k = 0$) gets `null` and a field `reasonX`/`reasonZ`. `--ref` adds `refX`/`refZ` from the
  unpruned reference search and exits with code 3 on disagreement. Errors (bad input, non-CSS, failed witness
  recheck) exit with code 1 and a message on stderr. The wrapper compiles `bin/graphlike_B` on first use.
* `board2inst.py BOARD.json PREFIX`: board entry to `PREFIX.hx`, `PREFIX.hz` (plus `PREFIX.board.json` with the
  board's claims). `gen_instances.py`: toric, planar, rotated and random test codes.
* Cross-checks: `bruteforce.py`, `mcb_check.py`, `circuits_check.py` (section 7). Driver: `validate.py`,
  `run_validation.sh`; edge cases: `tests/edge/run_edge.sh`. Results: RESULTS.md.

## 0. Conventions

A CSS code is a pair $H_X \in \mathbb{F}_2^{m_X \times n}$, $H_Z \in \mathbb{F}_2^{m_Z \times n}$ with
$H_X H_Z^{T} = 0$. Vectors $x \in \mathbb{F}_2^n$ are identified with subsets of $\{0,\dots,n-1\}$,
$|x|$ is the Hamming weight and $\langle x, y\rangle = \sum_e x_e y_e$.

$$d_Z = \min\{|x| : x \in \ker H_X \setminus \operatorname{rowsp} H_Z\},\qquad
  d_X = \min\{|x| : x \in \ker H_Z \setminus \operatorname{rowsp} H_X\},\qquad
  k = n - \operatorname{rk} H_X - \operatorname{rk} H_Z .$$

The program works with a pair $(A, B)$ with $A B^{T} = 0$ and computes
$$\delta(A,B) = \min\{|x| : x \in \ker A \setminus \operatorname{rowsp} B\}\qquad(\min\emptyset = \infty)$$
whenever every column of $A$ has weight at most 2. Then $d_Z = \delta(H_X, H_Z)$ and $d_X = \delta(H_Z, H_X)$.
Note $\operatorname{rowsp} B \subseteq \ker A$, because $A B^{T} = 0$.

## 1. The graph of $A$

**Definition 1.** Let $A$ have rows $0,\dots,m-1$ and columns of weight at most 2. The multigraph $G = G(A)$
has vertex set $V = \{0,\dots,m-1\} \cup \{\infty\}$ and edge set $E = \{0,\dots,n-1\}$; the edge $e$ joins
the two rows containing column $e$ if the column has weight 2, joins its row to $\infty$ if it has weight 1,
and is a loop at $\infty$ if it has weight 0. For $x \subseteq E$ and $v \in V$ let $\deg_x(v)$ be the
number of edge ends of $x$ at $v$ (a loop contributes 2). The cycle space is
$\mathcal{Z}(G) = \{x \subseteq E : \deg_x(v) \text{ even for all } v\}$.

**Lemma 1.** $\ker A = \mathcal{Z}(G)$.

*Proof.* For a row $i$ there is no loop at $i$, so $(Ax)_i = \sum_{e \ni i} x_e \equiv \deg_x(i) \pmod 2$.
Hence $x \in \mathcal{Z}(G)$ implies $Ax = 0$. Conversely, if $Ax = 0$ then all $\deg_x(i)$ with
$i \neq \infty$ are even, and $\deg_x(\infty) = 2|x| - \sum_{i \neq \infty} \deg_x(i)$ is even (handshake). $\square$

## 2. The nontriviality test

**Lemma 2.** Let $A B^{T} = 0$ and let $y_1, \dots, y_k \in \ker B$ with
$\operatorname{rowsp} A + \langle y_1,\dots,y_k\rangle = \ker B$. For $x \in \ker A$:
$x \in \operatorname{rowsp} B \iff \langle x, y_j\rangle = 0$ for all $j$.

*Proof.* The standard form on $\mathbb{F}_2^n$ is nondegenerate, so $U^{\perp\perp} = U$ for every subspace $U$;
with $\ker B = (\operatorname{rowsp} B)^{\perp}$ this gives $\operatorname{rowsp} B = (\ker B)^{\perp}$.
So $x \in \operatorname{rowsp} B$ iff $x \perp \ker B$. Since $x \in \ker A = (\operatorname{rowsp} A)^{\perp}$,
$x$ is orthogonal to $\operatorname{rowsp} A$, and $x \perp \ker B$ iff $x \perp y_j$ for all $j$. $\square$

**Lemma 3 (computation of the $y_j$, function `logical_basis`).** The program (i) inserts the rows of $A$ into
an echelon table (pivot = lowest set bit; inserting means reducing against the table and storing a nonzero
remainder), (ii) computes the basis $z_1,\dots,z_N$ of $\ker B$ attached to the free columns of the reduced
row echelon form of $B$, $N = n - \operatorname{rk} B$, (iii) inserts $z_1,\dots,z_N$ in turn and records
every nonzero remainder as some $y_j$. Then every $y_j \in \ker B$,
$\operatorname{rowsp} A + \langle y_j\rangle = \ker B$, and the number of $y_j$ is
$n - \operatorname{rk} A - \operatorname{rk} B$.

*Proof.* Rows of $A$ lie in $\ker B$ (as $B A^{T} = (A B^{T})^{T} = 0$), and by induction every stored
vector lies in $\ker B$, since a remainder is $z_t$ plus a sum of stored vectors. After step $t$ the table
spans $\operatorname{rowsp} A + \langle z_1,\dots,z_t\rangle$; at the end this is $\ker B$. Stored vectors
have distinct pivots, hence are independent; the table has $\dim \ker B = n - \operatorname{rk} B$ elements,
$\operatorname{rk} A$ of which come from step (i). The free-column vectors $z_f$ ($z_f$ has a 1 in the free
column $f$ and $R_{i f}$ in the pivot column of row $i$ of the RREF $R$) satisfy $R z_f = R_{if} + R_{if} = 0$
row by row and are independent (distinct free columns), so they form a basis of $\ker B$. $\square$

At run time the program checks $B y_j = 0$ for every $j$ and that the count equals
$n - \operatorname{rk} A - \operatorname{rk} B$ with both ranks computed by a separate elimination
(`rank_rows`).

**Corollary 4.** Put $S_j = \operatorname{supp} y_j \subseteq E$ and, for $S \subseteq E$,
$$D(S) = \min\{|x| : x \in \mathcal{Z}(G),\ |x \cap S| \text{ odd}\}.$$
Then $\delta(A,B) = \min_j D(S_j)$, and $\delta(A,B) = \infty$ if $k = 0$.

*Proof.* $\langle x, y_j\rangle \equiv |x \cap S_j|$. By Lemmas 1 and 2, $x \in \ker A \setminus \operatorname{rowsp} B$
iff $x \in \mathcal{Z}(G)$ and $|x \cap S_j|$ is odd for some $j$. Take minima. $\square$

## 3. Double covers

A walk $W = (v_0, e_1, v_1, \dots, e_L, v_L)$ has length $L$ and $S$-parity
$\pi_S(W) = \#\{i : e_i \in S\} \bmod 2$; it is closed if $v_L = v_0$.

**Definition 5.** For $S \subseteq E$ the double cover $\widetilde{G}_S$ has vertex set $V \times \mathbb{F}_2$;
every edge $e$ with ends $u, w$ (possibly $u = w$) has the two lifts joining $(u,a)$ and $(w, a + [e \in S])$,
$a \in \mathbb{F}_2$. The projection $(u,a) \mapsto u$ maps walks to walks of the same length.

**Lemma 6 (unique lifting).** For a walk $W$ in $G$ starting at $v_0$ and $a_0 \in \mathbb{F}_2$ there is a
unique walk in $\widetilde{G}_S$ starting at $(v_0, a_0)$ that projects to $W$; it ends at
$(v_L, a_0 + \pi_S(W))$. Every walk in $\widetilde{G}_S$ is the lift of its projection.

*Proof.* Induction on $L$: the traversal of $e_i$ from $v_{i-1}$ to $v_i$ has exactly one lift starting at
$(v_{i-1}, a_{i-1})$, and it ends at $(v_i, a_{i-1} + [e_i \in S])$. $\square$

So walks from $(v,0)$ to $(v,1)$ in $\widetilde{G}_S$ correspond, length-preservingly, to closed walks at $v$
in $G$ with odd $S$-parity.

**Lemma 7 (walks and cycles).**
(a) Let $W$ be a closed walk of length $L$ with $\pi_S(W) = 1$ and let $x(W)$ be the set of edges traversed an odd
number of times. Then $x(W) \in \mathcal{Z}(G)$, $|x(W)| \le L$ and $|x(W) \cap S|$ is odd.
(b) Every $x \in \mathcal{Z}(G)$ with $|x \cap S|$ odd contains a circuit $C$ (a closed walk without repeated
vertices other than the endpoint; loops and pairs of parallel edges included) with $|C \cap S|$ odd, $|C| \le |x|$.

*Proof.* (a) Let $\mu(e)$ be the number of traversals of $e$ and $\epsilon(e,v) \in \{0,1,2\}$ the number of ends
of $e$ at $v$. In a closed walk every visit of $v$ uses two edge ends, so $\sum_e \mu(e)\epsilon(e,v)$ is even, and
$\deg_{x(W)}(v) \equiv \sum_e \mu(e)\epsilon(e,v) \pmod 2$. Further $|x(W)| \le \sum_e \mu(e) = L$ and
$|x(W) \cap S| \equiv \sum_{e \in S} \mu(e) = \pi_S(W) \pmod 2$.
(b) Veblen's theorem: the edge set of a graph with all degrees even is a disjoint union of circuits (if $x \ne \emptyset$,
either $x$ contains a loop, or every vertex of positive degree has degree $\ge 2$ and $x$ contains a circuit;
removing it keeps all degrees even; induct on $|x|$). As $|x \cap S| = \sum_i |C_i \cap S|$ is odd, some $C_i$
has odd intersection, and $|C_i| \le |x|$. $\square$

**Proposition 8.** $D(S) = \min_{v \in V} \operatorname{dist}_{\widetilde{G}_S}((v,0),(v,1))$.

*Proof.* "$\le$": a shortest path from $(v,0)$ to $(v,1)$ of length $L$ projects to a closed walk with odd
$S$-parity (Lemma 6), so Lemma 7(a) gives $x \in \mathcal{Z}(G)$ with $|x \cap S|$ odd and $|x| \le L$.
"$\ge$": let $x$ realise $D(S)$; Lemma 7(b) gives a circuit $C$ with odd parity and $|C| \le D(S)$; traversing $C$
from one of its vertices $v$ and lifting from $(v,0)$ ends at $(v,1)$ (Lemma 6). $\square$

**Proposition 9 (vertex elimination).** Let $v_1,\dots,v_N$ be any ordering of $V$, let $G_i$ be the subgraph
induced on $\{v_i,\dots,v_N\}$ and $\widetilde{G}_{S,i}$ its double cover. Then
$$D(S) = \min_i \operatorname{dist}_{\widetilde{G}_{S,i}}((v_i,0),(v_i,1)).$$

*Proof.* A walk in $\widetilde{G}_{S,i}$ is a walk in $\widetilde{G}_S$, so each term is $\ge D(S)$ by
Proposition 8. Conversely take $C$ as in the proof of Proposition 8 and let $i$ be minimal with $v_i \in C$.
Then $C \subseteq G_i$, and the lift of the traversal of $C$ from $v_i$ is a walk of length $|C| \le D(S)$ from
$(v_i,0)$ to $(v_i,1)$ in $\widetilde{G}_{S,i}$. $\square$

Together with Corollary 4:
$$\delta(A,B) = \min_{i}\ \min_{j}\ \delta_{ij},\qquad
  \delta_{ij} = \operatorname{dist}_{\widetilde{G}_{S_j,i}}((v_i,0),(v_i,1)).\tag{$*$}$$
The program orders $V$ by decreasing degree (any order is correct; removing high-degree vertices first makes
later searches cheaper).

## 4. The bit-sliced search (function `search_min`)

For each root $v_i$ the program runs the level-synchronous breadth-first search from $(v_i,0)$ in all $k$ covers
$\widetilde{G}_{S_j,i}$ simultaneously. For $u \in V$ and $a \in \mathbb{F}_2$ it keeps bit vectors
$R_a[u], F_a[u] \in \mathbb{F}_2^k$, and $s_e \in \mathbb{F}_2^k$ with $(s_e)_j = [e \in S_j]$.

**Lemma 10 (invariant).** After level $t$, bit $j$ of $R_a[u]$ is $[\operatorname{dist}_j((v_i,0),(u,a)) \le t]$
and bit $j$ of $F_a[u]$ is $[\operatorname{dist}_j((v_i,0),(u,a)) = t]$, distances taken in $\widetilde{G}_{S_j,i}$.

*Proof.* Level 0: $R_0[v_i] = F_0[v_i] = \mathbf{1}$, everything else 0. Step $t \to t+1$: the program forms, over the
adjacency entries "$u$ to $w$ through $e$" with $u$ in the frontier and $w$ not deleted,
$$N_0[w] \mathrel{|}= (F_0[u] \wedge \neg s_e) \vee (F_1[u] \wedge s_e),\qquad
  N_1[w] \mathrel{|}= (F_1[u] \wedge \neg s_e) \vee (F_0[u] \wedge s_e),$$
so bit $j$ of $N_a[w]$ says that some lift of an edge at $w$ joins $(w,a)$ to a vertex $(u, a + (s_e)_j)$ at
distance exactly $t$ (each non-loop edge is stored in both directions, a loop once, which suffices because the
lifts of a loop at $u$ join $(u,b)$ and $(u, b + (s_e)_j)$ symmetrically). Then $F' = N \wedge \neg R$ and
$R' = R \vee F'$: a vertex is at distance exactly $t+1$ iff it is not within distance $t$ and has a neighbour at
distance exactly $t$. $\square$

**Pruning.** With the current minimum `best` (initially $\infty$) the search from $v_i$ stops at the first level
$t$ with $R_1[v_i] \ne 0$ (then $\min_j \delta_{ij} = t$ by Lemma 10, and `best` becomes $t$ since only levels
$t < $ `best` are run), or when $t$ reaches `best` (then $\min_j \delta_{ij} \ge$ `best`), or when the frontier
is empty (then $(v_i,1)$ is unreachable in every cover). Hence after the last root `best` equals the right-hand
side of $(*)$, that is $\delta(A,B)$. If $k > 0$ then $\delta(A,B) < \infty$, and the program aborts if the search
returns $\infty$.

## 5. Witness and independent recheck

Let $(i^*, j^*)$ be the pair at which `best` was last lowered, $t^* = \delta(A,B)$. The program runs an ordinary
breadth-first search in $\widetilde{G}_{S_{j^*},i^*}$ with parent pointers (function `cover_bfs`) and asserts that
the path to $(v_{i^*},1)$ has length $t^*$. By Lemmas 6 and 7(a) the set $x$ of edges used an odd number of times
lies in $\ker A$, has $\langle x, y_{j^*}\rangle = 1$, hence $x \notin \operatorname{rowsp} B$ (Lemma 2), and
$|x| \le t^*$. Since $x \in \ker A \setminus \operatorname{rowsp} B$, also $|x| \ge \delta(A,B) = t^*$, so
$|x| = t^*$.

The recheck does not use the $y_j$: it counts $|x|$, evaluates every row of $A$ on $x$, and compares
$\operatorname{rk}\binom{B}{x}$ with $\operatorname{rk} B$ by a separate forward elimination (`rank_rows`).
Any failure aborts with a message. The pairing $\langle x, y_{j^*}\rangle = 1$ is asserted as well.

## 6. Cost

Linear algebra: $O(n^2 (m_X + m_Z)/64)$ word operations. Search: at most $\min(\delta, \text{eccentricity})$
levels per root, each at most $2|E|\lceil k/64\rceil$ word operations, so $O(|V|\,\delta\,|E|\,k/64)$ in total.
Memory: $O(|V| k/64 + n(m_X+m_Z+n)/64)$ words. Observed (RESULTS.md): the nine board instances of the brief take
at most 0.06 s each, instances with $n \approx 2000$ (toric $L = 31$, $d = 31$; random, $k = 300$) at most 0.08 s.

## 7. Cross-checks based on other principles

* `--ref` (in `graphlike_B.c`): for every $j$ and every vertex an unpruned scalar breadth-first search in the full
  cover $\widetilde{G}_{S_j}$ (Proposition 8 directly, no bit slicing, no elimination order, no pruning).
  The program exits with code 3 if it disagrees with the main search.
* `bruteforce.py`: the definition. Kernel basis by reducing $[A^{T} \mid I]$, enumeration of all of $\ker A$
  in Gray-code order, membership in $\operatorname{rowsp} B$ by elimination. No graph, no logical basis.
  Feasible for $\dim \ker A \le 24$.
* `mcb_check.py`: a minimum cycle basis computed by networkx (de Pina's algorithm) and Lemma 11. Edges are
  subdivided so that the networkx graph is simple and the edge set of a basis cycle can be read off from its
  subdivision vertices; the script checks that the recovered vectors are cycles and form a basis of
  $\mathcal{Z}(G)$ (count and rank). Minimality of the basis is taken from networkx. Shares no code with
  `graphlike_B.c`.

**Lemma 11 (a minimum basis contains a minimum nontrivial element).** Let $Z$ be a finite-dimensional
$\mathbb{F}_2$-vector space, $w : Z \to \mathbb{R}$ any function, $U \subsetneq Z$ a subspace and $\mathcal{B}$ a
basis of $Z$ minimising $\sum_{b \in \mathcal{B}} w(b)$. Then
$\min\{w(z) : z \in Z \setminus U\} = \min\{w(b) : b \in \mathcal{B} \setminus U\}$.

*Proof.* $\mathcal{B} \setminus U \neq \emptyset$ because $\mathcal{B}$ spans $Z \neq U$, and "$\le$" holds because
$\mathcal{B} \setminus U \subseteq Z \setminus U$. Let $z \in Z \setminus U$ realise the minimum and write
$z = \sum_{b \in T} b$ with $T \subseteq \mathcal{B}$. If some $b_0 \in T$ had $w(b_0) > w(z)$, then
$(\mathcal{B} \setminus \{b_0\}) \cup \{z\}$ would be a basis (it spans $b_0 = z + \sum_{T \setminus b_0} b$) of
smaller total weight. So $w(b) \le w(z)$ for all $b \in T$, and $T \not\subseteq U$ since $z \notin U$. $\square$

With $Z = \mathcal{Z}(G) = \ker A$, $w = |\cdot|$ and $U = \operatorname{rowsp} B$ this gives
$\delta(A,B) = \min\{|C| : C \in \mathcal{B},\ C \notin \operatorname{rowsp} B\}$.

* `circuits_check.py`: exhaustive enumeration of all circuits of length at most a bound $\ell$, each tested for
  membership in $\operatorname{rowsp} B$ (fully reduced echelon basis of $B$; $v$ lies in the span iff $v$ equals the
  sum of the basis vectors at the pivot positions in $\operatorname{supp} v$). By Lemma 12, running it with
  $\ell = d$ proves $\delta(A,B) \ge d$ when no nontrivial circuit of length $< d$ exists, finds a nontrivial
  circuit of length $d$, and counts the minimum-weight logical operators of this type. No covers, no logical basis.

**Lemma 12.** Every $x \in \ker A \setminus \operatorname{rowsp} B$ with $|x| = \delta(A,B)$ is a circuit. Hence
$\delta(A,B) = \min\{|C| : C \text{ a circuit of } G,\ C \notin \operatorname{rowsp} B\}$, and the minimum-weight
elements of $\ker A \setminus \operatorname{rowsp} B$ are exactly the nontrivial circuits of length $\delta(A,B)$.

*Proof.* By Lemma 1 and Veblen's theorem (proof of Lemma 7(b)), $x = C_1 \sqcup \dots \sqcup C_r$ with circuits $C_i$,
each in $\ker A$. As $x \notin \operatorname{rowsp} B$, some $C_i \notin \operatorname{rowsp} B$, so
$|C_i| \ge \delta(A,B) = |x| = \sum_j |C_j|$, whence $r = 1$. The converse is clear. $\square$

*Enumeration.* A circuit of length $L \ge 2$ has a unique smallest vertex $s$ and exactly two edges at $s$ (for
$L = 2$: two parallel edges). The depth-first search from $s$ through vertices larger than $s$, extending simple paths
to at most $\ell - 1$ edges and closing with an unused edge back to $s$, meets it exactly twice, once per direction;
it is counted only when the closing edge has a larger index than the first edge. Loops ($L = 1$) are listed directly.

## 8. References (bibliographic details from memory, not checked at the source)

* J. C. de Pina, *Applications of shortest path methods*, PhD thesis, University of Amsterdam, 1995.
  The double-cover ("signed graph") shortest odd cycle subroutine of section 3 is the one used there.
* T. Kavitha, K. Mehlhorn, D. Michail, K. E. Paluch, *An $\tilde O(m^2 n)$ algorithm for minimum cycle basis of
  graphs*, Algorithmica 52 (2008); cited by the networkx implementation used in `mcb_check.py`.
* J. D. Horton, *A polynomial-time algorithm to find the shortest cycle basis of a graph*, SIAM J. Comput. 16 (1987).
  (Not used; mentioned because the brief asks for a different principle.)
* O. Veblen, *An application of modular equations in analysis situs*, Ann. of Math. (2) 14 (1912/13).
  (Decomposition of even graphs into circuits, Lemma 7(b); the short proof is given above.)
