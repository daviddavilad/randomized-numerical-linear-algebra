# Randomized algorithms for low-rank matrix and tensor decompositions

**Pearce, K. J., & Martinsson, P.-G. (2026).** Randomized algorithms for low-rank matrix and tensor decompositions. In *Numerical Analysis and Scientific Computing: an Anthology* (pp. 197–252). Springer Nature Switzerland. [doi](https://doi.org/10.1007/978-3-032-33428-2_7) · [arXiv preprint](https://arxiv.org/abs/2512.05286)

Read: 2026-09-12 (pass 2) (Skipped section 4 since it covers randomized algorithms for low-rank matrix decompositions, can implement in the future if interested. It's basically a "how to carry out the different problems presented earlier in the reading" section) (Also skipped section 5 & 6 on tensors, because I do not quite understand them yet and want to give them their due time)

## Problem

The survey synthesizes the developments of numerical linear algebra (NLA) methods for low-rank decompositions of matrices and tensors. Methods covered include randomized dimension reduction, such as the well-known singular value decomposition (SVD) and interpolative (ID) and CUR decompositions.

Other approaches covered include new methods for fast matrix sketching and sampling techniques, along with the extension of all these methods applied to tensors.

Throughout its history, the main focus of the NLA field has been driven by two directions mainly:

- Accuracy, trying to keep the approximated matrix as close as possible to the original matrix in terms of information captured.

- Efficiency, trying to make those computations faster by using several approximation methodologies.

The textbook approach then becomes trying to find the right balance between both. Making an algorithm computationally cheap, while giving up the least accuracy possible.

Randomization was introduced into NLA (thus naming it RNLA) because it gives us the ability to reduce the computational cost. By randomizing the SVD problem (rSVD), we are able to capture the majority of the information about the matrix. Recall the classical SVD algorithm:

$$A = U\Sigma V^T,$$

We can throw random directions through A to capture most of the information about A in the following way:

Generate a random matrix

$$\Omega \in \mathbb{R}^{n \times \ell},$$

where usually

$$\ell = k + p$$

and $p$ is a small oversampling amount, perhaps $5$-$10$. Then compute

$$\boxed{Y = A\Omega}.$$

The columns of $Y$ span approximately the right subspace, but they're messy and not orthogonal. So perform QR:

$$Y = QR.$$

Keep

$$\boxed{Q = \mathrm{orth}(Y)}.$$

Now compute

$$\boxed{B = Q^T A}.$$

Dimensions:

$$\underbrace{Q^T}_{\ell \times m}\ \underbrace{A}_{m \times n} = \underbrace{B}_{\ell \times n}.$$

Compute

$$B = \widehat{U}\Sigma V^T.$$

Because $B$ is small, this SVD is much cheaper.

As a result, by introducing randomness, we are able to formulate a much smaller and computationally efficient problem.

## Approach

The paper first reviews the classical numerical linear algebra decompositions that later serve as building blocks for randomized algorithms. It then introduces randomized dimensionality reduction through **sketching** and **sampling**, whose goal is to replace a large matrix problem with a smaller problem that approximately preserves the important geometry of the original matrix.

### Singular Value Decomposition

For a matrix

$$A \in \mathbb{F}^{m \times n},$$

the singular value decomposition is

$$A = U \Sigma V^*,$$

where $U$ and $V$ have orthonormal columns and

$$\Sigma = \mathrm{diag}(\sigma_1,\ldots,\sigma_r), \qquad \sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_r \geq 0.$$

The rank-$k$ truncated SVD is

$$A_k = \sum_{i=1}^{k} \sigma_i u_i v_i^*.$$

By the Eckart--Young theorem, $`A_k`$ is the optimal rank-$k$ approximation to $A$. Its errors under the spectral and Frobenius norms are

$$\|A-A_k\|_2 = \sigma_{k+1},$$

and

$$\|A-A_k\|_F = \left(\sum_{j=k+1}^{r} \sigma_j^2\right)^{1/2}.$$

Thus, the singular values directly quantify how much information is lost when the matrix is compressed to rank $k$.

---

### QR Decomposition

Every matrix admits a factorization

$$A = QR,$$

where $Q$ has orthonormal columns and $R$ is upper triangular.

The paper is particularly interested in **column-pivoted QR (CPQR)**,

$$AP = QR,$$

where $P$ permutes the columns of $A$.

The purpose of pivoting is to move informative or linearly independent columns toward the front of the matrix, making it possible to truncate the decomposition after $k$ steps and obtain a useful low-rank approximation.

The columns of $Q$ form an orthonormal basis for the range of $A$. The paper therefore uses the notation

$$Q = \mathrm{col}(A)$$

when only this basis is required.

CPQR is also useful for selecting skeleton columns for ID and CUR decompositions. At each step, it greedily chooses the remaining column with the largest Euclidean norm in the active submatrix.

---

### LU Decomposition

The LU decomposition is written as

$$PA = LU,$$

where $P$ is a permutation matrix, $L$ is lower triangular, and $U$ is upper triangular.

As with QR, the factorization can be stopped after $k$ steps to produce a partial LU decomposition. Pivoting determines which rows are selected during the process.

LU therefore provides another method for identifying informative rows or columns that may subsequently be used as skeletons in low-rank decompositions.

---

### Interpolative and CUR Decompositions

Unlike the SVD, which represents a matrix through new orthogonal basis vectors, the **Interpolative Decomposition (ID)** and **CUR decomposition** approximate the matrix using actual rows and/or columns of $A$.

These selected rows and columns are called **skeletons**.

#### Interpolative Decomposition

Let

$$R = A(I,:)$$

contain selected rows of $A$, where $I$ is the corresponding row index set. A row ID has the form

$$A \approx AR^\dagger R.$$

Similarly, let

$$C = A(:,J)$$

contain selected columns indexed by $J$. A column ID has the form

$$A \approx CC^\dagger A.$$

If $A$ has exactly rank $k$ and the selected $k$ rows or columns span the corresponding row or column space, these become exact decompositions.

For low-rank approximation, the quality of the ID depends strongly on the choice of skeletons. The paper gives the guarantee

$$\|A - CC^\dagger A\|_F \leq \sqrt{k+1}\,\|A-A_k\|_F$$

for an appropriately constructed column ID.

Skeletons can be selected through pivoting methods based on QR, LU, or the SVD.

One method discussed is **Osinsky's method**, which starts from one rank-$k$ truncated SVD

$$A_k = U \Sigma V^*.$$

It uses information contained in the dominant right-singular subspace $V$ to iteratively select columns of $A$.

Defining

$$\widetilde{A} = A - AVV^*,$$

the method treats $\widetilde{A}$ as the residual left outside the dominant rank-$k$ subspace and chooses skeletons so that the resulting interpolation error remains controlled.

The full iterative derivation is deferred to a later pass.

#### CUR Decomposition

The CUR decomposition simultaneously selects rows and columns:

$$A \approx CUR,$$

where

$$C = A(:,J), \qquad R = A(I,:).$$

For fixed skeletons $C$ and $R$, the optimal middle matrix in Frobenius norm is

$$U = C^\dagger A R^\dagger.$$

Thus,

$$A \approx C\left(C^\dagger A R^\dagger\right)R.$$

CUR is attractive because the approximation remains expressed in terms of actual rows and columns of the original data, which can make the factors more interpretable than singular vectors.

However, directly forming the optimal $U$ still requires accessing the full matrix $A$, motivating faster randomized methods for selecting and manipulating the skeletons.

---

### Randomized Dimensionality Reduction

The authors then introduce **dimensionality reduction maps (DRMs)**. These maps transform a large row or column space into a much smaller space while attempting to preserve its geometry.

For

$$A \in \mathbb{F}^{m \times n},$$

a map

$$\Gamma : \mathbb{F}^{m} \rightarrow \mathbb{F}^{d}, \qquad k \leq d \ll m,$$

is a subspace embedding with distortion $\varepsilon$ if

$$(1-\varepsilon)\|Ax\| \leq \|\Gamma Ax\| \leq (1+\varepsilon)\|Ax\|$$

for every

$$x \in \mathbb{F}^{n}.$$

Thus, although $\Gamma A$ is much smaller than $A$, distances within the relevant subspace are approximately preserved.

If $\Gamma$ is drawn randomly from a probability distribution that satisfies this property with high probability, it is a **randomized dimensionality reduction map**.

The paper distinguishes two major approaches:

$$\boxed{\text{Sketching: mix coordinates}}$$

and

$$\boxed{\text{Sampling: select coordinates}}.$$

---

#### Randomized Matrix Sketching

A row sketch is formed as

$$\Gamma A,$$

while a column sketch can be formed as

$$A\Omega.$$

Rather than selecting individual rows or columns, sketching generally forms **random linear combinations** of them. The goal is to compress the matrix while approximately retaining the geometry of its important subspaces.

The survey discusses several important sketching matrices.

**Gaussian sketch.**

The entries of

$$\Gamma \in \mathbb{F}^{d \times m}$$

are independently sampled as

$$\Gamma_{ij} \sim \mathcal{N}\left(0,\frac{1}{d}\right).$$

Gaussian sketches are mathematically convenient and provide strong theoretical guarantees, although multiplying by a dense Gaussian matrix may itself be expensive.

**Subsampled Randomized Trigonometric Transform (SRTT).**

These methods combine random sign changes, permutations, Fourier-type transforms, and subsampling. A representative construction is

$$\Gamma = \sqrt{\frac{m}{d}}\,\Pi_{m \rightarrow d}F\Phi\Pi_m.$$

The purpose is to first mix information throughout the matrix and then safely subsample it.

**Sparse sign matrices.**

These replace dense Gaussian matrices with sparse random matrices whose nonzero entries are typically Rademacher random variables,

$$\rho_i \in \{-1,+1\}.$$

Their sparsity can substantially reduce the cost of forming the sketch.

**CountSketch.**

CountSketch is an extremely sparse special case of a sparse sign matrix in which each column contains only one nonzero random signed entry.

This makes sketching particularly inexpensive, although a larger embedding dimension can be required.

**SparseStack.**

SparseStack combines several extremely sparse randomized embeddings.

The survey reports that recent variants can obtain accuracy comparable to Gaussian sketches while retaining sparse and computationally efficient matrix operations.

**Sparse randomized trigonometric transforms (SparseRTT).**

These combine sparse random transformations with fast trigonometric transforms, attempting to obtain both good mixing and computational efficiency.

The central trade-off among these methods is therefore the cost of constructing and applying the sketch versus the embedding dimension $d$ required to retain the important geometry of $A$.

---

#### Randomized Matrix Sampling

Randomized sampling differs from sketching because it selects **actual rows or columns** rather than random linear combinations of them.

Suppose rows of

$$A \in \mathbb{F}^{m \times n}$$

are sampled according to probabilities

$$p_1,\ldots,p_m.$$

The paper discusses several probability distributions.

**Uniform sampling**

$$p_i = \frac{1}{m}.$$

This is essentially free to construct, but it works best when information is distributed relatively evenly throughout the matrix.

Its success depends on the **coherence** of the matrix. High coherence indicates that important information is concentrated in only a small number of coordinates, making uniform sampling more likely to miss important rows or columns.

---

**Squared-norm sampling**

Rows can instead be sampled according to

$$p_i = \frac{\|a_i\|_2^2}{\|A\|_F^2},$$

where $`a_i^*`$ is the $i$-th row of $A$.

Rows with larger Euclidean norms are therefore assigned greater sampling probability. This incorporates information about the matrix into the sampling distribution rather than treating every row equally.

Computing all of the row norms costs approximately

$$O(mn).$$

---

**Leverage-score sampling**

Another approach is to sample according to the statistical leverage scores,

$$p_i = \frac{\ell_i}{n},$$

where

$$\ell_i = \|a_i\|_{(A^*A)^{-1}}^2.$$

Equivalently,

$$\ell_i = a_i^*(A^*A)^{-1}a_i.$$

Leverage scores measure how strongly an individual row contributes to the geometry of the column space, rather than measuring only its raw Euclidean magnitude.

This gives stronger theoretical embedding guarantees than uniform or squared-norm sampling in many settings.

However, computing the exact leverage scores can be as expensive as computing

$$\mathrm{svd}(A),$$

which largely defeats their purpose for large-scale problems.

The survey therefore points toward methods for **approximating leverage scores efficiently**.

For the sampling distributions summarized by the paper, the approximate costs of constructing the distributions are

| Sampling method | Cost of distribution |
|---|---:|
| Uniform | $`O(1)`$ |
| Squared norms | $`O(mn)`$ |
| Leverage scores | $`O(mn^2)`$ |

This illustrates a direct trade-off between the amount of information used to construct the sampling distribution and its computational cost.

---

**Determinantal Point Process (DPP) sampling**

The paper also discusses determinantal point processes, where subsets of rows are selected jointly rather than treating each draw independently.

For a subset $S$,

$$p_S = \frac{\det\left((AA^*)(S,S)\right)}{\det(I+AA^*)}.$$

The determinant encourages selection of rows that span diverse directions rather than repeatedly selecting highly similar vectors.

DPP sampling can have strong theoretical properties, but the corresponding probability distribution is considerably more expensive to construct and sample from.

---

### Randomized Low-Rank Matrix Decompositions

The randomized dimensionality-reduction tools above can be inserted into classical algorithms such as the SVD, QR, ID, and CUR decompositions.

The general strategy is

$$\boxed{\text{large matrix} \rightarrow \text{randomly compressed representation} \rightarrow \text{classical deterministic algorithm on the smaller problem}}.$$

For example, randomized SVD first constructs a low-dimensional approximation to the dominant range of $A$,

$$Y = A\Omega, \qquad Q = \mathrm{orth}(Y),$$

then compresses

$$B = Q^*A$$

and computes the expensive SVD only on the much smaller matrix $B$.

I skipped the detailed randomized matrix algorithms in Section 4 on this pass because they primarily show how the sketching and sampling tools introduced above are inserted into the classical decompositions.

These methods can be revisited and implemented in a later pass.

---

### Randomized Low-Rank Tensor Decompositions

For a future read

## Results (Taxonomy)

The authors classify the current literature into the following subsections:

- Naive approaches of decomposing matrices, such as the SVD, QR, LU, ID, CUR. The main drawback of these approaches is that the computational cost is equivalent to processing the whole matrix, as these approaches are *exact* decompositions.

- Randomized methods and dimensionality reduction algorithms, such as sketching and sampling. The drawback of these approaches is that they are an *approximation*, not exact. However, the most advanced methods in the space tend to resemble the original matrix quite accurately with the advantage of being much cheaper computationally. Not only that, they are also much more scalable.

- Randomized low-rank matrix algorithms: Randomized sketches and samples can then be incorporated into classical algorithms such as SVD, QR, ID, and CUR. The general structure is

$$\text{compress/sample} \rightarrow \text{solve a smaller deterministic problem} \rightarrow \text{lift the result back to the original space}.$$

These methods trade exact processing of the full matrix for a controlled approximation, often with probabilistic error guarantees.

- Randomized Low-Rank Tensor Decompositions

A useful way to interpret the overall taxonomy is through the trade-off

$$\boxed{\text{accuracy} \quad\leftrightarrow\quad \text{computational cost and scalability}}$$

Classical methods prioritize direct and highly accurate matrix factorizations, while randomized methods introduce controlled approximation in exchange for less computational cost and higher scalability.

## Questions and gaps

Things I have not covered yet:

- Operator norm $$\|A\|_2 = \sup_{\|x\|=1}\|Ax\|_2$$
- Frobenius norm $$\|A\|_F = \left(\sum_{i,j}|A(i,j)|^2\right)^{1/2}$$
- The (Hermitian) transpose of $A$
- Moore-Penrose pseudoinverse of $A$
- Isotropic vectors
- Rademacher distribution and complex Rademacher distribution

I do not quite understand section 2.5.3 (ID/CUR from one truncated SVD, Osinsky's method), I get lost in the notation in the iteration section.

Questions:

- In randomized matrix sampling, why most of the methods sample with probabilities that are proportional to the norm of the vectors $`\|a_i\|_2^2`$. Why does this work mathematically? (Couldnt get to the derivation) Is it because vectors that have higher euclidean norms carry more information about the matrix?

- In randomized matrix sampling, is the row leverage scores method related to the effective dimension of the matrix, $`d_{\mathrm{eff}}`$, in some way?

Is this correct?

Alaoui & Mahoney (2015) claim that sampling landmarks proportional to these ridge leverage scores gives much better approximation at the same $m$, with *guarantees*. There is an important identity,

$$\sum_i \ell_i(\lambda) = \mathrm{tr}\left(K\left(K + \lambda n I\right)^{-1}\right) = d_{\mathrm{eff}}(\lambda).$$

- I have heard about the possibility of rotating a matrix in the direction that contains the most data about A. How does that work mathematically and how can you use it to design more efficient algorithms?

## Open problems / extensions of the work

In the conclusion, the authors claim that the literature in the methods matrix case is reaching a "certain state of maturity", while the tensor environment is still evolving rapidly.

I couldn't find any potential extensions of the work explicitly mentioned, but throughout the paper, the authors suggest that there are open fields/questions in the application problems, as described in the introduction.

## Relevance to my work

Some of the citations seem rather interesting for my work on scalable kernel learning and scalable adaptive sparse Gaussian Processes (Adaptive SGPs):

- [159] Y. Panagakis, J. Kossaifi, G. Chrysos, J. Oldfield, T. Patti, M. Nicolaou, A. Anandkumar, and S. Zafeiriou. Tensor methods in deep learning, pages 1009–1048. Academic Press, 01 2024.

- [5] A. Aidini, G. Tsagkatakis, N. D. Sidiropoulos, and P. Tsakalides. Few-shot classification using tensor completion. In Proceedings of the 57th Annual Asilomar Conference on Signals, Systems, and Computers, pages 1283–1287. IEEE, 2023.

- [58] M. Dereziński and M. W. Mahoney. Recent and upcoming developments in randomized numerical linear algebra for machine learning. In Proceedings of the 30th ACM SIGKDD Conference on Knowledge Discovery and Data Mining, KDD '24, page 6470–6479, New York, NY, USA, 2024. Association for Computing Machinery.

Page 16 argues the following: "The coherence is a measure of how evenly information in the matrix is distributed; high coherence indicates that its key information is very localized and that uniform sampling may perform poorly as a result; see, e.g., [104]. In spite of this drawback, uniform sampling is a popular choice for large-scale linear algebra computations as a cost-effective OSE. However, it is unclear how to best choose the embedding dimension $d$ in practice since the coherence is typically unavailable, cf. [141, Section 9.6.4]." This could be a potential direction to explore, especially applied to scalable kernel learning. There might be a way to keep kernel computations and learning cheap if we know how information is distributed in the matrix and how concentrated it is. If we know that information is highly concentrated within one specific section of the matrix, we could attack a much smaller problem that retains most of the information about the data. And the algorithm could be dynamic as the coherence measure changes. Similar idea for scalable adaptive sparse Gaussian Processes.

There's a big drawback though. The authors cite that computing leverage scores are as expensive as svd(A), which kind of defeats the point of scalability. Need to think on this.... The authors also claim there exist fast methods to approximate them, cf. [6, 47, 71, 149]."

Uniform sampling is $`O(1)`$, squared norms is $`O(mn)`$, and leverage scores is $`O(mn^2)`$. Leverage scores seems rather expensive for a simple sampling task, there must be some way of efficiently solving this problem.