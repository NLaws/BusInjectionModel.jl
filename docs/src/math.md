
# Single Phase Bus Injection Model (Unrelaxed)

## Rectangular voltage variables
Notation:
- ``s_j`` net complex power injection at node ``j``
- ``\mathcal{N}`` set of all nodes in network
- ``v_j`` voltage at node ``j``

```math
\begin{aligned}
&s_j = \sum_{k: j \sim k} Y[j,k]^* ( |v_j|^2 - v_j v_k^*) \ \forall j \in \mathcal{N} \\
&v_{\text{substation bus}} = v_0
\end{aligned}
```

## Polar voltage variables
Notation:
- ``G_{ij}`` real entry of admittance matrix at row ``i``, column ``j``
- ``B_{ij}`` imaginary entry of admittance matrix at row ``i``, column ``j``
- ``|v_j|`` voltage magnitude at bus ``j``
- ``\angle v_j`` voltage angle at bus ``j``
- ``p_j`` real power injection at bus ``j``
- ``q_j`` reactive power injection at bus ``j``


```math
\begin{aligned}
&p_j =  |v_j| \sum_{i \in 1\dots\mathcal{N}} |v_i| \left[
     G_{ij} \cos(\angle v_j - \angle v_i) + B_{ij} \sin(\angle v_j - \angle v_i)
     \right]  \\
&q_j =  |v_j| \sum_{i \in 1\dots\mathcal{N}} |v_i| \left[
     G_{ij} \sin(\angle v_j - \angle v_i) - B_{ij} \cos(\angle v_j - \angle v_i)
     \right]  \\
&|v_{\text{substation bus}}| = v_0 \\
&\angle v_{\text{substation bus}} = 0
\end{aligned}
```

# Multi-Phase Bus Injection Model (Unrelaxed)

## Rectangular voltage variables
Notation:
- ``\boldsymbol s_j`` net complex power injection at bus ``j``, vector of phases
- ``\Phi_{j}`` phases connected to bus ``j`` (take sub-set of vector)
- ``\boldsymbol v_j^{\Phi_{jk}}`` complex voltage vector at bus ``j`` for the phases connected to
  bus ``k``
- ``k: j \sim k`` set of busses connected to bus ``j``
- ``H`` conjugate transpose
- ``\boldsymbol Y_{jk}^H`` phase admittance matrix for busses ``j`` and ``k``


```math
\begin{aligned}
\boldsymbol s_j^{\Phi_{j}} = \sum_{k: j \sim k} \text{diag} \left(  
     \boldsymbol v_j^{\Phi_{jk}} \left[ \boldsymbol v_j^{\Phi_{jk}} - \boldsymbol v_k^{\Phi_{jk}}  \right]^H \boldsymbol Y_{jk}^H
\right)
\quad \forall j \in \mathcal{N}
\end{aligned}
```