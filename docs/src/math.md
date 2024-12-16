
# Single Phase Bus Injection Model (Unrelaxed)

Notation:
- ``s_j`` net complex power injection at node ``j``
- ``\mathcal{N}` set of all nodes in network
- ``v_j`` voltage at node ``j``

```math
\begin{aligned}
s_j = \sum_{k: j \sim k} Y[j,k]^* ( |v_j|^2 - v_j v_k^*) \ \forall j \in \mathcal{N} \\
v_{\text{substation bus}} = v_0
\end{aligned}
```
