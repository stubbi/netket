import nqs
import networkx as nx
import numpy as np

operators = {}

# Ising 1D
g = nqs.graph.Hypercube(length=20, n_dim=1, pbc=True)
hi = nqs.hilbert.Spin(s=0.5, graph=g)
operators["Ising 1D"] = nqs.operator.Ising(h=1.321, hilbert=hi)

# Heisenberg 1D
g = nqs.graph.Hypercube(length=20, n_dim=1, pbc=True)
hi = nqs.hilbert.Spin(s=0.5, total_sz=0, graph=g)
operators["Heisenberg 1D"] = nqs.operator.Heisenberg(hilbert=hi)

# Bose Hubbard
g = nqs.graph.Hypercube(length=3, n_dim=2, pbc=True)
hi = nqs.hilbert.Boson(n_max=3, n_bosons=6, graph=g)
operators["Bose Hubbard"] = nqs.operator.BoseHubbard(U=4.0, hilbert=hi)

# Graph Hamiltonian
sigmax = [[0, 1], [1, 0]]
mszsz = [[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1]]
edges = [
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 8],
    [8, 9],
    [9, 10],
    [10, 11],
    [11, 12],
    [12, 13],
    [13, 14],
    [14, 15],
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 0],
]

g = nqs.graph.CustomGraph(edges=edges)
hi = nqs.hilbert.CustomHilbert(local_states=[-1, 1], graph=g)
ha = nqs.operator.GraphOperator(
    hi, siteops=[sigmax], bondops=[mszsz], bondops_colors=[0]
)
operators["Graph Hamiltonian"] = ha

# Custom Hamiltonian
sx = [[0, 1], [1, 0]]
sy = [[0, 1.0j], [-1.0j, 0]]
sz = [[1, 0], [0, -1]]
g = nqs.graph.CustomGraph(edges=[[i, i + 1] for i in range(20)])
hi = nqs.hilbert.CustomHilbert(local_states=[1, -1], graph=g)

sx_hat = nqs.operator.LocalOperator(hi, [sx] * 3, [[0], [1], [5]])
sy_hat = nqs.operator.LocalOperator(hi, [sy] * 4, [[2], [3], [4], [9]])
szsz_hat = nqs.operator.LocalOperator(hi, sz, [0]) * nqs.operator.LocalOperator(
    hi, sz, [1]
)
szsz_hat += nqs.operator.LocalOperator(hi, sz, [4]) * nqs.operator.LocalOperator(
    hi, sz, [5]
)
szsz_hat += nqs.operator.LocalOperator(hi, sz, [6]) * nqs.operator.LocalOperator(
    hi, sz, [8]
)
szsz_hat += nqs.operator.LocalOperator(hi, sz, [7]) * nqs.operator.LocalOperator(
    hi, sz, [0]
)

operators["Custom Hamiltonian"] = sx_hat + sy_hat + szsz_hat
operators["Custom Hamiltonian Prod"] = sx_hat * 1.5 + (2.0 * sy_hat)

rg = nqs.utils.RandomEngine(seed=1234)


def test_produce_elements_in_hilbert():
    for name, ha in operators.items():
        hi = ha.hilbert
        print(name, hi)
        assert len(hi.local_states) == hi.local_size
        assert hi.size > 0
        rstate = np.zeros(hi.size)

        local_states = hi.local_states

        for i in range(1000):
            hi.random_vals(rstate, rg)

            conns = ha.get_conn(rstate)

            for connector, newconf in zip(conns[1], conns[2]):
                rstatet = np.array(rstate)
                hi.update_conf(rstatet, connector, newconf)

                for rs in rstatet:
                    assert rs in local_states


def test_operator_is_hermitean():
    for name, ha in operators.items():
        hi = ha.hilbert
        print(name, hi)
        assert len(hi.local_states) == hi.local_size

        rstate = np.zeros(hi.size)

        local_states = hi.local_states

        for i in range(100):
            hi.random_vals(rstate, rg)
            conns = ha.get_conn(rstate)

            for mel, connector, newconf in zip(conns[0], conns[1], conns[2]):
                rstatet = np.array(rstate)
                hi.update_conf(rstatet, connector, newconf)

                conns1 = ha.get_conn(rstatet)
                foundinv = False
                for meli, connectori, newconfi in zip(conns1[0], conns1[1], conns1[2]):
                    rstatei = np.array(rstatet)
                    hi.update_conf(rstatei, connectori, newconfi)
                    if np.array_equal(rstatei, rstate):
                        foundinv = True
                        assert meli == np.conj(mel)
                assert foundinv

def test_no_segfault():
    g = nqs.graph.Hypercube(8, 1)
    hi = nqs.hilbert.Spin(g, 0.5)

    lo = nqs.operator.LocalOperator(hi, [[1,0],[0,1]], [0])
    lo = lo.transpose()

    hi = None

    lo = lo * lo

    assert True
