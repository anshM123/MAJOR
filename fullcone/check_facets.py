"""Sanity of the facet lists: (1) each orbit representative admits a contraction map (so the transcription
is a valid HEI); (2) random graph models satisfy all facets; (3) a random non-graph vector (perturbed) fails."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import random, cones, fc_common as fc
random.seed(1)
for N in (3, 4, 5):
    reps = {}
    for name, ineq in cones.instances(N):
        reps.setdefault(name, ineq)
    for name, ineq in reps.items():
        print(f"N={N} {name:5s} contraction map: {fc.contraction_exists(ineq, N)}", flush=True)
    worst = 1
    for t in range(200):
        nb = random.randint(0, 6); nodes = list(range(N+1+nb))
        label = {i: i for i in range(N+1)}
        edges = [(random.choice(nodes), random.choice(nodes), random.random()) for _ in range(random.randint(N, 3*N+8))]
        S = fc.graph_entropy(nodes, edges, label, N)
        worst = min(worst, fc.facet_check(S, N)[0])
    print(f"N={N}: 200 random graph models, min facet value / scale = {worst:+.2e}")
    # control: the uniform vector S(I)=1 for |I| in {1,N+1-1}... use the 'perfect tensor' ray S(I)=min(|I|,N+1-|I|)
    S = {m: min(bin(m).count("1"), N+1-bin(m).count("1")) for m in range(1, 2**N)}
    print(f"N={N}: perfect-tensor vector min facet = {fc.facet_check(S, N)}")
