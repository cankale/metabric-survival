def kaplan_meier(zeit, ereignisse):
    S_t = 1.0
    S_hat = []
    N = len(zeit)
    pairs = sorted(zip(zeit, ereignisse), key=lambda p: (p[0], -p[1]))
    for i, (t, d) in enumerate(pairs):
        n = N - i
        S_t *= (1.0 - d / n)
        S_hat.append((t,S_t))

    return S_hat