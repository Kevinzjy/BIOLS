def calcNormFactors(df):
    f75 = df.apply(lambda x: np.quantile(x, .75), axis=0)
    ref_idx = abs(f75-np.mean(f75)).idxmin()

    f = df.apply(lambda x: calcFactorTMM(x, df[ref_idx]), axis=0)
    
    f = f / np.exp(np.mean(np.log(f)))
    return f


def calcFactorTMM(obs, ref, libsize_obs=None, libsize_ref=None, logratioTrim=.3, sumTrim=0.05, doWeighting=True, Acutoff=-1e10):
    # TMM between two libraries
    obs = obs.astype(float)
    ref = ref.astype(float)

    nO = obs.sum() if libsize_obs is None else libsize_obs
    nR = obs.sum() if libsize_ref is None else libsize_ref

    logR = np.log2((obs/nO)/(ref/nR))               # log ratio of expression, accounting for library size
    absE = (np.log2(obs/nO) + np.log2(ref/nR)) / 2  # absolute expression
    v = (nO-obs)/nO/obs + (nR-ref)/nR/ref           # estimated asymptotic variance

    # remove infinite values, cutoff based on A
    fin = np.isfinite(logR) & np.isfinite(absE) & (absE > Acutoff)

    logR = logR[fin]
    absE = absE[fin]
    v = v[fin]

    # if max(abs(logR)) < 1e-6:
    #     return 1

    # taken from the original mean() function
    n = logR.shape[0]
    loL = np.floor(n * logratioTrim) + 1
    hiL = n + 1 - loL
    loS = np.floor(n * sumTrim) + 1
    hiS = n + 1 - loS

    # keep <- (rank(logR) %in% loL:hiL) & (rank(absE) %in% loS:hiS)
    # a fix from leonardo ivan almonacid cardenas, since rank() can return
    # non-integer values when there are a lot of ties
    keep = set(logR.sort_values().index[int(loL):int(hiL)]) & set(absE.sort_values().index[int(loS):int(hiS)])

    if doWeighting:
        f = sum(logR[keep]/v[keep]) / sum(1/v[keep])
    else:
        f = np.mean(logR[keep])

    if np.isnan(f):
        f = 0

    return 2 ** f