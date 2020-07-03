class WeightedRandomizer:
    def __init__ (self, weights):
        self.__max = .0
        self.__weights = []
        for value, weight in weights.items ():
            self.__max += weight
            self.__weights.append ( (self.__max, value) )

    def random (self):
        import random
        r = random.random() * self.__max
        for ceil, value in self.__weights:
            if ceil > r: 
                return value


def grid_KFold_RF(X, y, title, k=5, n_iter=100, n_estimators=100, random_state=None, verbose=True, plot=False):
    """Random forest model with k-fold validation
    
    Paramters
    ---------
    X : pandas.DataFrame,
        input features
    y : pandas.Series,
        label of rows
    k : int, default 5
        k-fold validation,
    n_iter : int, default 100
        number of iteration
    n_estimators : int, default 100
        number of estimator in Random Forest model
    random_state : default None
        random state in Random Forest model
    verbose : int, default 1
        if output in verbose mode (with progress bar and ROC-curve)

    Returns
    -------
    if verbose == 1:
        title, mean_fpr, tprs, mean_auc, std_auc
    else:
        title, mean_auc, std_auc
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import interp
    from sklearn.model_selection import KFold
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_curve, auc
    from logger import ProgressBar
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Arial"

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    rf = RandomForestClassifier(n_jobs=6,
                                n_estimators=n_estimators, 
                                random_state=random_state,
                                oob_score=True)

    if verbose:
        prog = logger.ProgressBar()
        prog.update(0)

    for i in range(n_iter):
        tmp_train = []
        tmp_test = []
        tmp_tprs = []
        tmp_aucs = []

        # get valid split data
        is_valid = 1
        cv = KFold(k, shuffle=True)
        for train, test in cv.split(X, y):
            # remove invalid split
            if not (0 < np.mean(y[train]) < 1 and 0 < np.mean(y[test]) < 1):
                is_valid = 0
            tmp_train.append(train)
            tmp_test.append(test)

        if is_valid == 0:
            continue

        # Train
        for train, test in zip(tmp_train, tmp_test):
            rf = RandomForestClassifier(n_estimators=n_estimators, random_state=random_state, oob_score=True)
            probas_ = rf.fit(X.iloc[train], y[train]).predict_proba(X.iloc[test])

            # Compute ROC curve and area the curve
            fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
            roc_auc = auc(fpr, tpr)

            tmp_tprs.append(interp(mean_fpr, fpr, tpr))        
            tmp_tprs[-1][0] = 0.0
            tmp_aucs.append(roc_auc)

        tmp_mean_tpr = np.mean(tmp_tprs, axis=0)
        tmp_mean_tpr[-1] = 1.0
        tprs.append(tmp_mean_tpr)
        aucs.append(auc(mean_fpr, tmp_mean_tpr))

        if verbose:
            prog.update(100 * i // n_iter)
        if plot:
            plt.plot(mean_fpr, tmp_mean_tpr, lw=1, alpha=0.3, )


    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)

    if verbose:
        prog.update(100)
    
    if plot:
        plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
                label='Chance', alpha=.8)
        plt.plot(mean_fpr, mean_tpr, color='b',
                label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
                lw=2, alpha=.8)
        plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                        label=r'$\pm$ 1 std. dev.')
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic example')
        plt.legend(loc="lower right")
        plt.show()
        
    return title, mean_fpr, tprs, mean_auc, std_auc


def calcNormFactors(df):
    import numpy as np
    f75 = df.apply(lambda x: np.quantile(x, .75), axis=0)
    ref_idx = abs(f75-np.mean(f75)).idxmin()

    f = df.apply(lambda x: calcFactorTMM(x, df[ref_idx]), axis=0)
    
    f = f / np.exp(np.mean(np.log(f)))
    return f


def calcFactorTMM(obs, ref, libsize_obs=None, libsize_ref=None, logratioTrim=.3, sumTrim=0.05, doWeighting=True, Acutoff=-1e10):
    import numpy as np
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