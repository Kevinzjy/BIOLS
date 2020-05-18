import random
import numpy as np
import utils
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['font.sans-serif'] = "Arial"


class WeightedRandomizer:
    def __init__ (self, weights):
        self.__max = .0
        self.__weights = []
        for value, weight in weights.items ():
            self.__max += weight
            self.__weights.append ( (self.__max, value) )

    def random (self):
        r = random.random() * self.__max
        for ceil, value in self.__weights:
            if ceil > r: 
                return value


def get_n50(sequence_lengths):
    sequence_lengths = sorted(sequence_lengths, reverse=True)
    total_bases = sum(sequence_lengths)
    target_bases = total_bases * 0.5
    bases_so_far = 0
    for sequence_length in sequence_lengths:
        bases_so_far += sequence_length
        if bases_so_far >= target_bases:
            return sequence_length
    return 0


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
    from scipy import interp
    from sklearn.model_selection import KFold
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_curve, auc

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    rf = RandomForestClassifier(n_jobs=6,
                                n_estimators=n_estimators, 
                                random_state=random_state,
                                oob_score=True)

    if verbose:
        prog = utils.ProgressBar()
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
