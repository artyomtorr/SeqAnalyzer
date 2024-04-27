from multiprocessing import Pool

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(
        self, n_estimators=10, max_depth=None, max_features=None, random_state=42
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, X, y, n_jobs=1): 
        self.classes_ = sorted(np.unique(y))
        if n_jobs == 1:
            for i in range(self.n_estimators):
                np.random.seed(self.random_state+i)
                feat_ids = np.random.choice(X.shape[1], size=self.max_features, replace=False)
                self.feat_ids_by_tree.append(feat_ids.tolist())
                boot_ids = np.random.choice(X.shape[0], size=X.shape[0], replace=True)
                x_boot = X[boot_ids]
                y_boot = y[boot_ids]
                tree = DecisionTreeClassifier(
                    max_depth=self.max_depth,
                    max_features=self.max_features,
                    random_state=self.random_state
                )
                tree.fit(x_boot[:, feat_ids], y_boot)
                self.trees.append(tree)
        else:
            with Pool(n_jobs) as p:
                results = p.starmap(self._fit_tree, [(X, y, i) for i in range(self.n_estimators)])
            self.trees, self.feat_ids_by_tree = zip(*results)
        return self

    def _fit_tree(self, X, y, i):
        np.random.seed(self.random_state+i)
        feat_ids = np.random.choice(X.shape[1], size=self.max_features, replace=False)
        boot_ids = np.random.choice(X.shape[0], size=X.shape[0], replace=True)
        x_boot = X[boot_ids]
        y_boot = y[boot_ids]
        tree = DecisionTreeClassifier(
            max_depth=self.max_depth,
            max_features=self.max_features,
            random_state=self.random_state
        )
        tree.fit(x_boot[:, feat_ids], y_boot)
        return tree, feat_ids.tolist()

    def predict_proba(self, X, n_jobs=1):
        if n_jobs == 1:
            y_pred_proba = np.zeros(shape=(X.shape[0], len(self.classes_)))
            for i in range(self.n_estimators):
                feat_ids = self.feat_ids_by_tree[i]
                y_pred_proba += self.trees[i].predict_proba(X[:, feat_ids])
            return y_pred_proba / self.n_estimators
        else:
            with Pool(n_jobs) as p:
                results = p.starmap(self._predict_proba_tree, [(X, i) for i in range(self.n_estimators)])
            y_pred_proba_sum = sum(list(results))
            return y_pred_proba_sum / self.n_estimators

    def _predict_proba_tree(self, X, i):
        feat_ids = self.feat_ids_by_tree[i]
        y_pred_proba = self.trees[i].predict_proba(X[:, feat_ids])
        return y_pred_proba

    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)
        return predictions