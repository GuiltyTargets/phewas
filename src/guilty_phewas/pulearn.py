# -*- coding: utf-8 -*-

from collections import defaultdict

from GAT2VEC.evaluation.classification import Classification
from GAT2VEC import parsers, paths
import pandas as pd
from sklearn import svm
from sklearn.metrics import f1_score, accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

# TODO move to GAT2VEC (forked)?


class PULearn(Classification):
    """ """

    def evaluate(self, model, label=False, evaluation_scheme="tr", cost_p: float = None, cost_n: float = None):
        """Evaluates the model according to the given evaluation scheme.

        :param model:
        :param label:
        :param evaluation_scheme:
        :param cost_p: The cost given to the positive values # TODO cost or weight??
        :param cost_n: The cost given to the negative values # TODO cost or weight??
        :return:
        """
        embedding = 0
        clf = self.get_classifier()

        if not label:
            embedding = parsers.get_embeddingDF(model)

        if evaluation_scheme == "bsvm":
            if cost_n and cost_p and cost_n > 0 and cost_p > 0:
                results = self.evaluate_bsvm(embedding, 5, cost_p=cost_p, cost_n=cost_n)
            else:
                raise ValueError("Biased SVM requires both negative and positive weights.")
        if evaluation_scheme == "cv":
            results = self.evaluate_cv(clf, embedding, 5)
        elif evaluation_scheme == "tr" or label:
            results = defaultdict(list)
            for tr in self.TR:
                print("TR ... ", tr)
                if label:
                    model = paths.get_embedding_path_wl(self.dataset_dir, self.output_dir, tr)
                    if isinstance(model, str):
                        embedding = parsers.get_embeddingDF(model)
                results.update(self.evaluate_tr(clf, embedding,  tr))

        print("Training Finished")

        df = pd.DataFrame(results)
        return df.groupby(axis=0, by="TR").mean()

    def evaluate_bsvm(self, embedding, n_splits, cost_p=.5, cost_n=200.):
        """Implement a biased SVM classifier.

        :param clf: Classifier object
        :param embedding:
        :param n_splits:
        :param cost_p:
        :param cost_n:
        :return: Dictionary containing numerical results of the classification.
        """

        embedding = embedding[self.label_ind, :]
        results = defaultdict(list)
        clf = svm.SVC(probability=True)
        for i in range(10):
            rskf = StratifiedKFold(n_splits=n_splits, shuffle=True)
            for train_idx, test_idx in rskf.split(embedding, self.labels):
                X_train, X_test, Y_train, Y_test = self._get_split(embedding, test_idx, train_idx)
                pred, probs = self.get_biased_predictions(clf, X_train, X_test, Y_train, cost_p, cost_n)
                # print(f"TR {i} { {0: len([1 for x in pred if x == 0]), 1: len([1 for x in pred if x == 1])}}")
                results["TR"].append(i)
                results["accuracy"].append(accuracy_score(Y_test, pred))
                results["f1micro"].append(f1_score(Y_test, pred, average='micro'))
                results["f1macro"].append(f1_score(Y_test, pred, average='macro'))
                if self.label_count == 2:
                    results["auc"].append(roc_auc_score(Y_test, probs[:, 1]))
                else:
                    results["auc"].append(0)
                print(f"TR {i}")
        return results

    @staticmethod
    def get_biased_predictions(clf, X_train, X_test, Y_train, cost_p, cost_n):
        sample_weight = [cost_p if x == 1 else cost_n for x in Y_train]
        clf.fit(X_train, Y_train, sample_weight)
        return clf.predict(X_test), clf.predict_proba(X_test)


if __name__ == '__main__':
    from GAT2VEC.gat2vec import Gat2Vec

    dir_ = "C:/Users/Mauricio/Thesis/bel_data/alzh"
    g2v = Gat2Vec(dir_, dir_, label=False, tr=[0.1, 0.3, 0.5])
    walk_length = 4
    num_walks = 30
    window_size = 10
    dimension = 128

    model = g2v.train_gat2vec(
        num_walks,
        walk_length,
        dimension,
        window_size,
        output=True,
    )
    print('***** model obtained *****.')
    pul = PULearn(dir_, dir_, tr=[0.1, 0.3, 0.5])
    print('PU created.')
    auc_df = pul.evaluate(model, label=False, evaluation_scheme="bsvm")
    print(auc_df)
