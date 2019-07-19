# -*- coding: utf-8 -*-

from guiltytargets.constants import gat2vec_config
from guiltytargets.gat2vec import Classification, Gat2Vec


def main():
    dir_ = "C:/Users/Mauricio/Thesis/bel_data/alzh"  # windows
    # dir_ = "C:\\Users\\Mauricio\\Thesis\\git\\reproduction\\data\\lc"

    classifier = Classification(dir_, dir_, tr=gat2vec_config.training_ratio)

    g2v = Gat2Vec(dir_, dir_, label=False, tr=gat2vec_config.training_ratio)
    model = g2v.train_gat2vec(
        gat2vec_config.num_walks,
        gat2vec_config.walk_length,
        gat2vec_config.dimension,
        gat2vec_config.window_size,
        output=True,
    )

    classifier = Classification(dir_, dir_, tr=gat2vec_config.training_ratio)

    auc_df = classifier.evaluate(model, label=False, evaluation_scheme="cv")


if __name__ == '__main__':
    main()
