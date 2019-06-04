from GAT2VEC.evaluation.classification import Classification
from GAT2VEC.gat2vec import Gat2Vec
from guiltytargets.constants import gat2vec_config

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

