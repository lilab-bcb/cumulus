import sys
import h5py
import pegasusio as io

from scipy.sparse import csr_matrix


def convert_to_10x_h5(mtx_folder, out_name, genome):
    data = io.read_input(mtx_folder, genome=genome)
    n_obs, n_feature = data.shape

    with h5py.File(f"{out_name}.h5", "w") as f:
        grp = f.create_group("matrix")

        grp.create_dataset(name="barcodes", data=data.obs_names.map(lambda s: f"{s}-1".encode("ascii", "ignore")).values)

        if not isinstance(data.X, csr_matrix):
            data.X = csr_matrix(data.X)
        grp.create_dataset(name="data", data=data.X.data)
        grp.create_dataset(name="indices", data=data.X.indices)
        grp.create_dataset(name="indptr", data=data.X.indptr)
        grp.create_dataset(name="shape", data=(n_feature, n_obs)) # feature-by-barcode

        feature_grp = grp.create_group("features")
        feature_grp.create_dataset(name="_all_tag_keys", data=[b"genome"])
        feature_grp.create_dataset(name="feature_type", data=data.var['featuretype'].apply(lambda s: s.encode("ascii", "ignore")).values)
        feature_grp.create_dataset(name="genome", shape=(n_feature,), dtype=f"S{len(genome)}", fillvalue=genome.encode("ascii", "ignore"))
        feature_grp.create_dataset(name="id", data=data.var['featureid'].apply(lambda s: s.encode("ascii", "ignore")).values)
        feature_grp.create_dataset(name="name", data=data.var_names.map(lambda s: s.encode("ascii", "ignore")).values)


if __name__ == "__main__":
    mtx_folder = sys.argv[1]
    out_name = sys.argv[2]
    genome = "Unknown" if len(sys.argv) <= 3 else sys.argv[3]

    convert_to_10x_h5(mtx_folder, out_name, genome)
