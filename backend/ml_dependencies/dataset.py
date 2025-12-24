# ml_dependencies/dataset.py
import os
import numpy as np
import torch
from typing import Optional

# DGL: compat v2+ and legacy
try:
    from dgl import graph as _dgl_graph
    def _graph(edges, num_nodes: int):
        src, dst = edges
        if not torch.is_tensor(src):
            src = torch.as_tensor(src, dtype=torch.int64)
        else:
            src = src.to(torch.int64)
        if not torch.is_tensor(dst):
            dst = torch.as_tensor(dst, dtype=torch.int64)
        else:
            dst = dst.to(torch.int64)
        return _dgl_graph((src, dst), num_nodes=int(num_nodes))
except Exception:
    from dgl.convert import graph as _dgl_graph
    def _graph(edges, num_nodes: int):
        src, dst = edges
        src = torch.as_tensor(src, dtype=torch.int64)
        dst = torch.as_tensor(dst, dtype=torch.int64)
        return _dgl_graph((src, dst), num_nodes=int(num_nodes))


def _resolve_dataset_path(graph_representation: str, target: str, data_dir: Optional[str] = None) -> str:
    """
        Return the absolute path of nmrshiftdb2_graph_{graph_representation}_{target}.npz
        using, in order of priority:
            1) provided data_dir (if not None)
            2) NMR_DATA_DIR env var
            3) ./data next to ml_dependencies/
    """
    fname = "nmrshiftdb2_graph_{rep}_{tgt}.npz".format(rep=graph_representation, tgt=target)
        # 1) explicit parameter
    if data_dir:
        p = os.path.join(data_dir, fname)
        if os.path.exists(p):
            return os.path.abspath(p)
        # 2) env var
    env_dir = os.getenv("NMR_DATA_DIR", "").strip()
    if env_dir:
        p = os.path.join(env_dir, fname)
        if os.path.exists(p):
            return os.path.abspath(p)
        # 3) neighboring ./data
    base = os.path.dirname(__file__)  # ml_dependencies/
    p = os.path.abspath(os.path.join(base, "data", fname))
    return p


class GraphDataset(object):
    def __init__(self, target: str, graph_representation: str, data_dir: Optional[str] = None):
        self.target = target
        self.graph_representation = graph_representation
        self.split = None
        self._data_path = _resolve_dataset_path(graph_representation, target, data_dir)
        self.load()

    def load(self):
        if not os.path.exists(self._data_path):
            raise FileNotFoundError(
                "Dataset not found: {p}\n"
                "-> Place the file in the correct location or set NMR_DATA_DIR to the folder containing:\n"
                "   - nmrshiftdb2_graph_{rep}_1H.npz\n"
                "   - nmrshiftdb2_graph_{rep}_13C.npz".format(
                    p=self._data_path, rep=self.graph_representation
                )
            )

        [mol_dict] = np.load(self._data_path, allow_pickle=True)["data"]

        self.n_node = mol_dict["n_node"]
        self.n_edge = mol_dict["n_edge"]
        self.node_attr = mol_dict["node_attr"]
        self.edge_attr = mol_dict["edge_attr"]
        self.src = mol_dict["src"]
        self.dst = mol_dict["dst"]

        self.shift = mol_dict["shift"]
        self.mask = mol_dict["mask"]
        self.smi = mol_dict["smi"]

        self.n_csum = np.concatenate([[0], np.cumsum(self.n_node)])
        self.e_csum = np.concatenate([[0], np.cumsum(self.n_edge)])

    def __getitem__(self, idx):
        g = _graph(
            (self.src[self.e_csum[idx]: self.e_csum[idx + 1]],
             self.dst[self.e_csum[idx]: self.e_csum[idx + 1]]),
            num_nodes=int(self.n_node[idx]),
        )
        g.ndata["node_attr"] = torch.from_numpy(
            self.node_attr[self.n_csum[idx]: self.n_csum[idx + 1]]
        ).float()
        g.edata["edge_attr"] = torch.from_numpy(
            self.edge_attr[self.e_csum[idx]: self.e_csum[idx + 1]]
        ).float()

        n_node = self.n_node[idx: idx + 1].astype(int)
        numHshifts = np.zeros(n_node)
        shift = self.shift[self.n_csum[idx]: self.n_csum[idx + 1]]
        shift_test = shift
        mask = self.mask[self.n_csum[idx]: self.n_csum[idx + 1]].astype(bool)

        if self.target == "1H":
            # same logic as the original repo (group mean, etc.)
            shift = np.hstack([np.mean(s) for s in self.shift[idx]])
            numHshifts = np.hstack([len(s) for s in self.shift[idx][mask]])
            shift_test = np.hstack([np.hstack(s) for s in self.shift[idx][mask]])

        return g, n_node, numHshifts, shift_test, shift, mask

    def __len__(self):
        return int(self.n_node.shape[0])
