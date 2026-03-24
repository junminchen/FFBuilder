import torch
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List, Tuple

def generate_fingerprints(smiles_list: List[str], n_bits: int = 2048) -> torch.Tensor:
    """
    Generates Morgan fingerprints for a list of SMILES.
    """
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fps.append(np.zeros(n_bits))
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
        arr = np.zeros((0,), dtype=np.int8)
        Chem.DataStructs.ConvertToNumpyArray(fp, arr)
        fps.append(arr)
    return torch.tensor(np.array(fps), dtype=torch.float32)

def simple_pca(X: torch.Tensor, n_components: int = 50) -> torch.Tensor:
    """
    Simple PCA implementation using SVD in torch.
    """
    X_mean = torch.mean(X, dim=0)
    X_centered = X - X_mean
    U, S, V = torch.pca_lowrank(X_centered, q=n_components)
    return torch.matmul(X_centered, V[:, :n_components])

class SimpleGP:
    """
    A simple Gaussian Process implementation in PyTorch.
    Uses RBF kernel.
    """
    def __init__(self, length_scale=1.0, noise=1e-4):
        self.length_scale = length_scale
        self.noise = noise
        self.X_train = None
        self.y_train = None
        self.K_inv = None

    def _rbf_kernel(self, X1, X2):
        dist_sq = torch.cdist(X1, X2).pow(2)
        return torch.exp(-0.5 * dist_sq / (self.length_scale**2))

    def fit(self, X, y):
        self.X_train = X
        self.y_train = y
        K = self._rbf_kernel(X, X) + self.noise * torch.eye(X.shape[0])
        self.K_inv = torch.inverse(K)

    def predict(self, X_test) -> Tuple[torch.Tensor, torch.Tensor]:
        K_s = self._rbf_kernel(self.X_train, X_test)
        K_ss = self._rbf_kernel(X_test, X_test) + self.noise * torch.eye(X_test.shape[0])
        
        mu = torch.matmul(K_s.t(), torch.matmul(self.K_inv, self.y_train))
        cov = K_ss - torch.matmul(K_s.t(), torch.matmul(self.K_inv, K_s))
        return mu, torch.diag(cov).clamp(min=1e-9).sqrt()

def expected_improvement(mu, sigma, best_f, maximize=True):
    """
    Calculates Expected Improvement.
    """
    from torch.distributions import Normal
    dist = Normal(torch.tensor([0.0]), torch.tensor([1.0]))
    
    if maximize:
        improvement = mu - best_f
    else:
        improvement = best_f - mu
        
    Z = improvement / sigma
    ei = improvement * dist.cdf(Z) + sigma * torch.exp(dist.log_prob(Z))
    return ei

def suggest_next(candidate_smiles, train_smiles, train_properties, n_suggest=5, maximize=True):
    """
    Mimics HiTPoly's screening loop to suggest next molecules to simulate.
    """
    if not train_smiles:
        return candidate_smiles[:n_suggest]
    
    all_smiles = train_smiles + candidate_smiles
    features = generate_fingerprints(all_smiles)
    # Reduce dimensionality if many features
    if features.shape[1] > 50:
        features = simple_pca(features, n_components=min(50, len(all_smiles)))
    
    X_train = features[:len(train_smiles)]
    X_cand = features[len(train_smiles):]
    y_train = torch.tensor(train_properties, dtype=torch.float32).view(-1, 1)
    
    gp = SimpleGP()
    gp.fit(X_train, y_train)
    
    mu, sigma = gp.predict(X_cand)
    best_f = y_train.max() if maximize else y_train.min()
    
    ei = expected_improvement(mu, sigma, best_f, maximize=maximize)
    
    # Get top indices
    _, top_indices = torch.topk(ei.view(-1), min(n_suggest, len(candidate_smiles)))
    
    return [candidate_smiles[i] for i in top_indices.tolist()]

if __name__ == "__main__":
    # Test suggest_next
    candidates = ["CCCC", "CCCCCO", "c1ccccc1", "C1CCCCC1", "CC(=O)O"]
    train_smi = ["CCO", "CCC"]
    train_prop = [0.5, 0.6]
    suggestions = suggest_next(candidates, train_smi, train_prop, n_suggest=2)
    print(f"Suggestions: {suggestions}")
