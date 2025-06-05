import numpy as np
from scipy.optimize import minimize
from itertools import combinations

def simulate_allele_frequencies(fst_matrix, num_variants=100, random_seed=None):
    """
    Simulate population allele frequencies based on a pairwise FST matrix.
    
    Parameters:
    -----------
    fst_matrix : numpy.ndarray
        Square matrix of pairwise FST values between populations
    num_variants : int
        Number of genetic variants to simulate
    random_seed : int, optional
        Seed for random number generation
        
    Returns:
    --------
    numpy.ndarray
        Matrix of allele frequencies with shape (num_variants, num_populations)
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    # Number of populations
    num_pops = fst_matrix.shape[0]
    
    # Generate allele frequencies
    allele_freqs = np.zeros((num_variants, num_pops))
    
    # For each variant
    for i in range(num_variants):
        # Start with a random ancestral allele frequency
        p_ancestral = np.random.beta(2, 2)  # Beta distribution centered around 0.5
        
        # Initial allele frequencies - start all populations with the same frequency
        # but with small random perturbations
        p = np.array([p_ancestral + np.random.normal(0, 0.01) for _ in range(num_pops)])
        p = np.clip(p, 0.05, 0.95)  # Ensure frequencies are not too extreme
        
        # Apply differentiation between populations based on FST
        # We'll iteratively adjust population pairs to match their target FST
        for _ in range(3):  # Multiple passes to improve convergence
            for j in range(num_pops):
                for k in range(j+1, num_pops):
                    target_fst = fst_matrix[j, k]
                    
                    if target_fst <= 0:
                        continue  # Skip if target FST is zero or negative
                    
                    # Current allele frequencies
                    p_j = p[j]
                    p_k = p[k]
                    
                    # Calculate current FST between these two populations
                    p_mean = (p_j + p_k) / 2
                    h_t = 2 * p_mean * (1 - p_mean)
                    h_s = (2 * p_j * (1 - p_j) + 2 * p_k * (1 - p_k)) / 2
                    
                    if h_t <= 0:
                        continue  # Skip if h_t is zero or negative
                    
                    current_fst = max(0, (h_t - h_s) / h_t)
                    
                    # If current FST is too low, increase differentiation
                    if current_fst < target_fst:
                        # Determine direction of adjustment (which population gets higher/lower)
                        if p_j >= p_k:
                            # Increase p_j, decrease p_k
                            adjust = min(0.05, (target_fst - current_fst))
                            p[j] = min(0.95, p_j + adjust * (1 - p_j))
                            p[k] = max(0.05, p_k - adjust * p_k)
                        else:
                            # Decrease p_j, increase p_k
                            adjust = min(0.05, (target_fst - current_fst))
                            p[j] = max(0.05, p_j - adjust * p_j)
                            p[k] = min(0.95, p_k + adjust * (1 - p_k))
                    
                    # If current FST is too high, decrease differentiation
                    elif current_fst > target_fst * 1.2:  # Allow some tolerance
                        # Move frequencies closer to their mean
                        adjust = min(0.05, (current_fst - target_fst))
                        p[j] = p_j + adjust * (p_mean - p_j)
                        p[k] = p_k + adjust * (p_mean - p_k)
        
        # Store the frequencies
        allele_freqs[i] = p
    
    # Fine-tune allele frequencies to better match target FST
    tuned_freqs = fine_tune_frequencies(allele_freqs, fst_matrix)
    
    return tuned_freqs

def calculate_fst(allele_freqs):
    """
    Calculate pairwise FST values from population allele frequencies.
    
    Parameters:
    -----------
    allele_freqs : numpy.ndarray
        Matrix of allele frequencies with shape (num_variants, num_populations)
        
    Returns:
    --------
    numpy.ndarray
        Matrix of pairwise FST values
    """
    num_variants, num_pops = allele_freqs.shape
    fst_matrix = np.zeros((num_pops, num_pops))
    
    # Calculate pairwise FST for different populations
    for i, j in combinations(range(num_pops), 2):
        p_i = allele_freqs[:, i]
        p_j = allele_freqs[:, j]
        
        # Calculate FST for each variant
        fst_values = []
        
        for v in range(num_variants):
            # Mean allele frequency between the two populations
            p_mean = (p_i[v] + p_j[v]) / 2
            
            # Total heterozygosity
            h_t = 2 * p_mean * (1 - p_mean)
            
            # Average within-population heterozygosity
            h_s = (2 * p_i[v] * (1 - p_i[v]) + 2 * p_j[v] * (1 - p_j[v])) / 2
            
            # Calculate FST only when h_t > 0 to avoid division by zero or negative values
            if h_t > 0:
                fst = (h_t - h_s) / h_t
                # Ensure FST is between 0 and 1
                fst = max(0, min(1, fst))
                fst_values.append(fst)
        
        # Average FST across all variants
        if fst_values:
            fst_matrix[i, j] = fst_matrix[j, i] = np.mean(fst_values)
    
    return fst_matrix

def fine_tune_frequencies(allele_freqs, target_fst_matrix, max_iter=5000):
    """
    Fine-tune allele frequencies to better match target FST values.
    
    Parameters:
    -----------
    allele_freqs : numpy.ndarray
        Initial matrix of allele frequencies
    target_fst_matrix : numpy.ndarray
        Target matrix of pairwise FST values
    max_iter : int
        Maximum number of iterations for tuning
        
    Returns:
    --------
    numpy.ndarray
        Tuned matrix of allele frequencies
    """
    num_variants, num_pops = allele_freqs.shape
    tuned_freqs = allele_freqs.copy()
    
    for iteration in range(max_iter):
        # Calculate current FST matrix
        current_fst = calculate_fst(tuned_freqs)
        
        # Calculate mean absolute error
        all_errors = []
        for i, j in combinations(range(num_pops), 2):
            all_errors.append(abs(target_fst_matrix[i, j] - current_fst[i, j]))
        
        mean_error = np.mean(all_errors) if all_errors else 0
        print(f"Iteration {iteration+1}, Mean absolute error: {mean_error:.4f}")
        
        if mean_error < 0.001:
            print("Target accuracy reached, stopping iterations")
            break
        
        # For each population pair
        for i, j in combinations(range(num_pops), 2):
            target = target_fst_matrix[i, j]
            current = current_fst[i, j]
            
            if abs(target - current) < 0.01:
                continue  # Skip if already close enough
            
            # Randomly select a subset of variants to adjust
            num_to_adjust = max(5, int(num_variants * 0.2))
            adjust_variants = np.random.choice(num_variants, size=num_to_adjust, replace=False)
            
            for var_idx in adjust_variants:
                p_i = tuned_freqs[var_idx, i]
                p_j = tuned_freqs[var_idx, j]
                p_mean = (p_i + p_j) / 2
                
                # Need to increase FST (increase differentiation)
                if current < target:
                    # Determine which frequency to increase and which to decrease
                    if p_i >= p_j:
                        # Move them further apart
                        adjust = min(0.05, (target - current) / 2)
                        tuned_freqs[var_idx, i] = min(0.95, p_i + adjust * (1 - p_i))
                        tuned_freqs[var_idx, j] = max(0.05, p_j - adjust * p_j)
                    else:
                        # Move them further apart in the opposite direction
                        adjust = min(0.05, (target - current) / 2)
                        tuned_freqs[var_idx, i] = max(0.05, p_i - adjust * p_i)
                        tuned_freqs[var_idx, j] = min(0.95, p_j + adjust * (1 - p_j))
                
                # Need to decrease FST (decrease differentiation)
                else:
                    # Move frequencies closer to their mean
                    adjust = min(0.05, (current - target) / 2)
                    tuned_freqs[var_idx, i] = p_i + adjust * (p_mean - p_i)
                    tuned_freqs[var_idx, j] = p_j + adjust * (p_mean - p_j)
                
                # Ensure frequencies stay within bounds
                tuned_freqs[var_idx, i] = max(0.05, min(0.95, tuned_freqs[var_idx, i]))
                tuned_freqs[var_idx, j] = max(0.05, min(0.95, tuned_freqs[var_idx, j]))
    
    return tuned_freqs

def run_simulation_example():
    """Run an example simulation and report the results"""
    # Define a simple FST matrix
    fst_matrix = np.array([
        [0.0, 0.05, 0.15],
        [0.05, 0.0, 0.10],
        [0.15, 0.10, 0.0]
    ])
    
    print("Target FST matrix:")
    print(fst_matrix)
    
    # Simulate allele frequencies
    print("\nSimulating allele frequencies...")
    allele_freqs = simulate_allele_frequencies(fst_matrix, num_variants=100, random_seed=42)
    
    # Calculate the resulting FST matrix
    result_fst = calculate_fst(allele_freqs)
    print("\nCalculated FST from simulated frequencies:")
    print(result_fst)
    
    # Calculate errors
    error = np.abs(fst_matrix - result_fst)
    print("\nAbsolute error by element:")
    print(error)
    
    # Calculate mean absolute error for off-diagonal elements
    off_diag_indices = np.where(~np.eye(fst_matrix.shape[0], dtype=bool))
    mean_abs_error = np.mean(error[off_diag_indices])
    print(f"\nMean absolute error: {mean_abs_error:.4f}")
    
    np.savetxt("allele_freqs.txt", allele_freqs)
    return allele_freqs, result_fst

if __name__ == "__main__":
    run_simulation_example()
