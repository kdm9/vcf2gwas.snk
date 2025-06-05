import numpy as np
from scipy.optimize import minimize
from itertools import combinations
import pandas as pd

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
        p = np.clip(p, 0.02, 0.98)  # Ensure frequencies are not too extreme
        
        # Apply differentiation between populations based on FST
        # We'll iteratively adjust population pairs to match their target FST
        for _ in range(20):  # Multiple passes to improve convergence
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
                            adjust = min(0.01, (target_fst - current_fst))
                            p[j] = min(0.98, p_j + adjust * (1 - p_j))
                            p[k] = max(0.02, p_k - adjust * p_k)
                        else:
                            # Decrease p_j, increase p_k
                            adjust = min(0.01, (target_fst - current_fst))
                            p[j] = max(0.02, p_j - adjust * p_j)
                            p[k] = min(0.98, p_k + adjust * (1 - p_k))
                    
                    # If current FST is too high, decrease differentiation
                    elif current_fst > target_fst * 1.1:  # Allow some tolerance
                        # Move frequencies closer to their mean
                        adjust = min(0.01, (current_fst - target_fst))
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

def fine_tune_frequencies(allele_freqs, target_fst_matrix, max_iter=500, verbose=True):
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
    verbose : bool
        Whether to print progress information
        
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
        if verbose:
            print(f"Iteration {iteration+1}, Mean absolute error: {mean_error:.4f}")
        
        if mean_error < 0.001:
            if verbose:
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

def simulate_individual_genotypes(pop_allele_freqs, pop_sizes, random_seed=None):
    """
    Simulate individual genotypes based on population allele frequencies.
    
    Parameters:
    -----------
    pop_allele_freqs : numpy.ndarray
        Matrix of allele frequencies with shape (num_variants, num_populations)
    pop_sizes : list or numpy.ndarray
        Number of individuals to simulate for each population
    random_seed : int, optional
        Seed for random number generation
        
    Returns:
    --------
    tuple
        (genotypes, population_labels, variant_ids)
        - genotypes: numpy.ndarray with shape (num_variants, total_individuals)
          Values are 0, 1, or 2 representing the count of alternative alleles
        - population_labels: list with population ID for each individual
        - variant_ids: list of variant IDs
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    num_variants, num_pops = pop_allele_freqs.shape
    
    # Validate pop_sizes
    if len(pop_sizes) != num_pops:
        raise ValueError(f"Length of pop_sizes ({len(pop_sizes)}) must match number of populations ({num_pops})")
    
    # Create variant IDs
    variant_ids = [f"var_{i+1}" for i in range(num_variants)]
    
    # Calculate total number of individuals
    total_individuals = sum(pop_sizes)
    
    # Initialize genotype matrix and population labels
    genotypes = np.zeros((num_variants, total_individuals), dtype=int)
    population_labels = []
    
    # Track individual index
    ind_idx = 0
    
    # For each population
    for pop_idx in range(num_pops):
        pop_size = pop_sizes[pop_idx]
        pop_name = f"pop_{pop_idx+1}"
        
        # Add population labels
        population_labels.extend([pop_name] * pop_size)
        
        # For each variant
        for var_idx in range(num_variants):
            # Get allele frequency for this population and variant
            p = pop_allele_freqs[var_idx, pop_idx]
            
            # Generate genotypes according to Hardy-Weinberg equilibrium
            # 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternative
            # Probabilities under HWE: [(1-p)^2, 2p(1-p), p^2]
            
            # Calculate HWE probabilities
            probs = [(1-p)**2, 2*p*(1-p), p**2]
            
            # Generate genotypes for this variant and population
            genotypes[var_idx, ind_idx:ind_idx+pop_size] = np.random.choice(
                [0, 1, 2], size=pop_size, p=probs
            )
        
        # Update individual index
        ind_idx += pop_size
    
    return genotypes, population_labels, variant_ids

def create_genotype_dataframe(genotypes, population_labels, variant_ids):
    """
    Create a pandas DataFrame from simulated genotype data.
    
    Parameters:
    -----------
    genotypes : numpy.ndarray
        Matrix of genotypes with shape (num_variants, total_individuals)
    population_labels : list
        Population assignment for each individual
    variant_ids : list
        IDs for each variant
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with variants as rows and individuals as columns
    """
    # Create individual IDs
    individual_ids = [f"{pop}_ind_{i+1}" for i, pop in enumerate(population_labels)]
    
    # Create DataFrame
    df = pd.DataFrame(
        data=genotypes,
        index=variant_ids,
        columns=individual_ids
    )
    
    # Add population information
    pop_info = pd.DataFrame({'population': population_labels}, index=individual_ids)
    
    return df, pop_info

def run_full_simulation(fst_matrix, pop_sizes, num_variants=100, random_seed=None):
    """
    Run a complete simulation from FST matrix to individual genotypes.
    
    Parameters:
    -----------
    fst_matrix : numpy.ndarray
        Square matrix of pairwise FST values between populations
    pop_sizes : list or numpy.ndarray
        Number of individuals to simulate for each population
    num_variants : int
        Number of genetic variants to simulate
    random_seed : int, optional
        Seed for random number generation
        
    Returns:
    --------
    dict
        Dictionary containing:
        - 'allele_frequencies': simulated population allele frequencies
        - 'calculated_fst': calculated FST matrix from simulated frequencies
        - 'genotypes': raw genotype matrix
        - 'genotype_df': pandas DataFrame of genotypes
        - 'population_info': DataFrame with population assignments
        - 'population_labels': list of population labels for each individual
        - 'variant_ids': list of variant IDs
    """
    if random_seed is not None:
        np.random.seed(random_seed)
    
    # Simulate population allele frequencies
    print("Simulating population allele frequencies...")
    allele_freqs = simulate_allele_frequencies(
        fst_matrix, num_variants=num_variants, random_seed=random_seed
    )
    
    # Calculate resulting FST matrix
    calculated_fst = calculate_fst(allele_freqs)
    
    # Simulate individual genotypes
    print(f"Simulating genotypes for individuals ({sum(pop_sizes)} total)...")
    genotypes, pop_labels, var_ids = simulate_individual_genotypes(
        allele_freqs, pop_sizes, random_seed=random_seed
    )
    
    # Create DataFrame representation
    genotype_df, pop_info = create_genotype_dataframe(genotypes, pop_labels, var_ids)
    
    # Calculate error in FST
    error = np.abs(fst_matrix - calculated_fst)
    off_diag_indices = np.where(~np.eye(fst_matrix.shape[0], dtype=bool))
    mean_abs_error = np.mean(error[off_diag_indices])
    print(f"Mean absolute error in FST: {mean_abs_error:.4f}")
    
    return {
        'allele_frequencies': allele_freqs,
        'calculated_fst': calculated_fst,
        'genotypes': genotypes,
        'genotype_df': genotype_df,
        'population_info': pop_info,
        'population_labels': pop_labels,
        'variant_ids': var_ids
    }

def write_simulation_to_files(simulation_results, output_prefix="simulation"):
    """
    Write simulation results to files.
    
    Parameters:
    -----------
    simulation_results : dict
        Dictionary of simulation results from run_full_simulation
    output_prefix : str
        Prefix for output files
        
    Returns:
    --------
    list
        List of file paths created
    """
    files_created = []
    
    # Write allele frequencies
    freq_file = f"{output_prefix}_allele_frequencies.csv"
    pd.DataFrame(
        simulation_results['allele_frequencies'],
        index=[f"var_{i+1}" for i in range(len(simulation_results['allele_frequencies']))],
        columns=[f"pop_{i+1}" for i in range(simulation_results['allele_frequencies'].shape[1])]
    ).to_csv(freq_file)
    files_created.append(freq_file)
    
    # Write calculated FST matrix
    fst_file = f"{output_prefix}_calculated_fst.csv"
    pd.DataFrame(
        simulation_results['calculated_fst'],
        index=[f"pop_{i+1}" for i in range(simulation_results['calculated_fst'].shape[0])],
        columns=[f"pop_{i+1}" for i in range(simulation_results['calculated_fst'].shape[1])]
    ).to_csv(fst_file)
    files_created.append(fst_file)
    
    # Write genotypes
    geno_file = f"{output_prefix}_genotypes.csv"
    simulation_results['genotype_df'].to_csv(geno_file)
    files_created.append(geno_file)
    
    # Write population information
    pop_file = f"{output_prefix}_population_info.csv"
    simulation_results['population_info'].to_csv(pop_file)
    files_created.append(pop_file)
    
    print(f"Simulation results written to {len(files_created)} files with prefix '{output_prefix}'")
    return files_created

def parse_fst_matrix(fst_file, n_pops):
    """
    Parse FST matrix from file.
    
    Parameters:
    -----------
    fst_file : str
        Path to FST matrix file
    n_pops : int
        Number of populations
        
    Returns:
    --------
    numpy.ndarray
        FST matrix
    """
    if fst_file is None:
        # Generate random FST matrix
        fst_matrix = np.zeros((n_pops, n_pops))
        for i in range(n_pops):
            for j in range(i+1, n_pops):
                fst = np.random.uniform(0.01, 0.2)  # Random FST between 0.01 and 0.2
                fst_matrix[i, j] = fst
                fst_matrix[j, i] = fst
    else:
        # Read FST matrix from file
        fst_matrix = np.loadtxt(fst_file)
        
        # Check dimensions
        if fst_matrix.shape != (n_pops, n_pops):
            raise ValueError(f"FST matrix dimensions ({fst_matrix.shape}) do not match the number of populations ({n_pops})")
    
    return fst_matrix



def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Simulate genetic variants in VCF format with population structure')
    parser.add_argument('--n-variants', type=int, default=1000, help='Number of variants to simulate')
    parser.add_argument('--n-pops', type=int, default=3, help='Number of populations')
    parser.add_argument('--samples-per-pop', type=str, default='50,50,50', help='Comma-separated list of samples per population')
    parser.add_argument('--fst-file', type=str, default=None, help='File with FST matrix (if not provided, random values will be used)')
    parser.add_argument('--output', type=str, default='simulated_variants.vcf', help='Output VCF file')
    parser.add_argument('--seed', type=int, default=None, help='Random seed')
    parser.add_argument('--check-fst', action='store_true', help='Calculate and print realized FST from simulated data')
    parser.add_argument('--max-iterations', type=int, default=5, help='Maximum number of iterations to try to match target FST')
    args = parser.parse_args()
    
    # Set random seed
    if args.seed is not None:
        np.random.seed(args.seed)
    
    # Parse samples per population
    samples_per_pop = [int(x) for x in args.samples_per_pop.split(',')]
    
    # Check if the number of specified populations matches the expected number
    if len(samples_per_pop) != args.n_pops:
        raise ValueError(f"Number of samples per population ({len(samples_per_pop)}) does not match the number of populations ({args.n_pops})")
    
    # Parse FST matrix
    target_fst_matrix = parse_fst_matrix(args.fst_file, args.n_pops)
    
    # Print simulation parameters
    print(f"Simulating {args.n_variants} variants across {args.n_pops} populations")
    print(f"Samples per population: {samples_per_pop}")
    print(f"Total samples: {sum(samples_per_pop)}")
    print(f"Target FST matrix:\n{target_fst_matrix}")
    print(f"Output file: {args.output}")
    

    # Create VCF metadata
    print("Creating VCF metadata...")
    metadata = create_vcf_metadata(sum(samples_per_pop), samples_per_pop)
    
    # Write VCF file
    print("Writing VCF file...")
    genotypes_to_vcf(genotypes, metadata, args.output)
    
    print("Simulation complete!")

# Example usage:
def example_run():
    """Example run of the complete simulation framework"""
    # Define a simple FST matrix for 3 populations
    fst_matrix = np.array([
        [0.0, 0.05, 0.15],
        [0.05, 0.0, 0.10],
        [0.15, 0.10, 0.0]
    ])
    
    # Define the number of individuals in each population
    pop_sizes = [50, 40, 30]  # Different sizes for each population
    
    print("Target FST matrix:")
    print(fst_matrix)
    print("\nPopulation sizes:", pop_sizes)
    
    # Run the full simulation
    results = run_full_simulation(fst_matrix, pop_sizes, num_variants=100, random_seed=42)
    
    print("\nCalculated FST matrix:")
    print(results['calculated_fst'])
    
    # Get some summary statistics on the genotypes
    geno_df = results['genotype_df']
    print("\nGenotype matrix shape:", results['genotypes'].shape)
    print("First 5 variants × 5 individuals:")
    print(geno_df.iloc[:5, :5])
    
    # Calculate allele frequencies from the genotypes and compare to the input frequencies
    pop_info = results['population_info']
    
    # Calculate mean genotype by population for the first 5 variants
    print("\nMean genotype by population (first 5 variants):")
    for pop in np.unique(results['population_labels']):
        ind_in_pop = pop_info[pop_info['population'] == pop].index
        pop_mean = geno_df.loc[results['variant_ids'][:5], ind_in_pop].mean(axis=1) / 2
        print(f"{pop}: {pop_mean.values}")
    
    return results

if __name__ == "__main__":
    example_run()
