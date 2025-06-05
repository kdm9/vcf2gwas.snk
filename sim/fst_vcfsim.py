#!/usr/bin/env python
import numpy as np
from scipy.optimize import minimize
from itertools import combinations
import pandas as pd
import argparse
import sys
import os
import datetime

def simulate_allele_frequencies(fst_matrix, num_variants=100, random_seed=None, verbose=False):
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
    
    num_pops = fst_matrix.shape[0]
    allele_freqs = np.zeros((num_variants, num_pops))
    
    for i in range(num_variants):
        # Start with a random ancestral allele frequency, beta dist (heavily U-shaped) centered around 0.5
        p_ancestral = np.random.beta(0.5, 0.5)  
        
        # Initial allele frequencies - start all populations with the same frequency
        # but with small random perturbations
        p = np.array([p_ancestral + np.random.normal(0, 0.01) for _ in range(num_pops)])
        p = np.clip(p, 0.002, 0.998)  # Ensure frequencies are not too extreme
        
        # Apply differentiation between populations based on FST
        # We'll iteratively adjust population pairs to match their target FST
        for _ in range(10):  # Multiple passes to improve convergence
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
                            adjust = 0.01 #min(0.02, (target_fst - current_fst))
                            p[j] = min(0.98, p_j + adjust * (1 - p_j))
                            p[k] = max(0.02, p_k - adjust * p_k)
                        else:
                            # Decrease p_j, increase p_k
                            adjust = 0.01 #min(0.02, (target_fst - current_fst))
                            p[j] = max(0.02, p_j - adjust * p_j)
                            p[k] = min(0.98, p_k + adjust * (1 - p_k))
                    
                    # If current FST is too high, decrease differentiation
                    elif current_fst > target_fst * 1.2:  # Allow some tolerance
                        # Move frequencies closer to their mean
                        adjust = 0.01# min(0.02, (current_fst - target_fst))
                        p[j] = p_j + adjust * (p_mean - p_j)
                        p[k] = p_k + adjust * (p_mean - p_k)
        
        # Store the frequencies
        allele_freqs[i] = p
    
    if verbose:
        print("Initial pop. freqs:")
        print(allele_freqs)
    # Fine-tune allele frequencies to better match target FST
    tuned_freqs = fine_tune_frequencies(allele_freqs, fst_matrix, verbose=verbose)
    
    return tuned_freqs


def print_hist(histo, indent=""):
    counts, bins = histo
    binrest = bins[1:]
    print(indent, f"{bins[0]:0.2f}".ljust(6))
    percs = counts / counts.sum() * 100
    for b, c in zip(binrest, percs):
        print(indent, f"{b:0.2f}".ljust(6), "#" * int(c))

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

def fine_tune_frequencies(allele_freqs, target_fst_matrix, max_iter=5000, verbose=False):
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
    
    last_mean_error = 1000000
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
        if last_mean_error - mean_error < 1e-6:
            if verbose:
                print("Mean error converged, stopping iterations")
            break
        last_mean_error = mean_error
        
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

def run_full_simulation(fst_matrix, pop_sizes, num_variants=100, random_seed=None, verbose=False):
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
    verbose : bool
        Whether to print progress information
        
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
    if verbose:
        print("Simulating population allele frequencies...")
    allele_freqs = simulate_allele_frequencies(
        fst_matrix, num_variants=num_variants, random_seed=random_seed, verbose=verbose,
    )
    
    # Calculate resulting FST matrix
    calculated_fst = calculate_fst(allele_freqs)

    if verbose:
        print(f"Simulated pop freqs yield this FST matrix:\n{calculated_fst}")
        print(f"c.f. target FST matrix:\n{fst_matrix}")
    
    # Simulate individual genotypes
    if verbose:
        print(f"Simulating genotypes for individuals ({sum(pop_sizes)} total)...")
    genotypes, pop_labels, var_ids = simulate_individual_genotypes(
        allele_freqs, pop_sizes, random_seed=random_seed
    )
    af = genotypes.sum(1)/(genotypes.shape[1]*2)
    af_hist = np.histogram(af, bins=np.arange(0, 1.1, 0.1))
    if verbose:
        print("Allele frequency histogram:")
        print_hist(af_hist, indent="  ")
        print()
    
    # Create DataFrame representation
    genotype_df, pop_info = create_genotype_dataframe(genotypes, pop_labels, var_ids)
    
    # Calculate error in FST
    error = np.abs(fst_matrix - calculated_fst)
    off_diag_indices = np.where(~np.eye(fst_matrix.shape[0], dtype=bool))
    mean_abs_error = np.mean(error[off_diag_indices])
    if verbose:
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

def write_simulation_to_files(simulation_results, output_prefix="simulation", verbose=False):
    """
    Write simulation results to files.
    
    Parameters:
    -----------
    simulation_results : dict
        Dictionary of simulation results from run_full_simulation
    output_prefix : str
        Prefix for output files
    verbose : bool
        Whether to print progress information
        
    Returns:
    --------
    list
        List of file paths created
    """
    files_created = []
    
    freq_file = f"{output_prefix}_allele_frequencies.csv"
    pd.DataFrame(
        simulation_results['allele_frequencies'],
        index=[f"var_{i+1}" for i in range(len(simulation_results['allele_frequencies']))],
        columns=[f"pop_{i+1}" for i in range(simulation_results['allele_frequencies'].shape[1])]
    ).to_csv(freq_file)
    files_created.append(freq_file)
    
    fst_file = f"{output_prefix}_calculated_fst.csv"
    pd.DataFrame(
        simulation_results['calculated_fst'],
        index=[f"pop_{i+1}" for i in range(simulation_results['calculated_fst'].shape[0])],
        columns=[f"pop_{i+1}" for i in range(simulation_results['calculated_fst'].shape[1])]
    ).to_csv(fst_file)
    files_created.append(fst_file)
    
    pop_file = f"{output_prefix}_population_info.csv"
    simulation_results['population_info'].to_csv(pop_file)
    files_created.append(pop_file)
    
    vcf_file = f"{output_prefix}.vcf"
    write_vcf(simulation_results, vcf_file)
    files_created.append(vcf_file)
    
    if verbose:
        print(f"Simulation results written to {len(files_created)} files with prefix '{output_prefix}'")
    return files_created

def write_vcf(simulation_results, output_file):
    """
    Write simulation results to a VCF file.
    
    Parameters:
    -----------
    simulation_results : dict
        Dictionary of simulation results from run_full_simulation
    output_file : str
        Path to output VCF file
    """
    genotypes = simulation_results['genotypes']
    variant_ids = simulation_results['variant_ids']
    population_labels = simulation_results['population_labels']
    
    # Create individual IDs without population prefix for VCF
    individual_ids = [f"ind_{i+1}" for i in range(len(population_labels))]
    
    with open(output_file, 'w') as f:
        # Write VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.datetime.now().strftime('%Y%m%d')}\n")
        f.write("##source=PopulationSimulator\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        
        # Add contig information (using chromosome 1 for all variants)
        f.write("##contig=<ID=1,length=249250621>\n")

        # Add population information
        for i, pop_id in enumerate(set(population_labels)):
            f.write(f"##POPULATION=<ID={pop_id},Description=\"Simulated population {i+1}\">\n")
            
        # Add sample information
        for i, (ind, pop) in enumerate(zip(individual_ids, population_labels)):
            f.write(f"##SAMPLE=<ID={ind},Population={pop}>\n")
        
        # Write header line
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + individual_ids
        f.write("\t".join(header) + "\n")
        
        # Write variant information
        for i, var_id in enumerate(variant_ids):
            # Create random position for variant (since we're using artificial data)
            pos = i * 100 + 1000  # arbitrary spacing
            
            # Create genotype strings
            gt_strings = []
            for gt in genotypes[i]:
                if gt == 0:
                    gt_string = "0|0"  # Homozygous reference
                elif gt == 1:
                    gt_string = "0|1"  # Heterozygous
                else:  # gt == 2
                    gt_string = "1|1"  # Homozygous alternate
                gt_strings.append(gt_string)
            
            # Create line
            line = [
                "1",                # CHROM (using chromosome 1 for all variants)
                str(pos),           # POS
                var_id,             # ID
                "A",                # REF (arbitrary reference allele)
                "T",                # ALT (arbitrary alternate allele)
                ".",                # QUAL (no quality score)
                "PASS",             # FILTER
                ".",                # INFO (no additional info)
                "GT",               # FORMAT (only genotype)
            ] + gt_strings
            
            f.write("\t".join(line) + "\n")

def parse_fst_matrix(fst_file=None):
    """
    Parse a FST matrix file
    
    Parameters:
    -----------
    fst_file : str
        Path to FST matrix file
        
    Returns:
    --------
    numpy.ndarray
        FST matrix as a numpy array
    """
    try:
        if fst_file is None:
            # Generate random FST matrix
            matrix = np.zeros((n_pops, n_pops))
            for i in range(n_pops):
                for j in range(i+1, n_pops):
                    fst = np.random.uniform(0.01, 0.2)  # Random FST between 0.01 and 0.2
                    matrix[i, j] = fst
                    matrix[j, i] = fst
        else:
            # Read FST matrix from file
            matrix = np.loadtxt(fst_file)
        
        # Check if matrix is square
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("FST matrix must be square")
        
        # Check if matrix is symmetric
        if not np.allclose(matrix, matrix.T):
            print("Warning: FST matrix is not symmetric. Using (matrix + matrix.T)/2")
            matrix = (matrix + matrix.T) / 2
        
        # Check values are between 0 and 1
        if np.any(matrix < 0) or np.any(matrix > 1):
            raise ValueError("FST values must be between 0 and 1")
        
        # Check diagonal is all zeros
        np.fill_diagonal(matrix, 0)
        
        return matrix
    
    except Exception as e:
        raise ValueError(f"Error parsing FST matrix: {str(e)}")

def parse_pop_sizes(sizes_str):
    """
    Parse a string representation of population sizes.
    
    Parameters:
    -----------
    sizes_str : str
        String representation of population sizes, e.g., "50,40,30"
        
    Returns:
    --------
    list
        List of population sizes
    """
    try:
        sizes = [int(x) for x in sizes_str.strip().split(',')]
    except Exception as e:
        raise ValueError(f"Error parsing population sizes: {str(e)}")
    if any(size <= 0 for size in sizes):
        raise ValueError("Population sizes must be positive integers")
    return sizes

def main():
    parser = argparse.ArgumentParser(
        description="Simulate genetic data for multiple populations based on FST values."
    )
    
    # Required arguments
    parser.add_argument(
        "--fst", "-f", required=True,
        help="FST matrix file"
    )
    parser.add_argument(
        "--pop-sizes", "-p", required=True,
        help="Population sizes as a comma-separated string, e.g., '50,40,30'"
    )
    
    # Optional arguments
    parser.add_argument(
        "--num-variants", "-n", type=int, default=10000,
        help="Number of genetic variants to simulate (default: 100)"
    )
    parser.add_argument(
        "--output-prefix", "-o", type=str, default="simulation",
        help="Prefix for output files (default: 'simulation')"
    )
    parser.add_argument(
        "--seed", "-s", type=int, default=None,
        help="Random seed for reproducibility (default: None)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print detailed progress information"
    )
    
    args = parser.parse_args()
    
    # Parse FST matrix and population sizes
    fst_matrix = parse_fst_matrix(args.fst)
    pop_sizes = parse_pop_sizes(args.pop_sizes)
    
    # Check if number of populations matches between FST matrix and population sizes
    if len(pop_sizes) != fst_matrix.shape[0]:
        raise ValueError(f"Number of populations in FST matrix ({fst_matrix.shape[0]}) "
                       f"does not match number of population sizes ({len(pop_sizes)})")
    
    if args.verbose:
        print(f"FST matrix ({fst_matrix.shape[0]} populations):")
        print(fst_matrix)
        print(f"Population sizes: {pop_sizes}")
        print(f"Number of variants: {args.num_variants}")
        print(f"Output prefix: {args.output_prefix}")
        print(f"Random seed: {args.seed}")
    
    # Run simulation
    simulation_results = run_full_simulation(
        fst_matrix=fst_matrix,
        pop_sizes=pop_sizes,
        num_variants=args.num_variants,
        random_seed=args.seed,
        verbose=args.verbose
    )
    
    # Write results to files
    files_created = write_simulation_to_files(
        simulation_results=simulation_results,
        output_prefix=args.output_prefix,
        verbose=args.verbose
    )
    print("Simulation completed successfully.")
        
if __name__ == "__main__":
    main()
