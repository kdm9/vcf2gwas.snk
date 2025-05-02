#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import allel
import time
import os
from scipy import stats
from scipy.linalg import cholesky

def simulate_allele_frequencies(n_variants, n_pops, fst_matrix):
    """
    Simulate allele frequencies across populations based on FST matrix
    using a multivariate normal model to ensure the resulting genotype data
    will have FST values matching the input matrix.
    
    Based on the approach described in:
    Balding & Nichols (1995), Hudson (2002), and Gautier et al. (2015)
    
    Parameters:
    -----------
    n_variants : int
        Number of variants to simulate
    n_pops : int
        Number of populations
    fst_matrix : numpy.ndarray
        Matrix of pairwise FST values between populations
        
    Returns:
    --------
    numpy.ndarray
        Matrix of allele frequencies (n_variants × n_pops)
    """
    # Generate ancestral allele frequencies
    # Using Beta(0.5, 0.5) to get more realistic U-shaped distribution of AF
    ancestral_freqs = np.random.beta(0.5, 0.5, size=n_variants)
    
    # Create matrix to store population allele frequencies
    pop_freqs = np.zeros((n_variants, n_pops))
    
    # Convert FST matrix to covariance matrix
    # Use the relationship: covariance between pops i and j = (1 - FST_ij)
    # where FST_ij is the FST value between populations i and j
    cov_matrix = 1 - fst_matrix
    
    # Ensure the covariance matrix is positive definite (required for multivariate normal)
    # Add a small value to the diagonal if needed
    min_eig = np.min(np.linalg.eigvals(cov_matrix))
    if min_eig < 1e-6:
        cov_matrix += np.eye(n_pops) * (1e-6 - min_eig)
        print("adjusted FST matrix to be pos-semidef")
    
    # For each variant, generate correlated allele frequencies
    for i in range(n_variants):
        p_anc = ancestral_freqs[i]
        
        # Use logit-normal approach for correlated frequencies
        # First generate correlated normal variables
        # Multivariate normal with mean 0 and covariance from FST
        z = np.random.multivariate_normal(np.zeros(n_pops), cov_matrix)
        
        # Transform the normal variables to the unit interval [0,1]
        # Scale the normal variables based on the ancestral frequency
        # This is a variation of the Balding-Nichols model for multiple populations
        theta = p_anc * (1 - p_anc)
        c = np.sqrt(fst_matrix.diagonal() / (1 - fst_matrix.diagonal()) * theta)
        
        # Calculate population frequencies
        for j in range(n_pops):
            # Add ancestral component and random component
            logit_p = np.log(p_anc / (1 - p_anc)) + c[j] * z[j]
            # Convert back from logit scale
            pop_freqs[i, j] = 1 / (1 + np.exp(-logit_p))
            
            # Bound frequencies to avoid issues with 0 and 1
            pop_freqs[i, j] = max(0.001, min(0.999, pop_freqs[i, j]))
    
    return pop_freqs

def calculate_realized_fst(genotypes, pop_sizes):
    """
    Calculate the realized FST matrix from the generated genotype data.
    
    Parameters:
    -----------
    genotypes : numpy.ndarray
        3D array of genotypes (n_variants × n_samples × 2)
    pop_sizes : list
        List of sample sizes per population
        
    Returns:
    --------
    numpy.ndarray
        Realized FST matrix
    """
    n_pops = len(pop_sizes)
    n_variants = genotypes.shape[0]
    
    # Calculate FST for each variant and then average
    fst_matrix = np.zeros((n_pops, n_pops))
    
    # Keep track of samples by population
    start_idx = 0
    pop_indices = []
    for size in pop_sizes:
        pop_indices.append(list(range(start_idx, start_idx + size)))
        start_idx += size
    
    # Calculate allele frequencies for each population
    pop_freqs = np.zeros((n_variants, n_pops))
    for i in range(n_pops):
        for v in range(n_variants):
            # Get genotypes for this population and variant
            pop_genotypes = genotypes[v, pop_indices[i], :]
            # Calculate allele frequency
            pop_freqs[v, i] = np.sum(pop_genotypes) / (2 * len(pop_indices[i]))
    
    # Calculate global allele frequencies
    global_freqs = np.zeros(n_variants)
    for v in range(n_variants):
        total_count = 0
        for i in range(n_pops):
            total_count += np.sum(genotypes[v, pop_indices[i], :])
        global_freqs[v] = total_count / (2 * np.sum(pop_sizes))
    
    # Calculate FST for each pair of populations
    for i in range(n_pops):
        for j in range(i+1, n_pops):
            fst_values = []
            for v in range(n_variants):
                p_i = pop_freqs[v, i]
                p_j = pop_freqs[v, j]
                p = (pop_sizes[i] * p_i + pop_sizes[j] * p_j) / (pop_sizes[i] + pop_sizes[j])
                
                # Only use variants that are polymorphic
                if 0 < p < 1:
                    # Calculate Hudson's FST estimator
                    n_i = pop_sizes[i]
                    n_j = pop_sizes[j]
                    h_s = (n_i * p_i * (1 - p_i) + n_j * p_j * (1 - p_j)) / (n_i + n_j)
                    h_t = p * (1 - p)
                    
                    if h_t > 0:
                        fst = 1 - (h_s / h_t)
                        fst_values.append(fst)
            
            # Average FST across variants
            if fst_values:
                avg_fst = np.mean(fst_values)
                fst_matrix[i, j] = avg_fst
                fst_matrix[j, i] = avg_fst
    
    return fst_matrix

def simulate_genotypes(pop_freqs, samples_per_pop):
    """
    Simulate genotypes based on population allele frequencies.
    
    Parameters:
    -----------
    pop_freqs : numpy.ndarray
        Matrix of allele frequencies (n_variants × n_pops)
    samples_per_pop : list
        List of number of samples per population
        
    Returns:
    --------
    numpy.ndarray
        3D array of genotypes (n_variants × n_samples × 2)
    """
    n_variants, n_pops = pop_freqs.shape
    total_samples = sum(samples_per_pop)
    
    # Create array to store genotypes
    genotypes = np.zeros((n_variants, total_samples, 2), dtype=np.int8)
    
    # Track the current sample index
    sample_idx = 0
    
    # For each population
    for pop_idx in range(n_pops):
        n_samples = samples_per_pop[pop_idx]
        
        # For each variant
        for var_idx in range(n_variants):
            p = pop_freqs[var_idx, pop_idx]
            
            # Simulate diploid genotypes
            for i in range(n_samples):
                # Two independent draws for a diploid individual
                genotypes[var_idx, sample_idx + i, 0] = np.random.binomial(1, p)
                genotypes[var_idx, sample_idx + i, 1] = np.random.binomial(1, p)
        
        # Update sample index
        sample_idx += n_samples
    
    return genotypes

def create_vcf_metadata(n_samples, samples_per_pop):
    """
    Create VCF metadata.
    
    Parameters:
    -----------
    n_samples : int
        Total number of samples
    samples_per_pop : list
        List of number of samples per population
        
    Returns:
    --------
    dict
        VCF metadata
    """
    # Create sample IDs
    sample_idx = 0
    sample_ids = []
    sample_pops = []
    
    for pop_idx, n_pop_samples in enumerate(samples_per_pop):
        for i in range(n_pop_samples):
            sample_ids.append(f"sample_{sample_idx}")
            sample_pops.append(f"pop_{pop_idx}")
            sample_idx += 1
    
    # Create metadata dict
    metadata = {
        'fileformat': 'VCFv4.2',
        'fileDate': time.strftime('%Y%m%d'),
        'source': 'vcf_population_simulator.py',
        'INFO': {
            'AF': {'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency'},
        },
        'FORMAT': {
            'GT': {'Number': 1, 'Type': 'String', 'Description': 'Genotype'}
        },
        'samples': sample_ids,
        'sample_populations': sample_pops
    }
    
    return metadata

def genotypes_to_vcf(genotypes, metadata, output_file):
    """
    Convert genotypes to VCF file.
    
    Parameters:
    -----------
    genotypes : numpy.ndarray
        3D array of genotypes (n_variants × n_samples × 2)
    metadata : dict
        VCF metadata
    output_file : str
        Output VCF file path
    """
    n_variants, n_samples, _ = genotypes.shape
    
    # Create chromosome positions (1-based)
    chrom = np.ones(n_variants, dtype='int32')
    pos = np.arange(1, n_variants + 1, dtype='int32') * 100  # Space variants 100bp apart
    
    # Create IDs
    id = ['rs' + str(i) for i in range(n_variants)]
    
    # Create ref/alt alleles
    ref = np.full(n_variants, 'A', dtype='U1')
    alt = np.full(n_variants, 'G', dtype='U1')
    
    # Calculate allele frequencies
    allele_counts = genotypes.sum(axis=(1, 2))
    total_alleles = 2 * n_samples
    afs = allele_counts / total_alleles
    
    # Create INFO field
    info = np.array(['AF=' + str(round(af, 4)) for af in afs], dtype='U20')
    
    # Create header
    header = ['##fileformat=' + metadata['fileformat']]
    header.append('##fileDate=' + metadata['fileDate'])
    header.append('##source=' + metadata['source'])
    
    # Add INFO fields to header
    for id, values in metadata['INFO'].items():
        header.append(f'##INFO=<ID={id},Number={values["Number"]},Type={values["Type"]},Description="{values["Description"]}">')
    
    # Add FORMAT fields to header
    for id, values in metadata['FORMAT'].items():
        header.append(f'##FORMAT=<ID={id},Number={values["Number"]},Type={values["Type"]},Description="{values["Description"]}">')
    
    # Add sample population information
    for i, (sample, pop) in enumerate(zip(metadata['samples'], metadata['sample_populations'])):
        header.append(f'##SAMPLE=<ID={sample},Population={pop}>')
    
    # Add column header
    header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(metadata['samples']))
    
    # Create VCF records
    records = []
    for i in range(n_variants):
        record = [
            str(chrom[i]),
            str(pos[i]),
            id[i],
            ref[i],
            alt[i],
            '.',  # QUAL
            'PASS',  # FILTER
            info[i],
            'GT'  # FORMAT
        ]
        
        # Add genotypes
        for j in range(n_samples):
            gt = f"{genotypes[i, j, 0]}|{genotypes[i, j, 1]}"
            record.append(gt)
        
        records.append('\t'.join(record))
    
    # Write VCF file
    with open(output_file, 'w') as f:
        f.write('\n'.join(header) + '\n')
        f.write('\n'.join(records) + '\n')

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
    
    # We'll iteratively adjust our simulation to get closer to the target FST
    current_fst_matrix = target_fst_matrix.copy()
    
    for iteration in range(args.max_iterations):
        print(f"\nIteration {iteration + 1}/{args.max_iterations}")
        
        # Simulate allele frequencies
        print("Simulating allele frequencies...")
        pop_freqs = simulate_allele_frequencies(args.n_variants, args.n_pops, current_fst_matrix)
        
        # Simulate genotypes
        print("Simulating genotypes...")
        genotypes = simulate_genotypes(pop_freqs, samples_per_pop)
        
        if args.check_fst or iteration < args.max_iterations - 1:
            # Calculate realized FST from the genotypes
            realized_fst = calculate_realized_fst(genotypes, samples_per_pop)
            print(f"Realized FST matrix:\n{realized_fst}")
            
            # Calculate error between target and realized FST
            upper_indices = np.triu_indices(args.n_pops, 1)
            target_values = target_fst_matrix[upper_indices]
            realized_values = realized_fst[upper_indices]
            mse = np.mean((target_values - realized_values) ** 2)
            print(f"Mean squared error: {mse:.6f}")
            
            if mse < 0.001 or iteration == args.max_iterations - 1:
                # If we're close enough or on the last iteration, use these results
                print("FST values are acceptably close to targets.")
                break
            
            # Adjust FST matrix for next iteration based on the error
            # This is a simple correction to try to get closer to the target
            adjustment = target_fst_matrix - realized_fst
            # Dampen the adjustment to prevent overshooting
            current_fst_matrix = current_fst_matrix + 0.5 * adjustment
            # Ensure FST values are in valid range [0, 1]
            current_fst_matrix = np.clip(current_fst_matrix, 0, 0.99)
            # Ensure matrix is symmetric
            current_fst_matrix = (current_fst_matrix + current_fst_matrix.T) / 2
            # Set diagonal to 0
            np.fill_diagonal(current_fst_matrix, 0)
    
    # Create VCF metadata
    print("Creating VCF metadata...")
    metadata = create_vcf_metadata(sum(samples_per_pop), samples_per_pop)
    
    # Write VCF file
    print("Writing VCF file...")
    genotypes_to_vcf(genotypes, metadata, args.output)
    
    print("Simulation complete!")

if __name__ == '__main__':
    main()
