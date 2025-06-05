import numpy as np
import scipy.optimize as optimize

class GRMSimulator:
    def __init__(self, target_grm):
        """
        Initialize simulator with target Genotype Relationship Matrix.
        
        Parameters:
        - target_grm (numpy.ndarray): Target symmetric GRM
        """
        self.target_grm = target_grm
        self.num_individuals, self.num_samples = target_grm.shape
        
    def _calculate_grm(self, genotypes):
        """
        Calculate GRM from genotype matrix.
        
        Parameters:
        - genotypes (numpy.ndarray): Genotype matrix
        
        Returns:
        numpy.ndarray: Calculated GRM
        """
        # Center genotypes
        centered_genotypes = genotypes - np.mean(genotypes, axis=0)
        
        # Normalize and compute relationship matrix
        normalized_genotypes = np.nan_to_num(centered_genotypes / np.std(centered_genotypes, axis=0))
        return np.dot(normalized_genotypes, normalized_genotypes.T) / genotypes.shape[1]
    
    def _objective_function(self, flat_genotypes):
        """
        Compute difference between target and simulated GRM.
        
        Parameters:
        - flat_genotypes (numpy.ndarray): Flattened genotype matrix
        
        Returns:
        float: Frobenius norm of GRM difference
        """
        # Reshape flattened genotypes
        genotypes = flat_genotypes.reshape(self.num_individuals, -1)
        
        # Ensure binary genotypes (0, 1, 2)
        genotypes = np.round(np.clip(genotypes, 0, 2)).astype(int)
        
        # Calculate GRM
        simulated_grm = self._calculate_grm(genotypes)
        
        # Compute difference
        return np.linalg.norm(self.target_grm - simulated_grm)
    
    def simulate_variants(self, num_variants, initial_guess=None, max_iterations=1000):
        """
        Simulate genetic variants to match target GRM.
        
        Parameters:
        - num_variants (int): Number of variants to simulate
        - initial_guess (numpy.ndarray, optional): Initial genotype matrix
        - max_iterations (int): Maximum optimization iterations
        
        Returns:
        dict: Simulation results
        """
        # Initialize genotypes if no initial guess
        if initial_guess is None:
            initial_guess = np.random.randint(0, 3, 
                size=(self.num_individuals, num_variants)).astype(float)
        
        # Flatten for optimization
        initial_flat = initial_guess.flatten()
        
        # Bounds to ensure genotypes are 0, 1, or 2
        bounds = [(0, 2)] * (self.num_individuals * num_variants)
        
        # Optimize genotypes
        result = optimize.minimize(
            self._objective_function, 
            initial_flat, 
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': max_iterations}
        )
        
        # Reshape and round optimized genotypes
        optimized_genotypes = np.round(
            result.x.reshape(self.num_individuals, num_variants)
        ).astype(int)
        
        # Calculate final GRM
        final_grm = self._calculate_grm(optimized_genotypes)
        
        return {
            'genotypes': optimized_genotypes,
            'simulated_grm': final_grm,
            'grm_difference': result.fun,
            'optimization_success': result.success
        }

# Example usage
def main():
    # Generate a sample target GRM
    np.random.seed(42)
    sample_target_grm = np.random.rand(10, 10)
    sample_target_grm = (sample_target_grm + sample_target_grm.T) / 2
    sample_target_grm = np.dot(sample_target_grm, sample_target_grm.T)
    sample_target_grm /= np.trace(sample_target_grm)
    
    # Initialize simulator
    simulator = GRMSimulator(sample_target_grm)
    
    # Simulate variants
    simulation_results = simulator.simulate_variants(num_variants=500)
    
    # Print results
    print("GRM Difference:", simulation_results['grm_difference'])
    print("Optimization Success:", simulation_results['optimization_success'])

    print(sample_target_grm)
    print(simulation_results['simulated_grm'])

if __name__ == "__main__":
    main()
