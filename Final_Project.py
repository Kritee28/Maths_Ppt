import numpy as np

# Hypothetical parameters
sigma_max = 300  # Maximum allowable stress in MPa
delta_max = 0.005  # Maximum allowable displacement in meters
P_cr = 10000  # Critical buckling load in N
f_min = 5.0  # Minimum allowable natural frequency in Hz

# Hypothetical material and geometric properties for initial design
material_density = 7800  # kg/m³ for steel
initial_design = [80,.05 , .02]  # [length (m), width (m), thickness (m)]

# Objective function: minimize weight
def objective(x):
    length, width, thickness = x
    volume = length * width * thickness
    weight = volume * material_density
    return weight

# Constraint functions
def stress_constraint(x):
    stress = calculate_stress(x)  # Hypothetical stress function
    return stress - sigma_max    # Constraint: σ - σ_max <= 0

def displacement_constraint(x):
    displacement = calculate_displacement(x)  # Hypothetical displacement function
    return displacement - delta_max    # Constraint: δ - δ_max <= 0

def buckling_constraint(x):
    load = calculate_buckling_load(x)  # Hypothetical buckling load function
    return P_cr - load  # Constraint: P_cr - P <= 0

def frequency_constraint(x):
    frequency = calculate_frequency(x)  # Hypothetical natural frequency function
    return f_min - frequency   # Constraint: f_min - f <= 0

# Hypothetical functions for stress, displacement, buckling, and frequency
def calculate_stress(x):
    return 100 + 20 * x[1]  # Example stress formula

def calculate_displacement(x):
    return 0.003 + 0.002 * x[2]  # Example displacement formula

def calculate_buckling_load(x):
    return 8000 + 1000 * x[2]  # Example buckling load formula

def calculate_frequency(x):
    return 6.0 - 0.1 * x[0]  # Example frequency formula

# Barrier objective function
def barrier_objective(x, mu):
    weight = objective(x)
    
    # Adding the barrier terms for constraints
    constraints = [
        stress_constraint(x),
        displacement_constraint(x),
        buckling_constraint(x),
        frequency_constraint(x)
    ]
    
    barrier_term = 0.0
    for g in constraints:
        if g < 0:
            barrier_term -= np.log(-g)  # Only include if g < 0, ensuring feasibility

    return weight + mu * barrier_term

# Perform the optimization using the barrier method
def optimize_with_barrier(initial_design, mu=1, mu_decay=0.1, tol=1e-6, max_iter=100):
    current_design = np.array(initial_design)
    min_design_bounds = np.array([10, .02, .005])  # Minimum bounds for [length, width, thickness]
    max_design_bounds = np.array([1e5, 1e5, 1e5])     # Maximum bounds for [length, width, thickness]
    alpha = 0.05  # Step size for gradient descent

    for iteration in range(max_iter):
        # Compute gradient numerically for barrier objective
        
        grad=np.array([
             objective(current_design)/current_design[0] + mu*(1/(frequency_constraint(current_design))),
             objective(current_design)/current_design[1] + mu*(1/(stress_constraint(current_design))),
             objective(current_design)/current_design[2] + mu*((1/buckling_constraint(current_design))+ 1/displacement_constraint(current_design))   
       ])
        
        print ("current_design is  ",current_design)
        # Update design variables
        current_design -= alpha * grad
        current_design = np.clip(current_design, min_design_bounds, max_design_bounds)
        
        # Check for convergence
        grad_norm = np.linalg.norm(grad)
        if grad_norm < tol:
            print(f"Converged in {iteration + 1} iterations")
            break
        print("gradient norm is  ",grad_norm,end="\n")

        # Decrease the barrier parameter
        mu *= mu_decay
    
    return current_design, objective(current_design)

# Run the optimization
optimized_design, minimized_weight = optimize_with_barrier(initial_design)
print("Optimized Design Variables:", optimized_design)
print("Minimized Weight:", minimized_weight)