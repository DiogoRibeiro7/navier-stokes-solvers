#!/usr/bin/env python3
"""
Visualization script for Navier-Stokes solver results
Author: Diogo Ribeiro (dfr@esmad.ipp.pt)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import os

def load_solution(filename):
    """Load solution data from file."""
    try:
        data = np.loadtxt(filename, comments='#')
        return data
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        sys.exit(1)

def plot_velocity_field(data, output_file='velocity_field.png'):
    """Plot velocity magnitude and streamlines."""
    x = data[:, 0]
    y = data[:, 1]
    u = data[:, 2]
    v = data[:, 3]
    
    # Determine grid dimensions
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    
    # Reshape data
    X = x.reshape((ny, nx))
    Y = y.reshape((ny, nx))
    U = u.reshape((ny, nx))
    V = v.reshape((ny, nx))
    
    # Compute velocity magnitude
    vel_mag = np.sqrt(U**2 + V**2)
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Velocity magnitude
    im1 = ax1.contourf(X, Y, vel_mag, levels=20, cmap='viridis')
    ax1.quiver(X[::2, ::2], Y[::2, ::2], U[::2, ::2], V[::2, ::2], 
               color='white', alpha=0.6)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title('Velocity Magnitude')
    ax1.set_aspect('equal')
    plt.colorbar(im1, ax=ax1, label='|u|')
    
    # Plot 2: Streamlines
    ax2.streamplot(X, Y, U, V, color=vel_mag, cmap='plasma', 
                   density=2, linewidth=1)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title('Streamlines')
    ax2.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_vorticity(data, output_file='vorticity.png'):
    """Plot vorticity field."""
    if data.shape[1] < 5:
        print("Warning: No vorticity data in file")
        return
    
    x = data[:, 0]
    y = data[:, 1]
    omega = data[:, 4]  # Assuming vorticity is 5th column
    
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    
    X = x.reshape((ny, nx))
    Y = y.reshape((ny, nx))
    Omega = omega.reshape((ny, nx))
    
    plt.figure(figsize=(8, 7))
    levels = np.linspace(Omega.min(), Omega.max(), 30)
    im = plt.contourf(X, Y, Omega, levels=levels, cmap='RdBu_r')
    plt.colorbar(im, label='ω')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Vorticity Field')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_centerline_profiles(data, output_file='profiles.png'):
    """Plot velocity profiles along centerlines."""
    x = data[:, 0]
    y = data[:, 1]
    u = data[:, 2]
    v = data[:, 3]
    
    nx = len(np.unique(x))
    ny = len(np.unique(y))
    
    X = x.reshape((ny, nx))
    Y = y.reshape((ny, nx))
    U = u.reshape((ny, nx))
    V = v.reshape((ny, nx))
    
    # Get centerline indices
    mid_x = nx // 2
    mid_y = ny // 2
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Vertical centerline (u velocity)
    ax1.plot(U[:, mid_x], Y[:, mid_x], 'b-', linewidth=2, label='u velocity')
    ax1.set_xlabel('u')
    ax1.set_ylabel('y')
    ax1.set_title('Vertical Centerline Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Horizontal centerline (v velocity)
    ax2.plot(X[mid_y, :], V[mid_y, :], 'r-', linewidth=2, label='v velocity')
    ax2.set_xlabel('x')
    ax2.set_ylabel('v')
    ax2.set_title('Horizontal Centerline Profile')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_solution.py <solution_file.dat> [output_prefix]")
        print("\nExample:")
        print("  python plot_solution.py fd_final_solution.dat")
        print("  python plot_solution.py spectral_final_solution.dat results/spectral")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(input_file)[0]
    
    print(f"Loading solution from: {input_file}")
    data = load_solution(input_file)
    print(f"Data shape: {data.shape}")
    print(f"Columns: x, y, u, v, [omega], [psi], [p]")
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate plots
    print("\nGenerating plots...")
    plot_velocity_field(data, f"{output_prefix}_velocity.png")
    plot_centerline_profiles(data, f"{output_prefix}_profiles.png")
    
    # Plot vorticity if available
    if data.shape[1] >= 5:
        plot_vorticity(data, f"{output_prefix}_vorticity.png")
    
    print("\n✓ All plots generated successfully!")
    print(f"\nView results:")
    print(f"  {output_prefix}_velocity.png")
    print(f"  {output_prefix}_profiles.png")
    if data.shape[1] >= 5:
        print(f"  {output_prefix}_vorticity.png")

if __name__ == "__main__":
    main()
