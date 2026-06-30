import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── Load force data ───────────────────────────────────────────
df = pd.read_csv('forces.txt')

# ── Plot Cl and Cd vs iteration ───────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)

ax1.plot(df['iter'], df['Cl'], color='royalblue', linewidth=1.5)
ax1.set_ylabel('Cl', fontweight='bold')
ax1.set_title('NACA 0012 — Force Coefficients vs Iteration', fontweight='bold')
ax1.grid(True, linestyle='--', alpha=0.6)

ax2.plot(df['iter'], df['Cd'], color='firebrick', linewidth=1.5)
ax2.set_xlabel('Iteration', fontweight='bold')
ax2.set_ylabel('Cd', fontweight='bold')
ax2.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig('forces_history.png', dpi=150, bbox_inches='tight')
plt.show()

# ── Converged values (mean of last 20% of data) ───────────────
last = df.iloc[int(0.8 * len(df)):]
Cl_mean = last['Cl'].mean()
Cd_mean = last['Cd'].mean()

print(f"\n=== Converged Results ===")
print(f"Cl        = {Cl_mean:.4f}")
print(f"Cd        = {Cd_mean:.4f}")
print(f"Cl/Cd     = {Cl_mean/Cd_mean:.2f}")

# ── Thin airfoil theory comparison ───────────────────────────
AoA_deg = 5.0
Cl_theory = 2.0 * np.pi * np.sin(np.radians(AoA_deg))
print(f"\n=== Validation ===")
print(f"Thin airfoil theory Cl = {Cl_theory:.4f}")
print(f"LBM Cl                 = {Cl_mean:.4f}")
print(f"Difference             = {abs(Cl_mean - Cl_theory):.4f}")

# ── Plot airfoil geometry (run this separately to verify geometry) ──
import pandas as pd
import matplotlib.pyplot as plt

geom = pd.read_csv('airfoil_geometry.csv')

fig, ax = plt.subplots(figsize=(10, 3))
ax.plot(geom['x_upper'], geom['y_upper'], 'b-', label='Upper surface')
ax.plot(geom[' y_lower'], geom[' y_lower'], 'r-', label='Lower surface')
ax.set_aspect('equal')
ax.legend()
ax.set_title('NACA 0012 Geometry (Lattice Units)', fontweight='bold')
ax.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig('airfoil_geometry.png', dpi=150)
plt.show()