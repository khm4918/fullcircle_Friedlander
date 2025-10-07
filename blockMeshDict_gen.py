import numpy as np
import matplotlib.pyplot as plt
import os

params = {
   'axes.labelsize': 20,
   'font.size': 20,
   'legend.fontsize': 12,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'figure.figsize': [10,10]
   } 
plt.rcParams.update(params)

#Plot figure of vertices
showfig = True

# Parameters
r0 = 0.5        # starting radius (mm)
r_max = 18     # max radius (mm)
theta_deg = 90  # wedge angle in degrees
n_radial = 1000  # radial cells
n_theta = 30    # angular divisions
n_axial = 1     # axial divisions
z_min = 0.0
z_max = 1.0

# Derived parameters
theta_rad = np.deg2rad(theta_deg)
delta_theta = theta_rad / n_theta
delta_z = (z_max - z_min) / n_axial

# Generate radii so that cells are approx square
radii = [r0]
for _ in range(n_radial):
    r_next = radii[-1] * (1 + delta_theta)
    if r_next > r_max:
        break
    radii.append(r_next)
n_r = len(radii) - 1

# Create vertex list
vertices = []
for z_i in range(n_axial + 1):
    z = z_min + z_i * delta_z
    for i_theta in range(n_theta + 1):
        theta = i_theta * delta_theta
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        for r in radii:
            x = r * cos_t
            y = r * sin_t
            vertices.append((x, y, z))

if showfig == True:
    x = []
    y = []
    
    for j in range(len(vertices)):
        x.append(vertices[j][0])
        y.append(vertices[j][1])
        
    x = np.array(x)
    y = np.array(y)
    
    plt.plot(x,y,'x')
    plt.xlim(-1, 20)
    plt.ylim(-1, 20)
    plt.show()

# Index helper
def idx(z, t, r):
    return z * (n_theta + 1) * (n_r + 1) + t * (n_r + 1) + r

# Write blockMeshDict
os.makedirs("system", exist_ok=True)
with open("system/blockMeshDict", "w") as f:
    f.write("FoamFile\n{\n    version 2.0;\n    format ascii;\n    class dictionary;\n    object blockMeshDict;\n}\n\n")
    f.write("convertToMeters 0.001;\n\n")

    f.write("vertices\n(\n")
    for v in vertices:
        f.write(f"    ({v[0]:.6f} {v[1]:.6f} {v[2]:.6f})\n")
    f.write(");\n\n")

    f.write("blocks\n(\n")
    for i_theta in range(n_theta):
        for i_r in range(n_r):
            v0 = idx(0, i_theta, i_r)
            v1 = idx(0, i_theta+1, i_r)
            v2 = idx(0, i_theta+1, i_r+1)
            v3 = idx(0, i_theta, i_r+1)
            v4 = idx(1, i_theta, i_r)
            v5 = idx(1, i_theta+1, i_r)
            v6 = idx(1, i_theta+1, i_r+1)
            v7 = idx(1, i_theta, i_r+1)

            # ✅ Corrected vertex order to prevent inside-out error
            f.write(f"    hex ({v0} {v3} {v2} {v1} {v4} {v7} {v6} {v5}) "
                    f"({1} {1} {n_axial}) simpleGrading (1 1 1)\n")
    f.write(");\n\n")

    f.write("edges\n(\n);\n\n")

    # Define boundaries
    f.write("boundary\n(\n")
    
    # symmetryPlane along theta = 0
    f.write("    symmetryX\n    {\n        type symmetryPlane;\n        faces\n        (\n")
    for i_r in range(n_r):
        v0 = idx(0, 0, i_r)
        v1 = idx(0, 0, i_r+1)
        v2 = idx(1, 0, i_r+1)
        v3 = idx(1, 0, i_r)
        f.write(f"            ({v0} {v1} {v2} {v3})\n")
    f.write("        );\n    }\n")
    
    # symmetryPlane along theta = 90 deg
    f.write("    symmetryY\n    {\n        type symmetryPlane;\n        faces\n        (\n")
    for i_r in range(n_r):
        v0 = idx(0, n_theta, i_r)
        v1 = idx(0, n_theta, i_r+1)
        v2 = idx(1, n_theta, i_r+1)
        v3 = idx(1, n_theta, i_r)
        f.write(f"            ({v0} {v1} {v2} {v3})\n")
    f.write("        );\n    }\n")
    
    # empty front and back patches
    f.write("    front\n    {\n        type empty;\n        faces\n        (\n")
    for i_theta in range(n_theta):
        for i_r in range(n_r):
            v0 = idx(1, i_theta, i_r)
            v1 = idx(1, i_theta+1, i_r)
            v2 = idx(1, i_theta+1, i_r+1)
            v3 = idx(1, i_theta, i_r+1)
            f.write(f"            ({v0} {v1} {v2} {v3})\n")
    f.write("        );\n    }\n")
    
    f.write("    back\n    {\n        type empty;\n        faces\n        (\n")
    for i_theta in range(n_theta):
        for i_r in range(n_r):
            v0 = idx(0, i_theta, i_r)
            v1 = idx(0, i_theta+1, i_r)
            v2 = idx(0, i_theta+1, i_r+1)
            v3 = idx(0, i_theta, i_r+1)
            f.write(f"            ({v0} {v1} {v2} {v3})\n")
    f.write("        );\n    }\n")
    
    # Wall at inner arc
    f.write("    innerWall\n    {\n        type wall;\n        faces\n        (\n")
    for i_theta in range(n_theta):
        v0 = idx(0, i_theta, 0)
        v1 = idx(0, i_theta+1, 0)
        v2 = idx(1, i_theta+1, 0)
        v3 = idx(1, i_theta, 0)
        f.write(f"            ({v0} {v1} {v2} {v3})\n")
    f.write("        );\n    }\n")
    
    # Wall at outer arc
    f.write("    outerWall\n    {\n        type wall;\n        faces\n        (\n")
    for i_theta in range(n_theta):
        v0 = idx(0, i_theta, n_r)
        v1 = idx(0, i_theta+1, n_r)
        v2 = idx(1, i_theta+1, n_r)
        v3 = idx(1, i_theta, n_r)
        f.write(f"            ({v0} {v1} {v2} {v3})\n")
    f.write("        );\n    }\n")
    
    f.write(");\n\n")
    f.write("mergePatchPairs\n(\n);\n")

print("✅ blockMeshDict written to system/blockMeshDict")