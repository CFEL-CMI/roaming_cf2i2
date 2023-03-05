from coordinates import init_coord, change_coord, writeXYZ, atomLabels
import numpy as np
import itertools


if __name__ == "__main__":

    (rCF_ref, rCI_ref, aFCF_ref, aICI_ref), xyz_ref = init_coord()
    r = rCI_ref
    t = 0
    p = 0
    writeXYZ("cf2i2_3d_grid1.xyz", atomLabels, xyz_ref, comment=f"reference geometry: r={r}, theta={t}, phi={p}")

    r_grid = (np.linspace(1.0, 5.0, 30), np.linspace(5.0, 10.0, 20))
    t_grid = (np.linspace(0.0, np.pi, 30), np.linspace(0.0, np.pi, 20))
    p_grid = (np.linspace(0.0, 2*np.pi, 30), np.linspace(0.0, 2*np.pi, 20))

    rtp_grid = [[elem for elem in itertools.product(r, t, p)]\
                for r, t, p in zip(r_grid, t_grid, p_grid)]

    for rtp in rtp_grid:
        for (r, t, p) in rtp:
            xyz = change_coord(xyz_ref, r, t, p)
            writeXYZ("cf2i2_3d_grid1.xyz", atomLabels, xyz, comment=f"r={r}, theta={t}, phi={p}", append=True)
