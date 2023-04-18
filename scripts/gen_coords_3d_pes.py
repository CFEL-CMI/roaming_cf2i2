from coordinates import init_coord, change_coord, writeXYZ, atomLabels
import numpy as np
import itertools


if __name__ == "__main__":

    # grid_filename = "cf2i2_3d_grid1.xyz"
    grid_filename = "cf2i2_3d_grid2.xyz"

    (rCF_ref, rCI_ref, aFCF_ref, aICI_ref), xyz_ref = init_coord()
    r = rCI_ref
    t = 0
    p = 0
    writeXYZ(grid_filename, atomLabels, xyz_ref, comment=f"reference geometry: r={r}, theta={t}, phi={p}")

    # this grid is too large (cf2i2_3d_grid1.xyz)
    # r_grid = (np.linspace(1.0, 5.0, 30), np.linspace(5.0, 10.0, 20))
    # t_grid = (np.linspace(0.0, np.pi, 30), np.linspace(0.0, np.pi, 20))
    # p_grid = (np.linspace(0.0, 2*np.pi, 30), np.linspace(0.0, 2*np.pi, 20))

    # generate smaller grid (cf2i2_3d_grid2.xyz)
    r_grid = (np.linspace(1.0, 5.0, 20), np.linspace(5.0, 10.0, 10))
    t_grid = (np.linspace(0.0, np.pi, 20), np.linspace(0.0, np.pi, 10))
    p_grid = (np.linspace(0.0, 2*np.pi, 20), np.linspace(0.0, 2*np.pi, 10))

    rtp_grid = [[elem for elem in itertools.product(r, t, p)]\
                for r, t, p in zip(r_grid, t_grid, p_grid)]

    print(f"grid size: {[len(elem) for elem in rtp_grid]}")

    for rtp in rtp_grid:
        for (r, t, p) in rtp:
            xyz = change_coord(xyz_ref, r, t, p)
            writeXYZ(grid_filename, atomLabels, xyz, comment=f"r={r}, theta={t}, phi={p}", append=True)
