from numpy.typing import NDArray
from typing import Tuple, List
import numpy as np 
from numpy.linalg import norm
from numpy import radians as rad
from scipy.spatial.transform import Rotation as R
def init_coord():
    FCF = 109.5  #degrees
    ICI = 112.5  #degrees
    CF = 1.34
    CI = 2.15
    atoms = ['C', 'F', 'F', 'I', 'I']
    coord = [[0,0,0],
            [-CF* np.sin(rad(FCF/2)), 0, -CF*np.cos(rad(FCF/2))],
            [CF*np.sin(rad(FCF/2)), 0, -CF*np.cos(rad(FCF/2))],
            [0, CI*np.sin(rad(ICI/2)), CI*np.cos(rad(ICI/2))],
            [0, -CI*np.sin(rad(ICI/2)), CI*np.cos(rad(ICI/2))]]
    return np.array(coord)
    
def change_coord(coord, r, theta, phi):
    vec1, vec2 = np.array(coord[1]), np.array(coord[2])
    zaxis = np.abs((vec1 + vec2)/norm(vec1 + vec2))
    xaxis = np.abs((vec1 - vec2)/norm(vec1 - vec2))
    yaxis = np.abs(np.cross(vec1, vec2)/norm(np.cross(vec1, vec2)))
    coord[4] = [0, -r*np.sin(rad(theta)), r*np.cos(rad(theta))]
    rot_vec = zaxis * rad(phi)
    rotation = R.from_rotvec(rot_vec)
    coord[4] = rotation.apply(coord[4])
    return coord


def writeXYZ(fileName: str, atomLabels: List[str], atomCoords: NDArray[np.float_],
             comment: str = "", append: bool = False
             ) -> None:
    with open(fileName, 'a' if append else 'w') as fl:
        fl.write(str(len(atomCoords)) + "\n")
        fl.write(comment + "\n")
        fl.write("\n".join(atom + "     " + "   ".join("%18.12f"%x for x in coords)
                     for atom, coords in zip(atomLabels, atomCoords) ))
        fl.write("\n")


def test_coords(atomLabels):
    coo0 = init_coord()
    r0 = 1.44
    theta0 = 0
    phi0 = 0
    writeXYZ('test_coords.xyz', atomLabels, coo0)

    r = np.linspace(r0, 10, 100)
    theta = [theta0]*len(r)
    phi = [phi0]*len(r)
    for x in zip(r, theta, phi):
        coo = change_coord(coo0, *x)
        writeXYZ('test_coords.xyz', atomLabels, coo, append=True)

    theta = np.linspace(0, np.pi, 100)
    r = [r0]*len(theta)
    phi = [phi0]*len(theta)
    for x in zip(r, theta, phi):
        coo = change_coord(coo0, *x)
        writeXYZ('test_coords.xyz', atomLabels, coo, append=True)

    phi = np.linspace(0, 2*np.pi, 100)
    r = [r0]*len(phi)
    theta = [theta0]*len(phi)
    for x in zip(r, theta, phi):
        coo = change_coord(coo0, *x)
        writeXYZ('test_coords.xyz', atomLabels, coo, append=True)


if __name__ == "__main__":
    atomLabels = ['C', 'F', 'F', 'I', 'I']
    test_coords(atomLabels)
    # new_coord = change_coord(coord = init_coord(), r = 1.44, theta = 120/2 , phi = 90)
    # writeXYZ(fileName = 'cf2i2_1.xyz', atomLabels = ['C', 'F', 'F', 'I', 'I'], atomCoords = new_coord)
    
