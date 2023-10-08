from numpy.typing import NDArray
from typing import Tuple, List
import numpy as np 
from numpy.linalg import norm
from numpy import radians as rad
from scipy.spatial.transform import Rotation as R
import random
import itertools

atomLabels = ['C', 'F', 'F', 'I', 'I']

def init_coord():
    global FCF, ICI, CF, CI
    FCF = 109.5  #degrees
    ICI = 112.5  #degrees
    CF = 1.34 #angstrom
    CI = 2.15 #angstrom
    atoms = ['C', 'F', 'F', 'I', 'I']
    coord = [[0,0,0],
            [-CF* np.sin(rad(FCF/2)), 0, -CF*np.cos(rad(FCF/2))],
            [CF*np.sin(rad(FCF/2)), 0, -CF*np.cos(rad(FCF/2))],
            [0, CI*np.sin(rad(ICI/2)), CI*np.cos(rad(ICI/2))],
            [0, -CI*np.sin(rad(ICI/2)), CI*np.cos(rad(ICI/2))]]
    return (CF, CI, FCF, ICI), np.array(coord)
 

def init_axes():
    _, coord = init_coord()
    vec1, vec2 = np.array(coord[1]), np.array(coord[2])
    zaxis = np.abs((vec1 + vec2)/norm(vec1 + vec2))
    xaxis = np.abs((vec1 - vec2)/norm(vec1 - vec2))
    yaxis = np.abs(np.cross(vec1, vec2)/norm(np.cross(vec1, vec2)))
    return(xaxis, yaxis, zaxis)
 

def change_coord(coord, r, theta, phi):
    x, y, z = init_axes()
    
    if theta == rad(ICI/2) and phi == 0 :
        theta = rad(ICI/2)
        coord[4] = [0, -r*np.sin(theta), r*np.cos(theta)]
    elif theta != rad(ICI/2) and phi == 0 :
        coord[4] = [0, -r*np.sin(theta), r*np.cos(theta)]
    elif phi != 0 and theta == rad(ICI/2):
        rot_vec = z * phi
        rotation = R.from_rotvec(rot_vec)
        theta = rad(ICI/2)
        coord[4] = [0, -r*np.sin(theta), r*np.cos(theta)]
        coord[4] = rotation.apply(coord[4])
    else:
        rot_vec = z * phi
        rotation = R.from_rotvec(rot_vec)
        coord[4] = [0, -r*np.sin(theta), r*np.cos(theta)]
        coord[4] = rotation.apply(coord[4])
    return coord


def rotate_frag(coord, r, angle, axis):
    x, y, z = init_axes()
    axes = {'x': x, 'y':y, 'z':z}
    assert r > CI
    
    theta = rad(ICI/2)
    coord[4] = [0, -r*np.sin(theta), r*np.cos(theta)]
    rot_vec = axes[axis] * angle
    rotation = R.from_rotvec(rot_vec)
    coord[0 : 4] = rotation.apply(coord[ 0 : 4 ])
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
    _, coo0 = init_coord()
    r0 = 2.15   #  distance of carbon iodine
    theta0 = 0
    phi0 = 0
    writeXYZ('test_coords.xyz', atomLabels, coo0)
    writeXYZ('test_Frag_coords.xyz', atomLabels, coo0)

    r = np.linspace(r0, 10, 100)
    theta = [theta0]*len(r)    # theta and phi is given in radians
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
        
    rot_angle = np.random.random(300)
    rot_angle *= np.pi
    r = [3]*len(rot_angle)
    axes = [random.choice(['x', 'y', 'z'])]*len(rot_angle)
    for x in zip(r, rot_angle, axes):
        coo = rotate_frag(coo0, *x)
        writeXYZ('test_Frag_coords.xyz', atomLabels, coo, append=True)
        

def test_coords2(atomLabels):
    (rCF, rCI, aFCF, aICI), coo0 = init_coord()
    r0 = rCI * 2
    writeXYZ('test_coords2.xyz', atomLabels, coo0)

    theta = np.linspace(0.001, np.pi-0.001, 10)
    phi = np.linspace(0, 2*np.pi, 100)
    theta_phi = [elem for elem in itertools.product(theta, phi)]
    for (theta, phi) in theta_phi:
        print(theta, phi)
        coo = change_coord(coo0, r0, theta, phi)
        writeXYZ('test_coords2.xyz', atomLabels, coo, append=True)

            
            
if __name__ == "__main__":
    atomLabels = ['C', 'F', 'F', 'I', 'I']
    # test_coords(atomLabels)
    test_coords2(atomLabels)
    # new_coord = change_coord(coord = init_coord(), r = 1.44, theta = 120/2 , phi = 90)
    # writeXYZ(fileName = 'cf2i2_1.xyz', atomLabels = ['C', 'F', 'F', 'I', 'I'], atomCoords = new_coord)
    
