import numpy as np
from scipy import interpolate 


def readPDBcoords(fname, resids):
	coords = {}
	sel = 'CA'
	with open(fname) as f:
		for line in f:
			if 'ATOM' in line:
				lines = line.split()
				resid = int(lines[4])-1
				atom_type = lines[2]
				if resid in resids and atom_type == sel:
					coords[resid] = [float(lines[5]),
									float(lines[6]),
									float(lines[7])]
	return coords

def import_(name="paths.txt"):
	paths = []
	with open(name) as f:
		for line in f:
			lines = line.split()
			paths.append([int(x) for x in lines[:-1]])
	spaths = np.hstack(paths)
	nodes = np.unique(spaths)
	return paths, nodes
	
def intrpl_(path, coords, smoothness=0.01):
	x_vals = []
	y_vals = []
	z_vals = []
	for res in path:
		coor = coords[res]
		x_vals.append(coor[0])
		y_vals.append(coor[1])
		z_vals.append(coor[2])
	degree = len(x_vals) - 1
	if degree > 3:
		degree = 3
	tck, _ = interpolate.splprep([x_vals, y_vals, z_vals], s=0, k=degree)
	unew = np.arange(0, 1.01, smoothness)
	out = interpolate.splev(unew, tck)
	return out

def draw_paths(paths, coords, smoothness=0.01, output="test.bild"):
	t = open(output, "w")
	for p in paths:
		ipath = intrpl_(p, coords)
		for c in range(len(ipath[0]) - 1):
			x1 = ipath[0][c]
			y1 = ipath[1][c]
			z1 = ipath[2][c]
			x2 = ipath[0][c+1]
			y2 = ipath[1][c+1]
			z2 = ipath[2][c+1]
			t.write(".color red\n")
			clt = ".cylinder {} {} {} {} {} {} {}\n".format(round(x1,3),
							round(y1,3), round(z1,3), round(x2,3), round(y2,3),
							round(z2,3), 0.1)
			t.write(clt)
	t.close()

paths, nodes = import_("paths-PIC.txt")
coords = readPDBcoords("network-all.pdb", nodes)
draw_paths(paths, coords, smoothness=0.05)
