import numpy as np
from heapq import heappush, heappop
import networkx as nx
import collections
from operator import itemgetter
from itertools import count, islice
import time
import argparse

class PathBuffer(object):
	"""
		This is a heap-like Priority Queue, modified for paths
		determined with Yen's algorithm. This data structure
		is similar to the one used in Networkx
	"""
	def __init__(self):
		self.paths = set()
		self.sortedpaths = list()
		self.counter = count()
		
	def __len__(self):
		return len(self.sortedpaths)
		
	def push(self, cost, path):
		hashable_path = tuple(path)
		if hashable_path not in self.paths:
			heappush(self.sortedpaths, (cost, next(self.counter), path))
			self.paths.add(hashable_path)
			
	def pop(self):
		(cost, num, path) = heappop(self.sortedpaths)
		hashable_path = tuple(path)
		self.paths.remove(hashable_path)
		return path

class Subgraph(PathBuffer):
	def __init__(self, A, s, t):
		self.A = []
		with open(A) as f:
			for line in f:
				lines = line.split()
				self.A.append([float(x) for x in lines])
		self.A = np.asarray(self.A, dtype=np.float64)
		self.G = nx.from_numpy_matrix(self.A)
		self.s = s
		self.t = t
		
	def dijkstra(self, A, s, t, ignore_nodes=None, ignore_edges=None):
		"""
		This is a variant of Dijkstra's algorithm that uses a
		bidirectional search to find the shortest path. The code
		structure is similar to that of Networkx with the exception
		that it has beeen modified to work explicity on Numpy arrays.
		Additionally, the parameter ignore_edges_nodes allows one to 
		employ this algorithm with Yen's algorithm. 
		"""
		if s == t:
			return (0, [s])
		
		push = heappush
		pop = heappop
		dists = [{}, {}]
		paths = [{s: [s]}, {t: [t]}]
		fringe = [[], []]
		seen = [{s: 0}, {t: 0}]
		c = count()
		push(fringe[0], (0, next(c), s))
		push(fringe[1], (0, next(c), t))
		finalpath = []
		dir = 1
		while fringe[0] and fringe[1]:
			# choose direction, 0 is forward
			dir = 1 - dir
			(d, _, current_node) = pop(fringe[dir])
			neighbors = np.where(A[:,current_node] != 0.0)[0]
			if ignore_nodes:
				#print current_node, ignore_edges_nodes, neighbors
				def filter_nodes(neigh):
					for n in neigh:
						if n not in ignore_nodes:
							yield n
				neighbors = filter_nodes(neighbors)
				neighbors = list(neighbors)
			if ignore_edges:
				def filter_edges(neigh):
					for n in neigh:
						if (current_node, n) not in ignore_edges and (n, current_node) not in ignore_edges:
							yield n
				neighbors = filter_edges(neighbors)
				neighbors = list(neighbors)
			if current_node in dists[dir]:
				continue
			dists[dir][current_node] = d
			if current_node in dists[1-dir]:
				return (finaldist, finalpath) # We are done
			if len(neighbors) == 0:
				raise KeyError
			for n in neighbors:
				cost = A[current_node][n]
				cn_dist = dists[dir][current_node] + cost
				if n in dists[dir]:
					if cn_dist < dists[dir][n]:
						raise ValueError("Negative Weights?")
				elif n not in seen[dir] or cn_dist < seen[dir][n]:
					seen[dir][n] = cn_dist
					push(fringe[dir], (cn_dist, next(c), n))
					paths[dir][n] = paths[dir][current_node] + [n]
					if n in seen[0] and n in seen[1]:
						totaldist = seen[0][n] + seen[1][n]
						if finalpath == [] or finaldist > totaldist:
							finaldist = totaldist 
							revpath = paths[1][n][:]
							revpath.reverse()
							finalpath = paths[0][n] + revpath[1:]
		if finalpath == []:
			raise KeyError

	def multisource_dijkstra(self, A, sources, t):
		paths = {source: [source] for source in sources}
		push = heappush
		pop = heappop
		dist = {}
		seen = {}
		c = count()
		fringe = []
		for source in sources:
			seen[source] = 0
			push(fringe, (0, next(c), source))
		while fringe:
			(d, _, v) = pop(fringe)
			if v in dist:
				continue
			dist[v] = d
			if v == target:
				break
			neighbors = np.where(A[:,v] != 0.0)[0]
			for n in neighbors:
				cost = A[v][n]
				vn_dist = dist[v] + cost
				if n in dist:
					if vn_dist < dist[n]:
						raise ValueError("Negative Weights?")
				elif n not in seen or vn_dist < seen[n]:
					seen[n] = vn_dist
					push(fringe, (vn_dist, next(c), n))
					paths[n] = paths[v] + [n]
		return (dist[target], paths[target])
				
	
	def find_opt_path(self):
		if len(self.s) > 1:
			opt_path_len, opt_path = self.multisource_dijkstra(self.A, self.s, self.t)
		else:
			opt_path_len, opt_path = self.dijkstra(self.A, self.s[0], self.t)
		self.opt_path = opt_path
		self.opt_path_len = opt_path_len
		self.s = self.opt_path[0]
	
	def _neighbors(self, node):
		neighbors = np.where(self.A[:,node] != 0.0)[0]
		neighbors = [x for x in neighbors if x != node]
		return neighbors
	
	def get_all_neighbors(self, path):
		all_neighbors = []
		for node in path:
			all_neighbors.append(self._neighbors(node))
		all_neighbors = np.hstack(all_neighbors)
		return list(set(all_neighbors))

	def get_level(self):
		if self.level == 1:
			return self.get_all_neighbors(self.opt_path)
		else:
			master = []
			master.append(self.get_all_neighbors(self.opt_path))
			count = 0
			while count < level:
				master.append(self.get_all_neighbors(master[count]))
				count += 1
			master = np.hstack(master)
			return list(set(master))
	
	def create_subgraph(self, level=1):
		print("Expanding {} neighbor(s) away from optimal path".format(level))
		self.level = level
		neighborslist = self.get_level()
		paths = [self.opt_path]
		for node in neighborslist:
			node_path = nx.all_shortest_paths(self.G, source=self.s,
							target=node,
							weight='weight')
			target_path = nx.all_shortest_paths(self.G, source=node,
							target=self.t,
							weight='weight')
			for p in node_path:
				for j in target_path:
					try:
						if j[1] not in p: # Avoids paths that fold back to the previous node
							path = p[:-1] + j
							if path not in paths:
								paths.append(path)
					except IndexError:
						pass
		nodelist = np.hstack(paths)
		nodelist = list(set(nodelist))
		self.S = self.G.subgraph(nodelist)
		unmapped_paths = paths
		mapping = {list(self.S.nodes())[n]: n for n in range(len(nodelist))} # Mapping from G to S
		reversemapping = {n: list(self.S.nodes())[n] for n in range(len(nodelist))} # Mapping from S to G
		SadjM = nx.to_numpy_array(self.S, dtype=np.float64)
		rtarget = mapping[self.t]
		return SadjM, mapping, reversemapping, unmapped_paths, rtarget

	def map_paths_to_graph(self, paths, mapping):
		mapped_paths = []
		for p in paths:
			mapped_paths.append([mapping[x] for x in p])
		return mapped_paths
		
	def _dist(self, A, path):
		return np.sum([A[path[x]][path[x+1]] for x in range(len(path)-1)])
	
	def find_paths_from_graph(self, S, k, t, somepath):
		"""
		This is Yen's algorithm for finding k-shortest paths.
		It uses Dijkstra's algorithm to find the shortest path
		between a spur and the target. Additionally, we ignore
		nodes and edges from already discovered paths.
		"""
		print("Searching for {} paths".format(k))
		found_paths = [somepath]
		pathlist = [somepath]
		prev_path = somepath
		Qpath = PathBuffer()
		EXIT = False
		while len(found_paths) < k:
			#roots = [prev_path[:x] for x in range(1, len(prev_path))]
			ignore_nodes = set()
			ignore_edges = set()
			for i in range(1, len(prev_path)):
				root = prev_path[:i]	
				root_len = self._dist(S, root)
				for p in pathlist:
					if p[:i] == root:
						ignore_edges.add((p[i-1], p[i]))
				try:
					length, spur = self.dijkstra(S, root[-1], t,
								ignore_nodes=ignore_nodes,
								ignore_edges=ignore_edges)
					path = root[:-1] + spur
					path_len = root_len + length
					Qpath.push(path_len, path)
				except KeyError:
					pass
				except TypeError: # We've run out of nodes and need to expand further away
					EXIT = True 
					break
				ignore_nodes.add(root[-1])
			if EXIT == True: # This breaks the main loop if we've run out of nodes
				break 
			if Qpath:
				new_path = Qpath.pop()
				pathlist.append(new_path)
				if new_path not in found_paths:
					found_paths.append(new_path)
				prev_path = new_path
			else:
				break
		return found_paths
		
	def output_paths(self, paths, npathstowrite, outname="paths.txt"):
		# First, we need to ensure that our paths are sorted by length
		unsorted = [p + [self._dist(self.A, p)] for p in paths]
		sorted_paths = sorted(unsorted, key=itemgetter(-1))
		t = open(outname, "w")
		c = 0
		for p in sorted_paths:
			s = " ".join([str(x) for x in p[:-1]])
			t.write(s+" {}\n".format(p[-1]))
			c += 1
			if c > npathstowrite:
				break
		t.close()


# Main execution			
parser = argparse.ArgumentParser()
parser.add_argument("-i", dest="input", type=str, default="network.dat", help="path to adjacency matrix input file, default='network.dat'")
parser.add_argument("-s", dest="source", type=str, default="1", help="source node (can be multiple sources), default=1")
parser.add_argument("-t", dest="target", type=int, default=1, help="target node, default=1")
parser.add_argument("-l", dest="level", type=int, default=1, help="level of neighbors to expand to, default=1")
parser.add_argument("-k", dest="npaths", type=int, default=1000, help="number of paths to find, default=1000")
parser.add_argument("-o", dest="output", type=str, default="paths.txt", help="output file name, default='paths.txt'")
parser.add_argument("-m", dest="maxlevel", type=int, default=5, help="Max level to expand to (default=5), this should only be changed if the calculation fails")
args = parser.parse_args()
adjM = args.input
ofile = args.output
try:
	source = [int(x)-1 for x in args.source.split(",")]
except:
	source = [int(args.source)-1]
target = args.target-1
npaths = args.npaths
level = args.level
maxlevel = args.maxlevel


if source == target:
	raise ValueError("Source {} == target {}, exiting".format(source, target))

	
S = Subgraph(adjM, source, target)
start = time.time()
S.find_opt_path()
source = S.opt_path[0]
print("Calculating suboptimal paths between {} and {}\n".format(source+1, target+1))
SadjM, mapping, reversemapping, unmapped_paths, mapped_t = S.create_subgraph(level)
paths = S.find_paths_from_graph(SadjM, npaths, mapped_t, [mapping[x] for x in unmapped_paths[0]]) 
mapped_paths = []
for p in paths:
	mapped_paths.append([reversemapping[x] for x in p])
print("Done!\nResults written to {}".format(ofile))
end = time.time() - start
print("Total time to find {} good paths is {} seconds".format(len(paths), end))
S.output_paths(mapped_paths, npaths, outname=ofile)

