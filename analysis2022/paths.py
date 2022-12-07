from collections import defaultdict as _defaultdict
from numpy.random import binomial as _binomial
from numpy.random import random as _random



class Paths(object):
	"""docstring for Paths"""
	def __init__(self):
		self._paths =_defaultdict(lambda: _defaultdict(int))
		self._number_paths_length = _defaultdict(int)
		self._tot_number_of_paths = 0
		self._nodes = set()

	def add_path(self, path, number_of_instances_of_path=1):
		if number_of_instances_of_path != 0:
			self._paths[len(path)-1][path] += number_of_instances_of_path
			self._tot_number_of_paths += number_of_instances_of_path
			self._number_paths_length[len(path)-1] += number_of_instances_of_path
			self._nodes.update(path)

	def add_paths(self, paths):
		for l in paths.lengths():
			for p in paths.paths_of_length(l):
				self.add_path(p, number_of_instances_of_path = paths.path_counter(p))

	def update_from_ngram(self, filename, splitter = ",", last_entry_frequency = False, print_one = False):
		with open(filename, "r") as f:
			for line in f:
				path = tuple(line[:-1].split(splitter))
				if last_entry_frequency:
					self.add_path(path[:-1], float(path[-1]))
				else:
					self.add_path(path)
				if print_one == True:
					if last_entry_frequency:
						print(path[:-1], path[-1])
					else:
						print(path, 1)
					print(line)
					print_one = False

	def update_from_list(self, list_of_paths):
		for path in list_of_paths:
			self.add_path(path)

	def total_number(self):
		return self._tot_number_of_paths
		
	def path_counter(self, path):
		return self._paths[len(path)-1][path]

	def lengths(self):
		return set(self._paths.keys())

	def paths_of_length(self, path_length):
		return set(self._paths[path_length].keys())

	def iter_paths(self):
		return {path for path in self._paths[length] for length in self._paths}
	
	def print_stats(self):
		print("total number of paths: {}".format(self._tot_number_of_paths))
		print("length | total #,  # uniques | min #, max # |")
		print("=============================================")
		for l in range(max( self._paths)+1):
			if l in self._paths:
				print("{leng:6d} | {tot:7.2f},  {uniq:9.2f} | {mini:5.2f}, {maxi:5.2f} |".format(
					leng = l, tot = sum(self._paths[l].values()), uniq = len(self._paths[l].keys()), mini = min(self._paths[l].values()), maxi = max(self._paths[l].values())))
	
	def train_test_split(self, train_percentage=0.8):
		"""
			Splits the Paths object into two Paths objects, at random.
			Train set will on average get the share of paths equal to 
			`train_percentage`.
			returns two Paths objects
		"""
		train = Paths()
		test = Paths()
		for l in self.lengths():
			for path in self.paths_of_length(l):
				n = self.path_counter(path)
				i = _binomial(n, train_percentage)
				if i>0:
					train.add_path(path, i)
				if n-i>0:
					test.add_path(path, n-i)
		return train, test

	def remove_path(self,path, n=1):
		leng = len(path)-1
		assert leng in self._paths
		assert path in self._paths[leng]
		assert n < self._paths[leng][path]

		self._paths[len(path)-1][path] -= n
		self._tot_number_of_paths -= n
		self._number_paths_length[leng] -= n
		

		

	def get_visited_nodes(self):
		return self._nodes
