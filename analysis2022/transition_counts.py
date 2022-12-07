# import numpy as _np
from collections import defaultdict as _defaultdict
from copy import deepcopy as _deepcopy
# from scipy.special import loggamma as _loggamma
from constraint import Constraint_V2 as _Constraint

class _k_Transitions():
	"""
		_k_Transitions object tracks the  history-future 
		counts of a given length k. Parameter k determines 
		the lenght of the histories that are tracked.
	"""
	def __init__(self, nodes,  k, constraint = None):
		assert isinstance(k, int), "K has to be type int."
		
		self._nodes = set(nodes)
		self._n_nodes = len(self._nodes)
		assert self._n_nodes > 0, "At least one node is required."
		
		self._k = k
		
		if constraint is None:
			constraint = _Constraint(nodes = self._nodes)
		assert isinstance(constraint, _Constraint), "constraint has to be an instance of `Constraint` class"
		self._constraint = constraint

		self._counts = _defaultdict(lambda:_defaultdict(int))
	
	def rescale_counters(self, factor):
		for hist in self._counts:
			for future in self._counts[hist]:
				self._counts[hist][future] *= factor

	def add_transition(self, hist, future, number_of_instances=1):
		"""
			hist: touple of nodes, history of an event `future`
			future: element of nodes
			number_of_instance: float, represents the number of observations of the hist/future pair
			and_subhistories: boolean, lets user deside whether to compute subhistories or not,
				might speed up the calculation if you will never use the subhistories counts.
		"""
		assert self._constraint.check(hist, future), "Constraint not satisfied. " + str(hist) + " " + str(future)
		assert len(hist) == self._k, "History should have length k.\nHistory: {};\nk: {}".format(hist, self._k)
		self._counts[hist][future] += number_of_instances

	def relabel_nodes(self, node_to_label):
		"""
			takes counts of each path and transforms it, relabeling 
			each node with its group. this function will speed up the 
			fitting of the group dynamcis a lot. instead of computing 
			the counts all the time from paths, I can compute them 
			once, and then just transform them using different group maps
			
			returns new PathStatistics object

		"""
		assert False, "Not implemented yet, because self._constraint.relabel_nodes() does not exist."
		# check whether the map acts on the whole set of nodes
		assert set(node_to_label.keys()).issuperset(self._nodes), "Some nodes cannot be mapped with this map."

		# transform hist is a function that sends a history to a group space, 
		# using node_to_label on each node in the history.
		def transform_hist(hist): return tuple(node_to_label(node) for node in hist)
		
		def transform_future(future):
			return node_to_label[x] 
			

		# start fresh:
		new_constraint = self._constraint.relabel_nodes(node_to_label)
		new_transitions = _k_Transitions(k = self._k, 
										 nodes = set(node_to_label[node] for node in self._nodes),
										 constraint = self._constraint)

		for hist in self._counts[order]:
			for future in self._counts[order][hist]:
				new_transitions.add_transition(transform_hist(hist), transform_future(future), 
											   number_of_instances = self._counts[hist][future]) 
		return new_transitions

	def export_all_counts(self,):
		"""
			Return the number of counts of history being 
			followed by a future as a dictionary of dictionaries,
			for a given history order
		"""
		return _deepcopy(self._counts)

	def get_histories(self):
		""" returns observed histories"""
		return list(self._counts.keys())

	def get_history_counter(self, hist):
		"""
			returns the number of counts of history
		"""
		if hist in self._counts:
			return sum(self._counts[hist].values())
		return 0

	def get_transition_counter(self, hist, future):
		"""
			returns the number of counts of a future continuing hist
		"""
		if hist in self._counts:
			if future in self._counts[hist]:
				return self._counts[hist][future]
		return 0

	def get_transitions_given_hist(self, hist):
		"""
			returns a dictionary of transitions that are observed,
			and their counts. 
		"""
		if hist in self._counts:
			return dict(self._counts[hist])
		
	def get_k(self):
		"""
			returns max order.
		"""
		return self._k

	def get_nodes(self,):
		"""
			labels of all the elements mentioned in the path statistics
		"""
		return _deepcopy(self._nodes)

	def get_n_nodes(self,):
		"""
			returns the number of nodes
		"""
		return self._n_nodes

	def get_constraint(self,):
		"""
			returns the constraint, which is a constraint object
		"""
		return _deepcopy(self._constraint)

	def iter_hist(self,):
		"""
			This will work for now.
			Returns a list of histories.
			Hopefully it will be an actual iterator once.
		"""
		return iter(self._counts.keys())

	def iter_future_given_hist(self, hist):
		"""
			This will work for now.
			Returns a list of futures given hist.
			Hopefully it will be an actual iterator once.
		"""
		return self._counts[hist].keys()
	
	# def iter_hist_future_pairs(self,):
		# assert False, "Not implemented yet."

	def __add__(self, k_transitions):
		"""
			Combines two _k_transition objects.
			Assumes the constraints are the same, 
			but it will work as long as the transitions 
			that are added are respecting constraint of 
			this object.
		"""

		assert self._k == k_transitions.get_k(), "Transition order k is not the same."
		assert self._nodes == k_transitions.get_nodes(), "Sets of nodes not the same."
		
		combined = _deepcopy(self)

		#this should go with iter_hist_future_pairs
		for hist in k_transitions.iter_hist():
			for future in k_transitions.iter_future_given_hist(hist):
				combined.add_transition(hist, future, number_of_instances=k_transitions.get_transition_counter(hist, future))
		return combined

	def plot_observed_hist_counts_distribution(self, figsize = (10,10), filter_histories = None):
		import matplotlib.pyplot as plt
		import numpy as np

		histories = self.get_histories()
		if filter_histories is None:
			pass
		elif isinstance(filter_histories, int):
			histories = [hist for hist in histories if self.get_history_counter(hist)>filter_histories]
		else:
			histories = filter_histories(histories)
		
		y_of_x = _defaultdict(int)
		for history in histories:
			y_of_x[self.get_history_counter(history)] +=1
		
		fig, ax = plt.subplots(figsize = figsize)
		
		x = range(1,max(y_of_x.keys()))
		
		obs_y = np.array([y_of_x[xi] for xi in x])
		obs_y = obs_y/np.sum(obs_y)

		ax.bar(x, height = obs_y, color = (0,0.5,0), alpha = 1)
		ax.set_xticks(x, x)
		ax.set_xlabel("history counter")
		ax.set_ylabel("number of histories with counter x")
		ax.legend([], title = "{} observations".format(sum(y_of_x.values())))
		
		return ax

	def plot_transition_counts(self, figsize = (10,20), filter_histories = None):
		import matplotlib.pyplot as plt
		import numpy as np
		histories = self.get_histories()
		if filter_histories is None:
			pass
		elif isinstance(filter_histories, int):
			histories = [hist for hist in histories if self.get_history_counter(hist)>filter_histories]
		else:
			histories = filter_histories(histories)

		histories.sort(key = lambda hist: -self.get_history_counter(hist))
		
		n_hist = len(histories)
		
		# scale the numbers of plots per x and y axis
		# so that the plots look more or less square.
		if figsize[1]<figsize[0]:
			product_times_ratio = n_hist *figsize[1]/figsize[0]
			n = int(np.ceil(np.sqrt(product_times_ratio))) 
			m = int(np.ceil(int(n_hist/n))) 
		else:
			product_times_ratio = n_hist *figsize[0]/figsize[1]
			m = int(np.ceil(np.sqrt(product_times_ratio)) )
			n = int(np.ceil(n_hist/m) )
			
		ij_from_hist = {hist:(int(e/m), e%m) for e,hist in enumerate(histories)}
		fig, ax = plt.subplots(n, m, figsize = figsize)
		
		if n == 1:
			ax = [ax]
		if m == 1:
			ax = [[axi] for axi in ax]

		for hist in histories:
			
			i,j = ij_from_hist[hist]
			
			x = [str(xi) for xi in self._constraint.future_of(hist)]
			x.sort(key = lambda i:-self.get_transition_counter(hist, i))
			obs_y = np.array([self.get_transition_counter(hist, xi) for xi in x])
			obs_y = obs_y/np.sum(obs_y)

			ax[i][j].bar(x, height = obs_y, color = (0,0.5,0), alpha = 1)
			ax[i][j].set_xticks(x, x)
			ax[i][j].set_xticklabels(x, rotation = "vertical")
			title = "c(x|"+hist[0]
			for node in hist[1:]:
				title+= ", "+str(node)
			title += ")"
			ax[i][j].set_title(title)
			ax[i][j].legend([],title = "{} observations".format(self.get_history_counter(hist)))
		
		return ax

class Transitions(object):
	"""
		The class that injests the paths and counts the
		transitions. The transitions are then stored in 
		_k_Transitions objects.
	"""
	def __init__(self, nodes,  max_order, stationary = None, constraint = None):
		# self._nodes = set(nodes)
		# self._n_nodes = len(self._nodes)
		# assert self._n_nodes > 0, "At least one node is required."
		
		assert isinstance(max_order, int), "max_order has to be type int."
		assert max_order >= 0, "max_order has to be positive"
		self._K = max_order
		

		# defines at which transition we assume that the paths become stationary
		if stationary is None:
			self._stationary = float("inf")
		else:
			self._stationary = stationary

		if constraint is None:
			constraint = _Constraint(nodes = nodes)
		assert isinstance(constraint, _Constraint), "constraint has to be an instance of `Constraint` class"

		self._transitions = {k:_k_Transitions(nodes,
											  k,  
											  constraint = constraint) 
							 for k in range(max_order+1)}
	
	def rescale_counters(self, factor):
		for k in self._transitions:
			self._transitions[k].rescale_counters(factor)

	def add_transition(self, history, future, number_of_instances = 1):
		"""
			adds a transition.
		"""
		self._transitions[len(history)].add_transition(history, future, number_of_instances=number_of_instances)

	def add_path(self, path, number_of_instances = 1):
		"""
			All histories are saved only insofar they occur without
			anything before it, but the one of order K are saved also
			if there is something before them.

			Example:

				K = 3
				path = (1,2,3,4,5,6,7)
				histories and futures extracted are going to be:
					() 		: 1
					(1)		: 2
					(1,2)	: 3
					(1,2,3)	: 4
					(2,3,4)	: 5
					(3,4,5)	: 6
					(4,5,6)	: 7
		"""
		for e, future in enumerate(path):
			if e<=self._K:
				history = path[:e]
				order = e
			else:
				history = path[(e-self._K):e]
				assert len(history) == self._K, "Path {path};\nhist: {hist};\nK: {K};\ne: {e};".format(path = path, hist = history, K = self._K, e = e)
				order = self._K

			self._transitions[order].add_transition(history, future, number_of_instances = number_of_instances)
			
			# allws self.stationarity<self._K
			# in other words, it allows sub-paths
			if e >= self._stationary:
				# start at 1, because I added all regular paths anyway
				# until order + 1 because I want to pick them all up, including empty history ()
				for i in range(1,order+1):
					shorter_hist = history[i:]
					assert len(shorter_hist) == order-i, "shorter_hist {}, e {}, i {}".format(shorter_hist, order, i)
					self._transitions[order-i].add_transition(shorter_hist, future, number_of_instances = number_of_instances)
				
	def add_paths(self, paths):
		"""
			injests paths object.
		"""
		# this should go with paths.iter_paths() but it would be too slow now.
		for length in paths.lengths():
			for path in paths.paths_of_length(length):
				self.add_path(path, number_of_instances = paths.path_counter(path))

	def orders(self,):
		"""
			returns the k's that are used here
		"""
		return list(range(self._K + 1))

	def iter_orders(self,):
		""" iterates through history lengths."""
		return iter(self._transitions.keys())

	def iter_hist_of_order(self,order):
		"""iterate through all histories of a specific order"""
		return self._transitions[order].iter_hist()

	def transitions_of_order(self,k):
		"""
			returns transitions of order k
		"""
		return self._transitions[k]
		
	def get_histories(self):
		histories = []
		for k in self._transitions:
			histories = histories + self._transitions[k].get_histories()
		return histories

	def get_history_counter(self, hist):
		"""
			returns the number of counts of history
		"""
		order =	len(hist)
		assert self.is_order(order), "Not an allowed order: {}, max_order = {}".format(order, self._K)
		return self._transitions[order].get_history_counter(hist)
	
	def get_transition_counter(self, hist, future):
		"""
			returns the number of counts of a future continuing hist
		"""
		order =	len(hist) 
		assert self.is_order(order), "Not an allowed order: {}, max_order = {}".format(order, self._K)
		return self._transitions[order].get_transition_counter(hist, future)

	def get_transitions_given_hist(self,hist):
		"""returns dictionary[future] = count"""
		order =	len(hist)
		assert self.is_order(order), "Not an allowed order: {}, max_order = {}".format(order, self._K)
		return self._transitions[len(hist)].get_transitions_given_hist(hist)

	def get_max_order(self):
		"""
			Maximum order of transitions that is saved.
			
		"""
		return self._K

	def get_nodes(self,):
		"""
			labels of all the elements mentioned in the path statistics
		"""
		return _deepcopy(self._transitions[0]._nodes)

	def get_n_nodes(self,):
		"""
			returns the number of nodes
		"""
		return self._transitions[0]._n_nodes

	def get_constraint(self,):
		"""
			returns the constraint, which is a constraint object
		"""
		return _deepcopy(self._transitions[0]._constraint)

	def is_order(self, order):
		"""checks whether something is a positive integer number smaller than K+1"""
		return isinstance(order, int) and order <= self._K and order >= 0

	def relabel_nodes(self, node_to_label):
		"""
			takes counts of each path and transforms it, relabeling 
			each node with its group. this function will speed up the 
			fitting of the group dynamcis a lot. instead of computing 
			the counts all the time from paths, I can compute them 
			once, and then just transform them using different group maps
			
			returns new PathStatistics object

		"""
		for k in self._transitions:
			self._transitions[k].relabel_nodes(node_to_label)

	def __add__(self, transitions):
		"""
			Combines two transition objects.
			Assumes the constraints are the same.
		"""
		combined = _deepcopy(self)
		for k in transitions.orders():
			combined._transitions[k] += transitions.transitions_of_order(k)
		return combined
	